#!/usr/bin/env python3
"""Run isolated LTE parameter sweeps from a saved index fit.

Examples:
    ./lte_sweep.py amo
    ./lte_sweep.py amo --target 0.9 --timeout 300
    ./lte_sweep.py amo --drop-low-modulation --year-offsets=-.001,0,.001
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import time
import uuid


ROOT = Path(__file__).resolve().parent
RUNNER = ROOT / "lte_run.py"


def parse_offsets(value: str) -> list[float]:
    offsets = [float(item) for item in value.split(",") if item.strip()]
    if not offsets:
        raise argparse.ArgumentTypeError("at least one comma-separated offset is required")
    return sorted(set(offsets), key=abs)


def correlation(results_path: Path) -> float:
    model: list[float] = []
    data: list[float] = []
    with results_path.open(newline="", encoding="utf-8") as results_file:
        for row in csv.reader(results_file):
            if len(row) < 3:
                continue
            try:
                model.append(float(row[1]))
                data.append(float(row[2]))
            except ValueError:
                continue

    if len(model) < 2:
        raise ValueError(f"{results_path} has fewer than two numeric model/data rows")
    model_mean = sum(model) / len(model)
    data_mean = sum(data) / len(data)
    covariance = sum(
        (model_value - model_mean) * (data_value - data_mean)
        for model_value, data_value in zip(model, data)
    )
    model_variance = sum((value - model_mean) ** 2 for value in model)
    data_variance = sum((value - data_mean) ** 2 for value in data)
    if model_variance == 0.0 or data_variance == 0.0:
        raise ValueError(f"{results_path} contains a constant series")
    return covariance / math.sqrt(model_variance * data_variance)


def set_response_option(path: Path, name: str, value: int) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    prefix = f"{name} "
    for index, line in enumerate(lines):
        if line.lstrip().startswith(prefix):
            lines[index] = f"{name:<30}{value}"
            break
    else:
        lines.append(f"{name:<30}{value}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def prepare_trial(
    source_dir: Path,
    working_dir: Path,
    trial_index: str,
    year: float,
    init: float,
    drop_low_modulation: bool,
) -> None:
    working_dir.mkdir()
    for name in (f"{source_dir.name}.dat", "lt.exe.p", "lt.exe.resp"):
        shutil.copy2(source_dir / name, working_dir / name)
    shutil.move(working_dir / f"{source_dir.name}.dat", working_dir / f"{trial_index}.dat")

    parameters_path = working_dir / "lt.exe.p"
    parameters = json.loads(parameters_path.read_text(encoding="utf-8"))
    parameters["year"] = year
    parameters["init"] = init
    if drop_low_modulation:
        modulations = parameters.get("ltep")
        if not isinstance(modulations, list) or len(modulations) < 2:
            raise ValueError("lt.exe.p must contain at least two ltep values")
        low_index = min(range(len(modulations)), key=lambda index: abs(float(modulations[index])))
        modulations.pop(low_index)
        set_response_option(working_dir / "lt.exe.resp", "NM", len(modulations))
    parameters_path.write_text(json.dumps(parameters, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("index", help="source index directory, for example amo")
    parser.add_argument("--cv", nargs=2, metavar=("BEGIN", "END"), default=("1880", "1885"))
    parser.add_argument("--timeout", type=float, default=300.0)
    parser.add_argument("--target", type=float, default=0.9, help="stop at this full-series correlation")
    parser.add_argument("--max-trials", type=int, default=25)
    parser.add_argument(
        "--year-offsets",
        type=parse_offsets,
        default=parse_offsets("-0.001,-0.0005,0,0.0005,0.001"),
        help="comma-separated offsets around saved year (default: -0.001,-0.0005,0,0.0005,0.001)",
    )
    parser.add_argument(
        "--init-offsets",
        type=parse_offsets,
        default=parse_offsets("-1,-0.5,0,0.5,1"),
        help="comma-separated offsets around saved init (default: -1,-0.5,0,0.5,1)",
    )
    parser.add_argument(
        "--drop-low-modulation",
        action="store_true",
        help="remove the smallest-magnitude ltep entry and update NM for every trial",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="directory for preserved trial outputs (default: lte_sweeps/<index>-<timestamp>)",
    )
    args = parser.parse_args()

    if not RUNNER.is_file():
        raise FileNotFoundError(f"Missing runner: {RUNNER}")
    source_dir = ROOT / args.index
    required = (source_dir / f"{args.index}.dat", source_dir / "lt.exe.p", source_dir / "lt.exe.resp")
    if not source_dir.is_dir() or not all(path.is_file() for path in required):
        raise FileNotFoundError(f"{source_dir} must contain {args.index}.dat, lt.exe.p, and lt.exe.resp")
    if args.timeout <= 0.0 or args.max_trials < 1:
        raise ValueError("--timeout and --max-trials must be positive")

    timestamp = time.strftime("%Y%m%d-%H%M%S")
    output_dir = args.output_dir or ROOT / "lte_sweeps" / f"{args.index}-{timestamp}"
    output_dir.mkdir(parents=True, exist_ok=False)
    parameters = json.loads((source_dir / "lt.exe.p").read_text(encoding="utf-8"))
    seed_year = float(parameters["year"])
    seed_init = float(parameters["init"])
    candidates = itertools.islice(
        itertools.product(args.year_offsets, args.init_offsets), args.max_trials
    )
    summary_path = output_dir / "summary.csv"
    best_score = -math.inf

    with summary_path.open("w", newline="", encoding="utf-8") as summary_file:
        summary = csv.DictWriter(
            summary_file,
            fieldnames=("trial", "year", "init", "correlation", "returncode"),
        )
        summary.writeheader()
        for number, (year_offset, init_offset) in enumerate(candidates, start=1):
            trial_index = f".lte-sweep-{uuid.uuid4().hex}"
            working_dir = ROOT / trial_index
            trial_dir = output_dir / f"trial-{number:03d}"
            year = seed_year + year_offset
            init = seed_init + init_offset
            try:
                prepare_trial(
                    source_dir,
                    working_dir,
                    trial_index,
                    year,
                    init,
                    args.drop_low_modulation,
                )
                command = [
                    sys.executable,
                    str(RUNNER),
                    trial_index,
                    "--cv",
                    *args.cv,
                    "--timeout",
                    str(args.timeout),
                ]
                with (working_dir / "run.log").open("w", encoding="utf-8") as log_file:
                    completed = subprocess.run(
                        command,
                        cwd=ROOT,
                        stdout=log_file,
                        stderr=subprocess.STDOUT,
                        check=False,
                    )
                score = (
                    correlation(working_dir / "lte_results.csv")
                    if completed.returncode == 0
                    else math.nan
                )
                shutil.move(str(working_dir), str(trial_dir))
                summary.writerow(
                    {
                        "trial": number,
                        "year": year,
                        "init": init,
                        "correlation": score,
                        "returncode": completed.returncode,
                    }
                )
                summary_file.flush()
                print(f"trial {number:03d}: correlation={score:.6f} year={year:.9g} init={init:.9g}")
                if math.isfinite(score) and score > best_score:
                    best_score = score
                if best_score >= args.target:
                    print(f"target reached: {best_score:.6f}")
                    break
            finally:
                if working_dir.exists():
                    shutil.rmtree(working_dir)

    print(f"results: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
