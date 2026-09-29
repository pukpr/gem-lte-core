#!/usr/bin/env python3
"""Cheap, no-search "take the best of what already exists" pass
(2026-09-29, per user request): for every grid cell, compare the
Feb2026 and Sep2026 results directly from their own lte_results.csv
(Pearson CC between Model and Data columns -- an apples-to-apples,
zero-Ada-invocation comparison). Where Feb2026 is meaningfully better,
trial-copy its lt.exe.p/lt.exe.resp/.dat into the Sep2026 directory and
gate the swap on two fast checks before keeping it:
  1. Manifold/backbone range (sweep.backbone_ok, instant -- just reads
     lt.exe.p's ltep[nm-1]).
  2. A live dLOD reading via TEST_ONLY=TRUE -- Ada's native "evaluate
     the currently-loaded parameters exactly once, no search, then
     save and exit" mode (Counter=0, Spread=0, exits before any Markov
     call -- found in gem-lte-primitives-solution.adb line ~1820).
     Takes a few seconds, not a real optimization run.
Gated on DLOD_FLOOR (0.994, this project's established floor) and
BACKBONE_MIN/MAX (0.1/2.5, same as sweep.py's own gates).

This deliberately does NOT run any random-descent search -- it only
ever compares/copies already-computed results, so it's fast: the goal
per the user's direction is quick, cheap progress and a better
collective starting point for whatever optimization comes next, not a
fresh fit.
"""
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, "/home/paul/eval/gem-lte-core/experiments/Sep2026")
import sweep

import numpy as np

DLOD_RE = re.compile(r"([\d.]+):dLOD:")
IMPROVEMENT_MARGIN = 0.02  # Feb2026 must beat Sep2026 by at least this much


def cc_from_results(csv_path: Path) -> float | None:
    if not csv_path.is_file():
        return None
    try:
        d = np.loadtxt(csv_path, delimiter=",")
    except Exception:
        return None
    model, data = d[:, 1], d[:, 2]
    if np.std(model) == 0 or np.std(data) == 0:
        return None
    return float(np.corrcoef(model, data)[0, 1])


def quick_dlod(cell_dir: Path, cell: str, timeout_s: int = 30) -> float | None:
    """TEST_ONLY=TRUE: evaluate the currently-loaded lt.exe.p exactly
    once (no search), parse the live dLOD reading from stdout."""
    overrides = dict(sweep.BASE_OVERRIDES)
    overrides["CLIMATE_INDEX"] = f"{cell}.dat"
    overrides["TEST_ONLY"] = "true"
    overrides["TIMEOUT"] = str(timeout_s)
    sweep.write_lte_run_sh(cell_dir, overrides)
    sweep.clear_stale_checkpoint(cell_dir, cell)
    try:
        proc = subprocess.run(
            ["bash", "-c", "ulimit -s 65536; ./lte_run.sh"],
            cwd=str(cell_dir), capture_output=True, text=True,
            timeout=timeout_s + 20,
        )
    except subprocess.TimeoutExpired:
        return None
    out = proc.stdout + proc.stderr
    matches = DLOD_RE.findall(out)
    if not matches:
        return None
    return float(matches[-1])


def copy_cell(src_dir: Path, dst_dir: Path, cell: str) -> None:
    for fname in ("lt.exe.p", "lt.exe.resp", f"{cell}.dat"):
        src = src_dir / fname
        if src.is_file():
            shutil.copy2(src, dst_dir / fname)


def backup(cell_dir: Path) -> dict:
    saved = {}
    for fname in ("lt.exe.p", "lt.exe.resp", "lte_run.sh", "lte_results.csv",
                  "lt.exe.windings.json"):
        p = cell_dir / fname
        if p.exists():
            saved[fname] = p.read_bytes()
    return saved


def restore(cell_dir: Path, saved: dict) -> None:
    for fname, content in saved.items():
        (cell_dir / fname).write_bytes(content)


def min_score(entry) -> float:
    scores = (entry or {}).get("scores") or {}
    if "train" in scores:
        return min(scores["train"], scores["validate"], scores["test"])
    if "pair" in scores:
        return min(scores["pair"])
    return -999.0


def main() -> None:
    names = sorted(p.name for p in sweep.HERE.iterdir()
                    if p.is_dir() and sweep.parse_grid_name(p.name))
    by_cell = sweep.latest_entry_per_cell(sweep.load_sweep_ledger())
    results = {}

    for i, cell in enumerate(names):
        sep_dir = sweep.HERE / cell
        feb_dir = sweep.FEB / cell
        sep_cc = cc_from_results(sep_dir / "lte_results.csv")
        feb_cc = cc_from_results(feb_dir / "lte_results.csv")
        print(f"[{i+1}/{len(names)}] {cell}: sep_cc={sep_cc} feb_cc={feb_cc}",
              flush=True)

        if feb_cc is None:
            results[cell] = dict(action="skipped_no_feb_data",
                                  sep_cc=sep_cc, feb_cc=feb_cc)
            continue
        # sep_cc is None whenever Sep2026 has no result of its own yet
        # (e.g. one of the 29 cells that lacked even a .dat file until
        # today) -- treat that as "always worse", not "skip", so Feb2026
        # fills the gap by default (still subject to the same manifold/
        # dLOD gates below).
        effective_sep_cc = sep_cc if sep_cc is not None else -999.0
        sep_disp = f"{sep_cc:.3f}" if sep_cc is not None else "none"
        if feb_cc <= effective_sep_cc + IMPROVEMENT_MARGIN:
            results[cell] = dict(action="kept_sep", sep_cc=sep_cc, feb_cc=feb_cc)
            continue

        bb = sweep.read_backbone(feb_dir)
        if bb is None or not (sweep.BACKBONE_MIN <= abs(bb) <= sweep.BACKBONE_MAX):
            print(f"  Feb2026 wins on CC ({feb_cc:.3f} vs {sep_disp}) but "
                  f"backbone={bb} fails manifold range -- skipped", flush=True)
            results[cell] = dict(action="rejected_manifold", sep_cc=sep_cc,
                                  feb_cc=feb_cc, backbone=bb)
            continue

        saved = backup(sep_dir)
        copy_cell(feb_dir, sep_dir, cell)
        dlod = quick_dlod(sep_dir, cell)
        if dlod is None or dlod < sweep.DLOD_FLOOR:
            print(f"  Feb2026 wins on CC ({feb_cc:.3f} vs {sep_disp}), "
                  f"backbone OK ({bb:.3f}), but dLOD={dlod} fails floor "
                  f"{sweep.DLOD_FLOOR} -- rolled back", flush=True)
            restore(sep_dir, saved)
            results[cell] = dict(action="rejected_dlod", sep_cc=sep_cc,
                                  feb_cc=feb_cc, backbone=bb, dlod=dlod)
            continue

        entry = dict(cell=cell, accepted=True, attempt_tag="feb2026_import",
                     validate_mode=False, metric="CC",
                     scores=dict(metric="CC", pair=[feb_cc, feb_cc]),
                     dlod=dlod, n_attempts=0,
                     timestamp=time.strftime("%Y-%m-%dT%H:%M:%S"),
                     backbone=bb,
                     note=(f"Imported Feb2026's result over Sep2026's "
                           f"(CC {sep_disp} -> {feb_cc:.3f}), no search "
                           f"-- passed manifold ({bb:.3f}) and dLOD "
                           f"({dlod:.4f}) gates."))
        sweep.append_ledger(entry)
        print(f"  IMPORTED from Feb2026: CC {sep_disp} -> {feb_cc:.3f}, "
              f"backbone={bb:.3f}, dLOD={dlod:.4f}", flush=True)
        results[cell] = dict(action="imported", sep_cc=sep_cc, feb_cc=feb_cc,
                              backbone=bb, dlod=dlod)

    print("\n=== SUMMARY ===", flush=True)
    counts = {}
    for r in results.values():
        counts[r["action"]] = counts.get(r["action"], 0) + 1
    for action, n in sorted(counts.items()):
        print(f"  {action}: {n}", flush=True)
    for cell, r in results.items():
        print(f"  {cell}: {r}", flush=True)
    (sweep.HERE / "PICK_BEST_DONE").write_text(
        f"{time.strftime('%Y-%m-%dT%H:%M:%S')}\n" +
        "\n".join(f"{c}: {r}" for c, r in results.items()))


if __name__ == "__main__":
    main()
