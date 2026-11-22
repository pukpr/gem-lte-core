#!/usr/bin/env python3
"""Set up (and optionally start) a multi-leg lt.exe run: N records fitted on one shared
manifold and set of windings, each leg regressing its own amplitudes (see LEGS.md).

    python3 make_leg_run.py sam_ks040w050 --leg sam --leg kS040_W050 --seed pdo_sam/result.p --seed-test --run 1800 8

Creates <subdir>/ next to this script with
    seed.p, lt.exe.p          copy of --seed (shared manifold, windings, harmonics)
    lt.exe.resp               --template's settings with LEGS set to the number of legs
    lt.exe.resp.1 .. .N       CLIMATE_INDEX (and SAMPLING when not monthly) per leg
    <leg data>.dat            copied, so the run directory is self-contained
    lte_run.sh                held-out window, metric, run time and threads; draws legs.png after the search
    dlod3.dat                 link to experiments/Feb2026/dlod3.dat
A leg is an index name (experiments/Feb2026/<name>/<name>.dat) or a path to a .dat file;
append :LABEL to set the plot label (default: the name in capitals) and @SAMPLING for
non-monthly data, e.g.  --leg ../QBO/era_qbo_u50_point0E_5day.dat:ERA50@73
--seed-test   scores the seed without searching (TEST_ONLY): per-leg train / held-out correlation
--run S T     starts the search in the background for S seconds on T threads (log: <subdir>/run.out)
An existing subdir is never overwritten unless --force is given. --root DIR builds it under DIR instead of here.
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
FEB = HERE.parents[1] / "Feb2026"
TEMPLATE = HERE / "amo_pdo_nino4_iode" / "lt.exe.resp"   # common-manifold settings (INIT_DATE 1880, MAXH 9, ...)


def parse_leg(spec):
    sampling = None
    if "@" in spec:
        spec, sampling = spec.rsplit("@", 1)
    label = None
    if ":" in spec:
        spec, label = spec.rsplit(":", 1)
    p = Path(spec)
    if p.suffix != ".dat":
        p = FEB / spec / f"{spec}.dat"
    elif not p.is_absolute():
        p = (Path.cwd() / p).resolve()
    if not p.exists():
        sys.exit(f"leg data not found: {p}")
    return p, (label or p.stem.upper()), sampling


def write_resp(template, n_legs, out):
    text = Path(template).read_text()
    text, n = re.subn(r"(?m)^LEGS\s+.*$", f"LEGS                           {n_legs}", text)
    if n == 0:
        text = text.rstrip("\n") + f"\nLEGS                           {n_legs}\n"
    out.write_text(text)


def scores(run_dir, n_legs, ts, te):
    rows = []
    for k in range(1, n_legs + 1):
        a = np.loadtxt(run_dir / f"lte_results.{k}.csv", delimiter=",")
        d, m, y = a[:, 0], a[:, 1], a[:, 2]
        ok = y != 0
        h = (d >= ts) & (d <= te)
        r = lambda s: float(np.corrcoef(m[s], y[s])[0, 1]) if s.sum() > 3 else float("nan")
        rows.append((r(ok & ~h), r(ok & h)))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("subdir")
    ap.add_argument("--leg", action="append", required=True, help="index name or .dat path [:LABEL][@SAMPLING]")
    ap.add_argument("--seed", required=True, help="lt.exe.p to start from (relative to this directory or absolute)")
    ap.add_argument("--template", default=str(TEMPLATE), help="lt.exe.resp to copy settings from")
    ap.add_argument("--held", nargs=2, type=float, default=(2000, 2005), metavar=("START", "END"))
    ap.add_argument("--metric", default="CC")
    ap.add_argument("--extra-env", default="", help="extra exports for lte_run.sh, e.g. 'TIDE_LEGACY=TRUE'")
    ap.add_argument("--seed-test", action="store_true")
    ap.add_argument("--run", nargs=2, type=int, metavar=("SECONDS", "THREADS"))
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--root", default=str(HERE), help="directory to build the run in (needs lt.exe and leg_scores2.py there or linked)")
    a = ap.parse_args()

    run_dir = Path(a.root).resolve() / a.subdir
    if run_dir.exists() and not a.force:
        sys.exit(f"{run_dir} exists; use --force to rebuild it (results there would be overwritten)")
    legs = [parse_leg(s) for s in a.leg]
    seed = Path(a.seed) if Path(a.seed).is_absolute() else (HERE / a.seed)
    if not seed.exists():
        sys.exit(f"seed not found: {seed}")
    run_dir.mkdir(exist_ok=True)

    shutil.copy(seed, run_dir / "seed.p")
    shutil.copy(seed, run_dir / "lt.exe.p")
    write_resp(a.template, len(legs), run_dir / "lt.exe.resp")
    for k, (dat, label, sampling) in enumerate(legs, 1):
        shutil.copy(dat, run_dir / dat.name)
        lines = [f"CLIMATE_INDEX                  {dat.name}"]
        if sampling:
            lines.append(f"SAMPLING                       {float(sampling)}")
        (run_dir / f"lt.exe.resp.{k}").write_text("\n".join(lines) + "\n")
    labels = " ".join(label for _, label, _ in legs)
    ts, te = a.held
    names = ", ".join(f"{label} (leg {k})" for k, (_, label, _) in enumerate(legs, 1))
    sh = run_dir / "lte_run.sh"
    sh.write_text(f"""#!/usr/bin/env bash
# {len(legs)} legs on one manifold: {names}; seed {seed.relative_to(HERE) if seed.is_relative_to(HERE) else seed}.
# {ts:g}-{te:g} is held out on every leg. Made by make_leg_run.py; see LEGS.md.
# usage: ./lte_run.sh [seconds, default 1800] [threads, default 16]       TEST_ONLY=TRUE ./lte_run.sh  scores only
cd "$(dirname "$0")"
export DLOD_REF=TRUE EXCLUDE=true TRAIN_START={ts:g} TRAIN_END={te:g} F9=1 LOCKA=FALSE LOCKT=FALSE METRIC={a.metric} TREND=true {a.extra_env}
export TIMEOUT=${{1:-1800}} NUMBER_OF_PROCESSORS=${{2:-16}} TEST_ONLY=${{TEST_ONLY:-FALSE}}
ulimit -s 65536
../lt.exe
python3 ../leg_scores2.py . {labels} --held {ts:g} {te:g}
""")
    sh.chmod(0o755)
    link = run_dir / "dlod3.dat"
    if not link.exists():
        link.symlink_to(os.path.relpath(FEB / "dlod3.dat", run_dir))
    print(f"set up {run_dir}: legs {names}")

    if a.seed_test:
        env = {**os.environ, "TEST_ONLY": "TRUE"}
        subprocess.run(["bash", "-c", "ulimit -s 65536; ../lt.exe"], cwd=run_dir, env={
            **env, **dict(kv.split("=", 1) for kv in
                          f"DLOD_REF=TRUE EXCLUDE=true TRAIN_START={ts:g} TRAIN_END={te:g} F9=1 LOCKA=FALSE "
                          f"LOCKT=FALSE METRIC={a.metric} TREND=true TIMEOUT=60 NUMBER_OF_PROCESSORS={len(legs)} "
                          f"{a.extra_env}".split())},
            stdout=open(run_dir / "seed_test.out", "w"), stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL)
        for (_, label, _), (tr, ho) in zip(legs, scores(run_dir, len(legs), ts, te)):
            print(f"  seed, {label:12s} train {tr:+.3f}  held out {ts:g}-{te:g} {ho:+.3f}")

    if a.run:
        secs, threads = a.run
        subprocess.Popen(["nohup", "bash", "lte_run.sh", str(secs), str(threads)], cwd=run_dir,
                         stdout=open(run_dir / "run.out", "w"), stderr=subprocess.STDOUT,
                         stdin=subprocess.DEVNULL, start_new_session=True)
        print(f"  search started: {secs} s on {threads} threads; log {run_dir}/run.out, "
              f"then legs.png and result scores (copy lt.exe.p to result.p to keep it)")


if __name__ == "__main__":
    main()
