#!/usr/bin/env python3
"""Full 89-quad basin-seeded fast-refine Kaplan sweep (Oct2026_quad5_besk).

Each of the 89 Kaplan quads in all89_basins.json (30 Atlantic, 41 Pacific, 18 Indian) is the
centre of a 5-leg run with its neighbours.json neighbours, seeded from our refined climate-index
fits:
  - Atlantic <- atlantic_family + AMO (besk 0.2014)
  - Pacific  <- indopacific_family + NINO4 (besk 0.2306)
  - Indian   <- indopacific_family + IODE (besk 0.2306)

For each quad:
  1. Set up 5-leg directory (template.resp, LEG_WINDINGS/LEG_P off, 2000-2005 held out,
     ACCEL=TRUE, MANIFOLD_PRE=TRUE).
  2. Run native Ada REFINE=TRUE (seed-anchored lpap prior, lam=1.5, 3 rounds x 8 iters,
     shared Block B + conservative harmonic screening max_swaps=2 min_gain=0.002).
  3. Run TEST_ONLY=TRUE (F9=1) to write lte_results.1..5.csv, run.out, legs.png, result.p.
  4. Immediately rescore the centre leg unfiltered (F9=-1, ACCEL=TRUE, MANIFOLD_PRE=TRUE)
     into _global_sst/nofilt/<q>/ for the global monthly SST reconstruction.
  5. Write DONE_89 marker (resumable).
"""
import json
import os
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
E = HERE.parent
DLOD = E / "Oct2026_layer" if (E / "Oct2026_layer").exists() else E / "Oct2026_dlod"
MAKE = E / "Oct2026_generator/legs/make_leg_run.py"
NB_PATH = HERE / "neighbours.json" if (HERE / "neighbours.json").exists() else E / "Oct2026_quad5/neighbours.json"
BASIN_PATH = HERE / "all89_basins.json" if (HERE / "all89_basins.json").exists() else HERE / "sample.json"
NOFILT = HERE / "_global_sst/nofilt"

REFINE_LAM, REFINE_ROUNDS, REFINE_ITERS = 1.5, 3, 8
THREADS, JOBS = 4, 4

basins = json.loads(BASIN_PATH.read_text())
nb = json.loads(NB_PATH.read_text())


def log(msg):
    with open(HERE / "sweep89.log", "a") as f:
        f.write(time.strftime("%Y-%m-%d %H:%M:%S ") + msg + "\n")


def rescore_nofilt(q, d):
    out = NOFILT / q
    out.mkdir(parents=True, exist_ok=True)
    shutil.copy(d / "result.p", out / "lt.exe.p")
    shutil.copy(d / f"{q}.dat", out / f"{q}.dat")
    lines = [
        l for l in (d / "lt.exe.resp").read_text().splitlines()
        if not re.match(r"^(LEGS|-?LEG_)", l)
    ]
    lines.append(f"CLIMATE_INDEX                  {q}.dat")
    (out / "lt.exe.resp").write_text("\n".join(lines) + "\n")
    lt_link = out / "lt.exe"
    if not lt_link.exists():
        lt_link.symlink_to(E / "Feb2026/enso_opt")
    env = dict(
        os.environ,
        DLOD_REF="TRUE",
        EXCLUDE="true",
        TRAIN_START="2000",
        TRAIN_END="2005",
        F9="-1",
        ACCEL="TRUE",
        MANIFOLD_PRE="TRUE",
        LOCKA="FALSE",
        LOCKT="FALSE",
        METRIC="CC",
        TREND="true",
        TEST_ONLY="TRUE",
        TIMEOUT="60",
        NUMBER_OF_PROCESSORS="2",
    )
    with open(out / "run.out", "w") as f:
        r = subprocess.run(
            ["bash", "-c", "ulimit -s 65536; ./lt.exe"],
            cwd=out,
            stdout=f,
            stderr=subprocess.STDOUT,
            stdin=subprocess.DEVNULL,
            env=env,
        )
    return r.returncode == 0 and (out / "lte_results.csv").exists()


def run(q):
    d = HERE / q
    if (d / "DONE_89").exists() and (NOFILT / q / "lte_results.csv").exists():
        log(f"{q}: done already, skipped")
        return
    t0 = time.time()
    args = [
        sys.executable, str(MAKE), q,
        "--root", str(HERE),
        "--seed", str(HERE / f"seed_{basins[q]}.p"),
        "--template", str(HERE / "template.resp"),
        "--extra-env", "ACCEL=TRUE MANIFOLD_PRE=TRUE",
        "--leg", f"{DLOD / q / (q + '.dat')}:C_{q[1:]}",
    ]
    for n in nb[q]:
        args += ["--leg", f"{DLOD / n['quad'] / (n['quad'] + '.dat')}:{n['side']}_{n['quad'][1:]}"]
    if d.exists():
        args.append("--force")
    out = subprocess.run(args, capture_output=True, text=True)
    if out.returncode:
        log(f"{q} ({basins[q]}): setup FAILED\n{out.stdout}\n{out.stderr}")
        return

    ref_env = dict(
        os.environ,
        REFINE="TRUE",
        REFINE_HOMOTOPY="FALSE",
        REFINE_ANCHOR_SEED="TRUE",
        REFINE_HARM_SCREEN="TRUE",
        REFINE_HARM_MAX_SWAPS="2",
        REFINE_HARM_MIN_GAIN="0.002",
        REFINE_PER_LEG="TRUE",
        REFINE_LAM=str(REFINE_LAM),
        REFINE_ROUNDS=str(REFINE_ROUNDS),
        REFINE_ITERS=str(REFINE_ITERS),
        REFINE_APPLY="TRUE",
        DLOD_DAT="../dlod3.dat",
        ACCEL="TRUE",
        MANIFOLD_PRE="TRUE",
        DLOD_REF="TRUE",
        EXCLUDE="true",
        TRAIN_START="2000",
        TRAIN_END="2005",
        F9="1",
        LOCKA="FALSE",
        LOCKT="FALSE",
        METRIC="CC",
        TREND="true",
        NUMBER_OF_PROCESSORS=str(THREADS),
    )
    with open(d / "refine.out", "w") as f:
        subprocess.run(
            ["bash", "-c", "ulimit -s 65536; ../lt.exe"],
            cwd=d, stdout=f, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, env=ref_env,
        )
    with open(d / "run.out", "w") as f:
        subprocess.run(
            ["bash", "lte_run.sh", "60", str(THREADS)],
            cwd=d, stdout=f, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
            env=dict(os.environ, TEST_ONLY="TRUE"),
        )
    if (d / "lt.exe.p").exists():
        (d / "result.p").write_bytes((d / "lt.exe.p").read_bytes())

    nf_ok = rescore_nofilt(q, d)
    if (d / "legs.png").exists() and nf_ok:
        (d / "DONE_89").write_text(time.strftime("%Y-%m-%d %H:%M:%S\n"))
    summary = [l for l in (d / "refine.out").read_text().splitlines() if l.startswith("REFINE summary:")]
    log(
        f"{q} ({basins[q]}): {time.time() - t0:.1f} s; "
        f"{summary[-1] if summary else 'no summary'}; "
        f"legs.png={'ok' if (d / 'legs.png').exists() else 'MISS'} nofilt={'ok' if nf_ok else 'MISS'}"
    )


if __name__ == "__main__":
    NOFILT.mkdir(parents=True, exist_ok=True)
    for f, src in (
        ("lt.exe", E / "Feb2026/enso_opt"),
        ("leg_scores2.py", E / "Oct2026_generator/legs/leg_scores2.py"),
        ("dlod3.dat", E / "Feb2026/dlod3.dat"),
    ):
        if not (HERE / f).exists():
            (HERE / f).symlink_to(src)
    dlod_nf = NOFILT / "dlod3.dat"
    if not dlod_nf.exists():
        dlod_nf.symlink_to(E / "Feb2026/dlod3.dat")

    qs = sys.argv[1:] or sorted(basins)
    log(
        f"89-quad fast-refine start: {len(qs)} quads, lam={REFINE_LAM} rounds={REFINE_ROUNDS} "
        f"iters={REFINE_ITERS} harm_swaps=2(>=0.002), {THREADS} threads x {JOBS} jobs"
    )
    with ThreadPoolExecutor(JOBS) as ex:
        list(ex.map(run, qs))
    if len(qs) == len(basins):
        (HERE / "SWEEP89_DONE").write_text(time.strftime("%Y-%m-%d %H:%M:%S\n"))
    log("89-quad fast-refine done")
