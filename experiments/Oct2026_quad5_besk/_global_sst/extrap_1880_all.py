#!/usr/bin/env python3
"""Extrapolate all 89 fitted Oct2026_quad5_besk quads backward to 1880 (unfiltered, F9=-1)
using native lt.exe with EXCLUDE=FALSE TRAIN_START=1950.0 TRAIN_END=2023.0.

Pre-1950 rows (1880.0000..1949.9167) are padded with 0.0 and strictly outside TRAIN_START..TRAIN_END,
so 1880-1949 (840 months) is a pure non-fit out-of-bag (OOB) backward validation."""
import json
import os
import re
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SWEEP = HERE.parent
EXP = SWEEP.parent
OUTROOT = HERE / "extrap1880"
OUTROOT.mkdir(parents=True, exist_ok=True)

dlod_link = OUTROOT / "dlod3.dat"
if not dlod_link.exists():
    dlod_link.symlink_to(EXP / "Feb2026/dlod3.dat")

basins = json.loads((SWEEP / "all89_basins.json").read_text())
pre_t = np.array([y + m / 12.0 for y in range(1880, 1950) for m in range(12)])
pre_rows = np.column_stack([pre_t, np.zeros_like(pre_t)])


def run_one(q):
    src = SWEEP / q
    out = OUTROOT / q
    out.mkdir(parents=True, exist_ok=True)
    shutil.copy(src / "result.p", out / "lt.exe.p")
    orig = np.loadtxt(src / f"{q}.dat")
    np.savetxt(out / f"{q}.dat", np.vstack([pre_rows, orig]), fmt="%.6f\t%.6f")
    lines = [
        l for l in (src / "lt.exe.resp").read_text().splitlines()
        if not re.match(r"^(LEGS|-?LEG_)", l)
    ]
    lines.append(f"CLIMATE_INDEX                  {q}.dat")
    (out / "lt.exe.resp").write_text("\n".join(lines) + "\n")
    lt_link = out / "lt.exe"
    if not lt_link.exists():
        lt_link.symlink_to(EXP / "Feb2026/enso_opt")
    env = dict(
        os.environ,
        DLOD_REF="TRUE",
        EXCLUDE="FALSE",
        TRAIN_START="1950.0",
        TRAIN_END="2023.0",
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
    ok = (out / "lte_results.csv").exists()
    return q, r.returncode, ok


t0 = time.time()
with ThreadPoolExecutor(8) as ex:
    res = list(ex.map(run_one, sorted(basins)))
ok_n = sum(1 for _, rc, ok in res if rc == 0 and ok)
print(f"Extrapolated {ok_n}/{len(res)} quads to 1880-2023 (F9=-1) in {time.time() - t0:.1f} s")
for q, rc, ok in res:
    if rc != 0 or not ok:
        print(f"  FAILED: {q} rc={rc} ok={ok}")
