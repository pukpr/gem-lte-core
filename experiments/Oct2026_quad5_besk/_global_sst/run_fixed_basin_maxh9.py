#!/usr/bin/env python3
"""Run all 89 Kaplan quads on the FIXED basin lpap manifolds (no per-quad lpap drift),
with ACCEL=FALSE, MAXH=9, NH=5 (all harmonics <= 9), F9=-1, extrapolated backward to 1880
(EXCLUDE=FALSE, TRAIN_START=1950.0, TRAIN_END=2023.0 -> pure 1880-1949 OOB validation)."""
import json
import os
import re
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SWEEP = HERE.parent
EXP = SWEEP.parent
OUTROOT = HERE / "fixed_basin_maxh9"
OUTROOT.mkdir(parents=True, exist_ok=True)

dlod_link = OUTROOT / "dlod3.dat"
if not dlod_link.exists():
    dlod_link.symlink_to(EXP / "Feb2026/dlod3.dat")

basins = json.loads((SWEEP / "all89_basins.json").read_text())
pre_t = np.array([y + m / 12.0 for y in range(1880, 1950) for m in range(12)])
pre_rows = np.column_stack([pre_t, np.zeros_like(pre_t)])

# Build MAXH=9 fixed basin seed .p files (preserving each basin's exact lpap, besk, ltep, IR,
# and keeping only harmonics <= 9 in the first 5 slots)
HARM_MAXH9 = {
    "atlantic": [4.0, 5.0, 6.0, 7.0, 9.0],
    "pacific":  [2.0, 3.0, 5.0, 6.0, 7.0],
    "indian":   [3.0, 5.0, 6.0, 7.0, 9.0],
}
seed_paths = {}
for b, hs in HARM_MAXH9.items():
    p = json.loads((SWEEP / f"seed_{b}.p").read_text())
    p["harm"] = hs + [0.0] * (len(p["harm"]) - len(hs))
    sp = OUTROOT / f"seed_{b}_maxh9.p"
    sp.write_text(json.dumps(p, indent=2) + "\n")
    seed_paths[b] = sp

# Build template.resp with ACCEL=FALSE, MAXH=9, NH="1 1 1 1 1", no LEGS
base_resp = []
for l in (SWEEP / "template.resp").read_text().splitlines():
    if re.match(r"^(LEGS|-?LEG_|MAXH|NH\s)", l):
        continue
    base_resp.append(l)
base_resp += [
    "MAXH                           9",
    'NH                             "1 1 1 1 1"',
    "ACCEL                          FALSE",
]


def run_one(q):
    src = SWEEP / q
    out = OUTROOT / q
    out.mkdir(parents=True, exist_ok=True)
    b = basins[q]
    (out / "lt.exe.p").write_bytes(seed_paths[b].read_bytes())
    orig = np.loadtxt(src / f"{q}.dat")
    np.savetxt(out / f"{q}.dat", np.vstack([pre_rows, orig]), fmt="%.6f\t%.6f")
    resp_lines = list(base_resp) + [f"CLIMATE_INDEX                  {q}.dat"]
    (out / "lt.exe.resp").write_text("\n".join(resp_lines) + "\n")
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
        ACCEL="FALSE",
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
print(f"Fixed-basin MAXH=9 ACCEL=FALSE: {ok_n}/{len(res)} quads done in {time.time() - t0:.1f} s")
for q, rc, ok in res:
    if rc != 0 or not ok:
        print(f"  FAILED: {q} rc={rc} ok={ok}")
