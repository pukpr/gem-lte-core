#!/usr/bin/env python3
"""Rescore all 30 fitted Oct2026_quad5_besk centre quads unfiltered (F9=-1) with TEST_ONLY=TRUE,
preserving ACCEL=TRUE and MANIFOLD_PRE=TRUE from the sweep."""
import json
import os
import re
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
SWEEP = HERE.parent
EXP = SWEEP.parent
NOFILT = HERE / "nofilt"
NOFILT.mkdir(parents=True, exist_ok=True)

dlod_link = NOFILT / "dlod3.dat"
if not dlod_link.exists():
    dlod_link.symlink_to(EXP / "Feb2026/dlod3.dat")

sample = json.loads((SWEEP / "sample.json").read_text())


def rescore_one(q):
    src = SWEEP / q
    out = NOFILT / q
    out.mkdir(parents=True, exist_ok=True)
    shutil.copy(src / "result.p", out / "lt.exe.p")
    shutil.copy(src / f"{q}.dat", out / f"{q}.dat")
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
    ok = (out / "lte_results.csv").exists()
    return q, r.returncode, ok


t0 = time.time()
with ThreadPoolExecutor(8) as ex:
    res = list(ex.map(rescore_one, sorted(sample)))
ok_n = sum(1 for _, rc, ok in res if rc == 0 and ok)
print(f"Rescored {ok_n}/{len(res)} quads unfiltered (F9=-1) in {time.time() - t0:.1f} s")
for q, rc, ok in res:
    if rc != 0 or not ok:
        print(f"  FAILED: {q} rc={rc} ok={ok}")
