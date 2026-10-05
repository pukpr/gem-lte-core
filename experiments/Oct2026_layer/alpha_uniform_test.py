#!/usr/bin/env python3
"""Manifold-INDEPENDENT positive gain on the LAYER=1 quads (2026-10-03): ALPHA=100 makes the
ALPHA window cover the whole (layered) manifold, so the summed winding term is multiplied by
ALPHA_GAIN whenever it is positive. TEST_ONLY on copies of _backcast (1880-2023 data, MLR fitted
on 1950-2023), so one run gives the fit period and the 1880-1949 back-cast.
Usage: alpha_uniform_test.py <gain> [<gain> ...]   (results in _alpha_g<gain>/)"""
import json, re, shutil, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
L = Path(__file__).resolve().parent; E = L.parent
sys.path.insert(0, str(E / "Sep2026")); import sweep, sweep_cc_refine as r
def run(job):
    gain, q = job; src = L / "_backcast" / q; d = L / f"_alpha_g{gain}" / q
    if d.exists(): shutil.rmtree(d)
    d.mkdir(parents=True)
    for f in ("lt.exe.p", "lt.exe.resp", "dlod_ref.dat", f"{q}.dat"): shutil.copy(src / f, d / f)
    ov = dict(sweep.BASE_OVERRIDES); ov.update(CLIMATE_INDEX=f"{q}.dat", METRIC="CC", F9="1", ACCEL="FALSE", EXCLUDE="false", TRAIN_START="1950", TRAIN_END="2023",
                                               TEST_ONLY="true", TIMEOUT="30", NUMBER_OF_PROCESSORS="2", ALPHA="100", ALPHA_GAIN=str(gain))
    try:
        r.run_lt(d, ov, True, 30); c = json.loads((d / "lt.exe.windings.json").read_text())["secular_context"]
        return dict(q=q, gain=gain, alpha=c.get("alpha"), alpha_gain=c.get("alpha_gain"))
    except Exception as exc: return dict(q=q, gain=gain, error=repr(exc)[:150])
for g in sys.argv[1:]:
    root = L / f"_alpha_g{g}"; root.mkdir(exist_ok=True)
    for f, tgt in (("lt.exe", L / "lt.exe"), ("dlod3.dat", E / "Feb2026" / "dlod3.dat")):
        if not (root / f).exists(): (root / f).symlink_to(tgt)
quads = sorted(p.name for p in (L / "_backcast").iterdir() if p.is_dir() and sweep.parse_grid_name(p.name))
with ThreadPoolExecutor(8) as ex: recs = list(ex.map(run, [(g, q) for g in sys.argv[1:] for q in quads]))
bad = [x for x in recs if "error" in x or x.get("alpha") != 100.0]
print(len(recs), "runs;", len(bad), "problems", bad[:3])
