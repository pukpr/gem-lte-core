#!/usr/bin/env python3
"""LAYER=1 seed scan for AMO: the first-stage winding is seeded with the phase the MLR
itself found for the lowest winding (0.515 rad, i.e. impA = A cos(phi), impB = A sin(phi))
and a range of amplitudes A (both signs), since the 600 s search from the old impA/impB
never moved them. 300 s searches, each re-scored with TEST_ONLY."""
import json, re, shutil, sys, math
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
B = Path(__file__).resolve().parent; FEB = B.parent; E = FEB.parent
sys.path.insert(0, str(E / "Sep2026")); import sweep, sweep_cc_refine as r
idx = "amo"; src = FEB / idx; TIMEOUT = 300
j0 = json.loads((src / "lt.exe.windings.json").read_text())["k_amp_phase"][0]; PHI = j0[2]
ov0 = dict(sweep.BASE_OVERRIDES); ov0.update(dict(re.findall(r'^export (\w+)=(\S+)', (src / "lte_run.sh").read_text(), re.M)))
def run(A):
    d = B / f"{idx}__seedA{A:+.0f}"
    if d.exists(): shutil.rmtree(d)
    d.mkdir()
    for f in ("lt.exe.resp", "dlod_ref.dat"): shutil.copy(src / f, d / f)
    shutil.copy(src / f"{idx}.dat", d / f"{d.name}.dat")
    p = json.loads((src / "lt.exe.p").read_text()); p["impA"] = A * math.cos(PHI); p["impB"] = A * math.sin(PHI)
    (d / "lt.exe.p").write_text(json.dumps(p, indent=2))
    ov = {**ov0, "CLIMATE_INDEX": f"{d.name}.dat", "LAYER": "1"}; rec = dict(A=A)
    try:
        for stage in ("start", "end"):
            if stage == "end": r.run_lt(d, {**ov, "TEST_ONLY": "false", "TIMEOUT": str(TIMEOUT), "NUMBER_OF_PROCESSORS": "5"}, False, TIMEOUT)
            out = r.parse(r.run_lt(d, {**ov, "TEST_ONLY": "true", "TIMEOUT": "30", "NUMBER_OF_PROCESSORS": "2"}, True, 30))
            a = np.loadtxt(d / "lte_results.csv", delimiter=","); pp = json.loads((d / "lt.exe.p").read_text())
            rec[stage] = dict(triplet=out["triplet"], dlod=out["dlod"], cc_all=float(np.corrcoef(a[:, 1], a[:, 2])[0, 1]), F=[float(a[:, 3].min()), float(a[:, 3].max())], impA=pp["impA"], impB=pp["impB"], ltep=pp["ltep"])
    except Exception as exc: rec["error"] = repr(exc)[:300]
    print(json.dumps(rec), flush=True); return rec
with ThreadPoolExecutor(3) as ex: recs = list(ex.map(run, [5, -5, 15, -15, 30, -30, 60, -60, 120]))
json.dump(recs, open(B / "seed_scan_amo.json", "w"), indent=1); print("DONE")
