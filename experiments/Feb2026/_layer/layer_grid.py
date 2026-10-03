#!/usr/bin/env python3
"""Brute-force profile of the first-layer (lowest winding) parameters for AMO, LAYER=1:
amplitude x phase x k1 grid, every other parameter fixed, MLR re-solved (TEST_ONLY).
The random search perturbs one of ~150 parameters at a time by small relative steps and
left impA/impB/k1 essentially unmoved, so the layer is profiled directly here."""
import json, re, shutil, sys, math, itertools
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
B = Path(__file__).resolve().parent; FEB = B.parent; E = FEB.parent
sys.path.insert(0, str(E / "Sep2026")); import sweep, sweep_cc_refine as r
idx = sys.argv[1] if len(sys.argv) > 1 else "amo"; src = FEB / idx; G = B / f"grid_{idx}"; G.mkdir(exist_ok=True)
for f, tgt in (("lt.exe", FEB / "enso_opt"), ("dlod3.dat", FEB / "dlod3.dat")):
    if not (G / f).exists(): (G / f).symlink_to(tgt)
p0 = json.loads((src / "lt.exe.p").read_text()); k0 = p0["ltep"][0]
ov0 = dict(sweep.BASE_OVERRIDES); ov0.update(dict(re.findall(r'^export (\w+)=(\S+)', (src / "lte_run.sh").read_text(), re.M)))
def run(job):
    i, (A, ph, kf) = job; d = G / f"g{i:03d}"
    if d.exists(): shutil.rmtree(d)
    d.mkdir()
    for f in ("lt.exe.resp", "dlod_ref.dat"): shutil.copy(src / f, d / f)
    shutil.copy(src / f"{idx}.dat", d / f"{d.name}.dat")
    p = json.loads(json.dumps(p0)); p["impA"] = A * math.cos(ph); p["impB"] = A * math.sin(ph); p["ltep"][0] = k0 * kf
    (d / "lt.exe.p").write_text(json.dumps(p, indent=2)); rec = dict(A=A, phase=ph, kf=kf, dir=d.name)
    try:
        out = r.parse(r.run_lt(d, {**ov0, "CLIMATE_INDEX": f"{d.name}.dat", "LAYER": "1", "TEST_ONLY": "true", "TIMEOUT": "30", "NUMBER_OF_PROCESSORS": "2"}, True, 30))
        rec.update(triplet=out["triplet"], dlod=out["dlod"])
    except Exception as exc: rec["error"] = repr(exc)[:200]
    return rec
grid = list(itertools.product((3, 6, 10, 15, 22, 30, 45), [i * math.pi / 6 for i in range(12)], (0.5, 1.0, 1.5, 2.0)))
with ThreadPoolExecutor(8) as ex: recs = list(ex.map(run, enumerate(grid)))
json.dump(recs, open(B / f"grid_{idx}.json", "w"), indent=1); print("DONE", len(recs), sum("error" in x for x in recs))
