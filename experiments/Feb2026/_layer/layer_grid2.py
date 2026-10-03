#!/usr/bin/env python3
"""LAYER v2 (2026-10-02 18:16 build): impA/impB stay the first-order Bessel correction of the
tidal manifold; the lowest-winding layer's amplitude pair is offs/bg. For one index:
 1. LAYER=0 TEST_ONLY on the new binary (must reproduce the current fit)
 2. grid: amplitude x phase (0..150 deg; +180 is the same layer up to sign, absorbed by the
    MLR) x k1 factor, LAYER=1, everything else fixed, MLR re-solved (TEST_ONLY)
 3. 300 s searches from the 3 best grid points, plus a LAYER=0 control search, re-scored.
Usage: layer_grid2.py <index>"""
import json, re, shutil, sys, math, itertools
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
B = Path(__file__).resolve().parent; FEB = B.parent; E = FEB.parent
sys.path.insert(0, str(E / "Sep2026")); import sweep, sweep_cc_refine as r
idx = sys.argv[1]; src = FEB / idx; G = B / f"grid2_{idx}"; G.mkdir(exist_ok=True)
for f, tgt in (("lt.exe", FEB / "enso_opt"), ("dlod3.dat", FEB / "dlod3.dat")):
    if not (G / f).exists(): (G / f).symlink_to(tgt)
p0 = json.loads((src / "lt.exe.p").read_text()); k0 = p0["ltep"][0]
ov0 = dict(sweep.BASE_OVERRIDES); ov0.update(dict(re.findall(r'^export (\w+)=(\S+)', (src / "lte_run.sh").read_text(), re.M)))
gs, ge = float(ov0.get("TRAIN_START", 2000)), float(ov0.get("TRAIN_END", 2005))
def prep(name, p):
    d = G / name
    if d.exists(): shutil.rmtree(d)
    d.mkdir()
    for f in ("lt.exe.resp", "dlod_ref.dat"): shutil.copy(src / f, d / f)
    shutil.copy(src / f"{idx}.dat", d / f"{d.name}.dat"); (d / "lt.exe.p").write_text(json.dumps(p, indent=2)); return d
def score(d, layer):
    out = r.parse(r.run_lt(d, {**ov0, "CLIMATE_INDEX": f"{d.name}.dat", "LAYER": str(layer), "TEST_ONLY": "true", "TIMEOUT": "30", "NUMBER_OF_PROCESSORS": "2"}, True, 30))
    a = np.loadtxt(d / "lte_results.csv", delimiter=","); t, m, x = a[:, 0], a[:, 1], a[:, 2]; gap = (t >= gs) & (t < ge)
    pp = json.loads((d / "lt.exe.p").read_text())
    return dict(triplet=out["triplet"], dlod=out["dlod"], cc_out=float(np.corrcoef(m[~gap], x[~gap])[0, 1]), cc_gap=float(np.corrcoef(m[gap], x[gap])[0, 1]),
                F=[float(a[:, 3].min()), float(a[:, 3].max())], offs=pp["offs"], bg=pp["bg"], k1=pp["ltep"][0])
def grid_job(job):
    i, (A, ph, kf) = job; p = json.loads(json.dumps(p0)); p["offs"] = A * math.cos(ph); p["bg"] = A * math.sin(ph); p["ltep"][0] = k0 * kf
    rec = dict(A=A, phase=ph, kf=kf)
    try: rec.update(score(prep(f"g{i:03d}", p), 1))
    except Exception as exc: rec["error"] = repr(exc)[:200]
    return rec
def search_job(job):
    name, p, layer = job; d = prep(name, p); rec = dict(name=name, layer=layer)
    try:
        rec["start"] = score(d, layer)
        r.run_lt(d, {**ov0, "CLIMATE_INDEX": f"{d.name}.dat", "LAYER": str(layer), "TEST_ONLY": "false", "TIMEOUT": "300", "NUMBER_OF_PROCESSORS": "5"}, False, 300)
        rec["end"] = score(d, layer)
    except Exception as exc: rec["error"] = repr(exc)[:200]
    print(json.dumps(rec), flush=True); return rec
res = dict(baseline=score(prep("baseline_L0", p0), 0)); print("baseline", json.dumps(res["baseline"]), flush=True)
grid = list(itertools.product((1, 3, 6, 10, 15, 22, 30), [i * math.pi / 6 for i in range(6)], (0.75, 1.0, 1.5, 2.0)))
with ThreadPoolExecutor(8) as ex: res["grid"] = list(ex.map(grid_job, enumerate(grid)))
ok = sorted([g for g in res["grid"] if "error" not in g and g.get("triplet")], key=lambda g: -g["cc_out"])
print("grid done:", len(ok), "of", len(grid), "| best:", json.dumps(ok[:3]), flush=True)
jobs = [("search_L0_control", p0, 0)]
for n, g in enumerate(ok[:3]):
    p = json.loads(json.dumps(p0)); p["offs"] = g["A"] * math.cos(g["phase"]); p["bg"] = g["A"] * math.sin(g["phase"]); p["ltep"][0] = k0 * g["kf"]; jobs.append((f"search_L1_best{n+1}", p, 1))
with ThreadPoolExecutor(3) as ex: res["search"] = list(ex.map(search_job, jobs))
json.dump(res, open(B / f"grid2_{idx}.json", "w"), indent=1); print("DONE")
