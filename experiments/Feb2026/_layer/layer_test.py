#!/usr/bin/env python3
"""LAYER experiment (2026-10-02): two-stage manifold, lowest winding as the first layer.
Scratch copies of Feb2026/<index>; for each: baseline score (TEST_ONLY, LAYER=0), then
searches of TIMEOUT s for LAYER=0 (control), 1 (replace) and 2 (additive), each re-scored
with TEST_ONLY. Usage: layer_test.py <timeout> <index> [<index> ...]"""
import json, re, shutil, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
B = Path(__file__).resolve().parent; FEB = B.parent; E = FEB.parent
sys.path.insert(0, str(E / "Sep2026")); import sweep, sweep_cc_refine as r
TIMEOUT = int(sys.argv[1]); IDX = sys.argv[2:]
for f, tgt in (("lt.exe", FEB / "enso_opt"), ("dlod3.dat", FEB / "dlod3.dat")):
    if not (B / f).exists(): (B / f).symlink_to(tgt)
def env_of(src):
    ov = dict(sweep.BASE_OVERRIDES)
    if (src / "lte_run.sh").exists(): ov.update(dict(re.findall(r'^export (\w+)=(\S+)', (src / "lte_run.sh").read_text(), re.M)))
    return ov
def score(d, ov):
    out = r.parse(r.run_lt(d, {**ov, "TEST_ONLY": "true", "TIMEOUT": "30", "NUMBER_OF_PROCESSORS": "2"}, True, 30))
    a = np.loadtxt(d / "lte_results.csv", delimiter=","); t, m, x = a[:, 0], a[:, 1], a[:, 2]
    ts, te = float(ov.get("TRAIN_START", 2000)), float(ov.get("TRAIN_END", 2005)); gap = (t >= ts) & (t < te)
    j = json.loads((d / "lt.exe.windings.json").read_text()); p = json.loads((d / "lt.exe.p").read_text())
    return dict(triplet=out["triplet"], dlod=out["dlod"], cc_all=float(np.corrcoef(m, x)[0, 1]), cc_train=float(np.corrcoef(m[~gap], x[~gap])[0, 1]), cc_test=float(np.corrcoef(m[gap], x[gap])[0, 1]),
                F_range=[float(a[:, 3].min()), float(a[:, 3].max())], k=[round(v[0], 4) for v in j["k_amp_phase"]], amp=[round(v[1], 3) for v in j["k_amp_phase"]], impA=p["impA"], impB=p["impB"], ltep=p["ltep"])
def run(job):
    idx, layer, search = job; src = FEB / idx; d = B / f"{idx}__L{layer}{'s' if search else 'b'}"
    if d.exists(): shutil.rmtree(d)
    d.mkdir()
    for f in ("lt.exe.p", "lt.exe.resp", "dlod_ref.dat"):
        if (src / f).exists(): shutil.copy(src / f, d / f)
    shutil.copy(src / f"{idx}.dat", d / f"{d.name}.dat")
    ov = env_of(src); ov.update(CLIMATE_INDEX=f"{d.name}.dat", LAYER=str(layer))
    rec = dict(index=idx, layer=layer, search=search)
    try:
        if search: r.run_lt(d, {**ov, "TEST_ONLY": "false", "TIMEOUT": str(TIMEOUT), "NUMBER_OF_PROCESSORS": "5"}, False, TIMEOUT)
        rec.update(score(d, ov))
    except Exception as exc: rec["error"] = repr(exc)[:300]
    print(json.dumps(rec), flush=True); return rec
jobs = [(i, 0, False) for i in IDX] + [(i, L, True) for i in IDX for L in (0, 1, 2)]
with ThreadPoolExecutor(3) as ex: recs = list(ex.map(run, jobs))
json.dump(recs, open(B / f"results_{'_'.join(IDX)}_{TIMEOUT}.json", "w"), indent=1); print("DONE")
