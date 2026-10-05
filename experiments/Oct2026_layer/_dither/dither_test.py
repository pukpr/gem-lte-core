#!/usr/bin/env python3
"""ANNUAL_DITHER test (build 2026-10-03 12:31). Fixed-parameter grid (regression re-solved) on
amo, pdo, kN000_W090/W110/W130: stage 1 (after layer) a in A1 x 4 phases, stage 0 (tidal manifold)
a in A0 x 2 phases, controls a = 0 (exact reproduction) and a = 1e-6 (climatology removal alone).
Optional: 'search' arg runs 300 s searches from each target's best grid point and from a=1e-6."""
import json, math, re, shutil, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
B = Path(__file__).resolve().parent; L = B.parent; E = L.parent; FEB = E / "Feb2026"
sys.path.insert(0, str(E / "Sep2026")); import sweep, sweep_cc_refine as r
for f, tgt in (("lt.exe", FEB / "enso_opt"), ("dlod3.dat", FEB / "dlod3.dat")):
    if not (B / f).exists(): (B / f).symlink_to(tgt)
TARGETS = {"amo": FEB / "amo", "pdo": FEB / "pdo", **{q: L / q for q in ("kN000_W090", "kN000_W110", "kN000_W130")}}
def make(name, tgt, a, ph, st):
    src = TARGETS[tgt]; d = B / name
    if d.exists(): shutil.rmtree(d)
    d.mkdir()
    for f in ("lt.exe.resp", "lt.exe.p", "dlod_ref.dat"): shutil.copy(src / f, d / f)
    shutil.copy(src / f"{tgt}.dat", d / f"{name}.dat")
    for k, v in (("ANNUAL_DITHER", repr(a)), ("ANNUAL_DITHER_PHASE", repr(ph)), ("ANNUAL_DITHER_STAGE", str(st))): sweep.set_resp_key(d / "lt.exe.resp", k, v)
    return d
def score(tgt, d, search=0):
    env = {**dict(re.findall(r'^export (\w+)=(\S+)', (TARGETS[tgt] / "lte_run.sh").read_text(), re.M)), "CLIMATE_INDEX": f"{d.name}.dat"}
    if search: r.run_lt(d, {**env, "TEST_ONLY": "false", "TIMEOUT": str(search), "NUMBER_OF_PROCESSORS": "4"}, False, search)
    out = r.parse(r.run_lt(d, {**env, "TEST_ONLY": "true", "TIMEOUT": "30", "NUMBER_OF_PROCESSORS": "2"}, True, 30))
    a = np.loadtxt(d / "lte_results.csv", delimiter=","); t, m, x = a[:, 0], a[:, 1], a[:, 2]
    k = (x != 0) & (t >= (1950 if tgt.startswith("k") else 1880)); tt = t[k] - t[k].mean()
    D = x[k] - np.polyval(np.polyfit(tt, x[k], 1), tt); M = m[k] - np.polyval(np.polyfit(tt, m[k], 1), tt)
    hi, lo = D >= np.quantile(D, .9), D <= np.quantile(D, .1)
    ctx = json.loads((d / "lt.exe.windings.json").read_text())["secular_context"]
    return dict(name=d.name, tgt=tgt, a=ctx.get("annual_dither"), ph=ctx.get("annual_dither_phase"), st=ctx.get("annual_dither_stage"), dlod=out["dlod"],
                **{kk: out["triplet"][kk] for kk in ("train", "validate", "test")}, warm=float(M[hi].mean() / D[hi].mean()), cold=float(M[lo].mean() / D[lo].mean()))
if __name__ == "__main__":
    A1 = (0.01, 0.02, 0.04, 0.08); A0 = (2.0, 8.0); PH = [i * math.pi / 2 for i in range(4)]
    jobs = [(t, 0.0, 0.0, 1) for t in TARGETS] + [(t, 1e-6, 0.0, 1) for t in TARGETS] + [(t, a, p, 1) for t in TARGETS for a in A1 for p in PH] + [(t, a, p, 0) for t in TARGETS for a in A0 for p in PH[:2]]
    with ThreadPoolExecutor(4) as ex:
        G = list(ex.map(lambda j: score(j[0], make(f"g_{j[0]}_s{j[3]}_a{j[1]:g}_p{j[2]:.2f}", *j)), jobs))
    for t in TARGETS:
        gs = [x for x in G if x["tgt"] == t]; b0 = next(x for x in gs if x["a"] == 0.0); b1 = next(x for x in gs if 0 < x["a"] < 1e-5)
        f = lambda x: f"{x['train']:.3f}/{x['validate']:.3f}/{x['test']:+.3f} w{x['warm']:.2f}c{x['cold']:.2f}"
        print(f"{t}: a=0 {f(b0)} | deseason only {f(b1)}")
        for st, As in ((1, A1), (0, A0)):
            for a in As:
                print(f"   stage {st} a={a:<5g} " + "  ".join(f"ph{x['ph']:.1f}:{f(x)}" for x in gs if x["st"] == st and abs(x["a"] - a) < 1e-9), flush=True)
    json.dump(dict(grid=G), open(B / "dither_grid.json", "w"), indent=1); print("DONE")
