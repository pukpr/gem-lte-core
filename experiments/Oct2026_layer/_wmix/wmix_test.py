#!/usr/bin/env python3
"""LAYER=1 with impC recast as the mixing weight w of the tidal manifold (F2 = w*Fb + layer),
build 2026-10-03 10:10. Targets: Feb2026 amo, pdo; Oct2026_layer kN000_W090/W110/W130.
 1. grid: TEST_ONLY with impC = w for w in WGRID (w = 0 is the pure layer)
 2. 300 s searches from w = 0 and from the best grid w (by train), re-scored
Reports train/validate/test, dLOD, and warm/cold tail capture (top/bottom decile of the
linearly detrended data, model mean / data mean) over the scored record."""
import json, re, shutil, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
B = Path(__file__).resolve().parent; L = B.parent; E = L.parent; FEB = E / "Feb2026"
sys.path.insert(0, str(E / "Sep2026")); import sweep_cc_refine as r
for f, tgt in (("lt.exe", FEB / "enso_opt"), ("dlod3.dat", FEB / "dlod3.dat")):
    if not (B / f).exists(): (B / f).symlink_to(tgt)
TARGETS = {"amo": FEB / "amo", "pdo": FEB / "pdo", **{q: L / q for q in ("kN000_W090", "kN000_W110", "kN000_W130")}}
WGRID = (-0.01, -0.005, -0.002, 0.0, 0.002, 0.005, 0.01)
def make(name, tgt, w, p=None):
    src = TARGETS[tgt]; d = B / name
    if d.exists(): shutil.rmtree(d)
    d.mkdir()
    for f in ("lt.exe.resp", "dlod_ref.dat"): shutil.copy(src / f, d / f)
    shutil.copy(src / f"{tgt}.dat", d / f"{name}.dat")
    p = p or json.loads((src / "lt.exe.p").read_text()); p["impC"] = w; (d / "lt.exe.p").write_text(json.dumps(p, indent=2)); return d
def env(tgt, d): return {**dict(re.findall(r'^export (\w+)=(\S+)', (TARGETS[tgt] / "lte_run.sh").read_text(), re.M)), "CLIMATE_INDEX": f"{d.name}.dat"}
def score(tgt, d, search=0):
    if search: r.run_lt(d, {**env(tgt, d), "TEST_ONLY": "false", "TIMEOUT": str(search), "NUMBER_OF_PROCESSORS": "4"}, False, search)
    out = r.parse(r.run_lt(d, {**env(tgt, d), "TEST_ONLY": "true", "TIMEOUT": "30", "NUMBER_OF_PROCESSORS": "2"}, True, 30))
    a = np.loadtxt(d / "lte_results.csv", delimiter=","); t, m, x, F = a[:, 0], a[:, 1], a[:, 2], a[:, 3]
    k = (x != 0) & (t >= 1950 if tgt.startswith("k") else t >= 1880); tt = t[k] - t[k].mean()
    D = x[k] - np.polyval(np.polyfit(tt, x[k], 1), tt); M = m[k] - np.polyval(np.polyfit(tt, m[k], 1), tt)
    hi, lo = D >= np.quantile(D, .9), D <= np.quantile(D, .1); p = json.loads((d / "lt.exe.p").read_text())
    return dict(name=d.name, tgt=tgt, w=p["impC"], dlod=out["dlod"], **{kk: out["triplet"][kk] for kk in ("train", "validate", "test")},
                warm=float(M[hi].mean() / D[hi].mean()), cold=float(M[lo].mean() / D[lo].mean()), F2=[float(F.min()), float(F.max())],
                ltep=p["ltep"], p=p)
def show(x): return (f"{x['name']:26s} w {x['w']:+.4f}: {x['train']:.4f}/{x['validate']:.4f}/{x['test']:.4f} dLOD {x['dlod']:.4f} | "
                     f"tail capture warm {x['warm']:.2f} cold {x['cold']:.2f} | F2 {x['F2'][0]:+.2f}..{x['F2'][1]:+.2f}")
if __name__ == "__main__":
    with ThreadPoolExecutor(8) as ex:
        G = list(ex.map(lambda j: score(j[0], make(f"g_{j[0]}_w{j[1]:+g}", j[0], j[1])), [(t, w) for t in TARGETS for w in WGRID]))
    print("GRID (fixed parameters, regression re-solved)")
    for x in G: print("  " + show(x))
    jobs = []
    for t in TARGETS:
        gs = [x for x in G if x["tgt"] == t]; best = max(gs, key=lambda x: x["train"])
        jobs += [(t, 0.0), (t, best["w"])] if best["w"] != 0.0 else [(t, 0.0), (t, 0.005 if best["w"] == 0 else best["w"])]
    with ThreadPoolExecutor(4) as ex:
        S = list(ex.map(lambda j: score(j[0], make(f"s_{j[0]}_w{j[1]:+g}", j[0], j[1]), search=300), jobs))
    print("SEARCHES (300 s), start w -> final")
    for x, j in zip(S, jobs): print(f"  start w {j[1]:+.4f} -> " + show(x))
    json.dump(dict(grid=[{k: v for k, v in x.items() if k != 'p'} for x in G], search=[{k: v for k, v in x.items() if k != 'p'} for x in S]), open(B / "wmix.json", "w"), indent=1)
    print("DONE")
