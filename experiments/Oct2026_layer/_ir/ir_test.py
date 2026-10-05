#!/usr/bin/env python3
"""Symmetric IR (12-month delay differential) test, build 2026-10-03 11:02 (impC = annual impulse;
LAYER_MIX off). IR = 0 in every current fit, and the search's relative steps cannot move a parameter
off 0, so: (1) TEST_ONLY grid of IR, (2) 300 s searches seeded at IR = 0 (control) and IR = 0.1.
Targets: Feb2026 amo, pdo; Oct2026_layer Pacific quads. Reports train/validate/test, dLOD, final IR,
and warm/cold tail capture (top/bottom decile of linearly detrended data; model mean / data mean)."""
import json, re, shutil, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
B = Path(__file__).resolve().parent; L = B.parent; E = L.parent; FEB = E / "Feb2026"
sys.path.insert(0, str(E / "Sep2026")); import sweep_cc_refine as r
for f, tgt in (("lt.exe", FEB / "enso_opt"), ("dlod3.dat", FEB / "dlod3.dat")):
    if not (B / f).exists(): (B / f).symlink_to(tgt)
QUADS = ("kN000_W090", "kN000_W110", "kN000_W130", "kN000_W150", "kN000_W170", "kN000_E170", "kN020_W110", "kS020_W090")
TARGETS = {"amo": FEB / "amo", "pdo": FEB / "pdo", **{q: L / q for q in QUADS}}
def make(name, tgt, ir):
    src = TARGETS[tgt]; d = B / name
    if d.exists(): shutil.rmtree(d)
    d.mkdir()
    for f in ("lt.exe.resp", "dlod_ref.dat"): shutil.copy(src / f, d / f)
    shutil.copy(src / f"{tgt}.dat", d / f"{name}.dat")
    p = json.loads((src / "lt.exe.p").read_text()); p["IR"] = ir; (d / "lt.exe.p").write_text(json.dumps(p, indent=2)); return d
def env(tgt, d): return {**dict(re.findall(r'^export (\w+)=(\S+)', (TARGETS[tgt] / "lte_run.sh").read_text(), re.M)), "CLIMATE_INDEX": f"{d.name}.dat"}
def score(tgt, d, search=0):
    if search: r.run_lt(d, {**env(tgt, d), "TEST_ONLY": "false", "TIMEOUT": str(search), "NUMBER_OF_PROCESSORS": "4"}, False, search)
    out = r.parse(r.run_lt(d, {**env(tgt, d), "TEST_ONLY": "true", "TIMEOUT": "30", "NUMBER_OF_PROCESSORS": "2"}, True, 30))
    a = np.loadtxt(d / "lte_results.csv", delimiter=","); t, m, x = a[:, 0], a[:, 1], a[:, 2]
    k = (x != 0) & (t >= (1950 if tgt.startswith("k") else 1880)); tt = t[k] - t[k].mean()
    D = x[k] - np.polyval(np.polyfit(tt, x[k], 1), tt); M = m[k] - np.polyval(np.polyfit(tt, m[k], 1), tt)
    hi, lo = D >= np.quantile(D, .9), D <= np.quantile(D, .1); p = json.loads((d / "lt.exe.p").read_text())
    return dict(name=d.name, tgt=tgt, IR=p["IR"], dlod=out["dlod"], **{kk: out["triplet"][kk] for kk in ("train", "validate", "test")},
                warm=float(M[hi].mean() / D[hi].mean()), cold=float(M[lo].mean() / D[lo].mean()), amp=float(M.std() / D.std()))
if __name__ == "__main__":
    IRS = (-0.2, -0.1, -0.05, 0.0, 0.05, 0.1, 0.2)
    with ThreadPoolExecutor(8) as ex:
        G = list(ex.map(lambda j: score(j[0], make(f"g_{j[0]}_ir{j[1]:+g}", j[0], j[1])), [(t, v) for t in TARGETS for v in IRS]))
    print("GRID (fixed parameters, regression re-solved): IR: train/validate warm/cold")
    for t in TARGETS: print(f"  {t:11s} " + "  ".join(f"{x['IR']:+.2f}:{x['train']:.3f}/{x['validate']:.3f} w{x['warm']:.2f}c{x['cold']:.2f}" for x in G if x["tgt"] == t), flush=True)
    with ThreadPoolExecutor(4) as ex:
        S = list(ex.map(lambda j: score(j[0], make(f"s_{j[0]}_ir{j[1]:+g}", j[0], j[1]), search=300), [(t, v) for t in TARGETS for v in (0.0, 0.1)]))
    print("SEARCHES 300 s: seed IR=0 (control) | seed IR=0.1")
    for t in TARGETS:
        a, b = [x for x in S if x["tgt"] == t]
        f = lambda x: f"{x['train']:.3f}/{x['validate']:.3f}/{x['test']:+.3f} IR {x['IR']:+.3f} warm {x['warm']:.2f} cold {x['cold']:.2f} amp {x['amp']:.2f}"
        print(f"  {t:11s} {f(a)} | {f(b)}")
    json.dump(dict(grid=G, search=S), open(B / "ir.json", "w"), indent=1); print("DONE")
