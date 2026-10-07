#!/usr/bin/env python3
"""DECLARATION_YEAR.md: YEAR-detuning ("wrong clock") test. Usage: year_test.py y1 | y2"""
import json, re, shutil, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
H = Path(__file__).resolve().parent; E = H.parent; Y = H / "_year"
sys.path.insert(0, str(E / "Sep2026")); import sweep, sweep_cc_refine as r
sys.argv_mode = sys.argv[1]; sys.argv = [sys.argv[0]]
from levels_test import eta2, dbic
BASE = 0.00405305; IDX = ("amo", "pdo", "nino4")
for f, tgt in (("lt.exe", E / "Oct2026_layer_ir" / "lt.exe"), ("dlod3.dat", E / "Feb2026" / "dlod3.dat")):
    try: (Y / f).symlink_to(tgt)
    except FileExistsError: pass
def run(i, dl, search):
    name = f"{i}_{'s' if search else 'f'}{dl:+.4f}"; d = Y / name
    if d.exists(): shutil.rmtree(d)
    d.mkdir(); src = H / "seeds" / i
    for f in ("lt.exe.p", "lt.exe.resp", "dlod_ref.dat"): shutil.copy(src / f, d / f)
    shutil.copy(src / f"{i}.dat", d / f"{name}.dat"); sweep.set_resp_key(d / "lt.exe.resp", "YEAR", repr(BASE + dl))
    env = {**dict(re.findall(r'^export (\w+)=(\S+)', (src / "lte_run.sh").read_text(), re.M)), "CLIMATE_INDEX": f"{name}.dat"}
    if search: r.run_lt(d, {**env, "TEST_ONLY": "false", "TIMEOUT": "300", "NUMBER_OF_PROCESSORS": "4"}, False, 300)
    out = r.parse(r.run_lt(d, {**env, "TEST_ONLY": "true", "TIMEOUT": "30", "NUMBER_OF_PROCESSORS": "2"}, True, 30))
    a = np.loadtxt(d / "lte_results.csv", delimiter=","); t, x, F2 = a[:, 0], a[:, 2], a[:, 3]
    t3 = out["triplet"]
    return dict(index=i, delta=dl, search=search, P=eta2(x, F2), dBIC=dbic(x, F2, t)[0], k=dbic(x, F2, t)[1], dlod=out["dlod"], **{k: t3[k] for k in ("train", "validate", "test")})
if __name__ == "__main__":
    mode = sys.argv_mode
    D = (0.0, -0.0005, 0.0005, -0.001, 0.001, -0.002, 0.002, -0.005, 0.005, -0.01, 0.01, -0.02, 0.02, -0.05, 0.05) if mode == "y1" else (0.0, -0.01, 0.01, -0.05, 0.05)
    with ThreadPoolExecutor(8 if mode == "y1" else 4) as ex: R = list(ex.map(lambda j: run(*j), [(i, d, mode == "y2") for i in IDX for d in D]))
    json.dump(R, open(H / f"year_{mode}.json", "w"), indent=1)
    for i in IDX:
        rs = sorted([x for x in R if x["index"] == i], key=lambda x: x["delta"]); z = next(x for x in rs if x["delta"] == 0)
        print(f"== {i} ({mode}): delta | P (level eta2_adj) | train/validate/test | dLOD | dBIC (k)")
        for x in rs: print(f"   {x['delta']:+.4f} | {x['P']:.3f}{' *' if x['P'] >= max(y['P'] for y in rs) else '  '} | {x['train']:.3f}/{x['validate']:.3f}/{x['test']:+.3f} | {x['dlod']:.4f}{'!' if x['dlod'] < 0.99 else ' '} | {x['dBIC']:+.0f} ({x['k']})")
