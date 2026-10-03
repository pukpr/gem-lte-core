#!/usr/bin/env python3
"""AMO with LAYER=1 reduced to a SINGLE regressed winding (NM 6 -> 2: ltep = [k1 layer, k],
no harmonics), from the user's 19:39 Feb2026/amo .p/.resp. Scans k with everything else
fixed (TEST_ONLY), user's lte_run.sh env. Usage: single_k.py k1 k2 ..."""
import json, re, shutil, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
B = Path(__file__).resolve().parent; FEB = B.parent.parent; E = FEB.parent
sys.path.insert(0, str(E / "Sep2026")); import sweep, sweep_cc_refine as r
src = FEB / "amo"
for f, tgt in (("lt.exe", FEB / "enso_opt"), ("dlod3.dat", FEB / "dlod3.dat")):
    if not (B / f).exists(): (B / f).symlink_to(tgt)
ov0 = dict(re.findall(r'^export (\w+)=(\S+)', (src / "lte_run.sh").read_text(), re.M))
def make(name, k, p=None):
    d = B / name
    if d.exists(): shutil.rmtree(d)
    d.mkdir()
    for f in ("lt.exe.resp", "dlod_ref.dat"): shutil.copy(src / f, d / f)
    shutil.copy(src / "amo.dat", d / f"{name}.dat")
    p = p or json.loads((src / "lt.exe.p").read_text()); p["ltep"] = [p["ltep"][0], k]
    (d / "lt.exe.p").write_text(json.dumps(p, indent=2)); sweep.set_resp_key(d / "lt.exe.resp", "NM", "2")
    return d
def score(d, search=0):
    ov = {**ov0, "CLIMATE_INDEX": f"{d.name}.dat"}
    if search: r.run_lt(d, {**ov, "TEST_ONLY": "false", "TIMEOUT": str(search), "NUMBER_OF_PROCESSORS": "5"}, False, search)
    out = r.parse(r.run_lt(d, {**ov, "TEST_ONLY": "true", "TIMEOUT": "30", "NUMBER_OF_PROCESSORS": "2"}, True, 30))
    j = json.loads((d / "lt.exe.windings.json").read_text()); p = json.loads((d / "lt.exe.p").read_text())
    return dict(name=d.name, triplet=out["triplet"], dlod=out["dlod"], kap=j["k_amp_phase"], ltep=p["ltep"], offs=p["offs"], bg=p["bg"])
if __name__ == "__main__":
    ks = [float(x) for x in sys.argv[1:]]
    with ThreadPoolExecutor(8) as ex: recs = list(ex.map(lambda k: score(make(f"k{k:g}", k)), ks))
    for x in recs:
        t = x["triplet"]; print(f"k {x['ltep'][1]:5.3f}: train/val/test {t['train']:.4f}/{t.get('validate',0):.4f}/{t['test']:.4f} dLOD {x['dlod']:.4f} amp {x['kap'][1][1]:.3f} phase {x['kap'][1][2]:+.2f}")
