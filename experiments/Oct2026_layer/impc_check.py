#!/usr/bin/env python3
"""Does impC (once-a-year delta added to the model output) matter? TEST_ONLY of each fit as is
and with impC = 0, same settings; amo/pdo/nino4 seeds + the 89 LAYER quads."""
import json, re, shutil, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
L = Path(__file__).resolve().parent; E = L.parent; FEB = E / "Feb2026"; B = L / "_impc0"
sys.path.insert(0, str(E / "Sep2026")); import sweep, sweep_cc_refine as r
for f, tgt in (("lt.exe", L / "lt.exe"), ("dlod3.dat", FEB / "dlod3.dat")):
    if not (B / f).exists(): (B / f).symlink_to(tgt)
def run(job):
    name, src, datf, zero = job; d = B / f"{name}__{'z' if zero else 'a'}"
    if d.exists(): shutil.rmtree(d)
    d.mkdir()
    for f in ("lt.exe.resp", "dlod_ref.dat"): shutil.copy(src / f, d / f)
    shutil.copy(src / datf, d / f"{d.name}.dat"); p = json.loads((src / "lt.exe.p").read_text()); imp = p["impC"]
    if zero: p["impC"] = 0.0
    (d / "lt.exe.p").write_text(json.dumps(p, indent=2))
    ov = dict(re.findall(r'^export (\w+)=(\S+)', (src / "lte_run.sh").read_text(), re.M)); ov.update(CLIMATE_INDEX=f"{d.name}.dat", TEST_ONLY="true", TIMEOUT="30", NUMBER_OF_PROCESSORS="2")
    out = r.parse(r.run_lt(d, ov, True, 30)); return dict(name=name, zero=zero, impC=imp, **out["triplet"])
jobs = [(i, FEB / i, f"{i}.dat", z) for i in ("amo", "pdo", "nino4") for z in (False, True)]
quads = sorted(p.name for p in L.iterdir() if p.is_dir() and sweep.parse_grid_name(p.name))
jobs += [(q, L / q, f"{q}.dat", z) for q in quads for z in (False, True)]
with ThreadPoolExecutor(8) as ex: recs = list(ex.map(run, jobs))
R = {(x["name"], x["zero"]): x for x in recs}
print(f"{'fit':11s} {'impC':>8s} | train / validate / test as is  ->  with impC=0")
for i in ("amo", "pdo", "nino4"):
    a, z = R[(i, False)], R[(i, True)]; print(f"{i:11s} {a['impC']:+8.4f} | {a['train']:.4f}/{a['validate']:.4f}/{a['test']:.4f} -> {z['train']:.4f}/{z['validate']:.4f}/{z['test']:.4f}")
d = np.array([[R[(q, True)][k] - R[(q, False)][k] for k in ("train", "validate", "test")] for q in quads]); ic = np.array([R[(q, False)]["impC"] for q in quads])
print(f"89 LAYER quads: |impC| median {np.median(abs(ic)):.4f}, max {abs(ic).max():.4f}; change when zeroed: train median {np.median(d[:,0]):+.4f} (worst {d[:,0].min():+.4f}), validate median {np.median(d[:,1]):+.4f} (worst {d[:,1].min():+.4f}), test median {np.median(d[:,2]):+.4f} (range {d[:,2].min():+.4f}..{d[:,2].max():+.4f})")
print("largest train losses:", ", ".join(f"{q} {d[i,0]:+.4f} (impC {ic[i]:+.3f})" for i, q in sorted(enumerate(quads), key=lambda z: d[z[0], 0])[:5]))
