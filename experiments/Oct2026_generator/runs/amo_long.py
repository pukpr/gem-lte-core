#!/usr/bin/env python3
"""AMO only, TIDES=ZONAL v2 (2026-10-06): 3 independent 1200 s searches from Feb2026/_rr2/start/amo,
plus a 1200 s TIDES=TABLE control. TEST_ONLY re-score, level-structure and penalized-winding stats."""
import sys, re, json, shutil, numpy as np
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
E=Path("/home/paul/eval/gem-lte-core/experiments"); H=E/"Oct2026_generator"/"runs"; sys.path.insert(0,str(E/"Sep2026")); import sweep_cc_refine as r
sys.path.insert(0,str(E/"Oct2026_surrogate")); sys.argv=[sys.argv[0]]
from levels_test import eta2, dbic
src=E/"Feb2026/_rr2/start/amo"; T=1200
def run(name, tides):
    d=H/name
    if d.exists(): shutil.rmtree(d)
    d.mkdir()
    for f in ("lt.exe.p","lt.exe.resp","dlod_ref.dat","amo.dat"): shutil.copy(src/f,d/f)
    env=dict(re.findall(r'^export (\w+)=(\S+)',(src/"lte_run.sh").read_text(),re.M)); env.update(CLIMATE_INDEX="amo.dat",TIDES=tides)
    r.run_lt(d,{**env,"TEST_ONLY":"false","TIMEOUT":str(T),"NUMBER_OF_PROCESSORS":"4"},False,T)
    p=r.parse(r.run_lt(d,{**env,"TEST_ONLY":"true","TIMEOUT":"20","NUMBER_OF_PROCESSORS":"2"},True,20)); t3=p["triplet"]
    a=np.loadtxt(d/"lte_results.csv",delimiter=","); b=dbic(a[:,2],a[:,3],a[:,0]); pp=json.load(open(d/"lt.exe.p"))
    keys=("delA","delB","asym","ma","mp","init","impA","impB","offs","bg")
    return dict(name=name,tides=tides,train=t3["train"],validate=t3["validate"],test=t3["test"],dlod=p["dlod"],P=float(eta2(a[:,2],a[:,3])),dBIC=float(b[0]),k=float(b[1]),
                ltep=[round(x,4) for x in pp["ltep"]],shared={k:round(pp[k],4) for k in keys})
jobs=[("long_amo_zonal_1","ZONAL"),("long_amo_zonal_2","ZONAL"),("long_amo_zonal_3","ZONAL"),("long_amo_table","TABLE")]
with ThreadPoolExecutor(4) as ex: R=list(ex.map(lambda j: run(*j), jobs))
json.dump(R,open(H/"amo_long.json","w"),indent=1)
st=json.load(open(src/"lt.exe.p")); print("start:", {k:round(st[k],4) for k in ("delA","delB","asym","ma","mp","init","impA","impB","offs","bg")}, "ltep", [round(x,4) for x in st["ltep"]])
for x in R: print(f"{x['name']:18s} {x['tides']}: {x['train']:.4f}/{x['validate']:.4f}/{x['test']:+.4f} dLOD {x['dlod']:.4f} | level {x['P']:.3f} | dBIC {x['dBIC']:+.0f} at k={x['k']} | ltep {x['ltep']} | {x['shared']}")
