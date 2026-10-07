#!/usr/bin/env python3
"""Step B: constituents fixed at the LOD-only regression (DLOD_REF=FALSE loads them; LOCKT=TRUE
freezes them); everything else is searched for 300 s from each index's current fit. Control: the
same search with the climate-fit constituents frozen (LOCKT=TRUE, DLOD_REF default), for equal
freedom."""
import sys, re, json, shutil, numpy as np
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
E=Path("/home/paul/eval/gem-lte-core/experiments"); H=E/"Oct2026_lodonly"; sys.path.insert(0,str(E/"Sep2026")); import sweep_cc_refine as r
sys.path.insert(0,str(E/"Oct2026_surrogate")); sys.argv=[sys.argv[0]]
from levels_test import eta2, dbic
def run(i, mode):
    src=E/"Oct2026_surrogate/seeds"/i; d=H/f"stepB_{i}_{mode}"
    if d.exists(): shutil.rmtree(d)
    d.mkdir()
    for f in ("lt.exe.p","lt.exe.resp","dlod_ref.dat",f"{i}.dat"): shutil.copy(src/f,d/f)
    env=dict(re.findall(r'^export (\w+)=(\S+)',(src/"lte_run.sh").read_text(),re.M)); env.update(CLIMATE_INDEX=f"{i}.dat",LOCKT="TRUE")
    if mode=="lod": env.update(DLOD_REF="FALSE")
    r.run_lt(d,{**env,"TEST_ONLY":"false","TIMEOUT":"300","NUMBER_OF_PROCESSORS":"5"},False,300)
    p=r.parse(r.run_lt(d,{**env,"TEST_ONLY":"true","TIMEOUT":"20","NUMBER_OF_PROCESSORS":"2"},True,20)); t3=p["triplet"]
    a=np.loadtxt(d/"lte_results.csv",delimiter=","); b=dbic(a[:,2],a[:,3],a[:,0])
    pp=json.load(open(d/"lt.exe.p"))
    return dict(i=i,mode=mode,train=t3["train"],validate=t3["validate"],test=t3["test"],dlod=p["dlod"],P=float(eta2(a[:,2],a[:,3])),dBIC=float(b[0]),k=float(b[1]),ltep=[round(x,4) for x in pp["ltep"]])
with ThreadPoolExecutor(3) as ex: R=list(ex.map(lambda j: run(*j), [(i,m) for m in ("lod","fit") for i in ("amo","pdo","nino4")]))
json.dump(R,open(H/"stepB.json","w"),indent=1)
for x in R: print(f"{x['i']:5s} {'LOD-only constituents ' if x['mode']=='lod' else 'climate-fit constituents'} (frozen, 300 s search): train/val/test {x['train']:.3f}/{x['validate']:.3f}/{x['test']:+.3f}  dLOD {x['dlod']:.4f}  level-structure {x['P']:.3f}  penalized winding dBIC {x['dBIC']:+.0f} at k={x['k']}  ltep {x['ltep']}")
