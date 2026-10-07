#!/usr/bin/env python3
"""First ZONAL-mode searches (2026-10-06): AMO, PDO, NINO4 from their current fits (Feb2026/_rr2/start
for amo/pdo, Oct2026_surrogate/seeds for nino4). TIDES=ZONAL vs TIDES=TABLE (control), 300 s each,
same settings; TEST_ONLY re-score; level-structure and penalized-winding statistics."""
import sys, re, json, shutil, numpy as np
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
E=Path("/home/paul/eval/gem-lte-core/experiments"); H=E/"Oct2026_generator"/"runs"; sys.path.insert(0,str(E/"Sep2026")); import sweep_cc_refine as r
sys.path.insert(0,str(E/"Oct2026_surrogate")); sys.argv=[sys.argv[0]]
from levels_test import eta2, dbic
SRC={"amo":E/"Feb2026/_rr2/start/amo","pdo":E/"Feb2026/_rr2/start/pdo","nino4":E/"Oct2026_surrogate/seeds/nino4"}
def run(i, mode, search):
    src=SRC[i]; d=H/f"{'s' if search else 'f'}_{i}_{mode}"
    if d.exists(): shutil.rmtree(d)
    d.mkdir()
    for f in ("lt.exe.p","lt.exe.resp","dlod_ref.dat",f"{i}.dat"): shutil.copy(src/f,d/f)
    env=dict(re.findall(r'^export (\w+)=(\S+)',(src/"lte_run.sh").read_text(),re.M)); env.update(CLIMATE_INDEX=f"{i}.dat",TIDES=mode)
    if search: r.run_lt(d,{**env,"TEST_ONLY":"false","TIMEOUT":"300","NUMBER_OF_PROCESSORS":"5"},False,300)
    p=r.parse(r.run_lt(d,{**env,"TEST_ONLY":"true","TIMEOUT":"20","NUMBER_OF_PROCESSORS":"2"},True,20)); t3=p["triplet"]
    a=np.loadtxt(d/"lte_results.csv",delimiter=","); b=dbic(a[:,2],a[:,3],a[:,0]); pp=json.load(open(d/"lt.exe.p"))
    return dict(i=i,mode=mode,search=search,metric=t3["metric"],train=t3["train"],validate=t3["validate"],test=t3["test"],dlod=p["dlod"],P=float(eta2(a[:,2],a[:,3])),dBIC=float(b[0]),k=float(b[1]),ltep=[round(x,4) for x in pp["ltep"]])
jobs=[(i,m,False) for i in SRC for m in ("TABLE","ZONAL")]+[(i,m,True) for i in SRC for m in ("ZONAL","TABLE")]
with ThreadPoolExecutor(3) as ex: R=list(ex.map(lambda j: run(*j), jobs))
json.dump(R,open(H/"zonal_search.json","w"),indent=1)
for x in R: print(f"{x['i']:5s} {x['mode']:5s} {'300s search' if x['search'] else 'no refit   '}: {x['metric']} {x['train']:.4f}/{x['validate']:.4f}/{x['test']:+.4f} dLOD {x['dlod']:.4f} | level {x['P']:.3f} | winding dBIC {x['dBIC']:+.0f} at k={x['k']} | ltep {x['ltep']}")
