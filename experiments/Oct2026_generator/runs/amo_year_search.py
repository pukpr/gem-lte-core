#!/usr/bin/env python3
"""AMO, TIDES=ZONAL v3 (comb year = 365.2422484 + YEAR, anchored 1990.5): 1200 s searches at YEAR
values giving Mt comb aliases of ~120, 127 (production), 135, 143 yr. TEST_ONLY re-score + stats."""
import sys, re, json, shutil, numpy as np
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
E=Path("/home/paul/eval/gem-lte-core/experiments"); H=E/"Oct2026_generator"/"runs"; sys.path.insert(0,str(E/"Sep2026")); import sweep_cc_refine as r, sweep
sys.path.insert(0,str(E/"Oct2026_surrogate")); sys.argv=[sys.argv[0]]
from levels_test import eta2, dbic
src=E/"Feb2026/_rr2/start/amo"; T=1200
alias=lambda y: abs(1/((365.2422484+y)/9.1329504-40))
def run(year):
    d=H/f"ys_amo_{year:+.5f}"
    if d.exists(): shutil.rmtree(d)
    d.mkdir()
    for f in ("lt.exe.p","lt.exe.resp","dlod_ref.dat","amo.dat"): shutil.copy(src/f,d/f)
    sweep.set_resp_key(d/"lt.exe.resp","YEAR",repr(year))
    env=dict(re.findall(r'^export (\w+)=(\S+)',(src/"lte_run.sh").read_text(),re.M)); env.update(CLIMATE_INDEX="amo.dat",TIDES="ZONAL")
    r.run_lt(d,{**env,"TEST_ONLY":"false","TIMEOUT":str(T),"NUMBER_OF_PROCESSORS":"4"},False,T)
    q=r.parse(r.run_lt(d,{**env,"TEST_ONLY":"true","TIMEOUT":"20","NUMBER_OF_PROCESSORS":"2"},True,20)); t3=q["triplet"]
    a=np.loadtxt(d/"lte_results.csv",delimiter=","); b=dbic(a[:,2],a[:,3],a[:,0]); pp=json.load(open(d/"lt.exe.p"))
    return dict(year=year,alias=alias(year),train=t3["train"],validate=t3["validate"],test=t3["test"],dlod=q["dlod"],P=float(eta2(a[:,2],a[:,3])),dBIC=float(b[0]),k=float(b[1]),ltep=[round(x,4) for x in pp["ltep"]],year_p=pp["year"])
with ThreadPoolExecutor(4) as ex: R=list(ex.map(run,(0.0,0.00405305,0.008,0.012)))
json.dump(R,open(H/"amo_year_search.json","w"),indent=1)
for x in R: print(f"YEAR {x['year']:+.5f} (Mt alias {x['alias']:.1f} yr): {x['train']:.4f}/{x['validate']:.4f}/{x['test']:+.4f} dLOD {x['dlod']:.4f} | level {x['P']:.3f} | dBIC {x['dBIC']:+.0f} at k={x['k']} | ltep {x['ltep']} | p.year {x['year_p']:.2e}")
