import sys, re, json, shutil, numpy as np
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
E=Path("/home/paul/eval/gem-lte-core/experiments"); H=E/"Oct2026_generator/savebug"; sys.path.insert(0,str(E/"Sep2026")); import sweep_cc_refine as r, sweep
src=H/"w4"; ENV=dict(re.findall(r'^export (\w+)=(\S+)',(E/"Oct2026_generator/runs/f_amo_ZONAL/lte_run.sh").read_text(),re.M)); ENV["METRIC"]="CC"   # pinned: do not follow the live folder
P=json.loads((src/"lt.exe.p").read_text())
def run(name, mod, dela0=False, search=0):
    d=H/name
    if d.exists(): shutil.rmtree(d)
    d.mkdir()
    for f in ("lt.exe.resp","dlod_ref.dat"): shutil.copy(src/f,d/f)
    shutil.copy(src/"w4.dat",d/f"{name}.dat"); p=json.loads(json.dumps(P))
    if dela0: p["delA"]=0.0
    (d/"lt.exe.p").write_text(json.dumps(p,indent=2)); sweep.set_resp_key(d/"lt.exe.resp","IMPULSE_MOD12","TRUE" if mod else "FALSE")
    env={**ENV,"CLIMATE_INDEX":f"{name}.dat"}
    if search: r.run_lt(d,{**env,"TEST_ONLY":"false","TIMEOUT":str(search),"NUMBER_OF_PROCESSORS":"8"},False,search)
    q=r.parse(r.run_lt(d,{**env,"TEST_ONLY":"true","TIMEOUT":"20","NUMBER_OF_PROCESSORS":"2"},True,20)); pp=json.loads((d/"lt.exe.p").read_text())
    return name,q["triplet"],{k:round(pp[k],4) for k in ("delA","delB","asym")}
mode=sys.argv[1]
jobs=[("imp_off",False),("imp_on",True),("imp_on_delA0",True,True)] if mode=="fixed" else [("imp_off_s",False,False,300),("imp_on_s",True,False,300)]
with ThreadPoolExecutor(2) as ex: R=list(ex.map(lambda j: run(*j), jobs))
for n,t,k in R: print(f"{n:14s}: train/validate/test {t['train']:.4f}/{t['validate']:.4f}/{t['test']:+.4f} | {k}")
