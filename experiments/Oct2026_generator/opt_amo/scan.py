#!/usr/bin/env python3
"""Scans on a frozen snapshot of the user's AMO ZONAL fit (start/): YEAR, impulse month x asym, and
(IMPULSE_MOD12 on) the second-impulse ratio delA/asym. TEST_ONLY, regression re-solved."""
import sys, re, json, shutil, numpy as np
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
E=Path("/home/paul/eval/gem-lte-core/experiments"); H=E/"Oct2026_generator/opt_amo"; S=H/"start"
sys.path.insert(0,str(E/"Sep2026")); import sweep_cc_refine as r, sweep
ENV=dict(re.findall(r'^export (\w+)=(\S+)',(S/"lte_run.sh").read_text(),re.M)); P=json.loads((S/"lt.exe.p").read_text())
def run(name, p=None, resp=None, search=0):
    d=H/"w"/name
    if d.exists(): shutil.rmtree(d)
    d.mkdir(parents=True)
    for f in ("lt.exe.resp","dlod_ref.dat"): shutil.copy(S/f,d/f)
    shutil.copy(S/"f_amo_ZONAL.dat",d/f"{name}.dat"); (d/"lt.exe.p").write_text(json.dumps(p or P,indent=2))
    for k,v in (resp or {}).items(): sweep.set_resp_key(d/"lt.exe.resp",k,v)
    env={**ENV,"CLIMATE_INDEX":f"{name}.dat"}
    if search: r.run_lt(d,{**env,"TEST_ONLY":"false","TIMEOUT":str(search),"NUMBER_OF_PROCESSORS":"5"},False,search)
    out=r.run_lt(d,{**env,"TEST_ONLY":"true","TIMEOUT":"20","NUMBER_OF_PROCESSORS":"2"},True,20)
    m=re.findall(r"([-\d.]+):dLOD:",out); a=np.loadtxt(d/"lte_results.csv",delimiter=","); t=a[:,0]; g=(t>=2000)&(t<2005)
    cc=lambda k: float(np.corrcoef(a[k,1],a[k,2])[0,1])
    return dict(name=name,train=cc(~g),validate=float("nan"),test=cc(g),dlod=float(m[-1]) if m else float("nan"),p=json.loads((d/"lt.exe.p").read_text()))
def par(jobs, n=8):
    with ThreadPoolExecutor(n) as ex: return list(ex.map(lambda j: run(*j), jobs))
def show(x, extra=""): print(f"  {x['name']:22s} CC outside 2000-05 {x['train']:.4f} | 2000-05 {x['test']:+.4f} | dLOD {x['dlod']:.4f} {extra}")
alias=lambda y: abs(1/((365.2422484+y)/9.1329504-40))
if __name__=="__main__":
    base=run("base"); print("BASELINE (your fit, TEST_ONLY):"); show(base)
    Y=[round(0.00405305+d,6) for d in np.arange(-0.008,0.0161,0.001)]
    R=par([(f"year_{y:+.6f}",None,{"YEAR":repr(y)}) for y in Y])
    print("YEAR scan (Mt comb alias):")
    for y,x in zip(Y,R): show(x, f"YEAR {y:+.6f}  Mt alias {alias(y):6.1f} yr" + ("  <- current" if abs(y-0.00405305)<1e-9 else ""))
    by=max(zip(Y,R),key=lambda z:z[1]["train"])[0]
    jobs=[]
    for m in range(12):
        for sc in (0.5,0.75,1.0,1.25,1.5):
            p=json.loads(json.dumps(P)); p["delB"]=((m-6)%12+18)/12.0; p["asym"]=P["asym"]*sc; jobs.append((f"imp_m{m:02d}_s{sc}",p,{"YEAR":repr(by)}))
    R2=par(jobs); best=sorted(R2,key=lambda x:-x["train"])[:8]
    print(f"IMPULSE month x asym scale (at the best YEAR {by:+.6f}); top 8 by train:")
    for x in best: show(x)
    mbest=max(R2,key=lambda x:x["train"]); pb=mbest["p"]
    jobs=[(f"ratio_{q:+.2f}",dict(pb,delA=q*pb["asym"]),{"YEAR":repr(by),"IMPULSE_MOD12":"TRUE"}) for q in (-1.0,-0.5,-0.25,-0.1,0.0,0.1,0.25,0.5,1.0)]
    R3=par(jobs); print(f"SECOND IMPULSE (IMPULSE_MOD12 on) at the best month ({mbest['name']}): delA = ratio x asym")
    for x in R3: show(x)
    json.dump(dict(base={k:v for k,v in base.items() if k!='p'},best_year=by,best_imp=mbest["name"],best_imp_p=pb,
                   ratio=[{k:v for k,v in x.items() if k!='p'} for x in R3]),open(H/"scan.json","w"),indent=1)
