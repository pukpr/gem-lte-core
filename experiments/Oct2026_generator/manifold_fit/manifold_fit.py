#!/usr/bin/env python3
"""Intermediate fit (2026-10-06): fit the TIDES=ZONAL manifold (column 4) to the manifold of the
user's good TABLE-mode AMO fit (runs/f_amo_ZONAL, 13:51, CC 0.810), using lt.exe's FORCING=TRUE +
MLR=TRUE mode (the model output IS the layered manifold), so the inner stages (comb, integrator,
Bessel ~0.2 winding, layer) are steered by the manifold shape, not by the outer AMO windings.
 1. harness check: TABLE mode, TEST_ONLY, same parameters -> CC must be ~1
 2. ZONAL TEST_ONLY from the same parameters (starting CC)
 3. two 1200 s ZONAL searches (production YEAR)
 4. carry each fitted parameter set back to the AMO data (ZONAL, outer windings re-solved, TEST_ONLY)"""
import sys, re, json, shutil, numpy as np
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
E=Path("/home/paul/eval/gem-lte-core/experiments"); H=E/"Oct2026_generator"/"manifold_fit"; R0=E/"Oct2026_generator"/"manifold_fit"/"v2_src"   # frozen copy: TABLE-mode fit from the user's 13:53 backup
sys.path.insert(0,str(E/"Sep2026")); import sweep_cc_refine as r, sweep as r_sweep
sys.path.insert(0,str(E/"Oct2026_surrogate")); sys.argv=[sys.argv[0]]
from levels_test import eta2, dbic
for f,tgt in (("lt.exe",E/"Feb2026/enso_opt"),("dlod3.dat",E/"Feb2026/dlod3.dat")):
    try: (H/f).symlink_to(tgt)
    except FileExistsError: pass

ENV=dict(re.findall(r'^export (\w+)=(\S+)',(R0/"lte_run.sh").read_text(),re.M))
def neutral(p):
    """outer terms neutralized: no IR, no annual impulse, no harmonics, one dummy winding at k=200"""
    q=json.loads(json.dumps(p)); q["IR"]=0.0; q["impC"]=0.0; q["ltep"]=q["ltep"][:2]+[200.0]; q["harm"]=[]; return q
def prep(name, p, datafile, neutralize=False):
    if neutralize: p=neutral(p)
    d=H/name
    if d.exists(): shutil.rmtree(d)
    d.mkdir()
    for f in ("lt.exe.resp","dlod_ref.dat"): shutil.copy(R0/f,d/f)
    if neutralize:
        rp=d/"lt.exe.resp"; t=re.sub(r'^NH\s.*\n','',rp.read_text(),flags=re.M); rp.write_text(t+'NH                             ""\n')
        r_sweep.set_resp_key(rp,"VALIDATE","FALSE")
    shutil.copy(datafile,d/f"{name}.dat"); (d/"lt.exe.p").write_text(json.dumps(p,indent=2)); return d
def go(d, tides, manifold, search=0):
    r_sweep.set_resp_key(d/"lt.exe.resp","TIDES",tides)
    env={**ENV,"CLIMATE_INDEX":f"{d.name}.dat","TIDES":tides,"METRIC":"CC","EXCLUDE":"false","VALIDATE":"FALSE"}
    if manifold: env.update(TRAIN_START="1880.0",TRAIN_END="2023.0")   # score the whole record: the manifold is flat within 2000-05
    if search: r.run_lt(d,{**env,"TEST_ONLY":"false","TIMEOUT":str(search),"NUMBER_OF_PROCESSORS":"5"},False,search)
    q=r.parse(r.run_lt(d,{**env,"TEST_ONLY":"true","TIMEOUT":"20","NUMBER_OF_PROCESSORS":"2"},True,20))
    return q, json.loads((d/"lt.exe.p").read_text())
def cc_col(d, c1, c2):
    b=np.loadtxt(d/"lte_results.csv",delimiter=","); return float(np.corrcoef(b[:,c1],b[:,c2])[0,1])
P0=json.loads((R0/"lt.exe.p").read_text()); TGT=H/"target_manifold_table.dat"
d=prep("table_target",P0,R0/"f_amo_ZONAL.dat"); q,_=go(d,"TABLE",False); b=np.loadtxt(d/"lte_results.csv",delimiter=",")
np.savetxt(TGT,b[:,[0,3]],fmt="%.6f",delimiter="\t")
print(f"0. TABLE target regenerated from the backup params: AMO CC {np.corrcoef(b[:,1],b[:,2])[0,1]:.4f} (user's 13:51 run: 0.8102); tides in JSON: {json.loads((d/'lt.exe.windings.json').read_text())['secular_context'].get('tides')}")
# step 1 dropped: in TABLE mode the target is the run's own manifold, the regression fits it exactly
# and is flagged singular, so no candidate is accepted (degenerate by construction)
d=prep("start_zonal",P0,TGT,True); go(d,"ZONAL",True); print(f"2. ZONAL from the TABLE params (outer neutralized): CC(model, target) = {cc_col(d,1,2):.5f}; CC(ZONAL manifold, target) = {cc_col(d,3,2):.5f}")
KEYS=("delA","delB","asym","ma","mp","init","impA","impB","offs","bg","year")
def search(k):
    d=prep(f"zfit_{k}",P0,TGT,True); q,p=go(d,"ZONAL",True,1200)
    return dict(k=k,cc=cc_col(d,1,2),mcc=cc_col(d,3,2),dlod=q["dlod"],p=p)
with ThreadPoolExecutor(2) as ex: S=list(ex.map(search,(1,2)))
print("   start inner:", {k:round(P0[k],4) for k in KEYS}, "ltep", [round(x,4) for x in P0["ltep"]])
for s in S: print(f"3. ZONAL manifold fit #{s['k']}: CC(model, target) = {s['cc']:.5f}, CC(manifold, target) = {s['mcc']:.5f}, dLOD {s['dlod']:.4f} | ltep {[round(x,4) for x in s['p']['ltep']]} | " + str({k:round(s['p'][k],4) for k in KEYS}))
amo=R0/"f_amo_ZONAL.dat"
for s in S:
    # inner (shared) parameters from the manifold fit, outer terms (windings, harmonics, IR, impC) from AMO's own fit
    pc=json.loads(json.dumps(P0)); fp=s["p"]
    for key in ("lpap","delA","delB","asym","ma","mp","shfT","init","year","impA","impB","offs","bg"): pc[key]=fp[key]
    pc["ltep"]=fp["ltep"][:2]+P0["ltep"][2:]
    d=prep(f"amo_from_zfit_{s['k']}",pc,amo); q,_=go(d,"ZONAL",False); t3=q["triplet"]
    b=np.loadtxt(d/"lte_results.csv",delimiter=","); bb=dbic(b[:,2],b[:,3],b[:,0])
    print(f"4. AMO data with manifold-fit #{s['k']} params (ZONAL, outer re-solved): CC train/test {t3['train']:.4f}/{t3['test']:+.4f} | level {eta2(b[:,2],b[:,3]):.3f} | dBIC {bb[0]:+.0f} at k={bb[1]} | manifold CC vs target {np.corrcoef(b[:,3],np.loadtxt(TGT)[:,1])[0,1]:.4f}")
json.dump([{k:v for k,v in s.items()} for s in S],open(H/"manifold_fit.json","w"),indent=1); print("DONE")
