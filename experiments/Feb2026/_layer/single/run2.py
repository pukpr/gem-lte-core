import json, shutil, sys
from pathlib import Path
sys.argv = [sys.argv[0]]; from single_k import *
# reference: the user's NM=6 layered fit, scored the same way
d = B / "ref_nm6"
if d.exists(): shutil.rmtree(d)
d.mkdir()
for f in ("lt.exe.resp", "lt.exe.p", "dlod_ref.dat"): shutil.copy(src / f, d / f)
shutil.copy(src / "amo.dat", d / "ref_nm6.dat")
out = {"ref_nm6": score(d)}; print(json.dumps(out["ref_nm6"]), flush=True)
for k in (2.5, 2.6):
    out[f"search_k{k}"] = score(make(f"search_k{k:g}", k), search=300); print(json.dumps(out[f"search_k{k}"]), flush=True)
json.dump(out, open(B / "run2.json", "w"), indent=1); print("DONE")
