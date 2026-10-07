#!/usr/bin/env python3
"""AMO/PDO common-manifold round-robin v2 -- see DECLARATION.md (noise floor, sharing penalty on
VALIDATE, success / plateau / budget exits). Usage: nohup setsid python3 round_robin2.py > rr2.log 2>&1 &"""
import json, re, shutil, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
B = Path(__file__).resolve().parent; FEB = B.parent; E = FEB.parent; W = B / "work"; START = B / "start"
sys.path.insert(0, str(E / "Sep2026")); import sweep_cc_refine as r
IDX = ("amo", "pdo"); DLOD_MIN = 0.985; ROUNDS = 12; T = 300
SHARED = ("lpap", "delA", "delB", "asym", "ma", "mp", "shfT", "init", "year", "impA", "impB", "offs", "bg")
def log(m): print(f"{time.strftime('%H:%M:%S')} {m}", flush=True)
for d in (W, START): d.mkdir(exist_ok=True)
for f, tgt in (("lt.exe", E / "Oct2026_layer_ir" / "lt.exe"), ("dlod3.dat", FEB / "dlod3.dat")):
    try: (W / f).symlink_to(tgt)
    except FileExistsError: pass
for i in IDX:
    (START / i).mkdir(exist_ok=True)
    for f in ("lt.exe.p", "lt.exe.resp", "lte_run.sh", "dlod_ref.dat", f"{i}.dat"):
        if not (START / i / f).exists(): shutil.copy(FEB / i / f, START / i / f)
ENV = {i: dict(re.findall(r'^export (\w+)=(\S+)', (START / i / "lte_run.sh").read_text(), re.M)) for i in IDX}
def shared_of(p): return {**{k: p[k] for k in SHARED}, "ltep01": p["ltep"][:2]}
def with_shared(p, s):
    q = json.loads(json.dumps(p))
    for k in SHARED: q[k] = s[k]
    q["ltep"] = list(s["ltep01"]) + q["ltep"][2:]; return q
def rundir(name, i, p):
    d = W / name
    if d.exists(): shutil.rmtree(d)
    d.mkdir()
    for f in ("lt.exe.resp", "dlod_ref.dat"): shutil.copy(START / i / f, d / f)
    shutil.copy(START / i / f"{i}.dat", d / f"{name}.dat"); (d / "lt.exe.p").write_text(json.dumps(p, indent=2)); return d
def score(d, i):
    out = r.parse(r.run_lt(d, {**ENV[i], "CLIMATE_INDEX": f"{d.name}.dat", "TEST_ONLY": "true", "TIMEOUT": "30", "NUMBER_OF_PROCESSORS": "2"}, True, 30))
    a = np.loadtxt(d / "lte_results.csv", delimiter=","); t3 = out["triplet"]
    return dict(train=t3["train"], validate=t3["validate"], test=t3["test"], dlod=out["dlod"], col4=a[:, [0, 3]].copy(), dir=d.name)
def search(name, i, p, ref=None):
    d = rundir(name, i, p); env = {**ENV[i], "CLIMATE_INDEX": f"{name}.dat", "TEST_ONLY": "false", "TIMEOUT": str(T), "NUMBER_OF_PROCESSORS": "5"}
    if ref is not None:
        np.savetxt(d / "manifold_ref.dat", ref, fmt="%.6f", delimiter="\t"); env["MANIFOLD"] = str(d / "manifold_ref.dat")
    r.run_lt(d, env, False, T); return json.loads((d / "lt.exe.p").read_text()), d
def nice(s): return f"{s['train']:.4f}/{s['validate']:.4f}/{s['test']:+.4f} dLOD {s['dlod']:.4f}"
state = {}
P0 = {i: json.loads((START / i / "lt.exe.p").read_text()) for i in IDX}
own = {i: score(rundir(f"own_{i}", i, P0[i]), i) for i in IDX}
for i in IDX: log(f"own fit {i}: {nice(own[i])}")
# Phase 1: noise floor
def rep(job):
    i, k = job; p, d = search(f"rep_{i}_{k}", i, P0[i]); return i, score(d, i)
with ThreadPoolExecutor(3) as ex: R = list(ex.map(rep, [(i, k) for i in IDX for k in range(3)]))
eps_i = {i: max(abs(s["validate"] - own[i]["validate"]) for j, s in R if j == i) for i in IDX}
for j, s in R: log(f"replicate {j}: {nice(s)}")
eps = max(max(eps_i.values()), 0.005); log(f"noise floor eps_i {eps_i} -> eps {eps:.4f}")
# Phase 2: round-robin
common = shared_of(P0["amo"]); P = {i: with_shared(P0[i], common) for i in IDX}
ev = {i: score(rundir(f"r00_{i}", i, P[i]), i) for i in IDX}
pen = {i: own[i]["validate"] - ev[i]["validate"] for i in IDX}; best = max(pen.values())
log(f"round 0 (common = AMO's manifold): " + " | ".join(f"{i} {nice(ev[i])} Delta {pen[i]:+.4f}" for i in IDX) + f" -> max Delta {best:+.4f}")
hist = [dict(round=0, accepted=True, pen=pen, ev={i: {k: v for k, v in ev[i].items() if k != 'col4'} for i in IDX})]; rejects = 0; exit_reason = "BUDGET"
if best <= eps and all(ev[i]["dlod"] >= DLOD_MIN for i in IDX): exit_reason = "SUCCESS"
for rd in range(1, ROUNDS + 1):
    if exit_reason == "SUCCESS": break
    X = IDX[(rd - 1) % 2]
    pX, _ = search(f"r{rd:02d}_search_{X}", X, P[X], ev[X]["col4"])
    cand = shared_of(pX); Pc = {i: (pX if i == X else with_shared(P[i], cand)) for i in IDX}
    evc = {i: score(rundir(f"r{rd:02d}_{i}", i, Pc[i]), i) for i in IDX}
    penc = {i: own[i]["validate"] - evc[i]["validate"] for i in IDX}; bc = max(penc.values())
    ok = all(evc[i]["dlod"] >= DLOD_MIN for i in IDX) and bc < best - eps
    log(f"round {rd} (search {X}): " + " | ".join(f"{i} {nice(evc[i])} Delta {penc[i]:+.4f}" for i in IDX) + f" -> max Delta {bc:+.4f} vs {best:+.4f}: {'ACCEPT' if ok else 'reject'}")
    hist.append(dict(round=rd, search=X, accepted=ok, pen=penc, ev={i: {k: v for k, v in evc[i].items() if k != 'col4'} for i in IDX}))
    if ok:
        P, ev, best, rejects = Pc, evc, bc, 0
        if best <= eps: exit_reason = "SUCCESS"
    else:
        rejects += 1
        if rejects >= 2: exit_reason = "PLATEAU"; break
    json.dump(dict(eps=eps, eps_i=eps_i, own={i: {k: v for k, v in own[i].items() if k != 'col4'} for i in IDX}, history=hist), open(B / "rr2_state.json", "w"), indent=1)
for i in IDX:
    (B / f"final_{i}").mkdir(exist_ok=True); src = W / ev[i]["dir"]
    for f in ("lt.exe.p", "lt.exe.resp", "lte_results.csv", "lt.exe.windings.json"):
        if (src / f).exists(): shutil.copy(src / f, B / f"final_{i}" / f)
json.dump(dict(eps=eps, eps_i=eps_i, exit=exit_reason, best_max_delta=best, own={i: {k: v for k, v in own[i].items() if k != 'col4'} for i in IDX}, history=hist), open(B / "rr2_state.json", "w"), indent=1)
log(f"EXIT {exit_reason}: max Delta {best:+.4f} (eps {eps:.4f}); final " + " | ".join(f"{i} {nice(ev[i])}" for i in IDX))
