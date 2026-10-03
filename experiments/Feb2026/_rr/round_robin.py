#!/usr/bin/env python3
"""Round-robin AMO <-> PDO with a COMMON manifold (2026-10-02, LAYER=1 build).

Shared (everything that sets lte_results.csv col 4 and the dLOD): lpap, delA, delB, asym,
ma, mp, shfT, init, year, impA, impB, offs, bg, ltep[0] (layer k1), ltep[1] (Bessel k).
Per index: ltep[2..] (regressed windings) + harm, impC, IR, annual terms (regressed).
Round i: search index X (TIMEOUT s) from the current common manifold + X's own private
parameters, with MANIFOLD=<current common col 4> as a soft anchor; take X's shared
parameters, re-score BOTH indices with them (TEST_ONLY, MLR re-solved); accept if the
joint score (mean training CC) improves and dLOD >= DLOD_MIN for both, else keep the old
common manifold (X's private parameters are kept from its search only on accept).
Works on copies in _rr/work; Feb2026/amo and Feb2026/pdo are never modified.
Usage: round_robin.py [rounds] [timeout]"""
import json, re, shutil, sys, time
from pathlib import Path
import numpy as np
B = Path(__file__).resolve().parent; FEB = B.parent; E = FEB.parent; W = B / "work"
sys.path.insert(0, str(E / "Sep2026")); import sweep, sweep_cc_refine as r
ROUNDS = int(sys.argv[1]) if len(sys.argv) > 1 else 10; TIMEOUT = int(sys.argv[2]) if len(sys.argv) > 2 else 300
DLOD_MIN = 0.994; IDX = ("amo", "pdo")
SHARED = ("lpap", "delA", "delB", "asym", "ma", "mp", "shfT", "init", "year", "impA", "impB", "offs", "bg")
for f, tgt in (("lt.exe", FEB / "enso_opt"), ("dlod3.dat", FEB / "dlod3.dat")):
    if not (W / f).exists(): W.mkdir(exist_ok=True); (W / f).symlink_to(tgt)
ENV = {i: dict(re.findall(r'^export (\w+)=(\S+)', (FEB / i / "lte_run.sh").read_text(), re.M)) for i in IDX}
def log(m):
    print(f"{time.strftime('%H:%M:%S')} {m}", flush=True)
def shared_of(p): return {**{k: p[k] for k in SHARED}, "ltep01": p["ltep"][:2]}
def with_shared(p, s):
    q = json.loads(json.dumps(p))
    for k in SHARED: q[k] = s[k]
    q["ltep"] = list(s["ltep01"]) + q["ltep"][2:]; return q
def rundir(name, idx, p):
    d = W / name
    if d.exists(): shutil.rmtree(d)
    d.mkdir()
    for f in ("lt.exe.resp", "dlod_ref.dat"): shutil.copy(FEB / idx / f, d / f)
    shutil.copy(FEB / idx / f"{idx}.dat", d / f"{name}.dat"); (d / "lt.exe.p").write_text(json.dumps(p, indent=2)); return d
def evaluate(name, idx, p):
    d = rundir(name, idx, p)
    out = r.parse(r.run_lt(d, {**ENV[idx], "CLIMATE_INDEX": f"{name}.dat", "TEST_ONLY": "true", "TIMEOUT": "30", "NUMBER_OF_PROCESSORS": "2"}, True, 30))
    a = np.loadtxt(d / "lte_results.csv", delimiter=",")
    return dict(triplet=out["triplet"], dlod=out["dlod"], col4=a[:, :4:3].copy(), dir=d.name)
def search(name, idx, p, ref_col4):
    d = rundir(name, idx, p); ref = d / "manifold_ref.dat"
    np.savetxt(ref, ref_col4, fmt="%.6f", delimiter="\t")
    r.run_lt(d, {**ENV[idx], "CLIMATE_INDEX": f"{name}.dat", "MANIFOLD": str(ref), "TEST_ONLY": "false", "TIMEOUT": str(TIMEOUT), "NUMBER_OF_PROCESSORS": "8"}, False, TIMEOUT)
    return json.loads((d / "lt.exe.p").read_text())
def joint(ev): return float(np.mean([ev[i]["triplet"]["train"] for i in IDX]))
def fmt(ev): return " | ".join(f"{i} {ev[i]['triplet']['train']:.4f}/{ev[i]['triplet'].get('validate',0):.4f}/{ev[i]['triplet']['test']:.4f} dLOD {ev[i]['dlod']:.4f}" for i in IDX)
P = {i: json.loads((FEB / i / "lt.exe.p").read_text()) for i in IDX}
hist = []
ev0 = {i: evaluate(f"orig_{i}", i, P[i]) for i in IDX}; log(f"own manifolds (as seeded): {fmt(ev0)}")
common = shared_of(P["amo"]); P = {i: with_shared(P[i], common) for i in IDX}
ev = {i: evaluate(f"r00_{i}", i, P[i]) for i in IDX}; best = joint(ev)
log(f"round 0, common = AMO's manifold: {fmt(ev)}  joint {best:.4f}"); hist.append(dict(round=0, src="amo-seed", accepted=True, joint=best, ev={i: {k: v for k, v in ev[i].items() if k != 'col4'} for i in IDX}))
for rd in range(1, ROUNDS + 1):
    X = IDX[(rd - 1) % 2]
    pX = search(f"r{rd:02d}_search_{X}", X, P[X], ev[X]["col4"])
    cand = shared_of(pX); Pc = {X: pX, **{i: with_shared(P[i], cand) for i in IDX if i != X}}
    evc = {i: evaluate(f"r{rd:02d}_{i}", i, Pc[i]) for i in IDX}; jc = joint(evc)
    ca, cp = evc["amo"]["col4"], evc["pdo"]["col4"]; n = min(len(ca), len(cp))   # same start (1880), PDO ends earlier
    m4 = float(np.corrcoef(ca[:n, 1], cp[:n, 1])[0, 1])
    ok = jc > best and all(evc[i]["dlod"] >= DLOD_MIN for i in IDX)
    log(f"round {rd} (search {X}): {fmt(evc)}  joint {jc:.4f} vs {best:.4f}  col4 CC amo/pdo {m4:.5f} -> {'ACCEPT' if ok else 'reject'}")
    hist.append(dict(round=rd, src=X, accepted=ok, joint=jc, ev={i: {k: v for k, v in evc[i].items() if k != 'col4'} for i in IDX}))
    if ok: P, ev, best, common = Pc, evc, jc, cand
    json.dump(dict(history=hist, best_joint=best, best_params=P), open(B / "rr_state.json", "w"), indent=1)
for i in IDX:
    (B / f"best_{i}").mkdir(exist_ok=True); src = W / ev[i]["dir"]
    for f in ("lt.exe.p", "lt.exe.resp", "lte_results.csv", "lt.exe.windings.json"): shutil.copy(src / f, B / f"best_{i}" / f)
log(f"DONE best joint {best:.4f}: {fmt(ev)}")
