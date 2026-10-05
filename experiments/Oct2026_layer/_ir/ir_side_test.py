#!/usr/bin/env python3
"""Warm-side-only IR (IR_SIDE=+1), build 2026-10-03 11:47, same 10 targets as ir_test.py.
Grid: IR_SIDE=+1 with IR in IRS; references IR=-0.2 symmetric (IR_SIDE 0) and cold-only (-1).
Searches: 300 s seeded IR=-0.2, IR_SIDE=+1; compared with the IR=0 control searches in ir.json."""
import json, sys
sys.argv = [sys.argv[0]]
from ir_test import *
sys.path.insert(0, str(E / "Sep2026")); import sweep
def mk(name, t, ir, side):
    d = make(name, t, ir); sweep.set_resp_key(d / "lt.exe.resp", "IR_SIDE", str(side)); return d
IRS = (-0.4, -0.3, -0.2, -0.1, 0.0)
jobs = [(t, v, 1) for t in TARGETS for v in IRS] + [(t, -0.2, s) for t in TARGETS for s in (0, -1)]
with ThreadPoolExecutor(8) as ex:
    G = list(ex.map(lambda j: dict(side=j[2], **score(j[0], mk(f"gs_{j[0]}_ir{j[1]:+g}_s{j[2]}", j[0], j[1], j[2]))), jobs))
print("GRID (fixed params): warm-only IR: train/validate warm/cold   || IR=-0.2 symmetric | cold-only")
for t in TARGETS:
    w = [x for x in G if x["tgt"] == t and x["side"] == 1]; s0 = next(x for x in G if x["tgt"] == t and x["side"] == 0); s1 = next(x for x in G if x["tgt"] == t and x["side"] == -1)
    f = lambda x: f"{x['train']:.3f}/{x['validate']:.3f} w{x['warm']:.2f}c{x['cold']:.2f}"
    print(f"  {t:11s} " + "  ".join(f"{x['IR']:+.1f}:{f(x)}" for x in w) + f"  || sym {f(s0)} | cold {f(s1)}", flush=True)
with ThreadPoolExecutor(4) as ex:
    S = list(ex.map(lambda t: dict(side=1, **score(t, mk(f"ss_{t}_warm", t, -0.2, 1), search=300)), list(TARGETS)))
ctrl = {x["tgt"]: x for x in json.load(open(B / "ir.json"))["search"] if x["name"].endswith("_ir+0")}
print("SEARCHES 300 s: IR=0 control | warm-only seeded IR=-0.2")
f = lambda x: f"{x['train']:.3f}/{x['validate']:.3f}/{x['test']:+.3f} IR {x['IR']:+.3f} warm {x['warm']:.2f} cold {x['cold']:.2f} amp {x['amp']:.2f}"
for x in S: print(f"  {x['tgt']:11s} {f(ctrl[x['tgt']])} | {f(x)}")
json.dump(dict(grid=G, search=S), open(B / "ir_side.json", "w"), indent=1); print("DONE")
