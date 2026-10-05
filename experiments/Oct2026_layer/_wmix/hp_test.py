#!/usr/bin/env python3
"""LAYER=1 + high-passed mixing (LAYER_HP months), build 2026-10-03 10:35: F2 = w*(Fb - runmean_N(Fb)) + layer,
w = impC. Grid TEST_ONLY over N and w on amo, pdo, kN000_W090/W110/W130; then a 300 s search from each target's
best grid point (by mean of train and validate), compared with the earlier w=0 searches (s_*_w+0)."""
import json, sys
sys.argv = [sys.argv[0]]
from wmix_test import *
sys.path.insert(0, str(E / "Sep2026")); import sweep
NS = (13, 61); WS = (-0.05, -0.02, -0.01, -0.005, -0.002, 0.0, 0.002, 0.005, 0.01, 0.02, 0.05)
def mk(name, t, w, n):
    d = make(name, t, w); sweep.set_resp_key(d / "lt.exe.resp", "LAYER_HP", str(n)); return d
with ThreadPoolExecutor(8) as ex:
    G = list(ex.map(lambda j: dict(N=j[2], **score(j[0], mk(f"h_{j[0]}_N{j[2]}_w{j[1]:+g}", j[0], j[1], j[2]))), [(t, w, n) for t in TARGETS for n in NS for w in WS]))
for t in TARGETS:
    for n in NS:
        print(f"{t:11s} N={n:2d}: " + "  ".join(f"{x['w']:+.3f}:{x['train']:.3f}/{x['validate']:.3f} w{x['warm']:.2f}c{x['cold']:.2f}" for x in G if x["tgt"] == t and x["N"] == n), flush=True)
jobs = []
for t in TARGETS:
    best = max([x for x in G if x["tgt"] == t and x["w"] != 0.0], key=lambda x: (x["train"] + x["validate"]) / 2); jobs.append((t, best["w"], best["N"]))
with ThreadPoolExecutor(4) as ex:
    S = list(ex.map(lambda j: dict(N=j[2], **score(j[0], mk(f"hs_{j[0]}_N{j[2]}_w{j[1]:+g}", j[0], j[1], j[2]), search=300)), jobs))
print("SEARCHES 300 s from best nonzero grid point (earlier w=0 search for reference):")
ref = {x["tgt"]: x for x in json.load(open(B / "wmix.json"))["search"] if x["name"].endswith("_w+0")}
for x, j in zip(S, jobs):
    r0 = ref[x["tgt"]]
    print(f"  {x['tgt']:11s} N={j[2]} start w {j[1]:+.3f} -> w {x['w']:+.4f}: {x['train']:.4f}/{x['validate']:.4f}/{x['test']:.4f} dLOD {x['dlod']:.4f} warm {x['warm']:.2f} cold {x['cold']:.2f}"
          f"  | w=0 search: {r0['train']:.4f}/{r0['validate']:.4f}/{r0['test']:.4f} warm {r0['warm']:.2f} cold {r0['cold']:.2f}")
json.dump(dict(grid=[{k: v for k, v in x.items() if k != 'p'} for x in G], search=[{k: v for k, v in x.items() if k != 'p'} for x in S]), open(B / "hp.json", "w"), indent=1)
print("DONE")
