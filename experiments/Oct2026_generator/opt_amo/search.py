#!/usr/bin/env python3
"""1200 s searches from the snapshot (start/) and from the most promising scan alternatives, equal budget:
 current fit (continuation) | impulse month index 5, asym x0.5 | impulse month index 8, asym x1.0 |
 second impulse delA = -0.1 x asym (IMPULSE_MOD12 on). YEAR = current (0.00405305) throughout."""
import sys, json
sys.argv=[sys.argv[0]]
from scan import run, P, show, par
def imp(label_m, sc): p=json.loads(json.dumps(P)); p["delB"]=((label_m-6)%12+18)/12.0; p["asym"]=P["asym"]*sc; return p
jobs=[("s_current",None,None,1200),("s_month5_half",imp(11,0.5),None,1200),("s_month8",imp(2,1.0),None,1200),
      ("s_second_m0.1",dict(P,delA=-0.1*P["asym"]),{"IMPULSE_MOD12":"TRUE"},1200)]
R=par(jobs,4)
for x in R: show(x, f"| ltep {[round(v,4) for v in x['p']['ltep']]} delA {x['p']['delA']:.3f} delB {x['p']['delB']:.4f} asym {x['p']['asym']:.3f}")
json.dump([{k:v for k,v in x.items()} for x in R],open("search.json","w"),indent=1); print("DONE")
