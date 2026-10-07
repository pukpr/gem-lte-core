#!/usr/bin/env python3
"""Equal-budget (1200 s) refits under the corrected tidal sum (sin throughout, jerk x Amplitude):
TABLE from the 0.833 TABLE fit, ZONAL from the user's fit. Same resp/run settings (snapshot start/)."""
import sys, json
sys.argv=[sys.argv[0]]
from scan import run, show, par
TP=json.load(open("../savebug/imp_off/lt.exe.p"))
jobs=[("fix_table",TP,{"TIDES":"TABLE","TIDE_LEGACY":"FALSE"},1200),("fix_zonal",None,{"TIDES":"ZONAL","TIDE_LEGACY":"FALSE"},1200)]
R=par(jobs,2)
for x in R: show(x, f"| ltep {[round(v,4) for v in x['p']['ltep']]} shfT {x['p']['shfT']:.5f}")
json.dump(R,open("refit.json","w"),indent=1); print("DONE")
