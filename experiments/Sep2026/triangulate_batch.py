#!/usr/bin/env python3
"""Triangulated (k=3, inverse-distance-weighted) re-attempts for the 8
stragglers left after the overnight run, plus a re-check of kS020_W090
(accepted but with suspicious validate=0.0/test=-1.87 scores -- excluded
from its own donor pool this time, so it can't seed from itself)."""
import sys, time
sys.path.insert(0, "/home/paul/eval/gem-lte-core/experiments/Sep2026")
import sweep

STRAGGLERS = ["kN060_W170", "kN040_W070", "kN020_E150", "kN000_W150",
              "kS020_E170", "kN000_E050", "kS020_E050", "kN040_W090"]
RECHECK = ["kS020_W090"]
TIMEOUT = 120
K = 3
RETRIES = 1


def run_triangulated(cell: str, exclude: set[str]) -> dict:
    result = None
    for i in range(RETRIES + 1):
        result = sweep.solve_cell_triangulated(cell, TIMEOUT, k=K, exclude=exclude)
        result["retry_index"] = i
        if result.get("accepted"):
            return result
    return result


for cell in STRAGGLERS + RECHECK:
    exclude = {cell} if cell in RECHECK else frozenset()
    neighbors = sweep.find_k_nearest_solved(cell, k=K, exclude=exclude)
    print(f"\n=== {cell} (triangulated, k={K}: {[(n, round(d)) for n, d in neighbors]}) ===",
          flush=True)
    t0 = time.time()
    result = run_triangulated(cell, exclude)
    result["timestamp"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    result["elapsed_s"] = round(time.time() - t0, 1)
    result["seed_from_path"] = "triangulated:" + ",".join(result.get("triangulated_from", []))
    sweep.append_ledger(result)
    if result.get("accepted"):
        print(f"[triangulate] SOLVED {cell}: {result['attempt_tag']} "
              f"dLOD={result['dlod']:.6f} scores={result['scores']}", flush=True)
    else:
        print(f"[triangulate] STILL UNSOLVED {cell}: {result.get('reason')}", flush=True)

print("\n[triangulate] done.")
