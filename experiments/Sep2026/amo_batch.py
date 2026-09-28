#!/usr/bin/env python3
"""Batch: seed Atlantic grid cells near AMO from AMO's own manifold.
Per the user's guidance: prioritize breadth (get as many solved as
possible) over perfecting any one stuck cell -- unsolved ones are left
for later triangulation/kriging from whichever neighbors DID solve."""
import sys, json, time
sys.path.insert(0, "/home/paul/eval/gem-lte-core/experiments/Sep2026")
import sweep

CELLS = ["kN040_W030", "kN040_W050", "kN060_W030", "kN020_W050",
         "kN040_W010", "kN040_W070"]
donor_dir = sweep.FEB / "amo"
TIMEOUT = 120

for cell in CELLS:
    print(f"\n=== {cell} (seed: amo) ===", flush=True)
    t0 = time.time()
    result = solve = sweep.solve_cell_with_retry(cell, donor_dir, TIMEOUT, retries=1)
    result["seed_from_path"] = str(donor_dir)
    result["timestamp"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    result["elapsed_s"] = round(time.time() - t0, 1)
    sweep.append_ledger(result)
    if result.get("accepted"):
        print(f"[amo_batch] SOLVED {cell}: {result['attempt_tag']} "
              f"dLOD={result['dlod']:.6f} scores={result['scores']}", flush=True)
    else:
        print(f"[amo_batch] UNSOLVED {cell}: {result.get('reason')} "
              f"(after {result.get('retry_index', 0) + 1} full tries)", flush=True)

print("\n[amo_batch] done.")
