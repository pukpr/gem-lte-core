#!/usr/bin/env python3
"""Fix the 3 known near-miss/wrong-donor cells before continuing the
outward flood-fill: re-seed the two pna-orphaned cells from the nearer,
oceanic kN040_W170; retry kN060_E010 (pure stochastic near-miss last
time) from baltic again."""
import sys, time
sys.path.insert(0, "/home/paul/eval/gem-lte-core/experiments/Sep2026")
import sweep

JOBS = [
    ("kN060_W170", sweep.HERE / "kN040_W170"),
    ("kN040_W150", sweep.HERE / "kN040_W170"),
    ("kN060_E010", sweep.FEB / "baltic"),
]
TIMEOUT = 120

for cell, donor_dir in JOBS:
    print(f"\n=== {cell} (re-seed: {donor_dir.name}) ===", flush=True)
    t0 = time.time()
    result = sweep.solve_cell_with_retry(cell, donor_dir, TIMEOUT, retries=1)
    result["seed_from_path"] = str(donor_dir)
    result["timestamp"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    result["elapsed_s"] = round(time.time() - t0, 1)
    sweep.append_ledger(result)
    if result.get("accepted"):
        print(f"[fixups] SOLVED {cell}: {result['attempt_tag']} "
              f"dLOD={result['dlod']:.6f} scores={result['scores']}", flush=True)
    else:
        print(f"[fixups] STILL UNSOLVED {cell}: {result.get('reason')}", flush=True)

print("\n[fixups] done.")
