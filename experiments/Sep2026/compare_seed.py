#!/usr/bin/env python3
"""One-off guided comparison: nino4 as donor for kN000_W150, varying
METRIC (CC/DTW), VALIDATE (TRUE/FALSE), and IDATE (1880 matching nino4's
own anchor with a verbatim init copy, vs 1950 -- kN000_W150's own natural
start -- with nino4's own value AT 1950, derived via Bessel-inversion of
its live lte_results.csv). Reuses sweep.py's attempt()/run_attempt()
machinery directly. Short (90s) internal timeouts for a guided first
pass -- per the user's own framing, the goal is a reasonable comprehensive
answer, not a perfectly converged single-cell result."""
import sys, json, time
sys.path.insert(0, "/home/paul/eval/gem-lte-core/experiments/Sep2026")
import sweep

HERE = sweep.HERE
CELLS = {"A": HERE / "kN000_W150_A", "B": HERE / "kN000_W150_B"}

results = []
for label, cell_dir in CELLS.items():
    for metric in ("DTW", "CC"):
        for validate in (True, False):
            tag = f"{label}_{metric}_{'T' if validate else 'F'}"
            print(f"\n=== {tag} ===", flush=True)
            t0 = time.time()
            r = sweep.attempt(cell_dir, validate=validate, metric=metric, timeout_s=90)
            r["tag"] = tag
            r["elapsed_s"] = round(time.time() - t0, 1)
            print(json.dumps(r, indent=2), flush=True)
            results.append(r)

print("\n\n=== SUMMARY ===")
for r in results:
    t = r.get("triplet") or r.get("pair") or {}
    print(f"{r['tag']:>14s}  deadlocked={r['deadlocked']!s:5s}  "
          f"dlod={r.get('dlod')}  scores={t}")

(HERE / "compare_seed_results.json").write_text(json.dumps(results, indent=2))
