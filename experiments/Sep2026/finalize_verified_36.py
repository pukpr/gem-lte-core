#!/usr/bin/env python3
"""Finalize the 34 (of 36) never-attempted cells that verify_untouched_36.py
found usable from their own prior Feb2026 parameter set: write the winning
metric into lte_run.sh (TEST_ONLY=false), sync lt.exe.p/resp/lte_run.sh to
Sep2026/<cell>/, and record a ledger entry for each."""
import json
import re
import sys
import time

sys.path.insert(0, "/home/paul/eval/gem-lte-core/experiments/Sep2026")
import sweep

report = json.load(open(sweep.HERE / "verify_untouched_36_report.json"))

finalized = []
for cell, e in sorted(report.items()):
    if e.get("needs_donor_seed"):
        continue
    winner = e["winner"]
    res = e[winner.lower()]
    cell_dir = sweep.FEB / cell

    # Finalize lte_run.sh: winning metric, TEST_ONLY=false
    run_sh = cell_dir / "lte_run.sh"
    text = run_sh.read_text()
    text = text.replace("export TEST_ONLY=true", "export TEST_ONLY=false")
    text = re.sub(r"export METRIC=\S+", f"export METRIC={winner}", text)
    run_sh.write_text(text)

    # Sync to Sep2026/<cell>/
    target = sweep.HERE / cell
    target.mkdir(exist_ok=True)
    for fname in ("lt.exe.p", "lt.exe.resp", "lte_run.sh"):
        src = cell_dir / fname
        if src.exists():
            (target / fname).write_bytes(src.read_bytes())

    # Build ledger entry
    if res.get("triplet"):
        t = res["triplet"]
        scores = dict(metric=winner, train=t["train"], validate=t["validate"], test=t["test"])
        validate_mode = True
    else:
        p = res["pair"]
        scores = dict(metric=winner, pair=[p["val1"], p["val2"]])
        validate_mode = False

    sweep.append_ledger(dict(
        cell=cell, accepted=True, attempt_tag="manual-gold-standard",
        validate_mode=validate_mode, metric=winner, scores=scores,
        dlod=res["dlod"], backbone=e["backbone"], n_attempts=None,
        seed_from_path=(
            f"manual: reused this cell's OWN prior Feb2026 parameter set "
            f"(never touched by sweep.py before, no ledger entry existed) "
            f"per user instruction (2026-09-28) to try prior work before "
            f"donor-seeding. Batch-verified DTW vs CC via "
            f"verify_untouched_36.py; {winner} won (quality={e['winner_quality']})."
        ),
        timestamp=time.strftime("%Y-%m-%dT%H:%M:%S"),
        note="Part of the 34/36 batch that succeeded from prior work alone, no search needed.",
    ))
    finalized.append(cell)
    print(f"finalized {cell}: {winner} {scores}")

print(f"\n=== finalized {len(finalized)} cells ===")
by_cell = sweep.latest_entry_per_cell(sweep.load_sweep_ledger())
solved = [c for c, e in by_cell.items() if e.get("accepted")]
print(f"total solved now: {len(solved)}")
