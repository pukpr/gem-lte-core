#!/usr/bin/env python3
"""Pilot test (2026-09-29, per user request): does seeding ImpC=0.01 (the
Annual_Impulse term, structurally stuck at 0.0 by gem-random_descent.adb's
Markov procedure in every prior fit) improve the fit once the search can
actually explore it? Re-runs the normal DTW-first cascade in place on a
small, diverse pilot set (one already-negative-test-score cell, one weak
pair-metric cell, one already-strong cell) and reports before/after scores
plus where ImpC ended up (near its 0.01 seed = probably not useful; grown
meaningfully in either direction = probably real signal).

Commits (via sweep.append_ledger) only if the new result passes the normal
gate_ok/backbone_ok checks AND its min score is >= the old one -- this is
an exploratory test, not a known-corrupted-cell fix, so unlike
resweep_corrupted.py it does NOT commit regardless of score comparison.
"""
import json
import sys
import time

sys.path.insert(0, "/home/paul/eval/gem-lte-core/experiments/Sep2026")
import sweep

PILOT_CELLS = ['kN060_W050', 'kN000_E110']
TIMEOUT = 400


def min_score(entry) -> float:
    scores = (entry or {}).get("scores") or {}
    if "train" in scores:
        return min(scores["train"], scores["validate"], scores["test"])
    if "pair" in scores:
        return min(scores["pair"])
    return -999.0


def backup(cell: str) -> dict:
    d = sweep.HERE / cell
    saved = {}
    for fname in ("lt.exe.p", "lt.exe.resp", "lte_run.sh", "lte_results.csv",
                  "lt.exe.windings.json"):
        p = d / fname
        if p.exists():
            saved[fname] = p.read_bytes()
    return saved


def restore(cell: str, saved: dict) -> None:
    d = sweep.HERE / cell
    for fname, content in saved.items():
        (d / fname).write_bytes(content)


def main() -> None:
    by_cell = sweep.latest_entry_per_cell(sweep.load_sweep_ledger())
    results = {}

    for cell in PILOT_CELLS:
        print(f"\n=== {cell} ===", flush=True)
        cell_dir = sweep.HERE / cell
        old_entry = by_cell.get(cell)
        old_min = min_score(old_entry)
        impc_before = json.loads((cell_dir / "lt.exe.p").read_text()).get("impC")
        print(f"  old min score: {old_min:.3f}  impC before: {impc_before}",
              flush=True)

        saved = backup(cell)
        t0 = time.time()
        try:
            result = sweep._run_cascade(cell, cell_dir, TIMEOUT)
        except Exception as exc:
            print(f"  EXCEPTION: {exc}", flush=True)
            restore(cell, saved)
            results[cell] = dict(action="exception", old_min=old_min)
            continue
        elapsed = time.time() - t0

        if not result.get("accepted"):
            print(f"  cascade FAILED ({result.get('reason')}) after "
                  f"{elapsed:.0f}s -- rolling back", flush=True)
            restore(cell, saved)
            results[cell] = dict(action="rolled_back_cascade_failed",
                                  old_min=old_min)
            continue

        new_min = min_score(result)
        impc_after = json.loads((cell_dir / "lt.exe.p").read_text()).get("impC")
        print(f"  cascade succeeded after {elapsed:.0f}s: {result['attempt_tag']} "
              f"scores={result['scores']}", flush=True)
        print(f"  new min score: {new_min:.3f}  (old: {old_min:.3f})  "
              f"impC after: {impc_after}", flush=True)

        if new_min >= old_min:
            result["timestamp"] = time.strftime("%Y-%m-%dT%H:%M:%S")
            result["impC"] = impc_after
            result["note"] = (
                f"ImpC=0.01 seeding pilot (2026-09-29): min score "
                f"{old_min:.3f} -> {new_min:.3f}, impC {impc_before} -> "
                f"{impc_after}."
            )
            sweep.append_ledger(result)
            print(f"  COMMITTED (delta {new_min - old_min:+.3f})", flush=True)
            results[cell] = dict(action="committed", old_min=old_min,
                                  new_min=new_min, impc_before=impc_before,
                                  impc_after=impc_after)
        else:
            print(f"  did NOT improve -- rolling back", flush=True)
            restore(cell, saved)
            results[cell] = dict(action="rolled_back_no_improvement",
                                  old_min=old_min, new_min=new_min,
                                  impc_before=impc_before,
                                  impc_after=impc_after)

    print("\n=== SUMMARY ===", flush=True)
    for cell, r in results.items():
        print(f"  {cell}: {r}", flush=True)
    (sweep.HERE / "PILOT_IMPC_DONE").write_text(
        f"{time.strftime('%Y-%m-%dT%H:%M:%S')}\n" +
        "\n".join(f"{c}: {r}" for c, r in results.items()))


if __name__ == "__main__":
    main()
