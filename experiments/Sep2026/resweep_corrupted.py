#!/usr/bin/env python3
"""Re-sweep the 8 SEVERELY-corrupted VALIDATE=TRUE cells found by
audit_validate_bug.py (2026-09-29) after fixing the Ada Forcing/harm
cross-thread resync bug in gem-lte-primitives-solution.adb.

Unlike sweep_lockw_batch.py, this does NOT reseed from a donor --
solve_cell()/setup_cell_dir() would discard the cell's own already-mostly-
correct fitted state (LPAP/LT frequencies were never corrupted by the bug,
only Forcing and the harm multipliers were). Instead it re-runs the normal
DTW-first cascade (sweep._run_cascade) IN PLACE, starting from the cell's
own current (corrupted-save) lt.exe.p, now with the fixed binary --
letting the search re-converge to a clean, self-consistent save.

Acceptance is NOT gated on beating the old ledger score: for a cell the
audit flagged as severely corrupted, the old recorded score is exactly
the untrustworthy one the bug produced (the "reported 0.6, replayed 0.3"
symptom) -- keeping the old, self-inconsistent checkpoint just because its
own (likely fictitious) reported number was higher would defeat the point.
Any cascade result that passes the normal gate_ok/backbone_ok checks is
committed; only a failed cascade (deadlock/no gate pass) rolls back to the
original corrupted state, which is preserved either way as a .prebugfix
backup.
"""
import sys
import time

sys.path.insert(0, "/home/paul/eval/gem-lte-core/experiments/Sep2026")
import sweep

# kN000_E150, kN040_W050, kS020_E150, and kS040_E170 already re-swept and
# committed (2026-09-29). kS020_E150 needed a THIRD pass after a deeper
# Ada bug was found: Forcing was missing entirely from
# Monitor.Report_Final/Winner (the VALIDATE=TRUE cross-thread lockbox),
# so the earlier KeepForcing fix only resynced within a single thread,
# not across threads -- now fixed and verified (kS020_E150 audit
# rel=9.3e-10, corr=1.0). These 4 are the remaining cells from the
# original "severe 8" not yet re-swept with the fully-fixed binary.
TARGET_CELLS = ['kS040_E070', 'kS040_E050', 'kN020_E150', 'kS040_E150']

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


def audit_one(cell: str) -> dict | None:
    """Re-check this one cell's post-cascade Forcing-vs-lt.exe.p agreement,
    reusing audit_validate_bug's own logic."""
    import audit_validate_bug as A
    return A.audit_cell(sweep.HERE / cell)


def main() -> None:
    by_cell = sweep.latest_entry_per_cell(sweep.load_sweep_ledger())
    results = {}

    for cell in TARGET_CELLS:
        print(f"\n=== {cell} ===", flush=True)
        cell_dir = sweep.HERE / cell
        old_entry = by_cell.get(cell)
        old_min = min_score(old_entry)
        print(f"  old ledger min score: {old_min:.3f} (may be inflated -- "
              f"this is one of the audit-confirmed severely-corrupted cells)",
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
                  f"{elapsed:.0f}s -- rolling back to pre-bugfix state",
                  flush=True)
            restore(cell, saved)
            results[cell] = dict(action="rolled_back_cascade_failed",
                                  old_min=old_min)
            continue

        new_min = min_score(result)
        bb = sweep.read_backbone(cell_dir)
        audit_after = audit_one(cell)
        print(f"  cascade succeeded after {elapsed:.0f}s: {result['attempt_tag']} "
              f"scores={result['scores']}  backbone={bb}", flush=True)
        print(f"  new min score: {new_min:.3f}  (old ledger: {old_min:.3f})",
              flush=True)
        if audit_after:
            print(f"  post-fix audit: rel={audit_after.get('rel'):.2e} "
                  f"corr={audit_after.get('corr'):.9f}", flush=True)

        result["timestamp"] = time.strftime("%Y-%m-%dT%H:%M:%S")
        result["note"] = (
            f"Re-swept in place after fixing the Ada Forcing/harm cross-thread "
            f"resync bug (audit_validate_bug.py flagged this cell as severely "
            f"corrupted, corr={((by_cell.get(cell) or {}).get('audit_corr'))}). "
            f"Old (untrustworthy) ledger min={old_min:.3f}, new min={new_min:.3f}. "
            f"backbone={bb}. Committed regardless of the old/new score "
            f"comparison since the old score is exactly the bug's own "
            f"fabricated one, not a real baseline to protect."
        )
        sweep.append_ledger(result)
        print(f"  COMMITTED", flush=True)
        results[cell] = dict(action="committed", old_min=old_min,
                              new_min=new_min, backbone=bb,
                              post_fix_audit=audit_after)

    print("\n=== SUMMARY ===", flush=True)
    for cell, r in results.items():
        print(f"  {cell}: {r}", flush=True)
    (sweep.HERE / "RESWEEP_CORRUPTED_DONE").write_text(
        f"{time.strftime('%Y-%m-%dT%H:%M:%S')}\n" +
        "\n".join(f"{c}: {r}" for c, r in results.items()))


if __name__ == "__main__":
    main()
