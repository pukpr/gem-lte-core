#!/usr/bin/env python3
"""Re-sweep the weakest-scoring currently-accepted cells using the new
LOCKW canonical-winding lock (2026-09-28), donor-seeded from each cell's
nearest already-solved neighbor. Per the validated usage pattern (see
HANDOFF.md): LOCKW works as a drift-prevention guardrail on a reasonable
starting point, not a cold-start forcer -- donor-seeding supplies that
starting point, LOCKW keeps the search from drifting off it to a
spurious alias during refinement.

For each cell: back up current state, run the normal DTW-first cascade
(sweep._run_cascade) but with LOCKW=TRUE added to BASE_OVERRIDES for the
duration, judge success by whether the new result's minimum score beats
the old one. Roll back if not.
"""
import sys
import time

sys.path.insert(0, "/home/paul/eval/gem-lte-core/experiments/Sep2026")
import sweep

TARGET_CELLS = ['kS020_W090', 'kN040_W010', 'kN020_W090', 'kN020_W030',
                 'kS040_W010', 'kS040_W070', 'kN060_W130', 'kS040_E170']

TIMEOUT = 400  # longer than the usual 120s -- LOCKW needs more iterations
                # to find a good fit around a locked winding, per this
                # session's direct finding


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
    for fname in ("lt.exe.p", "lt.exe.resp", "lte_run.sh"):
        p = d / fname
        if p.exists():
            saved[fname] = p.read_bytes()
    return saved


def restore(cell: str, saved: dict) -> None:
    d = sweep.HERE / cell
    for fname, content in saved.items():
        (d / fname).write_bytes(content)


def main() -> None:
    sweep.BASE_OVERRIDES["LOCKW"] = "TRUE"
    sweep.BASE_OVERRIDES["CANON_BACKBONE"] = "0.207"
    sweep.BASE_OVERRIDES["CANON_HARMONICS"] = \
        "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 18 19 23 27"
    sweep.BASE_OVERRIDES["CANON_SUBHARMONIC"] = "11"

    by_cell = sweep.latest_entry_per_cell(sweep.load_sweep_ledger())
    results = {}

    for cell in TARGET_CELLS:
        print(f"\n=== {cell} ===", flush=True)
        old_entry = by_cell.get(cell)
        old_min = min_score(old_entry)
        print(f"  current min score: {old_min:.3f}", flush=True)

        lat, lon = sweep.parse_grid_name(cell)
        nodes = sweep.current_solved_nodes()
        donor_name, donor_dist = min(
            ((n, sweep.haversine_km(lat, lon, *co)) for n, co in nodes.items()
             if n != cell), key=lambda x: x[1])
        donor_dir = sweep.FEB / donor_name if donor_name in sweep.FLAGSHIPS else sweep.HERE / donor_name
        print(f"  donor: {donor_name} ({donor_dist:.0f} km)", flush=True)

        saved = backup(cell)
        t0 = time.time()
        try:
            result = sweep.solve_cell_with_retry(cell, donor_dir, TIMEOUT, retries=1)
        except Exception as exc:
            print(f"  EXCEPTION: {exc}", flush=True)
            restore(cell, saved)
            results[cell] = dict(action="exception", old_min=old_min)
            continue
        elapsed = time.time() - t0

        if not result.get("accepted"):
            print(f"  cascade FAILED ({result.get('reason')}) after {elapsed:.0f}s -- rolling back", flush=True)
            restore(cell, saved)
            results[cell] = dict(action="rolled_back_cascade_failed", old_min=old_min)
            continue

        new_min = min_score(result)
        bb = sweep.read_backbone(sweep.HERE / cell)
        print(f"  cascade succeeded after {elapsed:.0f}s: {result['attempt_tag']} "
              f"scores={result['scores']}  backbone={bb}", flush=True)
        print(f"  new min score: {new_min:.3f}  (old: {old_min:.3f})", flush=True)

        if new_min > old_min:
            result["seed_from_path"] = str(donor_dir)
            result["timestamp"] = time.strftime("%Y-%m-%dT%H:%M:%S")
            result["note"] = (
                f"LOCKW canonical-backbone re-sweep (donor-seeded from {donor_name}, "
                f"{donor_dist:.0f}km). Min score improved {old_min:.3f} -> {new_min:.3f}. "
                f"backbone={bb}."
            )
            sweep.append_ledger(result)
            print(f"  COMMITTED (improved by {new_min - old_min:.3f})", flush=True)
            results[cell] = dict(action="committed", old_min=old_min, new_min=new_min, backbone=bb)
        else:
            print(f"  did NOT improve -- rolling back, keeping original", flush=True)
            restore(cell, saved)
            results[cell] = dict(action="rolled_back_no_improvement", old_min=old_min, new_min=new_min)

    print("\n=== SUMMARY ===", flush=True)
    for cell, r in results.items():
        print(f"  {cell}: {r}", flush=True)
    (sweep.HERE / "LOCKW_BATCH_DONE").write_text(
        f"{time.strftime('%Y-%m-%dT%H:%M:%S')}\n" + "\n".join(f"{c}: {r}" for c, r in results.items()))


if __name__ == "__main__":
    main()
