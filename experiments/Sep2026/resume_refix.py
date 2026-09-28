#!/usr/bin/env python3
"""Resume the backbone-invalidation refix at N=2 processors.

invalidate_and_refix.py's own find_bad_cells() only catches cells whose
LATEST ledger entry is still accepted=True with a bad backbone -- once a
cell has been invalidated (latest entry accepted=False, reason
"invalidated-bad-backbone" or a subsequent failed refix attempt
"unsolved-needs-review"), a bare re-run of that script finds zero bad
cells and does nothing. This script instead resumes the specific 24-cell
backlog directly (recomputed from the ledger by "still needs refix":
latest entry accepted=False with one of those two reasons), reusing the
exact same nearest-donor flood-fill loop, now correctly at
NUMBER_OF_PROCESSORS=2 (sweep.py's BASE_OVERRIDES was reverted from 1
back to 2 -- see HANDOFF.md)."""
import sys, time, traceback
sys.path.insert(0, "/home/paul/eval/gem-lte-core/experiments/Sep2026")
import sweep

TIMEOUT = 120
RETRIES = 1
LOG = sweep.HERE / "refix.log"

REFIX_REASONS = {"invalidated-bad-backbone", "unsolved-needs-review"}


def log(msg: str) -> None:
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)


def find_still_bad_cells() -> list[str]:
    by_cell = sweep.latest_entry_per_cell(sweep.load_sweep_ledger())
    return [c for c, e in by_cell.items()
            if not e.get("accepted") and e.get("reason") in REFIX_REASONS]


def main() -> None:
    bad_cells = find_still_bad_cells()
    log(f"=== resuming refix at N=2: {len(bad_cells)} cells still bad: {bad_cells} ===")

    remaining = set(bad_cells)
    total = len(bad_cells)
    processed = 0
    while remaining:
        nodes = sweep.current_solved_nodes()

        def dist(c):
            lat, lon = sweep.parse_grid_name(c)
            return min(sweep.haversine_km(lat, lon, *co) for co in nodes.values())

        cell = min(remaining, key=dist)
        remaining.discard(cell)
        lat, lon = sweep.parse_grid_name(cell)
        donor_name, donor_dist = min(
            ((n, sweep.haversine_km(lat, lon, *co)) for n, co in nodes.items()),
            key=lambda x: x[1])
        donor_dir = (sweep.FEB / donor_name if donor_name in sweep.FLAGSHIPS
                     else sweep.HERE / donor_name)

        log(f"--- [{processed+1}/{total}] {cell}  "
            f"seed={donor_name} ({donor_dist:.0f} km, single-donor) ---")
        t0 = time.time()
        try:
            result = sweep.solve_cell_with_retry(cell, donor_dir, TIMEOUT, retries=RETRIES)
            if not result.get("accepted"):
                log(f"    single-donor failed, trying triangulated (k=3)...")
                result = sweep.solve_cell_triangulated_with_retry(
                    cell, TIMEOUT, k=3, retries=RETRIES)
                donor_dir = "triangulated:" + ",".join(result.get("triangulated_from", []))
        except Exception as exc:
            log(f"EXCEPTION on {cell}: {exc}")
            log(traceback.format_exc())
            result = dict(cell=cell, accepted=False, reason=f"exception: {exc}")

        result["seed_from_path"] = str(donor_dir)
        result["timestamp"] = time.strftime("%Y-%m-%dT%H:%M:%S")
        result["elapsed_s"] = round(time.time() - t0, 1)
        sweep.append_ledger(result)

        if result.get("accepted"):
            bb = sweep.read_backbone(sweep.HERE / cell)
            log(f"REFIXED {cell}: {result['attempt_tag']} dLOD={result['dlod']:.6f} "
                f"backbone={bb} scores={result.get('scores')}")
        else:
            log(f"STILL BAD {cell}: {result.get('reason')}")
        processed += 1

    by_cell = sweep.latest_entry_per_cell(sweep.load_sweep_ledger())
    solved = [c for c, e in by_cell.items() if e.get("accepted")]
    still_bad = find_still_bad_cells()
    log(f"=== resume-refix complete: {len(solved)} total solved, "
        f"{len(still_bad)} still bad: {still_bad} ===")
    (sweep.HERE / "RESUME_REFIX_DONE").write_text(
        f"{len(solved)} solved, {len(still_bad)} still bad as of "
        f"{time.strftime('%Y-%m-%dT%H:%M:%S')}\n")


if __name__ == "__main__":
    main()
