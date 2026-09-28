#!/usr/bin/env python3
"""Step 1: invalidate every accepted cell whose backbone winding is a
spurious high-frequency alias (see sweep.backbone_ok) -- confirmed this
session by comparing the user's manually-tuned kS020_E050 "gold standard"
(backbone 0.207, matching the shared project-wide tidal backbone) against
its automated neighbors (backbone 18.68), which passed CC/dLOD gates by
coincidence despite being physically nonsensical, then propagated via
flood-fill donor inheritance.

Step 2: re-fit every invalidated cell using ONLY the confirmed-good
donor pool, nearest-to-current-good-solved-set first (recomputed each
iteration so fixed cells become valid donors for their own neighbors,
same flood-fill logic as run_overnight.py), now with the new
backbone_ok gate active so a repeat drift gets rejected automatically
instead of silently accepted again. Falls back to triangulated (k=3)
seeding if a single nearest donor doesn't resolve a cell."""
import sys, json, time, traceback
sys.path.insert(0, "/home/paul/eval/gem-lte-core/experiments/Sep2026")
import sweep

TIMEOUT = 120
RETRIES = 1
LOG = sweep.HERE / "refix.log"


def log(msg: str) -> None:
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)


def find_bad_cells() -> list[str]:
    by_cell = sweep.latest_entry_per_cell(sweep.load_sweep_ledger())
    bad = []
    for cell, e in by_cell.items():
        if not e.get("accepted"):
            continue
        cell_dir = sweep.HERE / cell
        if cell_dir.exists() and not sweep.backbone_ok(cell_dir):
            bad.append(cell)
    return bad


def main() -> None:
    bad_cells = find_bad_cells()
    log(f"=== invalidating {len(bad_cells)} bad-backbone cells: {bad_cells} ===")
    for cell in bad_cells:
        sweep.append_ledger(dict(
            cell=cell, accepted=False, reason="invalidated-bad-backbone",
            timestamp=time.strftime("%Y-%m-%dT%H:%M:%S"),
            seed_from_path="invalidation",
        ))

    remaining = set(bad_cells)
    processed = 0
    while remaining:
        nodes = sweep.current_solved_nodes()  # now correctly excludes invalidated cells

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

        log(f"--- [{processed+1}/{len(bad_cells)}] {cell}  "
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
    still_bad = find_bad_cells()
    log(f"=== refix complete: {len(solved)} total solved, "
        f"{len(still_bad)} still bad-backbone: {still_bad} ===")
    (sweep.HERE / "REFIX_DONE").write_text(
        f"{len(solved)} solved, {len(still_bad)} still bad as of "
        f"{time.strftime('%Y-%m-%dT%H:%M:%S')}\n")


if __name__ == "__main__":
    main()
