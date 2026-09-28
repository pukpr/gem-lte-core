#!/usr/bin/env python3
"""Overnight run: work through every remaining poor-holdout SST grid cell,
flood-filling from the current (continuously growing) solved set, which
now includes nino34/nao/tna/iode/iodw as additional donor nodes alongside
the original 6 flagships -- the first real coverage of the Indian Ocean
basin. Each cell is wrapped in its own try/except so one failure (a bug,
a missing file, a stuck search) can never take down the whole run.
Writes a DONE marker + summary when the queue is exhausted."""
import sys, json, time, traceback
sys.path.insert(0, "/home/paul/eval/gem-lte-core/experiments/Sep2026")
import sweep

TIMEOUT = 120
RETRIES = 1
LOG = sweep.HERE / "overnight.log"


def log(msg: str) -> None:
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    with LOG.open("a") as f:
        f.write(line + "\n")


def main() -> None:
    log(f"=== overnight run starting, flagships={sweep.FLAGSHIPS} ===")
    processed = 0
    while True:
        remaining = sweep.pick_pilot_cells(n=1000)  # everything left, nearest-first
        if not remaining:
            log("queue empty -- all poor cells attempted.")
            break

        nodes = sweep.current_solved_nodes()
        # pick_pilot_cells already sorts by distance-to-current-solved-set,
        # but the set grows as we go, so re-sort against the LATEST nodes
        # each iteration rather than trusting a stale snapshot.
        def dist(c):
            lat, lon = sweep.parse_grid_name(c)
            return min(sweep.haversine_km(lat, lon, *co) for co in nodes.values())
        cell = min(remaining, key=dist)

        lat, lon = sweep.parse_grid_name(cell)
        donor_name, donor_dist = min(
            ((n, sweep.haversine_km(lat, lon, *co)) for n, co in nodes.items()),
            key=lambda x: x[1])
        donor_dir = (sweep.FEB / donor_name if donor_name in sweep.FLAGSHIPS
                     else sweep.HERE / donor_name)

        log(f"--- [{processed+1}/{processed+len(remaining)}] {cell}  "
            f"seed={donor_name} ({donor_dist:.0f} km) ---")
        t0 = time.time()
        try:
            result = sweep.solve_cell_with_retry(cell, donor_dir, TIMEOUT, retries=RETRIES)
        except Exception as exc:
            log(f"EXCEPTION on {cell}: {exc}")
            log(traceback.format_exc())
            result = dict(cell=cell, accepted=False, reason=f"exception: {exc}")
        result["seed_from_path"] = str(donor_dir)
        result["timestamp"] = time.strftime("%Y-%m-%dT%H:%M:%S")
        result["elapsed_s"] = round(time.time() - t0, 1)
        sweep.append_ledger(result)

        if result.get("accepted"):
            log(f"SOLVED {cell}: {result['attempt_tag']} "
                f"dLOD={result['dlod']:.6f} scores={result.get('scores')}")
        else:
            log(f"UNSOLVED {cell}: {result.get('reason')}")
        processed += 1

    entries = sweep.load_sweep_ledger()
    by_cell = {}
    for e in entries:
        by_cell[e["cell"]] = e
    solved = [e for e in by_cell.values() if e.get("accepted")]
    log(f"=== overnight run complete: {len(solved)}/{len(by_cell)} cells solved "
        f"(cumulative, all attempts) ===")
    (sweep.HERE / "OVERNIGHT_DONE").write_text(
        f"{len(solved)}/{len(by_cell)} solved as of {time.strftime('%Y-%m-%dT%H:%M:%S')}\n")


if __name__ == "__main__":
    main()
