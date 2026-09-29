#!/usr/bin/env python3
"""Full 89-quad ImpC=0.01 exploration sweep (2026-09-29, per user request:
"consider the ImpC term in the next complete sweep, and report the
improvement"). Validated first on a 3-cell pilot (pilot_impc.py) -- all
three improved (kN040_W010 -0.200->0.817 [already-unstuck impC, not
directly attributable], kN060_W050 0.235->0.666 with impC growing 6x
from its 0.01 seed, kN000_E110 0.718->0.796 with impC barely moving) --
strong, consistent evidence the seeding is worthwhile project-wide.

For every one of the 89 grid cells: seed impC=0.01 via seed_impc.py's
logic (skipped if already nonzero -- a handful of cells, like
kN040_W010, already escaped 0.0 organically at some point in this
project's history), back up, re-run the normal DTW-first cascade in
place, and commit ONLY if the new result's min score is >= the old one
(this is an exploratory hypothesis test across an already-good dataset,
not a known-corruption fix like resweep_corrupted.py -- cells where
ImpC doesn't help should just keep their existing good state, not be
overwritten with a same-or-worse one).

Tracks impC before/after for every committed cell, to distinguish real
signal (moved meaningfully from the seed) from noise (stayed near 0.01
or snapped back toward 0).
"""
import json
import re
import sys
import time

sys.path.insert(0, "/home/paul/eval/gem-lte-core/experiments/Sep2026")
import sweep

SEED_VALUE = 0.01
IMPC_RE = re.compile(r'("impC"\s*:\s*)([^,\n]+)(,?)')
# Short runs, per the user's direction (2026-09-29): now that the stale
# per-CLIMATE_INDEX checkpoint bug is fixed at the Ada source, progress
# genuinely persists run-to-run, so short bursts compound instead of
# risking a regression to a worse state -- no need for a long single
# search per cell.
TIMEOUT = 90


def seed_impc(cell_dir) -> tuple[float, bool]:
    """Returns (impC value now in the file, whether this call changed it)."""
    p_path = cell_dir / "lt.exe.p"
    text = p_path.read_text()
    params = json.loads(text)
    current = params.get("impC", 0.0)
    if current != 0.0:
        return current, False
    new_text, n = IMPC_RE.subn(rf"\g<1>{SEED_VALUE}\g<3>", text, count=1)
    if n != 1:
        raise RuntimeError(f"could not locate impC line in {p_path}")
    p_path.write_text(new_text)
    return SEED_VALUE, True


def min_score(entry) -> float:
    scores = (entry or {}).get("scores") or {}
    if "train" in scores:
        return min(scores["train"], scores["validate"], scores["test"])
    if "pair" in scores:
        return min(scores["pair"])
    return -999.0


def backup(cell_dir) -> dict:
    saved = {}
    for fname in ("lt.exe.p", "lt.exe.resp", "lte_run.sh", "lte_results.csv",
                  "lt.exe.windings.json"):
        p = cell_dir / fname
        if p.exists():
            saved[fname] = p.read_bytes()
    return saved


def restore(cell_dir, saved: dict) -> None:
    for fname, content in saved.items():
        (cell_dir / fname).write_bytes(content)


def main() -> None:
    # Only target cells whose impC is STILL exactly 0.0 -- everything else
    # has already been seeded/explored (the 3 pilot cells, plus whatever
    # escaped 0.0 organically, plus the Feb2026-vs-Sep2026 import pass
    # incidentally seeding a bunch more). Computed dynamically rather than
    # a hardcoded exclusion list, since that set has grown since this
    # script was first written.
    all_cells = sorted(p.name for p in sweep.HERE.iterdir()
                        if p.is_dir() and sweep.parse_grid_name(p.name))
    names = []
    for c in all_cells:
        p_path = sweep.HERE / c / "lt.exe.p"
        if not p_path.is_file():
            continue
        if json.loads(p_path.read_text()).get("impC", 0.0) == 0.0:
            names.append(c)
    print(f"{len(names)} cells still have impC==0.0, targeting those", flush=True)
    by_cell = sweep.latest_entry_per_cell(sweep.load_sweep_ledger())
    results = {}

    for i, cell in enumerate(names):
        print(f"\n=== [{i+1}/{len(names)}] {cell} ===", flush=True)
        cell_dir = sweep.HERE / cell
        old_entry = by_cell.get(cell)
        old_min = min_score(old_entry)

        impc_before, was_seeded = seed_impc(cell_dir)
        print(f"  old min score: {old_min:.3f}  impC: {impc_before}"
              f"{' (freshly seeded)' if was_seeded else ' (already nonzero)'}",
              flush=True)

        saved = backup(cell_dir)
        t0 = time.time()
        try:
            result = sweep._run_cascade(cell, cell_dir, TIMEOUT)
        except Exception as exc:
            print(f"  EXCEPTION: {exc}", flush=True)
            restore(cell_dir, saved)
            results[cell] = dict(action="exception", old_min=old_min)
            continue
        elapsed = time.time() - t0

        if not result.get("accepted"):
            print(f"  cascade FAILED ({result.get('reason')}) after "
                  f"{elapsed:.0f}s -- rolling back", flush=True)
            restore(cell_dir, saved)
            results[cell] = dict(action="rolled_back_cascade_failed",
                                  old_min=old_min, impc_before=impc_before)
            continue

        new_min = min_score(result)
        impc_after = json.loads((cell_dir / "lt.exe.p").read_text()).get("impC")
        print(f"  cascade succeeded after {elapsed:.0f}s: {result['attempt_tag']} "
              f"scores={result['scores']}", flush=True)
        print(f"  new min score: {new_min:.3f}  (old: {old_min:.3f})  "
              f"impC after: {impc_after}", flush=True)

        if new_min >= old_min:
            result["timestamp"] = time.strftime("%Y-%m-%dT%H:%M:%S")
            result["note"] = (
                f"ImpC=0.01 exploration sweep (2026-09-29): min score "
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
            restore(cell_dir, saved)
            results[cell] = dict(action="rolled_back_no_improvement",
                                  old_min=old_min, new_min=new_min,
                                  impc_before=impc_before,
                                  impc_after=impc_after)

    print("\n=== SUMMARY ===", flush=True)
    improved = [c for c, r in results.items() if r["action"] == "committed"]
    print(f"  {len(improved)}/{len(names)} cells improved and committed",
          flush=True)
    for cell, r in results.items():
        print(f"  {cell}: {r}", flush=True)
    (sweep.HERE / "SWEEP_IMPC_DONE").write_text(
        f"{time.strftime('%Y-%m-%dT%H:%M:%S')}\n" +
        "\n".join(f"{c}: {r}" for c, r in results.items()))


if __name__ == "__main__":
    main()
