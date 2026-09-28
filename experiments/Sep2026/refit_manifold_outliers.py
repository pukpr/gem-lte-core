#!/usr/bin/env python3
"""Refit the 6 manifold-coherence outliers found by
plot_hidden_latent_forcing_all_quads.py (kS040_W030, kS020_E150,
kN040_E150, kN000_W150, kS020_E010, kN020_E150) via donor-seeding from a
clean, geographically-proximal, non-outlier quad -- per the user's
judgment (2026-09-28) that these 6 show nothing physically distinct from
their neighbors, so the divergence is more likely a fit artifact (wrong
harmonic branch / local optimum) than real physics.

For each cell: back up its current (accepted) state, run the normal
sweep.py donor-seeded cascade, then judge success on BOTH criteria --
not just CC/DTW score (which was already fine before!) but whether the
resulting manifold's detrended correlation to the cross-quad ensemble
mean actually improves. Roll back to the original state if either the
cascade fails outright or the manifold coherence doesn't improve.
"""
import shutil
import sys
import time
from pathlib import Path

sys.path.insert(0, "/home/paul/eval/gem-lte-core/experiments/Sep2026")
import sweep
import plot_hidden_latent_forcing_all_quads as m
import numpy as np

OUTLIERS = ['kS040_W030', 'kS020_E150', 'kN040_E150', 'kN000_W150',
            'kS020_E010', 'kN020_E150']

DONORS = {
    'kS040_W030': 'kS040_W050',
    'kS020_E150': 'kS020_E170',
    'kN040_E150': 'kN040_E130',
    'kN000_W150': 'nino34',
    'kS020_E010': 'kS020_E030',
    'kN020_E150': 'kN020_E170',
}


def ensemble_mean_excluding(exclude: set[str]) -> tuple[list[float], list[float]]:
    """Detrended, unit-RMS ensemble mean manifold built from every quad
    NOT in `exclude` -- the reference to compare a refit candidate
    against, so a cell can't just correlate with its own stale copy."""
    cells = sorted((p.name for p in sweep.FEB.iterdir()
                     if p.is_dir() and m.NAME_RE.match(p.name)))
    series = []
    for cell in cells:
        if cell in exclude:
            continue
        path = m.freshest_results_csv(cell)
        if path is None:
            continue
        times, forcing = m.read_time_and_forcing(path)
        if len(times) < 2:
            continue
        series.append((times, m.detrend_and_renormalize(times, forcing)))
    start = max(t[0] for t, _ in series)
    stop = min(t[-1] for t, _ in series)
    n_months = int((stop - start) * 12) + 1
    grid_t = [start + i / 12 for i in range(n_months)]
    stack = np.full((len(series), n_months), np.nan)
    for row, (t, f) in enumerate(series):
        stack[row, :] = np.interp(grid_t, t, f, left=np.nan, right=np.nan)
    return grid_t, list(np.nanmean(stack, axis=0))


def coherence(cell: str, ref_t: list[float], ref_v: list[float]) -> float:
    path = m.freshest_results_csv(cell)
    times, forcing = m.read_time_and_forcing(path)
    dt = m.detrend_and_renormalize(times, forcing)
    return m.pearson_correlation(dt, np.interp(times, ref_t, ref_v))


def backup(cell: str) -> dict:
    d = sweep.HERE / cell
    saved = {}
    for fname in ("lt.exe.p", "lt.exe.resp", "lte_run.sh", "lte_results.csv"):
        p = d / fname
        if p.exists():
            saved[fname] = p.read_bytes()
    return saved


def restore(cell: str, saved: dict) -> None:
    d = sweep.HERE / cell
    for fname, content in saved.items():
        (d / fname).write_bytes(content)


def main() -> None:
    print("Computing baseline ensemble reference (excluding all 6 outliers)...")
    ref_t, ref_v = ensemble_mean_excluding(set(OUTLIERS))

    results = {}
    for cell in OUTLIERS:
        print(f"\n=== {cell} (donor: {DONORS[cell]}) ===")
        before_coh = coherence(cell, ref_t, ref_v)
        print(f"  coherence BEFORE refit: {before_coh:.3f}")

        saved = backup(cell)
        by_cell = sweep.latest_entry_per_cell(sweep.load_sweep_ledger())
        prev_entry = by_cell.get(cell)

        donor_name = DONORS[cell]
        donor_dir = sweep.FEB / donor_name if donor_name in sweep.FLAGSHIPS else sweep.HERE / donor_name
        result = sweep.solve_cell_with_retry(cell, donor_dir, 120, retries=1)

        if not result.get("accepted"):
            print(f"  cascade FAILED ({result.get('reason')}) -- rolling back, keeping original")
            restore(cell, saved)
            results[cell] = dict(action="rolled_back_cascade_failed", before=before_coh)
            continue

        after_coh = coherence(cell, ref_t, ref_v)
        print(f"  cascade succeeded: {result['attempt_tag']} scores={result['scores']}")
        print(f"  coherence AFTER refit: {after_coh:.3f}")

        if after_coh > before_coh + 0.05:
            result["seed_from_path"] = str(donor_dir)
            result["timestamp"] = time.strftime("%Y-%m-%dT%H:%M:%S")
            result["note"] = (
                f"Refit to fix manifold-coherence outlier status (detrended "
                f"correlation to ensemble improved {before_coh:.3f} -> "
                f"{after_coh:.3f}). Superseded previous entry: "
                f"{prev_entry.get('attempt_tag') if prev_entry else None} "
                f"scores={prev_entry.get('scores') if prev_entry else None}."
            )
            sweep.append_ledger(result)
            print(f"  COMMITTED (coherence improved by {after_coh - before_coh:.3f})")
            results[cell] = dict(action="committed", before=before_coh, after=after_coh)
        else:
            print(f"  coherence did NOT meaningfully improve -- rolling back, keeping original")
            restore(cell, saved)
            results[cell] = dict(action="rolled_back_no_improvement", before=before_coh, after=after_coh)

    print("\n=== SUMMARY ===")
    for cell, r in results.items():
        print(f"  {cell}: {r}")


if __name__ == "__main__":
    main()
