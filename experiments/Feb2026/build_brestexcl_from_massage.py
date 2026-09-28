#!/usr/bin/env python3
"""build_brestexcl_from_massage.py — replaces brestexcl/brestexcl.dat with
the massage.py-processed reconstruction of the full PSMSL Brest record
(psmsl_brest_1.rlrdata, 1807-2025), per direct instruction.

Gap rule (exact, per direct instruction -- rebuilt from scratch): fill
ONLY single isolated -99999 months by linear interpolation from their
immediate neighbors. ANY run of 2 or more consecutive missing months is
left alone entirely -- not interpolated, dropped from the working
sequence -- regardless of length (this replaces an earlier, wrong
length-THRESHOLD heuristic that still interpolated some multi-month
runs).

Two earlier bugs, found by directly inspecting the output rather than
trusting the algorithm description, still apply and are avoided here:

1. Interpolating a multi-month gap creates an artificially smooth,
   noise-free ramp; if the 13-month boxcar window sits mostly/entirely
   inside it, the output reflects the ramp's own local endpoint level
   instead of real data's natural noise (caught directly: a spurious
   +175 spike at the 1915 gap, 12 months, confirmed absent from the raw
   data itself).
2. NaN-padding excluded gaps into a FIXED monthly grid is worse for
   gaps shorter than the 13-month window: no window is ever fully NaN,
   so every gap month still gets an output, just computed from 1-2 real
   points instead of 13 -- produced an even larger spurious spike
   (+314.8) right at the edge of the same 1915 gap.

Fix for both: excluded (>=2-month) gaps are DROPPED from the working
sequence entirely (never inserted as NaN placeholders), so every boxcar
window always draws exactly window_size real/interpolated
consecutive-in-SEQUENCE samples -- exactly how the original 1943-1953
gap was always handled, generalized to every multi-month gap.

Processing on the resulting compressed sequence: 13-point boxcar smooth,
then a single GLOBAL quadratic (order=2) detrend via polyfit -- matching
~/github/pukpr/GEM-LTE/experiments/Feb2026/massage.py (invoked there as
`python3 massage.py <rlrdata> 2`).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
IDX_DIR = ROOT / "brestexcl"
PSMSL_FILE = IDX_DIR / "psmsl_brest_1.rlrdata"
OUT_FILE = IDX_DIR / "brestexcl.dat"


def find_missing_runs(missing: np.ndarray) -> list[tuple[int, int]]:
    runs = []
    i = 0
    n = len(missing)
    while i < n:
        if missing[i]:
            j = i
            while j < n and missing[j]:
                j += 1
            runs.append((i, j))
            i = j
        else:
            i += 1
    return runs


def boxcar_filter(values: np.ndarray, window_size: int = 13) -> np.ndarray:
    """Verbatim massage.py boxcar_filter -- run on a plain, gap-free
    (compressed) sequence."""
    n = len(values)
    filtered = np.zeros(n)
    half = window_size // 2
    for i in range(n):
        start, end = max(0, i - half), min(n, i + half + 1)
        filtered[i] = np.mean(values[start:end])
    return filtered


def main() -> int:
    dates_full, raw_full = [], []
    with open(PSMSL_FILE) as f:
        for line in f:
            parts = line.strip().split(";")
            dates_full.append(float(parts[0]))
            raw_full.append(float(parts[1]))
    raw_full = np.array(raw_full)
    n = len(raw_full)
    proj_dates_full = 1807.0 + np.arange(n) / 12.0  # validated row-index convention
    missing = raw_full <= -99998.0

    runs = find_missing_runs(missing)
    single_gaps = [r for r in runs if (r[1] - r[0]) == 1]
    multi_gaps = [r for r in runs if (r[1] - r[0]) > 1]
    print(f"[gaps] {len(runs)} missing runs found: {len(single_gaps)} single "
          f"isolated months (interpolated) and {len(multi_gaps)} multi-month "
          f"runs (dropped from the sequence entirely):")
    for s, e in multi_gaps:
        print(f"    {proj_dates_full[s]:.3f}-{proj_dates_full[e-1]:.3f} ({e-s} months)")

    # --- interpolate ONLY single isolated missing months -----------------
    filled = raw_full.copy()
    valid_idx = np.nonzero(~missing)[0]
    for (s, e) in single_gaps:
        before = valid_idx[valid_idx < s]
        after = valid_idx[valid_idx >= e]
        if len(before) == 0 or len(after) == 0:
            continue
        i0, i1 = before[-1], after[0]
        filled[s:e] = np.interp(np.arange(s, e), [i0, i1], [raw_full[i0], raw_full[i1]])
    print(f"[interpolate] filled {len(single_gaps)} single isolated months; "
          f"dropping {sum(e - s for s, e in multi_gaps)} months across "
          f"{len(multi_gaps)} multi-month runs")

    # --- drop every multi-month run from the sequence entirely ----------
    drop_mask = np.zeros(n, dtype=bool)
    for s, e in multi_gaps:
        drop_mask[s:e] = True
    keep = ~drop_mask
    seq_dates = proj_dates_full[keep]
    seq_values = filled[keep]

    # --- exact massage.py algorithm on the compressed sequence -----------
    filtered = boxcar_filter(seq_values, window_size=13)
    p = np.polyfit(seq_dates, filtered, 2)
    result = filtered - np.polyval(p, seq_dates)

    with open(OUT_FILE, "w") as f:
        for t, v in zip(seq_dates, result):
            f.write(f"{t:.6f}\t{v:.8f}\n")

    print(f"[written] {OUT_FILE}: {len(seq_dates)} rows, "
          f"{seq_dates[0]:.4f} - {seq_dates[-1]:.4f}")
    print(f"[stats] mean={result.mean():+.3f}  std={result.std():.3f}  "
          f"min={result.min():.3f}  max={result.max():.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
