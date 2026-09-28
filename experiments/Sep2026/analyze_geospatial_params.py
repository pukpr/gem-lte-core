#!/usr/bin/env python3
"""Statistical (non-plotting) check for geospatial dependence in the
fitted PARAMETERS themselves -- tidal-constituent (LPAP) amplitudes,
initial condition ("init"), and every other scalar in lt.exe.p -- across
all 89 k_sst quads. Complements plot_winding_regions.py (which handles
the derived k_amp_phase winding spectrum and is being kept as the
canonical geospatial clustering per the user, 2026-09-28); this script
covers what's left: the raw fit parameters, with no plot, just numbers.

Method: for each scalar variable (17 core scalars + 42 LPAP amplitudes =
59 tests), fit a linear model against geography ([lat, sin(lon)*90,
cos(lon)*90], degrees->comparable-scale, matching plot_winding_regions.py's
convention) and get R^2. Significance is judged by PERMUTATION test (200
shuffles of the geography labels, refit each time) rather than a
parametric p-value, since 59 tests demand real multiple-comparisons
discipline -- Bonferroni-correct alpha=0.05/59 approx 0.00085 (permutation
test resolution capped at 1/201 approx 0.005, so anything at the
permutation floor is flagged as "possibly significant, needs more
permutations to confirm" rather than claimed outright).
"""
import json
import sys

sys.path.insert(0, "/home/paul/eval/gem-lte-core/experiments/Sep2026")
import sweep

import numpy as np

SCALAR_FIELDS = ["offs", "bg", "impA", "impB", "impC", "delA", "delB", "asym",
                  "ann1", "ann2", "sem1", "sem2", "year", "IR", "ma", "mp",
                  "shfT", "init"]

N_PERM = 5000
BONFERRONI_ALPHA = 0.05


def freshest_p(cell: str) -> dict:
    feb = sweep.FEB / cell / "lt.exe.p"
    sep = sweep.HERE / cell / "lt.exe.p"
    candidates = [p for p in (feb, sep) if p.exists()]
    path = max(candidates, key=lambda p: p.stat().st_mtime)
    return json.loads(path.read_text())


def geo_design(coords: np.ndarray) -> np.ndarray:
    lat = coords[:, 0]
    lon_rad = np.radians(coords[:, 1])
    return np.column_stack([np.ones(len(lat)), lat,
                             np.sin(lon_rad) * 90, np.cos(lon_rad) * 90])


def r_squared(y: np.ndarray, X: np.ndarray) -> float:
    coeffs, *_ = np.linalg.lstsq(X, y, rcond=None)
    pred = X @ coeffs
    ss_res = np.sum((y - pred) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1 - ss_res / ss_tot if ss_tot > 0 else 0.0


def permutation_p(y: np.ndarray, X_geo_cols: np.ndarray, ones: np.ndarray,
                   observed_r2: float, rng: np.random.Generator) -> float:
    count = 0
    n = len(y)
    for _ in range(N_PERM):
        perm = rng.permutation(n)
        X_shuffled = np.column_stack([ones, X_geo_cols[perm]])
        r2 = r_squared(y, X_shuffled)
        if r2 >= observed_r2:
            count += 1
    return (count + 1) / (N_PERM + 1)


def main() -> None:
    cells = sorted((p.name for p in sweep.FEB.iterdir()
                     if p.is_dir() and sweep.parse_grid_name(p.name)))
    coords = np.array([sweep.parse_grid_name(c) for c in cells])
    X = geo_design(coords)
    ones = X[:, :1]
    geo_cols = X[:, 1:]

    n_periods = None
    all_tests = []
    rng = np.random.default_rng(0)

    for field in SCALAR_FIELDS:
        y = np.array([freshest_p(c)[field] for c in cells])
        r2 = r_squared(y, X)
        p = permutation_p(y, geo_cols, ones, r2, rng)
        all_tests.append((f"scalar:{field}", r2, p))

    # LPAP amplitudes (index-wise, since periods are confirmed identical
    # across all 89 quads)
    first_p = freshest_p(cells[0])
    n_periods = len(first_p["lpap"])
    periods = [t[0] for t in first_p["lpap"]]
    for i in range(n_periods):
        y = np.array([freshest_p(c)["lpap"][i][1] for c in cells])  # amplitude
        r2 = r_squared(y, X)
        p = permutation_p(y, geo_cols, ones, r2, rng)
        all_tests.append((f"lpap_amp[period={periods[i]:.3f}d]", r2, p))

    bonf_thresh = BONFERRONI_ALPHA / len(all_tests)
    print(f"Total tests: {len(all_tests)}  Bonferroni-corrected alpha: {bonf_thresh:.5f}")
    print(f"(permutation test floor: {1/(N_PERM+1):.4f} -- p-values at that floor "
          f"need more permutations to resolve further, flagged separately)\n")

    all_tests.sort(key=lambda t: t[2])
    print(f"{'variable':<32} {'R^2':>8} {'perm p':>10}  verdict")
    print("-" * 70)
    for name, r2, p in all_tests:
        if p < bonf_thresh:
            verdict = "SIGNIFICANT (survives Bonferroni)"
        elif p <= 1 / (N_PERM + 1):
            verdict = "at permutation floor -- needs more perms to confirm"
        elif p < 0.05:
            verdict = "nominal p<0.05, does NOT survive correction"
        else:
            verdict = "not significant"
        print(f"{name:<32} {r2:>8.3f} {p:>10.4f}  {verdict}")


if __name__ == "__main__":
    main()
