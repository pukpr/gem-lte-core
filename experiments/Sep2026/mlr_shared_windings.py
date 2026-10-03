#!/usr/bin/env python3
"""Joint multi-quad MLR test of the "common eigenvalue" hypothesis
(2026-09-28 discussion with the user): is there a small SET of winding
values (k) shared across all 89 quads, with each quad's own contribution
captured purely by per-winding amplitude+phase -- and per the user's
follow-up, are those amplitudes/phases themselves smooth (not
independently free) functions of geographic position (lat, lon)?

Design: for quad i (coords lat_i, lon_i) at time t (t' = t - quad's own
t0):
  manifold(i, t) = level_i + trend_i * t'
                   + sum_k [ A_k(lat_i,lon_i) * cos(k t') + B_k(lat_i,lon_i) * sin(k t') ]

level_i, trend_i stay per-quad free parameters (a long-term SST
warming/offset is not what's being tested here). A_k, B_k are each a
linear combination of a low-order (degree-1, 4-term) real spatial basis
-- shared GLOBALLY across all 89 quads -- rather than 89 independent
values per k. This turns the whole multi-quad, multi-winding fit into
ONE linear regression (all unknowns enter linearly), solved via a single
lstsq call.

Compared against two baselines to isolate what's doing the work:
  (a) per-quad-independent amplitude/phase for the SAME shared k-set
      (no geographic smoothness constraint) -- tests whether the
      geographic parametrization costs much fit quality.
  (b) per-quad's OWN already-established best fit (train/validate/test
      from the ledger) -- a different metric (CC/DTW against real SST,
      not R^2 against the manifold) but the closest available reference
      point for "how good is good".
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, "/home/paul/eval/gem-lte-core/experiments/Sep2026")
import sweep
import plot_hidden_latent_forcing_all_quads as m

import numpy as np

BACKBONE = 0.207
HARMONICS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 18, 19, 23, 27]
SUBHARMONIC = BACKBONE / 11
CANDIDATE_KS = [h * BACKBONE for h in HARMONICS] + [SUBHARMONIC]


def spatial_basis(lat_deg: float, lon_deg: float) -> np.ndarray:
    """Degree-1 real spherical-harmonic-like basis: [1, sin(lat),
    cos(lat)sin(lon), cos(lat)cos(lon)] -- a smooth global dipole-order
    expansion, 4 free coefficients per winding per (cos/sin) rather than
    89 independent ones."""
    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)
    return np.array([1.0, np.sin(lat), np.cos(lat) * np.sin(lon), np.cos(lat) * np.cos(lon)])


N_SPATIAL_BASIS = 4
IR_LAG = 12  # months -- matches Ada's Model(I) -= IR*Model(I-12)


def apply_ir_filter(raw: np.ndarray, ir: float) -> np.ndarray:
    """Recursive 12-month-lag IIR feedback, matching
    gem-lte-primitives-solution.adb's Model(I) -= IR*Model(I-12) exactly
    (unconditional whenever IR /= 0, per that source -- not gated by
    UNCOMPENSATED, which only affects the separate Data_Records-side
    transform). IR is per-quad and already fixed (read from that quad's
    own lt.exe.p, not re-estimated here), so this is a known LINEAR
    operator -- applying it to each raw basis column before stacking
    keeps the whole multi-quad fit a single linear regression, per the
    user's framing of IR as a secondary factor that reshapes the
    manifold rather than a new free winding parameter."""
    if ir == 0.0:
        return raw
    out = raw.copy()
    for i in range(IR_LAG, len(out)):
        out[i] = raw[i] - ir * out[i - IR_LAG]
    return out


def load_quads(root: Path | None = None):
    """Default: each quad's freshest lte_results.csv/lt.exe.p across
    Feb2026 and Sep2026 (original behaviour). With root: that folder's
    quads only, e.g. --root ../Feb2026_sweep."""
    base = root if root is not None else sweep.FEB
    cells = sorted((p.name for p in base.iterdir()
                     if p.is_dir() and sweep.parse_grid_name(p.name)
                     and (root is None or (p / "lte_results.csv").exists())))
    data = []
    for cell in cells:
        if root is not None:
            path = root / cell / "lte_results.csv"
            p_path = root / cell / "lt.exe.p"
        else:
            path = m.freshest_results_csv(cell)
            p_path = max((sweep.FEB / cell / "lt.exe.p", sweep.HERE / cell / "lt.exe.p"),
                         key=lambda p: p.stat().st_mtime if p.exists() else -1)
        times, forcing = m.read_time_and_forcing(path)
        lat, lon = sweep.parse_grid_name(cell)
        ir = json.loads(p_path.read_text()).get("IR", 0.0)
        data.append(dict(cell=cell, times=np.array(times), forcing=np.array(forcing),
                          lat=lat, lon=lon, ir=ir))
    return data


def fit_geo_smooth(data):
    """Joint regression: per-quad level+trend (free), per-winding
    amp/phase constrained to a shared degree-1 spatial basis."""
    n_quads = len(data)
    n_k = len(CANDIDATE_KS)
    n_quad_params = 2 * n_quads  # level_i, trend_i
    n_shared_params = 2 * n_k * N_SPATIAL_BASIS  # A_k,B_k each x 4 basis funcs

    rows = []
    targets = []
    row_owner = []  # which quad each row belongs to, for per-quad R^2 later

    for qi, d in enumerate(data):
        t0 = d["times"][0]
        tt = d["times"] - t0
        sbasis = spatial_basis(d["lat"], d["lon"])  # (4,)
        n_t = len(tt)

        # Build each RAW basis column over this quad's full time series,
        # then IR-filter each one (linear op, see apply_ir_filter) before
        # stacking into the design matrix -- level_i and trend_i*t' get
        # filtered too, matching the Ada Model(I) -= IR*Model(I-12) which
        # applies to the whole assembled Model, not just the winding part.
        quad_cols = np.zeros((n_quad_params + n_shared_params, n_t))
        quad_cols[2 * qi, :] = apply_ir_filter(np.ones(n_t), d["ir"])
        quad_cols[2 * qi + 1, :] = apply_ir_filter(tt.copy(), d["ir"])
        col = n_quad_params
        for k in CANDIDATE_KS:
            cos_col = apply_ir_filter(np.cos(k * tt), d["ir"])
            sin_col = apply_ir_filter(np.sin(k * tt), d["ir"])
            for j in range(N_SPATIAL_BASIS):
                quad_cols[col + j, :] = sbasis[j] * cos_col
            col += N_SPATIAL_BASIS
            for j in range(N_SPATIAL_BASIS):
                quad_cols[col + j, :] = sbasis[j] * sin_col
            col += N_SPATIAL_BASIS

        for ti in range(n_t):
            rows.append(quad_cols[:, ti])
            targets.append(d["forcing"][ti])
            row_owner.append(qi)

    X = np.array(rows)
    y = np.array(targets)
    row_owner = np.array(row_owner)
    print(f"Joint geo-smooth design matrix: {X.shape[0]} rows x {X.shape[1]} cols "
          f"({n_quad_params} per-quad + {n_shared_params} shared-spatial)")

    coeffs, *_ = np.linalg.lstsq(X, y, rcond=None)
    pred = X @ coeffs
    overall_r2 = 1 - np.sum((y - pred) ** 2) / np.sum((y - y.mean()) ** 2)
    print(f"OVERALL R^2 (geo-smooth shared basis): {overall_r2:.4f}")

    per_quad_r2 = {}
    for qi, d in enumerate(data):
        mask = row_owner == qi
        yq, pq = y[mask], pred[mask]
        r2 = 1 - np.sum((yq - pq) ** 2) / np.sum((yq - yq.mean()) ** 2)
        per_quad_r2[d["cell"]] = r2
    return overall_r2, per_quad_r2, coeffs


def fit_independent_baseline(data):
    """Same shared k-set, but each quad gets its OWN independent
    amp/phase (no geographic constraint) -- isolates how much the
    geo-smoothness parametrization costs vs. full per-quad freedom."""
    per_quad_r2 = {}
    r2s = []
    for d in data:
        t0 = d["times"][0]
        tt = d["times"] - t0
        cols = [apply_ir_filter(np.ones_like(tt), d["ir"]),
                apply_ir_filter(tt.copy(), d["ir"])]
        for k in CANDIDATE_KS:
            cols.append(apply_ir_filter(np.cos(k * tt), d["ir"]))
            cols.append(apply_ir_filter(np.sin(k * tt), d["ir"]))
        X = np.column_stack(cols)
        coeffs, *_ = np.linalg.lstsq(X, d["forcing"], rcond=None)
        pred = X @ coeffs
        r2 = 1 - np.sum((d["forcing"] - pred) ** 2) / np.sum((d["forcing"] - d["forcing"].mean()) ** 2)
        per_quad_r2[d["cell"]] = r2
        r2s.append(r2)
    print(f"Per-quad-independent baseline: mean R^2 = {np.mean(r2s):.4f}, "
          f"min = {np.min(r2s):.4f}, max = {np.max(r2s):.4f}")
    return per_quad_r2


def main() -> None:
    ap = argparse.ArgumentParser(description="Joint multi-quad shared-winding MLR (collective R^2).")
    ap.add_argument("--root", type=Path, default=None,
                    help="read every quad from this folder only (default: freshest of Feb2026/Sep2026)")
    args = ap.parse_args()
    if args.root is not None:
        print(f"Quads from: {args.root.resolve()}")
    print(f"Candidate shared k-set ({len(CANDIDATE_KS)} windings): "
          f"{[round(k,4) for k in CANDIDATE_KS]}")
    data = load_quads(args.root.resolve() if args.root is not None else None)
    print(f"Loaded {len(data)} quads\n")

    print("=== Baseline: fully independent per-quad amp/phase (same k-set) ===")
    indep_r2 = fit_independent_baseline(data)

    print("\n=== Geo-smooth: amp/phase constrained to degree-1 spatial basis ===")
    overall_r2, geo_r2, coeffs = fit_geo_smooth(data)

    print("\n=== Per-quad comparison (worst 15 geo-smooth R^2) ===")
    ranked = sorted(geo_r2.items(), key=lambda kv: kv[1])
    for cell, r2 in ranked[:15]:
        print(f"  {cell}: geo-smooth R^2={r2:.3f}  independent R^2={indep_r2[cell]:.3f}  "
              f"gap={indep_r2[cell]-r2:.3f}")

    print("\n=== Summary ===")
    geo_vals = np.array(list(geo_r2.values()))
    indep_vals = np.array([indep_r2[c] for c in geo_r2])
    print(f"Geo-smooth per-quad R^2: mean={geo_vals.mean():.3f} min={geo_vals.min():.3f} max={geo_vals.max():.3f}")
    print(f"Independent per-quad R^2: mean={indep_vals.mean():.3f} min={indep_vals.min():.3f} max={indep_vals.max():.3f}")
    print(f"Mean gap (independent - geo-smooth): {(indep_vals-geo_vals).mean():.3f}")
    n_params_geo = 2*len(data) + 2*len(CANDIDATE_KS)*N_SPATIAL_BASIS
    n_params_indep = len(data) * (2 + 2*len(CANDIDATE_KS))
    print(f"Parameter count: geo-smooth={n_params_geo}  independent={n_params_indep}  "
          f"(geo-smooth uses {100*n_params_geo/n_params_indep:.1f}% as many parameters)")


if __name__ == "__main__":
    main()
