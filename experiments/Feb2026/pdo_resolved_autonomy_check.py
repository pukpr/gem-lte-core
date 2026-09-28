#!/usr/bin/env python3
"""pdo_resolved_autonomy_check.py — does an autonomous (self-referential)
feature creep into the FINAL STAGE of the PDO resolved-fluid pipeline
(pdo_resolved_shallow_water.py's fixed/moving gauge output), specifically
a delayed-differential-style term such as a 12-month lag, or any other
lag this project's own AR(1)-nulled discipline would flag as genuine
memory rather than ordinary noise persistence?

This directly reuses the ALREADY-VALIDATED sindy_discovery.py machinery
(build_library/sparse_fit for "does a lagged-self term survive sparse
selection alongside the Forcing/winding library", fit_nonauto_only +
autonomy_significance for "is any such lag-structure more than an
AR(1)-matched null would give by chance") -- built and synthetic-
validated earlier this session specifically for this question, now
pointed at the resolved-PDE model's own gauge output instead of the
production Ada model's regression output (sindy_discovery.py real pdo
already answers the question for THAT stage; see printed comparison).

Why this is a different question for a PDE than for the ODE-only
production model: the resolved model is LITERALLY a partial differential
equation with real spatial derivatives and finite wave-transit time
(the QBO vertical model's emergent 4-7 month lag came from EXACTLY this
mechanism, with no explicit lag parameter anywhere). A PDE's own genuine
memory is still driven entirely by the external Forcing history (it is
NOT a self-sustained autonomous oscillator -- with zero forcing it just
decays, per every sanity check this session has run), but from the
OUTSIDE, at a single sampling point, that memory can look exactly like
an autonomous lag term to a fitting procedure that only sees the
gauge output and the INSTANTANEOUS forcing value, not the forcing's full
history. The moving gauge is a second, sharper way this can happen: its
sampling position depends on the CURRENT forcing value, but what it
samples there was shaped by PAST forcing through wave propagation -- so
its output is a genuinely path-dependent (memory-full) functional of the
forcing history, not an instantaneous transform of it.

Usage
-----
    ./pdo_resolved_autonomy_check.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from baltic_resolved_shallow_water import run_solver, G                    # noqa: E402
from pdo_resolved_shallow_water import (build_grid, build_forcing_at_yl,   # noqa: E402
                                        production_yl, diagonal_tilt_field,
                                        full_fcor)
from sindy_discovery import (build_library, sparse_fit, fit_nonauto_only,  # noqa: E402
                             autonomy_significance, AR_LAGS_MONTHS)

ROOT = Path(__file__).resolve().parent


def regenerate(idx: str = "pdo", years: float = 140.0, dx_km: float = 300.0,
               drho_rho: float = 0.005, depth: float = 700.0,
               r_drag: float = 1e-7, nu_scale: float = 0.02,
               tide_amp: float = 0.02, lat_south: float = 20.0,
               lat_north: float = 60.0, alpha: float = 0.5,
               k_shift: float = 0.1):
    """Re-runs exactly pdo_resolved_shallow_water.run()'s solve step
    (same defaults), returning the arrays that script only prints."""
    _, mask, open_mask, dx, nx, ny, lat = build_grid(dx_km, lat_south, lat_north)
    H = np.full((ny, nx), depth)
    f_cor = full_fcor(lat, nx)
    g_reduced = G * drho_rho
    c = np.sqrt(g_reduced * depth)
    dt = 0.4 * dx / c
    nu = nu_scale * dx ** 2 / dt

    yl = production_yl(idx)
    dates, forcing, obs = build_forcing_at_yl(idx, yl)
    t_years = dates - dates[0]

    def forcing_fn(t_sec: float) -> float:
        t_yr = t_sec / (365.2422 * 86400.0)
        return float(np.interp(t_yr, t_years, forcing))

    dt_out_target = 1.0 / 12.0 * 365.2422 * 86400.0
    save_every = max(1, round(dt_out_target / dt))
    n_steps = int(years * 365.2422 * 86400.0 / dt)
    n_steps -= n_steps % save_every

    gauge_j, gauge_i = ny // 2, nx // 2
    tilt = diagonal_tilt_field(mask, alpha)
    gauge, max_abs, _, row_hist = run_solver(
        H, mask, open_mask, f_cor, dx, g_reduced, r_drag, nu, tide_amp,
        forcing_fn, n_steps, dt, save_every, gauge_j, gauge_i, ny, nx,
        tilt_field=tilt, save_row=True)

    gauge_dates = dates[0] + np.arange(len(gauge)) / 12.0
    gauge_forcing = np.interp(gauge_dates, dates, forcing)

    x_cols = np.arange(nx)
    virtual_i = np.clip(gauge_i + k_shift * gauge_forcing, 1.0, nx - 2.0)
    moving = np.array([np.interp(vi, x_cols, row_hist[t])
                       for t, vi in enumerate(virtual_i)])

    return gauge_dates, gauge_forcing, gauge, moving


def check(label: str, forcing: np.ndarray, data: np.ndarray, m_max: float,
         dm: float, threshold: float) -> None:
    m_grid = np.arange(dm, m_max + dm / 2, dm)

    Theta, names, kinds, y = build_library(forcing, data, m_grid, AR_LAGS_MONTHS)
    result = sparse_fit(Theta, y, names, kinds, threshold)
    print(f"\n-- {label} --")
    print(f"  sparse fit (Forcing/winding library + lag-{AR_LAGS_MONTHS} "
          f"self-terms): R^2={result['r2']:.3f}, {result['n_survivors']} "
          f"survivors")
    print(f"  raw variance share (uncorrected, see nulled test below): "
          f"nonauto={result['nonauto_share']:.1%}  auto={result['auto_share']:.1%}")
    for name, kind, c in result["survivors"][:10]:
        print(f"    [{kind:7s}] {name:24s} {c:+.3f}")

    resid = fit_nonauto_only(forcing, data, m_grid, threshold)
    r2_nonauto_only = 1.0 - np.var(resid) / np.var(data)
    print(f"  non-autonomous-only fit R^2={r2_nonauto_only:.3f} "
          f"(Forcing/winding library alone, no lag terms available at all)")

    for lags, tag in ([1], "lag-1 (ordinary AR memory)"), \
                     ([12], "lag-12 ONLY (literal 12mo delayed-differential)"), \
                     ([1, 12], "lag-1 AND lag-12 jointly"):
        sig = autonomy_significance(resid, lags, n_surrogates=200, seed=0)
        verdict = ("SIGNIFICANT beyond AR(1) null" if sig["significant"]
                  else "consistent with ordinary AR(1) noise persistence")
        print(f"  [{tag}] rho={sig['rho']:.2f}  real R^2={sig['r2_real']:.3f} "
              f"vs null mean={sig['null_mean']:.3f} p95={sig['null_p95']:.3f} "
              f"(percentile {sig['percentile']:.0f}) -> {verdict}")


def main() -> int:
    print("Regenerating PDO resolved-model gauges (fixed + moving, "
          "k_shift=0.1)...")
    gauge_dates, gauge_forcing, gauge, moving = regenerate()
    print(f"  {len(gauge_dates)} monthly samples, "
          f"{gauge_dates[0]:.0f}-{gauge_dates[-1]:.0f}")

    check("FIXED gauge (resolved PDE, before any moving-gauge pullback)",
         gauge_forcing, gauge, m_max=1.0, dm=0.01, threshold=0.08)
    check("MOVING gauge (k_shift=0.1, material-coordinate pullback)",
         gauge_forcing, moving, m_max=1.0, dm=0.01, threshold=0.08)
    return 0


if __name__ == "__main__":
    sys.exit(main())
