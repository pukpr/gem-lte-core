#!/usr/bin/env python3
"""nao_winding_rank_check.py — does the 12-month delayed-difference that
suppresses NAO's 60/120yr band (nao_resolved_shallow_water.py) also
degrade ridge COHERENCE, not just ridge power? Runs winding_rank.py's
generic rank_series() (AR(1)-floor-nulled peak/sharpness/continuity
ridge test, the same discipline as winding_scalogram.py) on the NAO
resolved model's FIXED (Iceland) gauge, both raw and 12mo-differenced,
against the SAME production-manifold-derived Forcing in both cases.

Usage
-----
    ./nao_winding_rank_check.py
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from baltic_resolved_shallow_water import run_solver, G                    # noqa: E402
from amo_resolved_shallow_water import full_fcor                           # noqa: E402
from nao_resolved_shallow_water import build_grid                         # noqa: E402
from pdo_resolved_shallow_water import (build_forcing_at_yl, production_yl,  # noqa: E402
                                        diagonal_tilt_field)
from signal_operators import delayed_difference                            # noqa: E402
from winding_rank import rank_series                                       # noqa: E402

ROOT = Path(__file__).resolve().parent


def regenerate_fixed(idx: str = "nao", years: float = 140.0, dx_km: float = 300.0,
                     g_eq: float = 9.81, depth: float = 300.0,
                     r_drag: float = 1e-7, nu_scale: float = 0.02,
                     tide_amp: float = 0.02, lat_south: float = 25.0,
                     lat_north: float = 70.0, alpha: float = 0.5):
    _, mask, open_mask, dx, nx, ny, lat = build_grid(dx_km, lat_south, lat_north)
    H = np.full((ny, nx), depth)
    f_cor = full_fcor(lat, nx)
    c = np.sqrt(g_eq * depth)
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

    j_iceland = int(np.argmin(np.abs(lat - 65.0)))
    i0 = nx // 2
    tilt = diagonal_tilt_field(mask, alpha)

    t0 = time.time()
    gauge, max_abs, _, _ = run_solver(
        H, mask, open_mask, f_cor, dx, g_eq, r_drag, nu, tide_amp,
        forcing_fn, n_steps, dt, save_every, j_iceland, i0, ny, nx,
        tilt_field=tilt, save_row=False)
    print(f"  regenerated in {time.time()-t0:.1f}s, max|zeta|={max_abs:.4f}")

    gauge_dates = dates[0] + np.arange(len(gauge)) / 12.0
    gauge_forcing = np.interp(gauge_dates, dates, forcing)
    return gauge_dates, gauge_forcing, gauge


def main() -> int:
    print("Regenerating NAO resolved-model FIXED (Iceland) gauge...")
    t, forcing, gauge = regenerate_fixed()

    print("\n-- winding_rank on RAW fixed gauge vs Forcing --")
    ranked_raw = rank_series(t, gauge, forcing, m_max=1.0, dm=0.005, sigma=15.0,
                             t0_step=5.0)
    for r in ranked_raw:
        print(f"  M={r['M']:.3f}  peak={r['peak_bits']:.2f} bits  "
              f"fwhm={r['fwhm_M']:.3f}  continuity={r['continuity']:.2f}  "
              f"PASS={r['pass']}")
    if not ranked_raw:
        print("  (no candidates even reached the peak-power screen)")

    lag = 12
    diff = delayed_difference(gauge, lag)
    t_diff = t[lag:]
    forcing_diff = forcing[lag:]

    print("\n-- winding_rank on 12mo-DELAYED-DIFFERENCE of the fixed gauge "
          "vs the SAME (trimmed) Forcing --")
    ranked_diff = rank_series(t_diff, diff, forcing_diff, m_max=1.0, dm=0.005,
                              sigma=15.0, t0_step=5.0)
    for r in ranked_diff:
        print(f"  M={r['M']:.3f}  peak={r['peak_bits']:.2f} bits  "
              f"fwhm={r['fwhm_M']:.3f}  continuity={r['continuity']:.2f}  "
              f"PASS={r['pass']}")
    if not ranked_diff:
        print("  (no candidates even reached the peak-power screen)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
