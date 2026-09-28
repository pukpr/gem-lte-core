#!/usr/bin/env python3
"""nao_baltic_check.py — applies the NAO atmospheric shallow-water model
(nao_resolved_shallow_water.py) to the BALTIC time series instead of NAO
itself, motivated by the documented real correlation between NAO and
Baltic coastal MSL (wind/inverse-barometer forcing -- the same
atmospheric circulation NAO indexes also sets up the winds and pressure
patterns that raise/lower Baltic sea level), and by baltic's own
strongest independently-established winding ridge (winding_rank.py
--index baltic: M=1.250, peak=4.21 bits, fwhm=0.040, continuity=1.00,
PASS -- see project memory [[cv-tools]]) landing almost exactly on
6x this project's shared 0.2075 cross-index backbone (0.2075*6=1.245),
which NAO's own fitted top winding (0.2077) also sits on. If the
atmospheric circulation the NAO model represents is genuinely an
intermediate step between the shared tidal manifold and Baltic's own
response, some gauge extraction from this SAME model might reveal M=1.25
directly, or fit real Baltic data better than it fit real NAO.

Grounded first: real NAO vs real Baltic MSL (1880-2023 overlap, monthly),
zero-lag r=-0.290, and this IS the best-aligned lag (no lag out to +/-24
months improves it) -- a real, immediate, moderate relationship, matching
the literature's own framing (Baltic responds to the CURRENT wind/
pressure pattern, not a delayed one).

Usage
-----
    ./nao_baltic_check.py
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from baltic_resolved_shallow_water import run_solver, G                    # noqa: E402
from amo_resolved_shallow_water import full_fcor                           # noqa: E402
from nao_resolved_shallow_water import build_grid                          # noqa: E402
from pdo_resolved_shallow_water import (build_forcing_at_yl, production_yl,  # noqa: E402
                                        diagonal_tilt_field)
from mode_evidence import sweep_M                                          # noqa: E402
from signal_operators import delayed_difference, dipole                    # noqa: E402
from winding_rank import rank_series                                       # noqa: E402
from qbo_compensating_dtw import standardize, dtw_distance, null_baseline  # noqa: E402

ROOT = Path(__file__).resolve().parent


def regenerate_all(idx: str = "nao", years: float = 140.0, dx_km: float = 300.0,
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
    j_azores = int(np.argmin(np.abs(lat - 38.0)))
    i0 = nx // 2
    tilt = diagonal_tilt_field(mask, alpha)

    t0 = time.time()
    gauge_iceland, max_abs, _, row_iceland = run_solver(
        H, mask, open_mask, f_cor, dx, g_eq, r_drag, nu, tide_amp,
        forcing_fn, n_steps, dt, save_every, j_iceland, i0, ny, nx,
        tilt_field=tilt, save_row=True)
    gauge_azores, _, _, _ = run_solver(
        H, mask, open_mask, f_cor, dx, g_eq, r_drag, nu, tide_amp,
        forcing_fn, n_steps, dt, save_every, j_azores, i0, ny, nx,
        tilt_field=tilt, save_row=False)
    print(f"  regenerated in {time.time()-t0:.1f}s, max|zeta|={max_abs:.4f}")

    gauge_dates = dates[0] + np.arange(len(gauge_iceland)) / 12.0
    gauge_forcing = np.interp(gauge_dates, dates, forcing)

    x_cols = np.arange(nx)
    k_shift = 0.1
    virtual_i = np.clip(i0 + k_shift * gauge_forcing, 1.0, nx - 2.0)
    moving = np.array([np.interp(vi, x_cols, row_iceland[t])
                       for t, vi in enumerate(virtual_i)])
    dipole_series = dipole(gauge_azores, gauge_iceland)

    return gauge_dates, gauge_forcing, gauge_iceland, moving, dipole_series


def smooth(x: np.ndarray, win_months: int) -> np.ndarray:
    k = np.ones(win_months) / win_months
    return np.convolve(x, k, mode="same")


def corr_dtw(name, series, obs_ref, seed_offset=0):
    r_raw = float(np.corrcoef(series, obs_ref)[0, 1])
    sign = -1.0 if r_raw < 0 else 1.0
    print(f"\n-- {name} vs real Baltic --")
    print(f"  Pearson r (sign-matched): {sign * r_raw:+.3f}")

    model_s = standardize(series)
    obs_s = standardize(obs_ref)
    model_signed = sign * model_s
    d_true, _ = dtw_distance(model_signed, obs_s)
    rng = np.random.default_rng(seed_offset)
    null_mean, null_std = null_baseline(rng, model_signed, obs_s, n_trials=16)
    z_raw = (null_mean - d_true) / null_std if null_std > 0 else float("nan")

    win = 3 * 12
    m3 = standardize(smooth(series, win))
    o3 = standardize(smooth(obs_ref, win))
    r3 = float(np.corrcoef(m3, o3)[0, 1])
    sign3 = -1.0 if r3 < 0 else 1.0
    d3, _ = dtw_distance(sign3 * m3, o3)
    rng2 = np.random.default_rng(seed_offset + 1)
    null_mean3, null_std3 = null_baseline(rng2, sign3 * m3, o3, n_trials=16)
    z3 = (null_mean3 - d3) / null_std3 if null_std3 > 0 else float("nan")

    print(f"  DTW raw-monthly z-score: {z_raw:+.2f}   "
          f"DTW 3yr-smoothed z-score: {z3:+.2f} (r_3yr={sign3*r3:+.3f})")


def ridge_check(name, series, dates_local, forcing_local, forcing_label,
                m_max, dm):
    ranked = rank_series(dates_local, series, forcing_local, m_max=m_max,
                         dm=dm, sigma=15.0, t0_step=5.0)
    near_125 = [r for r in ranked if abs(r["M"] - 1.25) < 0.03]
    print(f"  [{name}] winding_rank vs {forcing_label} (m_max={m_max}): "
          f"{len(ranked)} candidate(s)")
    for r in ranked[:5]:
        flag = "  <-- near baltic's real M=1.25" if abs(r["M"] - 1.25) < 0.03 else ""
        print(f"    M={r['M']:.3f}  peak={r['peak_bits']:.2f}bits  "
              f"fwhm={r['fwhm_M']:.3f}  cont={r['continuity']:.2f}  "
              f"pass={r['pass']}{flag}")
    if not ranked:
        print("    (no candidates reached the peak-power screen)")
    elif not near_125:
        print("    (no candidate landed near M=1.25)")


def main() -> int:
    cache = ROOT / "nao" / "nao_baltic_check_cache.npz"
    if cache.exists():
        print(f"Loading cached NAO atmospheric model output from {cache}...")
        d = np.load(cache)
        gauge_dates, gauge_forcing, fixed, moving, dipole_series = (
            d["gauge_dates"], d["gauge_forcing"], d["fixed"], d["moving"],
            d["dipole_series"])
    else:
        print("Regenerating NAO atmospheric model (fixed/moving/dipole "
              "gauges, driven by nao's own manifold)...")
        gauge_dates, gauge_forcing, fixed, moving, dipole_series = regenerate_all()
        np.savez(cache, gauge_dates=gauge_dates, gauge_forcing=gauge_forcing,
                fixed=fixed, moving=moving, dipole_series=dipole_series)
        print(f"  cached to {cache} for reuse")

    diff = delayed_difference(fixed, 12)
    diff_dates = gauge_dates[12:]
    diff_forcing = gauge_forcing[12:]

    baltic_dat = np.loadtxt(ROOT / "baltic" / "baltic.dat")
    dates_baltic, obs_baltic = baltic_dat[:, 0].copy(), baltic_dat[:, 1].copy()
    obs_on_gauge = np.interp(gauge_dates, dates_baltic, obs_baltic)
    obs_on_diff = np.interp(diff_dates, dates_baltic, obs_baltic)

    # BALTIC's own real production manifold -- the coordinate its own
    # M=1.25 ridge was established in (winding_rank --index baltic) --
    # the fairer test of "does this model's output carry baltic's own
    # winding signature" than scanning against NAO's driving forcing.
    yl_baltic = production_yl("baltic")
    dates_b, forcing_baltic_full, _ = build_forcing_at_yl("baltic", yl_baltic)
    forcing_baltic_on_gauge = np.interp(gauge_dates, dates_b, forcing_baltic_full)
    forcing_baltic_on_diff = np.interp(diff_dates, dates_b, forcing_baltic_full)

    variants = [
        ("FIXED single point (Iceland)", fixed, gauge_dates, obs_on_gauge,
         gauge_forcing, forcing_baltic_on_gauge),
        ("MOVING single point (k_shift=0.1)", moving, gauge_dates, obs_on_gauge,
         gauge_forcing, forcing_baltic_on_gauge),
        ("SPATIAL DIPOLE (Azores-Iceland)", dipole_series, gauge_dates,
         obs_on_gauge, gauge_forcing, forcing_baltic_on_gauge),
        ("FIXED, 12mo delayed-difference", diff, diff_dates, obs_on_diff,
         diff_forcing, forcing_baltic_on_diff),
    ]

    print("\n================ correlation + DTW vs real Baltic MSL "
          "================")
    for i, (name, series, dts, obs, _, _) in enumerate(variants):
        corr_dtw(name, series, obs, seed_offset=2 * i)

    print("\n================ ridge scan vs NAO's OWN driving Forcing "
          "================")
    for name, series, dts, _, forc_nao, _ in variants:
        ridge_check(name, series, dts, forc_nao, "NAO's own Forcing",
                   m_max=1.5, dm=0.005)

    print("\n============= ridge scan vs BALTIC's OWN production Forcing "
          "(the fairer M=1.25 test) =============")
    for name, series, dts, _, _, forc_balt in variants:
        ridge_check(name, series, dts, forc_balt, "BALTIC's own Forcing",
                   m_max=1.5, dm=0.005)

    print("\n(reference: real NAO vs real Baltic MSL, zero-lag, best "
          "alignment: r=-0.290; baltic's own independently-fitted "
          "winding_rank ridge: M=1.250, peak=4.21 bits, fwhm=0.040, "
          "continuity=1.00, PASS)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
