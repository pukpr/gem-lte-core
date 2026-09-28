#!/usr/bin/env python3
"""baltic_atmospheric_check.py — drives the SAME atmospheric shallow-
water geometry built for NAO (nao_resolved_shallow_water.py: equivalent
depth ~300m, c~54m/s, full Coriolis, 25N-70N North Atlantic-European
sector) with BALTIC's OWN real production manifold instead of NAO's --
the direct next step after nao_baltic_check.py found that NAO's manifold,
run through this basin, does not reproduce baltic's own M=1.25 ridge or
outperform the real (and sign-meaningful, inverse-barometer-consistent)
NAO-Baltic correlation.

The question this asks is different and complementary: baltic's real
sea-level variability is documented to be substantially ATMOSPHERICALLY
forced (wind stress + inverse-barometer loading), not purely a slow
oceanic gravity-wave process -- so does treating BALTIC's OWN tidal
manifold through this FAST, atmospheric-equivalent-depth basin (rather
than baltic_resolved_shallow_water.py's slow OCEAN reduced-gravity
basin) produce a genuinely different, perhaps sharper, winding structure
against real Baltic data?

Gauges: FIXED and MOVING at a real Baltic-relevant latitude (58N, near
the Stockholm long tide-gauge record), plus the same Azores/Iceland
SPATIAL DIPOLE as the NAO run (does an NAO-like dipole structure emerge
from BALTIC's own manifold on its own?), plus the 12mo delayed-
difference of the Baltic-latitude fixed gauge.

Usage
-----
    ./baltic_atmospheric_check.py
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

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


def smooth(x: np.ndarray, win_months: int) -> np.ndarray:
    k = np.ones(win_months) / win_months
    return np.convolve(x, k, mode="same")


def regenerate(idx: str = "baltic", years: float = 140.0, dx_km: float = 300.0,
              g_eq: float = 9.81, depth: float = 300.0, r_drag: float = 1e-7,
              nu_scale: float = 0.02, tide_amp: float = 0.02,
              lat_south: float = 25.0, lat_north: float = 70.0,
              alpha: float = 0.5, baltic_lat: float = 58.0):
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

    j_baltic = int(np.argmin(np.abs(lat - baltic_lat)))
    j_iceland = int(np.argmin(np.abs(lat - 65.0)))
    j_azores = int(np.argmin(np.abs(lat - 38.0)))
    i0 = nx // 2
    tilt = diagonal_tilt_field(mask, alpha)

    print(f"  driving with {idx}'s own manifold (yl={yl:.4f}d); gauge rows: "
          f"baltic-lat j={j_baltic} ({lat[j_baltic]:.1f}N), "
          f"iceland j={j_iceland} ({lat[j_iceland]:.1f}N), "
          f"azores j={j_azores} ({lat[j_azores]:.1f}N)")

    t0 = time.time()
    gauge_baltic, max_abs, _, row_baltic = run_solver(
        H, mask, open_mask, f_cor, dx, g_eq, r_drag, nu, tide_amp,
        forcing_fn, n_steps, dt, save_every, j_baltic, i0, ny, nx,
        tilt_field=tilt, save_row=True)
    gauge_iceland, _, _, _ = run_solver(
        H, mask, open_mask, f_cor, dx, g_eq, r_drag, nu, tide_amp,
        forcing_fn, n_steps, dt, save_every, j_iceland, i0, ny, nx,
        tilt_field=tilt, save_row=False)
    gauge_azores, _, _, _ = run_solver(
        H, mask, open_mask, f_cor, dx, g_eq, r_drag, nu, tide_amp,
        forcing_fn, n_steps, dt, save_every, j_azores, i0, ny, nx,
        tilt_field=tilt, save_row=False)
    print(f"  regenerated in {time.time()-t0:.1f}s, max|zeta|={max_abs:.4f}")

    gauge_dates = dates[0] + np.arange(len(gauge_baltic)) / 12.0
    gauge_forcing = np.interp(gauge_dates, dates, forcing)

    x_cols = np.arange(nx)
    k_shift = 0.1
    virtual_i = np.clip(i0 + k_shift * gauge_forcing, 1.0, nx - 2.0)
    moving = np.array([np.interp(vi, x_cols, row_baltic[t])
                       for t, vi in enumerate(virtual_i)])
    dipole_series = dipole(gauge_azores, gauge_iceland)

    return gauge_dates, gauge_forcing, gauge_baltic, moving, dipole_series


def robust_z(model_signed, obs_s, seeds, n_trials=32):
    """A single null_baseline draw at n_trials=16 was found (this same
    session) to overstate z by ~0.6-0.9 for this data -- average several
    independent draws at a larger n_trials for an honest number, rather
    than report whatever the first seed happens to give."""
    d_true, _ = dtw_distance(model_signed, obs_s)
    zs = []
    for seed in seeds:
        rng = np.random.default_rng(seed)
        null_mean, null_std = null_baseline(rng, model_signed, obs_s,
                                            n_trials=n_trials)
        zs.append((null_mean - d_true) / null_std if null_std > 0 else float("nan"))
    return float(np.mean(zs)), float(np.std(zs))


def corr_dtw(name, series, obs_ref, seed_offset=0):
    r_raw = float(np.corrcoef(series, obs_ref)[0, 1])
    sign = -1.0 if r_raw < 0 else 1.0
    print(f"\n-- {name} vs real Baltic --")
    print(f"  Pearson r: raw={r_raw:+.3f}  sign-matched={sign*r_raw:+.3f}")

    model_s = standardize(series)
    obs_s = standardize(obs_ref)
    model_signed = sign * model_s
    seeds = [seed_offset, seed_offset + 100, seed_offset + 200]
    z_raw, z_raw_std = robust_z(model_signed, obs_s, seeds)

    win = 3 * 12
    m3 = standardize(smooth(series, win))
    o3 = standardize(smooth(obs_ref, win))
    r3 = float(np.corrcoef(m3, o3)[0, 1])
    sign3 = -1.0 if r3 < 0 else 1.0
    z3, z3_std = robust_z(sign3 * m3, o3, seeds)

    print(f"  DTW raw-monthly z-score: {z_raw:+.2f} (+/-{z_raw_std:.2f} "
          f"across seeds)   DTW 3yr-smoothed z-score: {z3:+.2f} "
          f"(+/-{z3_std:.2f}) (r_3yr={sign3*r3:+.3f})")
    return dict(r_raw=r_raw, sign=sign, z_raw=z_raw, z3=z3, r3=sign3 * r3)


def ridge_check(name, series, dates_local, forcing_local, m_max, dm):
    ranked = rank_series(dates_local, series, forcing_local, m_max=m_max,
                         dm=dm, sigma=15.0, t0_step=5.0)
    near_125 = [r for r in ranked if abs(r["M"] - 1.25) < 0.03]
    print(f"  [{name}] winding_rank vs BALTIC's own Forcing (m_max={m_max}): "
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
    return ranked


def plot_result(gauge_dates, fixed, moving, dipole_series, variants,
                dtw_results, out_path: Path, m_max=1.5, dm=0.005) -> None:
    """Mirrors nao_resolved_shallow_water.py's plot style: a basin map,
    the raw undetrended gauges, one sign-matched real-data overlay panel
    per variant, and a ridge-scan panel (mode_evidence.sweep_M on the
    linear-Forcing residual, same convention as nao/amo's own ridge
    panels) for the headline 12mo-diff finding, marking both baltic's
    real M=1.250 and this run's own M=1.29 near-miss."""
    _, mask, _, _, _, _, _ = build_grid(300.0, 25.0, 70.0)
    H = np.full(mask.shape, 300.0)

    n_panels = 2 + len(variants) + 1
    fig, axes = plt.subplots(n_panels, 1, figsize=(11, 3.1 * n_panels))
    depth_show = np.where(mask, H, np.nan)
    im = axes[0].imshow(depth_show, origin="lower", cmap="Blues", aspect="auto")
    axes[0].set_title("atmospheric basin driven by BALTIC's own manifold "
                     "(equiv. depth 300m, full Coriolis)", fontsize=9)
    fig.colorbar(im, ax=axes[0], shrink=0.7, orientation="horizontal", pad=0.3)

    axes[1].plot(gauge_dates, fixed, linewidth=0.5, label="FIXED (58N)")
    axes[1].plot(gauge_dates, moving, linewidth=0.5, alpha=0.8, label="MOVING")
    axes[1].plot(gauge_dates, dipole_series, linewidth=0.5, alpha=0.8,
                label="DIPOLE (Azores-Iceland)")
    axes[1].set_title("raw (undetrended) gauges", fontsize=9)
    axes[1].set_xlabel("year")
    axes[1].legend(fontsize=7)

    def zscore(x):
        return (x - np.mean(x)) / np.std(x)

    diff_series = None
    diff_dates_local = None
    for ax, (name, series, dts, obs, forc), res in zip(
            axes[2:2 + len(variants)], variants, dtw_results):
        sign = res["sign"]
        r_overlay = sign * res["r_raw"]
        ax.plot(dts, zscore(obs), color="0.55", linewidth=0.8,
               label="real Baltic (standardized)")
        ax.plot(dts, sign * zscore(series), color="tab:blue", linewidth=1.0,
               alpha=0.85, label=f"{name} (r={r_overlay:+.3f}, "
                                  f"DTW_3yr z={res['z3']:+.2f})")
        ax.axhline(0, color="k", linewidth=0.4)
        ax.set_title(name, fontsize=9)
        ax.legend(fontsize=7)
        ax.set_xlabel("year")
        if "12mo" in name:
            diff_series, diff_dates_local, diff_forc = series, dts, forc

    ax = axes[-1]
    A = np.column_stack([np.ones_like(diff_forc), diff_forc])
    coef, *_ = np.linalg.lstsq(A, diff_series, rcond=None)
    resid = diff_series - A @ coef
    m_grid = np.arange(dm, m_max + dm / 2, dm)
    m_grid_out, curve = sweep_M(diff_forc, resid, (dm, m_max), len(m_grid))
    ax.plot(m_grid_out, np.abs(curve), color="tab:green")
    ylim_top = ax.get_ylim()[1]
    for m, lbl, c, y_frac in [(1.250, "baltic real M=1.25", "tab:red", 0.95),
                              (1.290, "this run's M=1.29", "tab:purple", 0.80)]:
        ax.axvline(m, color=c, linestyle=":", alpha=0.8)
        ax.text(m + 0.015, ylim_top * y_frac, lbl, color=c, fontsize=7)
    ax.set_xlabel("winding number M")
    ax.set_ylabel("|ridge scan cc|")
    ax.set_title("FIXED, 12mo delayed-difference: ridge scan vs BALTIC's "
                "own Forcing", fontsize=9)

    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"\n  saved {out_path}")


def main() -> int:
    cache = ROOT / "baltic" / "baltic_atmospheric_check_cache.npz"
    if cache.exists():
        print(f"Loading cached atmospheric-basin-driven-by-baltic output "
              f"from {cache}...")
        d = np.load(cache)
        gauge_dates, gauge_forcing, fixed, moving, dipole_series = (
            d["gauge_dates"], d["gauge_forcing"], d["fixed"], d["moving"],
            d["dipole_series"])
    else:
        print("Regenerating atmospheric basin driven by BALTIC's own "
              "manifold...")
        gauge_dates, gauge_forcing, fixed, moving, dipole_series = regenerate()
        np.savez(cache, gauge_dates=gauge_dates, gauge_forcing=gauge_forcing,
                fixed=fixed, moving=moving, dipole_series=dipole_series)
        print(f"  cached to {cache}")

    diff = delayed_difference(fixed, 12)
    diff_dates = gauge_dates[12:]
    diff_forcing = gauge_forcing[12:]

    baltic_dat = np.loadtxt(ROOT / "baltic" / "baltic.dat")
    obs_on_gauge = np.interp(gauge_dates, baltic_dat[:, 0], baltic_dat[:, 1])
    obs_on_diff = np.interp(diff_dates, baltic_dat[:, 0], baltic_dat[:, 1])

    variants = [
        ("FIXED, Baltic-latitude (58N)", fixed, gauge_dates, obs_on_gauge, gauge_forcing),
        ("MOVING (k_shift=0.1)", moving, gauge_dates, obs_on_gauge, gauge_forcing),
        ("SPATIAL DIPOLE (Azores-Iceland)", dipole_series, gauge_dates, obs_on_gauge, gauge_forcing),
        ("FIXED, 12mo delayed-difference", diff, diff_dates, obs_on_diff, diff_forcing),
    ]

    print("\n================ correlation + DTW vs real Baltic MSL "
          "================")
    dtw_results = []
    for i, (name, series, dts, obs, _) in enumerate(variants):
        dtw_results.append(corr_dtw(name, series, obs, seed_offset=10 * i))

    print("\n================ ridge scan vs BALTIC's OWN Forcing "
          "(m_max=1.5, testing for M=1.25) ================")
    for name, series, dts, _, forc in variants:
        ridge_check(name, series, dts, forc, m_max=1.5, dm=0.005)

    print("\n(reference: baltic's own independently-fitted winding_rank "
          "ridge: M=1.250, peak=4.21 bits, fwhm=0.040, continuity=1.00, "
          "PASS; real NAO vs real Baltic, raw: r=-0.279)")

    plot_result(gauge_dates, fixed, moving, dipole_series, variants,
               dtw_results, ROOT / "baltic" / "baltic_atmospheric_check.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
