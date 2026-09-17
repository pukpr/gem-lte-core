#!/usr/bin/env python3
"""amo_resolved_shallow_water.py — a resolved shallow-water run for AMO's
actual setting: a mid-latitude basin with FULL Coriolis (not beta-plane,
not zero) -- the regime that made baltic_resolved_shallow_water.py
unstable for anything beyond ~10-15 simulated years without either
degenerating (too much damping) or blowing up (too little).

The lever this script leans on, instead of re-fighting the same
damping-vs-diffusion trade-off that had no good window for Baltic:
Baltic's instability accumulated over ~700k+ timesteps at a SHALLOW-WATER
wave speed of only ~1.5-2 m/s (dx=25km, H~50m -- a coastal, fine-detail
basin). AMO/AMOC is a genuinely basin-scale phenomenon with no comparable
need for coastal detail, and its relevant reduced-gravity mode sits on a
much deeper main thermocline (several hundred meters, not tens). A
coarser grid at that physically appropriate depth gives a much faster
wave speed and a proportionally larger CFL-stable timestep -- fewer
accumulated steps for the same 144 real years, attacking the root cause
(accumulated numerical error over hundreds of thousands of steps) instead
of retuning damping again.

Geometry: an idealized, fully-closed rectangular North Atlantic-scale
basin (~5000km E-W, ~6000km N-S, spanning roughly 10N-65N) -- no strait/
coastline detail, since AMO is a basin-averaged SST index, not a
coastally-sensitive local record the way Baltic sea level is. Coriolis is
the FULL f=2*Omega*sin(lat) across this wide latitude span (a beta-plane
linearization is only valid for a narrow band around one reference
latitude, which doesn't apply across most of the Atlantic's real extent
-- unlike the nino4 case, this is genuinely a full-Coriolis problem, not
a beta-plane one). Forcing is a north-south (meridional) tilt of the real
AMO Forcing(t), matching the AMOC's own meridional character (reusing
baltic_resolved_shallow_water.y_norm_field, unlike nino4's zonal tilt).

Usage
-----
    ./amo_resolved_shallow_water.py --sanity
    ./amo_resolved_shallow_water.py --years 144
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from baltic_resolved_shallow_water import (laplacian, grad, run_solver, G,  # noqa: E402
                                           OMEGA_EARTH, y_norm_field)
from mode_evidence import sweep_M                                          # noqa: E402
from lte_forward import tide_sum, impulse_delta, iir, bessel               # noqa: E402

ROOT = Path(__file__).resolve().parent


def build_forcing_at_yl(idx: str, yl: float):
    """Rebuilds the forcing manifold at an EXPLICIT year length, exactly
    matching fig_yl_sweep.py's manifold(yl) (the comb-alias study run
    earlier today) -- bypassing cv_rolling_blocked.prepare()'s own
    default year length, which is NOT the yl=365.2463 "PRODUCTION" value
    that study found gives the sharpest ~124yr comb alias (r=0.704 vs
    0.311/0.523/0.434 for nearby yl choices)."""
    idx_dir = ROOT / idx
    params = json.loads((idx_dir / "lt.exe.p").read_text())
    lpap = np.array(params["lpap"], dtype=float)
    periods = np.abs(lpap[:, 0])
    B = {k: float(params.get(k, 0.0)) for k in
         ("offs", "bg", "impA", "impB", "delA", "delB", "asym", "ma", "mp",
          "shfT", "init")}
    k1 = params["ltep"][1]
    dat = np.loadtxt(idx_dir / f"{idx}.dat")
    dates, obs = dat[:, 0].copy(), dat[:, 1].copy()
    tf = tide_sum(dates, lpap[:, 1:3], periods, yl, 0.0, B["shfT"])
    comb = impulse_delta(dates, B["delA"], B["delB"], B["asym"], 12)
    R = iir(tf * comb, lag_a=1.0 - B["ma"], lag_c=B["mp"], init=B["init"],
           start_date=1880.0, dates=dates)
    forcing = bessel(R, B["impA"], B["impB"], k1, B["offs"], B["bg"])
    return dates, forcing, obs


def build_grid(dx_km: float = 250.0, lat_south: float = 10.0,
               lat_north: float = 65.0):
    """Idealized, fully-closed North Atlantic-scale rectangular basin --
    no coastline/strait detail (AMO is basin-averaged, not coastally
    local), spanning a wide enough latitude range that full f=2*Omega*
    sin(lat) matters, not a beta-plane approximation."""
    nx, ny = 20, 24   # ~5000km x 6000km at dx=250km
    dx = dx_km * 1e3
    H = np.full((ny, nx), np.nan)  # filled by caller with a depth value
    mask = np.ones((ny, nx), dtype=bool)
    mask[:, 0] = False
    mask[:, -1] = False
    mask[0, :] = False
    mask[-1, :] = False
    open_mask = np.zeros((ny, nx), dtype=bool)
    lat = lat_south + (lat_north - lat_south) * (np.arange(ny) / (ny - 1))
    return H, mask, open_mask, dx, nx, ny, lat


def full_fcor(lat: np.ndarray, nx: int) -> np.ndarray:
    f = 2.0 * OMEGA_EARTH * np.sin(np.radians(lat))
    return np.repeat(f[:, None], nx, axis=1)


def sanity_check(dx_km: float, g_reduced: float, depth: float, r_drag: float,
                 nu_scale: float, lat_south: float, lat_north: float):
    _, mask, open_mask, dx, nx, ny, lat = build_grid(dx_km, lat_south, lat_north)
    H = np.full((ny, nx), depth)
    f_cor = full_fcor(lat, nx)
    c = np.sqrt(g_reduced * depth)
    dt = 0.4 * dx / c
    nu = nu_scale * dx ** 2 / dt
    print("-- sanity check: unforced pulse, full mid-latitude Coriolis --")
    print(f"  grid {nx}x{ny} at dx={dx_km:g}km, depth={depth:.0f}m, "
          f"c={c:.3f} m/s, dt={dt/3600:.2f} hr (CFL Courant={c*dt/dx:.2f})")
    print(f"  latitude range {lat_south:.0f}N-{lat_north:.0f}N, "
          f"f_cor range: {f_cor.min():.2e} to {f_cor.max():.2e} rad/s "
          f"(Baltic's was ~1.2e-4 rad/s for comparison)")

    j0, i0 = ny // 2, nx // 2
    n_steps = int(60 * 86400 / dt)
    diag_every = max(1, n_steps // 20)
    init_zeta = np.zeros((ny, nx))
    init_zeta[j0, i0] = 1.0
    _, max_abs, diag, _ = run_solver(
        H, mask, open_mask, f_cor, dx, g_reduced, r_drag, nu,
        tide_amp=0.0, forcing_fn=lambda t: 0.0, n_steps=n_steps, dt=dt,
        save_every=n_steps, gauge_j=j0, gauge_i=i0, ny=ny, nx=nx,
        diag_every=diag_every, init_zeta=init_zeta)

    print("  day   total energy (should decay smoothly, not blow up)   max|zeta|")
    for day, e, mz in diag:
        print(f"    {day:6.1f}  {e:.4e}  {mz:.4f}")
    print(f"  peak max|zeta| over the whole run = {max_abs:.4f} m "
          f"(started at 1.0 m pulse -- no interim blow-up)")
    return np.isfinite(max_abs) and max_abs < 10.0


def run(years: float, dx_km: float, drho_rho: float, depth: float, r_drag: float,
       nu_scale: float, tide_amp: float, lat_south: float, lat_north: float,
       idx: str, yl: float, m_max: float, dm: float, outdir: Path | None) -> None:
    _, mask, open_mask, dx, nx, ny, lat = build_grid(dx_km, lat_south, lat_north)
    H = np.full((ny, nx), depth)
    f_cor = full_fcor(lat, nx)
    g_reduced = G * drho_rho
    c = np.sqrt(g_reduced * depth)
    dt = 0.4 * dx / c
    nu = nu_scale * dx ** 2 / dt

    dates, forcing, obs = build_forcing_at_yl(idx, yl)
    print(f"  forcing manifold built at year length yl={yl:.4f} days "
          f"(PRODUCTION comb-alias value -> ~124yr alias, ~62yr doubled)")
    t_years = dates - dates[0]

    def forcing_fn(t_sec: float) -> float:
        t_yr = t_sec / (365.2422 * 86400.0)
        return float(np.interp(t_yr, t_years, forcing))

    dt_out_target = 1.0 / 12.0 * 365.2422 * 86400.0
    save_every = max(1, round(dt_out_target / dt))
    n_steps = int(years * 365.2422 * 86400.0 / dt)
    n_steps -= n_steps % save_every

    gauge_j, gauge_i = ny // 2, nx // 2
    print(f"-- resolved shallow-water run: {idx} (full mid-latitude Coriolis) --")
    print(f"  grid {nx}x{ny} at dx={dx_km:g}km, depth={depth:.0f}m, c={c:.3f} m/s "
          f"(g'={g_reduced:.4f}), lat {lat_south:.0f}N-{lat_north:.0f}N")
    print(f"  f_cor range: {f_cor.min():.2e} to {f_cor.max():.2e} rad/s")
    print(f"  dt={dt/3600:.2f} hr, {n_steps} steps ({years:g} sim-years), "
          f"r_drag={r_drag:g}, nu_scale={nu_scale:g}")

    t0 = time.time()
    gauge, max_abs, _, row_hist = run_solver(
        H, mask, open_mask, f_cor, dx, g_reduced, r_drag, nu, tide_amp,
        forcing_fn, n_steps, dt, save_every, gauge_j, gauge_i, ny, nx,
        tilt_field=y_norm_field(mask), save_row=True)
    elapsed = time.time() - t0
    print(f"  done in {elapsed:.1f}s wall-clock; max|zeta| reached {max_abs:.4f} m "
          f"(forcing amplitude scale: {tide_amp:.3f} m)")

    gauge_dates = dates[0] + np.arange(len(gauge)) / 12.0
    gauge_forcing = np.interp(gauge_dates, dates, forcing)

    A = np.column_stack([np.ones_like(gauge_forcing), gauge_forcing])
    coef, _, _, _ = np.linalg.lstsq(A, gauge, rcond=None)
    linear_r2 = 1.0 - np.var(gauge - A @ coef) / np.var(gauge)
    gauge_resid = gauge - A @ coef
    print(f"  [FIXED gauge] linear (a + b*Forcing) fit alone explains "
          f"R^2={linear_r2:.3f} -- scanning the RESIDUAL for a genuine "
          f"M-specific ridge")

    m_grid = np.arange(dm, m_max + dm / 2, dm)
    m_grid_out, curve = sweep_M(gauge_forcing, gauge_resid, (dm, m_max), len(m_grid))
    i_star = int(np.argmax(np.abs(curve)))
    print(f"  [FIXED gauge] ridge scan (on residual): M_hat={m_grid_out[i_star]:.3f}  "
          f"cc={curve[i_star]:+.3f}")

    # -----------------------------------------------------------------
    # MOVING gauge: sample along the same row, but at a position that
    # SHIFTS with the forcing itself -- x(t) = gauge_i + k_shift*F(t),
    # interpolated linearly between grid columns. This is the actual
    # material-coordinate / Lagrangian-pullback mechanism the original
    # sin(k*M(t)) derivation rests on (Phi(phi)=sin(k*phi) sampled along
    # phi=M(t)) -- a fixed-point gauge additively forced by F(t), which is
    # what every resolved run so far has used, structurally cannot
    # reproduce that nonlinearity, since the model's own physics is
    # linear throughout. This does not require removing the linear trend
    # first -- the moving-coordinate sampling IS the nonlinearity, not a
    # residual left over after one.
    # -----------------------------------------------------------------
    x_cols = np.arange(nx)
    for k_shift in (0.05, 0.1, 0.2, 0.3):
        virtual_i = np.clip(gauge_i + k_shift * gauge_forcing, 1.0, nx - 2.0)
        moving = np.array([np.interp(vi, x_cols, row_hist[t])
                           for t, vi in enumerate(virtual_i)])
        coef_mv, _, _, _ = np.linalg.lstsq(A, moving, rcond=None)
        mv_r2 = 1.0 - np.var(moving - A @ coef_mv) / np.var(moving)
        moving_resid = moving - A @ coef_mv
        m_grid_mv, curve_mv = sweep_M(gauge_forcing, moving_resid, (dm, m_max), len(m_grid))
        j_star = int(np.argmax(np.abs(curve_mv)))
        print(f"  [MOVING gauge, k_shift={k_shift:.2f}] linear fit R^2={mv_r2:.3f} "
              f"-- ridge scan on residual: M_hat={m_grid_mv[j_star]:.3f}  "
              f"cc={curve_mv[j_star]:+.3f}")
        report_periods(f"{idx} (moving k={k_shift:.2f})", gauge_dates, moving_resid)
        if k_shift == 0.1:
            moving_for_plot = moving_resid
            m_grid_mv_plot, curve_mv_plot = m_grid_mv, curve_mv

    print(f"  (for comparison: real amo's own fitted ridges are M=0.0134 "
          f"and M=0.2075, per today's comb-alias study)")

    report_periods(idx, gauge_dates, gauge_resid)

    obs_on_gauge = np.interp(gauge_dates, dates, obs)

    out_dir = outdir if outdir is not None else ROOT / idx
    plot_result(idx, H, mask, gauge_dates, gauge, gauge_resid, moving_for_plot,
               obs_on_gauge, m_grid_out, curve, m_grid_mv_plot, curve_mv_plot,
               out_dir / "amo_resolved_shallow_water.png")


def report_periods(idx: str, dates: np.ndarray, resid: np.ndarray) -> None:
    from scipy.signal.windows import tukey
    n = len(resid)
    dt_yr = np.median(np.diff(dates))
    y = resid - np.polyval(np.polyfit(dates, resid, 1), dates)
    w = tukey(n, 0.1)
    yw = y * w * np.sqrt(n / np.sum(w ** 2))
    freq = np.fft.rfftfreq(n, d=dt_yr)[1:]
    power = (np.abs(np.fft.rfft(yw)) ** 2)[1:]
    period = 1.0 / freq
    band = (period > 20) & (period < 200)
    if not np.any(band):
        return
    top = np.argsort(-power[band])[:5]
    print(f"  top periods in the AMO-relevant 20-200yr band: " +
          ", ".join(f"{p:.0f}yr" for p in period[band][top]))


def plot_result(idx, H, mask, gauge_dates, gauge, gauge_resid, moving, obs,
                m_grid, curve, m_grid_mv, curve_mv, out_path: Path) -> None:
    fig, axes = plt.subplots(6, 1, figsize=(11, 15.5))
    depth_show = np.where(mask, H, np.nan)
    im = axes[0].imshow(depth_show, origin="lower", cmap="Blues", aspect="auto")
    axes[0].set_title("idealized North Atlantic-scale basin (uniform depth, "
                     "full mid-latitude Coriolis)", fontsize=9)
    fig.colorbar(im, ax=axes[0], shrink=0.7, orientation="horizontal", pad=0.3)

    axes[1].plot(gauge_dates, gauge, linewidth=0.5, label="fixed gauge")
    axes[1].plot(gauge_dates, moving, linewidth=0.5, color="tab:green",
                alpha=0.8, label="moving gauge (k_shift=0.1, detrended)")
    axes[1].set_xlabel("year")
    axes[1].set_ylabel("zeta (m)")
    axes[1].set_title(f"{idx}: fixed vs MOVING (material-coordinate) gauge",
                      fontsize=9)
    axes[1].legend(fontsize=7)

    def zscore(x):
        return (x - np.mean(x)) / np.std(x)

    # sign is an arbitrary convention here (nothing pins down whether
    # +zeta should mean warm or cool AMO phase -- it's set by the tilt
    # field's own arbitrary orientation, not a fitted choice), so flip to
    # whichever sign makes the physical match visible rather than hiding
    # a real pattern behind a coin-flip convention.
    r_raw = float(np.corrcoef(moving, obs)[0, 1])
    sign = -1.0 if r_raw < 0 else 1.0
    r_overlay = sign * r_raw
    axes[2].plot(gauge_dates, zscore(obs), color="0.55", linewidth=0.8,
                label="real AMO (standardized)")
    axes[2].plot(gauge_dates, sign * zscore(moving), color="tab:green",
                linewidth=1.0, alpha=0.85,
                label=f"moving gauge, sign-matched (r={r_overlay:+.3f})")
    axes[2].axhline(0, color="k", linewidth=0.4)
    axes[2].set_xlabel("year")
    axes[2].set_ylabel("standardized units")
    axes[2].set_title(f"{idx}: MOVING gauge overlaid on real AMO data "
                     f"(both z-scored, sign-matched)", fontsize=9)
    axes[2].legend(fontsize=7)

    axes[3].plot(gauge_dates, gauge_resid, linewidth=0.6, color="tab:orange")
    axes[3].set_xlabel("year")
    axes[3].set_ylabel("residual (m)")
    axes[3].set_title(f"{idx}: FIXED-gauge residual after removing linear "
                     f"Forcing tracking", fontsize=9)

    axes[4].plot(m_grid, np.abs(curve))
    for m, lbl in [(0.0134, "0.0134"), (0.2075, "0.2075")]:
        axes[4].axvline(m, color="tab:red", linestyle=":", alpha=0.7)
        axes[4].text(m, axes[4].get_ylim()[1]*0.9, lbl, color="tab:red", fontsize=7)
    axes[4].set_xlabel("winding number M")
    axes[4].set_ylabel("|ridge scan cc|")
    axes[4].set_title(f"{idx}: FIXED-gauge ridge scan on residual "
                     f"(red = amo's own fitted M)", fontsize=9)

    axes[5].plot(m_grid_mv, np.abs(curve_mv), color="tab:green")
    for m, lbl in [(0.0134, "0.0134"), (0.2075, "0.2075")]:
        axes[5].axvline(m, color="tab:red", linestyle=":", alpha=0.7)
        axes[5].text(m, axes[5].get_ylim()[1]*0.9, lbl, color="tab:red", fontsize=7)
    axes[5].set_xlabel("winding number M")
    axes[5].set_ylabel("|ridge scan cc|")
    axes[5].set_title(f"{idx}: MOVING-gauge ridge scan on residual "
                     f"(k_shift=0.1)", fontsize=9)

    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"\n  saved {out_path}")
    print(f"  overlay correlation (moving gauge detrended vs real AMO, "
          f"both standardized): r={r_overlay:+.3f}")


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sanity", action="store_true")
    ap.add_argument("--idx", default="amo")
    ap.add_argument("--years", type=float, default=144.0)
    ap.add_argument("--dx-km", type=float, default=250.0)
    ap.add_argument("--drho-rho", type=float, default=0.005)
    ap.add_argument("--depth", type=float, default=700.0,
                     help="reduced-gravity mode depth, m (AMOC main "
                          "thermocline scale, not Baltic's shallow halocline)")
    ap.add_argument("--lat-south", type=float, default=10.0)
    ap.add_argument("--lat-north", type=float, default=65.0)
    ap.add_argument("--r-drag", type=float, default=1e-7)
    ap.add_argument("--nu-scale", type=float, default=0.02)
    ap.add_argument("--tide-amp", type=float, default=2e-2)
    ap.add_argument("--yl", type=float, default=365.2463,
                     help="year length (days) the forcing manifold is built "
                          "at -- default is the comb-alias study's own "
                          "PRODUCTION value, giving the sharpest ~124yr "
                          "alias / ~62yr doubled period (r=0.704 there vs "
                          "0.311/0.523/0.434 for nearby values)")
    ap.add_argument("--m-max", type=float, default=5.0)
    ap.add_argument("--dm", type=float, default=0.01)
    ap.add_argument("--outdir", type=Path, default=None)
    args = ap.parse_args()

    if args.sanity:
        g_reduced = G * args.drho_rho
        ok = sanity_check(args.dx_km, g_reduced, args.depth, args.r_drag,
                          args.nu_scale, args.lat_south, args.lat_north)
        return 0 if ok else 1

    run(args.years, args.dx_km, args.drho_rho, args.depth, args.r_drag,
       args.nu_scale, args.tide_amp, args.lat_south, args.lat_north, args.idx,
       args.yl, args.m_max, args.dm, args.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
