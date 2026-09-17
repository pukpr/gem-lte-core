#!/usr/bin/env python3
"""nino4_beta_plane.py — the revised resolved-fluid ENSO run set up from
the AMO/ENSO reconciliation: beta-plane Coriolis (f=beta*y, NOT the
literal f=0 of nino4_resolved_shallow_water.py) plus a real zonal
thermocline tilt (deep warm-pool west, shallow cold-tongue east), keeping
the linear (Rayleigh) damping that reconciliation concluded is the
physically right choice here (bottom friction has no channel to act
through at ENSO's ~100-200m thermocline depth, thousands of meters above
the actual seafloor -- unlike the Baltic's shallow, bottom-coupled
setting).

Why change away from zero Coriolis at all, given it was exactly what
fixed the earlier long-duration numerical instability (see
[[nino4-resolved-shallow-water]])? Because f=0 everywhere doesn't just
remove a numerical headache -- it also removes the ONE piece of physics
(the beta-effect, f varying with latitude) responsible for real ENSO
boundary dynamics: Kelvin waves reflecting into Rossby waves at the
eastern boundary and vice versa at the western boundary is the actual
"delay" in the delayed-oscillator picture of ENSO, and that conversion
requires a latitude-varying Coriolis parameter, not its absence. A
zero-Coriolis channel is numerically convenient but physically a
different (non-rotating wave-channel) system, not a stand-in for
equatorial wave dynamics.

The reconciliation's other conclusion carries over unchanged: don't
expect this run to reproduce anything AMO-like at 60-120 year
timescales. ENSO's own natural bandwidth (Kelvin transit ~2-3 months,
first Rossby mode ~6-12 months) is what the delayed-oscillator mechanism
runs on, and that's the band this run's own free oscillation periods
should be judged against -- 2-7 years, not decades.

Geometry: same overall channel as nino4_resolved_shallow_water.py
(~130E-80W, closed at both zonal ends), but H(x) now ramps from a
warm-pool-like 150m in the west to a cold-tongue-like 40m in the east
(a linear ramp -- a real profile has a sharper eastern transition, but a
linear ramp is the honest "modest compute" first pass, not claimed to be
survey-accurate). f(y) = beta*y, y measured in meters from the channel's
center row, zero exactly at the equator and small (not zero) at the
meridional edges of this narrow band.

Usage
-----
    ./nino4_beta_plane.py --sanity
    ./nino4_beta_plane.py --years 144
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from baltic_resolved_shallow_water import laplacian, grad, run_solver, G  # noqa: E402
from nino4_resolved_shallow_water import x_norm_field                    # noqa: E402
from cv_rolling_blocked import prepare                                   # noqa: E402
from mode_evidence import sweep_M                                        # noqa: E402

ROOT = Path(__file__).resolve().parent
OMEGA_EARTH = 7.2921e-5   # rad/s
R_EARTH = 6.371e6         # m
BETA = 2.0 * OMEGA_EARTH / R_EARTH  # df/dy at the equator, ~2.29e-11 /(m*s)


def build_grid(dx_km: float = 200.0, h_west: float = 150.0,
               h_east: float = 40.0):
    """Same channel footprint as nino4_resolved_shallow_water.py, but with
    a real zonal thermocline tilt: linear ramp from h_west (warm pool) to
    h_east (cold tongue) -- not survey-accurate (the real transition is
    sharper near the eastern boundary), but a defensible first pass."""
    nx, ny = 110, 9
    dx = dx_km * 1e3
    depth_profile = np.linspace(h_west, h_east, nx)
    H = np.tile(depth_profile, (ny, 1))
    mask = np.ones((ny, nx), dtype=bool)
    mask[:, 0] = False
    mask[:, -1] = False
    mask[0, :] = False
    mask[-1, :] = False
    open_mask = np.zeros((ny, nx), dtype=bool)
    return H, mask, open_mask, dx, nx, ny


def beta_plane_fcor(ny: int, nx: int, dx: float) -> np.ndarray:
    """f = beta*y, y in meters from the channel's center row (the
    equator) -- zero at the center, small but nonzero at the meridional
    edges of this narrow band, unlike f=0 everywhere."""
    j = np.arange(ny)
    y = (j - (ny - 1) / 2.0) * dx
    f = BETA * y
    return np.repeat(f[:, None], nx, axis=1)


def sanity_check(dx_km: float, g_reduced: float, r_drag: float, nu_scale: float,
                 h_west: float, h_east: float):
    H, mask, open_mask, dx, nx, ny = build_grid(dx_km, h_west, h_east)
    f_cor = beta_plane_fcor(ny, nx, dx)
    c_max = np.sqrt(g_reduced * h_west)  # fastest wave speed sets the CFL limit
    dt = 0.4 * dx / c_max
    nu = nu_scale * dx ** 2 / dt
    print("-- sanity check: unforced pulse, beta-plane Coriolis + zonal tilt --")
    print(f"  grid {nx}x{ny} at dx={dx_km:g}km, depth {h_west:.0f}m(W)->{h_east:.0f}m(E), "
          f"c_max={c_max:.3f} m/s, dt={dt/3600:.2f} hr (CFL Courant={c_max*dt/dx:.2f})")
    print(f"  f_cor range: {f_cor.min():.2e} to {f_cor.max():.2e} rad/s "
          f"(vs Baltic's full-latitude ~1.2e-4 rad/s -- "
          f"{1.2e-4/max(abs(f_cor.min()), abs(f_cor.max()), 1e-30):.0f}x smaller)")

    j0, i0 = ny // 2, 40
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


def run(years: float, dx_km: float, drho_rho: float, r_drag: float,
       nu_scale: float, tide_amp: float, h_west: float, h_east: float,
       idx: str, m_max: float, dm: float, outdir: Path | None) -> None:
    H, mask, open_mask, dx, nx, ny = build_grid(dx_km, h_west, h_east)
    f_cor = beta_plane_fcor(ny, nx, dx)
    g_reduced = G * drho_rho
    c_max = np.sqrt(g_reduced * h_west)
    dt = 0.4 * dx / c_max
    nu = nu_scale * dx ** 2 / dt

    prep = prepare(idx)
    dates, forcing = prep["dates"], prep["forcing"]
    t_years = dates - dates[0]

    def forcing_fn(t_sec: float) -> float:
        t_yr = t_sec / (365.2422 * 86400.0)
        return float(np.interp(t_yr, t_years, forcing))

    dt_out_target = 1.0 / 12.0 * 365.2422 * 86400.0
    save_every = max(1, round(dt_out_target / dt))
    n_steps = int(years * 365.2422 * 86400.0 / dt)
    n_steps -= n_steps % save_every

    gauge_j, gauge_i = ny // 2, 40
    print(f"-- resolved shallow-water run: {idx} (beta-plane Coriolis + zonal tilt) --")
    print(f"  grid {nx}x{ny} at dx={dx_km:g}km, depth {h_west:.0f}m(W)->{h_east:.0f}m(E), "
          f"c_max={c_max:.3f} m/s (g'={g_reduced:.4f})")
    print(f"  f_cor range: {f_cor.min():.2e} to {f_cor.max():.2e} rad/s")
    print(f"  dt={dt/3600:.2f} hr, {n_steps} steps ({years:g} sim-years), "
          f"r_drag={r_drag:g}, nu_scale={nu_scale:g}")

    t0 = time.time()
    gauge, max_abs, _, _ = run_solver(H, mask, open_mask, f_cor, dx, g_reduced,
                                   r_drag, nu, tide_amp, forcing_fn, n_steps, dt,
                                   save_every, gauge_j, gauge_i, ny, nx,
                                   tilt_field=x_norm_field(mask))
    elapsed = time.time() - t0
    print(f"  done in {elapsed:.1f}s wall-clock; max|zeta| reached {max_abs:.4f} m "
          f"(forcing amplitude scale: {tide_amp:.3f} m)")

    gauge_dates = dates[0] + np.arange(len(gauge)) / 12.0
    gauge_forcing = np.interp(gauge_dates, dates, forcing)

    A = np.column_stack([np.ones_like(gauge_forcing), gauge_forcing])
    coef, _, _, _ = np.linalg.lstsq(A, gauge, rcond=None)
    linear_r2 = 1.0 - np.var(gauge - A @ coef) / np.var(gauge)
    gauge_resid = gauge - A @ coef
    print(f"  linear (a + b*Forcing) fit alone explains R^2={linear_r2:.3f} "
          f"of the gauge -- scanning the RESIDUAL for a genuine M-specific ridge")

    m_grid = np.arange(dm, m_max + dm / 2, dm)
    m_grid_out, curve = sweep_M(gauge_forcing, gauge_resid, (dm, m_max), len(m_grid))
    i_star = int(np.argmax(np.abs(curve)))
    print(f"  virtual-gauge ridge scan (on residual): M_hat={m_grid_out[i_star]:.3f}  "
          f"cc={curve[i_star]:+.3f}")
    print(f"  (for comparison: real nino4's independently-found ridges are "
          f"M*=0.450, per cv_ridge_transfer.py -- though not AR1-significant there. "
          f"This run should be judged against ENSO's own 2-7yr band, not that value.")

    # also report the gauge's own dominant period(s) via a simple periodogram,
    # since that's the ENSO-relevant (2-7yr) diagnostic this run should
    # actually be judged against, not the AMO-derived M* comparison above
    report_periods(idx, gauge_dates, gauge_resid)

    out_dir = outdir if outdir is not None else ROOT / idx
    plot_result(idx, H, mask, gauge_dates, gauge, gauge_resid, m_grid_out, curve,
               out_dir / "nino4_beta_plane.png")


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
    band = (period > 1.0) & (period < 10.0)
    if not np.any(band):
        return
    top = np.argsort(-power[band])[:5]
    print(f"  top periods in the ENSO-relevant 1-10yr band: " +
          ", ".join(f"{p:.2f}yr" for p in period[band][top]))


def plot_result(idx, H, mask, gauge_dates, gauge, gauge_resid, m_grid, curve,
                out_path: Path) -> None:
    fig, axes = plt.subplots(4, 1, figsize=(11, 10.5))
    depth_show = np.where(mask, H, np.nan)
    im = axes[0].imshow(depth_show, origin="lower", cmap="Blues", aspect="auto")
    axes[0].set_title("idealized equatorial Pacific channel (depth, m): "
                     "warm-pool west -> cold-tongue east, beta-plane Coriolis",
                     fontsize=9)
    fig.colorbar(im, ax=axes[0], shrink=0.7, orientation="horizontal", pad=0.3)

    axes[1].plot(gauge_dates, gauge, linewidth=0.5)
    axes[1].set_xlabel("year")
    axes[1].set_ylabel("virtual gauge zeta (m)")
    axes[1].set_title(f"{idx}: simulated thermocline displacement", fontsize=9)

    axes[2].plot(gauge_dates, gauge_resid, linewidth=0.6, color="tab:orange")
    axes[2].set_xlabel("year")
    axes[2].set_ylabel("residual (m)")
    axes[2].set_title(f"{idx}: residual after removing linear Forcing tracking",
                      fontsize=9)

    axes[3].plot(m_grid, np.abs(curve))
    axes[3].axvline(0.45, color="tab:red", linestyle=":",
                    label="real nino4's M*=0.450 (AMO-scale comparison only)")
    axes[3].set_xlabel("winding number M")
    axes[3].set_ylabel("|ridge scan cc|")
    axes[3].set_title(f"{idx}: virtual-gauge ridge scan", fontsize=9)
    axes[3].legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"\n  saved {out_path}")


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sanity", action="store_true")
    ap.add_argument("--idx", default="nino4")
    ap.add_argument("--years", type=float, default=144.0)
    ap.add_argument("--dx-km", type=float, default=200.0)
    ap.add_argument("--drho-rho", type=float, default=0.005)
    ap.add_argument("--h-west", type=float, default=150.0,
                     help="warm-pool (western) thermocline depth, m")
    ap.add_argument("--h-east", type=float, default=40.0,
                     help="cold-tongue (eastern) thermocline depth, m")
    ap.add_argument("--r-drag", type=float, default=1e-7)
    ap.add_argument("--nu-scale", type=float, default=0.02)
    ap.add_argument("--tide-amp", type=float, default=2e-2)
    ap.add_argument("--m-max", type=float, default=5.0)
    ap.add_argument("--dm", type=float, default=0.01)
    ap.add_argument("--outdir", type=Path, default=None)
    args = ap.parse_args()

    if args.sanity:
        g_reduced = G * args.drho_rho
        ok = sanity_check(args.dx_km, g_reduced, args.r_drag, args.nu_scale,
                          args.h_west, args.h_east)
        return 0 if ok else 1

    run(args.years, args.dx_km, args.drho_rho, args.r_drag, args.nu_scale,
       args.tide_amp, args.h_west, args.h_east, args.idx, args.m_max, args.dm,
       args.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
