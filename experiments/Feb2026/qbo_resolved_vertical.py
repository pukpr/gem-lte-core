#!/usr/bin/env python3
"""qbo_resolved_vertical.py — a "shallow-water equivalent" for QBO,
formulated from first principles rather than copying the Baltic/AMO/ENSO
geometry, because QBO's own physics is structurally different in a way
that changes what the resolved domain should even be.

Why the domain is VERTICAL, not horizontal
-------------------------------------------
Every other resolved run this session (Baltic, AMO, both ENSO versions)
solved a 2D HORIZONTAL (longitude x latitude) problem, because those
indices' defining structure is spatial in the horizontal. QBO is the
opposite: it is a ZONAL-MEAN (wavenumber s=0) quantity by definition --
u(z,t) at the equator, with no longitude dependence at all. In the
Chapter-12 LTE reduction (12-7), the longitude/zonal-wavenumber term is
`SW(s)*zeta`, coming from separating d^2/dlambda^2 into -s^2 for a zonal
mode e^{is*lambda}. At s=0, SW(0)=0 IDENTICALLY -- the zonal restoring
term vanishes completely, not approximately. So there is no meaningful
"horizontal channel" version of this problem the way there was for ENSO's
zonal thermocline tilt. What QBO actually varies over, and what its
defining phenomenon (downward-propagating wind regimes) lives in, is
ALTITUDE. So this script reuses the exact same generic solver
(run_solver, zero Coriolis -- QBO is equatorial, so this is exact, not
approximate, the same condition already used for nino4) but relabels the
zonal grid axis as a VERTICAL (log-pressure altitude) axis instead of
longitude. No new PDE code -- same machinery, different physical meaning
of the coordinate, because the math (a 1D wave equation with no rotation)
is identical either way.

Forcing: the draconic (27.2122d) Doodson line ALONE, per the user's own
stated reasoning (wavenumber=0 leaves this as the dominant relevant
term) -- but NOT as a bare sinusoid. The first version of this script
fed tide_sum's raw ~27-day cosine straight to the solver, which produces
a clean fast oscillation with nothing resembling QBO's real ~28-month
cadence. The fix, pointed out directly: convolve/amplitude-modulate that
fast tone with the REAL Impulse_Delta (twice-yearly Dirac comb) and then
push it through the REAL IIR (memory + Coulomb-drag) stage -- exactly
the two pipeline steps every other index in this project already uses
(lte_forward.forward, Calc_Forcing) and exactly the same comb-alias /
stroboscopic-sampling mechanism already established for AMO earlier this
session (a fast tidal line aliased against a coarser sampling comb
produces an emergent SLOW, jagged beat). Matches qbo_isolated_iir.py's
own precedent of running an isolated Doodson subset through this same
real, nonlinear two-stage pipeline rather than comparing a bare cosine.
qbo30 and qbo50's own fitted draconic amplitude/phase are nearly
identical (-3.319 vs -3.384, 35.166 vs 35.167) -- confirming this really
is the same astronomical input felt at both levels; whatever difference
exists between the two records has to come from the RESOLVED VERTICAL
STRUCTURE (plus each level's own delA/delB/asym/ma/mp comb+IIR
calibration), not from a different raw forcing.

Real target to explain, established directly from data before building
anything (not assumed): qbo50 LAGS qbo30 by 4 months (best-lag cross-
correlation r=0.802, vs r=0.353 at zero lag). At the log-pressure altitude
separation between 30hPa (~24.5km) and 50hPa (~21.0km) -- about 3.5km --
a 4-month lag implies a descent rate close to the real, documented QBO
descent rate of order 1 km/month. This is the number the moving
(descending) gauge is tested against.

Two gauges, and the actual question
------------------------------------
1. FIXED gauges at the grid points nearest 30hPa and 50hPa altitude --
   the direct "does the resolved model's own vertical wave structure,
   with no shifting at all, already produce a lag between two fixed
   levels" test.
2. A DESCENDING moving gauge (position shifts downward at a tunable rate)
   -- tests whether treating "the thing you sample" as a single trajectory
   descending through the column (rather than two independent fixed
   points) gives a cleaner or physically more direct account of the
   30->50hPa relationship, echoing the same material-coordinate mechanism
   used for AMO/ENSO but re-interpreted vertically, which is what QBO's
   own real phenomenology (descending regimes) actually is.

Usage
-----
    ./qbo_resolved_vertical.py --sanity
    ./qbo_resolved_vertical.py --idx qbo30 --years 70
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
from baltic_resolved_shallow_water import laplacian, grad, run_solver, G  # noqa: E402
from lte_forward import (tide_sum, impulse_delta, iir,                    # noqa: E402
                         year_length as lte_year_length)

ROOT = Path(__file__).resolve().parent
DRACONIC_PERIOD = 27.212220815  # days
H_SCALE = 7.0  # km, standard atmospheric scale height for log-pressure altitude


def p_to_z(p_hpa: float, p0: float = 1000.0) -> float:
    """Log-pressure altitude (km): z = H*ln(p0/p)."""
    return H_SCALE * np.log(p0 / p_hpa)


def build_grid(dx_km: float = 0.5, p_top: float = 10.0, p_bot: float = 100.0):
    """1D-equivalent vertical channel: nx points spanning log-pressure
    altitude from p_bot (bottom, high pressure) to p_top (top, low
    pressure), a trivial ny=3 second dimension purely so the existing 2D
    solver runs unmodified -- no meridional structure is being resolved
    (wavenumber=0 already removed the zonal dimension; the meridional one
    is not the point of this test)."""
    z_bot, z_top = p_to_z(p_bot), p_to_z(p_top)
    nx = int((z_top - z_bot) / dx_km) + 1
    ny = 3
    dx = dx_km * 1e3
    H = np.full((ny, nx), 1.0)  # uniform "depth" -- see module docstring:
                                # this is an illustrative wave-speed proxy,
                                # not a literal fluid depth
    mask = np.ones((ny, nx), dtype=bool)
    mask[:, 0] = False   # rigid lid, bottom of modeled layer
    mask[:, -1] = False  # rigid lid, top of modeled layer
    mask[0, :] = False
    mask[-1, :] = False
    open_mask = np.zeros((ny, nx), dtype=bool)
    z = z_bot + dx_km * np.arange(nx)
    return H, mask, open_mask, dx, nx, ny, z


def z_norm_field(mask: np.ndarray) -> np.ndarray:
    ny, nx = mask.shape
    i = np.arange(nx)
    return np.repeat(((i / (nx - 1)) - 0.5)[None, :], ny, axis=0)


def draconic_forcing(idx: str, dates: np.ndarray, yl: float | None = None):
    """Isolated draconic (27.2122d) Doodson term, pushed through the REAL
    Impulse_Delta (twice-yearly Dirac comb) -> IIR (memory + Coulomb-drag)
    stages -- matches qbo_isolated_iir.py's precedent of running an
    isolated Doodson subset through this same real, nonlinear two-stage
    pipeline rather than comparing a bare cosine. This is what turns the
    raw ~27-day draconic tone into a genuinely slow, jagged, multi-year
    envelope (the same comb-alias/stroboscopic mechanism established for
    AMO earlier this session), not a fast sinusoid.

    `dates` MUST be the real, native (monthly) sample grid the index's
    own delA/delB/asym comb calibration was fit against -- Impulse_Delta
    identifies "this sample's month" from where each date's fractional
    year falls among `sampling` (12) bins, so evaluating it on a much
    finer synthetic grid would smear each instantaneous monthly impulse
    into a month-wide block, which is a different pipeline than the real
    one. Returns the monthly-cadence forcing; the caller interpolates it
    onto the solver's own finer timestep (exactly as amo_resolved_
    shallow_water.py's run() does with its own monthly forcing series).

    `yl` (year length, days): defaults to this index's own real
    production value (YEAR_IN_DAYS + resp YEAR + lt.exe.p "year"
    candidate, via lte_forward.year_length) rather than a generic
    365.2422 guess; still overridable for a future yl-sweep, following
    AMO's build_forcing_at_yl convention (amp/phase reused as-is, only
    the freq=yl/period mapping varies)."""
    p = json.loads((ROOT / idx / "lt.exe.p").read_text())
    resp = {}
    resp_path = ROOT / idx / "lt.exe.resp"
    if resp_path.is_file():
        toks = resp_path.read_text().split()
        for i in range(0, len(toks) - 1, 2):
            if toks[i].upper() == "YEAR":
                resp["YEAR"] = float(toks[i + 1])
    if yl is None:
        yl = lte_year_length(resp.get("YEAR", 0.0), float(p.get("year", 0.0)))
    lpap = np.array(p["lpap"], dtype=float)
    row = min(lpap, key=lambda r: abs(abs(r[0]) - DRACONIC_PERIOD))
    period, amp, phase = row
    tf = tide_sum(dates, np.array([[amp, phase]]), np.array([DRACONIC_PERIOD]),
                 yl, 0.0, 0.0)
    comb = impulse_delta(dates, p["delA"], p["delB"], p["asym"], sampling=12)
    forced = iir(tf * comb, lag_a=1.0 - p["ma"], lag_c=p["mp"],
                init=p["init"], start_date=0.0, dates=dates)
    return forced, yl


def sanity_check(dx_km: float, g_reduced: float, r_drag: float, nu_scale: float):
    H, mask, open_mask, dx, nx, ny, z = build_grid(dx_km)
    f_cor = np.zeros((ny, nx))
    c = np.sqrt(g_reduced * np.median(H[mask]))
    dt = 0.4 * dx / c
    nu = nu_scale * dx ** 2 / dt
    print("-- sanity check: unforced pulse, vertical channel --")
    print(f"  grid {nx}x{ny} at dz={dx_km:g}km (z={z[0]:.1f}-{z[-1]:.1f}km), "
          f"c={c:.3f} m/s, dt={dt/3600:.2f} hr (CFL Courant={c*dt/dx:.2f})")

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
    print(f"  peak max|zeta| over the whole run = {max_abs:.4f} "
          f"(started at 1.0 pulse -- no interim blow-up)")
    return np.isfinite(max_abs) and max_abs < 10.0


def run(years: float, dx_km: float, drho_rho: float, r_drag: float,
       nu_scale: float, tide_amp: float, idx: str, descent_km_per_yr: float,
       yl: float, outdir: Path | None) -> None:
    H, mask, open_mask, dx, nx, ny, z = build_grid(dx_km)
    f_cor = np.zeros((ny, nx))
    g_reduced = G * drho_rho
    c = np.sqrt(g_reduced * np.median(H[mask]))
    dt = 0.4 * dx / c
    nu = nu_scale * dx ** 2 / dt

    dat = np.loadtxt(ROOT / idx / f"{idx}.dat")
    dates_obs, obs = dat[:, 0].copy(), dat[:, 1].copy()
    t0_year = dates_obs[0]
    # forcing built on the REAL monthly dates_obs grid (see draconic_forcing
    # docstring: Impulse_Delta needs native monthly cadence, not a finer
    # synthetic grid), then interpolated onto the solver's own dt below.
    forcing, yl_used = draconic_forcing(idx, dates_obs, yl)
    t_years_grid = dates_obs - t0_year

    def forcing_fn(t_sec: float) -> float:
        t_yr = t_sec / (365.2422 * 86400.0)
        return float(np.interp(t_yr, t_years_grid, forcing))

    dt_out_target = 1.0 / 12.0 * 365.2422 * 86400.0
    save_every = max(1, round(dt_out_target / dt))
    n_steps = int(years * 365.2422 * 86400.0 / dt)
    n_steps -= n_steps % save_every

    z30, z50 = p_to_z(30.0), p_to_z(50.0)
    i30 = int(np.argmin(np.abs(z - z30)))
    i50 = int(np.argmin(np.abs(z - z50)))
    j0 = ny // 2
    print(f"-- resolved vertical run: {idx} (wavenumber=0, zero Coriolis) --")
    print(f"  grid {nx}x{ny} at dz={dx_km:g}km, z={z[0]:.1f}-{z[-1]:.1f}km, "
          f"c={c:.3f} m/s")
    print(f"  30hPa -> z={z30:.2f}km (grid i={i30}), 50hPa -> z={z50:.2f}km "
          f"(grid i={i50}), separation={abs(z[i30]-z[i50]):.2f}km")
    print(f"  dt={dt/3600:.2f} hr, {n_steps} steps ({years:g} sim-years), "
          f"forcing=draconic+comb+IIR (yl={yl_used:.4f}d), r_drag={r_drag:g}")

    t0 = time.time()
    _, max_abs, _, row_hist = run_solver(
        H, mask, open_mask, f_cor, dx, g_reduced, r_drag, nu, tide_amp,
        forcing_fn, n_steps, dt, save_every, j0, i30, ny, nx,
        tilt_field=z_norm_field(mask), save_row=True)
    elapsed = time.time() - t0
    print(f"  done in {elapsed:.1f}s wall-clock; max|zeta| reached {max_abs:.4f}")

    gauge_dates = t0_year + np.arange(len(row_hist)) / 12.0
    gauge30 = row_hist[:, i30]
    gauge50 = row_hist[:, i50]

    obs30 = np.interp(gauge_dates, dates_obs, obs) if idx == "qbo30" else None
    dat_other = "qbo50" if idx == "qbo30" else "qbo30"
    dat_other_arr = np.loadtxt(ROOT / dat_other / f"{dat_other}.dat")
    obs_other = np.interp(gauge_dates, dat_other_arr[:, 0], dat_other_arr[:, 1])
    obs_this = np.interp(gauge_dates, dates_obs, obs)

    def best_lag(a, b, max_lag=12):
        a = (a - a.mean()) / a.std()
        b = (b - b.mean()) / b.std()
        best_l, best_r = 0, -2.0
        for lag in range(-max_lag, max_lag + 1):
            r = np.corrcoef(a, np.roll(b, lag))[0, 1]
            if r > best_r:
                best_r, best_l = r, lag
        return best_l, best_r

    print(f"\n  FIXED gauges at 30hPa and 50hPa altitude, cc vs each other:")
    l, r = best_lag(gauge30, gauge50)
    print(f"    sim(30hPa) vs sim(50hPa): best lag={l:+d} mo, r={r:+.3f}  "
          f"(real data: lag=-4mo, r=0.802)")

    print(f"  FIXED gauges vs REAL {idx}.dat:")
    r_direct = np.corrcoef(gauge30, obs_this)[0, 1] if idx == "qbo30" else \
        np.corrcoef(gauge50, obs_this)[0, 1]
    print(f"    sim({'30' if idx=='qbo30' else '50'}hPa, matching level) vs "
          f"real {idx}: r={r_direct:+.3f}")

    # descending moving gauge: position shifts DOWNWARD (toward higher
    # pressure/lower altitude) linearly with time, at a tunable rate
    z_cols = z
    x_cols = np.arange(nx)
    t_yr_since_start = (gauge_dates - gauge_dates[0])
    for rate in (descent_km_per_yr * 0.5, descent_km_per_yr, descent_km_per_yr * 2.0):
        z_traj = z[i30] - rate * t_yr_since_start  # start at 30hPa, descend
        i_traj = np.interp(z_traj, z_cols[::-1], x_cols[::-1])
        i_traj = np.clip(i_traj, 1.0, nx - 2.0)
        moving = np.array([np.interp(vi, x_cols, row_hist[t])
                           for t, vi in enumerate(i_traj)])
        l1, r1 = best_lag(moving, obs_this)
        print(f"  DESCENDING gauge (rate={rate:.2f} km/yr, "
              f"{rate/12:.3f} km/mo): best lag vs real {idx}={l1:+d}mo r={r1:+.3f}")

    out_dir = outdir if outdir is not None else ROOT / idx
    plot_result(idx, z, mask.shape, gauge_dates, gauge30, gauge50, obs_this,
               obs_other, dat_other, out_dir / "qbo_resolved_vertical.png")


def plot_result(idx, z, shape, gauge_dates, gauge30, gauge50, obs_this,
                obs_other, dat_other, out_path: Path) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(11, 9))

    def zscore(x):
        return (x - np.mean(x)) / np.std(x)

    axes[0].plot(gauge_dates, zscore(gauge30), label="sim @ 30hPa altitude")
    axes[0].plot(gauge_dates, zscore(gauge50), label="sim @ 50hPa altitude",
                alpha=0.8)
    axes[0].set_title(f"{idx}: resolved vertical model, fixed gauges at "
                     f"30hPa vs 50hPa", fontsize=9)
    axes[0].legend(fontsize=7)
    axes[0].set_xlabel("year")

    this_level = "30" if idx == "qbo30" else "50"
    axes[1].plot(gauge_dates, zscore(obs_this), color="0.5",
                label=f"real {idx} ({this_level}hPa, standardized)")
    sim_this = gauge30 if idx == "qbo30" else gauge50
    axes[1].plot(gauge_dates, zscore(sim_this), color="tab:blue", alpha=0.8,
                label=f"sim @ {this_level}hPa (standardized)")
    axes[1].set_title(f"{idx}: sim at matching level vs real data", fontsize=9)
    axes[1].legend(fontsize=7)
    axes[1].set_xlabel("year")

    axes[2].plot(gauge_dates, zscore(obs_this), color="0.5",
                label=f"real {idx}")
    axes[2].plot(gauge_dates, zscore(obs_other), color="tab:orange",
                label=f"real {dat_other}")
    axes[2].set_title("for reference: the two REAL datasets together",
                     fontsize=9)
    axes[2].legend(fontsize=7)
    axes[2].set_xlabel("year")

    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"\n  saved {out_path}")


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sanity", action="store_true")
    ap.add_argument("--idx", default="qbo30", choices=["qbo30", "qbo50"])
    ap.add_argument("--years", type=float, default=70.0)
    ap.add_argument("--dx-km", type=float, default=0.5)
    ap.add_argument("--drho-rho", type=float, default=0.005)
    ap.add_argument("--r-drag", type=float, default=1e-7)
    ap.add_argument("--nu-scale", type=float, default=0.02)
    ap.add_argument("--tide-amp", type=float, default=2e-2)
    ap.add_argument("--descent-km-per-yr", type=float, default=12.0,
                     help="central descent rate to test, km/yr (default 12 "
                          "= 1 km/month, the real documented QBO rate)")
    ap.add_argument("--yl", type=float, default=365.2422,
                     help="year length (days) fed to tide_sum -- lpap "
                          "amp/phase is reused as-is (AMO convention)")
    ap.add_argument("--outdir", type=Path, default=None)
    args = ap.parse_args()

    if args.sanity:
        g_reduced = G * args.drho_rho
        ok = sanity_check(args.dx_km, g_reduced, args.r_drag, args.nu_scale)
        return 0 if ok else 1

    run(args.years, args.dx_km, args.drho_rho, args.r_drag, args.nu_scale,
       args.tide_amp, args.idx, args.descent_km_per_yr, args.yl, args.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
