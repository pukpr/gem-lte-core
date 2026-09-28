#!/usr/bin/env python3
"""pna_resolved_shallow_water.py — a resolved shallow-water run for the
Pacific-North American (PNA) teleconnection pattern, a genuine test of
whether the basin/waveguide geometry itself can be diagnosed from a
climate index's own winding character, rather than just another index
run through the same template.

Why PNA is a distinct, informative test case
-----------------------------------------------------------------------
- Real record starts ~1950 (n=917, 1950.0-2026.3, confirmed gap-free --
  detect_gaps('pna') returns none, uniform ~0.0833yr monthly steps
  throughout) -- roughly half the ~140yr span used for amo/pdo/nao/
  baltic/brestexcl, so this run is deliberately sized to that shorter
  span rather than padded out to 140 years: about half the timesteps,
  half the wall-clock cost, matching what a shorter record should cost.
  (IDATE=1880 in pna's own lt.exe.resp turns out to be moot here: since
  pna.dat itself only contains rows from 1950 onward, `iir`'s own
  start_date search already resolves to index 0 -- i.e. even the
  production Ada model's own forcing has no pre-1950 spin-up for this
  index, so this script's forcing needs none either. Verified before
  writing this rather than assumed.)
- PNA is unambiguously ATMOSPHERIC (a 500hPa geopotential height
  teleconnection pattern -- the reference NOAA/CPC definition, not an
  ocean or sea-level record), so this reuses the atmospheric equivalent-
  depth convention from nao_resolved_shallow_water.py (c~50m/s scale),
  not the ocean reduced-gravity one used for baltic/amo/pdo/brestexcl.
- Distinctive winding character, already established
  (winding_rank.py --index pna, re-verified live): THREE passing,
  unusually SHARP ridges -- M=1.460 (peak=2.68 bits, fwhm=0.030,
  cont=1.00 -- perfectly stationary), M=4.730 (2.71 bits, fwhm=0.030,
  cont=0.92), M=2.970 (2.31 bits, fwhm=0.030, cont=0.77). All three have
  fwhm=0.030 -- sharper than most other indices surveyed in
  WINDING_SCALOGRAM_FEASIBILITY.md (typically 0.04-0.09). M=1.460 is the
  order-7 harmonic of the shared 0.2076 backbone (0.2076*7=1.453,
  matching this index's own harm=[...,7,...] entry) -- "faster cycling"
  in the user's own words is this: PNA's dominant winding number sits at
  7x the shared backbone rate, well above nao's 0.83 (order ~4) or amo's
  -0.013 (order ~0), meaning PNA locks onto a much higher harmonic of
  the same shared clock than the other mid-latitude indices do.

The actual test: does the basin/waveguide (encoded here as the
equivalent depth -> wave speed -> what winding number a fixed-size basin
naturally locks onto) matter for reproducing that "faster cycling"
character? `--depth` is exposed specifically to sweep this without
adding any new degree of freedom to the MODEL itself -- it's the same
parameter every other resolved script already exposes, just tested
deliberately here across several plausible atmospheric-equivalent-depth
values instead of left at one default.

Usage
-----
    ./pna_resolved_shallow_water.py --sanity
    ./pna_resolved_shallow_water.py --years 76
    ./pna_resolved_shallow_water.py --years 76 --depth 150
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
from baltic_resolved_shallow_water import run_solver, G                    # noqa: E402
from amo_resolved_shallow_water import full_fcor                           # noqa: E402
from nao_resolved_shallow_water import build_grid                          # noqa: E402
from pdo_resolved_shallow_water import (build_forcing_at_yl, production_yl,  # noqa: E402
                                        diagonal_tilt_field)
from mode_evidence import sweep_M                                          # noqa: E402
from winding_rank import rank_series                                       # noqa: E402
from qbo_compensating_dtw import standardize, dtw_distance, null_baseline  # noqa: E402
from signal_operators import delayed_difference                            # noqa: E402

ROOT = Path(__file__).resolve().parent
IDX = "pna"
REAL_RIDGES = [1.460, 4.730, 2.970]


def robust_z(model_signed, obs_s, seeds, n_trials=32):
    d_true, _ = dtw_distance(model_signed, obs_s)
    zs = []
    for seed in seeds:
        rng = np.random.default_rng(seed)
        nm, ns = null_baseline(rng, model_signed, obs_s, n_trials=n_trials)
        zs.append((nm - d_true) / ns if ns > 0 else float("nan"))
    return float(np.mean(zs)), float(np.std(zs))


def sanity_check(dx_km, g_eq, depth, r_drag, nu_scale, lat_south, lat_north):
    _, mask, open_mask, dx, nx, ny, lat = build_grid(dx_km, lat_south, lat_north)
    H = np.full((ny, nx), depth)
    f_cor = full_fcor(lat, nx)
    c = np.sqrt(g_eq * depth)
    dt = 0.4 * dx / c
    nu = nu_scale * dx ** 2 / dt
    print("-- sanity check: unforced pulse, PNA (Pacific-North America) "
          "basin --")
    print(f"  grid {nx}x{ny} at dx={dx_km:g}km, equiv depth={depth:.0f}m, "
          f"c={c:.1f} m/s, dt={dt/3600:.2f} hr (CFL Courant={c*dt/dx:.2f})")
    print(f"  latitude range {lat_south:.0f}N-{lat_north:.0f}N")

    j0, i0 = ny // 2, nx // 2
    n_steps = int(20 * 86400 / dt)
    diag_every = max(1, n_steps // 20)
    init_zeta = np.zeros((ny, nx))
    init_zeta[j0, i0] = 1.0
    _, max_abs, diag, _ = run_solver(
        H, mask, open_mask, f_cor, dx, g_eq, r_drag, nu,
        tide_amp=0.0, forcing_fn=lambda t: 0.0, n_steps=n_steps, dt=dt,
        save_every=n_steps, gauge_j=j0, gauge_i=i0, ny=ny, nx=nx,
        diag_every=diag_every, init_zeta=init_zeta)

    print("  day   total energy (should decay smoothly, not blow up)   max|zeta|")
    for day, e, mz in diag:
        print(f"    {day:6.1f}  {e:.4e}  {mz:.4f}")
    print(f"  peak max|zeta| over the whole run = {max_abs:.4f} "
          f"(started at 1.0 pulse -- no interim blow-up)")
    return np.isfinite(max_abs) and max_abs < 10.0


def regenerate(years, dx_km=300.0, g_eq=9.81, depth=300.0, r_drag=1e-7,
              nu_scale=0.02, tide_amp=0.02, lat_south=20.0, lat_north=60.0,
              alpha=0.5, k_shift=0.1, smooth_sigma_months=0.0):
    _, mask, open_mask, dx, nx, ny, lat = build_grid(dx_km, lat_south, lat_north)
    H = np.full((ny, nx), depth)
    f_cor = full_fcor(lat, nx)
    c = np.sqrt(g_eq * depth)
    dt = 0.4 * dx / c
    nu = nu_scale * dx ** 2 / dt

    yl = production_yl(IDX)
    dates, forcing, obs = build_forcing_at_yl(IDX, yl)
    t_years = dates - dates[0]

    if smooth_sigma_months > 0:
        # Mass-conserving (Gaussian-kernel) smoothing of the monthly
        # forcing BEFORE it drives the PDE -- the twice-yearly
        # Impulse_Delta kicks are near-discontinuous at monthly
        # resolution (dF/dt ~80-90/yr at a kick vs ~0.005/yr on a
        # plateau), and an explicit time-stepper's local truncation
        # error is worst exactly where the forcing derivative is largest.
        # Plain linear interpolation (the previous behavior) already
        # smooths the kick somewhat, but ACCIDENTALLY and with a kinked
        # (discontinuous-derivative) result at every monthly node --
        # this instead deliberately spreads the SAME total phase change
        # (Gaussian convolution conserves the integral) over a chosen,
        # PDE-resolvable width, rather than leaving it to whatever
        # interpolation scheme happens to do with it.
        from scipy.ndimage import gaussian_filter1d
        forcing = gaussian_filter1d(forcing, sigma=smooth_sigma_months,
                                    mode="nearest")

    def forcing_fn(t_sec: float) -> float:
        t_yr = t_sec / (365.2422 * 86400.0)
        return float(np.interp(t_yr, t_years, forcing))

    dt_out_target = 1.0 / 12.0 * 365.2422 * 86400.0
    save_every = max(1, round(dt_out_target / dt))
    n_steps = int(years * 365.2422 * 86400.0 / dt)
    n_steps -= n_steps % save_every

    j0, i0 = ny // 2, nx // 2
    tilt = diagonal_tilt_field(mask, alpha)

    print(f"-- resolved shallow-water run: {IDX} (atmospheric, full "
          f"Coriolis) --")
    print(f"  grid {nx}x{ny} at dx={dx_km:g}km, equiv depth={depth:.0f}m, "
          f"c={c:.1f} m/s, lat {lat_south:.0f}N-{lat_north:.0f}N")
    print(f"  dt={dt/3600:.2f} hr, {n_steps} steps ({years:g} sim-years), "
          f"forcing yl={yl:.4f}d, r_drag={r_drag:g}")

    t0 = time.time()
    fixed, max_abs, _, row_hist = run_solver(
        H, mask, open_mask, f_cor, dx, g_eq, r_drag, nu, tide_amp,
        forcing_fn, n_steps, dt, save_every, j0, i0, ny, nx,
        tilt_field=tilt, save_row=True)
    print(f"  done in {time.time()-t0:.1f}s wall-clock; max|zeta| reached "
          f"{max_abs:.4f}")

    gauge_dates = dates[0] + np.arange(len(fixed)) / 12.0
    gauge_forcing = np.interp(gauge_dates, dates, forcing)

    x_cols = np.arange(nx)
    virtual_i = np.clip(i0 + k_shift * gauge_forcing, 1.0, nx - 2.0)
    moving = np.array([np.interp(vi, x_cols, row_hist[t])
                       for t, vi in enumerate(virtual_i)])

    return gauge_dates, gauge_forcing, fixed, moving


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sanity", action="store_true")
    ap.add_argument("--years", type=float, default=76.0)
    ap.add_argument("--dx-km", type=float, default=300.0)
    ap.add_argument("--g-eq", type=float, default=9.81)
    ap.add_argument("--depth", type=float, default=300.0,
                     help="atmospheric equivalent depth, m (the "
                          "waveguide/basin test knob -- sweep this to "
                          "see whether it shifts the emergent winding "
                          "structure toward PNA's own high M=1.46)")
    ap.add_argument("--lat-south", type=float, default=20.0)
    ap.add_argument("--lat-north", type=float, default=60.0)
    ap.add_argument("--r-drag", type=float, default=1e-7)
    ap.add_argument("--nu-scale", type=float, default=0.02)
    ap.add_argument("--tide-amp", type=float, default=0.02)
    ap.add_argument("--alpha", type=float, default=0.5)
    ap.add_argument("--m-max", type=float, default=5.0)
    ap.add_argument("--dm", type=float, default=0.01)
    ap.add_argument("--outdir", type=Path, default=None)
    ap.add_argument("--smooth-sigma-months", type=float, default=0.0,
                     help="mass-conserving Gaussian smoothing sigma (in "
                          "months) applied to the monthly forcing before "
                          "it drives the PDE, to deliberately spread the "
                          "near-discontinuous twice-yearly comb kicks "
                          "over a PDE-resolvable width instead of leaving "
                          "it to accidental linear-interpolation kinks")
    ap.add_argument("--tag", default="", help="suffix for cache/plot "
                     "filenames, e.g. for a depth sweep")
    args = ap.parse_args()

    if args.sanity:
        ok = sanity_check(args.dx_km, args.g_eq, args.depth, args.r_drag,
                          args.nu_scale, args.lat_south, args.lat_north)
        return 0 if ok else 1

    cache = ROOT / IDX / f"pna_resolved_check_cache{args.tag}.npz"
    if cache.exists():
        print(f"Loading cached output from {cache}...")
        d = np.load(cache)
        gauge_dates, gauge_forcing, fixed, moving = (
            d["gauge_dates"], d["gauge_forcing"], d["fixed"], d["moving"])
    else:
        gauge_dates, gauge_forcing, fixed, moving = regenerate(
            args.years, args.dx_km, args.g_eq, args.depth, args.r_drag,
            args.nu_scale, args.tide_amp, args.lat_south, args.lat_north,
            args.alpha, smooth_sigma_months=args.smooth_sigma_months)
        np.savez(cache, gauge_dates=gauge_dates, gauge_forcing=gauge_forcing,
                fixed=fixed, moving=moving)
        print(f"  cached to {cache}")

    dat = np.loadtxt(ROOT / IDX / f"{IDX}.dat")
    obs_on_gauge = np.interp(gauge_dates, dat[:, 0], dat[:, 1])

    diff = delayed_difference(fixed, 12)
    diff_dates = gauge_dates[12:]
    diff_forcing = gauge_forcing[12:]
    diff_obs = obs_on_gauge[12:]

    variants = [
        ("FIXED", fixed, gauge_dates, gauge_forcing, obs_on_gauge),
        ("MOVING (k_shift=0.1)", moving, gauge_dates, gauge_forcing, obs_on_gauge),
        ("FIXED, 12mo delayed-difference", diff, diff_dates, diff_forcing, diff_obs),
    ]
    results = {}
    print("\n================ correlation + DTW vs real pna "
          "================")
    for i, (name, series, dts, forc, obs_v) in enumerate(variants):
        r_raw = float(np.corrcoef(series, obs_v)[0, 1])
        sign = -1.0 if r_raw < 0 else 1.0
        model_s = standardize(series)
        obs_s = standardize(obs_v)
        z_raw, z_raw_std = robust_z(sign * model_s, obs_s,
                                    [10 * i, 10 * i + 1, 10 * i + 2])

        def smooth(x, win):
            k = np.ones(win) / win
            return np.convolve(x, k, mode="same")
        win = 3 * 12
        m3 = standardize(smooth(series, win))
        o3 = standardize(smooth(obs_v, win))
        r3 = float(np.corrcoef(m3, o3)[0, 1])
        sign3 = -1.0 if r3 < 0 else 1.0
        z3, z3_std = robust_z(sign3 * m3, o3, [10 * i, 10 * i + 1, 10 * i + 2])

        print(f"-- {name} vs real {IDX} --")
        print(f"  Pearson r: raw={r_raw:+.3f} sign-matched={sign*r_raw:+.3f}  "
              f"3yr-smoothed sign-matched={sign3*r3:+.3f}")
        print(f"  DTW z-score: raw={z_raw:+.2f}(+/-{z_raw_std:.2f})  "
              f"3yr-smoothed={z3:+.2f}(+/-{z3_std:.2f})")
        results[name] = dict(r_raw=r_raw, sign=sign, z_raw=z_raw, z3=z3,
                             r3=sign3 * r3)

    print(f"\n================ ridge scan vs {IDX}'s OWN Forcing "
          f"(m_max={args.m_max}, testing for the 3 established real "
          f"ridges: {REAL_RIDGES}) ================")
    ranked = {}
    for name, series, dts, forc, obs_v in variants:
        r = rank_series(dts, series, forc, m_max=args.m_max, dm=args.dm,
                        sigma=15.0, t0_step=5.0)
        ranked[name] = r
        print(f"  [{name}] {len(r)} candidate(s):")
        for row in r[:6]:
            near = min(REAL_RIDGES, key=lambda m: abs(m - row["M"]))
            flag = (f"  <-- near real M={near:.3f}"
                    if abs(near - row["M"]) < 0.05 else "")
            print(f"    M={row['M']:.3f}  peak={row['peak_bits']:.2f}bits  "
                  f"fwhm={row['fwhm_M']:.3f}  cont={row['continuity']:.2f}  "
                  f"pass={row['pass']}{flag}")
        if not r:
            print("    (no candidates reached the peak-power screen)")

    print(f"\n(reference: real {IDX}'s own winding_rank ridges (M/peak_bits/"
          f"cont): 1.460/2.68/1.00, 4.730/2.71/0.92, 2.970/2.31/0.77 -- "
          f"all PASS)")

    plot_result(gauge_dates, fixed, moving, variants, results, args.m_max,
               args.dm,
               (args.outdir if args.outdir is not None else ROOT / IDX) /
               f"pna_resolved_shallow_water{args.tag}.png")
    return 0


def plot_result(gauge_dates, fixed, moving, variants, results, m_max, dm,
                out_path: Path) -> None:
    n_panels = 2 + len(variants)
    fig, axes = plt.subplots(n_panels, 1, figsize=(11, 3.1 * n_panels))

    axes[0].plot(gauge_dates, fixed, linewidth=0.5, label="FIXED")
    axes[0].plot(gauge_dates, moving, linewidth=0.5, alpha=0.8, label="MOVING")
    axes[0].set_title(f"{IDX}: raw (undetrended) gauges", fontsize=9)
    axes[0].legend(fontsize=7)
    axes[0].set_xlabel("year")

    def zscore(x):
        return (x - np.mean(x)) / np.std(x)

    for ax, (name, series, dts, forc, obs_v) in zip(
            axes[1:1 + len(variants)], variants):
        res = results[name]
        sign = res["sign"]
        ax.plot(dts, zscore(obs_v), color="0.55", linewidth=0.8,
               label=f"real {IDX}")
        ax.plot(dts, sign * zscore(series), color="tab:blue", linewidth=0.9,
               alpha=0.85, label=f"{name} (r={sign*res['r_raw']:+.3f}, "
                                  f"3yr r={res['r3']:+.3f}, "
                                  f"DTW z={res['z_raw']:+.2f}/"
                                  f"{res['z3']:+.2f})")
        ax.axhline(0, color="k", linewidth=0.4)
        ax.set_title(f"{IDX}: {name}", fontsize=9)
        ax.legend(fontsize=7)
        ax.set_xlabel("year")

    ax = axes[-1]
    m_grid = np.arange(dm, m_max + dm / 2, dm)
    colors = ["tab:blue", "tab:green", "tab:orange"]
    for (name, series, dts, forc, obs_v), color in zip(variants, colors):
        Amat = np.column_stack([np.ones_like(forc), forc])
        coef, *_ = np.linalg.lstsq(Amat, series, rcond=None)
        resid = series - Amat @ coef
        m_grid_out, curve = sweep_M(forc, resid, (dm, m_max), len(m_grid))
        ax.plot(m_grid_out, np.abs(curve), color=color, alpha=0.8, label=name,
               linewidth=0.8)
    for m in REAL_RIDGES:
        ax.axvline(m, color="tab:red", linestyle=":", alpha=0.5)
    ax.set_xlabel("winding number M")
    ax.set_ylabel("|ridge scan cc|")
    ax.set_title(f"{IDX}: ridge scan (red dotted = the 3 established "
                f"real ridges)", fontsize=9)
    ax.legend(fontsize=7)

    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"\n  saved {out_path}")


if __name__ == "__main__":
    sys.exit(main())
