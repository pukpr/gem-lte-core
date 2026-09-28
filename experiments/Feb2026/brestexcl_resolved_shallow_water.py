#!/usr/bin/env python3
"""brestexcl_resolved_shallow_water.py — a resolved shallow-water run for
Brest (brestexcl), a real open-coast Atlantic tide-gauge record --
physically distinct from Baltic (nearly tideless, enclosed, wind/
pressure-driven) in the one respect that matters most for this model:
Brest sits on a genuinely macro-tidal coast (English Channel funnel
effect, spring range order several meters, among the largest in Europe)
where the astronomical tide is the DOMINANT physical signal, not a
weak proxy correlation the way it is for AMO/PDO/NAO/Baltic's own MSL
records. So unlike those, an ocean/reduced-gravity treatment here is
modeling something close to the real, literal mechanism, not an analogy.

Grounded first via winding_rank.py --index brestexcl (already-gap-aware
tooling, no new code needed for this step): FIVE independently-fitted,
PASSING, stationary ridges -- M=1.870 (peak=3.55 bits, cont=1.00),
M=0.955 (2.69 bits, 0.96), M=2.490 (2.62 bits, 0.81), M=0.410 (2.17 bits,
0.73), M=0.060 (2.83 bits, 0.73) -- a much richer harmonic ladder than
Baltic's single M=1.25 ridge, plausibly because strong real macro-tidal
forcing here excites genuine harmonic/compound-tide structure (shallow-
water tidal distortion) more readily than Baltic's weak, indirect
forcing does. Three of the five land within ~1-2% of small-integer
multiples of this index's own fitted backbone winding (ltep[4]=0.20759,
itself on this project's shared cross-index 0.2075 backbone):
0.410/0.20759=1.975 (~2x), 1.870/0.20759=9.008 (~9x),
2.490/0.20759=11.995 (~12x). winding_scalogram.png confirms these bands
persist continuously on BOTH sides of the real 1944.3-1954.3 data gap
(WWII-era blackout, already handled transparently by winding_scalogram.py's
own detect_gaps/valid_mask_from_gaps) -- not a gap artifact.

Domain: a modest Bay-of-Biscay/Channel-approach-scale idealized
rectangular basin (not Baltic's bespoke bathymetry-matched shape --
Brest's real dynamics involve continental-shelf/Channel-funnel geometry
this project has no bathymetry for, so this is an explicitly simplified
idealized basin, not a literal Channel model), centered on Brest's real
latitude (48.4N), full Coriolis, OCEAN reduced-gravity (same g'/depth
convention as amo/pdo, whose coarser-grid stability fix this reuses).

Gap-handling for the real-data comparison: brestexcl.dat's own
1944.3-1954.3 window is fabricated linear-interpolation filler (same gap
winding_scalogram.py already excludes) -- masked out of every
correlation/DTW comparison against real data below. The model's OWN
driving Forcing is unaffected (Impulse_Delta/IIR only need continuous
calendar dates, not real data values, so the isolated tidal manifold
itself has no gap).

Usage
-----
    ./brestexcl_resolved_shallow_water.py --sanity
    ./brestexcl_resolved_shallow_water.py --years 140
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
import json                                                                 # noqa: E402

from baltic_resolved_shallow_water import run_solver, G                    # noqa: E402
from amo_resolved_shallow_water import full_fcor                           # noqa: E402
from nao_resolved_shallow_water import build_grid                          # noqa: E402
from pdo_resolved_shallow_water import production_yl, diagonal_tilt_field  # noqa: E402
from mode_evidence import sweep_M                                          # noqa: E402
from winding_scalogram import detect_gaps                                  # noqa: E402
from winding_rank import rank_series                                       # noqa: E402
from qbo_compensating_dtw import standardize, dtw_distance, null_baseline  # noqa: E402
from lte_forward import tide_sum, impulse_delta, iir, bessel                # noqa: E402
from signal_operators import delayed_difference                            # noqa: E402

ROOT = Path(__file__).resolve().parent
IDX = "brestexcl"
BREST_LAT = 48.4

# established ridges, for reference throughout (winding_rank.py --index
# brestexcl, re-verified live at the top of this module's docstring)
REAL_RIDGES = [1.870, 0.955, 2.490, 0.410, 0.060]


def robust_z(model_signed, obs_s, seeds, n_trials=32):
    d_true, _ = dtw_distance(model_signed, obs_s)
    zs = []
    for seed in seeds:
        rng = np.random.default_rng(seed)
        nm, ns = null_baseline(rng, model_signed, obs_s, n_trials=n_trials)
        zs.append((nm - d_true) / ns if ns > 0 else float("nan"))
    return float(np.mean(zs)), float(np.std(zs))


def build_forcing_continuous(idx: str, yl: float, dates: np.ndarray) -> np.ndarray:
    """Rebuilds the SAME lpap->tide_sum->Impulse_Delta->IIR->Bessel
    pipeline as pdo_resolved_shallow_water.build_forcing_at_yl, but on an
    explicitly-supplied, GENUINELY CONTINUOUS date array rather than
    reading dates from the index's own .dat file.

    This matters specifically for brestexcl: its real .dat file does not
    contain fabricated linear-filler rows across its 1944.3-1954.3 gap
    (confirmed directly: lte_results.csv has the identical row count and
    the identical missing-date jump) -- it just SKIPS those ~120 monthly
    rows entirely. build_forcing_at_yl's own gauge_dates construction
    elsewhere in this project assumes uniform monthly spacing and then
    np.interp's the forcing onto that assumed grid FROM the real (gappy)
    dates/forcing arrays -- during a decade with no bracketing points
    less than ~10 years apart, np.interp silently draws a smooth linear
    RAMP across the whole missing decade instead of genuine oscillating
    tidal forcing. Caught directly: the resolved model's own gauge output
    was found to go nearly perfectly FLAT for ~15 years, aligned almost
    exactly with the real gap -- not a coincidence, and not a plotting
    artifact (verified against the raw cached array), but this exact
    interpolation bug. Fixed by recomputing the forcing pipeline on a
    synthetic, gap-free monthly `dates` array covering the full run --
    every stage (tide_sum, Impulse_Delta, IIR, Bessel) is a pure function
    of PARAMETERS and DATES, never of the real observed data values, so
    this recomputation is well-posed and requires no real data for the
    gap years at all."""
    idx_dir = ROOT / idx
    params = json.loads((idx_dir / "lt.exe.p").read_text())
    lpap = np.array(params["lpap"], dtype=float)
    periods = np.abs(lpap[:, 0])
    B = {k: float(params.get(k, 0.0)) for k in
         ("offs", "bg", "impA", "impB", "delA", "delB", "asym", "ma", "mp",
          "shfT", "init")}
    ltep = params["ltep"]
    k1 = ltep[len(ltep) - 1]
    tf = tide_sum(dates, lpap[:, 1:3], periods, yl, 0.0, B["shfT"])
    comb = impulse_delta(dates, B["delA"], B["delB"], B["asym"], 12)
    R = iir(tf * comb, lag_a=1.0 - B["ma"], lag_c=B["mp"], init=B["init"],
           start_date=1880.0, dates=dates)
    return bessel(R, B["impA"], B["impB"], k1, B["offs"], B["bg"])


def sanity_check(dx_km, g_reduced, depth, r_drag, nu_scale, lat_south, lat_north):
    _, mask, open_mask, dx, nx, ny, lat = build_grid(dx_km, lat_south, lat_north)
    H = np.full((ny, nx), depth)
    f_cor = full_fcor(lat, nx)
    c = np.sqrt(g_reduced * depth)
    dt = 0.4 * dx / c
    nu = nu_scale * dx ** 2 / dt
    print("-- sanity check: unforced pulse, Brest/Bay-of-Biscay-scale "
          "basin --")
    print(f"  grid {nx}x{ny} at dx={dx_km:g}km, depth={depth:.0f}m, "
          f"c={c:.3f} m/s, dt={dt/3600:.2f} hr (CFL Courant={c*dt/dx:.2f})")
    print(f"  latitude range {lat_south:.0f}N-{lat_north:.0f}N")

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


def regenerate(years=140.0, dx_km=150.0, drho_rho=0.005, depth=700.0,
              r_drag=1e-7, nu_scale=0.02, tide_amp=0.02, lat_south=40.0,
              lat_north=55.0, alpha=0.5, k_shift=0.1):
    _, mask, open_mask, dx, nx, ny, lat = build_grid(dx_km, lat_south, lat_north)
    H = np.full((ny, nx), depth)
    f_cor = full_fcor(lat, nx)
    g_reduced = G * drho_rho
    c = np.sqrt(g_reduced * depth)
    dt = 0.4 * dx / c
    nu = nu_scale * dx ** 2 / dt

    yl = production_yl(IDX)
    real_dat = np.loadtxt(ROOT / IDX / f"{IDX}.dat")
    t0_year = real_dat[0, 0]
    dates = t0_year + np.arange(0.0, years + 1.0, 1.0 / 12.0)
    forcing = build_forcing_continuous(IDX, yl, dates)
    t_years = dates - dates[0]

    def forcing_fn(t_sec: float) -> float:
        t_yr = t_sec / (365.2422 * 86400.0)
        return float(np.interp(t_yr, t_years, forcing))

    dt_out_target = 1.0 / 12.0 * 365.2422 * 86400.0
    save_every = max(1, round(dt_out_target / dt))
    n_steps = int(years * 365.2422 * 86400.0 / dt)
    n_steps -= n_steps % save_every

    j0 = int(np.argmin(np.abs(lat - BREST_LAT)))
    i0 = nx // 2
    tilt = diagonal_tilt_field(mask, alpha)

    print(f"-- resolved shallow-water run: {IDX} (full Coriolis, "
          f"OCEAN reduced-gravity) --")
    print(f"  grid {nx}x{ny} at dx={dx_km:g}km, depth={depth:.0f}m, "
          f"c={c:.3f} m/s (g'={g_reduced:.4f}), lat {lat_south:.0f}N-"
          f"{lat_north:.0f}N, gauge row j={j0} (lat={lat[j0]:.1f}N)")
    print(f"  dt={dt/3600:.2f} hr, {n_steps} steps ({years:g} sim-years), "
          f"forcing yl={yl:.4f}d, r_drag={r_drag:g}")

    t0 = time.time()
    fixed, max_abs, _, row_hist = run_solver(
        H, mask, open_mask, f_cor, dx, g_reduced, r_drag, nu, tide_amp,
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
    ap.add_argument("--years", type=float, default=140.0)
    ap.add_argument("--dx-km", type=float, default=150.0)
    ap.add_argument("--drho-rho", type=float, default=0.005)
    ap.add_argument("--depth", type=float, default=700.0)
    ap.add_argument("--lat-south", type=float, default=40.0)
    ap.add_argument("--lat-north", type=float, default=55.0)
    ap.add_argument("--r-drag", type=float, default=1e-7)
    ap.add_argument("--nu-scale", type=float, default=0.02)
    ap.add_argument("--tide-amp", type=float, default=0.02)
    ap.add_argument("--alpha", type=float, default=0.5)
    ap.add_argument("--m-max", type=float, default=3.0)
    ap.add_argument("--dm", type=float, default=0.01)
    ap.add_argument("--outdir", type=Path, default=None)
    args = ap.parse_args()

    if args.sanity:
        g_reduced = G * args.drho_rho
        ok = sanity_check(args.dx_km, g_reduced, args.depth, args.r_drag,
                          args.nu_scale, args.lat_south, args.lat_north)
        return 0 if ok else 1

    cache = ROOT / IDX / "brestexcl_resolved_check_cache.npz"
    if cache.exists():
        print(f"Loading cached output from {cache}...")
        d = np.load(cache)
        gauge_dates, gauge_forcing, fixed, moving = (
            d["gauge_dates"], d["gauge_forcing"], d["fixed"], d["moving"])
    else:
        gauge_dates, gauge_forcing, fixed, moving = regenerate(
            args.years, args.dx_km, args.drho_rho, args.depth, args.r_drag,
            args.nu_scale, args.tide_amp, args.lat_south, args.lat_north,
            args.alpha)
        np.savez(cache, gauge_dates=gauge_dates, gauge_forcing=gauge_forcing,
                fixed=fixed, moving=moving)
        print(f"  cached to {cache}")

    dat = np.loadtxt(ROOT / IDX / f"{IDX}.dat")
    obs_on_gauge = np.interp(gauge_dates, dat[:, 0], dat[:, 1])

    gap_start, gap_end = detect_gaps(IDX)[0]
    valid = (gauge_dates < gap_start) | (gauge_dates > gap_end)
    print(f"\nExcluding fabricated filler across the real gap "
          f"{gap_start:.2f}-{gap_end:.2f} from all real-data comparisons "
          f"below ({(~valid).sum()} of {len(valid)} monthly samples "
          f"excluded)")

    diff = delayed_difference(fixed, 12)
    diff_dates = gauge_dates[12:]
    diff_forcing = gauge_forcing[12:]
    diff_obs = obs_on_gauge[12:]
    diff_valid = valid[12:]

    variants = [
        ("FIXED", fixed, gauge_dates, gauge_forcing, obs_on_gauge, valid),
        ("MOVING (k_shift=0.1)", moving, gauge_dates, gauge_forcing,
         obs_on_gauge, valid),
        ("FIXED, 12mo delayed-difference", diff, diff_dates, diff_forcing,
         diff_obs, diff_valid),
    ]
    results = {}
    print("\n================ correlation + DTW vs real brestexcl "
          "(gap-excluded) ================")
    for i, (name, series, dts, forc, obs_v, val_v) in enumerate(variants):
        s_v, o_v = series[val_v], obs_v[val_v]
        r_raw = float(np.corrcoef(s_v, o_v)[0, 1])
        sign = -1.0 if r_raw < 0 else 1.0
        model_s = standardize(s_v)
        obs_s = standardize(o_v)
        z_raw, z_raw_std = robust_z(sign * model_s, obs_s,
                                    [10 * i, 10 * i + 1, 10 * i + 2])
        print(f"-- {name} vs real {IDX} --")
        print(f"  Pearson r: raw={r_raw:+.3f}  sign-matched={sign*r_raw:+.3f}")
        print(f"  DTW z-score (gap-excluded): {z_raw:+.2f} (+/-{z_raw_std:.2f})")
        results[name] = dict(r_raw=r_raw, sign=sign, z_raw=z_raw)

    print(f"\n================ ridge scan vs {IDX}'s OWN Forcing "
          f"(m_max={args.m_max}, testing for the 5 established real "
          f"ridges) ================")
    ranked = {}
    for name, series, dts, forc, obs_v, val_v in variants:
        r = rank_series(dts, series, forc, m_max=args.m_max,
                        dm=args.dm, sigma=15.0, t0_step=5.0)
        ranked[name] = r
        print(f"  [{name}] {len(r)} candidate(s):")
        for row in r[:6]:
            near = min(REAL_RIDGES, key=lambda m: abs(m - row["M"]))
            flag = (f"  <-- near real M={near:.3f}"
                    if abs(near - row["M"]) < 0.04 else "")
            print(f"    M={row['M']:.3f}  peak={row['peak_bits']:.2f}bits  "
                  f"fwhm={row['fwhm_M']:.3f}  cont={row['continuity']:.2f}  "
                  f"pass={row['pass']}{flag}")
        if not r:
            print("    (no candidates reached the peak-power screen)")

    print(f"\n(reference: real {IDX}'s own winding_rank ridges (M/peak_bits/"
          f"cont): 1.870/3.55/1.00, 0.955/2.69/0.96, 2.490/2.62/0.81, "
          f"0.410/2.17/0.73, 0.060/2.83/0.73 -- all PASS)")

    plot_result(gauge_dates, fixed, moving, variants, results, ranked,
               args.m_max, args.dm,
               (args.outdir if args.outdir is not None else ROOT / IDX) /
               "brestexcl_resolved_shallow_water.png")
    return 0


def plot_result(gauge_dates, fixed, moving, variants, results, ranked,
                m_max, dm, out_path: Path) -> None:
    n_panels = 2 + len(variants)
    fig, axes = plt.subplots(n_panels, 1, figsize=(11, 3.1 * n_panels))

    axes[0].plot(gauge_dates, fixed, linewidth=0.5, label="FIXED")
    axes[0].plot(gauge_dates, moving, linewidth=0.5, alpha=0.8, label="MOVING")
    valid0 = variants[0][5]
    gap_dates = gauge_dates[~valid0]
    if len(gap_dates):
        axes[0].axvspan(gap_dates.min(), gap_dates.max(), color="0.85",
                        label="real-data gap (excluded from comparisons)")
    axes[0].set_title(f"{IDX}: raw (undetrended) gauges", fontsize=9)
    axes[0].legend(fontsize=7)
    axes[0].set_xlabel("year")

    for ax, (name, series, dts, forc, obs_v, val_v) in zip(
            axes[1:1 + len(variants)], variants):
        res = results[name]
        sign = res["sign"]
        # z-scored using gap-excluded stats (matching what corr/DTW
        # actually used), but the MODEL is plotted CONTINUOUSLY across
        # its own full (already gap-free) record; only the OBS line is
        # NaN-broken across the real gap, so matplotlib doesn't draw a
        # misleading straight bridge line across missing years the way
        # plotting both on a gap-excluded x-axis alone did in an earlier
        # version of this plot.
        obs_mu, obs_sigma = obs_v[val_v].mean(), obs_v[val_v].std()
        series_mu, series_sigma = series[val_v].mean(), series[val_v].std()
        obs_z_full = (obs_v - obs_mu) / obs_sigma
        obs_z_full_broken = obs_z_full.copy()
        obs_z_full_broken[~val_v] = np.nan
        series_z_full = sign * (series - series_mu) / series_sigma

        ax.plot(dts, obs_z_full_broken, color="0.55", linewidth=0.7,
               label=f"real {IDX} (standardized, gap-excluded)")
        ax.plot(dts, series_z_full, color="tab:blue", linewidth=0.9,
               alpha=0.85, label=f"{name} (r={sign*res['r_raw']:+.3f}, "
                                  f"DTW z={res['z_raw']:+.2f})")
        ax.axhline(0, color="k", linewidth=0.4)
        ax.set_title(f"{IDX}: {name}", fontsize=9)
        ax.legend(fontsize=7)
        ax.set_xlabel("year")

    ax = axes[-1]
    m_grid = np.arange(dm, m_max + dm / 2, dm)
    colors = ["tab:blue", "tab:green", "tab:orange"]
    for (name, series, dts, forc, obs_v, val_v), color in zip(variants, colors):
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
    ax.set_title(f"{IDX}: ridge scan (red dotted = the 5 established real "
                f"ridges)", fontsize=9)
    ax.legend(fontsize=7)

    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"\n  saved {out_path}")


if __name__ == "__main__":
    sys.exit(main())
