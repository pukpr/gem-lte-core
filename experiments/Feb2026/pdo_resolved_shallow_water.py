#!/usr/bin/env python3
"""pdo_resolved_shallow_water.py — a resolved shallow-water run for PDO,
deliberately positioned BETWEEN the two existing resolved regimes rather
than reusing either one wholesale.

Why "between ENSO and AMO" is a real physical statement, not just a
convenient label
-----------------------------------------------------------------------
- nino4 (`nino4_resolved_shallow_water.py` / `nino4_beta_plane.py`):
  equatorial, zero/beta-plane Coriolis, ZONAL (east-west) thermocline
  tilt -- PDO's own dominant mode is well known to correlate with ENSO
  (the "atmospheric bridge"/tropical-extratropical teleconnection), and
  a large fraction of PDO's variance is literally ENSO re-expressed at
  North Pacific gyre scale.
- amo (`amo_resolved_shallow_water.py`): full mid-latitude Coriolis
  (f=2*Omega*sin(lat) across a wide 10N-65N span), MERIDIONAL
  (north-south) basin-tilt, driven by AMOC's own meridional character.
- PDO's real spatial fingerprint (Mantua et al. 1997's leading North
  Pacific SST EOF) is neither purely zonal nor purely meridional: a
  "horseshoe" pattern, one sign in the central/western North Pacific
  bordered by the opposite sign along the eastern rim -- a genuinely
  DIAGONAL/basin-perimeter structure. The domain here reuses AMO's own
  full-Coriolis machinery (a mid-latitude basin this wide cannot use a
  single-reference-latitude beta-plane, the same reasoning that ruled it
  out for AMO), but narrower and more equatorward (20N-60N, matching the
  real PDO index's own EOF domain -- computed poleward of 20N precisely
  so it isn't dominated by the raw tropical ENSO signal, while still
  being close enough to the tropics that Rossby-wave communication from
  the ENSO band is physically live) and forced through a TILT FIELD that
  is a genuine mix of the zonal (ENSO) and meridional (AMO) patterns --
  literally the two existing tilt fields blended, `alpha` controlling
  the mix (default 0.5, "half of each").

Forcing: the REAL, full 42-term lpap -> tide_sum -> Impulse_Delta -> IIR
-> Bessel pipeline (matching build_forcing_at_yl's own established
pattern in amo_resolved_shallow_water.py) -- unlike QBO, PDO has no
wavenumber=0 simplification, so there is no reason to isolate a single
term here. PDO's own fitted top winding number (`ltep[2]=0.207278`,
NM=3's Bessel-modulation frequency) lands almost exactly on the
0.2075 backbone this project has already found shared across
baltic/nino4/amo -- PDO is not an outlier calibration, it sits on the
same cross-index winding family.

Moving gauge: per the user's own framing, this is placed here
specifically because PDO's behavior is expected to need it MORE than a
single-tilt-axis basin would -- reuses the same material-coordinate
mechanism as AMO's script (a virtual sampling position that wobbles
along the saved row in proportion to the forcing itself), but riding on
top of the diagonal tilt field rather than a purely meridional one.

Usage
-----
    ./pdo_resolved_shallow_water.py --sanity
    ./pdo_resolved_shallow_water.py --years 140
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
from nino4_resolved_shallow_water import x_norm_field                      # noqa: E402
from amo_resolved_shallow_water import full_fcor                           # noqa: E402
from mode_evidence import sweep_M                                          # noqa: E402
from lte_forward import (tide_sum, impulse_delta, iir, bessel, read_resp,  # noqa: E402
                         year_length as lte_year_length)

ROOT = Path(__file__).resolve().parent


def build_forcing_at_yl(idx: str, yl: float):
    """Same pipeline as amo_resolved_shallow_water.build_forcing_at_yl
    (full 42-term lpap, real Impulse_Delta+IIR+Bessel) -- PDO has no
    wavenumber-0 shortcut, so the full pipeline is the right one, not a
    simplification."""
    idx_dir = ROOT / idx
    params = json.loads((idx_dir / "lt.exe.p").read_text())
    lpap = np.array(params["lpap"], dtype=float)
    periods = np.abs(lpap[:, 0])
    B = {k: float(params.get(k, 0.0)) for k in
         ("offs", "bg", "impA", "impB", "delA", "delB", "asym", "ma", "mp",
          "shfT", "init")}
    ltep = params["ltep"]
    k1 = ltep[len(ltep) - 1]     # matches Dipole_Model's m[nm-1] Bessel arg
    dat = np.loadtxt(idx_dir / f"{idx}.dat")
    dates, obs = dat[:, 0].copy(), dat[:, 1].copy()
    tf = tide_sum(dates, lpap[:, 1:3], periods, yl, 0.0, B["shfT"])
    comb = impulse_delta(dates, B["delA"], B["delB"], B["asym"], 12)
    R = iir(tf * comb, lag_a=1.0 - B["ma"], lag_c=B["mp"], init=B["init"],
           start_date=1880.0, dates=dates)
    forcing = bessel(R, B["impA"], B["impB"], k1, B["offs"], B["bg"])
    return dates, forcing, obs


def production_yl(idx: str) -> float:
    """This index's own real production year length (YEAR_IN_DAYS + resp
    YEAR + lt.exe.p 'year' candidate) -- unlike AMO's script, PDO has no
    prior comb-alias sweep establishing a special yl, so the honest
    default is the real calibrated value, not a borrowed constant.

    BUG FIX (2026-09-17): this used to parse lt.exe.resp with a naive
    `.read_text().split()`, which breaks on ANY quoted multi-word value
    appearing before the YEAR line (e.g. `NH "1 1 1 1 1"`, present in
    EVERY index's resp file this project has) -- the quotes get split
    apart into separate whitespace tokens, desynchronizing the assumed
    strict name/value alternation for everything after that point, so
    YEAR's real value was silently never found and resp_year silently
    stayed 0.0 for every index checked (baltic/nao/brestexcl/pdo/pna/
    qbo30/qbo50 all confirmed affected). Fixed by reusing lte_forward.
    read_resp, the same shlex-based (quote-aware) parser already used
    and validated for the real production forward() pipeline, instead of
    a second, buggy, ad hoc parser."""
    idx_dir = ROOT / idx
    params = json.loads((idx_dir / "lt.exe.p").read_text())
    resp = read_resp(idx_dir / "lt.exe.resp")
    resp_year = float(resp.get("YEAR", 0.0))
    return lte_year_length(resp_year, float(params.get("year", 0.0)))


def build_grid(dx_km: float = 300.0, lat_south: float = 20.0,
               lat_north: float = 60.0):
    """North Pacific-scale basin, narrower and more equatorward than
    AMO's Atlantic domain -- matches the real PDO index's own EOF domain
    (poleward of 20N), wide enough (9000km) for the basin's real
    east-west extent."""
    nx, ny = 30, 15   # ~9000km x 4500km at dx=300km
    dx = dx_km * 1e3
    H = np.full((ny, nx), np.nan)
    mask = np.ones((ny, nx), dtype=bool)
    mask[:, 0] = False
    mask[:, -1] = False
    mask[0, :] = False
    mask[-1, :] = False
    open_mask = np.zeros((ny, nx), dtype=bool)
    lat = lat_south + (lat_north - lat_south) * (np.arange(ny) / (ny - 1))
    return H, mask, open_mask, dx, nx, ny, lat


def diagonal_tilt_field(mask: np.ndarray, alpha: float = 0.5) -> np.ndarray:
    """Blend of the zonal (ENSO-style) and meridional (AMO-style) tilt
    fields -- `alpha` is the zonal fraction (0=pure AMO meridional tilt,
    1=pure ENSO zonal tilt, 0.5=equal mix). This is the concrete
    implementation of "PDO sits between ENSO and AMO": the real PDO
    horseshoe pattern (Mantua et al. 1997) is neither purely zonal nor
    purely meridional."""
    return alpha * x_norm_field(mask) + (1.0 - alpha) * y_norm_field(mask)


def sanity_check(dx_km, g_reduced, depth, r_drag, nu_scale, lat_south,
                 lat_north):
    _, mask, open_mask, dx, nx, ny, lat = build_grid(dx_km, lat_south, lat_north)
    H = np.full((ny, nx), depth)
    f_cor = full_fcor(lat, nx)
    c = np.sqrt(g_reduced * depth)
    dt = 0.4 * dx / c
    nu = nu_scale * dx ** 2 / dt
    print("-- sanity check: unforced pulse, North Pacific-scale basin --")
    print(f"  grid {nx}x{ny} at dx={dx_km:g}km, depth={depth:.0f}m, "
          f"c={c:.3f} m/s, dt={dt/3600:.2f} hr (CFL Courant={c*dt/dx:.2f})")
    print(f"  latitude range {lat_south:.0f}N-{lat_north:.0f}N, "
          f"f_cor range: {f_cor.min():.2e} to {f_cor.max():.2e} rad/s")

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
       idx: str, yl: float | None, alpha: float, m_max: float, dm: float,
       outdir: Path | None) -> None:
    _, mask, open_mask, dx, nx, ny, lat = build_grid(dx_km, lat_south, lat_north)
    H = np.full((ny, nx), depth)
    f_cor = full_fcor(lat, nx)
    g_reduced = G * drho_rho
    c = np.sqrt(g_reduced * depth)
    dt = 0.4 * dx / c
    nu = nu_scale * dx ** 2 / dt

    if yl is None:
        yl = production_yl(idx)
    dates, forcing, obs = build_forcing_at_yl(idx, yl)
    print(f"  forcing manifold built at year length yl={yl:.4f}d "
          f"(this index's own real production calibration)")
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
    print(f"-- resolved shallow-water run: {idx} (full Coriolis, "
          f"diagonal alpha={alpha:.2f} zonal/meridional tilt mix) --")
    print(f"  grid {nx}x{ny} at dx={dx_km:g}km, depth={depth:.0f}m, c={c:.3f} m/s "
          f"(g'={g_reduced:.4f}), lat {lat_south:.0f}N-{lat_north:.0f}N")
    print(f"  f_cor range: {f_cor.min():.2e} to {f_cor.max():.2e} rad/s")
    print(f"  dt={dt/3600:.2f} hr, {n_steps} steps ({years:g} sim-years), "
          f"r_drag={r_drag:g}, nu_scale={nu_scale:g}")

    t0 = time.time()
    gauge, max_abs, _, row_hist = run_solver(
        H, mask, open_mask, f_cor, dx, g_reduced, r_drag, nu, tide_amp,
        forcing_fn, n_steps, dt, save_every, gauge_j, gauge_i, ny, nx,
        tilt_field=tilt, save_row=True)
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

    def best_interior_ridge(m_grid, curve, floor=0.05):
        """argmax over the FULL grid is dominated by an M->0 boundary
        artifact (the scan's own edge, not a genuine winding resonance --
        visible in the ridge-scan plot as a spike that only decays away
        from M=0.002), so report the best peak away from that boundary
        separately from the raw global argmax."""
        mask = m_grid >= floor
        i = int(np.argmax(np.abs(curve[mask])))
        idx_full = np.flatnonzero(mask)[i]
        return m_grid[idx_full], curve[idx_full]

    m_grid = np.arange(dm, m_max + dm / 2, dm)
    m_grid_out, curve = sweep_M(gauge_forcing, gauge_resid, (dm, m_max), len(m_grid))
    i_star = int(np.argmax(np.abs(curve)))
    m_int, cc_int = best_interior_ridge(m_grid_out, curve)
    print(f"  [FIXED gauge] ridge scan (on residual): global argmax "
          f"M_hat={m_grid_out[i_star]:.3f} cc={curve[i_star]:+.3f} "
          f"(likely M->0 boundary artifact, not a real ridge) -- best "
          f"INTERIOR (M>=0.05) peak: M_hat={m_int:.3f}  cc={cc_int:+.3f}")

    x_cols = np.arange(nx)
    moving_for_plot, m_grid_mv_plot, curve_mv_plot = None, None, None
    best_r, best_k, best_moving = 0.0, None, None
    obs_on_gauge = np.interp(gauge_dates, dates, obs)
    for k_shift in (0.05, 0.1, 0.2, 0.3):
        virtual_i = np.clip(gauge_i + k_shift * gauge_forcing, 1.0, nx - 2.0)
        moving = np.array([np.interp(vi, x_cols, row_hist[t])
                           for t, vi in enumerate(virtual_i)])
        coef_mv, _, _, _ = np.linalg.lstsq(A, moving, rcond=None)
        mv_r2 = 1.0 - np.var(moving - A @ coef_mv) / np.var(moving)
        moving_resid = moving - A @ coef_mv
        m_grid_mv, curve_mv = sweep_M(gauge_forcing, moving_resid, (dm, m_max), len(m_grid))
        m_int_mv, cc_int_mv = best_interior_ridge(m_grid_mv, curve_mv)
        r_mv = abs(float(np.corrcoef(moving, obs_on_gauge)[0, 1]))
        print(f"  [MOVING gauge, k_shift={k_shift:.2f}] linear fit R^2={mv_r2:.3f} "
              f"-- best INTERIOR (M>=0.05) ridge: M_hat={m_int_mv:.3f}  "
              f"cc={cc_int_mv:+.3f} -- |r vs real PDO|={r_mv:.3f}")
        report_periods(f"{idx} (moving k={k_shift:.2f})", gauge_dates, moving_resid)
        if k_shift == 0.1:
            moving_for_plot = moving_resid
            m_grid_mv_plot, curve_mv_plot = m_grid_mv, curve_mv
        if r_mv > best_r:
            best_r, best_k, best_moving = r_mv, k_shift, moving

    print(f"  (for comparison: PDO's own fitted top winding is M=0.2073, "
          f"landing on this project's cross-index 0.2075 backbone shared "
          f"with baltic/nino4/amo)")
    print(f"  best moving-gauge k_shift by |r| vs real PDO: {best_k:.2f} "
          f"(|r|={best_r:.3f})")

    report_periods(idx, gauge_dates, gauge_resid)

    out_dir = outdir if outdir is not None else ROOT / idx
    plot_result(idx, H, mask, gauge_dates, gauge, gauge_resid, moving_for_plot,
               best_moving, obs_on_gauge, m_grid_out, curve, m_grid_mv_plot,
               curve_mv_plot, alpha, out_dir / "pdo_resolved_shallow_water.png")


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
    # PDO's own literature range: a shorter ENSO-adjacent band (2-7yr) and
    # a genuinely decadal band (10-60yr) -- report both, unlike AMO's
    # single 20-200yr window, since "between ENSO and AMO" means BOTH
    # timescales are plausible here.
    for lo, hi, tag in ((2, 7, "ENSO-band"), (10, 60, "decadal")):
        band = (period > lo) & (period < hi)
        if not np.any(band):
            continue
        top = np.argsort(-power[band])[:3]
        print(f"  top {tag} ({lo}-{hi}yr) periods: " +
              ", ".join(f"{p:.1f}yr" for p in period[band][top]))


def plot_result(idx, H, mask, gauge_dates, gauge, gauge_resid, moving,
                best_moving, obs, m_grid, curve, m_grid_mv, curve_mv, alpha,
                out_path: Path) -> None:
    fig, axes = plt.subplots(6, 1, figsize=(11, 15.5))
    depth_show = np.where(mask, H, np.nan)
    im = axes[0].imshow(depth_show, origin="lower", cmap="Blues", aspect="auto")
    axes[0].set_title(f"idealized North Pacific-scale basin (uniform depth, "
                     f"full Coriolis, diagonal tilt alpha={alpha:.2f})",
                     fontsize=9)
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

    r_raw = float(np.corrcoef(best_moving, obs)[0, 1])
    sign = -1.0 if r_raw < 0 else 1.0
    r_overlay = sign * r_raw
    axes[2].plot(gauge_dates, zscore(obs), color="0.55", linewidth=0.8,
                label="real PDO (standardized)")
    axes[2].plot(gauge_dates, sign * zscore(best_moving), color="tab:green",
                linewidth=1.0, alpha=0.85,
                label=f"best moving gauge, sign-matched (r={r_overlay:+.3f})")
    axes[2].axhline(0, color="k", linewidth=0.4)
    axes[2].set_xlabel("year")
    axes[2].set_ylabel("standardized units")
    axes[2].set_title(f"{idx}: best MOVING gauge overlaid on real PDO data "
                     f"(both z-scored, sign-matched)", fontsize=9)
    axes[2].legend(fontsize=7)

    axes[3].plot(gauge_dates, gauge_resid, linewidth=0.6, color="tab:orange")
    axes[3].set_xlabel("year")
    axes[3].set_ylabel("residual (m)")
    axes[3].set_title(f"{idx}: FIXED-gauge residual after removing linear "
                     f"Forcing tracking", fontsize=9)

    axes[4].plot(m_grid, np.abs(curve))
    axes[4].axvline(0.2073, color="tab:red", linestyle=":", alpha=0.7)
    axes[4].text(0.2073, axes[4].get_ylim()[1]*0.9, "0.2073", color="tab:red",
               fontsize=7)
    axes[4].set_xlabel("winding number M")
    axes[4].set_ylabel("|ridge scan cc|")
    axes[4].set_title(f"{idx}: FIXED-gauge ridge scan on residual "
                     f"(red = PDO's own fitted M)", fontsize=9)

    axes[5].plot(m_grid_mv, np.abs(curve_mv), color="tab:green")
    axes[5].axvline(0.2073, color="tab:red", linestyle=":", alpha=0.7)
    axes[5].text(0.2073, axes[5].get_ylim()[1]*0.9, "0.2073", color="tab:red",
               fontsize=7)
    axes[5].set_xlabel("winding number M")
    axes[5].set_ylabel("|ridge scan cc|")
    axes[5].set_title(f"{idx}: MOVING-gauge ridge scan on residual "
                     f"(k_shift=0.1)", fontsize=9)

    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"\n  saved {out_path}")
    print(f"  overlay correlation (best moving gauge vs real PDO, "
          f"both standardized): r={r_overlay:+.3f}")


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sanity", action="store_true")
    ap.add_argument("--idx", default="pdo")
    ap.add_argument("--years", type=float, default=140.0)
    ap.add_argument("--dx-km", type=float, default=300.0)
    ap.add_argument("--drho-rho", type=float, default=0.005)
    ap.add_argument("--depth", type=float, default=700.0)
    ap.add_argument("--lat-south", type=float, default=20.0)
    ap.add_argument("--lat-north", type=float, default=60.0)
    ap.add_argument("--r-drag", type=float, default=1e-7)
    ap.add_argument("--nu-scale", type=float, default=0.02)
    ap.add_argument("--tide-amp", type=float, default=2e-2)
    ap.add_argument("--yl", type=float, default=None,
                     help="year length (days); default is this index's "
                          "own real production calibration")
    ap.add_argument("--alpha", type=float, default=0.5,
                     help="tilt mix: 0=pure AMO meridional, 1=pure ENSO "
                          "zonal, 0.5=equal blend (default)")
    ap.add_argument("--m-max", type=float, default=5.0)
    ap.add_argument("--dm", type=float, default=0.002)
    ap.add_argument("--outdir", type=Path, default=None)
    args = ap.parse_args()

    if args.sanity:
        g_reduced = G * args.drho_rho
        ok = sanity_check(args.dx_km, g_reduced, args.depth, args.r_drag,
                          args.nu_scale, args.lat_south, args.lat_north)
        return 0 if ok else 1

    run(args.years, args.dx_km, args.drho_rho, args.depth, args.r_drag,
       args.nu_scale, args.tide_amp, args.lat_south, args.lat_north,
       args.idx, args.yl, args.alpha, args.m_max, args.dm, args.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
