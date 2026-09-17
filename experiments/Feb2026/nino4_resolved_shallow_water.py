#!/usr/bin/env python3
"""nino4_resolved_shallow_water.py — the same resolved shallow-water
experiment as baltic_resolved_shallow_water.py, but for the domain the
GEM_Chap12 derivation was actually built for: the equatorial Pacific
(NINO4), where Coriolis genuinely vanishes (sin(phi)=0 at phi=0), not just
approximately as in the Baltic's mid-latitude setting.

Why this matters numerically, not just physically: the Baltic run's
long-duration instability traced to the interaction between the Coriolis
sub-step and weak damping over hundreds of thousands of timesteps (see
baltic_resolved_shallow_water.py / the [[baltic-resolved-shallow-water]]
memory). At the equator, f_cor=0 is not an approximation to relax
carefully -- it is the exact, literal condition the original derivation
already assumes (Chapter 12's own equatorial reduction: "dropping the
Coriolis terms at the exact equator is legitimate"). Removing rotation
entirely removes that entire failure mode: the exact-rotation sub-step
becomes a no-op (cos(0)=1, sin(0)=0), and the remaining explicit forward-
backward gravity-wave scheme has a long, well-understood track record of
being stable under an ordinary CFL condition with NO delicate rotation-
damping trade-off -- so a full, real-length (144yr) resolved run may
actually be achievable here where it wasn't for the Baltic.

Geometry: unlike the Baltic's multi-basin coastline, the equatorial
Pacific needs none of that -- it's a single zonal channel between two
continental boundaries (~130E to ~80W, ~23,000 km), closed (solid-wall) at
both ends, narrow in latitude (no natural equatorial trapping scale exists
without the beta-effect that a nonzero-but-linearized Coriolis would give
-- consistent with reproducing THIS derivation's own zero-Coriolis
simplification, not the fuller Matsuno equatorial-wave theory that keeps
beta). Forcing is a ZONAL (east-west) tilt of the real nino4 Forcing(t) --
matching the classic wind-forced zonal thermocline-tilt picture of ENSO,
rather than the Baltic run's north-south tilt.

Depth/stratification: a single representative equatorial-Pacific
thermocline depth (not the real, well-known east-west shoaling tilt --
a known simplification, consistent with a first resolved pass) and the
same reduced-gravity convention as the Baltic run.

Usage
-----
    ./nino4_resolved_shallow_water.py --sanity
    ./nino4_resolved_shallow_water.py --years 144
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
from cv_rolling_blocked import prepare                                    # noqa: E402
from mode_evidence import sweep_M                                         # noqa: E402

ROOT = Path(__file__).resolve().parent


def build_grid(dx_km: float = 200.0):
    """Idealized equatorial Pacific channel: solid walls at both zonal
    ends (continental boundaries), narrow meridional band (no natural
    trapping scale without beta -- see module docstring), uniform depth."""
    nx, ny = 110, 9   # ~22,000 km x ~1,600 km at dx=200km
    dx = dx_km * 1e3
    H = np.full((ny, nx), 120.0)  # representative equatorial thermocline depth
    mask = np.ones((ny, nx), dtype=bool)
    mask[:, 0] = False    # western continental wall
    mask[:, -1] = False   # eastern continental wall
    mask[0, :] = False    # meridional channel walls
    mask[-1, :] = False
    open_mask = np.zeros((ny, nx), dtype=bool)  # closed basin: no open boundary
    return H, mask, open_mask, dx, nx, ny


def x_norm_field(mask: np.ndarray) -> np.ndarray:
    """East-west normalized coordinate in [-0.5, 0.5] for the zonal
    tilting tidal body force (matches ENSO's classic zonal-tilt picture,
    unlike the Baltic run's north-south tilt)."""
    ny, nx = mask.shape
    i = np.arange(nx)
    return np.repeat(((i / (nx - 1)) - 0.5)[None, :], ny, axis=0)


def sanity_check(dx_km: float, g_reduced: float, r_drag: float, nu_scale: float):
    H, mask, open_mask, dx, nx, ny = build_grid(dx_km)
    f_cor = np.zeros((ny, nx))  # exactly zero -- the equatorial condition
    c = np.sqrt(g_reduced * np.median(H[mask]))
    dt = 0.4 * dx / c
    nu = nu_scale * dx ** 2 / dt
    print("-- sanity check: unforced pulse, zero Coriolis --")
    print(f"  grid {nx}x{ny} at dx={dx_km:g}km, wave speed c={c:.3f} m/s, "
          f"dt={dt/3600:.2f} hr (CFL Courant={c*dt/dx:.2f})")

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


def run(years: float, dx_km: float, drho_rho: float, h1: float, r_drag: float,
       nu_scale: float, tide_amp: float, idx: str, m_max: float, dm: float,
       outdir: Path | None) -> None:
    H, mask, open_mask, dx, nx, ny = build_grid(dx_km)
    f_cor = np.zeros((ny, nx))
    g_reduced = G * drho_rho
    H_active = np.minimum(H, h1)
    c = np.sqrt(g_reduced * np.median(H_active[mask]))
    dt = 0.4 * dx / c
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

    gauge_j, gauge_i = ny // 2, 40  # within the NINO4 box (west-central basin)
    print(f"-- resolved shallow-water run: {idx} (equatorial, zero Coriolis) --")
    print(f"  grid {nx}x{ny} at dx={dx_km:g}km, active-layer c={c:.3f} m/s "
          f"(g'={g_reduced:.4f}, h1={h1:.0f}m capped depth)")
    print(f"  dt={dt/3600:.2f} hr, {n_steps} steps ({years:g} sim-years), "
          f"r_drag={r_drag:g}, nu_scale={nu_scale:g}")

    t0 = time.time()
    gauge, max_abs, _, _ = run_solver(H_active, mask, open_mask, f_cor, dx, g_reduced,
                                   r_drag, nu, tide_amp, forcing_fn, n_steps, dt,
                                   save_every, gauge_j, gauge_i, ny, nx,
                                   tilt_field=x_norm_field(mask))
    elapsed = time.time() - t0
    print(f"  done in {elapsed:.1f}s wall-clock; max|zeta| reached {max_abs:.4f} m "
          f"(forcing amplitude scale: {tide_amp:.3f} m)")

    gauge_dates = dates[0] + np.arange(len(gauge)) / 12.0
    gauge_forcing = np.interp(gauge_dates, dates, forcing)

    # fit_at_M's own design already has a linear k0*Forcing term, but if
    # the gauge is dominated by quasi-static tracking of Forcing's own
    # slow drift (gauge =~ c*Forcing(t)), that linear term alone soaks up
    # nearly all the variance at EVERY candidate M -- which is exactly
    # what a near-flat, ~0.99-everywhere scan means, not evidence there's
    # no ridge. Regress out the best-fit linear-in-Forcing component first
    # and scan the residual, so a genuine M-specific ridge (if any) isn't
    # masked by the dominant linear response.
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
          f"M*=0.450, per cv_ridge_transfer.py -- though not AR1-significant there)")

    out_dir = outdir if outdir is not None else ROOT / idx
    plot_result(idx, H, mask, gauge_dates, gauge, m_grid_out, curve,
               out_dir / "nino4_resolved_shallow_water.png")


def plot_result(idx, H, mask, gauge_dates, gauge, m_grid, curve, out_path: Path) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(11, 8))
    depth_show = np.where(mask, H, np.nan)
    im = axes[0].imshow(depth_show, origin="lower", cmap="Blues", aspect="auto")
    axes[0].set_title("idealized equatorial Pacific channel (depth, m), "
                     "zero Coriolis", fontsize=9)
    fig.colorbar(im, ax=axes[0], shrink=0.7, orientation="horizontal", pad=0.3)

    axes[1].plot(gauge_dates, gauge, linewidth=0.5)
    axes[1].set_xlabel("year")
    axes[1].set_ylabel("virtual gauge zeta (m)")
    axes[1].set_title(f"{idx}: simulated thermocline displacement", fontsize=9)

    axes[2].plot(m_grid, np.abs(curve))
    axes[2].axvline(0.45, color="tab:red", linestyle=":",
                    label="real nino4's M*=0.450")
    axes[2].set_xlabel("winding number M")
    axes[2].set_ylabel("|ridge scan cc|")
    axes[2].set_title(f"{idx}: virtual-gauge ridge scan", fontsize=9)
    axes[2].legend(fontsize=8)

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
    ap.add_argument("--h1", type=float, default=120.0)
    ap.add_argument("--r-drag", type=float, default=1e-7)
    ap.add_argument("--nu-scale", type=float, default=0.02)
    ap.add_argument("--tide-amp", type=float, default=2e-2)
    ap.add_argument("--m-max", type=float, default=5.0)
    ap.add_argument("--dm", type=float, default=0.01)
    ap.add_argument("--outdir", type=Path, default=None)
    args = ap.parse_args()

    if args.sanity:
        g_reduced = G * args.drho_rho
        ok = sanity_check(args.dx_km, g_reduced, args.r_drag, args.nu_scale)
        return 0 if ok else 1

    run(args.years, args.dx_km, args.drho_rho, args.h1, args.r_drag,
       args.nu_scale, args.tide_amp, args.idx, args.m_max, args.dm, args.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
