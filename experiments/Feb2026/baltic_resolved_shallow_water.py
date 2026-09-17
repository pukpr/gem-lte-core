#!/usr/bin/env python3
"""baltic_resolved_shallow_water.py — the "eventual resolved fluid-dynamics
run": a genuinely spatially-resolved, forced, rotating shallow-water PDE
solver on an idealized Baltic-shaped 2D grid, replacing the lumped 6-box
network (baltic_eigenmode_check.py) with something that actually resolves
real wave propagation, coastline geometry, and Coriolis effects — none of
which a lumped network can represent at all.

Scope, stated plainly: this is a coarse (25 km), LINEAR, single-active-
layer (reduced-gravity, i.e. the baroclinic/internal mode -- the one
implicated by every other check this session), rotating shallow-water
model on an idealized (not surveyed) Baltic-shaped domain, with linear
bottom drag and a small numerical diffusion for stability. It is NOT a
general circulation model: no nonlinear advection, no wind, no
thermodynamics, no real bathymetry. It IS a genuine PDE integration with
real 2D wave dynamics, which is the specific thing missing from every
other check built this session.

Method
------
Co-located (non-staggered) grid for simplicity and robustness -- a fully
staggered Arakawa C-grid is the more standard choice for a research model,
but a linear, lightly-diffused system on a co-located grid is far less
prone to computational-mode pathologies, and simplicity here is worth more
than that last increment of rigor for a first resolved pass:

    d(zeta)/dt = -H * (du/dx + dv/dy)
    du/dt = +f*v - g'*d(zeta)/dx - g'*d(Ftide)/dx - r*u + nu*laplacian(u)
    dv/dt = -f*u - g'*d(zeta)/dy - g'*d(Ftide)/dy - r*v + nu*laplacian(v)

Ftide(x,y,t) = tide_amp * Forcing(t) * y_norm(x,y) -- the real, already-
built lunisolar Forcing(t) manifold (cv_rolling_blocked.prepare()),
applied as a spatially-tilting body force (a linear potential in the
north-south coordinate), the simplest non-trivial way to inject an
external tidal-type forcing with a genuine, non-zero pressure gradient
(a spatially UNIFORM potential would have zero gradient and force
nothing) -- analogous to how the real degree-2 tidal potential tilts the
sea surface, without claiming to reproduce its exact spherical-harmonic
shape.

Time-stepping: forward-backward (update u,v from the current zeta, then
update zeta from the JUST-UPDATED u,v) -- explicit, simple, and stable for
the shallow-water gravity-wave part under the usual CFL condition, unlike
plain forward Euler (which is unconditionally unstable for an undamped
wave equation).

Geometry: an idealized (not surveyed) elongated main channel — shoaling
from Baltic-Proper-like depth in the south to Bothnian-Bay-like depth in
the north — with one side arm (Gulf-of-Finland-like) and a narrow,
shallow outlet (Danish-Straits-like) opening onto a fixed-level (zeta=0)
"North Sea" boundary strip. Depths reuse baltic_eigenmode_check.py's own
basin numbers, not new ones.

Output: a virtual gauge at the main basin's center, sampled monthly (same
cadence as the real .dat files), run through the SAME validated ridge-scan
machinery (mode_evidence.fit_at_M) used on real data and on the synthetic
ground truth -- an apples-to-apples comparison, not a new ad hoc metric.

Usage
-----
    ./baltic_resolved_shallow_water.py --sanity            # unforced pulse test
    ./baltic_resolved_shallow_water.py --years 145                 # full run
    ./baltic_resolved_shallow_water.py --years 145 --tide-amp 5e-3
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
from cv_rolling_blocked import prepare                      # noqa: E402
from mode_evidence import fit_at_M, sweep_M                 # noqa: E402

ROOT = Path(__file__).resolve().parent
G = 9.81
OMEGA_EARTH = 7.2921e-5  # rad/s


# ---------------------------------------------------------------------------
# Idealized geometry: reuses baltic_eigenmode_check.py's own basin depths
# ---------------------------------------------------------------------------

def build_grid(dx_km: float = 25.0):
    """Idealized Baltic-shaped domain on a co-located grid. Returns
    (H, mask, open_mask, lat, dx, nx, ny). H = depth (m), mask = True where
    wet, open_mask = True at the fixed-level North-Sea boundary strip."""
    nx, ny = 44, 64  # ~1100 km x 1600 km at dx=25km -- Baltic's rough extent
    dx = dx_km * 1e3
    H = np.zeros((ny, nx))
    mask = np.zeros((ny, nx), dtype=bool)

    # main channel: x in [14,28), shoaling south (Baltic Proper, 60m) to
    # north (Bothnian Bay, 40m) -- rows 0 (south) .. ny-1 (north)
    x0, x1 = 14, 28
    for j in range(ny):
        frac = j / (ny - 1)  # 0 = south, 1 = north
        if frac < 0.55:
            depth = 60.0       # Baltic Proper
        elif frac < 0.80:
            depth = 60.0       # Bothnian Sea
        else:
            depth = 40.0       # Bothnian Bay
        H[j, x0:x1] = depth
        mask[j, x0:x1] = True

    # Gulf-of-Finland-like side arm: east from the northern part of the
    # "Baltic Proper" section, rows around frac~0.45-0.55
    gf_j0, gf_j1 = int(0.42 * ny), int(0.55 * ny)
    H[gf_j0:gf_j1, x1:x1 + 10] = 40.0
    mask[gf_j0:gf_j1, x1:x1 + 10] = True

    # Gulf-of-Riga-like side arm: smaller, south of the Finland arm
    gr_j0, gr_j1 = int(0.30 * ny), int(0.38 * ny)
    H[gr_j0:gr_j1, x1:x1 + 6] = 25.0
    mask[gr_j0:gr_j1, x1:x1 + 6] = True

    # Danish-Straits-like narrow, shallow outlet at the south end,
    # narrowing in x, opening onto a fixed-level North Sea strip
    strait_j0, strait_j1 = 0, 6
    sx0, sx1 = 18, 22
    H[strait_j0:strait_j1, sx0:sx1] = 15.0
    mask[strait_j0:strait_j1, sx0:sx1] = True
    open_mask = np.zeros((ny, nx), dtype=bool)
    H[0, sx0:sx1] = 15.0
    mask[0, sx0:sx1] = True
    open_mask[0, sx0:sx1] = True

    # latitude field for the Coriolis parameter: Baltic spans ~54-66N
    lat = 54.0 + 12.0 * (np.arange(ny) / (ny - 1))
    lat2d = np.repeat(lat[:, None], nx, axis=1)
    return H, mask, open_mask, lat2d, dx, nx, ny


def y_norm_field(mask: np.ndarray) -> np.ndarray:
    """North-south normalized coordinate in [-0.5, 0.5], for the tilting
    tidal body force -- gives Ftide a genuine, non-zero gradient."""
    ny, nx = mask.shape
    j = np.arange(ny)
    return np.repeat(((j / (ny - 1)) - 0.5)[:, None], nx, axis=1)


# ---------------------------------------------------------------------------
# The solver
# ---------------------------------------------------------------------------

def laplacian(f: np.ndarray, dx: float) -> np.ndarray:
    lap = np.zeros_like(f)
    lap[1:-1, 1:-1] = (f[2:, 1:-1] + f[:-2, 1:-1] + f[1:-1, 2:] + f[1:-1, :-2]
                      - 4.0 * f[1:-1, 1:-1]) / dx ** 2
    return lap


def grad(f: np.ndarray, dx: float) -> tuple[np.ndarray, np.ndarray]:
    dfdy = np.zeros_like(f)
    dfdx = np.zeros_like(f)
    dfdy[1:-1, :] = (f[2:, :] - f[:-2, :]) / (2 * dx)
    dfdx[:, 1:-1] = (f[:, 2:] - f[:, :-2]) / (2 * dx)
    return dfdx, dfdy


def run_solver(H, mask, open_mask, f_cor, dx, g_reduced, r_drag, nu,
              tide_amp, forcing_fn, n_steps, dt, save_every, gauge_j, gauge_i,
              ny, nx, diag_every: int = 0, init_zeta: np.ndarray | None = None,
              tilt_field: np.ndarray | None = None, save_row: bool = False):
    """`save_row`: if True, additionally records the full zeta ROW at
    gauge_j (all nx columns) at every save step -- lets a caller sample a
    MOVING virtual gauge (a location that shifts with the forcing itself,
    the material-coordinate/Lagrangian-pullback picture the original
    sin(k*M(t)) derivation actually rests on) by post-hoc interpolating
    along that row, instead of only ever reading one fixed grid point."""
    """Coriolis is applied as an EXACT rotation sub-step, not a naive
    explicit forward-Euler term -- forward Euler on a pure rotation
    (du/dt=fv, dv/dt=-fu) has amplification factor sqrt(1+(f*dt)^2) > 1
    for ANY dt > 0, i.e. it is unconditionally unstable regardless of CFL.
    An exact rotation (cos/sin) has amplification exactly 1, so this
    operator-split step carries no spurious energy growth; the remaining
    pressure-gradient/friction/diffusion terms are still explicit
    forward-backward, which IS conditionally stable under the usual CFL
    limit on wave speed.

    `tilt_field` is the spatial pattern the tidal body force is applied
    through (Ftide = tide_amp*forcing_fn(t)*tilt_field) -- explicit rather
    than hardcoded, since different domains want different tilt directions
    (Baltic: north-south; the equatorial nino4 channel: zonal) and a
    caller's gauge point landing on a symmetry node of the WRONG tilt
    field silently reads exactly zero forever, which is exactly the bug
    this parameter replaces (see nino4_resolved_shallow_water.py's
    module docstring / [[baltic-resolved-shallow-water]] memory). Defaults
    to the north-south tilt for backward compatibility with existing
    Baltic callers."""
    ynorm = y_norm_field(mask) if tilt_field is None else tilt_field
    zeta = np.zeros((ny, nx)) if init_zeta is None else init_zeta.copy()
    u = np.zeros((ny, nx))
    v = np.zeros((ny, nx))
    cos_f, sin_f = np.cos(f_cor * dt), np.sin(f_cor * dt)
    gauge = np.empty(n_steps // save_every)
    row_history = np.empty((n_steps // save_every, nx)) if save_row else None
    diag = [] if diag_every else None
    g_idx = 0
    max_abs_zeta = 0.0
    for step in range(n_steps):
        t = step * dt
        Ftide = tide_amp * forcing_fn(t) * ynorm
        dzdx, dzdy = grad(zeta, dx)
        dFdx, dFdy = grad(Ftide, dx)

        u_rot = u * cos_f + v * sin_f
        v_rot = -u * sin_f + v * cos_f
        u_new = (u_rot - dt * g_reduced * (dzdx + dFdx) - dt * r_drag * u_rot
                + dt * nu * laplacian(u_rot, dx))
        v_new = (v_rot - dt * g_reduced * (dzdy + dFdy) - dt * r_drag * v_rot
                + dt * nu * laplacian(v_rot, dx))
        u_new[~mask] = 0.0
        v_new[~mask] = 0.0
        u, v = u_new, v_new

        dudx, _ = grad(u, dx)
        _, dvdy = grad(v, dx)
        zeta_new = zeta - dt * H * (dudx + dvdy) + dt * nu * laplacian(zeta, dx)
        zeta_new[open_mask] = 0.0
        zeta_new[~mask] = 0.0
        zeta = zeta_new

        max_abs_zeta = max(max_abs_zeta, float(np.max(np.abs(zeta))))
        if not np.all(np.isfinite(zeta)):
            raise FloatingPointError(f"blew up at step {step} (t={t/86400:.1f} d)")
        if diag_every and step % diag_every == 0:
            ke = 0.5 * np.sum(H[mask] * (u[mask] ** 2 + v[mask] ** 2))
            pe = 0.5 * g_reduced * np.sum(zeta[mask] ** 2)
            diag.append((t / 86400.0, ke + pe, float(np.max(np.abs(zeta)))))
        if (step + 1) % save_every == 0:
            gauge[g_idx] = zeta[gauge_j, gauge_i]
            if save_row:
                row_history[g_idx] = zeta[gauge_j, :]
            g_idx += 1
    row_out = row_history[:g_idx] if save_row else None
    return gauge[:g_idx], max_abs_zeta, diag, row_out


# ---------------------------------------------------------------------------
# Sanity check: unforced pulse, no forcing, check it propagates & decays
# ---------------------------------------------------------------------------

def sanity_check(dx_km: float, g_reduced: float, r_drag: float, nu_scale: float):
    """Unforced pulse test, via the SAME run_solver() the real run uses (no
    duplicated physics to drift out of sync) -- a 1 m pulse should decay
    smoothly under friction/diffusion, not blow up or oscillate wildly."""
    H, mask, open_mask, lat2d, dx, nx, ny = build_grid(dx_km)
    f_cor = 2 * OMEGA_EARTH * np.sin(np.radians(lat2d))
    c = np.sqrt(g_reduced * np.median(H[mask]))
    dt = 0.4 * dx / c
    nu = nu_scale * dx ** 2 / dt
    print("-- sanity check: unforced pulse --")
    print(f"  grid {nx}x{ny} at dx={dx_km:g}km, wave speed c={c:.3f} m/s, "
          f"dt={dt/3600:.2f} hr (CFL Courant={c*dt/dx:.2f})")

    j0, i0 = int(0.5 * ny), 21
    n_steps = int(60 * 86400 / dt)  # 60 days
    diag_every = max(1, n_steps // 20)

    init_zeta = np.zeros((ny, nx))
    init_zeta[j0, i0] = 1.0  # 1 m pulse, main channel center
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


# ---------------------------------------------------------------------------
# Full forced run
# ---------------------------------------------------------------------------

def run(years: float, dx_km: float, drho_rho: float, h1: float, r_drag: float,
       nu_scale: float, tide_amp: float, idx: str, m_max: float, dm: float,
       outdir: Path | None) -> None:
    H, mask, open_mask, lat2d, dx, nx, ny = build_grid(dx_km)
    f_cor = 2 * OMEGA_EARTH * np.sin(np.radians(lat2d))
    g_reduced = G * drho_rho
    H_active = np.minimum(H, h1)
    c = np.sqrt(g_reduced * np.median(H_active[mask]))
    dt = 0.4 * dx / c
    nu = nu_scale * dx ** 2 / dt  # numerical diffusion scaled to grid/dt

    prep = prepare(idx)
    dates, forcing = prep["dates"], prep["forcing"]
    t_years = dates - dates[0]

    def forcing_fn(t_sec: float) -> float:
        t_yr = t_sec / (365.2422 * 86400.0)
        return float(np.interp(t_yr, t_years, forcing))

    dt_out_target = 1.0 / 12.0 * 365.2422 * 86400.0  # monthly, seconds
    save_every = max(1, round(dt_out_target / dt))
    n_steps = int(years * 365.2422 * 86400.0 / dt)
    n_steps -= n_steps % save_every

    gauge_j, gauge_i = int(0.5 * ny), 21  # main basin center
    print(f"-- resolved shallow-water run: {idx} --")
    print(f"  grid {nx}x{ny} at dx={dx_km:g}km, active-layer c={c:.3f} m/s "
          f"(g'={g_reduced:.4f}, h1={h1:.0f}m capped depth)")
    print(f"  dt={dt/3600:.2f} hr, {n_steps} steps ({years:g} sim-years), "
          f"saving every {save_every} steps (~monthly)")

    t0 = time.time()
    gauge, max_abs, _, _ = run_solver(H_active, mask, open_mask, f_cor, dx, g_reduced,
                                   r_drag, nu, tide_amp, forcing_fn, n_steps, dt,
                                   save_every, gauge_j, gauge_i, ny, nx)
    elapsed = time.time() - t0
    print(f"  done in {elapsed:.1f}s wall-clock; max|zeta| reached {max_abs:.4f} m")

    gauge_dates = dates[0] + np.arange(len(gauge)) / 12.0
    gauge_forcing = np.interp(gauge_dates, dates, forcing)

    m_grid = np.arange(dm, m_max + dm / 2, dm)
    m_grid_out, curve = sweep_M(gauge_forcing, gauge, (dm, m_max), len(m_grid))
    i_star = int(np.argmax(np.abs(curve)))
    print(f"  virtual-gauge ridge scan (same method as real/synthetic data): "
          f"M_hat={m_grid_out[i_star]:.3f}  cc={curve[i_star]:+.3f}")
    print(f"  (for comparison: real baltic's independently-found ridge is "
          f"M*=1.250, per cv_ridge_transfer.py)")

    out_dir = outdir if outdir is not None else ROOT / idx
    plot_result(idx, H, mask, open_mask, gauge_dates, gauge, m_grid_out, curve,
               out_dir / "baltic_resolved_shallow_water.png")


def plot_result(idx, H, mask, open_mask, gauge_dates, gauge, m_grid, curve,
                out_path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
    depth_show = np.where(mask, H, np.nan)
    im = axes[0].imshow(depth_show, origin="lower", cmap="Blues")
    axes[0].imshow(np.where(open_mask, 1, np.nan), origin="lower", cmap="Reds",
                  vmin=0, vmax=1, alpha=0.8)
    axes[0].set_title("idealized geometry (depth, m)\nred = open North Sea boundary",
                      fontsize=9)
    fig.colorbar(im, ax=axes[0], shrink=0.8)
    axes[0].plot([21], [int(0.5 * mask.shape[0])], "r*", markersize=12)

    axes[1].plot(gauge_dates, gauge, linewidth=0.6)
    axes[1].set_xlabel("year")
    axes[1].set_ylabel("virtual gauge zeta (m)")
    axes[1].set_title(f"{idx}: simulated interface displacement", fontsize=9)

    axes[2].plot(m_grid, np.abs(curve))
    axes[2].axvline(1.25, color="tab:red", linestyle=":",
                    label="real baltic's M*=1.250")
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
    ap.add_argument("--sanity", action="store_true", help="run the unforced sanity check and exit")
    ap.add_argument("--idx", default="baltic")
    ap.add_argument("--years", type=float, default=145.0)
    ap.add_argument("--dx-km", type=float, default=25.0)
    ap.add_argument("--drho-rho", type=float, default=0.005)
    ap.add_argument("--h1", type=float, default=50.0)
    ap.add_argument("--r-drag", type=float, default=1e-6,
                     help="linear bottom-drag coefficient (1/s)")
    ap.add_argument("--nu-scale", type=float, default=0.02,
                     help="numerical diffusion, as a fraction of the "
                          "dx^2/dt stability scale")
    ap.add_argument("--tide-amp", type=float, default=2e-2,
                     help="tidal body-force amplitude (m, at |y_norm|=0.5)")
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
