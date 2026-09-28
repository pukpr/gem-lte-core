#!/usr/bin/env python3
"""pna_eigenmode_analysis.py — computes the ACTUAL linear normal modes
(eigenperiods) of the PNA resolved model's own discretized shallow-water
operator, INCLUDING Coriolis (full f=2*Omega*sin(lat), varying across
the 20N-60N domain -- a genuine beta effect, not a single reference
latitude), linear drag, and diffusion -- exactly the operator run_solver
integrates, not an idealized textbook rectangular-basin formula.

Motivation: the naive Coriolis-free rectangular-basin seiche formula
(2*L/c) gives eigenperiods of 0.6-3.1 DAYS for this grid/depth -- far
too fast to explain why Gaussian-smoothing the forcing at ~1-3 MONTH
widths (spreading the same phase change over a longer window, per
signal_operators/pna_resolved_shallow_water.py's --smooth-sigma-months)
is what reveals pna's real winding numbers. Coriolis fundamentally
changes the mode structure: with f varying across latitude (a real beta
effect), the basin supports both fast inertia-gravity (Poincare-like)
modes near the gravity-wave timescale AND much slower Rossby-wave-like
modes -- this computes the REAL spectrum for this exact configuration
rather than assuming which branch is relevant.

Method: the model is already exactly linear (no advective nonlinearity
anywhere in run_solver) -- the continuous-time operator d(state)/dt =
A @ state (state = [zeta, u, v] stacked over interior wet points) is
built via the impulse-response/probing method: apply the same
grad/laplacian discretization and Coriolis coupling used by run_solver
(in its dt->0 / continuous-time limit -- the exact-rotation sub-step is
the discrete-stable equivalent of the standard Coriolis term for small
dt) to each unit basis vector, one column of A at a time. Then
np.linalg.eig gives the exact discrete spectrum for this grid.

Usage
-----
    ./pna_eigenmode_analysis.py
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
from baltic_resolved_shallow_water import laplacian, grad, G                # noqa: E402
from amo_resolved_shallow_water import full_fcor                           # noqa: E402
from nao_resolved_shallow_water import build_grid                          # noqa: E402

ROOT = Path(__file__).resolve().parent


def build_operator(dx_km=300.0, lat_south=20.0, lat_north=60.0, g_eq=9.81,
                   depth=300.0, r_drag=1e-7, nu_scale=0.02):
    _, mask, open_mask, dx, nx, ny, lat = build_grid(dx_km, lat_south, lat_north)
    H = np.full((ny, nx), depth)
    f_cor = full_fcor(lat, nx)
    c = np.sqrt(g_eq * depth)
    dt_ref = 0.4 * dx / c
    nu = nu_scale * dx ** 2 / dt_ref

    wet = np.argwhere(mask)
    n_wet = len(wet)
    idx_of = -np.ones((ny, nx), dtype=int)
    for k, (j, i) in enumerate(wet):
        idx_of[j, i] = k
    print(f"grid {nx}x{ny}, {n_wet} wet interior points, state size "
          f"{3*n_wet} (zeta,u,v)")

    def rhs(zeta, u, v):
        dzdx, dzdy = grad(zeta, dx)
        dudx, _ = grad(u, dx)
        _, dvdy = grad(v, dx)
        zeta_dot = -H * (dudx + dvdy) + nu * laplacian(zeta, dx)
        u_dot = -g_eq * dzdx - r_drag * u + nu * laplacian(u, dx) + f_cor * v
        v_dot = -g_eq * dzdy - r_drag * v + nu * laplacian(v, dx) - f_cor * u
        zeta_dot[~mask] = 0.0
        u_dot[~mask] = 0.0
        v_dot[~mask] = 0.0
        return zeta_dot, u_dot, v_dot

    n = 3 * n_wet
    A = np.zeros((n, n))
    t0 = time.time()
    for k in range(n):
        field = k // n_wet
        wk = k % n_wet
        j, i = wet[wk]
        zeta = np.zeros((ny, nx)); u = np.zeros((ny, nx)); v = np.zeros((ny, nx))
        if field == 0:
            zeta[j, i] = 1.0
        elif field == 1:
            u[j, i] = 1.0
        else:
            v[j, i] = 1.0
        zeta_dot, u_dot, v_dot = rhs(zeta, u, v)
        col = np.concatenate([zeta_dot[mask], u_dot[mask], v_dot[mask]])
        A[:, k] = col
        if k % 200 == 0:
            print(f"  building operator: column {k}/{n} "
                  f"({time.time()-t0:.1f}s elapsed)")
    print(f"  operator built in {time.time()-t0:.1f}s")
    return A, n_wet, lat, wet, mask, dx, nx, ny


def main() -> int:
    A, n_wet, lat, wet, mask, dx, nx, ny = build_operator()

    print("computing eigenvalues (dense, may take a while for "
          f"{A.shape[0]}x{A.shape[1]})...")
    t0 = time.time()
    eigvals, eigvecs = np.linalg.eig(A)
    print(f"  done in {time.time()-t0:.1f}s")

    # eigenvalues of a damped oscillatory system are complex:
    # lambda = -decay_rate + i*omega -- omega gives the oscillation
    # period, -decay_rate (should be <=0 for a stable/dissipative system)
    # gives the e-folding decay time.
    omega = np.abs(eigvals.imag)
    decay = -eigvals.real
    with np.errstate(divide="ignore"):
        period_days = np.where(omega > 1e-12, 2 * np.pi / omega / 86400.0,
                               np.inf)
    with np.errstate(divide="ignore"):
        efold_days = np.where(decay > 1e-12, 1.0 / decay / 86400.0, np.inf)

    order = np.argsort(period_days)
    print("\nlongest-period genuine oscillatory modes (period, e-folding "
          "decay time, both days):")
    shown = 0
    for k in order[::-1]:
        if not np.isfinite(period_days[k]) or period_days[k] < 0.5:
            continue
        print(f"  period={period_days[k]:9.2f} days   "
              f"decay e-fold={efold_days[k]:9.2f} days   "
              f"(lambda={eigvals[k]:.3e})")
        shown += 1
        if shown >= 25:
            break

    print(f"\nmax stable (real part <= ~0): "
          f"{np.all(eigvals.real <= 1e-8)}")
    print(f"largest real part (should be near 0, not positive -- a "
          f"positive value means the discretized operator is itself "
          f"unstable): {eigvals.real.max():.4e}")

    # histogram of periods on a log scale
    finite_periods = period_days[np.isfinite(period_days)]
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.hist(np.log10(finite_periods[finite_periods > 0.01]), bins=60,
           color="tab:blue")
    for target, lbl in [(1/12, "1mo/12"), (1.0, "1mo"), (1.5, "1.5mo"),
                        (3.0, "3mo"), (30.0*1.0, "")]:
        pass
    for months, lbl in [(0.25, "1wk"), (1.0, "1mo"), (1.5, "1.5mo"),
                        (3.0, "3mo")]:
        days = months * 30.44
        ax.axvline(np.log10(days), color="tab:red", linestyle=":")
        ax.text(np.log10(days), ax.get_ylim()[1]*0.9, lbl, color="tab:red",
               fontsize=8, rotation=90, va="top")
    ax.set_xlabel("log10(period, days)")
    ax.set_ylabel("count of normal modes")
    ax.set_title("pna basin: real Coriolis-included linear normal mode "
                "period spectrum (red = smoothing sweet spots found "
                "empirically)")
    fig.tight_layout()
    out = ROOT / "pna" / "pna_eigenmode_spectrum.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"\nsaved {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
