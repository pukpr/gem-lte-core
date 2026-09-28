#!/usr/bin/env python3
"""pna_transfer_function.py — the actual forced (driven, NOT free-
oscillation) response of the pna resolved model's linear operator, as a
function of forcing frequency: |gauge response| / |forcing amplitude| at
each period, computed via the resolvent (i*omega*I - A)^-1 @ B of the
EXACT SAME linear operator built in pna_eigenmode_analysis.py.

Motivation: the free-eigenmode spectrum (pna_eigenmode_analysis.py) does
NOT support a genuine slow (Rossby-wave-timescale) resonance -- the only
high-"ringing-quality" (Q=efold/period > 1) modes are all sub-day, and
the slow (700-1300 day) modes decay in under 1.2 days, i.e. they are
overdamped relaxations, not real sloshing. So the 1-3 month smoothing
"sweet spot" found empirically (pna_resolved_shallow_water.py
--smooth-sigma-months) is NOT explained by a classical resonant
eigenmode. This computes the more directly relevant quantity instead:
the system's actual DRIVEN, forced-response transfer function |H(omega)|
-- for a linear, damped, driven system, the forced response can still
show a real local maximum at frequencies where a heavily-damped complex
eigenvalue's real part is small relative to how far omega sits from its
imaginary part, even without a distinct "ringing" eigenmode -- this is
the RIGOROUS way to check for a real resonance in the driven system,
rather than assuming the free-oscillation spectrum is the whole story.

Usage
-----
    ./pna_transfer_function.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pna_eigenmode_analysis import build_operator                          # noqa: E402
from baltic_resolved_shallow_water import grad                             # noqa: E402
from pdo_resolved_shallow_water import diagonal_tilt_field                 # noqa: E402

ROOT = Path(__file__).resolve().parent


def main() -> int:
    A, n_wet, lat, wet, mask, dx, nx, ny = build_operator()
    g_eq = 9.81
    tide_amp = 0.02
    alpha = 0.5

    tilt = diagonal_tilt_field(mask, alpha)
    dtiltdx, dtiltdy = grad(tilt, dx)
    Bu = (-g_eq * tide_amp * dtiltdx)[mask]
    Bv = (-g_eq * tide_amp * dtiltdy)[mask]
    Bzeta = np.zeros(n_wet)
    B = np.concatenate([Bzeta, Bu, Bv])

    j0, i0 = ny // 2, nx // 2
    gauge_wet_idx = None
    for k, (j, i) in enumerate(wet):
        if j == j0 and i == i0:
            gauge_wet_idx = k
            break
    assert gauge_wet_idx is not None, "gauge point not in wet interior"
    print(f"gauge at grid ({j0},{i0}) -> wet index {gauge_wet_idx} "
          f"(zeta component at state index {gauge_wet_idx})")

    n = A.shape[0]
    I = np.eye(n)

    periods_days = np.logspace(np.log10(0.2), np.log10(3650), 400)
    omega = 2 * np.pi / (periods_days * 86400.0)

    print(f"computing resolvent at {len(omega)} frequencies "
          f"(each a {n}x{n} linear solve)...")
    response = np.empty(len(omega), dtype=complex)
    for k, w in enumerate(omega):
        x = np.linalg.solve(1j * w * I - A, B)
        response[k] = x[gauge_wet_idx]
        if k % 50 == 0:
            print(f"  {k}/{len(omega)}")

    mag = np.abs(response)

    print("\ntop local maxima of |H(period)| (genuine resonant peaks in "
          "the DRIVEN response):")
    is_peak = (mag[1:-1] > mag[:-2]) & (mag[1:-1] > mag[2:])
    peak_idx = np.flatnonzero(is_peak) + 1
    order = peak_idx[np.argsort(-mag[peak_idx])]
    for k in order[:20]:
        print(f"  period={periods_days[k]:9.3f} days   |H|={mag[k]:.4e}")

    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.loglog(periods_days, mag, color="tab:blue")
    for months, lbl in [(0.25, "1wk"), (1.0, "1mo"), (1.5, "1.5mo"),
                        (3.0, "3mo"), (12.0, "1yr")]:
        days = months * 30.44
        ax.axvline(days, color="tab:red", linestyle=":", alpha=0.6)
        ax.text(days, ax.get_ylim()[1] if False else mag.max(), lbl,
               color="tab:red", fontsize=8, rotation=90, va="top")
    ax.set_xlabel("forcing period (days)")
    ax.set_ylabel("|gauge response| / |forcing amplitude|  (driven, "
                 "steady-state)")
    ax.set_title("pna basin: ACTUAL driven transfer function (not free "
                "eigenmodes) -- red = empirical smoothing sweet spots")
    fig.tight_layout()
    out = ROOT / "pna" / "pna_transfer_function.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"\nsaved {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
