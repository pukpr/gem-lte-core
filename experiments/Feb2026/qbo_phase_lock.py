#!/usr/bin/env python3
"""Empirical rotation-number / phase-residual diagnostic for qbo30.

The Ada model (gem-lte-primitives-solution.adb) locks the base winding
number to a single IIR-filter parameter whenever LOCKF is set:

    M(NM) = 1.0 / (Decay * mP)

For qbo30 (LOCKF TRUE, DECAY unset -> default 1.0, mp=0.0528558439226213)
this gives M(NM) = 18.9204 -- the dominant, essentially unbroken ridge in
winding_scalogram.png.

This script applies the standard circle-map / PLL diagnostic for phase
locking: build the continuous winding phase theta(t) = 2*pi*M(NM)*Forcing(t),
unwrap it, fit a straight line, and look at the RESIDUAL. A residual that
stays bounded (wanders but doesn't grow) is the signature of a locked
oscillator; a residual that grows systematically over time is the
signature of slipping (not locked).

As a second, independent, more physically direct check, it also tracks the
actual reversal-to-reversal (zero-crossing) period sequence of the real
Data and Model columns, and compares it against two reference rates:
  - mP's own "free-running" rate: (1/mP) IIR steps * dt (~1.58 yr) -- what
    the leaky-integrator-with-drag mechanism would produce on its own,
    absent any real external tidal entrainment.
  - the real lunar-vs-calendar-month beat (~2.37 yr): 1/(1/T_lunar -
    12/365.242), the same construction that gives the ordinary synodic
    month from the sidereal month, applied here to the draconic month
    (1/27.2122 - 1/365.242 = 1/29.403).

No new mechanism -- both diagnostics are computed directly from
qbo30/lte_results.csv and qbo30/lt.exe.p as they already exist.
"""
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from wavelet_scalogram import load_columns
from lte_forward import read_resp

ROOT = Path(__file__).resolve().parent

# 1/T_draconic - 1/T_year, the same construction that gives the ordinary
# synodic month from the sidereal month (1/27.3216 - 1/365.242 = 1/29.53).
T_DRACONIC = 27.212221
T_YEAR = 365.242
DRACONIC_SYNODIC_DAYS = 1.0 / (1.0 / T_DRACONIC - 1.0 / T_YEAR)
LUNAR_CAL_BEAT_DAYS = 1.0 / (1.0 / DRACONIC_SYNODIC_DAYS - 12.0 / T_YEAR)
LUNAR_CAL_BEAT_YEARS = LUNAR_CAL_BEAT_DAYS / T_YEAR


def locked_M(idx: str) -> tuple[float, float, float]:
    """M(NM) = 1/(Decay*mP), the same formula as the Ada Lock_Freq branch."""
    p = json.loads((ROOT / idx / "lt.exe.p").read_text())
    resp = read_resp(ROOT / idx / "lt.exe.resp")
    decay = float(resp.get("DECAY", 1.0))
    mP = p["mp"]
    return 1.0 / (decay * mP), decay, mP


def phase_residual(year, forcing, M):
    theta = 2.0 * np.pi * M * forcing
    theta_unwrapped = np.unwrap(theta)
    A = np.vstack([year, np.ones_like(year)]).T
    slope, intercept = np.linalg.lstsq(A, theta_unwrapped, rcond=None)[0]
    residual = theta_unwrapped - (slope * year + intercept)
    return residual, slope


def rolling_std(year, residual, window_years=15.0, step_years=2.0):
    centers = np.arange(year.min() + window_years / 2,
                         year.max() - window_years / 2, step_years)
    stds = []
    for c in centers:
        mask = np.abs(year - c) < window_years / 2
        stds.append(residual[mask].std() if mask.sum() > 5 else np.nan)
    return centers, np.array(stds)


def reversal_periods(year, x):
    """Half-cycle length (years) between successive zero-crossings of x,
    doubled to a full-cycle length -- same measure used earlier this
    session's zero-crossing check."""
    xc = x - x.mean()
    zc_idx = np.where(np.diff(np.sign(xc)) != 0)[0]
    zc_year = year[zc_idx]
    full_cycle = np.diff(zc_year) * 2.0
    centers = zc_year[:-1] + np.diff(zc_year) / 2.0
    return centers, full_cycle


def main():
    idx = "qbo30"
    year, dt, model, obs, forcing = load_columns(idx)
    M, decay, mP = locked_M(idx)
    natural_period_mP = (1.0 / mP) * dt

    print(f"[{idx}] Decay={decay:g}  mP={mP:.10f}  locked M(NM)={M:.4f}")
    print(f"  mP-only 'free-running' period: {natural_period_mP:.4f} yr "
          f"((1/mP)={1.0/mP:.3f} steps * dt={dt:.5f} yr/step)")
    print(f"  draconic-synodic (lunar-vs-calendar-month) beat: "
          f"{DRACONIC_SYNODIC_DAYS:.3f} d -> beat period "
          f"{LUNAR_CAL_BEAT_YEARS:.4f} yr")

    residual, slope = phase_residual(year, forcing, M)
    res_centers, res_std = rolling_std(year, residual)
    first_half = residual[year < np.median(year)]
    second_half = residual[year >= np.median(year)]
    print(f"  phase residual: overall std={residual.std():.3f} rad; "
          f"first-half std={first_half.std():.3f}, "
          f"second-half std={second_half.std():.3f} "
          f"({'bounded/locked' if second_half.std() < 1.5 * first_half.std() else 'GROWING (slipping)'})")

    cc_m, per_m = reversal_periods(year, model)
    cc_o, per_o = reversal_periods(year, obs)
    print(f"  realized full-cycle period: Model median={np.median(per_m):.3f} yr, "
          f"Data median={np.median(per_o):.3f} yr")

    fig, axes = plt.subplots(3, 1, figsize=(11, 11))

    axes[0].plot(year, residual, color="tab:blue", linewidth=0.6)
    axes[0].plot(res_centers, res_std, color="tab:red", linewidth=1.5,
                 label="rolling std (15yr window)")
    axes[0].plot(res_centers, -res_std, color="tab:red", linewidth=1.5)
    axes[0].axhline(0, color="black", linewidth=0.5)
    axes[0].set_ylabel("winding-phase residual (rad)")
    axes[0].set_title(f"{idx}: phase residual from linear fit to "
                       f"theta(t)=2*pi*M*Forcing(t), M={M:.4f}  "
                       f"(bounded = locked; growing = slipping)",
                       fontweight="bold")
    axes[0].legend(loc="upper right", fontsize=8)

    axes[1].plot(cc_m, per_m, "o-", color="tab:blue", markersize=3,
                 label="Model (col 2)")
    axes[1].plot(cc_o, per_o, "s-", color="tab:orange", markersize=3,
                 alpha=0.6, label="Data (col 3)")
    axes[1].axhline(natural_period_mP, color="purple", linestyle=":",
                     label=f"mP free-running rate ({natural_period_mP:.3f} yr)")
    axes[1].axhline(LUNAR_CAL_BEAT_YEARS, color="green", linestyle="--",
                     label=f"lunar-calendar beat ({LUNAR_CAL_BEAT_YEARS:.3f} yr)")
    axes[1].axhline(np.median(per_m), color="gray", linestyle="-",
                     linewidth=0.8,
                     label=f"Model median ({np.median(per_m):.3f} yr)")
    axes[1].set_ylabel("realized full-cycle period (yr)")
    axes[1].set_title("reversal-to-reversal period sequence vs. the two "
                       "'natural' reference rates", fontweight="bold")
    axes[1].legend(loc="upper right", fontsize=8)

    axes[2].hist(per_m, bins=20, alpha=0.6, color="tab:blue", label="Model")
    axes[2].hist(per_o, bins=20, alpha=0.6, color="tab:orange", label="Data")
    axes[2].axvline(natural_period_mP, color="purple", linestyle=":")
    axes[2].axvline(LUNAR_CAL_BEAT_YEARS, color="green", linestyle="--")
    axes[2].set_xlabel("full-cycle period (yr)")
    axes[2].set_ylabel("count")
    axes[2].legend(loc="upper right", fontsize=8)

    fig.tight_layout()
    out = ROOT / idx / "qbo_phase_lock.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"  saved {out}")


if __name__ == "__main__":
    main()
