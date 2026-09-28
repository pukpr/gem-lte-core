#!/usr/bin/env python3
"""pna_per_winding_check.py — per the user's reminder that real pna data
can be FILTERED to help guide the process: applies winding_transform
(the same local, Gaussian-windowed basis used throughout this project's
winding_scalogram/winding_rank tooling, unmodified) directly at PNA's
THREE established real winding numbers (1.460/2.970/4.730) -- not a
scan, a fixed-M local projection -- to BOTH real (filtered-by-
construction, since this IS a narrowband local projection) PNA data and
each smoothing-sweep variant's own best-performing gauge output, giving
a genuinely TIME-RESOLVED, per-winding comparison instead of a single
aggregate ridge-scan number.

Usage
-----
    ./pna_per_winding_check.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from winding_scalogram import standardize, winding_transform               # noqa: E402
from signal_operators import delayed_difference                           # noqa: E402

ROOT = Path(__file__).resolve().parent
IDX = "pna"
REAL_M = [1.460, 2.970, 4.730]
SIGMA = 15.0
T0_STEP = 5.0


def load(tag):
    d = np.load(ROOT / IDX / f"pna_resolved_check_cache{tag}.npz")
    return d["gauge_dates"], d["gauge_forcing"], d["fixed"], d["moving"]


def main() -> int:
    dat = np.loadtxt(ROOT / IDX / f"{IDX}.dat")

    # best-performing series found per winding, from the smoothing sweep
    gd15, gf15, fixed15, moving15 = load("_smooth1p5")
    gd3, gf3, fixed3, moving3 = load("_smooth3p0")
    diff3 = delayed_difference(fixed3, 12)
    diff3_dates = gd3[12:]
    diff3_forcing = gf3[12:]

    obs15 = np.interp(gd15, dat[:, 0], dat[:, 1])
    obs3diff = np.interp(diff3_dates, dat[:, 0], dat[:, 1])

    cases = [
        (1.460, "MOVING, sigma=1.5mo", gd15, moving15, gf15, obs15),
        (2.970, "MOVING, sigma=1.5mo", gd15, moving15, gf15, obs15),
        (4.730, "FIXED 12mo-diff, sigma=3.0mo", diff3_dates, diff3,
         diff3_forcing, obs3diff),
    ]

    fig, axes = plt.subplots(3, 1, figsize=(11, 10))
    for ax, (M, label, t, series, forcing, obs) in zip(axes, cases):
        t0_grid = np.arange(t[0] + SIGMA / 2, t[-1] - SIGMA / 2 + 1e-9,
                            T0_STEP)
        m_grid = np.array([M])
        G_obs, edge_mask, _ = winding_transform(t, standardize(obs), forcing,
                                                m_grid, t0_grid, SIGMA)
        G_mod, _, _ = winding_transform(t, standardize(series), forcing,
                                        m_grid, t0_grid, SIGMA)
        mag_obs = np.abs(G_obs[0])
        mag_mod = np.abs(G_mod[0])
        valid = ~edge_mask
        r = np.corrcoef(mag_obs[valid], mag_mod[valid])[0, 1]
        print(f"M={M}: corr(|G_obs|, |G_model|) over time, edge-excluded "
              f"= {r:+.3f}   [{label}]")

        ax2 = ax.twinx()
        ax.plot(t0_grid, mag_obs, color="0.4", label="real pna |G(t)|")
        ax2.plot(t0_grid, mag_mod, color="tab:red", alpha=0.8,
                label=f"model |G(t)| [{label}]")
        for j in np.where(edge_mask)[0]:
            ax.axvspan(t0_grid[j] - T0_STEP / 2, t0_grid[j] + T0_STEP / 2,
                      color="0.9", zorder=0)
        ax.set_title(f"M={M}: local winding amplitude over time "
                    f"(r={r:+.3f}, gray shading = edge-of-record)",
                    fontsize=9)
        ax.set_ylabel("real |G(t)|", color="0.4")
        ax2.set_ylabel("model |G(t)|", color="tab:red")
        ax.set_xlabel("year (window center)")
        lines1, labels1 = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax.legend(lines1 + lines2, labels1 + labels2, fontsize=7,
                 loc="upper left")

    fig.tight_layout()
    out = ROOT / IDX / "pna_per_winding_check.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"\nsaved {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
