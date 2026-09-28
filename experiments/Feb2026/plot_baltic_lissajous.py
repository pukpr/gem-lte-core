#!/usr/bin/env python3
"""plot_baltic_lissajous.py -- stacked Lissajous-style phase diagrams,
one row per shared winding K1..K7, comparing Baltic MSL (x-axis) against
Baltic Kaplan SST (y-axis).

Both series' reconstructed contribution from winding K_n is the same
sinusoid Amp*sin(theta + Phase), sharing the same theta (the winding's
own cycle) since it's the same physical frequency in both series. That
means Amp*sin(Phase) and Amp*cos(Phase) are just this one curve evaluated
at theta=0 and theta=pi/2 -- so both requested views are two marked
points on the SAME Lissajous ellipse, not two different curves. Each row
draws that ellipse in both panels; left panel highlights the theta=0
(sin-component) point, right panel highlights the theta=pi/2
(cos-component) point.

A straight diagonal line through the ellipse = in phase (or antiphase);
a circle/wide ellipse = quadrature; the tilt and eccentricity encode the
amplitude ratio and phase difference together, same as reading a
Lissajous figure off an oscilloscope in XY mode.
"""
from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# (name, M, MSL_amp, MSL_phase, SST_amp, SST_phase)
ROWS = [
    ("K1", -0.9247, 11.62297655729, -2.36338427942, 0.22577263392, -2.22944978289),
    ("K2",  0.1219,  7.81932176333,  1.76531972834, 0.14169038012, -2.24141362705),
    ("K3",  0.2076,  4.93930575688,  0.65176840196, 0.03680704287,  1.03135137822),
    ("K4",  2.0756,  4.93904939633, -2.22391915483, 0.12861208211,  1.63310242181),
    ("K5",  1.8681,  6.02491961827, -2.69284012571, 0.19895822950, -0.00547796482),
    ("K6",  1.2454, 31.64637364627, -2.60132542489, 0.36205225099,  0.77230536907),
    ("K7",  2.2832,  6.91947250481, -0.36796701026, 0.10539246400, -2.46314136735),
]

theta = np.linspace(0, 2 * np.pi, 400)

fig, axes = plt.subplots(len(ROWS), 2, figsize=(9, 3.1 * len(ROWS)))

for row, (name, M, ma, mp, sa, sp) in enumerate(ROWS):
    x = ma * np.sin(theta + mp)   # MSL reconstructed contribution
    y = sa * np.sin(theta + sp)   # SST reconstructed contribution

    x_sin, y_sin = ma * np.sin(mp), sa * np.sin(sp)          # theta=0
    x_cos, y_cos = ma * np.cos(mp), sa * np.cos(sp)          # theta=pi/2

    xlim = 1.15 * max(abs(x.min()), abs(x.max()))
    ylim = 1.15 * max(abs(y.min()), abs(y.max()))

    for col, (hx, hy, label, color) in enumerate([
        (x_sin, y_sin, "Amp·sin(Phase)", "tab:blue"),
        (x_cos, y_cos, "Amp·cos(Phase)", "tab:orange"),
    ]):
        ax = axes[row, col]
        ax.plot(x, y, color="0.55", lw=1.1)
        ax.axhline(0, color="0.85", lw=0.7, zorder=0)
        ax.axvline(0, color="0.85", lw=0.7, zorder=0)
        ax.annotate("", xy=(hx, hy), xytext=(0, 0),
                    arrowprops=dict(arrowstyle="-|>", color=color, lw=1.8))
        ax.scatter([hx], [hy], color=color, s=28, zorder=5)
        ax.set_xlim(-xlim, xlim)
        ax.set_ylim(-ylim, ylim)
        ax.set_title(f"{name}  (M={M:+.4f})   {label} = ({hx:.3g}, {hy:.3g})",
                     fontsize=9)
        if col == 0:
            ax.set_ylabel("SST", fontsize=9)
        if row == len(ROWS) - 1:
            ax.set_xlabel("MSL", fontsize=9)
        ax.tick_params(labelsize=7)

fig.suptitle("Baltic MSL vs. Kaplan SST: Lissajous phase diagrams per shared winding\n"
             "gray ellipse = full reconstructed-cycle trajectory; "
             "left = sin-component (θ=0), right = cos-component (θ=π/2) "
             "-- both points lie on the same curve",
             fontsize=11, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.965))

out = "/home/paul/eval/gem-lte-core/experiments/Feb2026/baltic_lissajous.png"
fig.savefig(out, dpi=140)
print(f"[saved] {out}")
