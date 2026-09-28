#!/usr/bin/env python3
"""plot_baltic_northsea_lissajous.py -- stacked Lissajous-style phase
diagrams for the windings actually SHARED between Baltic SST and North
Sea SST (the backbone and the dominant M~1.245 shelf-sea anchor).
X-axis = Baltic, Y-axis = North Sea, same construction as
plot_baltic_lissajous.py.

IMPORTANT: both amplitude/phase pairs below were extracted using a
SINGLE SHARED forcing (Baltic's own) for both regions' regressions --
NOT each region's own independently-fit forcing. An earlier version of
this script used each region's own forcing, which produced a spurious
~166 degree "antiphase" reading; that was a methodology bug (two
independently-fit forcings that are 99.66% correlated overall still
decorrelate almost completely once multiplied by M~1.245), not a real
finding. See BALTIC_ADMITTANCE_LISSAJOUS.md section 8 for the full
correction.
"""
from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# (name, M, Baltic_amp, Baltic_phase, NorthSea_amp, NorthSea_phase)
# -- phases both measured against Baltic's own shared forcing
ROWS = [
    ("backbone", 0.2078, 0.0206, np.radians(162.0), 0.0212, np.radians(169.2)),
    ("dominant (shelf-sea anchor)", 1.2467, 0.3411, np.radians(41.9), 0.2202, np.radians(37.8)),
]

theta = np.linspace(0, 2 * np.pi, 400)

fig, axes = plt.subplots(len(ROWS), 2, figsize=(9, 3.4 * len(ROWS)))

for row, (name, M, ma, mp, sa, sp) in enumerate(ROWS):
    x = ma * np.sin(theta + mp)   # Baltic reconstructed contribution
    y = sa * np.sin(theta + sp)   # North Sea reconstructed contribution

    x_sin, y_sin = ma * np.sin(mp), sa * np.sin(sp)          # theta=0
    x_cos, y_cos = ma * np.cos(mp), sa * np.cos(sp)          # theta=pi/2

    xlim = 1.15 * max(abs(x.min()), abs(x.max()))
    ylim = 1.15 * max(abs(y.min()), abs(y.max()))

    dphi = np.degrees((sp - mp + np.pi) % (2 * np.pi) - np.pi)

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
        ax.set_title(f"{name}  M={M:.4f}  phase diff={dphi:+.1f}°\n"
                     f"{label} = ({hx:.4f}, {hy:.4f})", fontsize=9)
        if col == 0:
            ax.set_ylabel("North Sea", fontsize=9)
        if row == len(ROWS) - 1:
            ax.set_xlabel("Baltic", fontsize=9)
        ax.tick_params(labelsize=7)

fig.suptitle("Baltic vs. North Sea SST: Lissajous phase diagrams for their shared windings\n"
             "gray ellipse = full reconstructed-cycle trajectory; "
             "left = sin-component (θ=0), right = cos-component (θ=π/2)\n"
             "Both phases measured against ONE shared forcing (Baltic's own) -- "
             "corrected from an earlier each-region's-own-forcing bug",
             fontsize=10, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.90))

out = "/home/paul/eval/gem-lte-core/experiments/Feb2026/baltic_northsea_lissajous.png"
fig.savefig(out, dpi=140)
print(f"[saved] {out}")
