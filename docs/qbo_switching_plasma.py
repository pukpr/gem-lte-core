#!/usr/bin/env python3
"""One-panel schematic: the QBO as a switched (AC) toroidal plasma (2026-10-22). See qbo_switching_plasma.md.
Time scale 1 yr -> 1 s. The probe trace is the comb device of qbo_tokamak.py (pulses -> integrator -> threshold).
    python3 qbo_switching_plasma.py  -> qbo_switching_plasma.png
"""
import json
from pathlib import Path

import matplotlib
import numpy as np
from matplotlib.patches import Ellipse, FancyArrowPatch, FancyBboxPatch, Rectangle

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
__file__ = str(HERE / "qbo_tokamak.py")
exec(open(HERE / "qbo_tokamak.py").read().split("BOUNDS = ")[0])
X = np.array(json.load(open(HERE / "qbo_tokamak.json"))["COMB|threshold"]["x"])
INK, MUTED, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
BLUE, ORANGE, RED, GOLD, GREEN = "#2a78d6", "#eb6834", "#c8322b", "#eda100", "#1baf7a"

fig = plt.figure(figsize=(16, 8.2), facecolor=SURF)
ax = fig.add_axes([0, 0.30, 1, 0.62]); ax.set_xlim(0, 16); ax.set_ylim(0, 6.2); ax.axis("off")
fig.text(0.012, 0.955, "The QBO as a switched toroidal plasma: a polarity-reversing pulsed supply drives a ring whose flow (field) reverses",
         fontsize=13.5, color=INK)


def box(x, y, w, h, title, body, fc):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05,rounding_size=0.12", fc=fc, ec=INK, lw=1.2))
    ax.text(x + w / 2, y + h - 0.32, title, ha="center", va="top", fontsize=11.5, color=INK, weight="bold")
    ax.text(x + w / 2, y + h - 0.85, body, ha="center", va="top", fontsize=9.5, color=INK, linespacing=1.35)


def arrow(a, b, col=INK):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=18, color=col, lw=2))


box(0.3, 3.3, 3.6, 2.6, "1  SWITCH  (the Sun)", "fires twice per cycle,\nflipping polarity each time\n(+, −, +, −, …)\n\n2 pulses per second", "#fde7c2")
box(0.3, 0.3, 3.6, 2.6, "2  MODULATOR  (the Moon)", "sets each pulse's height:\nMoon above / below the\nSun's plane (ecliptic)\n\n13.42 cycles per second", "#dbe9f8")
arrow((3.95, 1.6), (4.9, 3.0)); arrow((3.95, 4.6), (4.9, 3.6))
# capacitor + H-bridge glyph
ax.add_patch(FancyBboxPatch((4.95, 2.35), 1.9, 1.9, boxstyle="round,pad=0.05,rounding_size=0.1", fc="#eef0e6", ec=INK, lw=1.2))
ax.text(5.9, 3.95, "pulsed supply", ha="center", fontsize=10, color=INK, weight="bold")
for xx in (5.35, 6.45):
    ax.plot([xx, xx], [2.6, 3.65], color=INK, lw=1.4)
ax.plot([5.35, 6.45], [3.12, 3.12], color=ORANGE, lw=3)
for (xx, yy) in ((5.35, 3.45), (6.45, 3.45), (5.35, 2.8), (6.45, 2.8)):
    ax.add_patch(Rectangle((xx - 0.12, yy - 0.1), 0.24, 0.2, fc="#fbe0dd", ec=RED, lw=1))
ax.text(5.9, 2.45, "H-bridge switch", ha="center", fontsize=8.5, color=MUTED)
arrow((6.95, 3.3), (7.85, 3.3), ORANGE)
# torus
cx, cy = 10.0, 3.2
ax.add_patch(Ellipse((cx, cy), 4.6, 2.2, fc="#dbe9f8", ec=INK, lw=1.4))
ax.add_patch(Ellipse((cx, cy), 2.0, 0.8, fc=SURF, ec=INK, lw=1.4))
ax.add_patch(Rectangle((cx - 0.22, cy - 1.5), 0.44, 3.0, fc="#f6d9c7", ec=ORANGE, lw=1.5))
ax.text(cx, cy + 1.62, "drive coil", ha="center", fontsize=9, color=ORANGE)
th = np.linspace(-0.1 * np.pi, 0.9 * np.pi, 60)
ax.plot(cx + 1.75 * np.cos(th), cy + 0.82 * np.sin(th), color=BLUE, lw=2.6)
ax.add_patch(FancyArrowPatch((cx + 1.75 * np.cos(th[-3]), cy + 0.82 * np.sin(th[-3])), (cx + 1.75 * np.cos(th[-1]), cy + 0.82 * np.sin(th[-1])),
                             arrowstyle="-|>", mutation_scale=18, color=BLUE, lw=2.6))
th2 = th[::-1] + np.pi
ax.plot(cx + 2.05 * np.cos(th2), cy + 0.98 * np.sin(th2), color=RED, lw=2.6)
ax.add_patch(FancyArrowPatch((cx + 2.05 * np.cos(th2[-3]), cy + 0.98 * np.sin(th2[-3])), (cx + 2.05 * np.cos(th2[-1]), cy + 0.98 * np.sin(th2[-1])),
                             arrowstyle="-|>", mutation_scale=18, color=RED, lw=2.6))
ax.text(cx, 0.55, "3  RING PLASMA  (the equatorial stratosphere)\nflow around the ring reverses at the edge when the\ndrive passes a threshold; each reversal sinks inward",
        ha="center", fontsize=10, color=INK, linespacing=1.35)
ax.add_patch(plt.Circle((cx + 2.3, cy + 0.15), 0.12, fc=GREEN, ec=INK)); ax.text(cx + 2.45, cy + 0.42, "probe", fontsize=9, color=GREEN)
arrow((cx + 2.45, cy), (13.1, cy), GREEN)
ax.text(14.55, 4.7, "probe reads a\nreversing square wave", ha="center", fontsize=10, color=INK)
ax.text(14.55, 0.9, "period = 2.37 s\n(1 s = 1 yr: the QBO's\n28-month period)", ha="center", fontsize=10, color=INK, weight="bold")
ins = fig.add_axes([0.82, 0.52, 0.17, 0.17]); ins.axis("off")
u = model("COMB", X, "threshold")[2][30]; m = (TG >= 1990) & (TG < 2005)
ins.plot(TG[m] - 1990, u[m], color=GREEN, lw=1.8); ins.axhline(0, color=MUTED, lw=0.5)

# bottom strip: the arithmetic
ax2 = fig.add_axes([0.03, 0.02, 0.94, 0.25]); ax2.axis("off")
ax2.text(0, 0.95, "Why 2.37 s: stroboscopic beating, no tuning", fontsize=12, color=INK, weight="bold", va="top")
ax2.text(0, 0.72,
         "The modulator runs at 13.42 cycles per second. Flipping the sign at every pulse shifts it to 14.42; sampled twice a second, "
         "14.42 folds down to\n14.42 − 14 = 0.42 cycles per second, "
         "a 2.37-s beat (the wagon-wheel effect of a strobe). If the switch did not flip\npolarity the beat would be 1.73 s; "
         "if the modulator followed the Moon's declination instead of its latitude it would be 2.72 s. Only the latitude with a "
         "polarity-flipping\nswitch gives the observed QBO period.", fontsize=10, color=INK, va="top", linespacing=1.45)
ax2.text(0, 0.06, "Each ingredient exists in laboratory devices: polarity-switched (AC) tokamaks, reversed-field pinches whose edge field reverses "
         "by itself, oscillating-field current drive,\nand a rotating-annulus fluid experiment that reproduces descending QBO reversals. "
         "See qbo_switching_plasma.md.", fontsize=9.5, color=MUTED, va="bottom", linespacing=1.4)
fig.savefig(HERE / "qbo_switching_plasma.png", dpi=100, facecolor=SURF)
print("plot: qbo_switching_plasma.png")
