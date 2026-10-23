#!/usr/bin/env python3
"""Schematic of the modified tokamak (reversing poloidal drive) and the QBO record read as an MHD shot (2026-10-22).

Time scale: 1 year of the QBO record = 1 second of discharge (the 1953-2026 record is a 73-s shot; the QBO period
2.37 yr becomes 2.37 s; the semiannual comb becomes a 2-Hz pulsed supply; the draconic month becomes a 13.42-Hz
modulation). "Wind" becomes the edge toroidal field at three probe radii, 1 m/s -> 1 mT.
Uses the fitted comb/threshold device from qbo_tokamak.py (qbo_tokamak.json).
    python3 qbo_tokamak_device.py  -> qbo_tokamak_device.png (schematic + circuit), qbo_tokamak_shot.png (data as a shot)
"""
import json
from pathlib import Path

import matplotlib
import numpy as np
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, FancyBboxPatch, Rectangle, Wedge

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
__file__ = str(HERE / "qbo_tokamak.py")
exec(open(HERE / "qbo_tokamak.py").read().split("BOUNDS = ")[0])   # model, comb_drive, DATA, TG, ALIAS
X = np.array(json.load(open(HERE / "qbo_tokamak.json"))["COMB|threshold"]["x"])
INK, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
BLUE, ORANGE, GREEN, GOLD, PINK, RED = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#c8322b"
T0 = 1953.0


def arrow(ax, a, b, col=INK, lw=1.4, style="-|>", rad=0.0, ms=12):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle=style, mutation_scale=ms, color=col, lw=lw, connectionstyle=f"arc3,rad={rad}"))


def box(ax, xy, w, h, text, fc, ec=INK, fs=8.5, tc=INK):
    ax.add_patch(FancyBboxPatch(xy, w, h, boxstyle="round,pad=0.012,rounding_size=0.015", fc=fc, ec=ec, lw=1.1))
    ax.text(xy[0] + w / 2, xy[1] + h / 2, text, ha="center", va="center", fontsize=fs, color=tc, linespacing=1.25)


fig = plt.figure(figsize=(17, 12.5), facecolor=SURF)
gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.25], height_ratios=[1, 1.05], hspace=0.22, wspace=0.08)

# ---------------- (a) device overview: torus, TF coils, central solenoid, PF coils ----------------
ax = fig.add_subplot(gs[0, 0]); ax.set_xlim(-5, 5); ax.set_ylim(-4.3, 3.9); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("(a) Modified tokamak: axisymmetric torus, toroidal flow = QBO \"wind\"", loc="left", fontsize=11, color=INK)
R0, a0, tilt = 3.0, 1.05, 0.38
ax.add_patch(Ellipse((0, 0), 2 * (R0 + a0), 2 * (R0 + a0) * tilt, fc="#dbe9f8", ec=INK, lw=1.3, zorder=1))
ax.add_patch(Ellipse((0, 0), 2 * (R0 - a0), 2 * (R0 - a0) * tilt, fc=SURF, ec=INK, lw=1.3, zorder=2))
for ph in np.linspace(0, 2 * np.pi, 12, endpoint=False):           # TF coils as D-rings seen edge-on
    x, y = R0 * np.cos(ph), R0 * np.sin(ph) * tilt
    front = np.sin(ph) < 0
    ax.add_patch(Ellipse((x, y), 0.34 + 0.5 * abs(np.sin(ph)), 2.5 * a0 * (0.75 + 0.25 * abs(np.cos(ph))), fc="none",
                         ec=GOLD, lw=2.4 if front else 1.2, alpha=1 if front else 0.45, zorder=4 if front else 0))
ph_ = np.linspace(0, 2 * np.pi, 1400); th_ = 9 * ph_                   # helical (stellarator-like) winding, m > 0
xh = (R0 + 0.82 * a0 * np.cos(th_)) * np.cos(ph_); yh = (R0 + 0.82 * a0 * np.cos(th_)) * np.sin(ph_) * tilt + 0.82 * a0 * np.sin(th_) * 0.85
front_ = np.sin(ph_) < 0
ax.plot(np.where(front_, xh, np.nan), np.where(front_, yh, np.nan), color=PINK, lw=1.3, zorder=4)
ax.plot(np.where(~front_, xh, np.nan), np.where(~front_, yh, np.nan), color=PINK, lw=0.8, alpha=0.4, zorder=0)
ax.add_patch(Rectangle((-0.42, -2.4), 0.84, 4.8, fc="#f6d9c7", ec=ORANGE, lw=1.6, zorder=3))
for yy in np.linspace(-2.2, 2.2, 12):
    ax.plot([-0.42, 0.42], [yy, yy + 0.12], color=ORANGE, lw=0.9, zorder=3)
ax.text(-0.6, 2.75, "central solenoid (reversing poloidal drive)", ha="right", fontsize=8.5, color=ORANGE)
for yy, lab in ((2.05, "PF coil"), (-2.05, "PF coil")):
    ax.add_patch(Ellipse((0, yy), 2 * (R0 + 1.25), 2 * (R0 + 1.25) * tilt * 0.55, fc="none", ec=GREEN, lw=1.6, ls="--", zorder=0))
    ax.text(R0 + 1.32, yy, lab, fontsize=8, color=GREEN, va="center")
th = np.linspace(-0.25 * np.pi, 0.95 * np.pi, 80)                     # toroidal flow arrow (core +, edge reversed)
ax.plot((R0 - 0.35) * np.cos(th), (R0 - 0.35) * np.sin(th) * tilt, color=BLUE, lw=2.2, zorder=5)
arrow(ax, ((R0 - 0.35) * np.cos(th[-3]), (R0 - 0.35) * np.sin(th[-3]) * tilt), ((R0 - 0.35) * np.cos(th[-1]), (R0 - 0.35) * np.sin(th[-1]) * tilt), BLUE, 2.2, ms=16)
th2 = np.linspace(0.95 * np.pi, -0.25 * np.pi, 80)
ax.plot((R0 + 0.72) * np.cos(th2), (R0 + 0.72) * np.sin(th2) * tilt, color=RED, lw=2.2, zorder=5)
arrow(ax, ((R0 + 0.72) * np.cos(th2[-3]), (R0 + 0.72) * np.sin(th2[-3]) * tilt), ((R0 + 0.72) * np.cos(th2[-1]), (R0 + 0.72) * np.sin(th2[-1]) * tilt), RED, 2.2, ms=16)
ax.text(-4.9, -3.0, "blue: core toroidal field/flow (one sign)\nred: reversed edge layer (RFP-like), the newest QBO regime\n"
        "gold: toroidal-field coils (fixed, axisymmetric)\ngreen: poloidal-field coils   orange: drive solenoid\npink: helical coils (m > 0, high windings): fast ripple",
        fontsize=8.2, color=MUTED, va="top")
ax.text(4.95, 3.75, "toroidal = longitude around the equator (m = 0)\npoloidal = latitude-height plane\nedge = 10 hPa, core = lower stratosphere",
        fontsize=8.2, color=MUTED, ha="right", va="top")

# ---------------- (b) poloidal cross-section with reversal surfaces and probes ----------------
ax = fig.add_subplot(gs[1, 0]); ax.set_xlim(-1.55, 1.75); ax.set_ylim(-1.4, 1.45); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("(b) Poloidal cross-section: reversal surfaces move inward", loc="left", fontsize=11, color=INK)
layers = [(1.00, 0.86, RED), (0.86, 0.66, BLUE), (0.66, 0.45, RED), (0.45, 0.0, BLUE)]
for r1, r2, c in layers:
    ax.add_patch(Wedge((0, 0), r1, 0, 360, width=r1 - r2, fc=c, alpha=0.18, ec="none"))
for r in (1.0, 0.86, 0.66, 0.45):
    ax.add_patch(Circle((0, 0), r, fc="none", ec=INK if r == 1 else MUTED, lw=1.4 if r == 1 else 0.9, ls="-" if r == 1 else "--"))
for r, c in ((0.93, RED), (0.76, BLUE), (0.55, RED), (0.22, BLUE)):
    for ang in np.radians([100, 220, 340]):
        sym = "⊗" if c == RED else "⊙"
        ax.text(r * np.cos(ang), r * np.sin(ang), sym, ha="center", va="center", fontsize=12, color=c)
for r, L, c in ((0.90, 30, BLUE), (0.80, 50, GREEN), (0.75, 70, ORANGE)):
    ax.plot([r * np.cos(np.radians(20))], [r * np.sin(np.radians(20))], marker="o", ms=7, color=c, mec=INK, zorder=5)
    ax.annotate(f"probe r/a = {r:.2f}  (\"{L} hPa\")", (r * np.cos(np.radians(20)), r * np.sin(np.radians(20))),
                xytext=(1.08, 0.95 - 0.16 * [30, 50, 70].index(L)), fontsize=8.5, color=INK, arrowprops=dict(arrowstyle="-", color=c, lw=0.9))
for ang in np.radians([135, 150, 165]):
    arrow(ax, (1.18 * np.cos(ang), 1.18 * np.sin(ang)), (0.72 * np.cos(ang), 0.72 * np.sin(ang)), INK, 1.2, ms=10)
ax.text(-1.5, 1.32, "each new reversal surface is created at the edge\n(threshold on the drive) and penetrates inward:\n"
        "advection w (pinch) + diffusion η + damping τ_d\n= the QBO's downward propagation", fontsize=8.5, color=INK, va="top")
ax.text(-1.5, -1.05, "⊙ / ⊗  toroidal field (flow) out of / into the page\n---  reversal surfaces (B_φ = 0, \"zero-wind lines\")",
        fontsize=8.5, color=MUTED, va="top")

# ---------------- (c) reversing poloidal circuitry ----------------
ax = fig.add_subplot(gs[0, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
ax.set_title("(c) Reversing poloidal circuit (1 yr -> 1 s)", loc="left", fontsize=11, color=INK)
box(ax, (0.02, 0.70), 0.20, 0.20, "SOLAR CLOCK\n(strong carrier comb)\n2 Hz, alternating polarity\n= semiannual trigger,\nsign of solar declination", "#fde7c2")
box(ax, (0.02, 0.36), 0.20, 0.22, "LUNAR MODULATOR\n13.42 Hz (draconic)\nMoon's latitude above /\nbelow the clock's own\nplane (the ecliptic)", "#dbe9f8")
box(ax, (0.30, 0.52), 0.18, 0.24, "CAPACITOR BANK\ncharge per pulse\nV_n ∝ sin β_moon(t_n)\n× sin δ_sun(t_n)\n× (solar distance)³", "#eef0e6")
# H-bridge
hx, hy, hw, hh = 0.55, 0.47, 0.17, 0.34
ax.add_patch(Rectangle((hx, hy), hw, hh, fc="none", ec=INK, lw=1.1))
ax.text(hx + hw / 2, hy + hh + 0.025, "H-BRIDGE (thyristors)\nfires at each clock tick, polarity (−1)ⁿ", ha="center", fontsize=8.2, color=INK)
for (x0, y0), lab in (((hx, hy + hh * 0.72), "S1"), ((hx + hw, hy + hh * 0.72), "S2"), ((hx, hy + hh * 0.28), "S3"), ((hx + hw, hy + hh * 0.28), "S4")):
    ax.add_patch(Rectangle((x0 - 0.018, y0 - 0.03), 0.036, 0.06, fc="#fbe0dd", ec=RED, lw=1)); ax.text(x0, y0, lab, ha="center", va="center", fontsize=7.5)
ax.plot([hx + 0.03, hx + hw - 0.03], [hy + hh / 2, hy + hh / 2], color=ORANGE, lw=2)
ax.text(hx + hw / 2, hy + hh / 2 + 0.025, "solenoid primary", ha="center", fontsize=7.5, color=ORANGE)
box(ax, (0.79, 0.55), 0.19, 0.20, "L/R INTEGRATOR\nτ_i (leaky)\nloop voltage ->\nplasma current\n(poloidal field)", "#f6d9c7")
box(ax, (0.79, 0.22), 0.19, 0.20, "EDGE RELAXATION\n(RFP threshold ±θ,\nhysteresis)\nedge B_φ flips sign", "#fbe0dd")
box(ax, (0.50, 0.08), 0.26, 0.15, "INWARD PENETRATION\npinch w, diffusion η, damping τ_d\n-> probes at r/a = 0.90 / 0.80 / 0.75", "#e3f3ec")
box(ax, (0.25, 0.29), 0.21, 0.15, "HELICAL COIL DRIVER\n(stellarator-like, m > 0)\nhigh windings 8-26x of the\nmanifold, gated near the\nreversal surface", "#f9e1ec")
arrow(ax, (0.46, 0.365), (0.50, 0.365), PINK, 1.3)
box(ax, (0.50, 0.30), 0.22, 0.13, "PICKUP COILS +\nINTEGRATORS\nB_φ(t) at three radii\n= the \"QBO\" record", "#eef0e6")
arrow(ax, (0.22, 0.80), (0.30, 0.68)); arrow(ax, (0.22, 0.47), (0.30, 0.58)); arrow(ax, (0.48, 0.64), (0.55, 0.64))
arrow(ax, (0.12, 0.70), (0.12, 0.58), MUTED, 1.0, ms=9); ax.text(0.125, 0.635, "reference plane", fontsize=7.2, color=MUTED)
arrow(ax, (0.22, 0.86), (0.60, 0.835), RED, 1.0, rad=-0.15, ms=9); ax.text(0.38, 0.895, "gate + polarity", fontsize=7.5, color=RED)
arrow(ax, (0.72, 0.64), (0.79, 0.65)); arrow(ax, (0.885, 0.55), (0.885, 0.42)); arrow(ax, (0.79, 0.30), (0.76, 0.18)); arrow(ax, (0.61, 0.23), (0.61, 0.30))
# inset waveforms
ins = ax.inset_axes([0.0, 0.0, 0.40, 0.27]); ins.axis("off")
tt = np.arange(0, 12, 0.5); _, sds, dist, sb = astro(2000.0 + X[0] / 365.24 + tt)
pulses = sb * sds * dist
ins.vlines(tt, 0, pulses / np.abs(pulses).max(), color=BLUE, lw=1.8); ins.axhline(0, color=MUTED, lw=0.6)
D = comb_drive(X[0], X[6]); m = (TG >= 2000) & (TG < 2012)
ins.plot(TG[m] - 2000, D[m] / np.abs(D[m]).max(), color=ORANGE, lw=1.4)
ins.set_ylim(-1.5, 1.15)
ins.text(0, -1.25, "comb pulses (blue, 2 Hz, polarity (−1)ⁿ, height set by the lunar modulator)\nand the integrated drive (orange): the 2.37-s beat, 0–12 s", fontsize=7.5, color=MUTED, va="top")

# ---------------- (d) probe signals of the "shot" ----------------
ax = fig.add_subplot(gs[1, 1])
out = model("COMB", X, "threshold")[2]; s = TG >= T0
for k, (L, c, r) in enumerate(((30, BLUE, 0.90), (50, GREEN, 0.80), (70, ORANGE, 0.75))):
    t, y = DATA[L]; mm = (t >= T0) & (t < T0 + 24)
    sc = 1.0 if L != 50 else 100.0                               # qbo50.dat is stored in units of 100 m/s
    off = -60 * k
    ax.plot(t[mm] - T0, sc * y[mm] + off, color=c, lw=1.2, label=f"probe r/a = {r:.2f} (\"{L} hPa\")")
    mod = out[L][s][(TG[s] < T0 + 24)]; mt = TG[s][TG[s] < T0 + 24] - T0
    ax.plot(mt, 15 * (mod - mod.mean()) / mod.std() + off, color=INK, lw=0.9, ls="--", alpha=0.7)
    ax.text(24.3, off, f"r/a {r:.2f}", fontsize=8.5, color=c, va="center")
ax.set_xlim(0, 25.8); ax.set_xlabel("time in the shot (s)   [1 s = 1 yr of the 1953- record]", color=MUTED)
ax.set_ylabel("edge toroidal field B_φ (mT, offset per probe)", color=MUTED)
ax.set_title("(d) QBO record as probe signals (first 25 s): reversals every ~2.37 s, outer probe first",
             loc="left", fontsize=11, color=INK)
ax.plot([], [], color=INK, ls="--", lw=0.9, label="device (comb drive + threshold), scaled")
ax.legend(frameon=False, fontsize=8.5, loc="upper right", ncol=2, bbox_to_anchor=(1.0, 1.0))
ax.set_ylim(-170, 55)
ax.set_facecolor(SURF); ax.grid(True, color=GRID, lw=0.5); ax.tick_params(colors=MUTED)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
fig.savefig(HERE / "qbo_tokamak_device.png", dpi=95, facecolor=SURF, bbox_inches="tight")

# ---------------- second figure: whole shot, time-radius map ----------------
fig2, axs = plt.subplots(2, 1, figsize=(16, 7.5), facecolor=SURF, gridspec_kw=dict(height_ratios=[1, 1]))
tg = np.arange(T0, 2024.2, 1 / 12); rr = np.array([0.90, 0.80, 0.75])
obs = np.array([np.interp(tg, *DATA[L]) * (100 if L == 50 else 1) for L in (30, 50, 70)])
mod = np.array([np.interp(tg, TG, out[L]) for L in (30, 50, 70)]); mod = 25 * (mod - mod.mean(1, keepdims=True)) / mod.std(1, keepdims=True)
for a_, z, lab in ((axs[0], obs, "measured (KIT QBO winds read as B_φ, mT)"), (axs[1], mod, "device (comb + threshold), scaled to 25 mT rms")):
    im = a_.pcolormesh(tg - T0, rr, z, cmap="RdBu_r", vmin=-35, vmax=35, shading="gouraud")
    a_.set_xlim(0, 25); a_.set_ylim(0.75, 0.90); a_.set_ylabel("probe radius r/a", color=MUTED); a_.set_title(lab, loc="left", fontsize=10.5, color=INK)
    a_.set_facecolor(SURF); a_.tick_params(colors=MUTED)
    fig2.colorbar(im, ax=a_, pad=0.01, label="B_φ (mT)")
axs[1].set_xlabel("time in the shot (s)   [1 s = 1 yr, t = 0 at 1953]", color=MUTED)
fig2.suptitle("First 25 s of the 73-s shot: edge toroidal field between r/a = 0.90 and 0.75. Reversal fronts tilt to later times inward "
              "(inward penetration ≈ 0.6 s from 0.90 to 0.75)", x=0.01, ha="left", fontsize=11.5, color=INK)
fig2.tight_layout(); fig2.savefig(HERE / "qbo_tokamak_shot.png", dpi=95, facecolor=SURF)
print("plots: qbo_tokamak_device.png, qbo_tokamak_shot.png")
