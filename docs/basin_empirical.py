#!/usr/bin/env python3
"""Empirical counterpart of basin.py: Landreman-type field lines driven by the fitted GEM-LTE manifolds (2026-10-21).

Landreman (arXiv:2609.26742), family 1: along a field line, parametrized by zeta, the Cartesian coordinates obey
decoupled harmonic oscillators  x'' = -x, y'' = -y, z'' = -4z  (zeta-derivatives), so a field line is
    x = a cos(zeta),  y = b sin(zeta),  z = c cos(2 zeta + beta)
and the rotational transform is the frequency ratio, iota = 2.

GEM-LTE (Chapter 12, eq. 12-14): the response zeta(t) = sin(sqrt(A) M(t) + theta0) is, written in the manifold M
instead of time, exactly a harmonic oscillator in M:  d2(zeta)/dM2 = -A zeta.  The fitted model is a sum of such
oscillators, sin(2 pi k F(t) + phase), with windings k = h * k_base (integer harmonics h) plus slow windings, all on one
manifold F(t) (lte_results.csv column 4). Taking the manifold as the field-line parameter,
    zeta = 2 pi k_base F(t),
the base winding gives the (x, y) pair and a harmonic h gives z = cos(h zeta + phase): Landreman's family-1 form with
the 2:1 ratio replaced by the fitted h:1, i.e. iota = h. Amplitudes and phases are the fitted ones
(lt.exe.windings.N.json, k_amp_phase). Because h is an integer the curve closes (a rational surface); a slow winding
with a non-integer ratio to the base (k_slow/k_base = 0.581 for the common manifold) instead fills a torus
(an irrational surface), shown in the last 3D panel.

Note: the integer ratios are imposed by the fitted model (integer harmonics of one base), not discovered; what the data
choose is which harmonic dominates and the non-integer slow-winding ratio.

    python3 basin_empirical.py      -> basin_empirical.png
"""
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

HERE = Path(__file__).resolve().parent
LEGS = HERE.parent / "experiments/Oct2026_generator/legs"
INK, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
TIME = LinearSegmentedColormap.from_list("time", ["#cde2fb", "#3987e5", "#104281"])      # sequential, one hue

# (title, run directory, leg, harmonic h for z, manifold label)
CASES = [("NINO4, common manifold", "nino4_34_3_12", 1, 9, "current tides"),
         ("Baltic MSL, pre-fix manifold", "nao_baltic_legacy_balticseed_2000_2005", 1, 6, "pre-fix tides"),
         ("NAO, pre-fix manifold", "nao_baltic_legacy_naoseed_2000_2005", 2, 4, "pre-fix tides"),
         ("QBO 30 hPa (KIT), QBO manifold", "qbo30_point0E", 1, 8, "pre-fix tides")]
SLOW_CASE = ("NINO4, base and slow winding", "amo_pdo_nino4_iode", 3)


def load(run, leg):
    d = LEGS / run
    p = json.load(open(d / ("result.p" if (d / "result.p").exists() else "lt.exe.p")))
    a = np.loadtxt(d / f"lte_results.{leg}.csv", delimiter=",")
    W = json.load(open(d / f"lt.exe.windings.{leg}.json"))["k_amp_phase"]
    return p, a[:, 0], a[:, 3], W


def term(W, k):
    """amplitude and phase of the fitted term with |winding| closest to k"""
    w = min(W, key=lambda x: abs(abs(x[0]) - k))
    return abs(w[0]), w[1], w[2] * np.sign(w[0])


def style3d(ax, title):
    ax.set_title(title, fontsize=10, color=INK, loc="left")
    ax.set_facecolor(SURF); ax.tick_params(colors=MUTED, labelsize=7)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_pane_color((0.99, 0.99, 0.98, 1.0)); axis._axinfo["grid"]["color"] = GRID


fig = plt.figure(figsize=(17, 13), facecolor=SURF)
gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 0.55])

# (a) Landreman family 1, exact: iota = 2
ax = fig.add_subplot(gs[0, 0], projection="3d"); z = np.linspace(0, 2 * np.pi, 400)
ax.plot(np.cos(z), np.sin(z), 0.5 * np.cos(2 * z), color=INK, lw=1.5)
style3d(ax, "(a) Landreman family 1 field line:\nx=cos ζ, y=sin ζ, z=½cos 2ζ  (ι = 2, closed)")

summary = []
for i, (title, run, leg, h, man) in enumerate(CASES):
    p, t, F, W = load(run, leg); kb = abs(p["ltep"][-1])
    _, ab, pb = term(W, kb); kh, ah, ph = term(W, h * kb)
    zeta = 2 * np.pi * kb * F
    x, y, zz = np.cos(zeta + pb), np.sin(zeta + pb), (ah / ab) * np.cos(2 * np.pi * kh * F + ph)
    ax = fig.add_subplot(gs[(i + 1) // 3, (i + 1) % 3], projection="3d")
    sc = ax.scatter(x, y, zz, c=t, cmap=TIME, s=2, depthshade=False)
    style3d(ax, f"({'bcde'[i]}) {title}: z = {kh / kb:.0f}× harmonic\n(ι = {kh / kb:.0f}, closed; amplitude ratio {ah / ab:.2f}; {man})" + (f"\nmanifold spans only {kb * (F.max() - F.min()):.1f} turn: a segment of the line" if kb * (F.max() - F.min()) < 1 else ""))
    summary.append((title, t, F, kb, kh / kb, ah / ab))

# (f) irrational surface: base winding (toroidal) and the slow winding (poloidal) on a torus
p, t, F, W = load(*SLOW_CASE[1:]); kb = abs(p["ltep"][-1]); ks = abs(p["ltep"][0]); iota_slow = ks / kb
_, _, pb = term(W, kb); _, _, ps = term(W, ks)
phi, th = 2 * np.pi * kb * F + pb, 2 * np.pi * ks * F + ps; R0, r = 1.0, 0.4
ax = fig.add_subplot(gs[1, 2], projection="3d")
sc = ax.scatter((R0 + r * np.cos(th)) * np.cos(phi), (R0 + r * np.cos(th)) * np.sin(phi), r * np.sin(th), c=t, cmap=TIME, s=2, depthshade=False)
style3d(ax, f"(f) {SLOW_CASE[0]}: toroidal = base, poloidal = slow\n(ι = {ks / kb:.3f}, irrational: the line fills the torus)")
cb = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.08); cb.set_label("year", color=MUTED); cb.ax.tick_params(colors=MUTED)

# (g) the field-line parameters themselves: the fitted manifolds F(t)
ax = fig.add_subplot(gs[2, :]); ax.set_facecolor(SURF)
for (title, t, F, kb, h, ratio), col in zip(summary, ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]):
    ax.plot(t, kb * F, color=col, lw=0.9, label=f"{title}  (base {kb:.4f})")
ax.set_ylabel("k_base · F(t)  (field-line turns)", color=MUTED); ax.set_xlabel("year", color=MUTED)
ax.set_title("(g) The field-line parameter of each panel: the fitted manifold, in turns of the base winding "
             "(travel back and forth along the line is the non-autonomous, tidally paced part)", loc="left", fontsize=10, color=INK)
ax.grid(True, color=GRID, lw=0.6); ax.tick_params(colors=MUTED); ax.legend(frameon=False, fontsize=8, ncol=2, loc="lower left")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.suptitle("Landreman-type field lines driven by the fitted GEM-LTE manifolds (colour = year)", x=0.01, ha="left", fontsize=13, color=INK)
fig.tight_layout()
fig.savefig(HERE / "basin_empirical.png", dpi=100, facecolor=SURF)
print("plot: basin_empirical.png")
for title, t, F, kb, h, ratio in summary:
    print(f"  {title:32s} base {kb:.4f}, iota = {h:.0f}, z/x amplitude ratio {ratio:.2f}, manifold spans {kb * (F.max() - F.min()):.1f} turns")
print(f"  {SLOW_CASE[0]:32s} slow/base iota = {iota_slow:.4f}")
