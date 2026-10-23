#!/usr/bin/env python3
"""QBO emulator = axisymmetric device (slow, wavenumber 0) + helical "stellarator" coil set (fast, high windings,
wavenumber > 0), compared with ERA5 5-day winds at 30 hPa by characteristics, not point by point (2026-10-22).

Slow part   the saturating-core + gated-penetration device of qbo_shoulders.py (SATWARP), at 30 hPa
Fast part   the high-winding content of the two-leg QBO fit (legs/qbo30_point0E, leg 2 = ERA5 5-day, windings up to 26x
            the base on the common manifold): its model output minus a ~4-month running mean (periods under ~4 months;
            a 1-yr cut would carry the observed reversal edges into the emulator)
Gate        the helical ripple's amplitude depends on the slow state, as in a reversed-field pinch where the helical
            modes concentrate at the reversal surface:  G = 1 + g0 exp(-(u_slow/u0)^2) + g1 * max(du_slow/dt, 0)
            (g0, g1 set from the data's own state dependence, not fitted to the series)
Characteristics compared (data vs emulator, with and without the gate):
  1. power spectrum of the 5-day series (log-spectrum correlation, 0.1-30 cycles/yr)
  2. fast-ripple amplitude by state: westerly / easterly / transition / westerly-accelerating / easterly-accelerating
  3. reversal timing: mean period and the seasonal clustering of transitions
    python3 qbo_stellarator.py  -> qbo_stellarator.png, stdout
"""
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
_src = open(HERE / "qbo_shoulders.py").read().split("B = {")[0]
__file__ = str(HERE / "qbo_shoulders.py"); exec(_src)                       # run(), TG (5-day grid), DATA, astro ...
LEG = HERE.parent / "experiments/Oct2026_generator/legs/qbo30_point0E"
INK, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
BLUE, ORANGE, GREEN, PINK = "#2a78d6", "#eb6834", "#1baf7a", "#e87ba4"

a = np.loadtxt(LEG / "lte_results.2.csv", delimiter=","); t5, m5, y5 = a[:, 0], a[:, 1], a[:, 2]
n = 73; box = np.ones(n) / n
lp = lambda v: np.convolve(v, box, "same")  # noqa: E731
box3 = np.ones(25) / 25                                                      # ~4-month window
hp3 = lambda v: v - np.convolve(v, box3, "same")  # noqa: E731
fast_fit = hp3(m5)                                                           # helical (high-winding) content, periods < ~4 months
x = np.array(json.load(open(HERE / "qbo_shoulders.json"))["SATWARP"]["x"])
slow = np.interp(t5, TG, run("SATWARP", x)[1][30])
ys = lp(y5); slow = (slow - slow.mean()) / slow.std() * ys[n:-n].std() + ys[n:-n].mean()   # scale to the observed slow wind


def gate(us, g0, g1, u0=5.0):
    du = np.gradient(us, t5); du = du / np.percentile(np.abs(du), 90)
    G = 1 + g0 * np.exp(-(us / u0) ** 2) + g1 * np.maximum(du, 0); return G / G.mean()


# set the gate from the DATA's state dependence (fast amplitude by state, normalized)
yf = hp3(y5)
def states(us, f):
    du = np.gradient(us, t5); e = np.abs(f); k = slice(n, -n)
    S = {"westerly": us > 5, "easterly": us < -5, "transition": np.abs(us) < 5, "W-accelerating": du > 0, "E-accelerating": du < 0}
    return {s: e[k][v[k]].mean() for s, v in S.items()}


obs = states(ys, yf)
g0 = obs["transition"] / obs["westerly"] - 1; g1 = obs["W-accelerating"] / obs["E-accelerating"] - 1
EMU = {"slow only (axisymmetric)": slow + 0 * fast_fit,
       "slow + helical, no gate": slow + fast_fit,
       "slow + helical, gated": slow + fast_fit * gate(slow, g0, g1)}


def spectrum(v):
    v = v[n:-n] - v[n:-n].mean(); w = np.hanning(len(v)); P = np.abs(np.fft.rfft(v * w)) ** 2
    f = np.fft.rfftfreq(len(v), 1 / 73.0); P = np.convolve(P, np.ones(9) / 9, "same"); return f, P


fo, Po = spectrum(y5); band = (fo > 0.1) & (fo < 30)
print(f"gate set from data: g0 (transition boost) {g0:.2f}, g1 (westerly-acceleration boost) {g1:.2f}")
print(f"{'':28s} {'log-spectrum r':>14s} | fast |u'| by state (m/s): " + "  ".join(f"{s:>14s}" for s in obs))
print(f"{'observed':28s} {'':>14s} | " + "  ".join(f"{v:14.2f}" for v in obs.values()))
for k, e in EMU.items():
    f, P = spectrum(e); r = np.corrcoef(np.log(P[band] + 1e-9), np.log(Po[band]))[0, 1]
    st = states(lp(e), hp3(e))
    print(f"{k:28s} {r:14.3f} | " + "  ".join(f"{v:14.2f}" for v in st.values()))
for k, e in (("observed", y5), ("emulator (gated)", EMU["slow + helical, gated"])):
    s = np.sign(lp(e)[n:-n]); c = np.where(np.diff(s) != 0)[0]; tt = t5[n:-n][c]
    mo = np.round((tt % 1) * 12).astype(int) % 12; hist = np.bincount(mo, minlength=12)
    print(f"{k:18s}: mean reversal period {2 * np.mean(np.diff(tt)):.2f} yr; reversals by month J..D: {' '.join(str(h) for h in hist)}")

fig = plt.figure(figsize=(16, 10.5), facecolor=SURF); gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.05])
for row, (lab, e, col) in enumerate((("observed ERA5 30 hPa (5-day)", y5, INK), ("emulator: axisymmetric device + gated helical coils", EMU["slow + helical, gated"], BLUE))):
    ax = fig.add_subplot(gs[row, :]); m = (t5 >= 1995) & (t5 < 2006)
    ax.plot(t5[m], e[m], color=col, lw=0.9); ax.plot(t5[m], lp(e)[m], color=ORANGE, lw=0.3, alpha=0.9)
    ax.axhline(0, color=MUTED, lw=0.6); ax.set_xlim(1995, 2006); ax.set_ylabel("u (m/s)  [B_φ in mT]", color=MUTED)
    ax.set_title(f"{lab}, 1995-2005 (orange: 1-yr running mean = the axisymmetric part)", loc="left", fontsize=10.5, color=INK)
ax = fig.add_subplot(gs[2, 0])
ax.loglog(fo[band], Po[band] / Po[band].max(), color=INK, lw=1.4, label="observed")
for (k, e), col in zip(EMU.items(), (ORANGE, GREEN, BLUE)):
    f, P = spectrum(e); ax.loglog(f[band], P[band] / Po[band].max(), color=col, lw=1.1, alpha=0.85, label=k)
ax.set_xlabel("cycles per year  [= Hz in the shot]", color=MUTED); ax.set_title("power spectrum (5-day sampling)", loc="left", fontsize=10.5, color=INK)
ax.legend(frameon=False, fontsize=8)
ax = fig.add_subplot(gs[2, 1]); S = list(obs); xx = np.arange(len(S)); wdt = 0.2
ax.bar(xx - 1.5 * wdt, [obs[s] for s in S], wdt, color=INK, label="observed")
for j, ((k, e), col) in enumerate(zip(EMU.items(), (ORANGE, GREEN, BLUE))):
    st = states(lp(e), hp3(e)); ax.bar(xx + (j - 0.5) * wdt, [st[s] for s in S], wdt, color=col, label=k)
ax.set_xticks(xx); ax.set_xticklabels(S, fontsize=9); ax.set_ylabel("mean |fast u'| (m/s), periods < 4 months", color=MUTED)
ax.set_title("fast ripple by slow state: strongest at reversals and westerly onset", loc="left", fontsize=10.5, color=INK)
ax.legend(frameon=False, fontsize=8)
for a_ in fig.axes:
    a_.set_facecolor(SURF); a_.tick_params(colors=MUTED); a_.grid(True, color=GRID, lw=0.5)
    for sp in ("top", "right"):
        a_.spines[sp].set_visible(False)
fig.tight_layout(); fig.savefig(HERE / "qbo_stellarator.png", dpi=95, facecolor=SURF)
print("plot: qbo_stellarator.png")
