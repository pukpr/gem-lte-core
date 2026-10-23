#!/usr/bin/env python3
"""Absolute-phase test: do QBO reversals follow the astronomically fixed beat? (2026-10-22)

Drive, nothing fitted: impulses at the June and December solstices (day 172 and half a year later), heights
sin(beta_moon) * sin(dec_sun) * (solar distance)^3 from Meeus' mean elements (astro() in qbo_tokamak.py).
Its beat phase psi(t): complex demodulation of the impulse series at f0 = 0.4219 cycles/yr, Gaussian-smoothed over
~2 beat cycles, psi = 2 pi f0 t + theta(t). theta(t) wobbles with the lunar orbit terms.
Reversals: as in qbo_phase_lock.py (30/50/70 hPa). Phase of reversal k: psi(t_k) + pi for easterly onsets.
  1. implied lag: the circular mean offset between reversals and the drive's phase, in months (modulo the beat; the
     sign convention of the drive is not fixed physically, so the half-cycle alternative is also given)
  2. lag-independent test: R with psi(t) vs R with the plain 2.370-yr carrier; null = psi(t + Delta) for random
     Delta (3-40 yr either way), which keeps the wobble statistics but misaligns them (R ignores a constant offset)
    python3 qbo_abs_phase.py  -> qbo_abs_phase.png, stdout
"""
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
__file__ = str(HERE / "qbo_tokamak.py")
exec(open(HERE / "qbo_tokamak.py").read().split("DATA = {}")[0])            # astro()
_pl = open(HERE / "qbo_phase_lock.py").read()
exec(_pl.split("FS = np.arange")[0].split("INK, MUTED")[0])                 # FILES, F0
exec("def reversals" + _pl.split("def reversals")[1].split("def R(")[0])
INK, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
rng = np.random.default_rng(21)
P0 = 1 / F0

tn = np.concatenate([[y + 171.5 / 365.2422, y + (171.5 + 182.62) / 365.2422] for y in range(1880, 2070)])
_, sds, dist, sb = astro(tn); a = sb * sds * dist
tc = np.arange(1885.0, 2065.0, 1 / 12)
z = a * np.exp(-2j * np.pi * F0 * tn)
sig = 2 * P0 / 2.355
W = np.exp(-0.5 * ((tc[:, None] - tn[None, :]) / sig) ** 2)
zc = (W @ z) / W.sum(1)
theta = np.unwrap(np.angle(zc))
psi = lambda t, d=0.0: 2 * np.pi * F0 * t + np.interp(t + d, tc, theta)  # noqa: E731
print(f"drive beat-phase wobble theta(t): range {np.ptp(theta[(tc > 1953) & (tc < 2025)]) / (2 * np.pi):.3f} cycles over 1953-2024, "
      f"sd {np.std(theta[(tc > 1953) & (tc < 2025)]) / (2 * np.pi):.3f} cycles")
amp = np.abs(zc)
print(f"drive beat amplitude over 1953-2024: min/max {amp[(tc > 1953) & (tc < 2025)].min() / amp.mean():.2f} / {amp[(tc > 1953) & (tc < 2025)].max() / amp.mean():.2f} of mean")

R = lambda ph: np.mean(np.exp(1j * ph))  # noqa: E731
for L in FILES:
    t, s = reversals(*FILES[L]); half = np.where(s > 0, 0.0, np.pi)
    zl = R(2 * np.pi * F0 * t + half); za = R(psi(t) + half)
    off = (-np.angle(za) / (2 * np.pi)) % 1
    nul = []
    for _ in range(4000):
        d = rng.uniform(3, 40) * rng.choice([-1, 1]); nul.append(abs(R(psi(t, d) + half)))
    nul = np.array(nul)
    print(f"\n{L} hPa ({len(t)} reversals): R plain 2.370-yr carrier {abs(zl):.3f}; R astronomical phase {abs(za):.3f}; "
          f"time-shifted ephemeris {np.median(nul):.3f} (95th {np.percentile(nul, 95):.3f}), percentile of real {np.mean(nul < abs(za)) * 100:.1f}")
    print(f"   implied lag of westerly onsets behind the drive's beat (+ phase): {off * P0 * 12:5.1f} months "
          f"(or {((off + 0.5) % 1) * P0 * 12:5.1f} months with the drive's sign reversed); beat = {P0 * 12:.1f} months")
    for a0, b0 in ((1953, 1987), (1987, 2006), (2006, 2025)):
        m = (t >= a0) & (t < b0)
        print(f"   {a0}-{b0 - 1}: R plain {abs(R(2 * np.pi * F0 * t[m] + half[m])):.2f}  R astronomical {abs(R(psi(t[m]) + half[m])):.2f}")

fig, axs = plt.subplots(1, 2, figsize=(15, 4.8), facecolor=SURF)
ax = axs[0]; m = (tc > 1950) & (tc < 2026)
ax.plot(tc[m], (theta[m] - theta[m].mean()) / (2 * np.pi), color=INK, lw=1.6, label="drive phase wobble (ephemeris)")
cols = {30: "#2a78d6", 50: "#1baf7a", 70: "#eb6834"}
for L in FILES:
    t, s = reversals(*FILES[L]); ph = (F0 * t + np.where(s > 0, 0, 0.5)) % 1
    ph = np.angle(np.exp(2j * np.pi * (ph - np.angle(R(2 * np.pi * ph)) / (2 * np.pi)))) / (2 * np.pi)
    ax.scatter(t, -ph, s=14, color=cols[L], alpha=0.7, label=f"{L} hPa reversals (phase vs 2.370-yr carrier, centred)")
ax.set_ylabel("cycles (later = up)", color=MUTED); ax.set_title("does reversal timing follow the ephemeris wobble?", loc="left", fontsize=10.5, color=INK)
ax.legend(frameon=False, fontsize=8)
ax = axs[1]; ax.plot(tc[m], amp[m] / amp.mean(), color=INK, lw=1.4)
ax.set_title("drive beat amplitude (ephemeris), relative to mean", loc="left", fontsize=10.5, color=INK)
for a_ in axs:
    a_.set_facecolor(SURF); a_.tick_params(colors=MUTED); a_.grid(True, color=GRID, lw=0.5)
    for sp in ("top", "right"):
        a_.spines[sp].set_visible(False)
fig.tight_layout(); fig.savefig(HERE / "qbo_abs_phase.png", dpi=95, facecolor=SURF)
print("\nplot: qbo_abs_phase.png")
