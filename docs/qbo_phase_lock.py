#!/usr/bin/env python3
"""Do QBO reversals lock to the 2.370-yr draconic/semiannual beat? (2026-10-22)

Reversal dates: KIT monthly winds at 30/50/70 hPa, 3-month running mean, hysteresis +-3 m/s (a reversal is counted when
the wind passes the opposite threshold; its date is the zero crossing in between, interpolated).
Phase: phi = frac((t - 1953) * f) for westerly onsets, frac(... ) + 0.5 for easterly onsets; Rayleigh R = |mean e^{2 pi i phi}|.
Tests
  1. frequency scan 0.30-0.65 cycles/yr: where does R peak? predicted 0.4219 (2.370 yr); controls 1/1.73, 1/2.72
  2. interval-shuffle null: same set of half-cycle intervals in random order (same period spread, phase free to wander);
     R at 0.4219 and the scan maximum, real vs 5000 surrogates
  3. phase drift: the circular mean phase in 1953-1986, 1987-2005, 2006-2024
    python3 qbo_phase_lock.py  -> qbo_phase_lock.png, stdout
"""
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
LEGS = HERE.parent / "experiments/Oct2026_generator/legs"
FILES = {30: (LEGS / "qbo30_point0E/kit_qbo_u30.dat", 1.0), 50: (LEGS / "qbo50/qbo50.dat", 100.0), 70: (LEGS / "qbo70_point0E/kit_qbo_u70.dat", 1.0)}
F0 = 365.2422 / 27.212221 - 13 - 1 + 1    # placeholder, replaced below
F0 = (365.2422 / 27.212221) % 1.0         # 0.4219 cycles/yr
CTRL = {"2.370 yr (predicted)": F0, "1.73 yr (no polarity flip)": 1 - F0, "2.72 yr (declination)": 365.2422 / 27.321582 % 1.0}
INK, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
rng = np.random.default_rng(11)


def reversals(path, scale, thr=3.0):
    a = np.loadtxt(path); t, u = a[:, 0], a[:, 1] * scale
    m = (t >= 1953) & (t < 2024.3); t, u = t[m], np.convolve(u[m], np.ones(3) / 3, "same")
    out = []; state = np.sign(u[0]) or 1.0; last0 = None
    for i in range(1, len(u)):
        if np.sign(u[i]) != np.sign(u[i - 1]):
            last0 = t[i - 1] + (t[i] - t[i - 1]) * u[i - 1] / (u[i - 1] - u[i])
        if state < 0 and u[i] > thr:
            out.append((last0, +1)); state = 1.0
        elif state > 0 and u[i] < -thr:
            out.append((last0, -1)); state = -1.0
    return np.array([o[0] for o in out]), np.array([o[1] for o in out])


def R(t, s, f):
    ph = (t - 1953.0) * f + np.where(s > 0, 0.0, 0.5)
    z = np.mean(np.exp(2j * np.pi * ph)); return abs(z), np.angle(z) / (2 * np.pi) % 1


FS = np.arange(0.30, 0.65, 0.0005)
REV = {L: reversals(*FILES[L]) for L in FILES}
res = {}
for L, (t, s) in REV.items():
    d = np.diff(t)
    print(f"\n{L} hPa: {len(t)} reversals 1953-2024, half-cycle intervals {d.mean():.2f} +- {d.std():.2f} yr -> mean period {2 * d.mean():.3f} yr")
    scan = np.array([R(t, s, f)[0] for f in FS]); fpk = FS[np.argmax(scan)]
    print(f"   frequency scan: R peaks at {fpk:.4f} cycles/yr ({1 / fpk:.3f} yr), R = {scan.max():.3f}")
    for k, f in CTRL.items():
        print(f"   R at {k:28s} {R(t, s, f)[0]:.3f}")
    # interval-shuffle surrogates
    S0, Sm = [], []
    for _ in range(5000):
        tt = t[0] + np.concatenate([[0.0], np.cumsum(rng.permutation(d))])
        ss = s[0] * (-1.0) ** np.arange(len(tt))
        S0.append(R(tt, ss, F0)[0]); Sm.append(max(R(tt, ss, f)[0] for f in FS[::8]))
    S0, Sm = np.array(S0), np.array(Sm); r0 = R(t, s, F0)[0]
    print(f"   shuffle null: R(0.4219) real {r0:.3f} vs surrogates {np.median(S0):.3f} (95th {np.percentile(S0, 95):.3f}); p = {np.mean(S0 >= r0):.4f}")
    print(f"                 scan max real {scan.max():.3f} vs surrogate scan max {np.median(Sm):.3f} (95th {np.percentile(Sm, 95):.3f}); p = {np.mean(Sm >= scan.max()):.4f}")
    for a, b in ((1953, 1987), (1987, 2006), (2006, 2025)):
        m = (t >= a) & (t < b); r, ph = R(t[m], s[m], F0)
        print(f"   {a}-{b - 1}: {m.sum():2d} reversals, R {r:.2f}, mean phase {ph:.2f} cycles")
    res[L] = (t, s, scan, S0, r0)

fig, axs = plt.subplots(1, 3, figsize=(16, 4.8), facecolor=SURF)
cols = {30: "#2a78d6", 50: "#1baf7a", 70: "#eb6834"}
ax = axs[0]
for L, (t, s, scan, S0, r0) in res.items():
    ax.plot(FS, scan, color=cols[L], lw=1.4, label=f"{L} hPa")
for k, f in CTRL.items():
    ax.axvline(f, color=MUTED if "predicted" not in k else INK, lw=1.0, ls="-" if "predicted" in k else ":")
    ax.text(f, 1.0, k.split(" (")[0], rotation=90, va="top", ha="right", fontsize=8, color=MUTED)
ax.set_xlabel("beat frequency (cycles/yr)", color=MUTED); ax.set_ylabel("Rayleigh R of reversal phases", color=MUTED); ax.set_ylim(0, 1.02)
ax.set_title("phase clustering of reversals vs beat frequency", loc="left", fontsize=10.5, color=INK); ax.legend(frameon=False, fontsize=8.5)
ax = axs[1]
for L, (t, s, scan, S0, r0) in res.items():
    ax.hist(S0, bins=50, color=cols[L], alpha=0.35, density=True); ax.axvline(r0, color=cols[L], lw=2, label=f"{L} hPa real")
ax.set_xlabel("R at 0.4219 cycles/yr", color=MUTED); ax.set_title("real vs interval-shuffled reversal sequences", loc="left", fontsize=10.5, color=INK)
ax.legend(frameon=False, fontsize=8.5)
ax = axs[2]
for L, (t, s, scan, S0, r0) in res.items():
    ph = ((t - 1953.0) * F0 + np.where(s > 0, 0.0, 0.5)) % 1
    ax.scatter(t, ph, color=cols[L], s=18, label=f"{L} hPa")
ax.set_ylim(0, 1); ax.set_xlabel("year of reversal", color=MUTED); ax.set_ylabel("phase in the 2.370-yr beat (cycles)", color=MUTED)
ax.set_title("reversal phase over time: flat = locked", loc="left", fontsize=10.5, color=INK); ax.legend(frameon=False, fontsize=8.5)
for a in axs:
    a.set_facecolor(SURF); a.tick_params(colors=MUTED); a.grid(True, color=GRID, lw=0.5)
    for sp in ("top", "right"):
        a.spines[sp].set_visible(False)
fig.tight_layout(); fig.savefig(HERE / "qbo_phase_lock.png", dpi=95, facecolor=SURF)
print("\nplot: qbo_phase_lock.png")
