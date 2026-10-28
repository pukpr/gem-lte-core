#!/usr/bin/env python3
"""Spectral decomposition of the QBO: which tidal lines are present beyond the draconic beat? (2026-10-23)

Bookkeeping as in ChandlerWobble/docs (cw_data_tests.md): a lunar term of frequency f (cycles/yr) sampled by the annual /
semiannual comb folds to |n - f|; the Chandler wobble is 2 x draconic folded (0.8438) and Mm appears at 14 - anomalistic
(0.7449).
A. Folded lines in the monthly KIT winds, 30/50/70 hPa, 1953-2026 (monthly climatology removed).
B. Direct lunar lines (no folding) in daily ERA5 winds: 30 hPa (0 E), 50 hPa (zonal band), 70 hPa (to 2020); a smooth
   6-harmonic seasonal fit and a 1-yr running mean removed.
Each candidate: least-squares amplitude at its exact frequency; percentile among neighbouring frequencies (excluding
the line itself); phase in the first vs second half of the record (a forced line keeps its phase).
Then a free peak search listing the strongest peaks with the nearest candidate.
    python3 qbo_tidal_lines.py  -> qbo_tidal_lines.png, qbo_tidal_lines.json, stdout
"""
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
G = HERE.parent / "experiments/Oct2026_generator"
INK, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
Y = 365.2422
fD, fA, fT = Y / 27.212221, Y / 27.554550, Y / 27.321582            # draconic, anomalistic (Mm), tropical
fEv = Y / 31.811939                                                  # evection
NOD = 1 / 18.613


def fold(f, n): return abs(n - f)


CAND_A = {"QBO: draconic (13 - fD)": fold(fD, 13), "2 x draconic (CW line; also QBO 2nd harmonic)": fold(2 * fD, 27),
          "draconic + nodal": fold(fD, 13) + NOD, "Mm: anomalistic - 13": fold(fA, 13), "Mm: 14 - anomalistic": fold(fA, 14),
          "Mm, semiannual image (15 - fA)": fold(fA, 15), "tropical - 13": fold(fT, 13), "14 - tropical": fold(fT, 14),
          "Mf (2 x tropical) 27 - 2fT": fold(2 * fT, 27), "2Mm (2 x anomalistic) 27 - 2fA": fold(2 * fA, 27),
          "2Mm, 26 - ... ": fold(2 * fA, 26), "evection 11.48 -> 0.48": fold(fEv, 11), "evection 12 - fEv": fold(fEv, 12),
          "draconic semiannual image (15 - fD)": fold(fD, 15)}
CAND_B = {"Mm (anomalistic)": fA, "tropical month": fT, "draconic month": fD, "Mf (2 x tropical)": 2 * fT,
          "2 x draconic": 2 * fD, "2Mm (2 x anomalistic)": 2 * fA, "evection": fEv, "Msf (Mf - 2)": 2 * fT - 2,
          "Mt (3 x tropical - ...)": Y / 9.132948}


def prep(t, y, hp=False):
    m = np.isfinite(y); t, y = t[m], y[m]
    if hp:      # daily: smooth 6-harmonic seasonal fit (a stepwise monthly climatology leaves lines at integer cycles/yr)
        X = np.column_stack([np.ones_like(t)] + [g(2 * np.pi * k * t) for k in range(1, 7) for g in (np.cos, np.sin)])
        y = y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
    else:
        mon = np.floor((t % 1) * 12).astype(int) % 12
        y = y - np.array([y[mon == k].mean() for k in range(12)])[mon]
    if hp:
        n = int(round(1 / np.median(np.diff(t)))); y = y - np.convolve(y, np.ones(n) / n, "same"); t, y = t[n:-n], y[n:-n]
    return t, y - y.mean()


def amp_phase(t, y, f):
    X = np.column_stack([np.cos(2 * np.pi * f * t), np.sin(2 * np.pi * f * t)])
    c = np.linalg.lstsq(X, y, rcond=None)[0]; return np.hypot(*c), np.arctan2(-c[1], c[0])


def line_test(t, y, f, half_w, excl, nb=120):
    a, _ = amp_phase(t, y, f)
    fs = np.concatenate([np.linspace(f - half_w, f - excl, nb // 2), np.linspace(f + excl, f + half_w, nb // 2)])
    fs = fs[fs > 0.02]
    bg = np.array([amp_phase(t, y, g)[0] for g in fs])
    mid = 0.5 * (t[0] + t[-1]); h1, h2 = t < mid, t >= mid
    a1, p1 = amp_phase(t[h1], y[h1], f); a2, p2 = amp_phase(t[h2], y[h2], f)
    return a, np.mean(bg < a) * 100, np.degrees((p2 - p1 + np.pi) % (2 * np.pi) - np.pi), a2 / a1


def peaks(t, y, lo, hi, df, n=10):
    fs = np.arange(lo, hi, df); A = np.array([amp_phase(t, y, f)[0] for f in fs])
    loc = [i for i in range(1, len(A) - 1) if A[i] >= A[i - 1] and A[i] >= A[i + 1]]
    return fs, A, sorted(loc, key=lambda i: -A[i])[:n]


RES = {"A": {}, "B": {}}
print("A. Folded lines in the monthly KIT winds (1953-2026); pct = percentile among neighbours (+-0.12, excl. +-0.02); "
      "dphase = 2nd-half minus 1st-half phase; ratio = 2nd/1st amplitude")
KIT = {30: (G / "legs/qbo30_point0E/kit_qbo_u30.dat", 1), 50: (G / "legs/qbo50/qbo50.dat", 100), 70: (G / "legs/qbo70_point0E/kit_qbo_u70.dat", 1)}
specA = {}
for L, (f, sc) in KIT.items():
    a = np.loadtxt(f); t, y = prep(a[:, 0], a[:, 1] * sc); m = t >= 1953; t, y = t[m], y[m]
    print(f"\n  {L} hPa")
    for k, fr in sorted(CAND_A.items(), key=lambda kv: kv[1]):
        amp, pct, dph, rat = line_test(t, y, fr, 0.12, 0.02)
        RES["A"][f"{L}|{k}"] = dict(f=fr, amp=amp, pct=pct, dphase=dph, ratio=rat)
        flag = " <--" if pct >= 95 and abs(dph) < 45 else ""
        print(f"   {fr:7.4f} ({1 / fr:5.2f} yr) {k:42s} amp {amp:5.2f} m/s  pct {pct:5.1f}  dphase {dph:+5.0f}  ratio {rat:4.2f}{flag}")
    fs, A, top = peaks(t, y, 0.05, 2.6, 0.002)
    specA[L] = (fs, A)
    allc = {**CAND_A}
    print("   strongest peaks: " + "; ".join(f"{fs[i]:.4f} ({A[i]:.1f}) ~ {min(allc, key=lambda k: abs(allc[k] - fs[i]))} [{fs[i] - allc[min(allc, key=lambda k: abs(allc[k] - fs[i]))]:+.4f}]" for i in top[:8]))

print("\nB. Direct lunar lines in daily ERA5 (no folding); neighbours +-1.5, excl. +-0.06 cycles/yr")
DAY = {"30 hPa (0E)": G / "QBO/era_qbo_u30_point0E_daily.dat", "50 hPa": G / "QBO/era_qbo_u50_daily.dat", "70 hPa": G / "QBO/era_qbo_u70_daily.dat"}
specB = {}
for lab, f in DAY.items():
    a = np.loadtxt(f); t, y = prep(a[:, 0], a[:, 1], hp=True)
    print(f"\n  {lab} ({t[0]:.0f}-{t[-1]:.0f})")
    for k, fr in sorted(CAND_B.items(), key=lambda kv: kv[1]):
        amp, pct, dph, rat = line_test(t, y, fr, 1.5, 0.06)
        RES["B"][f"{lab}|{k}"] = dict(f=fr, amp=amp, pct=pct, dphase=dph, ratio=rat)
        flag = " <--" if pct >= 95 and abs(dph) < 45 else ""
        print(f"   {fr:8.4f} ({Y / fr:6.2f} d) {k:24s} amp {amp:5.3f} m/s  pct {pct:5.1f}  dphase {dph:+5.0f}  ratio {rat:4.2f}{flag}")
    fs, A, top = peaks(t, y, 2.5, 45.0, 0.01, 12)
    specB[lab] = (fs, A)
    print("   strongest peaks: " + "; ".join(f"{fs[i]:.2f} ({Y / fs[i]:.1f} d, {A[i]:.2f})" for i in top[:10]))
json.dump(RES, open(HERE / "qbo_tidal_lines.json", "w"), indent=1)

fig, axs = plt.subplots(2, 1, figsize=(16, 9), facecolor=SURF)
ax = axs[0]
for (L, (fs, A)), c in zip(specA.items(), ("#2a78d6", "#1baf7a", "#eb6834")):
    ax.plot(fs, A, color=c, lw=1.1, label=f"{L} hPa (KIT, monthly)")
for k, fr in CAND_A.items():
    ax.axvline(fr, color=MUTED, lw=0.6, ls=":")
    ax.text(fr, ax.get_ylim()[1] if False else 0.98, k.split(":")[0].split(" (")[0][:16], rotation=90, va="top", ha="right", fontsize=7, color=MUTED, transform=ax.get_xaxis_transform())
ax.set_yscale("log"); ax.set_xlim(0.05, 2.6); ax.set_xlabel("cycles per year", color=MUTED); ax.set_ylabel("amplitude (m/s)", color=MUTED)
ax.set_title("A. Folded tidal lines in the QBO winds (dotted: candidate positions)", loc="left", fontsize=11, color=INK); ax.legend(frameon=False, fontsize=8.5, loc="lower left")
ax = axs[1]
for (lab, (fs, A)), c in zip(specB.items(), ("#2a78d6", "#1baf7a", "#eb6834")):
    ax.plot(fs, A, color=c, lw=0.8, alpha=0.85, label=f"{lab} (ERA5 daily)")
for k, fr in CAND_B.items():
    if fr < 45:
        ax.axvline(fr, color=MUTED, lw=0.6, ls=":"); ax.text(fr, 0.98, k[:14], rotation=90, va="top", ha="right", fontsize=7, color=MUTED, transform=ax.get_xaxis_transform())
ax.set_yscale("log"); ax.set_xlim(2.5, 45); ax.set_xlabel("cycles per year", color=MUTED); ax.set_ylabel("amplitude (m/s)", color=MUTED)
ax.set_title("B. Direct lunar lines in the daily winds (no folding)", loc="left", fontsize=11, color=INK); ax.legend(frameon=False, fontsize=8.5, loc="lower left")
for a_ in axs:
    a_.set_facecolor(SURF); a_.tick_params(colors=MUTED); a_.grid(True, color=GRID, lw=0.5)
    for sp in ("top", "right"):
        a_.spines[sp].set_visible(False)
fig.tight_layout(); fig.savefig(HERE / "qbo_tidal_lines.png", dpi=95, facecolor=SURF)
print("\nplot: qbo_tidal_lines.png")
