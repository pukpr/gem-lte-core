#!/usr/bin/env python3
"""150 hPa central-Pacific wind fitted with and without the 5-day SOI as a second leg (2026-10-22).

Runs (experiments/Oct2026_generator/legs): wind150_1leg (wind alone) and wind150_soi5d (wind + 5-day SOI), both seeded
from the unified SOI fit soi_universal/soi_nino4_h20 (base 0.2016, harmonics to 20, delB on the monthly slot,
IMPULSE_EXACT), 2011-2016 held out. Compares per-leg scores, the shared windings each search kept, and each leg's
fitted amplitude per winding (lt.exe.windings*.json, k_amp_phase), including the share carried by windings >= 10x.
    python3 wind_soi_windings.py  -> wind_soi_windings.png, stdout
"""
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
L = HERE.parent / "experiments/Oct2026_generator/legs"
INK, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
H0, H1 = 2011.0, 2016.0
RUNS = {"seed (nothing searched)": (L / "wind150_soi5d", "seed.p", None),
        "wind alone": (L / "wind150_1leg", "lt.exe.p", {"WIND150": ("lte_results.csv", "lt.exe.windings.json")}),
        "wind + 5-day SOI": (L / "wind150_soi5d", "lt.exe.p", {"WIND150": ("lte_results.1.csv", "lt.exe.windings.1.json"),
                                                              "SOI5D": ("lte_results.2.csv", "lt.exe.windings.2.json")})}


def score(f):
    a = np.loadtxt(f, delimiter=","); t, m, y = a[:, 0], a[:, 1], a[:, 2]; ok = y != 0; h = (t >= H0) & (t < H1)
    return np.corrcoef(m[ok & ~h], y[ok & ~h])[0, 1], np.corrcoef(m[ok & h], y[ok & h])[0, 1]


def amps(f, base):
    w = json.load(open(f)); out = {}
    for k, a, ph in w["k_amp_phase"]:
        out[round(abs(k) / base, 2)] = abs(a)
    return out


AMP = {}
for name, (d, pfile, legs) in RUNS.items():
    p = json.load(open(d / pfile)); base = abs(p["ltep"][-1])
    print(f"\n{name}: base {base:.5f}, slow {abs(p['ltep'][0]) / base:.3f}x, harmonics {sorted(int(h) for h in p['harm'])}")
    if legs is None:
        continue
    for lab, (res, wj) in legs.items():
        tr, ho = score(d / res); A = amps(d / wj, base); AMP[(name, lab)] = A
        tot = sum(v ** 2 for v in A.values()); hi = sum(v ** 2 for h, v in A.items() if h >= 9.5)
        print(f"   {lab:8s} train {tr:+.3f}  held out {H0:.0f}-{H1:.0f} {ho:+.3f}  | share of fitted amplitude^2 in windings >= 10x: {hi / tot:.2f}")
        print("            amplitude by winding: " + "  ".join(f"{h:g}x {v:.2f}" for h, v in sorted(A.items())))
if ("wind alone", "WIND150") in AMP and ("wind + 5-day SOI", "WIND150") in AMP:
    a, b = AMP[("wind alone", "WIND150")], AMP[("wind + 5-day SOI", "WIND150")]
    ks = sorted(set(a) & set(b))
    if len(ks) > 2:
        print(f"\nwind-leg amplitudes, wind alone vs with SOI, on {len(ks)} shared windings: r = {np.corrcoef([a[k] for k in ks], [b[k] for k in ks])[0, 1]:+.2f}")
    s, w = AMP.get(("wind + 5-day SOI", "SOI5D")), b
    ks = sorted(set(s) & set(w))
    print(f"wind vs SOI leg amplitudes in the joint run, {len(ks)} windings: r = {np.corrcoef([w[k] for k in ks], [s[k] for k in ks])[0, 1]:+.2f}")

fig, ax = plt.subplots(figsize=(14, 4.8), facecolor=SURF); wdt = 0.27
allk = sorted(set().union(*[set(v) for v in AMP.values()]))
xx = np.arange(len(allk))
for j, ((name, lab), c) in enumerate(zip(AMP, ("#eb6834", "#2a78d6", "#1baf7a"))):
    A = AMP[(name, lab)]; tot = np.sqrt(sum(v ** 2 for v in A.values()))
    ax.bar(xx + (j - 1) * wdt, [A.get(k, 0) / tot for k in allk], wdt, color=c, label=f"{lab}, {name}")
ax.set_xticks(xx); ax.set_xticklabels([f"{k:g}x" for k in allk], fontsize=8.5); ax.set_ylabel("amplitude (share of leg total)", color=MUTED)
ax.set_title("Fitted amplitude per winding: 150 hPa wind alone, wind with the 5-day SOI leg, and the SOI leg", loc="left", fontsize=11, color=INK)
ax.legend(frameon=False, fontsize=9)
ax.set_facecolor(SURF); ax.tick_params(colors=MUTED); ax.grid(True, color=GRID, lw=0.5, axis="y")
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
fig.tight_layout(); fig.savefig(HERE / "wind_soi_windings.png", dpi=95, facecolor=SURF)
print("plot: wind_soi_windings.png")
