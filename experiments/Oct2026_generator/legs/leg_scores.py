#!/usr/bin/env python3
"""Correlations of the two saved leg fits (lte_results.1.csv, lte_results.2.csv) inside and outside 1992-2005, and a plot.
usage: leg_scores.py [run directory]"""
import sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
d = sys.argv[1] if len(sys.argv) > 1 else "."; cc = lambda a, b: float(np.corrcoef(a, b)[0, 1]); z = lambda v: (v - v.mean()) / v.std()
fig, ax = plt.subplots(2, 1, figsize=(17, 8), facecolor="#fcfcfb")
for i, (tag, name) in enumerate(((".1", "leg 1"), (".2", "leg 2"))):
    a = np.loadtxt(f"{d}/lte_results{tag}.csv", delimiter=","); t = a[:, 0]; w = (t >= 1992) & (t < 2005); af = t >= 2005; k = max(3, int(round(0.5 / np.median(np.diff(t)))) | 1); hi = lambda v: v - np.convolve(v, np.ones(k) / k, mode="same")
    print(f"{name}: {len(t)} samples {t[0]:.2f}-{t[-1]:.2f} | CC inside 1992-2005 {cc(a[w, 1], a[w, 2]):+.3f} | outside {cc(a[~w, 1], a[~w, 2]):+.3f} | after 2005 {cc(a[af, 1], a[af, 2]):+.3f} | faster than 6 months: inside {cc(hi(a[:, 1])[w], hi(a[:, 2])[w]):+.3f}, outside {cc(hi(a[:, 1])[~w], hi(a[:, 2])[~w]):+.3f}")
    A = ax[i]; A.set_facecolor("#fcfcfb"); A.axvspan(1992, 2005, color="#e4e3df", alpha=0.6, lw=0); A.plot(t, z(a[:, 2]), color="#2a78d6", lw=1.0, label="data (standardized)"); A.plot(t, z(a[:, 1]), color="#eb6834", lw=1.0, label="model (standardized)")
    A.set_title(f"{name}: inside 1992-2005 {cc(a[w, 1], a[w, 2]):+.3f}, outside {cc(a[~w, 1], a[~w, 2]):+.3f}", loc="left", fontsize=11); A.legend(frameon=False, fontsize=9, ncol=2); A.set_xlim(t[0], t[-1])
    for s in ("top", "right"): A.spines[s].set_visible(False)
fig.savefig(f"{d}/legs.png", dpi=100, bbox_inches="tight", facecolor="#fcfcfb"); print(f"plot: {d}/legs.png")
