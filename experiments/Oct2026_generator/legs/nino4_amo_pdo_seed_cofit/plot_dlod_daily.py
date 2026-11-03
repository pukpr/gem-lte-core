#!/usr/bin/env python3
"""Daily model tidal forcing (sign flipped, least-squares scaled) against dLOD over a 2-year window.
usage: python3 plot_dlod_daily.py [start end]   (dLOD in ../../../Feb2026/dlod3.dat begins 1962.0)"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

t0, t1 = (float(sys.argv[1]), float(sys.argv[2])) if len(sys.argv) > 2 else (1962.0, 1964.0)
c = np.loadtxt("dlod_compare.csv", delimiter=",")          # date, dLOD tidal fit, model tide sum
obs = np.loadtxt("../../../Feb2026/dlod3.dat")             # date, dLOD (ms)
t, model = c[:, 0], -c[:, 2]
assert np.allclose(t, obs[:, 0])
w = (t >= t0) & (t < t1)
k = np.polyfit(model, obs[:, 1], 1)                        # one scale + offset over the full record
r_all, r_win = np.corrcoef(obs[:, 1], model)[0, 1], np.corrcoef(obs[w, 1], model[w])[0, 1]
for y in range(1962, 2018, 2):
    g = (t >= y) & (t < y + 2); print(f"{y}-{y + 2}: r = {np.corrcoef(obs[g, 1], model[g])[0, 1]:+.3f}")

fig, ax = plt.subplots(figsize=(13, 4.2), facecolor="#fcfcfb"); ax.set_facecolor("#fcfcfb")
ax.plot(t[w], obs[w, 1], color="#2a78d6", lw=1.4, label="dLOD (ms)")
ax.plot(t[w], np.polyval(k, model[w]), color="#eb6834", lw=1.4, label="model tidal forcing, sign flipped and scaled")
ax.axhline(0, color="#d5d4cf", lw=0.8, zorder=0)
ax.set_xlim(t0, t1); ax.set_ylabel("ms"); ax.grid(axis="y", color="#ecebe7", lw=0.6)
ax.set_title(f"Daily dLOD vs sign-flipped model forcing, {t0:g}-{t1:g}: r = {r_win:+.2f} (full 1962-2019 record {r_all:+.2f})",
             loc="left", fontsize=11, color="#2b2a27")
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.legend(frameon=False, fontsize=9, ncol=2, loc="upper left")
out = f"dlod_daily_{t0:g}_{t1:g}.png"; fig.savefig(out, dpi=110, bbox_inches="tight", facecolor="#fcfcfb"); print("plot:", out)
