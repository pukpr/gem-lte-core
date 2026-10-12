#!/usr/bin/env python3
"""Two legs on two indices with 2000-2005 held out: train and held-out correlations per leg, dLOD-free summary of the
shared windings, the per-leg amplitudes by harmonic, and a plot. usage: leg_scores2.py [dir] [name1] [name2]"""
import json, sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
d = sys.argv[1] if len(sys.argv) > 1 else "."; names = sys.argv[2:4] if len(sys.argv) > 3 else ["leg 1", "leg 2"]; cc = lambda a, b: float(np.corrcoef(a, b)[0, 1]); z = lambda v: (v - v.mean()) / v.std()
fig, ax = plt.subplots(2, 1, figsize=(17, 8), facecolor="#fcfcfb"); F = []
for i, name in enumerate(names):
    a = np.loadtxt(f"{d}/lte_results.{i + 1}.csv", delimiter=","); a = a[np.abs(a[:, 2]) > 0] if (a[:60, 2] == 0).all() else a; t = a[:, 0]; g = (t >= 2000) & (t < 2005); F.append(a[:, [0, 3]])
    print(f"{name}: {t[0]:.1f}-{t[-1]:.1f} | train {cc(a[~g, 1], a[~g, 2]):+.3f} | held out 2000-05 {cc(a[g, 1], a[g, 2]):+.3f} | by 36-yr era: " + " ".join(f"{cc(a[(t >= y) & (t < y + 36), 1], a[(t >= y) & (t < y + 36), 2]):+.2f}" for y in range(1880, 2020, 36) if ((t >= y) & (t < y + 36)).sum() > 100))
    A = ax[i]; A.set_facecolor("#fcfcfb"); A.axvspan(2000, 2005, color="#e4e3df", alpha=0.7, lw=0); A.plot(t, z(a[:, 2]), color="#2a78d6", lw=1.1, label=f"{name} (standardized)"); A.plot(t, z(a[:, 1]), color="#eb6834", lw=1.1, label="model (standardized)")
    A.set_title(f"{name}: train {cc(a[~g, 1], a[~g, 2]):+.3f}, held out 2000-05 {cc(a[g, 1], a[g, 2]):+.3f}", loc="left", fontsize=11); A.legend(frameon=False, fontsize=9, ncol=2); A.set_xlim(1880, 2023)
    for s in ("top", "right"): A.spines[s].set_visible(False)
n = min(len(F[0]), len(F[1])); t0 = max(F[0][0, 0], F[1][0, 0]); m0 = np.interp(F[1][F[1][:, 0] >= t0, 0], F[0][:, 0], F[0][:, 1]); print(f"manifold of the two legs: CC {cc(m0, F[1][F[1][:, 0] >= t0, 1]):.5f}")
W = [json.load(open(f"{d}/lt.exe.windings.{i}.json")) for i in (1, 2)]; K = [k for k, _, _ in W[0]["k_amp_phase"]]; base = K[1]
print(f"shared windings: slow {K[0]:+.4f}, base {base:.4f}; amplitude (phase in degrees) by harmonic:")
for j in np.argsort(np.abs(K)): print(f"   {('slow' if j == 0 else 'h' + str(round(K[j] / base))):>5s}  {names[0]} {W[0]['k_amp_phase'][j][1]:6.3f} ({np.degrees(W[0]['k_amp_phase'][j][2]):+5.0f})   {names[1]} {W[1]['k_amp_phase'][j][1]:6.3f} ({np.degrees(W[1]['k_amp_phase'][j][2]):+5.0f})")
fig.savefig(f"{d}/legs.png", dpi=100, bbox_inches="tight", facecolor="#fcfcfb"); print(f"plot: {d}/legs.png")
