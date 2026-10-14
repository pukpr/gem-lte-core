#!/usr/bin/env python3
"""Two legs on two indices with 2000-2005 held out: train and held-out correlations per leg, dLOD-free summary of the
shared windings, the per-leg amplitudes by harmonic, and a plot. usage: leg_scores2.py [dir] [name1] [name2]"""
import json, sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
H0, H1 = 2000.0, 2005.0
X0 = 1880.0; OUTPNG = "legs.png"
if "--from" in sys.argv: i = sys.argv.index("--from"); X0 = float(sys.argv[i + 1]); del sys.argv[i:i + 2]
if "--out" in sys.argv: i = sys.argv.index("--out"); OUTPNG = sys.argv[i + 1]; del sys.argv[i:i + 2]
if "--held" in sys.argv: i = sys.argv.index("--held"); H0, H1 = float(sys.argv[i + 1]), float(sys.argv[i + 2]); del sys.argv[i:i + 3]
d = sys.argv[1] if len(sys.argv) > 1 else "."; names = sys.argv[2:] if len(sys.argv) > 3 else ["leg 1", "leg 2"]; NL = len(names); cc = lambda a, b: float(np.corrcoef(a, b)[0, 1]); z = lambda v: (v - v.mean()) / v.std()
fig, ax = plt.subplots(NL, 1, figsize=(17, 4 * NL), facecolor="#fcfcfb"); F = []
for i, name in enumerate(names):
    a = np.loadtxt(f"{d}/lte_results.{i + 1}.csv", delimiter=","); a = a[np.abs(a[:, 2]) > 0] if (a[:60, 2] == 0).all() else a; t = a[:, 0]; g = (t >= H0) & (t < H1); g2 = (t >= 2000) & (t < 2005); F.append(a[:, [0, 3]])
    dq = lambda v: v - np.polyval(np.polyfit(t[~g], v[~g], 2), t); trend_only = cc(np.polyval(np.polyfit(t[~g], a[~g, 2], 2), t)[~g], a[~g, 2])
    print(f"{name}: {t[0]:.1f}-{t[-1]:.1f} | train {cc(a[~g, 1], a[~g, 2]):+.3f} (trend alone {trend_only:+.3f}; with quadratic trend removed from both {cc(dq(a[:, 1])[~g], dq(a[:, 2])[~g]):+.3f}) | held out {H0:g}-{H1:g} {cc(a[g, 1], a[g, 2]):+.3f} | 2000-05 {cc(a[g2, 1], a[g2, 2]):+.3f} | by 36-yr era: " + " ".join(f"{cc(a[(t >= y) & (t < y + 36), 1], a[(t >= y) & (t < y + 36), 2]):+.2f}" for y in range(1880, 2020, 36) if ((t >= y) & (t < y + 36)).sum() > 100))
    A = ax[i]; A.set_facecolor("#fcfcfb"); A.axvspan(H0, H1, color="#e4e3df", alpha=0.7, lw=0); A.plot(t, z(a[:, 2]), color="#2a78d6", lw=1.1, label=f"{name} (standardized)"); A.plot(t, z(a[:, 1]), color="#eb6834", lw=1.1, label="model (standardized)")
    A.set_title(f"{name}: train {cc(a[~g, 1], a[~g, 2]):+.3f}, held out {H0:g}-{H1:g} {cc(a[g, 1], a[g, 2]):+.3f}", loc="left", fontsize=11); A.legend(frameon=False, fontsize=9, ncol=2); A.set_xlim(X0, max(2023, float(t[-1])))
    for s in ("top", "right"): A.spines[s].set_visible(False)
print("manifold of each leg against leg 1: " + ", ".join(f"{cc(np.interp(F[i][:, 0], F[0][:, 0], F[0][:, 1])[F[i][:, 0] <= F[0][-1, 0]], F[i][F[i][:, 0] <= F[0][-1, 0], 1]):.5f}" for i in range(1, NL)))
W = [json.load(open(f"{d}/lt.exe.windings.{i}.json")) for i in range(1, NL + 1)]; K = [k for k, _, _ in W[0]["k_amp_phase"]]; base = K[1]
print(f"shared windings: slow {K[0]:+.4f}, base {base:.4f}; amplitude (phase in degrees) by harmonic:")
for j in np.argsort(np.abs(K)): print(f"   {('slow' if j == 0 else 'h' + str(round(K[j] / base))):>5s}  " + "   ".join(f"{names[i]} {W[i]['k_amp_phase'][j][1]:6.3f} ({np.degrees(W[i]['k_amp_phase'][j][2]):+5.0f})" for i in range(NL)))
fig.savefig(f"{d}/{OUTPNG}", dpi=100, bbox_inches="tight", facecolor="#fcfcfb"); print(f"plot: {d}/{OUTPNG}")
