#!/usr/bin/env python3
"""Two legs on two indices with 2000-2005 held out: train and held-out correlations per leg, dLOD-free summary of the
shared windings, the per-leg amplitudes by harmonic, and a plot. usage: leg_scores2.py [dir] [name1] [name2]"""
import json, sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
H0, H1 = 2000.0, 2005.0
X0 = None; OUTPNG = "legs.png"   # X0 None: start at the decade of the earliest leg's first sample
if "--from" in sys.argv: i = sys.argv.index("--from"); X0 = float(sys.argv[i + 1]); del sys.argv[i:i + 2]
if "--out" in sys.argv: i = sys.argv.index("--out"); OUTPNG = sys.argv[i + 1]; del sys.argv[i:i + 2]
if "--held" in sys.argv: i = sys.argv.index("--held"); H0, H1 = float(sys.argv[i + 1]), float(sys.argv[i + 2]); del sys.argv[i:i + 3]
d = sys.argv[1] if len(sys.argv) > 1 else "."; names = sys.argv[2:] if len(sys.argv) > 3 else ["leg 1", "leg 2"]; NL = len(names); cc = lambda a, b: float(np.corrcoef(a, b)[0, 1]); z = lambda v: (v - v.mean()) / v.std()
fig, ax = plt.subplots(NL + 1, 1, figsize=(17, 4 * (NL + 1)), facecolor="#fcfcfb"); F = []   # last panel: every leg's manifold together
if X0 is None:
    def _start(k):
        a = np.loadtxt(f"{d}/lte_results.{k}.csv", delimiter=","); nz = a[a[:, 2] != 0, 0]; return float(nz[0]) if len(nz) else float(a[0, 0])
    X0 = 10.0 * np.floor(min(_start(k) for k in range(1, NL + 1)) / 10.0)
for i, name in enumerate(names):
    a = np.loadtxt(f"{d}/lte_results.{i + 1}.csv", delimiter=","); a = a[np.abs(a[:, 2]) > 0] if (a[:60, 2] == 0).all() else a; t = a[:, 0]; g = (t >= H0) & (t < H1); g2 = (t >= 2000) & (t < 2005); F.append(a[:, [0, 3]])
    dq = lambda v: v - np.polyval(np.polyfit(t[~g], v[~g], 2), t); trend_only = cc(np.polyval(np.polyfit(t[~g], a[~g, 2], 2), t)[~g], a[~g, 2])
    print(f"{name}: {t[0]:.1f}-{t[-1]:.1f} | train {cc(a[~g, 1], a[~g, 2]):+.3f} (trend alone {trend_only:+.3f}; with quadratic trend removed from both {cc(dq(a[:, 1])[~g], dq(a[:, 2])[~g]):+.3f}) | held out {H0:g}-{H1:g} {cc(a[g, 1], a[g, 2]):+.3f} | 2000-05 {cc(a[g2, 1], a[g2, 2]):+.3f} | by 36-yr era: " + " ".join(f"{cc(a[(t >= y) & (t < y + 36), 1], a[(t >= y) & (t < y + 36), 2]):+.2f}" for y in range(1880, 2020, 36) if ((t >= y) & (t < y + 36)).sum() > 100))
    A = ax[i]; A.set_facecolor("#fcfcfb"); A.axvspan(H0, H1, color="#e4e3df", alpha=0.7, lw=0); A.plot(t, z(a[:, 2]), color="#2a78d6", lw=1.1, label=f"{name} (standardized)"); A.plot(t, z(a[:, 1]), color="#eb6834", lw=1.1, label="model (standardized)")
    A.set_title(f"{name}: train {cc(a[~g, 1], a[~g, 2]):+.3f}, held out {H0:g}-{H1:g} {cc(a[g, 1], a[g, 2]):+.3f}", loc="left", fontsize=11); A.legend(frameon=False, fontsize=9, ncol=2); A.set_xlim(X0, max(2023, float(t[-1])))
    for s in ("top", "right"): A.spines[s].set_visible(False)
print("manifold of each leg against leg 1: " + ", ".join(f"{cc(np.interp(F[i][:, 0], F[0][:, 0], F[0][:, 1])[F[i][:, 0] <= F[0][-1, 0]], F[i][F[i][:, 0] <= F[0][-1, 0], 1]):.5f}" for i in range(1, NL)))
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]   # categorical order, fixed per leg
A = ax[NL]; A.set_facecolor("#fcfcfb"); A.axvspan(H0, H1, color="#e4e3df", alpha=0.7, lw=0); x1 = max(2023, max(float(f[-1, 0]) for f in F))
for i in range(NL):   # column 4 of lte_results.<leg>.csv, raw (same units on every leg), so offsets and amplitudes show as they are
    m = F[i][:, 0] <= F[0][-1, 0]; r = 1.0 if i == 0 else cc(np.interp(F[i][m, 0], F[0][:, 0], F[0][:, 1]), F[i][m, 1])
    A.plot(F[i][:, 0], F[i][:, 1], color=SERIES[i % len(SERIES)], lw=1.1, alpha=0.9, label=names[i] if i == 0 else f"{names[i]} (CC vs {names[0]} {r:.4f})")
lo, hi = A.get_ylim(); A.set_ylim(lo, hi + 0.18 * (hi - lo))   # headroom for the legend
if NL <= 4:   # end labels, pushed apart where the lines finish close together, each with a leader to its line's end
    ys = [float(f[-1, 1]) for f in F]; gap = 0.06 * (hi - lo); y = {}
    for i in sorted(range(NL), key=lambda j: ys[j]): y[i] = ys[i] if not y else max(ys[i], max(y.values()) + gap)
    for i in range(NL): A.annotate(names[i], (F[i][-1, 0], ys[i]), xytext=(F[i][-1, 0] + 1.5, y[i]), textcoords="data", va="center", fontsize=9, color="#52514d", arrowprops=dict(arrowstyle="-", color="#a3a29c", lw=0.6, shrinkA=0, shrinkB=1))
A.set_title(f"manifolds of all {NL} legs (lte_results column 4, not standardized)", loc="left", fontsize=11); A.legend(frameon=False, fontsize=9, ncol=min(NL, 4), loc="upper left"); A.set_xlim(X0, x1 + (6 if NL <= 4 else 0))
for s_ in ("top", "right"): A.spines[s_].set_visible(False)
W = [json.load(open(f"{d}/lt.exe.windings.{i}.json")) for i in range(1, NL + 1)]; K = [k for k, _, _ in W[0]["k_amp_phase"]]; base = K[1] if len(K) > 1 else K[0]   # no harmonics: slow winding only
KL = [[k for k, _, _ in w["k_amp_phase"]] for w in W]
if any(len(k) != len(K) or not np.allclose(k, K, rtol=0, atol=1e-12) for k in KL[1:]):   # LEG_WINDINGS: each leg its own ltep/harm
    print("per-leg windings (LEG_WINDINGS): wavenumber amplitude (phase in degrees)")
    for i in range(NL): print(f"   {names[i]:>8s}  " + "   ".join(f"{k:+.4f} {a_:6.3f} ({np.degrees(p_):+5.0f})" for k, a_, p_ in W[i]["k_amp_phase"]))
else: print(f"shared windings: slow {K[0]:+.4f}, base {base:.4f}; amplitude (phase in degrees) by harmonic:")
for j in (np.argsort(np.abs(K)) if not any(len(k) != len(K) or not np.allclose(k, K, rtol=0, atol=1e-12) for k in KL[1:]) else []): print(f"   {('slow' if j == 0 else f'{K[j] / K[0]:.0f}x'):>5s}  " + "   ".join((f"{names[i]} {W[i]['k_amp_phase'][j][1]:6.3f} ({np.degrees(W[i]['k_amp_phase'][j][2]):+5.0f})" if j < len(W[i]['k_amp_phase']) else f"{names[i]}      -       ") for i in range(NL)))   # a leg with fewer windings (shorter per-leg NH) shows -
fig.savefig(f"{d}/{OUTPNG}", dpi=100, bbox_inches="tight", facecolor="#fcfcfb"); print(f"plot: {d}/{OUTPNG}")
