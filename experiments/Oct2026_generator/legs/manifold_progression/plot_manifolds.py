#!/usr/bin/env python3
"""Shared forcing manifold (lte_results column 4) as legs are added to the NINO4 benchmark fit:
NINO4 alone -> +AMO -> +PDO -> +IODE (each leg-1 results file; the manifold is common to all legs)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

runs = [("NINO4", "../nino4_slow_then_harm/lte_results.csv", "#2a78d6"),
        ("+ AMO", "../nino4_amo_seed_nino4/lte_results.1.csv", "#eb6834"),
        ("+ AMO + PDO", "../nino4_amo_pdo_seed_cofit/lte_results.1.csv", "#1baf7a"),
        ("+ AMO + PDO + IODE", "../nino4_amo_pdo_iode_seed_cofit/lte_results.1.csv", "#eda100")]
M = [(n, np.loadtxt(f, delimiter=","), c) for n, f, c in runs]
t, base = M[0][1][:, 0], M[0][1][:, 3]

fig, (a1, a2) = plt.subplots(2, 1, figsize=(14, 8.5), sharex=True, facecolor="#fcfcfb", gridspec_kw={"height_ratios": [3, 2]})
for n, a, c in M:
    a1.plot(a[:, 0], a[:, 3], color=c, lw=1.0, label=n)
for n, a, c in M[1:]:
    d = a[:, 3] - base
    a2.plot(a[:, 0], d, color=c, lw=1.0, label=f"{n}  (r with NINO4 {np.corrcoef(a[:, 3], base)[0, 1]:.3f}, rms change {np.sqrt(np.mean(d ** 2)):.1f})")
for ax in (a1, a2):
    ax.set_facecolor("#fcfcfb"); ax.grid(axis="y", color="#ecebe7", lw=0.6); ax.axhline(0, color="#d5d4cf", lw=0.8, zorder=0)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.legend(frameon=False, fontsize=9, loc="upper left", ncol=2 if ax is a1 else 1)
a1.set_title("Shared forcing manifold (lte_results column 4) as legs join the NINO4 fit", loc="left", fontsize=12, color="#2b2a27")
a2.set_title("Change from NINO4 alone", loc="left", fontsize=11, color="#2b2a27")
a2.set_xlim(t[0], t[-1]); a2.set_xlabel("year")
fig.tight_layout(); fig.savefig("manifold_progression.png", dpi=105, facecolor="#fcfcfb"); print("plot: manifold_progression.png")
