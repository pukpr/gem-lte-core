#!/usr/bin/env python3
"""One fundamental cycle (1/k of the slow winding) of the winding superposition
G(F) = sum amp*sin(2*pi*k*F + phase), from each run's lt.exe.windings*.json.
Top: the four legs of the 4-leg fit. Bottom: NINO4 as legs are added (1- to 4-leg).
The linear k0*F term is left out; shading/bars mark the manifold range each fit actually visits."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
L4 = "../nino4_amo_pdo_iode_seed_cofit"
top = [(n, f"{L4}/lt.exe.windings.{i}.json") for i, n in enumerate(["NINO4", "AMO", "PDO", "IODE"], 1)]
prog = [("1-leg NINO4", "../nino4_slow_then_harm", "lt.exe.windings.json", "lte_results.csv"),
        ("2-leg + AMO", "../nino4_amo_seed_nino4", "lt.exe.windings.1.json", "lte_results.1.csv"),
        ("3-leg + PDO", "../nino4_amo_pdo_seed_cofit", "lt.exe.windings.1.json", "lte_results.1.csv"),
        ("4-leg + IODE", L4, "lt.exe.windings.1.json", "lte_results.1.csv")]


def G(w, F):
    return sum(a * np.sin(2 * np.pi * k * F + p) for k, a, p in w["k_amp_phase"])


def label(w):
    k0 = w["k_amp_phase"][0][0]
    return ", ".join("slow" if j == 0 else f"{k / k0:.0f}x" for j, (k, _, _) in enumerate(w["k_amp_phase"]))


rng = lambda d, f: (lambda m: (m.min(), m.max()))(np.loadtxt(f"{d}/{f}", delimiter=",")[:, 3])
k_slow = json.load(open(top[0][1]))["k_amp_phase"][0][0]
lo4, hi4 = rng(L4, "lte_results.1.csv")
F0 = np.floor(lo4) - 3; F = np.linspace(F0, F0 + 1 / k_slow, 1500)

fig, (a1, a2) = plt.subplots(2, 1, figsize=(12, 9.5), facecolor="#fcfcfb", sharex=True)
a1.axvspan(lo4, hi4, color="#e4e3df", alpha=0.6, lw=0)
for (n, f), c in zip(top, C):
    w = json.load(open(f)); a1.plot(F, G(w, F), color=c, lw=1.8, label=f"{n}  ({label(w)})")
a1.set_title(f"4-leg fit: winding superposition over one slow-winding cycle (1/k = {1 / k_slow:.1f} manifold units); shaded = manifold range 1880-2024",
             loc="left", fontsize=11, color="#2b2a27")
for i, ((n, d, wf, rf), c) in enumerate(zip(prog, C)):
    w = json.load(open(f"{d}/{wf}")); a2.plot(F, G(w, F), color=c, lw=1.8, label=f"{n}  ({label(w)}; slow 1/k = {1 / w['k_amp_phase'][0][0]:.1f})")
    lo, hi = rng(d, rf); a2.plot([lo, hi], [-2.35 - 0.12 * i] * 2, color=c, lw=3, solid_capstyle="butt")
a2.set_title("NINO4 leg as legs are added (bars at bottom: manifold range of each fit)", loc="left", fontsize=11, color="#2b2a27")
a2.set_xlabel("manifold F (lte_results column 4)")
for ax in (a1, a2):
    ax.set_facecolor("#fcfcfb"); ax.grid(axis="y", color="#ecebe7", lw=0.6); ax.axhline(0, color="#d5d4cf", lw=0.8, zorder=0)
    ax.set_ylabel("index units"); ax.legend(frameon=False, fontsize=9, loc="upper right" if ax is a1 else "upper left")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
a2.set_xlim(F[0], F[-1])
fig.tight_layout(); fig.savefig("winding_cycle.png", dpi=105, facecolor="#fcfcfb"); print("plot: winding_cycle.png")
