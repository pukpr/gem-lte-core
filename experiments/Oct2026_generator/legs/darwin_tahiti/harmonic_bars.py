#!/usr/bin/env python3
"""Bar chart of the fitted response: amplitude of every winding term (slow winding, base, and each harmonic of the base)
for the Darwin and Tahiti legs, from lt.exe.windings.1.json / .2.json, in hPa. The phase difference between the two
legs is written over each pair (180 degrees = exact see-saw)."""
import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
p = json.load(open("result.p")); kb = abs(p["ltep"][-1]); W = [json.load(open(f"lt.exe.windings.{i}.json"))["k_amp_phase"] for i in (1, 2)]
lab = lambda k: "slow" if abs(abs(k) - abs(p["ltep"][0])) < 1e-9 else f"h{round(abs(k) / kb)}"; order = sorted(range(len(W[0])), key=lambda j: (lab(W[0][j][0]) != "slow", abs(W[0][j][0])))
names = [lab(W[0][j][0]) for j in order]; A = [[W[i][j][1] for j in order] for i in (0, 1)]; dph = [((np.degrees(W[1][j][2] - W[0][j][2]) + 180) % 360) - 180 for j in order]
INK, MUTED, SURF, GRID, BL, OR = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df", "#2a78d6", "#eb6834"; x = np.arange(len(names))
fig, ax = plt.subplots(figsize=(16, 5.8), facecolor=SURF); ax.set_facecolor(SURF); ax.bar(x - 0.2, A[0], 0.38, color=BL, label="Darwin"); ax.bar(x + 0.2, A[1], 0.38, color=OR, label="Tahiti")
top = max(max(A[0]), max(A[1]))
for xi, a0, a1, d in zip(x, A[0], A[1], dph): ax.text(xi, max(a0, a1) + 0.02 * top, f"{d:+.0f}°", ha="center", va="bottom", fontsize=8, color=MUTED)
ax.set_xticks(x); ax.set_xticklabels([f"{n}\n(k={abs(W[0][j][0]):.3g})" for n, j in zip(names, order)], fontsize=8, color=MUTED); ax.set_ylabel("amplitude (hPa)", color=MUTED); ax.set_ylim(0, top * 1.15)
ax.set_title(f"Darwin and Tahiti: fitted response by winding term (base k = {kb:.4f}); numbers = Tahiti phase minus Darwin phase", loc="left", fontsize=11.5, color=INK)
ax.legend(frameon=False, fontsize=10, loc="upper right"); ax.grid(True, axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True); ax.tick_params(colors=MUTED)
for s in ("top", "right"): ax.spines[s].set_visible(False)
for s in ("left", "bottom"): ax.spines[s].set_color(GRID)
fig.savefig("harmonic_bars.png", dpi=100, bbox_inches="tight", facecolor=SURF)
for n, a0, a1, d in zip(names, A[0], A[1], dph): print(f"{n:>5s}: Darwin {a0:.3f}  Tahiti {a1:.3f} hPa   phase difference {d:+.0f}")
print("plot: harmonic_bars.png")
