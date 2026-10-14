#!/usr/bin/env python3
""""Harmonic spectrum of the zone-folded modulation" for every leg of a multi-leg run, computed exactly as in
Feb2026/plot.py (BZONE panel): the IR-matched data and the IR-decompensated, scaled model are each regressed on
1, F, t, t^2 and sin/cos(2 pi h k F) for h = 1..12 (k = base winding, F = manifold), and the amplitude of each
harmonic is plotted. usage: harmonic_spectrum.py DIR NAME1 NAME2 ..."""
import json, sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
d, names = sys.argv[1], sys.argv[2:]; p = json.load(open(f"{d}/lt.exe.p")); k = float(p["ltep"][-1]); NH = 12
INK, MUTED, SURF, GRID, BL, OR = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df", "#2a78d6", "#eb6834"
def decomp(v, ir):
    u = v.copy()
    if ir != 0 and abs(ir) < 1:
        for i in range(12, len(v)): u[i] = v[i] + ir * u[i - 12]
    return u
fig, axs = plt.subplots(2, (len(names) + 1) // 2, figsize=(8.5 * ((len(names) + 1) // 2), 9), facecolor=SURF, squeeze=False); print(f"base winding k = {k:.5f}; amplitudes by harmonic (data | model), each leg in its own units")
for i, name in enumerate(names):
    a = np.loadtxt(f"{d}/lte_results.{i + 1}.csv", delimiter=","); t, model, data, F = a[:, 0], a[:, 1], a[:, 2], a[:, 3]; ir = float(json.load(open(f"{d}/lt.exe.windings.{i + 1}.json"))["IR"])
    md = decomp(model, ir); dm = data.copy(); dm[12:] = dm[12:] + ir * data[:-12]; md = md * (np.dot(md, dm) / np.dot(md, md))
    tt = t - t[0]; X = np.column_stack([np.ones_like(t), F, tt, tt ** 2] + [g(2 * np.pi * h * k * F) for h in range(1, NH + 1) for g in (np.sin, np.cos)]); cd = np.linalg.lstsq(X, dm, rcond=None)[0]; cm = np.linalg.lstsq(X, md, rcond=None)[0]
    Ad = np.hypot(cd[4::2], cd[5::2]); Am = np.hypot(cm[4::2], cm[5::2]); A = axs[i // axs.shape[1], i % axs.shape[1]]; x = np.arange(1, NH + 1); A.set_facecolor(SURF)
    A.bar(x - 0.19, Ad, 0.36, color=BL, label="Data (IR-matched)"); A.bar(x + 0.19, Am, 0.36, color=OR, label="Model (IR-decompensated, scaled)"); A.set_xticks(x); A.set_xticklabels([f"h={h}\n(1/{h * k:.2g})" for h in x], fontsize=8, color=MUTED)
    A.set_title(f"{name}: harmonic spectrum of the zone-folded modulation   (dominant: data h={int(np.argmax(Ad)) + 1}, model h={int(np.argmax(Am)) + 1})", loc="left", fontsize=11, color=INK); A.set_xlabel("Harmonic of base winding k", color=MUTED); A.set_ylabel("Amplitude", color=MUTED)
    A.grid(True, axis="y", color=GRID, lw=0.6); A.set_axisbelow(True); A.tick_params(colors=MUTED); A.set_ylim(0, 1.2 * max(Ad.max(), Am.max())); A.legend(frameon=False, fontsize=9, loc="upper center", ncol=2)
    for s in ("top", "right"): A.spines[s].set_visible(False)
    for s in ("left", "bottom"): A.spines[s].set_color(GRID)
    print(f"{name:6s} IR {ir:+.3f} | " + "  ".join(f"h{h}: {Ad[h - 1]:.3f}|{Am[h - 1]:.3f}" for h in x) + f" | CC of the two spectra {np.corrcoef(Ad, Am)[0, 1]:+.2f}")
for j in range(len(names), axs.size): axs.flat[j].axis("off")
import os; fig.suptitle(f"{os.path.basename(os.path.abspath(d))}: one shared manifold, base winding k = {k:.4f}", x=0.02, ha="left", fontsize=13, color=INK); fig.tight_layout(rect=(0, 0, 1, 0.97)); out = f"{d}/harmonic_spectrum.png"; fig.savefig(out, dpi=100, facecolor=SURF); print("plot:", out)
