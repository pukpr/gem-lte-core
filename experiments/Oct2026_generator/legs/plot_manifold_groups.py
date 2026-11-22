#!/usr/bin/env python3
"""Plot two N-leg groups' forcing manifolds F(t) against each other.

    python3 plot_manifold_groups.py [RUN_A] [RUN_B] [--extra RUN] [-o manifold_groups.png]

F(t) is column 4 of lte_results.1.csv (shared by every leg of a run). Panels:
full record, a 10-year zoom of the monthly structure, F_A vs F_B scatter, and the
power spectrum of the monthly increments dF with the annual-impulse alias periods
that separate the groups marked (Mf 3.8 yr, the 2.0-2.15 yr cluster).
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
INK, INK2, GRID, SURF = "#2b2b29", "#6b6a63", "#e4e3de", "#fcfcfb"
C = ["#2a78d6", "#eb6834", "#1baf7a"]   # palette slots 1-3, validated (slot 3 dashed + labelled)


def manifold(run):
    a = np.loadtxt(run / "lte_results.1.csv", delimiter=",")
    return a[:, 0], a[:, 3]


def dlod(run):
    c = np.loadtxt(run / "dlod_compare.csv", delimiter=",")
    return np.corrcoef(c[:, 1], c[:, 2])[0, 1]


def label(run):
    legs = [l.split()[1].replace(".dat", "").upper()
            for f in sorted(run.glob("lt.exe.resp.[0-9]")) for l in f.read_text().splitlines()
            if l.startswith("CLIMATE_INDEX")]
    return "+".join(legs) if legs else run.name


def spectrum(f):
    d = np.diff(f)
    d = (d - d.mean()) * np.hanning(len(d))
    p = np.abs(np.fft.rfft(d)) ** 2
    fr = np.fft.rfftfreq(len(d), 1 / 12)
    return fr[1:], p[1:] / p[1:].sum()


def style(ax):
    ax.set_facecolor(SURF)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("a", nargs="?", default="amo_pdo_iode_nino4_feb")
    ap.add_argument("b", nargs="?", default="nao_baltic_lag0_w3")
    ap.add_argument("--extra", default="amo_pdo_iode_nino4_feb_overnight_cp")
    ap.add_argument("--zoom", nargs=2, type=float, default=(1995, 2005))
    ap.add_argument("-o", default="manifold_groups.png")
    a = ap.parse_args()
    runs = [HERE / a.a, HERE / a.b] + ([HERE / a.extra] if a.extra else [])
    data = [manifold(r) for r in runs]
    names = [f"{label(r)}  ({r.name}, dLOD {dlod(r):.3f})" for r in runs]
    short = [label(r) for r in runs]
    short[2:] = [f"{label(r)} best-scoring fit" for r in runs[2:]]

    fig = plt.figure(figsize=(13, 9.5), facecolor=SURF)
    gs = fig.add_gridspec(3, 2, height_ratios=[1.15, 1, 1.05], hspace=0.42, wspace=0.18)
    ax1 = fig.add_subplot(gs[0, :]); ax2 = fig.add_subplot(gs[1, :])
    ax3 = fig.add_subplot(gs[2, 0]); ax4 = fig.add_subplot(gs[2, 1])
    for ax in (ax1, ax2, ax3, ax4):
        style(ax)

    # 1. full record
    order = [2, 0, 1] if len(runs) > 2 else [0, 1]
    for i in order:
        t, f = data[i]
        ax1.plot(t, f, color=C[i], lw=0.7 if i == 2 else 1.1, ls="--" if i == 2 else "-",
                 alpha=0.9, label=names[i])
    for i in (0, 1):
        t, f = data[i]
        j = np.searchsorted(t, 1986)
        ax1.annotate(short[i], (t[j], f[j]), xytext=(0, 14 if i == 0 else -18), textcoords="offset points",
                     color=INK, fontsize=8.5, fontweight="bold")
    ax1.set_title("Forcing manifold F(t), monthly, as fitted (column 4 of lte_results)", loc="left", color=INK, fontsize=10.5)
    ax1.set_ylabel("F", color=INK2, fontsize=9)
    ax1.legend(loc="upper left", fontsize=7.5, frameon=False, labelcolor=INK)

    # 2. zoom
    z0, z1 = a.zoom
    for i in order:
        t, f = data[i]
        s = (t >= z0) & (t <= z1)
        ax2.plot(t[s], f[s], color=C[i], lw=0.8 if i == 2 else 1.3, ls="--" if i == 2 else "-",
                 marker=None if i == 2 else "o", ms=2.2, alpha=0.95)
    ax2.set_title(f"Zoom {z0:.0f}-{z1:.0f}: month-to-month structure (one impulse per year)", loc="left", color=INK, fontsize=10.5)
    ax2.set_ylabel("F", color=INK2, fontsize=9)

    # 3. scatter
    ka = {int(round(t * 12)): f for t, f in zip(*data[0])}
    kb = {int(round(t * 12)): f for t, f in zip(*data[1])}
    ks = sorted(set(ka) & set(kb))
    fa, fb = np.array([ka[k] for k in ks]), np.array([kb[k] for k in ks])
    yr = np.array(ks) / 12
    sc = ax3.scatter(fa, fb, c=yr, cmap="Blues", s=6, alpha=0.8, linewidths=0)
    cb = fig.colorbar(sc, ax=ax3, pad=0.01)
    cb.ax.tick_params(labelsize=7, colors=INK2); cb.outline.set_visible(False)
    cb.set_label("year", color=INK2, fontsize=8)
    ax3.set_xlabel(f"F, {short[0]}", color=INK2, fontsize=9)
    ax3.set_ylabel(f"F, {short[1]}", color=INK2, fontsize=9)
    ax3.set_title(f"F vs F: CC {np.corrcoef(fa, fb)[0, 1]:+.2f}, monthly increments {np.corrcoef(np.diff(fa), np.diff(fb))[0, 1]:+.2f}",
                  loc="left", color=INK, fontsize=10.5)

    # 4. spectrum of increments
    for i in (0, 1):
        fr, p = spectrum(data[i][1])
        per = 1 / fr
        s = (per >= 1.2) & (per <= 12)
        ax4.plot(per[s], p[s], color=C[i], lw=1.1, label=short[i])
    top = ax4.get_ylim()[1]
    for x, txt, y in ((2.15, " 2.0-2.15 yr alias cluster", 1.10), (3.80, " Mf alias 3.8 yr", 0.80)):
        ax4.axvline(x, color=INK2, lw=0.8, ls=":")
        ax4.text(x, top * y, txt, color=INK2, fontsize=7.5, ha="left", va="bottom")
    ax4.set_ylim(0, top * 1.2)
    ax4.axvspan(2.00, 2.15, color=GRID, alpha=0.6, lw=0)
    ax4.set_xscale("log")
    ax4.set_xticks([1.5, 2, 3, 4, 6, 8, 12]); ax4.set_xticklabels(["1.5", "2", "3", "4", "6", "8", "12"])
    ax4.set_xlabel("period of dF component (yr)", color=INK2, fontsize=9)
    ax4.set_ylabel("share of dF power", color=INK2, fontsize=9)
    ax4.set_title("Spectrum of the monthly increments dF, 1.2-12 yr", loc="left", color=INK, fontsize=10.5)
    ax4.legend(fontsize=8, frameon=False, labelcolor=INK)

    out = HERE / a.o
    fig.savefig(out, dpi=130, facecolor=SURF, bbox_inches="tight")
    print(out)


if __name__ == "__main__":
    main()
