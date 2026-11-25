#!/usr/bin/env python3
"""Plot the power spectrum of Baltic (or Stockholm) MSL alongside the
annual-impulse tidal aliasing ladder (Tropical T, Anomalistic P, Draconic D,
and Chandler Wobble P14) and compare Pre-Bessel dLOD F(t) vs LTE-wound output.

Usage:
    python3 plot_baltic_aliasing.py [run_dir] [--leg 2] [--out spectrum.png] [--fmax 1.15]
"""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

YR = 365.2422
F_T = YR / 27.32166155   # Tropical month (13.3682/yr)
F_D = YR / 27.21222082   # Draconic month (13.4220/yr)
F_P = YR / 27.55454988   # Anomalistic month (13.2552/yr)

def alias_yr(f):
    return abs(f - round(f))

# Annual-impulse aliasing ladder: (freq_cyc_yr, label, family_key, y_tier)
# family_key in {"D": Draconic, "T": Tropical/Mf, "P": Anomalistic/Mm, "M": Mixed/Solar}
LADDER = [
    (F_D - F_T,                "D-T\n18.6y",      "D", 0.94),
    (F_T - F_P,                "T-P\n8.85y",      "P", 0.78),
    (alias_yr(2 * F_D),        "2D\n6.41y",       "D", 0.94),
    (F_D - F_P,                "D-P\n6.00y",      "M", 0.76),
    (alias_yr(F_T + F_D),      "T+D\n4.77y",      "M", 0.94),
    (alias_yr(F_P),            "P (Mm)\n3.92y",   "P", 0.76),
    (alias_yr(2 * F_T),        "2T (Mf)\n3.79y",  "T", 0.94),
    (2 * (F_D - F_P),          "2(D-P)\n3.00y",   "M", 0.84),
    (alias_yr(F_T),            "T\n2.72y",        "T", 0.94),
    (alias_yr(F_D),            "D\n2.37y",        "D", 0.94),
    (1.0 - alias_yr(F_D),      "1-D\n1.73y",      "D", 0.94),
    (1.0 - alias_yr(F_T),      "1-T\n1.58y",      "T", 0.78),
    (1.0 - alias_yr(2 * F_T),  "1-2T\n1.36y",     "T", 0.94),
    (1.0 - alias_yr(F_T + F_D),"1-(T+D)\n1.27y",  "M", 0.78),
    (1.0 - alias_yr(2 * F_D),  "1-2D / P14\n14.2m","D", 0.95),
    (1.0000,                   "Ann\n1.00y",      "S", 0.80),
]

COLORS = {
    "D": "#c62828",  # Draconic (crimson)
    "T": "#1565c0",  # Tropical / Mf (blue)
    "P": "#2e7d32",  # Anomalistic / Mm (green)
    "M": "#6a1b9a",  # Mixed T-D / D-P (purple)
    "S": "#555555",  # Annual (gray)
}

def power_spec(x, dt=1.0 / 12.0, hann=True):
    x0 = x - np.mean(x)
    x0 = x0 / (np.std(x0) + 1e-15)
    if hann:
        w = np.hanning(len(x0))
        x0 = x0 * w
    f = np.fft.rfftfreq(len(x0), d=dt)
    P = np.abs(np.fft.rfft(x0)) ** 2
    P /= np.sum(P)
    return f, P

def band_ratio(f, P):
    p_d = P[(f >= 0.14) & (f <= 0.19)].sum()
    p_tp = P[(f >= 0.24) & (f <= 0.28)].sum()
    return p_d, p_tp, p_d / (p_tp + 1e-15)

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("dir", nargs="?",
                    default="/home/paul/gemini/gem-lte-core/experiments/Oct2026_generator/legs/atlantic_family",
                    help="Run directory containing lte_results.<leg>.csv")
    ap.add_argument("--leg", type=int, default=2, help="Leg index for Baltic (default: 2)")
    ap.add_argument("--name", default="Baltic MSL (baltic_3pt_na)", help="Series display name")
    ap.add_argument("--out", default="", help="Output PNG path (default: <dir>/baltic_aliasing_spectrum.png)")
    ap.add_argument("--fmax", type=float, default=1.12, help="Max frequency in cyc/yr (default: 1.12)")
    ap.add_argument("--log", action="store_true", help="Use log y-scale")
    args = ap.parse_args()

    d = Path(args.dir)
    csv_path = d / f"lte_results.{args.leg}.csv"
    arr = np.loadtxt(csv_path, delimiter=",")
    m = arr[:, 2] != 0
    t, y_mod, y_dat, F_pre = arr[m, 0], arr[m, 1], arr[m, 2], arr[m, 3]
    dt = float(np.median(np.diff(t)))

    f, P_dat = power_spec(y_dat, dt=dt, hann=True)
    _, P_mod = power_spec(y_mod, dt=dt, hann=True)
    _, P_pre = power_spec(F_pre, dt=dt, hann=True)

    fm = (f >= 0.02) & (f <= args.fmax)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8.2), sharex=True,
                                   gridspec_kw={"height_ratios": [1.35, 1.0], "hspace": 0.12})

    # --- Panel 1: Baltic Data vs Baltic LTE Model ---
    ax1.fill_between(f[fm], 0, P_dat[fm] * 100, color="#90a4ae", alpha=0.35, label=f"{args.name} (data)")
    ax1.plot(f[fm], P_dat[fm] * 100, color="#263238", lw=1.25)
    ax1.plot(f[fm], P_mod[fm] * 100, color="#d84315", lw=1.35, alpha=0.90, label=f"{args.name} LTE model (Bessel + wound)")

    ymax1 = max(float(P_dat[fm].max()), float(P_mod[fm].max())) * 100
    if args.log:
        ax1.set_yscale("log")
        ax1.set_ylim(1e-3, ymax1 * 3.5)
    else:
        ax1.set_ylim(0, ymax1 * 1.28)

    for freq, lbl, fam, ytier in LADDER:
        if freq < 0.02 or freq > args.fmax:
            continue
        col = COLORS[fam]
        ls = "-" if fam == "D" else "--"
        lw = 1.15 if fam == "D" else 0.85
        ax1.axvline(freq, color=col, ls=ls, lw=lw, alpha=0.55, zorder=1)
        y_pos = (ymax1 * (1.02 if ytier > 0.88 else 0.84)) if not args.log else (ymax1 * (2.2 if ytier > 0.88 else 0.95))
        ax1.text(freq, y_pos, lbl, color=col, fontsize=7.6, fontweight="bold" if fam == "D" else "normal",
                 ha="center", va="bottom",
                 bbox=dict(boxstyle="round,pad=0.14", fc="white", ec=col, alpha=0.82, lw=0.6))

    d_dat, tp_dat, r_dat = band_ratio(f, P_dat)
    d_mod, tp_mod, r_mod = band_ratio(f, P_mod)
    d_pre, tp_pre, r_pre = band_ratio(f, P_pre)

    ax1.set_ylabel("Normalized Power (%)", fontsize=10.5)
    ax1.set_title(
        f"{args.name} Power Spectrum & Annual-Impulse Tidal Aliasing Ladder (T, P, D → P14)\n"
        f"Draconic/Tropical ratio  2D/(2T,P):  Pre-Bessel dLOD F = {r_pre:.2f}  →  LTE Model = {r_mod:.2f}  |  Data = {r_dat:.2f}",
        fontsize=11.5, fontweight="bold", pad=10,
    )
    ax1.grid(True, alpha=0.22)

    # Legend with family color keys
    handles, labels = ax1.get_legend_handles_labels()
    handles += [
        Line2D([0], [0], color=COLORS["D"], lw=1.4, ls="-",  label="Draconic (D, 2D, D-T, P14)"),
        Line2D([0], [0], color=COLORS["T"], lw=1.2, ls="--", label="Tropical / Mf (T, 2T)"),
        Line2D([0], [0], color=COLORS["P"], lw=1.2, ls="--", label="Anomalistic / Mm (P, T-P)"),
        Line2D([0], [0], color=COLORS["M"], lw=1.2, ls="--", label="T-D / D-P intermodulation"),
    ]
    ax1.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.50, 0.60),
               ncol=3, fontsize=8.2, framealpha=0.92)

    # --- Panel 2: Pre-Bessel dLOD F(t) vs LTE Model (folding transfer) ---
    # Rescale Pre-Bessel F to its in-band [0.02..fmax] share so its T/P vs D structure is clearly visible
    P_pre_band = P_pre / (P_pre[fm].sum() + 1e-15) * 100.0
    P_mod_band = P_mod / (P_mod[fm].sum() + 1e-15) * 100.0
    ax2.fill_between(f[fm], 0, P_pre_band[fm], color="#1565c0", alpha=0.28, label="Pre-Bessel dLOD F(t) (in-band %)")
    ax2.plot(f[fm], P_pre_band[fm], color="#0d47a1", lw=1.35)
    ax2.plot(f[fm], P_mod_band[fm], color="#d84315", lw=1.35, alpha=0.90, label="Post-Bessel + LTE wound model (in-band %)")

    for freq, lbl, fam, _ in LADDER:
        if freq < 0.02 or freq > args.fmax:
            continue
        ax2.axvline(freq, color=COLORS[fam], ls="-" if fam == "D" else "--", lw=1.0 if fam == "D" else 0.7, alpha=0.45)

    # Annotate the Jacobi-Anger folding cascade 2T -> T+D -> 2D -> 1-2D (P14)
    ax2.annotate(
        "Bessel + LTE folding cascade:\n2T (Mf) → T+D → 2D (6.41y) → 1-2D (P14, 14.2m)",
        xy=(0.156, float(P_mod_band[(f >= 0.145) & (f <= 0.165)].max())),
        xytext=(0.31, max(float(P_pre_band[fm].max()), float(P_mod_band[fm].max())) * 0.72),
        arrowprops=dict(arrowstyle="->", color="#c62828", lw=1.4),
        fontsize=8.8, fontweight="bold", color="#b71c1c",
        bbox=dict(boxstyle="round,pad=0.3", fc="#fff8e1", ec="#c62828", alpha=0.92),
    )

    ax2.set_xlabel("Frequency (cycles / year)", fontsize=10.5)
    ax2.set_ylabel("In-Band Power (%)", fontsize=10.5)
    ax2.set_xlim(0.02, args.fmax)
    ax2.grid(True, alpha=0.22)
    ax2.legend(loc="upper right", fontsize=8.5, framealpha=0.90)

    # Secondary top period axis on Panel 1
    ax_top = ax1.twiny()
    ax_top.set_xlim(ax1.get_xlim())
    p_ticks = [20, 10, 6.4, 4.8, 3.8, 2.4, 1.73, 1.36, 1.186, 1.0]
    ax_top.set_xticks([1.0 / p for p in p_ticks if 0.02 <= 1.0 / p <= args.fmax])
    ax_top.set_xticklabels([f"{p:g}y" if p != 1.186 else "14.2m" for p in p_ticks if 0.02 <= 1.0 / p <= args.fmax], fontsize=8.2)
    ax_top.set_xlabel("Period", fontsize=9.5, labelpad=4)

    out_path = Path(args.out) if args.out else (d / "baltic_aliasing_spectrum.png")
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")
    print(f"  Draconic 2D (0.14-0.19) / (2T,P) (0.24-0.28): Pre-Bessel F={r_pre:.2f}, LTE Model={r_mod:.2f}, Baltic Data={r_dat:.2f}")

if __name__ == "__main__":
    main()
