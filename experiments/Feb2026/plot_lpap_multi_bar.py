#!/usr/bin/env python3
"""
plot_lpap_multi_bar.py

A command-line tool to make a horizontal multi-bar-chart of the absolute
amplitudes of each tidal factor (lpap in lt.exe.p) across a list of subdirectories.

--dlod [FILE] compares lpap with the dLOD-calibrated values: the dLOD set is added as one more dataset (default
FILE: experiments/wisdom_lpap/dlod_calibrated.lt.exe.p; with no subdirectories given, the current folder's lt.exe.p
is compared), and two panels are added beside the amplitudes: the amplitude ratio to dLOD (log2: +1 = doubled) and
the phase shift from dLOD (degrees). Negative amplitudes are folded into the phase (+180 deg) before comparing.
The comparison is also printed as a table.
"""

import os
import sys
import json
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

DLOD_DEFAULT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "wisdom_lpap",
                                            "dlod_calibrated.lt.exe.p"))

def load_lpap(target_dir: str) -> list:
    """
    Load the 'lpap' list from 'lt.exe.p' in *target_dir* (or from *target_dir* itself when it is a file).
    """
    path = target_dir if os.path.isfile(target_dir) else os.path.join(target_dir, "lt.exe.p")
    if not os.path.isfile(path):
        raise FileNotFoundError(f"'lt.exe.p' not found in: {target_dir}")

    with open(path, "r", encoding="utf-8") as fh:
        data = json.load(fh)

    if isinstance(data, dict):
        if "lpap" not in data:
            raise KeyError(f"Key 'lpap' not found in {path} JSON object.")
        lpap = data["lpap"]
    elif isinstance(data, list):
        lpap = data
    else:
        raise TypeError(f"{path} must contain a JSON object or a JSON array.")

    if not lpap:
        raise ValueError(f"The 'lpap' list is empty in {path}.")

    return lpap

def format_period(p: float) -> str:
    """
    Format a period in days, adding approximate month or year indicators
    for longer periods.
    """
    if p >= 365.25:
        years = p / 365.25
        return f"{p:.5g} d ({years:.1f} yr)"
    elif p >= 30.0:
        months = p / 30.4375
        return f"{p:.5g} d ({months:.1f} mo)"
    else:
        return f"{p:.5g} d"

def cluster_periods(periods: set, tol: float = 1e-5) -> dict:
    """
    Groups raw periods that are within relative tolerance of each other,
    mapping each raw period to the average period of its group.
    """
    if not periods:
        return {}
    sorted_p = sorted(list(periods))
    clusters = []
    current_cluster = [sorted_p[0]]
    
    for p in sorted_p[1:]:
        ref = current_cluster[0]
        denom = abs(ref)
        if denom == 0:
            diff = abs(p - ref)
        else:
            diff = abs(p - ref) / denom
            
        if diff <= tol:
            current_cluster.append(p)
        else:
            clusters.append(current_cluster)
            current_cluster = [p]
    if current_cluster:
        clusters.append(current_cluster)
        
    mapping = {}
    for cluster in clusters:
        rep = sum(cluster) / len(cluster)
        for p in cluster:
            mapping[p] = rep
    return mapping

def discover_subdirs() -> list:
    """
    Find all subdirectories under the current directory containing 'lt.exe.p',
    skipping standard non-data directories.
    """
    SKIP = {"locs", "scripts", "rlr_data", "__pycache__", ".git", ".qwen"}
    subdirs = []
    for name in sorted(os.listdir(".")):
        if os.path.isdir(name) and name not in SKIP and not name.startswith(".") and not name.startswith("qbo"):
            if os.path.isfile(os.path.join(name, "lt.exe.p")):
                subdirs.append(name)
    return subdirs

def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot horizontal multi-bar chart of absolute tidal amplitudes (lpap in lt.exe.p)."
    )
    parser.add_argument(
        "subdirs",
        nargs="*",
        help="Subdirectories containing 'lt.exe.p'. If omitted, all relevant subdirs in the current folder are discovered automatically."
    )
    parser.add_argument(
        "--mode",
        choices=["grouped", "overlay"],
        default="grouped",
        help="Plot mode: 'grouped' (side-by-side bars) or 'overlay' (outlined overlaid bars). (default: %(default)s)"
    )
    parser.add_argument(
        "--min-amp",
        type=float,
        default=1e-5,
        help="Minimum absolute amplitude to display. Tides with smaller amplitudes will be filtered out. (default: %(default)s)"
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=1e-4,
        help="Relative tolerance to align/merge close periods across different datasets. (default: %(default)s)"
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default="lpap_multi_bar_chart.png",
        help="Filename to save the plot. Set to 'none' to skip saving. (default: %(default)s)"
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Do not display the GUI plot window (useful for automated scripts or headless environments)."
    )
    parser.add_argument(
        "--figsize",
        type=float,
        nargs=2,
        metavar=("WIDTH", "HEIGHT"),
        help="Custom figure size (width height) in inches. If omitted, height scales automatically with number of periods."
    )
    parser.add_argument(
        "--dlod",
        nargs="?",
        const=DLOD_DEFAULT,
        default=None,
        metavar="FILE",
        help="Compare with the dLOD-calibrated lpap (default file: %(const)s). Adds amplitude-ratio and phase-shift "
             "panels and prints the comparison. With no subdirectories, the current folder's lt.exe.p is compared."
    )
    parser.add_argument(
        "--log",
        action="store_true",
        help="Logarithmic amplitude axis (the minor lines are otherwise invisible next to Mf)."
    )
    parser.add_argument(
        "--short-labels",
        action="store_true",
        help="Use short labels (directory base names) in the legend instead of full paths."
    )
    return parser.parse_args()

def main():
    args = parse_args()
    
    # 1. Discover subdirs if not provided
    subdirs = args.subdirs
    if args.dlod and not subdirs and os.path.isfile("lt.exe.p"):
        subdirs = ["."]
    if args.dlod and not os.path.isfile(args.dlod):
        print(f"Error: dLOD-calibrated lpap not found: {args.dlod}", file=sys.stderr)
        sys.exit(1)
    if args.dlod and args.output == "lpap_multi_bar_chart.png":
        args.output = "lpap_vs_dlod.png"
    if not subdirs:
        subdirs = discover_subdirs()
        if not subdirs:
            print("Error: No subdirectories with 'lt.exe.p' found in current directory.", file=sys.stderr)
            print("Please specify subdirectories explicitly or run from the directory containing them.", file=sys.stderr)
            sys.exit(1)
        print(f"Automatically discovered {len(subdirs)} subdirectories containing 'lt.exe.p'.")
    else:
        # Validate provided subdirs
        valid_subdirs = []
        for sd in subdirs:
            if os.path.isdir(sd):
                if os.path.isfile(os.path.join(sd, "lt.exe.p")):
                    valid_subdirs.append(sd)
                else:
                    print(f"Warning: 'lt.exe.p' not found in subdirectory '{sd}'. Skipping.", file=sys.stderr)
            else:
                print(f"Warning: '{sd}' is not a valid directory. Skipping.", file=sys.stderr)
        subdirs = valid_subdirs
        if not subdirs:
            print("Error: No valid subdirectories with 'lt.exe.p' provided.", file=sys.stderr)
            sys.exit(1)

    # 2. Load and filter lpap from each subdirectory
    subdirs_data = {}
    all_raw_periods = set()
    phases = {}   # dataset -> {raw period: folded phase}; used by --dlod
    DLOD_LABEL = "dLOD calibrated"
    sources = [(sd, sd) for sd in subdirs] + ([(DLOD_LABEL, args.dlod)] if args.dlod else [])

    for sd, src in sources:
        try:
            lpap = load_lpap(src)
            filtered = []
            phases[sd] = {}
            for triplet in lpap:
                if len(triplet) < 2:
                    continue
                period = float(triplet[0])
                amp = abs(float(triplet[1]))
                if len(triplet) >= 3:   # fold a negative amplitude into the phase
                    phases[sd][period] = (float(triplet[2]) + (np.pi if float(triplet[1]) < 0 else 0.0)) % (2 * np.pi)
                if amp >= args.min_amp:
                    filtered.append((period, amp))
                    all_raw_periods.add(period)
            if filtered:
                subdirs_data[sd] = filtered
            else:
                print(f"Note: All tidal factors in '{sd}' were below min-amp threshold ({args.min_amp}).")
        except Exception as e:
            print(f"Warning: Failed to load/parse {sd}: {e}", file=sys.stderr)
            
    if not subdirs_data:
        print("Error: No tidal amplitude data met the filtering criteria in any directory.", file=sys.stderr)
        sys.exit(1)
        
    # 3. Cluster close periods
    period_mapping = cluster_periods(all_raw_periods, args.tol)
    
    # 4. Aggregate amplitudes by clustered periods
    subdir_amplitudes = {}
    unique_clustered_periods = set()
    
    for sd, data in subdirs_data.items():
        sd_amps = {}
        for raw_p, amp in data:
            clustered_p = period_mapping[raw_p]
            unique_clustered_periods.add(clustered_p)
            # Take maximum amplitude if multiple raw periods map to the same cluster
            if clustered_p in sd_amps:
                sd_amps[clustered_p] = max(sd_amps[clustered_p], amp)
            else:
                sd_amps[clustered_p] = amp
        subdir_amplitudes[sd] = sd_amps
        
    sorted_periods = sorted(list(unique_clustered_periods))
    n_periods = len(sorted_periods)
    print(f"Plotting {n_periods} unique period bands across {len(subdir_amplitudes)} datasets.")
    
    # Map clustered periods to integer y-indices
    y_pos_map = {p: i for i, p in enumerate(sorted_periods)}
    y_pos = np.arange(n_periods)
    
    # 5. Set up the figure
    if args.figsize:
        figsize = tuple(args.figsize)
    else:
        # Scale figure height automatically based on number of periods
        figsize = (12, max(6, 0.28 * n_periods))
        
    if args.dlod:
        if not args.figsize:
            figsize = (17, figsize[1])
        fig, (ax, axr, axp) = plt.subplots(1, 3, figsize=figsize, sharey=True, gridspec_kw={"width_ratios": [3, 1.4, 1.2]})
    else:
        fig, ax = plt.subplots(figsize=figsize)

    # Color selection
    sd_list = list(subdir_amplitudes.keys())
    M = len(sd_list)
    if M <= 10:
        colors = plt.cm.tab10(np.linspace(0, 1, M))
    elif M <= 20:
        colors = plt.cm.tab20(np.linspace(0, 1, M))
    else:
        colors = plt.cm.viridis(np.linspace(0, 1, M))
        
    # 6. Plotting
    if args.mode == "grouped":
        # Multi-bar grouped side-by-side horizontal chart
        total_height = 0.8
        bar_height = total_height / M
        
        for j, sd in enumerate(sd_list):
            sd_amps = subdir_amplitudes[sd]
            ys = []
            amps = []
            for p in sorted_periods:
                if p in sd_amps:
                    y_idx = y_pos_map[p]
                    y_offset = (j - (M - 1) / 2.0) * bar_height
                    ys.append(y_idx + y_offset)
                    amps.append(sd_amps[p])
            
            label = os.path.basename(os.path.abspath(sd)) if (args.short_labels or sd == ".") else sd
            ax.barh(
                ys, amps,
                height=bar_height,
                align="center",
                color=colors[j],
                label=label,
                alpha=0.85
            )
    else:
        # Overlay mode: outlined bars with no fill
        for j, sd in enumerate(sd_list):
            sd_amps = subdir_amplitudes[sd]
            ys = []
            amps = []
            for p in sorted_periods:
                if p in sd_amps:
                    ys.append(y_pos_map[p])
                    amps.append(sd_amps[p])
            
            label = os.path.basename(sd) if args.short_labels else sd
            ax.barh(
                ys, amps,
                align="center",
                fill=False,
                edgecolor=colors[j],
                linewidth=1.2,
                label=label,
                alpha=0.8
            )
            
    # 7. Customize axes and styling
    labels = [format_period(p) for p in sorted_periods]
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()  # Smallest period at top, largest at bottom
    
    # Format X-axis with 5 significant figures/nice digits
    def _sig5(x, _pos):
        if x == 0:
            return "0"
        from math import log10, floor
        try:
            digits = 5 - 1 - floor(log10(abs(x)))
            digits = max(0, digits)
            return f"{x:.{digits}f}"
        except ValueError:
            return f"{x}"
            
    if args.log:
        ax.set_xscale("log")
    else:
        ax.xaxis.set_major_formatter(FuncFormatter(_sig5))

    # 7b. --dlod: amplitude ratio and phase shift of every other dataset against the dLOD calibration
    if args.dlod and DLOD_LABEL in subdir_amplitudes:
        ref_amp = subdir_amplitudes[DLOD_LABEL]
        ref_ph = {period_mapping[p]: v for p, v in phases[DLOD_LABEL].items() if p in period_mapping}
        others = [sd for sd in sd_list if sd != DLOD_LABEL]
        K = max(len(others), 1); hgt = 0.8 / K
        print(f"\n{'period (d)':>18s} {'dLOD amp':>10s} " + " ".join(f"{(os.path.basename(os.path.abspath(sd)) if args.short_labels or sd == '.' else sd)[:22]:>34s}" for sd in others))
        table = {p: [] for p in sorted_periods}
        for j, sd in enumerate(others):
            col = colors[sd_list.index(sd)]
            sph = {period_mapping[p]: v for p, v in phases[sd].items() if p in period_mapping}
            ys, rat, yp, dph = [], [], [], []
            for p in sorted_periods:
                y = y_pos_map[p] + (j - (K - 1) / 2.0) * hgt
                a, r = subdir_amplitudes[sd].get(p), ref_amp.get(p)
                if a and r:
                    ys.append(y); rat.append(np.log2(a / r))
                    if p in sph and p in ref_ph:
                        yp.append(y_pos_map[p]); dph.append(np.degrees(np.angle(np.exp(1j * (sph[p] - ref_ph[p])))))
                table[p].append((a, r, np.log2(a / r) if a and r else None, dph[-1] if (a and r and yp and yp[-1] == y_pos_map[p]) else None))
            axr.barh(ys, rat, height=hgt, color=col, alpha=0.85)
            axp.plot(dph, yp, "o", ms=4, color=col)
        for p, rows in table.items():
            r = ref_amp.get(p)
            print(f"{format_period(p):>18s} {r if r else 0:10.5f} " + " ".join(
                f"{(a or 0):10.5f} {('x%.2f' % 2 ** q) if q is not None else '   -':>9s} {('%+5.0f deg' % d) if d is not None else '':>11s}" for a, _, q, d in rows))
        axr.axvline(0, color="#52514e", lw=0.8); axr.grid(True, axis="x", ls=":", lw=0.5, alpha=0.7)
        axr.set_xlabel("amplitude / dLOD (log2: +1 = doubled)", fontsize=10, fontweight="bold")
        axp.axvline(0, color="#52514e", lw=0.8); axp.set_xlim(-180, 180); axp.set_xticks([-180, -90, 0, 90, 180])
        axp.grid(True, ls=":", lw=0.5, alpha=0.7); axp.set_xlabel("phase - dLOD (deg)", fontsize=10, fontweight="bold")
        axr.set_title("vs dLOD", fontsize=11, fontweight="bold"); axp.set_title("phase shift", fontsize=11, fontweight="bold")
    
    ax.set_xlabel("Absolute Amplitude (a.u.)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Tidal Period (days)", fontsize=11, fontweight="bold")
    
    title_mode = "Grouped" if args.mode == "grouped" else "Overlaid"
    ax.set_title(
        (f"Tidal amplitudes vs the dLOD calibration ({os.path.basename(args.dlod)})\n" if args.dlod else
         f"Tidal Amplitude Spectrum — {title_mode} Multi-Bar Chart ({len(subdir_amplitudes)} datasets)\n")
        + f"Filtered for Abs Amplitude >= {args.min_amp}",
        fontsize=12,
        fontweight="bold"
    )
    
    ax.grid(True, axis="x", ls=":", lw=0.5, alpha=0.7)
    
    # Legend placement
    ax.legend(loc="upper right", bbox_to_anchor=(1.0, 1.0), framealpha=0.95, fontsize=9)
    
    fig.tight_layout()
    
    # 8. Save output
    if args.output.lower() != "none":
        try:
            plt.savefig(args.output, dpi=300, bbox_inches="tight")
            print(f"Plot successfully saved to: {args.output}")
        except Exception as e:
            print(f"Error saving plot to '{args.output}': {e}", file=sys.stderr)
            
    # 9. Show plot
    if not args.no_show:
        # Check if running in a headless environment
        if 'DISPLAY' not in os.environ and os.name != 'nt':
            print("Notice: No DISPLAY environment variable detected. Skipping GUI display.", file=sys.stderr)
        else:
            try:
                plt.show()
            except Exception as e:
                print(f"Error displaying GUI plot: {e}", file=sys.stderr)

if __name__ == "__main__":
    main()
