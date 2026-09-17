#!/usr/bin/env python3
"""
plot_lpap_multi_bar.py

A command-line tool to make a horizontal multi-bar-chart of the absolute
amplitudes of each tidal factor (lpap in lt.exe.p) across a list of subdirectories.
"""

import os
import sys
import json
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

def load_lpap(target_dir: str) -> list:
    """
    Load the 'lpap' list from 'lt.exe.p' in *target_dir*.
    """
    path = os.path.join(target_dir, "lt.exe.p")
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
        "--short-labels",
        action="store_true",
        help="Use short labels (directory base names) in the legend instead of full paths."
    )
    return parser.parse_args()

def main():
    args = parse_args()
    
    # 1. Discover subdirs if not provided
    subdirs = args.subdirs
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
    
    for sd in subdirs:
        try:
            lpap = load_lpap(sd)
            filtered = []
            for triplet in lpap:
                if len(triplet) < 2:
                    continue
                period = float(triplet[0])
                amp = abs(float(triplet[1]))
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
            
            label = os.path.basename(sd) if args.short_labels else sd
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
            
    ax.xaxis.set_major_formatter(FuncFormatter(_sig5))
    
    ax.set_xlabel("Absolute Amplitude (a.u.)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Tidal Period (days)", fontsize=11, fontweight="bold")
    
    title_mode = "Grouped" if args.mode == "grouped" else "Overlaid"
    ax.set_title(
        f"Tidal Amplitude Spectrum — {title_mode} Multi-Bar Chart ({len(subdir_amplitudes)} datasets)\n"
        f"Filtered for Abs Amplitude >= {args.min_amp}",
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
