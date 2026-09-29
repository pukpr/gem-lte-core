#!/usr/bin/env python3
"""World Mercator status map of all 89 k_sst grid quads, colored by fit
status (validated status palette from the dataviz skill), each box
labeled with its train/validate/test scores (or pair values, or "?" if
never attempted by this sweep at all).

Categories (status palette, fixed roles):
  good      #0ca30c  -- accepted, no caveat, weakest score >= 0.3
  warning   #fab219  -- accepted, weak score in [0.0, 0.3) somewhere
  serious   #ec835a  -- accepted, a NEGATIVE score somewhere
  critical  #d03b3b  -- still bad / unsolved (accepted: false in ledger)
  gray      #9a9a9a  -- never attempted by this sweep (no ledger entry)
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, "/home/paul/eval/gem-lte-core/experiments/Sep2026")
import sweep

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import cartopy.crs as ccrs
import cartopy.feature as cfeature

FEB = sweep.FEB
NAME_RE = re.compile(r"^k([NS])(\d{3})_([EW])(\d{3})$")
REGION_DEG = 20.0

STATUS_COLORS = {
    "good": "#0ca30c",
    "warning": "#fab219",
    "serious": "#ec835a",
    "critical": "#d03b3b",
    "unattempted": "#9a9a9a",
}
STATUS_LABELS = {
    "good": "Solved — clean",
    "warning": "Solved — weak score (0–0.3)",
    "serious": "Solved — negative score",
    "critical": "Still unsolved",
    "unattempted": "Never attempted (no ledger entry)",
}


def decode(name: str):
    m = NAME_RE.match(name)
    if not m:
        return None
    ns, lat_s, ew, lon_s = m.groups()
    lat_c = float(lat_s) * (1 if ns == "N" else -1)
    lon_c = float(lon_s) * (1 if ew == "E" else -1)
    return lat_c, lon_c


def score_list(entry: dict) -> list[float]:
    """Flatten whatever score representation this ledger entry has
    (triplet train/validate/test, or a 2-element pair) into a plain
    list of floats for classification and display."""
    scores = entry.get("scores") or {}
    if "train" in scores and "validate" in scores and "test" in scores:
        return [scores["train"], scores["validate"], scores["test"]]
    if "pair" in scores:
        return list(scores["pair"])
    return []


def classify(entry: dict | None) -> str:
    if entry is None:
        return "unattempted"
    if not entry.get("accepted"):
        return "critical"
    vals = score_list(entry)
    if any(v < 0 for v in vals):
        return "serious"
    if any(0 <= v < 0.3 for v in vals):
        return "warning"
    return "good"


def label_text(name: str, entry: dict | None) -> str:
    if entry is None:
        return f"{name}\n(untouched)"
    if not entry.get("accepted"):
        return f"{name}\nSTILL BAD"
    vals = score_list(entry)
    if len(vals) == 3:
        return f"{name}\ntr={vals[0]:.2f} v={vals[1]:.2f} te={vals[2]:.2f}"
    if len(vals) == 2:
        return f"{name}\npair=({vals[0]:.2f},{vals[1]:.2f})"
    return f"{name}\n(no scores)"


def main() -> None:
    names = sorted(p.name for p in FEB.iterdir()
                    if p.is_dir() and NAME_RE.match(p.name))
    by_cell = sweep.latest_entry_per_cell(sweep.load_sweep_ledger())

    fig = plt.figure(figsize=(34, 20), dpi=150)
    ax = plt.axes(projection=ccrs.Mercator())
    ax.set_extent([-180, 180, -75, 80], crs=ccrs.PlateCarree())
    ax.coastlines(linewidth=0.6)
    ax.add_feature(cfeature.LAND, edgecolor="black", facecolor="#f5f4f0", linewidth=0.4)
    ax.add_feature(cfeature.OCEAN, facecolor="#eaf3fb")

    counts = {k: 0 for k in STATUS_COLORS}

    for name in names:
        decoded = decode(name)
        if decoded is None:
            continue
        lat_c, lon_c = decoded
        half = REGION_DEG / 2.0
        lat_lo, lat_hi = lat_c - half, lat_c + half
        lon_lo, lon_hi = lon_c - half, lon_c + half
        entry = by_cell.get(name)
        status = classify(entry)
        counts[status] += 1
        color = STATUS_COLORS[status]

        box_lons = [lon_lo, lon_hi, lon_hi, lon_lo, lon_lo]
        box_lats = [lat_lo, lat_lo, lat_hi, lat_hi, lat_lo]
        ax.plot(box_lons, box_lats, color=color, linewidth=1.2,
                transform=ccrs.PlateCarree())
        ax.fill(box_lons, box_lats, color=color, alpha=0.30,
                transform=ccrs.PlateCarree())

        ax.text(lon_c, lat_c, label_text(name, entry), fontsize=5.6,
                transform=ccrs.PlateCarree(), ha="center", va="center",
                color="#1a1a1a", linespacing=1.3,
                bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                          edgecolor=color, linewidth=0.8, alpha=0.85))

    legend_handles = [
        mpatches.Patch(facecolor=STATUS_COLORS[k], edgecolor=STATUS_COLORS[k],
                        alpha=0.6, label=f"{STATUS_LABELS[k]}  (n={counts[k]})")
        for k in ("good", "warning", "serious", "critical", "unattempted")
    ]
    ax.legend(handles=legend_handles, loc="lower left", fontsize=13,
              framealpha=0.95, title="Fit status (train/validate/test or pair)",
              title_fontsize=14)

    plt.title(
        f"GEM-LTE k_sst grid-quad fit status ({len(names)} quads, "
        f"{counts['good']+counts['warning']+counts['serious']} solved, "
        f"{counts['critical']} unsolved, {counts['unattempted']} untouched) "
        "— 2026-09-29",
        fontsize=20, pad=14)
    plt.tight_layout()
    out = sweep.HERE / "quad_status_map.png"
    plt.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")
    print("counts:", counts)


if __name__ == "__main__":
    main()
