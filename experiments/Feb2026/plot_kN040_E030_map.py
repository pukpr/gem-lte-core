#!/usr/bin/env python3
"""plot_kN040_E030_map.py -- geographic map of the standout Kaplan SST
region kN040_E030 (30-50N, 20-40E: Aegean/Black Sea/Eastern Mediterranean
approaches) that shows a strong, reproducible winding fit (val_r 0.60-0.70
across multiple record lengths/donors) despite not corresponding to any
named production climate index in this project.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

LAT_LO, LAT_HI = 30.0, 50.0
LON_LO, LON_HI = 20.0, 40.0

fig = plt.figure(figsize=(11, 7))

# panel 1: regional close-up
ax1 = fig.add_subplot(1, 2, 1, projection=ccrs.PlateCarree())
ax1.set_extent([LON_LO - 12, LON_HI + 10, LAT_LO - 8, LAT_HI + 8], crs=ccrs.PlateCarree())
ax1.add_feature(cfeature.LAND, facecolor="#ddd6c5", zorder=0)
ax1.add_feature(cfeature.OCEAN, facecolor="#bcd4e6", zorder=0)
ax1.add_feature(cfeature.COASTLINE, linewidth=0.6, edgecolor="0.3")
ax1.add_feature(cfeature.BORDERS, linewidth=0.3, edgecolor="0.5", linestyle=":")
ax1.gridlines(draw_labels=True, linewidth=0.3, color="0.6", alpha=0.6)
ax1.add_patch(plt.Rectangle((LON_LO, LAT_LO), LON_HI - LON_LO, LAT_HI - LAT_LO,
                             transform=ccrs.PlateCarree(), facecolor="none",
                             edgecolor="crimson", linewidth=2.2, zorder=5))
ax1.set_title("kN040_E030 region (30-50N, 20-40E)\nAegean / Sea of Marmara / "
              "Black Sea / Eastern Med approaches", fontsize=10)

# panel 2: global context
ax2 = fig.add_subplot(1, 2, 2, projection=ccrs.Robinson())
ax2.set_global()
ax2.add_feature(cfeature.LAND, facecolor="#ddd6c5", zorder=0)
ax2.add_feature(cfeature.OCEAN, facecolor="#bcd4e6", zorder=0)
ax2.add_feature(cfeature.COASTLINE, linewidth=0.4, edgecolor="0.4")
ax2.gridlines(linewidth=0.2, color="0.7", alpha=0.5)
ax2.add_patch(plt.Rectangle((LON_LO, LAT_LO), LON_HI - LON_LO, LAT_HI - LAT_LO,
                             transform=ccrs.PlateCarree(), facecolor="crimson",
                             edgecolor="crimson", linewidth=1.5, zorder=5, alpha=0.85))
ax2.set_title("global context", fontsize=10)

fig.suptitle("kN040_E030: strongest unclaimed winding candidate in the Kaplan "
             "SST sweep\nval_r = 0.60-0.70 across record lengths/donors, no "
             "overlap with any named production index (NINO4/AMO/PDO/TNA/IOD)",
             fontsize=11, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.92))

out = "/home/paul/eval/gem-lte-core/experiments/Feb2026/k_sst/kN040_E030_map.png"
fig.savefig(out, dpi=150)
print(f"[saved] {out}")
