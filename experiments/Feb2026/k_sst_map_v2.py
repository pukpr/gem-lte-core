#!/usr/bin/env python3
"""k_sst_map_v2.py — geographic map of the PROPERLY-GATED (winding_rank)
k_sst sweep results (k_sst_summary_v2.json), replacing the earlier
crude-method map. Two panels:
  (1) val_r (backbone + however many windings actually passed the gate,
      0-2 per region) by region, colored on a diverging scale.
  (2) number of passing region-specific windings (0, 1, or 2) — most
      regions are 0, so this is really a map of WHERE the shared
      backbone alone is insufficient and something extra was found.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
v2 = json.load(open(HERE / "k_sst" / "k_sst_summary_v2.json"))

lat = np.array([r["lat"] for r in v2])
lon = np.array([r["lon"] for r in v2])
val_r = np.array([r["val_r"] for r in v2])
n_pass = np.array([r["n_passing"] for r in v2])

fig, axes = plt.subplots(2, 1, figsize=(11, 9))

ax = axes[0]
sc = ax.scatter(lon, lat, c=val_r, cmap="RdBu_r", vmin=-0.4, vmax=0.7, s=260,
                 marker="s", edgecolor="0.3", linewidth=0.5)
for r in v2:
    if r["n_passing"] > 0:
        ax.text(r["lon"], r["lat"], "*", ha="center", va="center",
                fontsize=14, color="black", fontweight="bold")
cb = fig.colorbar(sc, ax=ax, label="val_r (holdout, shared-backbone model)")
n_star = int(np.sum(n_pass > 0))
n_ridges = int(np.sum(n_pass))
ax.set_title(f"Kaplan SST regions (AMO donor, >=1950): holdout r of the shared "
              f"backbone + any gated extra winding\n"
              f"'*' marks {n_star}/89 regions with >=1 passing ridge ({n_ridges} "
              f"total) -- NOT distinguishable from the calibrated null false-"
              f"positive rate (~27 expected by chance alone); no individual "
              f"ridge here should be trusted", fontsize=10)
ax.set_xlabel("lon"); ax.set_ylabel("lat")
ax.axhline(0, color="0.6", lw=0.5); ax.set_xlim(-190, 190); ax.set_ylim(-95, 95)
ax.set_facecolor("0.93")

ax = axes[1]
sc2 = ax.scatter(lon, lat, c=n_pass, cmap="viridis", vmin=0, vmax=2, s=260,
                  marker="s", edgecolor="0.3", linewidth=0.5)
fig.colorbar(sc2, ax=ax, label="# region-specific windings passing AR1+FWHM+continuity gate", ticks=[0, 1, 2])
ax.set_title("Where a genuine EXTRA standing-wave mode was found (90% of regions: none)", fontsize=10)
ax.set_xlabel("lon"); ax.set_ylabel("lat")
ax.axhline(0, color="0.6", lw=0.5); ax.set_xlim(-190, 190); ax.set_ylim(-95, 95)
ax.set_facecolor("0.93")

fig.tight_layout()
out = HERE / "k_sst" / "k_sst_map_v2.png"
fig.savefig(out, dpi=130)
print(f"[saved] {out}")

# also print the recurring-M cluster check
print("\n[recurring-M check] the 9 regions with n_passing>0:")
for r in sorted(v2, key=lambda r: -r["n_passing"]):
    if r["n_passing"] > 0:
        print(f"  {r['name']:12s} lat={r['lat']:+6.1f} lon={r['lon']:+6.1f}  "
              f"M={np.round(r['top_M'],3).tolist()}  val_r={r['val_r']:+.3f}")
