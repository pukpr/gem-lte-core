#!/usr/bin/env python3
"""Cluster the 89 k_sst quads by their winding amplitude spectrum
(lt.exe.windings.json's k_amp_phase triples: k = winding value, amp =
its amplitude, phase = its phase) and plot each resulting region's
aggregate amplitude-by-winding histogram as a mini bar chart inset on a
Mercator map, at that region's geographic centroid.

Per the user (2026-09-28): "(k,amp) values change gradually across the
globe" -- so region boundaries are found by CLUSTERING the actual
per-quad amplitude spectra (not hand-drawn boxes like "north pacific").
Binning uses |k| (winding values are treated as unsigned magnitudes,
consistent with how every backbone/winding value elsewhere in this
project is discussed -- e.g. ~0.207, ~1.25, ~4.73 -- sign is a phase-
convention artifact, not a distinct physical winding) in 0.1-wide
intervals per the user's instruction, with a final overflow bin folding
in everything above WINDING_CUTOFF (the vast majority of amplitude mass
sits under 5.0; a long thin tail runs out past 20).

Cluster count is chosen automatically via silhouette score, not
hardcoded, since the user asked to "estimate how best to separate the
regions" rather than assume a fixed number of named ocean basins.
"""
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, "/home/paul/eval/gem-lte-core/experiments/Sep2026")
import sweep

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import normalize, StandardScaler

GEO_WEIGHT = 2.5  # relative weight of geographic compactness vs winding-
                   # shape similarity in the joint clustering feature
                   # space (both standardized first) -- higher = more
                   # spatially compact "named ocean basin"-style regions,
                   # lower = purer winding-shape grouping regardless of
                   # location. Chosen empirically: weight 1.0 gave
                   # geographically-scattered clusters (confirmed
                   # visually this session -- same-color dots spread
                   # pole to pole), so spatial compactness needs to
                   # dominate for the output to resemble the "north
                   # pacific, equatorial pacific" style regions the user
                   # described, while the winding features still get a
                   # real vote in where the boundaries fall.

NAME_RE = re.compile(r"^k([NS])(\d{3})_([EW])(\d{3})$")
WINDING_CUTOFF = 5.0
BIN_WIDTH = 0.1
N_BINS = int(WINDING_CUTOFF / BIN_WIDTH) + 1  # last bin is the overflow bin

# Categorical palette (dataviz skill's validated default, light mode)
CLUSTER_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4",
                   "#008300", "#4a3aa7", "#e34948"]


def decode(name: str) -> tuple[float, float]:
    m = NAME_RE.match(name)
    ns, lat_s, ew, lon_s = m.groups()
    return float(lat_s) * (1 if ns == "N" else -1), float(lon_s) * (1 if ew == "E" else -1)


def freshest_windings(cell: str) -> Path | None:
    feb = sweep.FEB / cell / "lt.exe.windings.json"
    sep = sweep.HERE / cell / "lt.exe.windings.json"
    candidates = [p for p in (feb, sep) if p.exists()]
    if not candidates:
        return None
    return max(candidates, key=lambda p: p.stat().st_mtime)


def amplitude_histogram(triples: list[list[float]]) -> np.ndarray:
    hist = np.zeros(N_BINS)
    for k, amp, phase in triples:
        ak = abs(k)
        bin_idx = min(int(ak / BIN_WIDTH), N_BINS - 1)
        hist[bin_idx] += amp
    return hist


def geo_centroid(member_coords: np.ndarray) -> tuple[float, float]:
    """Mean lat/lon that handles longitude wraparound correctly (circular
    mean via atan2) -- naive averaging breaks for any region straddling
    the dateline, e.g. a cluster spanning both E170 and W170 would
    average to ~0 (India), not ~180 (the Pacific) where it belongs.
    Confirmed this was happening for the wide Pacific-spanning cluster
    before this fix."""
    lat = member_coords[:, 0].mean()
    lon_rad = np.radians(member_coords[:, 1])
    lon = np.degrees(np.arctan2(np.sin(lon_rad).mean(), np.cos(lon_rad).mean()))
    return lat, lon


def choose_k(features: np.ndarray, k_range: range) -> tuple[int, KMeans]:
    best_k, best_score, best_model = None, -1.0, None
    for k in k_range:
        model = KMeans(n_clusters=k, n_init=10, random_state=0).fit(features)
        if len(set(model.labels_)) < 2:
            continue
        score = silhouette_score(features, model.labels_)
        print(f"  k={k}: silhouette={score:.3f}")
        if score > best_score:
            best_k, best_score, best_model = k, score, model
    return best_k, best_model


def main() -> None:
    cells = sorted((p.name for p in sweep.FEB.iterdir()
                     if p.is_dir() and NAME_RE.match(p.name)))
    histograms, coords, used_cells = [], [], []
    for cell in cells:
        path = freshest_windings(cell)
        if path is None:
            continue
        data = json.loads(path.read_text())
        triples = data.get("k_amp_phase", [])
        if not triples:
            continue
        histograms.append(amplitude_histogram(triples))
        coords.append(decode(cell))
        used_cells.append(cell)

    H = np.array(histograms)
    coords = np.array(coords)
    print(f"loaded {len(used_cells)} quads' winding spectra")

    # Cluster on L1-normalized spectra (shape, not overall amplitude
    # magnitude, is what should define a "region" here -- two quads with
    # the same winding-mix but different total energy should cluster
    # together) COMBINED with geographic position, standardized and
    # weighted (see GEO_WEIGHT) -- pure winding-shape clustering alone
    # produced geographically scattered "regions" (confirmed visually:
    # same-color quads from the Arctic to the Southern Ocean), which
    # doesn't match what "region" means for this map.
    H_norm = normalize(H, norm="l1")
    wind_features = StandardScaler().fit_transform(H_norm)

    lat = coords[:, 0]
    lon_rad = np.radians(coords[:, 1])
    geo_raw = np.column_stack([lat, np.sin(lon_rad) * 90, np.cos(lon_rad) * 90])
    geo_features = StandardScaler().fit_transform(geo_raw) * GEO_WEIGHT

    joint_features = np.hstack([geo_features, wind_features])

    print("Choosing cluster count via silhouette score:")
    best_k, model = choose_k(joint_features, range(3, 10))
    labels = model.labels_
    print(f"chosen k={best_k}")

    fig = plt.figure(figsize=(26, 15))
    ax = plt.axes(projection=ccrs.Mercator())
    ax.set_extent([-180, 180, -70, 80], crs=ccrs.PlateCarree())
    ax.coastlines(linewidth=0.6)
    ax.add_feature(cfeature.LAND, edgecolor="black", facecolor="#f5f4f0", linewidth=0.4)
    ax.add_feature(cfeature.OCEAN, facecolor="#eaf3fb")

    for i, cell in enumerate(used_cells):
        lat, lon = coords[i]
        color = CLUSTER_COLORS[labels[i] % len(CLUSTER_COLORS)]
        ax.plot(lon, lat, "o", color=color, markersize=5,
                transform=ccrs.PlateCarree(), markeredgecolor="white",
                markeredgewidth=0.4)

    bin_centers = np.arange(N_BINS) * BIN_WIDTH + BIN_WIDTH / 2

    for cluster_id in range(best_k):
        members = labels == cluster_id
        if not members.any():
            continue
        region_hist = H[members].sum(axis=0)
        c_lat, c_lon = geo_centroid(coords[members])
        color = CLUSTER_COLORS[cluster_id % len(CLUSTER_COLORS)]

        px, py = ax.projection.transform_point(c_lon, c_lat, ccrs.PlateCarree())
        px_disp, py_disp = ax.transData.transform((px, py))
        px_fig, py_fig = fig.transFigure.inverted().transform((px_disp, py_disp))

        w, h = 0.14, 0.13
        inset = fig.add_axes((max(0.01, min(0.99 - w, px_fig - w / 2)),
                               max(0.01, min(0.99 - h, py_fig - h / 2)),
                               w, h))
        inset.bar(bin_centers[:-1], region_hist[:-1], width=BIN_WIDTH * 0.9,
                  color=color, edgecolor="none")
        if region_hist[-1] > 0:
            inset.bar([WINDING_CUTOFF + BIN_WIDTH], [region_hist[-1]],
                       width=BIN_WIDTH * 0.9, color=color, alpha=0.4)
        inset.set_title(f"region {cluster_id} (n={members.sum()})",
                         fontsize=7, color=color, fontweight="bold")
        inset.tick_params(labelsize=5)
        inset.set_xlim(0, WINDING_CUTOFF + 2 * BIN_WIDTH)
        inset.patch.set_alpha(0.9)
        for spine in inset.spines.values():
            spine.set_edgecolor(color)

    ax.set_title(
        f"Winding-amplitude spectrum by region — {len(used_cells)} k_sst quads, "
        f"clustered into {best_k} regions (KMeans on L1-normalized |k| histograms, "
        f"0.1-wide bins, silhouette-selected k)\nDot color = cluster membership; "
        f"bar charts = each region's total amplitude per winding value "
        f"(final bin = overflow, |k|>{WINDING_CUTOFF})",
        fontsize=15)

    plt.tight_layout()
    out = sweep.HERE / "winding_regions_map.png"
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")

    print("\nCluster membership:")
    for cluster_id in range(best_k):
        members = [used_cells[i] for i in range(len(used_cells)) if labels[i] == cluster_id]
        print(f"  region {cluster_id} (n={len(members)}): {members}")


if __name__ == "__main__":
    main()
