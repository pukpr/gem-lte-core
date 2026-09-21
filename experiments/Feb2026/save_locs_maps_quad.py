#!/usr/bin/env python3
"""save_locs_maps_quad.py -- generate locs/<region>_loc.png for every
k_sst grid-quadrangle region (kN040_W030, kS020_E110, ...), adapted from
~/github/pukpr/pukpr.github.io/results/save_locs_maps.py (which plots a
single point). These regions are real 20x20 degree boxes
(build_k_sst.py's REGION_DEG), not point stations, so the map draws the
actual quadrangle outline instead of a point -- the name itself encodes
the box center (kN040_W030 -> lat_c=+40, lon_c=-30; box = center +/- 10
degrees), matching build_k_sst.extract_regions()'s naming convention
exactly.

Usage: save_locs_maps_quad.py [region ...]   (default: all matching
dirs under the Feb2026 root that have an lt.exe.p)
"""
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

HERE = Path(__file__).resolve().parent
LOCS_DIR = HERE / "locs"
REGION_DEG = 20.0
NAME_RE = re.compile(r"^k([NS])(\d{3})_([EW])(\d{3})$")


def decode_region(name: str):
    m = NAME_RE.match(name)
    if not m:
        return None
    ns, lat_s, ew, lon_s = m.groups()
    lat_c = float(lat_s) * (1 if ns == "N" else -1)
    lon_c = float(lon_s) * (1 if ew == "E" else -1)
    return lat_c, lon_c


def plot_quadrangle(name: str, lat_c: float, lon_c: float) -> None:
    half = REGION_DEG / 2.0
    lat_lo, lat_hi = lat_c - half, lat_c + half
    lon_lo, lon_hi = lon_c - half, lon_c + half

    fig = plt.figure(figsize=(10, 5))
    ax = plt.axes(projection=ccrs.Robinson())
    ax.set_global()
    ax.coastlines()
    ax.add_feature(cfeature.LAND, edgecolor="black", facecolor="none")
    ax.add_feature(cfeature.OCEAN, facecolor="white")

    box_lons = [lon_lo, lon_hi, lon_hi, lon_lo, lon_lo]
    box_lats = [lat_lo, lat_lo, lat_hi, lat_hi, lat_lo]
    ax.plot(box_lons, box_lats, color="red", linewidth=1.5,
            transform=ccrs.PlateCarree())
    ax.fill(box_lons, box_lats, color="red", alpha=0.35,
            transform=ccrs.PlateCarree())
    ax.plot(lon_c, lat_c, "ro", markersize=4, transform=ccrs.PlateCarree())
    label_lat = min(lat_hi + 3, 88.0)
    ax.text(lon_c, label_lat, name, fontsize=10,
            transform=ccrs.PlateCarree(), ha="center", va="bottom",
            color="blue")

    plt.title(f"Location: {name}  "
              f"[{lat_lo:+.0f}°..{lat_hi:+.0f}° lat, "
              f"{lon_lo:+.0f}°..{lon_hi:+.0f}° lon]")
    plt.tight_layout()
    LOCS_DIR.mkdir(exist_ok=True)
    out = LOCS_DIR / f"{name}_loc.png"
    plt.savefig(out)
    plt.close(fig)
    print(f"wrote {out}")


def discover_regions() -> list[str]:
    SKIP = {"locs", "scripts", "rlr_data", "node_modules", "__pycache__", "checkpoints"}
    names = []
    for p in sorted(HERE.iterdir()):
        if p.is_dir() and p.name not in SKIP and (p / "lt.exe.p").exists():
            if NAME_RE.match(p.name):
                names.append(p.name)
    return names


def main() -> None:
    names = sys.argv[1:] if len(sys.argv) > 1 else discover_regions()
    if not names:
        print("No matching k_sst grid regions found.")
        return
    for name in names:
        decoded = decode_region(name)
        if decoded is None:
            print(f"skip {name}: doesn't match kNxxx_Exxx / kSxxx_Wxxx pattern")
            continue
        lat_c, lon_c = decoded
        plot_quadrangle(name, lat_c, lon_c)


if __name__ == "__main__":
    main()
