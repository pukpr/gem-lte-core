#!/usr/bin/env python3
"""extract_named_box.py -- pull one arbitrary (non-20-degree-grid-aligned)
lat/lon box out of the Kaplan SST NetCDF, cos(latitude)-weighted mean over
valid ocean cells, same convention as build_k_sst.py's extract_regions()
but for a box shaped to match a REAL named sub-basin instead of the coarse
20x20 macro-grid -- e.g. kN020_E050 (10-30N, 40-60E) lumps together the
Red Sea, the Persian Gulf, and a chunk of the Gulf of Aden/Arabian Sea in
one box; this pulls each one out on its own.

Usage:
    extract_named_box.py <name> <lat_lo> <lat_hi> <lon_lo> <lon_hi>
        [--min-year 1856]

Writes <name>/<name>.dat (year, value; project convention, monthly).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import netCDF4 as nc

HERE = Path(__file__).resolve().parent
NC_FILE = HERE / "k_sst" / "sst.mean.anom.nc"


def extract(name: str, lat_lo: float, lat_hi: float, lon_lo: float, lon_hi: float,
            min_year: float) -> None:
    ds = nc.Dataset(NC_FILE)
    lat = ds.variables["lat"][:].astype(float)
    lon = ds.variables["lon"][:].astype(float)
    time_var = ds.variables["time"]
    dates_cal = nc.num2date(time_var[:], time_var.units)
    proj_dates = np.array([d.year + (d.month - 1) / 12.0 for d in dates_cal])
    sst = ds.variables["sst"][:]

    coslat = np.cos(np.deg2rad(lat))
    lon_signed = np.where(lon > 180.0, lon - 360.0, lon)

    lat_sel = np.nonzero((lat >= lat_lo) & (lat < lat_hi))[0]
    lon_sel = np.nonzero((lon_signed >= lon_lo) & (lon_signed < lon_hi))[0]
    if len(lat_sel) == 0 or len(lon_sel) == 0:
        raise SystemExit(f"no grid cells in box lat[{lat_lo},{lat_hi}) "
                          f"lon[{lon_lo},{lon_hi})")

    box = sst[:, lat_sel[:, None], lon_sel[None, :]]
    mask = np.ma.getmaskarray(box)
    valid_frac = float((~mask).sum()) / mask.size
    w = coslat[lat_sel][:, None] * np.ones((1, len(lon_sel)))
    w = np.broadcast_to(w, box.shape[1:])

    series = np.full(len(proj_dates), np.nan)
    n_cells_used = np.full(len(proj_dates), 0)
    for t in range(box.shape[0]):
        bt = box[t]
        mt = ~np.ma.getmaskarray(bt)
        if mt.sum() == 0:
            continue
        series[t] = float(np.sum(np.where(mt, bt.filled(0.0) * w, 0.0)) /
                           np.sum(np.where(mt, w, 0.0)))
        n_cells_used[t] = int(mt.sum())

    ok = ~np.isnan(series) & (proj_dates >= min_year)
    print(f"[{name}] box lat[{lat_lo},{lat_hi}) lon[{lon_lo},{lon_hi})  "
          f"{len(lat_sel)}x{len(lon_sel)} cells, valid_frac={valid_frac:.2f}, "
          f"median cells/month used={int(np.median(n_cells_used[ok]))}, "
          f"n_months={ok.sum()}  span {proj_dates[ok].min():.1f}-{proj_dates[ok].max():.1f}")

    out_dir = HERE / name
    out_dir.mkdir(exist_ok=True)
    out = np.column_stack([proj_dates[ok], series[ok]])
    np.savetxt(out_dir / f"{name}.dat", out, fmt="%.6f\t%.6f")
    print(f"wrote {out_dir / (name + '.dat')}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("name")
    ap.add_argument("lat_lo", type=float)
    ap.add_argument("lat_hi", type=float)
    ap.add_argument("lon_lo", type=float)
    ap.add_argument("lon_hi", type=float)
    ap.add_argument("--min-year", type=float, default=1856.0)
    args = ap.parse_args()
    extract(args.name, args.lat_lo, args.lat_hi, args.lon_lo, args.lon_hi, args.min_year)


if __name__ == "__main__":
    main()
