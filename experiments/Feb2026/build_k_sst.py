#!/usr/bin/env python3
"""build_k_sst.py — chop the Kaplan Extended SST V2 global grid
(k_sst/sst.mean.anom.nc, downloaded from
https://downloads.psl.noaa.gov//Datasets/kaplan_sst/sst.mean.anom.nc)
into regional boxes, discover each region's own dominant winding
number(s) against a SHARED, borrowed forcing manifold (no per-region
Ada optimizer run available), write each region as a real ws.py-
compatible <region>/<region>.dat + <region>/lt.exe.p, then run ws.py's
own global-mode fit on every region and tabulate results.

Dataset: 5x5 degree global grid, monthly, 1856-01 to 2023-01 (2005
months), Kaplan optimal-analysis SST anomaly (already gap-filled by the
source product's own EOF-based reconstruction -- not fabricated here).

Region definition: REGION_DEG x REGION_DEG boxes (default 20x20),
cos(latitude)-weighted mean of all valid (non-land, non-ice-masked)
grid cells in the box, computed independently per month (so a box's
own series can still show a real missing month if literally nothing in
it was valid that month, though the source product rarely does this).
Boxes with fewer than MIN_OCEAN_FRAC valid cells on average are skipped
as insufficiently oceanic.

Shared manifold: every region reuses nino4's own tide/comb/IIR/Bessel
shape constants verbatim (offs, bg, impA, impB, delA, delB, asym, ma,
mp, shfT, init, year, the 42-row lpap table) -- justified by this
session's and the sibling repo's own finding that nearly all
independently-optimized indices converge on essentially ONE manifold
(manifold_overlay.py: 12/13 indices' own manifolds correlate >0.9976
with the mean). IR is set to 0 (no per-region delay-differential fit
attempted -- a real simplification, stated plainly, not hidden).

Winding discovery (since there is no per-region Ada optimizer run to
supply lt.exe.p's ltep): scan M over a grid, at each M fit a 4-column
regression [1, F, sin(2*pi*M*F), cos(2*pi*M*F)] (plus trend/seasonal)
against the region's own data, and compare its R^2 to the SAME
regression's R^2 on AR(1)-matched noise surrogates with the region's
own lag-1 autocorrelation and variance -- the region's own top winding
is the M with the largest excess (real R^2 minus surrogate mean R^2).
This is a simpler, whole-record version of winding_rank.py's own
AR(1)-floor discipline (not the full windowed peak/FWHM/continuity
triple test), appropriate for a first-pass sweep across many regions.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import netCDF4 as nc

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ws  # noqa: E402

K_SST_DIR = HERE / "k_sst"
NC_FILE = K_SST_DIR / "sst.mean.anom.nc"
DONOR_INDEX = "amo"
REGION_DEG = 20.0
MIN_OCEAN_FRAC = 0.30
MIN_YEAR = 1950.0  # Kaplan SST pre-1950 coverage is sparse/ship-track-only
                    # and not trustworthy at this grid resolution -- restrict
                    # the whole sweep (and the >=200-month coverage test) to
                    # the modern, better-observed record.
M_GRID = np.arange(0.03, 5.0001, 0.02)
N_AR1_SURROGATES = 16


def load_donor_params() -> dict:
    p = json.loads((HERE / DONOR_INDEX / "lt.exe.p").read_text())
    keep = ("offs", "bg", "impA", "impB", "impC", "delA", "delB", "asym",
            "ann1", "ann2", "year", "ma", "mp", "shfT", "init", "lpap")
    out = {k: p[k] for k in keep if k in p}
    # the shared manifold's own Bessel FM rate (its backbone winding) --
    # used BOTH to build F(t) for ridge discovery below AND appended as
    # the final ltep entry in every region's own written lt.exe.p, so the
    # SAME k1 drives the Bessel step in discovery and in the final fit --
    # avoids a circular/inconsistent-F bug (region-specific ltep can't be
    # known before discovery, but Bessel's k1 must be fixed to compute F
    # in the first place).
    out["backbone_k1"] = float(p["ltep"][-1])
    return out


def extract_regions() -> list[dict]:
    ds = nc.Dataset(NC_FILE)
    lat = ds.variables["lat"][:].astype(float)
    lon = ds.variables["lon"][:].astype(float)
    time_var = ds.variables["time"]
    dates_cal = nc.num2date(time_var[:], time_var.units)
    # project date convention: year + month_index/12 (month_index from
    # the calendar date's own month, 0-based) -- matches every other
    # index's .dat file in this project
    proj_dates = np.array([d.year + (d.month - 1) / 12.0 for d in dates_cal])
    sst = ds.variables["sst"][:]  # (time, lat, lon), masked array

    coslat = np.cos(np.deg2rad(lat))
    lat_edges = np.arange(-90.0, 90.0001, REGION_DEG)
    lon_edges = np.arange(0.0, 360.0001, REGION_DEG)

    regions = []
    for i in range(len(lat_edges) - 1):
        lat_lo, lat_hi = lat_edges[i], lat_edges[i + 1]
        lat_sel = np.nonzero((lat >= lat_lo) & (lat < lat_hi))[0]
        if len(lat_sel) == 0:
            continue
        for j in range(len(lon_edges) - 1):
            lon_lo, lon_hi = lon_edges[j], lon_edges[j + 1]
            lon_sel = np.nonzero((lon >= lon_lo) & (lon < lon_hi))[0]
            if len(lon_sel) == 0:
                continue
            box = sst[:, lat_sel[:, None], lon_sel[None, :]]  # (time, nlat, nlon)
            mask = np.ma.getmaskarray(box)
            valid_frac = float((~mask).sum()) / mask.size
            if valid_frac < MIN_OCEAN_FRAC:
                continue
            w = coslat[lat_sel][:, None] * np.ones((1, len(lon_sel)))
            w = np.broadcast_to(w, box.shape[1:])
            series = np.full(len(proj_dates), np.nan)
            for t in range(box.shape[0]):
                bt = box[t]
                mt = ~np.ma.getmaskarray(bt)
                if mt.sum() == 0:
                    continue
                series[t] = float(np.sum(np.where(mt, bt.filled(0.0) * w, 0.0)) /
                                   np.sum(np.where(mt, w, 0.0)))
            ok = ~np.isnan(series) & (proj_dates >= MIN_YEAR)
            if ok.sum() < 200:  # need a real amount of data to fit anything
                continue
            lat_c, lon_c = (lat_lo + lat_hi) / 2.0, (lon_lo + lon_hi) / 2.0
            lon_c_signed = lon_c - 360.0 if lon_c > 180.0 else lon_c
            name = f"k{'N' if lat_c >= 0 else 'S'}{abs(int(lat_c)):03d}_" \
                   f"{'E' if lon_c_signed >= 0 else 'W'}{abs(int(lon_c_signed)):03d}"
            regions.append(dict(name=name, lat_c=lat_c, lon_c=lon_c_signed,
                                 dates=proj_dates[ok], values=series[ok],
                                 valid_frac=valid_frac))
    return regions


def ar1_surrogate(x: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    x = x - x.mean()
    phi = np.corrcoef(x[:-1], x[1:])[0, 1]
    phi = np.clip(phi, -0.98, 0.98)
    sigma = x.std() * math.sqrt(max(1e-8, 1 - phi ** 2))
    n = len(x)
    s = np.empty(n)
    s[0] = rng.normal(0, x.std())
    for i in range(1, n):
        s[i] = phi * s[i - 1] + rng.normal(0, sigma)
    return s


def regression_r2(t, F, y, M, trend_on=True):
    # ANNUAL *and* SEMI-ANNUAL calendar harmonics are both included here,
    # matching ws.build_design's own convention -- without both, the M-scan
    # can mistake ordinary calendar-locked seasonal content (autonomous,
    # not tidal) for a genuine non-autonomous winding ridge. Raw SST
    # anomaly grids commonly still carry residual seasonal structure even
    # after climatology removal (the climatology itself has sampling
    # noise), so this is not optional for this dataset.
    cols = [np.ones_like(F), F, np.sin(2 * math.pi * M * F), np.cos(2 * math.pi * M * F)]
    if trend_on:
        t0 = t - t[0]
        cols += [np.cos(2 * math.pi * t), np.sin(2 * math.pi * t),
                 np.cos(4 * math.pi * t), np.sin(4 * math.pi * t), t0, t0 ** 2]
    A = np.column_stack(cols)
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ coef
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0


def discover_windings(t, F, y, n_top=3, seed=0) -> list[dict]:
    rng = np.random.default_rng(seed)
    real_r2 = np.array([regression_r2(t, F, y, M) for M in M_GRID])
    surr_r2 = np.zeros_like(M_GRID)
    for s in range(N_AR1_SURROGATES):
        ys = ar1_surrogate(y, rng)
        surr_r2 += np.array([regression_r2(t, F, ys, M) for M in M_GRID])
    surr_r2 /= N_AR1_SURROGATES
    excess = real_r2 - surr_r2

    # pick n_top local maxima of `excess`, at least 0.15 apart in M so we
    # don't just report the same ridge's shoulder three times
    order = np.argsort(-excess)
    picked = []
    for idx in order:
        M = M_GRID[idx]
        if all(abs(M - p["M"]) > 0.15 for p in picked):
            picked.append(dict(M=float(M), real_r2=float(real_r2[idx]),
                                surrogate_r2=float(surr_r2[idx]),
                                excess_r2=float(excess[idx])))
        if len(picked) >= n_top:
            break
    return picked


def main() -> int:
    if not NC_FILE.exists():
        print(f"missing {NC_FILE} -- fetch it first", file=sys.stderr)
        return 1
    donor = load_donor_params()
    print(f"[donor] borrowing shared manifold shape from {DONOR_INDEX}")

    print("[extract] chopping Kaplan SST into regions...")
    regions = extract_regions()
    print(f"[extract] {len(regions)} regions pass the >= {MIN_OCEAN_FRAC:.0%} "
          f"ocean-coverage / >=200-month threshold")

    for r in regions:
        idx_dir = K_SST_DIR / r["name"]
        idx_dir.mkdir(parents=True, exist_ok=True)
        with open(idx_dir / f"{r['name']}.dat", "w") as f:
            for t, v in zip(r["dates"], r["values"]):
                f.write(f"{t:.6f}\t{v:.6f}\n")
    print(f"[write] {len(regions)} <region>.dat files written under {K_SST_DIR}")

    results = []
    yl = ws.year_length(0.0, float(donor["year"]))
    for i, r in enumerate(regions):
        idx = r["name"]
        idx_dir = K_SST_DIR / idx
        # F(t) uses the SHARED backbone as Bessel's k1 -- fixed, not yet
        # dependent on this region's own (undiscovered) windings; see
        # load_donor_params for why this must be fixed before discovery.
        fake_prep = dict(params={**donor, "ltep": [donor["backbone_k1"]]},
                          yl=yl, dates=r["dates"], idate=r["dates"][0], nm=1,
                          year_startup=0.0, year_cand=float(donor["year"]))
        F = ws.build_forcing_at(fake_prep, r["dates"])
        picks = discover_windings(r["dates"], F, r["values"], n_top=3, seed=i)
        top_M = [p["M"] for p in picks]

        params = {k: v for k, v in donor.items() if k != "backbone_k1"}
        # backbone appended LAST so lt[nm-1] (Bessel's own k1 at fit time,
        # via ws.load_index) matches exactly what built F above
        params["ltep"] = top_M + [donor["backbone_k1"]]
        params["IR"] = 0.0
        (idx_dir / "lt.exe.p").write_text(json.dumps(params, indent=1))

        prep = ws.load_index(K_SST_DIR, idx)
        g = ws.fit_global(prep)
        results.append(dict(name=idx, lat=r["lat_c"], lon=r["lon_c"],
                             n_months=len(r["dates"]), top_M=top_M,
                             excess_r2=[p["excess_r2"] for p in picks],
                             train_r=g["train_r"], val_r=g["val_r"]))
        print(f"  [{i+1:3d}/{len(regions)}] {idx:12s} lat={r['lat_c']:+6.1f} "
              f"lon={r['lon_c']:+6.1f}  M={np.round(top_M,3).tolist()}  "
              f"train_r={g['train_r']:+.3f}  val_r={g['val_r']:+.3f}")

    out_json = K_SST_DIR / "k_sst_summary.json"
    with open(out_json, "w") as f:
        json.dump(results, f, indent=1)
    print(f"\n[written] {out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
