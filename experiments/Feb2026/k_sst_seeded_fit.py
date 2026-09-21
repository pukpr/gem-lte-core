#!/usr/bin/env python3
"""k_sst_seeded_fit.py -- seed each Kaplan SST grid cell's winding set
from its nearest already-validated named index (AMO / PDO / NINO4)
instead of blind per-cell discovery, which the previous pilot
(k_sst_winding_fit.py) showed has almost no statistical power once the
trend is properly removed at single-20-degree-cell resolution.

Framing: this is an information-theory exercise, not a from-scratch
climate-science discovery exercise. AMO and NINO4 are each, on their
own, complex enough to have decades of dedicated literature; if a SHARED
external manifold plus a HANDFUL of free numbers (the same few numbers,
borrowed wholesale from the nearest named index) is enough to get a
nominal, honestly-validated fit at nearby grid cells, that is itself the
finding worth reporting -- an extreme compression of "how much
information does this location's SST variability actually contain
beyond the shared tidal clock" to a handful of borrowed parameters.
Scope is deliberately narrow: nominal fits, post-1950 (Kaplan's
trustworthy window), good enough to test the kriging/spatial-coherence
idea -- not a rigorous per-cell discovery exercise.

Method per grid cell:
  1. Find the nearest of {amo, pdo, nino4} by great-circle distance.
  2. Take that donor's own already-fitted top windings (by amplitude,
     excluding the shared backbone) -- 2 for "core", 4 for "stretch" --
     AS FIXED CANDIDATE FREQUENCIES. No per-cell discovery/gating at all.
  3. Detrend the cell's own data (linear, fit on the training portion).
  4. Fit backbone + seeded windings via one direct OLS against the
     cell's own DETRENDED data on the training portion (>=1970); score
     R on both train and the held-out 1950-1970 validation portion.

This is a much cheaper, much more powerful test than blind discovery:
fitting a handful of PRE-SPECIFIED frequencies is an ordinary regression,
not a multiple-comparisons search across a whole M-grid, so it doesn't
need winding_rank's AR1-floor Monte Carlo at all -- cheap enough to run
across the whole grid, not just a pilot.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ws  # noqa: E402
import build_k_sst as bk  # noqa: E402
from geospatial_winding_analysis import haversine_km  # noqa: E402

DONOR_COORDS = {"amo": (40.0, -40.0), "pdo": (40.0, -180.0), "nino4": (0.0, -175.0)}
BACKBONE_TOL = 0.03


def donor_seed_Ms(idx: str, n: int) -> list[float]:
    d = json.loads((HERE / idx / "lt.exe.windings.json").read_text())
    backbone = json.loads((HERE / idx / "lt.exe.p").read_text())["ltep"][-1]
    extras = [(m, a) for m, a, p in d["k_amp_phase"] if abs(m - backbone) > BACKBONE_TOL]
    extras.sort(key=lambda x: -x[1])
    return [m for m, a in extras[:n]]


def nearest_donor(lat: float, lon: float) -> str:
    return min(DONOR_COORDS, key=lambda d: haversine_km(lat, lon, *DONOR_COORDS[d]))


def linear_detrend_fit(t, y):
    A = np.column_stack([np.ones_like(t), t - t[0]])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return coef


def fit_and_score(F_train, d_train, F_val, d_val, Ms):
    cols_tr = [np.ones_like(F_train)]
    for m in Ms:
        cols_tr += [np.sin(2 * np.pi * m * F_train), np.cos(2 * np.pi * m * F_train)]
    A_tr = np.column_stack(cols_tr)
    coef, *_ = np.linalg.lstsq(A_tr, d_train, rcond=None)
    train_r = float(np.corrcoef(A_tr @ coef, d_train)[0, 1]) if d_train.std() > 0 else 0.0

    cols_val = [np.ones_like(F_val)]
    for m in Ms:
        cols_val += [np.sin(2 * np.pi * m * F_val), np.cos(2 * np.pi * m * F_val)]
    A_val = np.column_stack(cols_val)
    val_pred = A_val @ coef
    val_r = float(np.corrcoef(val_pred, d_val)[0, 1]) if d_val.std() > 0 else 0.0
    return train_r, val_r


def main() -> None:
    donor = bk.load_donor_params()
    backbone = donor["backbone_k1"]
    seed_core = {i: donor_seed_Ms(i, 2) for i in DONOR_COORDS}
    seed_stretch = {i: donor_seed_Ms(i, 4) for i in DONOR_COORDS}
    print("[k_sst_seeded_fit] donor seed windings:")
    for i in DONOR_COORDS:
        print(f"  {i}: core={[f'{m:+.3f}' for m in seed_core[i]]}  "
              f"stretch={[f'{m:+.3f}' for m in seed_stretch[i]]}")

    yl = ws.year_length(0.0, float(donor["year"]))
    regions = bk.extract_regions()
    print(f"\n[k_sst_seeded_fit] {len(regions)} grid cells\n")

    results = []
    t0 = time.time()
    for i, r in enumerate(regions):
        name, lat, lon = r["name"], r["lat_c"], r["lon_c"]
        dates, values = r["dates"], r["values"]
        nz = values != 0.0
        donor_name = nearest_donor(lat, lon)

        fake_prep = dict(params={**donor, "ltep": [backbone]}, yl=yl, dates=dates,
                          idate=dates[0], nm=1, year_startup=0.0, year_cand=float(donor["year"]))
        F = ws.build_forcing_at(fake_prep, dates)

        split_year = 1970.0
        train_idx = np.nonzero((dates >= split_year) & nz)[0]
        val_idx = np.nonzero((dates < split_year) & nz)[0]
        if len(train_idx) < 50 or len(val_idx) < 20:
            results.append(dict(name=name, lat=lat, lon=lon, status="skipped (too little data)"))
            continue

        trend_coef = linear_detrend_fit(dates[train_idx], values[train_idx])
        A_tr_trend = np.column_stack([np.ones_like(dates[train_idx]), dates[train_idx] - dates[train_idx][0]])
        d_train = values[train_idx] - A_tr_trend @ trend_coef
        A_val_trend = np.column_stack([np.ones_like(dates[val_idx]), dates[val_idx] - dates[train_idx][0]])
        d_val = values[val_idx] - A_val_trend @ trend_coef

        F_train, F_val = F[train_idx], F[val_idx]
        row = dict(name=name, lat=lat, lon=lon, donor=donor_name)
        for label, seeds in [("core", seed_core[donor_name]), ("stretch", seed_stretch[donor_name])]:
            Ms = [backbone] + seeds
            train_r, val_r = fit_and_score(F_train, d_train, F_val, d_val, Ms)
            row[f"{label}_n_windings"] = len(Ms)
            row[f"{label}_train_r"] = train_r
            row[f"{label}_val_r"] = val_r
        results.append(row)
        print(f"  [{i+1:3d}/{len(regions)}] {name:12s} lat={lat:+5.1f} lon={lon:+6.1f}  "
              f"donor={donor_name:6s}  core: train_r={row['core_train_r']:+.3f} val_r={row['core_val_r']:+.3f}  "
              f"stretch: train_r={row['stretch_train_r']:+.3f} val_r={row['stretch_val_r']:+.3f}")

    ok = [r for r in results if "status" not in r]
    for label in ["core", "stretch"]:
        tr = np.array([r[f"{label}_train_r"] for r in ok])
        va = np.array([r[f"{label}_val_r"] for r in ok])
        print(f"\n[{label}] n={len(ok)}  train_r: mean={tr.mean():+.3f} median={np.median(tr):+.3f} "
              f">=0.5: {np.sum(tr>=0.5)}  >=0.7: {np.sum(tr>=0.7)}   |  "
              f"val_r: mean={va.mean():+.3f} median={np.median(va):+.3f} >=0.5: {np.sum(va>=0.5)}")

    for dn in DONOR_COORDS:
        sub = [r for r in ok if r["donor"] == dn]
        if sub:
            tr = np.array([r["core_train_r"] for r in sub])
            print(f"  by donor={dn}: n={len(sub)}  core train_r mean={tr.mean():+.3f} median={np.median(tr):+.3f}")

    out = HERE / "k_sst_seeded_fit.json"
    out.write_text(json.dumps(results, indent=1))
    print(f"\n[k_sst_seeded_fit] saved {out}  ({(time.time()-t0)/60:.1f} min)")


if __name__ == "__main__":
    main()
