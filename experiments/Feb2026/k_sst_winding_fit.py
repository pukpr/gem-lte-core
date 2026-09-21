#!/usr/bin/env python3
"""k_sst_winding_fit.py -- the "variant of ws.py" for the Kaplan SST
regional sweep, per the plan: fit backbone + a SMALL capped set of
extra windings (2-3 core, up to 5 stretch), train on the later part of
each region's own record, validate on an EARLIER held-out stretch,
DETREND before fitting/scoring (so R isn't inflated by matching a trend
rather than genuine tidal structure), and skip any accel/quadratic term
entirely (a real parabolic trend can be mistaken for the long comb-alias
"120yr" cycle -- see AMO_FORMULATION.md -- so the safest thing is to not
let the two compete for the same variance in the first place).

DATA-AVAILABILITY CONSTRAINT (checked before running anything): the bulk
k_sst/ regions built by build_k_sst.py are restricted to 1950-2023 only
(its own MIN_YEAR=1950, documented reason: pre-1950 Kaplan coverage is
sparse/ship-track-only and not trustworthy at 20-degree grid resolution).
Only the 8 top-level named regions already carried through a real Ada
optimization this session (kN_baltic, kN_northsea, kN_brest, kN040_E030,
kN040_E130, kS040_E130, kN_indostraits, kN_okhotsk) have real data back
to 1856. So the literal ">1880 train / <1880 validate" split from the
plan is used ONLY for those 8; the ~80+ bulk 1950-2023 regions use an
analogous held-out-earlier-stretch split within their own real range
(train >=1970, validate 1950-1970) rather than fabricate pre-1880 data
or discard the documented data-quality caveat.

Discovery reuses winding_rank.rank_series (the same AR1-floor + FWHM +
continuity gated test build_k_sst_v2.py already validated, including its
record-length-proportional sigma/t0_step scaling) -- run on the
TRAINING portion only, so no validation-period information leaks into
which windings get selected.

Rollout: starts with a PILOT (the 8 known-good named regions, plus the
top bulk regions by build_k_sst_v2.py's own existing val_r) rather than
the full 80+ sweep, per "start on regions that have decent fits and
radially split out from there".
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
import winding_rank as wr  # noqa: E402

K_SST_DIR = bk.K_SST_DIR
TOP_LEVEL_1856 = ["kN_baltic", "kN_northsea", "kN_brest", "kN040_E030",
                  "kN040_E130", "kS040_E130", "kN_indostraits", "kN_okhotsk"]


def linear_detrend(t: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    A = np.column_stack([np.ones_like(t), t - t[0]])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return y - A @ coef, coef


def fit_at_Ms(F: np.ndarray, y: np.ndarray, Ms: list[float]) -> tuple[np.ndarray, float]:
    cols = [np.ones_like(F)]
    for m in Ms:
        cols.append(np.sin(2 * np.pi * m * F))
        cols.append(np.cos(2 * np.pi * m * F))
    A = np.column_stack(cols)
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    fitted = A @ coef
    r = float(np.corrcoef(fitted, y)[0, 1]) if y.std() > 0 else 0.0
    return coef, r


def run_region(name: str, dates: np.ndarray, values: np.ndarray,
               donor: dict, split_year: float, max_core: int = 3,
               max_stretch: int = 5) -> dict:
    yl = ws.year_length(0.0, float(donor["year"]))
    fake_prep = dict(params={**donor, "ltep": [donor["backbone_k1"]]},
                      yl=yl, dates=dates, idate=dates[0], nm=1,
                      year_startup=0.0, year_cand=float(donor["year"]))
    F = ws.build_forcing_at(fake_prep, dates)
    backbone = donor["backbone_k1"]

    nz = values != 0.0
    train_idx = np.nonzero((dates >= split_year) & nz)[0]
    val_idx = np.nonzero((dates < split_year) & nz)[0]
    if len(train_idx) < 50 or len(val_idx) < 20:
        return dict(name=name, status="skipped (too little data on one side)")

    d_train, F_train = linear_detrend(dates[train_idx], values[train_idx])[0], F[train_idx]
    d_val = values[val_idx] - np.polyval(
        np.polyfit(dates[train_idx], values[train_idx], 1), dates[val_idx])
    F_val = F[val_idx]

    span_years = float(dates[train_idx][-1] - dates[train_idx][0])
    sigma_scaled = max(span_years * (15.0 / 167.0), 3.0)
    t0_step_scaled = max(span_years * (5.0 / 167.0), 1.0)
    try:
        ridges = wr.rank_series(dates[train_idx], d_train, F_train,
                                sigma=sigma_scaled, t0_step=t0_step_scaled)
    except SystemExit:
        ridges = []
    passing = [x for x in ridges if x["pass"] and abs(x["M"]) > 0.01
               and abs(x["M"] - backbone) > 0.03]
    passing.sort(key=lambda x: -x["peak_bits"])

    result = dict(name=name, n_candidates=len(passing),
                  candidates=[x["M"] for x in passing])
    for label, k in [("core", max_core), ("stretch", max_stretch)]:
        Ms = [backbone] + [x["M"] for x in passing[:k]]
        coef, train_r = fit_at_Ms(F_train, d_train, Ms)
        val_pred = np.column_stack(
            [np.ones_like(F_val)] +
            sum(([np.sin(2 * np.pi * m * F_val), np.cos(2 * np.pi * m * F_val)] for m in Ms), [])
        ) @ coef
        val_r = float(np.corrcoef(val_pred, d_val)[0, 1]) if d_val.std() > 0 else 0.0
        result[f"{label}_n_windings"] = len(Ms)
        result[f"{label}_M"] = Ms
        result[f"{label}_train_r"] = train_r
        result[f"{label}_val_r"] = val_r
    return result


def main() -> None:
    donor = bk.load_donor_params()
    print(f"[k_sst_winding_fit] donor backbone_k1={donor['backbone_k1']:.6f}")

    pilot = []
    for name in TOP_LEVEL_1856:
        dat = HERE / name / f"{name}.dat"
        if dat.exists():
            pilot.append((name, dat, 1880.0))

    summary_path = K_SST_DIR / "k_sst_summary_v2.json"
    if summary_path.exists():
        summ = json.loads(summary_path.read_text())
        summ = [s for s in summ if np.isfinite(s.get("val_r", float("nan")))]
        summ.sort(key=lambda s: -s["val_r"])
        for s in summ[:15]:
            dat = K_SST_DIR / s["name"] / f"{s['name']}.dat"
            if dat.exists():
                pilot.append((s["name"], dat, 1970.0))

    print(f"[k_sst_winding_fit] pilot: {len(pilot)} regions "
          f"({len(TOP_LEVEL_1856)} known-good 1856+ named, "
          f"{len(pilot)-len(TOP_LEVEL_1856)} top bulk-region by existing val_r)\n")

    results = []
    t0 = time.time()
    for i, (name, dat, split_year) in enumerate(pilot):
        raw = np.loadtxt(dat)
        dates, values = raw[:, 0], raw[:, 1]
        r = run_region(name, dates, values, donor, split_year)
        results.append(r)
        if "status" in r:
            print(f"  [{i+1:2d}/{len(pilot)}] {name:14s} {r['status']}")
            continue
        print(f"  [{i+1:2d}/{len(pilot)}] {name:14s} split={split_year:.0f}  "
              f"n_cand={r['n_candidates']:2d}  "
              f"core(n={r['core_n_windings']}): train_r={r['core_train_r']:+.3f} val_r={r['core_val_r']:+.3f}  "
              f"stretch(n={r['stretch_n_windings']}): train_r={r['stretch_train_r']:+.3f} val_r={r['stretch_val_r']:+.3f}")

    ok = [r for r in results if "status" not in r]
    core_val = np.array([r["core_val_r"] for r in ok])
    stretch_val = np.array([r["stretch_val_r"] for r in ok])
    core_train = np.array([r["core_train_r"] for r in ok])
    stretch_train = np.array([r["stretch_train_r"] for r in ok])
    print(f"\n[k_sst_winding_fit] {len(ok)}/{len(pilot)} regions fit")
    print(f"  CORE (backbone+<=3):    train_r mean={core_train.mean():+.3f}  "
          f">=0.70: {np.sum(core_train>=0.70)}/{len(ok)}   |  val_r mean={core_val.mean():+.3f}")
    print(f"  STRETCH (backbone+<=5): train_r mean={stretch_train.mean():+.3f}  "
          f">=0.75: {np.sum(stretch_train>=0.75)}/{len(ok)}   |  val_r mean={stretch_val.mean():+.3f}")

    out = HERE / "k_sst_winding_fit_pilot.json"
    out.write_text(json.dumps(results, indent=1))
    print(f"\n[k_sst_winding_fit] saved {out}  ({(time.time()-t0)/60:.1f} min)")


if __name__ == "__main__":
    main()
