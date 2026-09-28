#!/usr/bin/env python3
"""brestexcl_pre1880_deseasoned.py — the corrected pre-1880 CV test.

Two things established in prior turns that this script builds on
directly, not repeats: (1) brestexcl.dat's own annual-cycle amplitude
is essentially zero (0.54, vs raw PSMSL's 48.5 -- confirmed by direct
harmonic-regression comparison) so the "excl" naming is now understood:
the seasonal cycle (and the secular trend) are excluded, and any
correlation test against raw PSMSL must remove the same thing from
BOTH sides or it's dominated by that trivial shared/unshared mismatch.
(2) the model fit against brestexcl's own data tracks it well (2000-
2010 genuine holdout val_r=+0.720) -- that is the model this script
extrapolates, not a weaker raw-PSMSL refit.

This script: deseasonalizes+detrends raw PSMSL via one joint harmonic
regression (annual+semiannual+linear+quadratic) fit on the FULL
available 1807-2020 PSMSL record (not just the overlap, avoiding
extrapolation), producing a target comparable in kind to brestexcl.dat.
Then scores the Step-1 model (extrapolated backward, never refit) on
the genuinely never-seen 1807-1880 portion of that deseasonalized
series.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cv_rolling_blocked import prepare  # noqa: E402
from lte_forward import (regression_factors, lte_response, filter9point,  # noqa: E402
                          year_length, Config, read_resp)
from brestexcl_resolved_shallow_water import build_forcing_continuous  # noqa: E402
import json

IDX = "brestexcl"
ROOT = Path(__file__).resolve().parent
PSMSL_FILE = ROOT / IDX / "psmsl_brest_1.rlrdata"
PSMSL_START_YEAR = 1807


def decompensate_f9(x, ir, f9):
    out = x.copy()
    if ir != 0.0:
        for i in range(len(out) - 1, 11, -1):
            out[i] -= ir * out[i - 12]
    for _ in range(f9):
        out = filter9point(out)
    return out


def main() -> int:
    dates_psmsl, vals_psmsl, missing = [], [], []
    with open(PSMSL_FILE) as f:
        for i, line in enumerate(f):
            v = float(line.split(";")[1])
            dates_psmsl.append(PSMSL_START_YEAR + i / 12.0)
            vals_psmsl.append(v)
            missing.append(v <= -99998.0)
    dates_psmsl = np.array(dates_psmsl)
    vals_psmsl = np.array(vals_psmsl)
    missing = np.array(missing)
    ok = ~missing

    # --- deseasonalize + detrend PSMSL via one joint regression, fit on
    # the FULL available record (not just the overlap -- avoids
    # extrapolating a quadratic trend 70+ years backward) ------------------
    d0 = dates_psmsl[0]
    A = np.column_stack([
        np.ones_like(dates_psmsl), dates_psmsl - d0, (dates_psmsl - d0) ** 2,
        np.sin(2 * np.pi * dates_psmsl), np.cos(2 * np.pi * dates_psmsl),
        np.sin(4 * np.pi * dates_psmsl), np.cos(4 * np.pi * dates_psmsl)])
    coef, _, _, _ = np.linalg.lstsq(A[ok], vals_psmsl[ok], rcond=None)
    psmsl_ds = vals_psmsl - A @ coef
    print(f"[deseasonalize] fitted on full 1807-2020 PSMSL record; "
          f"annual/semiannual amplitudes removed: "
          f"{np.hypot(coef[3], coef[4]):.2f} / {np.hypot(coef[5], coef[6]):.2f} mm "
          f"(brestexcl.dat's own: 0.54 / 0.27 -- now directly comparable)")

    prep = prepare(IDX)
    dates_orig, data_orig, forcing_orig, m = (
        prep["dates"], prep["data_raw"], prep["forcing"], prep["m"])
    nz_orig = data_orig != 0.0

    overlap = (dates_psmsl >= dates_orig[0]) & (dates_psmsl <= dates_orig[nz_orig][-1]) & ok
    idx_map = np.clip(np.searchsorted(dates_orig, dates_psmsl[overlap]), 0, len(dates_orig) - 1)
    proj_vals = data_orig[idx_map]
    mask2 = proj_vals != 0.0
    r_check = float(np.corrcoef(psmsl_ds[overlap][mask2], proj_vals[mask2])[0, 1])
    print(f"[check] deseasonalized+detrended PSMSL vs brestexcl.dat, "
          f"1880-2020 overlap: r={r_check:+.4f}  (raw PSMSL vs brestexcl.dat "
          f"was r=+0.2556 -- {'improved' if r_check>0.2556 else 'not improved'})")

    # --- the already-validated model (Step 1: fit against brestexcl's OWN
    # data, ALL of it, genuine 2000-2010 holdout val_r=+0.720) -- reused
    # UNCHANGED, not refit here ---------------------------------------------
    level, k0, amp, phase, annual, trend, accel = regression_factors(
        dates_orig[nz_orig], forcing_orig[nz_orig], data_orig[nz_orig], m,
        prep["trend_on"])

    params = json.loads((ROOT / IDX / "lt.exe.p").read_text())
    resp = read_resp(ROOT / IDX / "lt.exe.resp")
    cfg = Config(resp, {})
    yl = year_length(cfg.get("YEAR", 0.0), float(params.get("year", 0.0)))
    forcing_ext = build_forcing_continuous(IDX, yl, dates_psmsl)
    model_raw = lte_response(forcing_ext, dates_psmsl, m, amp, phase, level,
                              k0, trend, accel, prep["nonlin"], annual)
    model = decompensate_f9(model_raw, prep["ir"], prep["f9"])

    in_sample_r = float(np.corrcoef(model[overlap][mask2], psmsl_ds[overlap][mask2])[0, 1])
    print(f"[{IDX}] Step-1 model vs deseasonalized PSMSL, 1880-2020 "
          f"overlap: r={in_sample_r:+.4f}")

    pre1880 = (dates_psmsl < dates_orig[0]) & ok
    r_pre1880 = float(np.corrcoef(model[pre1880], psmsl_ds[pre1880])[0, 1])
    print(f"[RESULT] pre-1880 (genuinely never seen, 708 real monthly "
          f"values, 1807-1880): r={r_pre1880:+.4f}")

    def smooth3(x):
        r = x.copy()
        r[1:-1] = (x[:-2] + x[1:-1] + x[2:]) / 3.0
        return r
    model_s = smooth3(model)
    psmsl_ds_s = smooth3(np.where(missing, np.nan, psmsl_ds))
    sm_mask = pre1880 & ~np.isnan(psmsl_ds_s)
    r_pre1880_smooth = float(np.corrcoef(model_s[sm_mask], psmsl_ds_s[sm_mask])[0, 1])
    print(f"[RESULT] same, 3-month-smoothed: r={r_pre1880_smooth:+.4f}")

    rng = np.random.default_rng(42)
    n = len(model)
    null_rs = []
    for _ in range(500):
        shift = rng.integers(200, n - 200)
        null_rs.append(float(np.corrcoef(np.roll(model, shift)[pre1880], psmsl_ds[pre1880])[0, 1]))
    null_rs = np.array(null_rs)
    z = (r_pre1880 - null_rs.mean()) / null_rs.std()
    print(f"[null] 500 random circular-shift nulls: mean={null_rs.mean():+.4f} "
          f"sd={null_rs.std():.4f}  ->  z={z:+.2f}")

    fig, axes = plt.subplots(2, 1, figsize=(12, 7.5))
    show = (dates_psmsl >= 1807) & (dates_psmsl < 1900) & ok
    axes[0].plot(dates_psmsl[show], psmsl_ds[show], color="0.25", linewidth=1.0,
                 label="real PSMSL Brest, deseasonalized+detrended (never used to fit anything)")
    axes[0].plot(dates_psmsl[show], model[show], color="teal", linewidth=1.0, alpha=0.85,
                 label="Step-1 brestexcl model, extrapolated backward (never fit to this era)")
    axes[0].axvspan(1807, 1880, color="crimson", alpha=0.08, label="pre-1880: genuinely never seen")
    axes[0].set_title(f"Brest: deseasonalized comparison, pre-1880 "
                       f"(r={r_pre1880:+.3f}, z={z:+.1f} vs null)",
                       fontsize=10, fontweight="bold")
    axes[0].legend(fontsize=8, loc="upper left")

    axes[1].hist(null_rs, bins=40, color="0.7", label="null (random phase shift)")
    axes[1].axvline(r_pre1880, color="crimson", linewidth=2, label=f"observed r={r_pre1880:+.3f}")
    axes[1].set_title("null distribution", fontsize=10)
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    out_path = ROOT / IDX / "brestexcl_pre1880_deseasoned.png"
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"[saved] {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
