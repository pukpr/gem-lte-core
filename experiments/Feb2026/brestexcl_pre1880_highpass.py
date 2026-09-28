#!/usr/bin/env python3
"""brestexcl_pre1880_highpass.py — the corrected pre-1880 CV test, using
a HIGH-PASS reconstruction (short boxcar smooth, then subtract a
several-YEAR-wide moving average) instead of the earlier boxcar+
quadratic-detrend guess at massage.py's algorithm.

Motivated directly by a visual observation: the earlier massage.py
reconstruction (13pt boxcar + GLOBAL quadratic detrend) matched
brestexcl.dat's TIMING reasonably well but was systematically damped in
AMPLITUDE (std ratio ~0.62-0.78 depending on window) -- a global
2-parameter quadratic can only remove one parabola shape over the whole
144-year record, so it can't track real DECADAL wander the way a
rolling several-year-wide average can. Grid search over (short smooth
width, long subtract-window width) confirms this: best found is an
~11-month short smooth followed by subtracting its own 5.5-year rolling
mean, giving r=+0.704 against brestexcl.dat over 1880-2020 -- a real,
substantial improvement over the earlier r=+0.552 (quadratic-detrend)
reconstruction, directly validating the high-pass-with-a-multi-year-
window structure.
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
SHORT_W = 11     # months, the short smoothing width
LONG_YEARS = 5.5  # years, the high-pass subtraction window


def decompensate_f9(x, ir, f9):
    out = x.copy()
    if ir != 0.0:
        for i in range(len(out) - 1, 11, -1):
            out[i] -= ir * out[i - 12]
    for _ in range(f9):
        out = filter9point(out)
    return out


def boxcar(values, window):
    n = len(values)
    out = np.zeros(n)
    half = window // 2
    for i in range(n):
        s, e = max(0, i - half), min(n, i + half + 1)
        out[i] = values[s:e].mean()
    return out


def main() -> int:
    times, values, full_idx = [], [], []
    with open(ROOT / IDX / "psmsl_brest_1.rlrdata") as f:
        for i, line in enumerate(f):
            t, v = float(line.split(";")[0]), float(line.split(";")[1])
            if v == -99999:
                continue
            times.append(t)
            values.append(v)
            full_idx.append(i)
    times, values, full_idx = np.array(times), np.array(values), np.array(full_idx)
    proj_dates = 1807.0 + full_idx / 12.0

    smoothed = boxcar(values, SHORT_W)
    long_window_months = int(round(LONG_YEARS * 12))
    long_ma = boxcar(smoothed, long_window_months)
    psmsl_hp = smoothed - long_ma

    prep = prepare(IDX)
    dates_orig, data_orig, forcing_orig, m = (
        prep["dates"], prep["data_raw"], prep["forcing"], prep["m"])
    nz_orig = data_orig != 0.0

    recon_by_month = {int(round(d * 12)): v for d, v in zip(proj_dates, psmsl_hp)}
    d_check, b_check, r_check_arr = [], [], []
    for do, bv in zip(dates_orig[nz_orig], data_orig[nz_orig]):
        key = int(round(do * 12))
        if key in recon_by_month:
            d_check.append(do)
            b_check.append(bv)
            r_check_arr.append(recon_by_month[key])
    d_check, b_check, r_check_arr = map(np.array, (d_check, b_check, r_check_arr))
    r_check = float(np.corrcoef(r_check_arr, b_check)[0, 1])
    print(f"[check] high-pass PSMSL ({SHORT_W}mo smooth, {LONG_YEARS}yr subtract) "
          f"vs brestexcl.dat, 1880-2020: r={r_check:+.4f}  "
          f"std ratio={r_check_arr.std()/b_check.std():.3f}  "
          f"(previous quadratic-detrend reconstruction: r=+0.552)")

    # --- Step-1 validated model, unchanged --------------------------------
    level, k0, amp, phase, annual, trend, accel = regression_factors(
        dates_orig[nz_orig], forcing_orig[nz_orig], data_orig[nz_orig], m,
        prep["trend_on"])
    params = json.loads((ROOT / IDX / "lt.exe.p").read_text())
    resp = read_resp(ROOT / IDX / "lt.exe.resp")
    cfg = Config(resp, {})
    yl = year_length(cfg.get("YEAR", 0.0), float(params.get("year", 0.0)))
    forcing_ext = build_forcing_continuous(IDX, yl, proj_dates)
    model_raw = lte_response(forcing_ext, proj_dates, m, amp, phase, level,
                              k0, trend, accel, prep["nonlin"], annual)
    model = decompensate_f9(model_raw, prep["ir"], prep["f9"])

    model_by_month = {int(round(d * 12)): v for d, v in zip(proj_dates, model)}
    mv = np.array([model_by_month[int(round(do * 12))] for do in d_check])
    in_sample_r = float(np.corrcoef(mv, r_check_arr)[0, 1])
    print(f"[{IDX}] Step-1 model vs high-pass-processed PSMSL, "
          f"1880-2020 overlap: r={in_sample_r:+.4f}")

    # --- the actual test: 1840-1880 -----------------------------------
    test_mask = (proj_dates >= 1840.0) & (proj_dates < dates_orig[0])
    d_test = proj_dates[test_mask]
    hp_test = psmsl_hp[test_mask]
    model_test = np.array([model_by_month.get(int(round(dt * 12)), np.nan) for dt in d_test])
    ok_test = ~np.isnan(model_test)
    r_test = float(np.corrcoef(model_test[ok_test], hp_test[ok_test])[0, 1])
    print(f"[RESULT] 1840-1880 (never-seen, n={ok_test.sum()}): r={r_test:+.4f}")

    rng = np.random.default_rng(42)
    n = len(model)
    null_rs = []
    for _ in range(500):
        shift = rng.integers(200, n - 200)
        shifted = np.roll(model, shift)
        shifted_by_month = {int(round(d * 12)): v for d, v in zip(proj_dates, shifted)}
        sv = np.array([shifted_by_month.get(int(round(dt * 12)), np.nan) for dt in d_test])
        okk = ~np.isnan(sv)
        null_rs.append(float(np.corrcoef(sv[okk], hp_test[okk])[0, 1]))
    null_rs = np.array(null_rs)
    z = (r_test - null_rs.mean()) / null_rs.std()
    print(f"[null] 500 random circular-shift nulls: mean={null_rs.mean():+.4f} "
          f"sd={null_rs.std():.4f}  ->  z={z:+.2f}")

    fig, axes = plt.subplots(3, 1, figsize=(13, 10))
    axes[0].plot(d_check, b_check, color="0.25", linewidth=0.8, label="real brestexcl.dat")
    axes[0].plot(d_check, r_check_arr, color="teal", linewidth=0.8, alpha=0.8,
                 label=f"high-pass reconstruction ({SHORT_W}mo smooth - {LONG_YEARS}yr MA)")
    axes[0].set_title(f"full overlap, 1880-2020 (r={r_check:+.3f})", fontsize=10, fontweight="bold")
    axes[0].legend(fontsize=8, loc="upper left")

    zmask = (d_check >= 2000) & (d_check < 2015)
    axes[1].plot(d_check[zmask], b_check[zmask], color="0.25", linewidth=1.2, marker="o",
                 markersize=3, label="real brestexcl.dat")
    axes[1].plot(d_check[zmask], r_check_arr[zmask], color="teal", linewidth=1.2, marker="o",
                 markersize=3, alpha=0.8, label="high-pass reconstruction")
    axes[1].set_title("zoom: 2000-2015", fontsize=10)
    axes[1].legend(fontsize=8, loc="upper left")

    show = (proj_dates >= 1835) & (proj_dates < 1900)
    axes[2].plot(proj_dates[show], psmsl_hp[show], color="0.25", linewidth=1.0,
                 label="real PSMSL Brest, high-pass-processed (never used to fit anything)")
    axes[2].plot(proj_dates[show], model[show], color="teal", linewidth=1.0, alpha=0.85,
                 label="Step-1 brestexcl model, extrapolated backward")
    axes[2].axvspan(1840, 1880, color="crimson", alpha=0.08, label="1840-1880: the test window")
    axes[2].axvspan(1835, 1840, color="0.85", label="excluded (pre-1840 step-change era)")
    axes[2].set_title(f"1840-1880 test (r={r_test:+.3f}, z={z:+.1f} vs null)",
                       fontsize=10, fontweight="bold")
    axes[2].legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    out = ROOT / IDX / "brestexcl_pre1880_highpass.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"[saved] {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
