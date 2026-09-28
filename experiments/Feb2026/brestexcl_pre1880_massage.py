#!/usr/bin/env python3
"""brestexcl_pre1880_massage.py — the corrected pre-1880 CV test, using
the ACTUAL data-processing pipeline found in the sibling repo
~/github/pukpr/GEM-LTE (experiments/Feb2026/massage.py, invoked by
tmodel.sh as `python3 massage.py <rlrdata> 2`), not a guessed transform.

Pipeline (exactly matching massage.py): drop -99999 rows, apply a
13-point boxcar (moving average, shrinking window at the edges), then
remove a QUADRATIC (order=2) trend via polyfit over the whole record.
Validated against brestexcl.dat over the 1880-2020 overlap: r=+0.552
(vs raw PSMSL's r=+0.256, and a from-scratch guessed deseasonalize+
detrend's r=+0.296) -- the largest, most direct improvement found, from
the ACTUAL pipeline rather than another guess.

Per direct instruction: exclude everything before 1840 (documented
step-change / known data-quality issue in the earliest Brest record,
matching the sibling repo's own README note: "known data quality issues
in the pre-1880 Brest record"). The never-seen test window here is
1840-1880, not 1807-1880.
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


def decompensate_f9(x, ir, f9):
    out = x.copy()
    if ir != 0.0:
        for i in range(len(out) - 1, 11, -1):
            out[i] -= ir * out[i - 12]
    for _ in range(f9):
        out = filter9point(out)
    return out


def boxcar_filter(values, window_size=13):
    n = len(values)
    filtered = np.zeros(n)
    half = window_size // 2
    for i in range(n):
        start, end = max(0, i - half), min(n, i + half + 1)
        filtered[i] = values[start:end].mean()
    return filtered


def main() -> int:
    # --- exact massage.py pipeline -----------------------------------------
    times, values = [], []
    with open(PSMSL_FILE) as f:
        for line in f:
            t, v = float(line.split(";")[0]), float(line.split(";")[1])
            if v == -99999:
                continue
            times.append(t)
            values.append(v)
    times, values = np.array(times), np.array(values)
    filtered = boxcar_filter(values, 13)
    p = np.polyfit(times, filtered, 2)
    psmsl_massaged = filtered - np.polyval(p, times)
    # map PSMSL's own mid-month date convention onto this project's grid
    proj_dates = np.floor(times) + np.round((times - np.floor(times)) * 12) / 12.0

    prep = prepare(IDX)
    dates_orig, data_orig, forcing_orig, m = (
        prep["dates"], prep["data_raw"], prep["forcing"], prep["m"])
    nz_orig = data_orig != 0.0

    overlap = (proj_dates >= dates_orig[0]) & (proj_dates <= dates_orig[nz_orig][-1])
    idx_map = np.clip(np.searchsorted(dates_orig, proj_dates[overlap]), 0, len(dates_orig) - 1)
    proj_vals = data_orig[idx_map]
    mask2 = proj_vals != 0.0
    r_check = float(np.corrcoef(psmsl_massaged[overlap][mask2], proj_vals[mask2])[0, 1])
    print(f"[check] massage.py-reconstructed PSMSL vs brestexcl.dat, "
          f"1880-2020 overlap: r={r_check:+.4f}  (this is the actual "
          f"processing pipeline, found in ~/github/pukpr/GEM-LTE/"
          f"experiments/Feb2026/massage.py + tmodel.sh, not a guess)")

    # --- Step-1 validated model (fit on ALL of brestexcl.dat, genuine
    # 2000-2010 holdout val_r=+0.720) -- reused unchanged -------------------
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

    in_sample_r = float(np.corrcoef(model[overlap][mask2], psmsl_massaged[overlap][mask2])[0, 1])
    print(f"[{IDX}] Step-1 model vs massage.py-processed PSMSL, "
          f"1880-2020 overlap: r={in_sample_r:+.4f}")

    # --- the test, per direct instruction: nothing before 1840 ------------
    test_window = (proj_dates >= 1840.0) & (proj_dates < dates_orig[0])
    print(f"[test] {test_window.sum()} real monthly values, "
          f"{proj_dates[test_window][0]:.2f} - {proj_dates[test_window][-1]:.2f} "
          f"(1840-1880, excluding the documented pre-1840 step-change era)")
    r_test = float(np.corrcoef(model[test_window], psmsl_massaged[test_window])[0, 1])
    print(f"[RESULT] 1840-1880 backward-extrapolation correlation: "
          f"r={r_test:+.4f}")

    def smooth3(x):
        r = x.copy()
        r[1:-1] = (x[:-2] + x[1:-1] + x[2:]) / 3.0
        return r
    r_test_smooth = float(np.corrcoef(smooth3(model)[test_window], smooth3(psmsl_massaged)[test_window])[0, 1])
    print(f"[RESULT] same, 3-month-smoothed: r={r_test_smooth:+.4f}")

    rng = np.random.default_rng(42)
    n = len(model)
    null_rs = []
    for _ in range(500):
        shift = rng.integers(200, n - 200)
        null_rs.append(float(np.corrcoef(np.roll(model, shift)[test_window], psmsl_massaged[test_window])[0, 1]))
    null_rs = np.array(null_rs)
    z = (r_test - null_rs.mean()) / null_rs.std()
    print(f"[null] 500 random circular-shift nulls: mean={null_rs.mean():+.4f} "
          f"sd={null_rs.std():.4f}  ->  z={z:+.2f}")

    fig, axes = plt.subplots(2, 1, figsize=(12, 7.5))
    show = (proj_dates >= 1835) & (proj_dates < 1900)
    axes[0].plot(proj_dates[show], psmsl_massaged[show], color="0.25", linewidth=1.0,
                 label="real PSMSL Brest, massage.py-processed (never used to fit anything)")
    axes[0].plot(proj_dates[show], model[show], color="teal", linewidth=1.0, alpha=0.85,
                 label="Step-1 brestexcl model, extrapolated backward (never fit to this era)")
    axes[0].axvspan(1840, 1880, color="crimson", alpha=0.08, label="1840-1880: the test window")
    axes[0].axvspan(1835, 1840, color="0.85", label="excluded (pre-1840 step-change era)")
    axes[0].set_title(f"Brest: massage.py-processed comparison, 1840-1880 "
                       f"(r={r_test:+.3f}, z={z:+.1f} vs null)",
                       fontsize=10, fontweight="bold")
    axes[0].legend(fontsize=8, loc="upper left")

    axes[1].hist(null_rs, bins=40, color="0.7", label="null (random phase shift)")
    axes[1].axvline(r_test, color="crimson", linewidth=2, label=f"observed r={r_test:+.3f}")
    axes[1].set_title("null distribution", fontsize=10)
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    out_path = ROOT / IDX / "brestexcl_pre1880_massage.png"
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"[saved] {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
