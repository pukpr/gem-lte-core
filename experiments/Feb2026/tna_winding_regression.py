#!/usr/bin/env python3
"""tna_winding_regression.py — closed-form non-autonomous winding-transfer
fit for TNA (Indian Ocean Dipole, eastern edge), following the same
method already validated for NINO4/TPI. TNA has a real, plausibly
accelerating secular rise (candidate cause: global warming) that must be
separated from the tidal/winding-driven variation rather than left to
inflate the apparent fit quality by coincidence -- the same trap already
caught once this session for PNA (a shared-trend artifact that collapsed
from DTW z=+5.13 to +1.19 after detrending).

TNA's resp has no UNCOMPENSATED override (standard mode, like nino4/pna)
and IR=0.028 is small/stable -- no TPI-style surprise here. The one
genuinely new piece of analysis this script adds: decompose the fitted
model into a TREND-ONLY component (the secular date/date^2 terms alone)
and a WINDINGS+SEASONAL-ONLY component (everything else), propagate the
SAME (IR-decompensate, then F9) linear post-processing to each piece
separately (valid because both steps are linear, so the pieces sum back
exactly to the full production-matching model), then score the
windings-only component against TREND-REMOVED real data -- an honest
test of whether the tidal/winding structure explains genuine variation,
not just riding a shared secular rise.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cv_rolling_blocked import prepare, blocked_cv, valid_range  # noqa: E402
from lte_forward import regression_factors, lte_response, filter9point  # noqa: E402
from signal_operators import delayed_difference  # noqa: E402

IDX = "tna"
ROOT = Path(__file__).resolve().parent


def decompensate_f9(x: np.ndarray, ir: float, f9: int) -> np.ndarray:
    """Same one-shot lag-12 decompensation + F9 pipeline fit_score applies
    to the full model -- factored out so it can be applied identically to
    sub-components and still sum back exactly (both steps are linear)."""
    out = x.copy()
    if ir != 0.0:
        for i in range(len(out) - 1, 11, -1):
            out[i] -= ir * out[i - 12]
    for _ in range(f9):
        out = filter9point(out)
    return out


def main() -> int:
    prep = prepare(IDX)
    dates, data, forcing, m = prep["dates"], prep["data"], prep["forcing"], prep["m"]
    print(f"[{IDX}] {len(m)} winding terms: {np.round(m, 4).tolist()}")
    print(f"[{IDX}] nonlin={prep['nonlin']}  trend_on={prep['trend_on']}  "
          f"f9={prep['f9']}  ir={prep['ir']:.4f}")
    t_min, t_max = valid_range(prep)
    print(f"[{IDX}] valid date range: {t_min:.2f} - {t_max:.2f}")

    # --- empirical check: is the claimed rising/accelerating trend real? --
    # (dates centered on the record's first valid date -- an uncentered
    # fit over ~1880-2022 makes the linear/quadratic terms numerically
    # collinear and uninterpretable)
    nz = data != 0.0
    d0_check = dates[nz][0]
    quad, lin, _ = np.polyfit(dates[nz] - d0_check, prep["data_raw"][nz], 2)
    trend_1880, trend_1960, trend_2022 = (
        quad * (yr - d0_check) ** 2 + lin * (yr - d0_check)
        for yr in (1880, 1960, 2022))
    print(f"[{IDX}] empirical quadratic fit to RAW data alone (centered): "
          f"linear={lin:+.5f}/yr  quadratic={quad:+.7f}/yr^2 -- implied "
          f"level at 1880/1960/2022: {trend_1880:+.3f}/{trend_1960:+.3f}/"
          f"{trend_2022:+.3f}")

    test_idx = np.nonzero((dates >= 2000.0) & (dates < 2010.0))[0]
    train_idx = np.nonzero((dates < 2000.0) | (dates >= 2010.0))[0]
    tr = train_idx[data[train_idx] != 0.0]
    te = test_idx[data[test_idx] != 0.0]

    level, k0, amp, phase, annual, trend, accel = regression_factors(
        dates[tr], forcing[tr], data[tr], m, prep["trend_on"])
    print(f"[{IDX}] fitted secular terms: trend={trend:+.6f}/yr  "
          f"accel={accel:+.8f}/yr^2  level={level:+.4f}  k0={k0:+.4f}")

    # full production-matching model
    model_full_raw = lte_response(forcing, dates, m, amp, phase, level, k0,
                                   trend, accel, prep["nonlin"], annual)
    model_full = decompensate_f9(model_full_raw, prep["ir"], prep["f9"])

    # trend-only sub-component: JUST trend*date + accel*(date-date0)^2
    # (level/k0*F excluded -- those are a constant offset and the tidal
    # linear term, not "secular rise", and left in the windings side)
    trend_component_raw = trend * dates + accel * (dates - dates[0]) ** 2.0
    trend_component = decompensate_f9(trend_component_raw, prep["ir"], prep["f9"])

    # windings+seasonal+level+k0F component = full - trend, but reconstruct
    # directly (not by subtraction) for numerical cleanliness
    rest_raw = lte_response(forcing, dates, m, amp, phase, level, k0,
                             0.0, 0.0, prep["nonlin"], annual)
    rest = decompensate_f9(rest_raw, prep["ir"], prep["f9"])

    assert np.allclose(model_full, trend_component + rest, atol=1e-8), \
        "linearity check failed -- decomposition doesn't sum back exactly"

    data_detrended = data - trend_component

    train_r_full = float(np.corrcoef(model_full[tr], data[tr])[0, 1])
    val_r_full = float(np.corrcoef(model_full[te], data[te])[0, 1])
    val_r_trend_only = float(np.corrcoef(trend_component[te], data[te])[0, 1])
    val_r_windings_vs_detrended = float(
        np.corrcoef(rest[te], data_detrended[te])[0, 1])
    print(f"[{IDX}] 2000-2010 holdout, FULL production-matching model:")
    print(f"    train_r={train_r_full:+.4f}   val_r={val_r_full:+.4f}")
    print(f"[{IDX}] TREND-ONLY component vs raw real data (how much of the "
          f"'skill' could be just the secular rise): val_r={val_r_trend_only:+.4f}")
    print(f"[{IDX}] WINDINGS+SEASONAL-ONLY component vs TREND-REMOVED real "
          f"data (the honest test of genuine oscillatory explanatory "
          f"power): val_r={val_r_windings_vs_detrended:+.4f}")

    lte_csv = ROOT / IDX / "lte_results.csv"
    if lte_csv.exists():
        raw = np.loadtxt(lte_csv, delimiter=",")
        d2, model_col = raw[:, 0], raw[:, 1]
        idx_map = np.clip(np.searchsorted(dates, d2), 0, len(dates) - 1)
        prod_test = np.nonzero((d2 >= 2000.0) & (d2 < 2010.0))[0]
        if len(prod_test) > 3:
            y = data[idx_map[prod_test]]
            x = model_col[prod_test]
            mask = y != 0.0
            if mask.sum() > 3 and np.std(x[mask]) > 0 and np.std(y[mask]) > 0:
                prod_r = float(np.corrcoef(x[mask], y[mask])[0, 1])
                print(f"[{IDX}] production lte_results.csv reference over "
                      f"the same window: {prod_r:+.4f}")

    # 6-block CV, both full-model and windings-only-vs-detrended
    print(f"[{IDX}] 6-block CV across full record:")
    edges = np.linspace(*valid_range(prep), 7)
    full_ccs, wind_ccs = [], []
    for k in range(6):
        b0, b1 = edges[k], edges[k + 1]
        te_k = np.nonzero((dates >= b0) & (dates < b1) & nz)[0]
        tr_k = np.nonzero(((dates < b0 - 0.5) | (dates >= b1 + 0.5)) & nz)[0]
        lvl, k0k, ampk, phk, annk, trk, acck = regression_factors(
            dates[tr_k], forcing[tr_k], data[tr_k], m, prep["trend_on"])
        mf_raw = lte_response(forcing, dates, m, ampk, phk, lvl, k0k, trk,
                               acck, prep["nonlin"], annk)
        mf = decompensate_f9(mf_raw, prep["ir"], prep["f9"])
        tc_raw = trk * dates + acck * (dates - dates[0]) ** 2.0
        tc = decompensate_f9(tc_raw, prep["ir"], prep["f9"])
        rs_raw = lte_response(forcing, dates, m, ampk, phk, lvl, k0k, 0.0,
                               0.0, prep["nonlin"], annk)
        rs = decompensate_f9(rs_raw, prep["ir"], prep["f9"])
        dd = data - tc
        cc_full = float(np.corrcoef(mf[te_k], data[te_k])[0, 1])
        cc_wind = float(np.corrcoef(rs[te_k], dd[te_k])[0, 1])
        full_ccs.append(cc_full); wind_ccs.append(cc_wind)
        print(f"    test={b0:.1f}-{b1:.1f}  full_cc={cc_full:+.4f}  "
              f"windings-vs-detrended_cc={cc_wind:+.4f}")
    full_ccs, wind_ccs = np.array(full_ccs), np.array(wind_ccs)
    print(f"    -> full: mean={full_ccs.mean():+.4f}  "
          f"median={np.median(full_ccs):+.4f}")
    print(f"    -> windings-vs-detrended: mean={wind_ccs.mean():+.4f}  "
          f"median={np.median(wind_ccs):+.4f}")

    # 12mo delayed-difference-of-forcing, for consistency with prior indices
    forcing_diff = np.concatenate([np.zeros(12), delayed_difference(forcing, 12)])
    lvl2, k02, amp2, ph2, ann2, tr2, ac2 = regression_factors(
        dates[tr], forcing_diff[tr], data[tr], m, prep["trend_on"])
    md_raw = lte_response(forcing_diff, dates, m, amp2, ph2, lvl2, k02, tr2,
                           ac2, prep["nonlin"], ann2)
    md = decompensate_f9(md_raw, prep["ir"], prep["f9"])
    diff_val_r = float(np.corrcoef(md[te], data[te])[0, 1])
    print(f"[{IDX}] 12mo delayed-difference-of-forcing variant: "
          f"val_r={diff_val_r:+.4f}  (plain-forcing val_r={val_r_full:+.4f} "
          f"for reference)")

    # --- plot: full fit, and the trend-vs-windings decomposition ----------
    fig, axes = plt.subplots(3, 1, figsize=(11, 10.5))
    axes[0].plot(dates[nz], data[nz], color="0.25", linewidth=0.8,
                 label="real TNA (F9-filtered)")
    axes[0].plot(dates[nz], model_full[nz], color="darkorange", linewidth=0.9,
                 alpha=0.85, label="full model (trend+windings, held out 2000-2010)")
    axes[0].axvspan(2000, 2010, color="gold", alpha=0.15, label="held-out window")
    axes[0].set_title(f"TNA: full closed-form fit  (holdout val_r={val_r_full:+.3f})",
                       fontsize=10, fontweight="bold")
    axes[0].legend(fontsize=8, loc="upper left")

    # cosmetic re-centering ONLY for display -- correlations above are
    # computed on the un-recentered arrays and are unaffected by an
    # additive constant; this just keeps the plot's y-axis interpretable
    # instead of sitting on trend*dates' large intrinsic offset.
    # the lag-12 decompensation applied to the trend-only piece produces a
    # real but purely transient boundary spike in its first ~12 months
    # (not enough history yet for the lag term) -- exclude that settle-in
    # window from the display mean/range so it doesn't dominate the axis;
    # it sums back exactly into the full model regardless (already asserted)
    # and is nowhere near the 2000-2010 window the numbers are scored on.
    settle_mask = nz & (dates >= dates[nz][0] + 2.0)
    trend_display = trend_component - np.mean(trend_component[settle_mask])
    data_detrended_display = data - trend_display

    axes[1].plot(dates[settle_mask], data[settle_mask], color="0.6", linewidth=0.7,
                 label="real TNA (raw, F9-filtered)")
    axes[1].plot(dates[settle_mask], trend_display[settle_mask], color="firebrick",
                 linewidth=1.8, label="fitted TREND-ONLY component (nuisance rise, re-centered for display)")
    axes[1].set_title("secular trend component fitted jointly with the "
                       "windings (compensated out, not just eyeballed)",
                       fontsize=10)
    axes[1].legend(fontsize=8, loc="upper left")

    zmask = (dates >= 1998) & (dates < 2012) & nz
    axes[2].plot(dates[zmask], data_detrended_display[zmask], color="0.25",
                 linewidth=1.3, marker="o", markersize=2,
                 label="real TNA, TREND REMOVED")
    rest_display = rest + np.mean(trend_component[settle_mask])
    axes[2].plot(dates[zmask], rest_display[zmask], color="darkorange", linewidth=1.3,
                 alpha=0.85, label="windings+seasonal-only model (held out here)")
    axes[2].axvspan(2000, 2010, color="gold", alpha=0.15)
    axes[2].set_title(f"zoom 1998-2012: detrended real data vs windings-only "
                       f"model  (val_r={val_r_windings_vs_detrended:+.3f})",
                       fontsize=10)
    axes[2].legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    out_path = ROOT / IDX / "tna_winding_regression.png"
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"[{IDX}] saved {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
