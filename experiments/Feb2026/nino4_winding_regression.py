#!/usr/bin/env python3
"""nino4_winding_regression.py — apply the closed-form, non-autonomous
winding/quantization transfer-function approach (the geoenergymath.com
"winding scalogram" framework, already validated exactly this way for
PNA this session) to NINO4, which has so far only been attempted via the
much more expensive resolved-PDE/beta-plane route.

Reuses cv_rolling_blocked.prepare()/fit_score() unmodified -- this
already encodes every lesson learned the hard way for PNA: the EXACT
production winding-number term list (prep['m'], not an approximate
nearest-neighbor match), the F9-FILTERED target (prep['data'], not
prep['data_raw']), the IR lag-12 decompensation step, and the ordinary
trend/seasonal regressors -- all in one already-tested tool.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cv_rolling_blocked import prepare, fit_score, blocked_cv, valid_range  # noqa: E402
from lte_forward import regression_factors, lte_response, filter9point  # noqa: E402
from signal_operators import delayed_difference  # noqa: E402

IDX = "nino4"


def fit_variant(prep: dict, forcing_variant: np.ndarray, train_idx: np.ndarray,
                 test_idx: np.ndarray) -> dict | None:
    """Same zero-new-DOF regression as fit_score, but with an arbitrary
    forcing array substituted in (e.g. the 12mo delayed-difference) instead
    of prep['forcing'] -- lets us test WHICH forcing variant the data
    prefers using the identical term count/structure, not add new terms."""
    dates, data, m = prep["dates"], prep["data"], prep["m"]
    train_idx = train_idx[data[train_idx] != 0.0]
    test_idx = test_idx[data[test_idx] != 0.0]
    if len(train_idx) < 2 * (2 * len(m) + 7) or len(test_idx) < 3:
        return None
    level, k0, amp, phase, annual, trend, accel = regression_factors(
        dates[train_idx], forcing_variant[train_idx], data[train_idx], m,
        prep["trend_on"])
    model = lte_response(forcing_variant, dates, m, amp, phase, level, k0,
                          trend, accel, prep["nonlin"], annual)
    if prep["ir"] != 0.0:
        model = model.copy()
        for i in range(len(model) - 1, 11, -1):
            model[i] -= prep["ir"] * model[i - 12]
    for _ in range(prep["f9"]):
        model = filter9point(model)
    x, y = model[test_idx], data[test_idx]
    xt, yt = model[train_idx], data[train_idx]
    if np.std(x) == 0.0 or np.std(y) == 0.0:
        return None
    return dict(train_r=float(np.corrcoef(xt, yt)[0, 1]),
                val_r=float(np.corrcoef(x, y)[0, 1]))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--unfiltered", action="store_true",
                     help="regress against the RAW (pre-F9) data instead of "
                          "the F9-filtered production target, and don't "
                          "apply F9 to the model output either -- an "
                          "apples-to-apples raw-vs-raw comparison.")
    args = ap.parse_args()

    prep = prepare(IDX)
    if args.unfiltered:
        prep = dict(prep)  # shallow copy -- don't mutate prepare()'s result
        prep["data"] = prep["data_raw"]
        prep["f9"] = 0
        label = "RAW (unfiltered)"
        out_name = "nino4_winding_regression_unfiltered.png"
    else:
        label = "F9-filtered"
        out_name = "nino4_winding_regression.png"
    dates, data, forcing, m = prep["dates"], prep["data"], prep["forcing"], prep["m"]
    print(f"[{IDX}] {len(m)} winding terms (exact production set): "
          f"{np.round(m, 4).tolist()}")
    print(f"[{IDX}] nonlin={prep['nonlin']}  trend_on={prep['trend_on']}  "
          f"f9={prep['f9']}  ir={prep['ir']:.4f}")
    t_min, t_max = valid_range(prep)
    print(f"[{IDX}] valid date range: {t_min:.2f} - {t_max:.2f}")

    # --- match PNA's own validated convention: holdout = 2000-2010 -------
    test_idx = np.nonzero((dates >= 2000.0) & (dates < 2010.0))[0]
    train_idx = np.nonzero((dates < 2000.0) | (dates >= 2010.0))[0]
    result = fit_score(prep, train_idx, test_idx)
    if result is None:
        print("[FAIL] fit_score returned None -- not enough usable rows "
              "in one side of the split.")
    else:
        # also report TRAIN correlation for the same fit, matching how the
        # PNA correction reported both train_r and val_r
        level, k0, amp, phase, annual, trend, accel = regression_factors(
            dates[train_idx], forcing[train_idx], data[train_idx], m,
            prep["trend_on"])
        model = lte_response(forcing, dates, m, amp, phase, level, k0, trend,
                              accel, prep["nonlin"], annual)
        if prep["ir"] != 0.0:
            model = model.copy()
            for i in range(len(model) - 1, 11, -1):
                model[i] -= prep["ir"] * model[i - 12]
        for _ in range(prep["f9"]):
            model = filter9point(model)
        tr = train_idx[data[train_idx] != 0.0]
        train_r = float(np.corrcoef(model[tr], data[tr])[0, 1])
        print(f"[{IDX}] ({label}) 2000-2010 holdout (matches "
              f"nino4site2000-2010.png convention):")
        print(f"    train_r={train_r:+.4f}   val_r={result['cc']:+.4f}   "
              f"n_train={result['n_train']}  n_test={result['n_test']}")

    # --- production's OWN reported skill on the identical window, for a
    # direct, apples-to-apples reference point -----------------------------
    lte_csv = Path(__file__).resolve().parent / IDX / "lte_results.csv"
    if lte_csv.exists():
        raw = np.loadtxt(lte_csv, delimiter=",")
        d2, model_col = raw[:, 0], raw[:, 1]
        # align production's own model column onto the SAME F9-filtered
        # `data`/`dates` arrays by nearest date (site outputs are on the
        # same monthly grid, but guard against any off-by-one padding)
        idx_map = np.searchsorted(dates, d2)
        idx_map = np.clip(idx_map, 0, len(dates) - 1)
        prod_test = np.nonzero((d2 >= 2000.0) & (d2 < 2010.0))[0]
        if len(prod_test) > 3:
            y = data[idx_map[prod_test]]
            x = model_col[prod_test]
            mask = y != 0.0
            if mask.sum() > 3 and np.std(x[mask]) > 0 and np.std(y[mask]) > 0:
                prod_r = float(np.corrcoef(x[mask], y[mask])[0, 1])
                print(f"[{IDX}] production lte_results.csv model-column "
                      f"correlation over the same 2000-2010 window "
                      f"(reference target): {prod_r:+.4f}")

    # --- broader robustness check: 6-block CV across the whole record ----
    folds = blocked_cv(prep, n_blocks=6, embargo_years=0.5)
    print(f"[{IDX}] 6-block CV across full record:")
    for i, f in enumerate(folds):
        print(f"    block {i}  test={f['test_span'][0]:.1f}-"
              f"{f['test_span'][1]:.1f}  n_test={f['n_test']:4d}  "
              f"cc={f['cc']:+.4f}")
    if folds:
        ccs = np.array([f["cc"] for f in folds])
        print(f"    -> mean={ccs.mean():+.4f}  median={np.median(ccs):+.4f}  "
              f"min={ccs.min():+.4f}  max={ccs.max():+.4f}")

    # --- 12-month delayed-difference check (per direct request) -----------
    # Same zero-new-DOF regression, same m/term count, only the forcing
    # array is swapped for its 12mo delayed difference (the same operator
    # that was NECESSARY for NAO/brestexcl to reveal their high-frequency
    # character) -- report honestly whether nino4 wants it, is indifferent,
    # or is hurt by it.
    forcing_diff = np.concatenate([np.zeros(12), delayed_difference(forcing, 12)])
    diff_result = fit_variant(prep, forcing_diff, train_idx, test_idx)
    print(f"[{IDX}] 12mo delayed-difference variant, SAME 2000-2010 holdout, "
          f"same term count (zero new DOF):")
    if diff_result is None:
        print("    (not enough usable rows)")
    else:
        print(f"    train_r={diff_result['train_r']:+.4f}   "
              f"val_r={diff_result['val_r']:+.4f}   "
              f"(plain-forcing val_r={result['cc']:+.4f} for reference)")

    # --- validation plot ---------------------------------------------------
    level, k0, amp, phase, annual, trend, accel = regression_factors(
        dates[train_idx], forcing[train_idx], data[train_idx], m,
        prep["trend_on"])
    model_full = lte_response(forcing, dates, m, amp, phase, level, k0, trend,
                               accel, prep["nonlin"], annual)
    if prep["ir"] != 0.0:
        model_full = model_full.copy()
        for i in range(len(model_full) - 1, 11, -1):
            model_full[i] -= prep["ir"] * model_full[i - 12]
    for _ in range(prep["f9"]):
        model_full = filter9point(model_full)

    fig, axes = plt.subplots(2, 1, figsize=(11, 7.5))
    nz = data != 0.0
    axes[0].plot(dates[nz], data[nz], color="0.25", linewidth=0.8,
                 label=f"real NINO4 ({label})")
    axes[0].plot(dates[nz], model_full[nz], color="crimson", linewidth=0.9,
                 alpha=0.85,
                 label="closed-form winding-transfer model "
                       "(trained outside 2000-2010)")
    axes[0].axvspan(2000, 2010, color="gold", alpha=0.15,
                     label="held-out window")
    axes[0].set_title(f"NINO4 ({label}): closed-form non-autonomous "
                       f"winding-transfer fit, full record  "
                       f"(holdout val_r={result['cc']:+.3f})",
                       fontsize=10, fontweight="bold")
    axes[0].legend(fontsize=8, loc="upper left")
    axes[0].set_xlim(1880, 2025)

    zmask = (dates >= 1998) & (dates < 2012) & nz
    axes[1].plot(dates[zmask], data[zmask], color="0.25", linewidth=1.3,
                 marker="o", markersize=2, label=f"real NINO4 ({label})")
    axes[1].plot(dates[zmask], model_full[zmask], color="crimson",
                 linewidth=1.3, alpha=0.85, label="model (held out here)")
    axes[1].axvspan(2000, 2010, color="gold", alpha=0.15)
    axes[1].set_title("zoom: 1998-2012 (2000-2010 = never seen by the fit)",
                       fontsize=10)
    axes[1].legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    out_path = Path(__file__).resolve().parent / IDX / out_name
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"[{IDX}] saved {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
