#!/usr/bin/env python3
"""tpi_winding_regression.py — closed-form non-autonomous winding/
quantization transfer-function fit for TPI (Pacific tripole index),
matching the approach already validated for PNA and NINO4.

IMPORTANT DIFFERENCE FROM NINO4/PNA, found by reading the resp file
before blindly reusing cv_rolling_blocked.fit_score(): tpi/lt.exe.resp
sets UNCOMPENSATED=TRUE. Per src/gem-lte-primitives.adb and
src/gem-lte-primitives-solution.adb, this is NOT a no-op flag -- it
switches the production regression to a structurally different (and,
per the Ada comments, "jointly optimal" vs "provably non-optimal")
scheme:

  standard mode (nino4/pna, UNCOMPENSATED unset/False):
    fit raw_modulation (the plain sin/cos+trend basis) against a
    PRE-EMPHASIZED target DR(t) = Data(t) + IR*Data(t-12), then
    reconstruct Model_final(t) = raw_modulation(t) - IR*raw_modulation(t-12)
    as a first-order approximate inverse. This is what
    cv_rolling_blocked.regression_factors()/fit_score() implements.

  uncompensated mode (tpi, UNCOMPENSATED=True):
    fold the SAME lag-12 operator directly into every column of the
    design matrix (X_L(t) = X(t) - IR*X(t-12)), and regress X_L
    directly against the RAW (untransformed) target. This is a single
    jointly-optimal regression (L is not orthogonal, so pre-emphasizing
    the target and fitting is a genuinely different, suboptimal
    solution) -- not implemented anywhere in this project's Python
    reimplementation, so it has to be built here by hand.

By linearity, Model_final(t) = raw_modulation(t) - IR*raw_modulation(t-12)
turns out to be the SAME reconstruction formula in both modes -- the
difference is entirely in which side of the regression the lag-12
operator is applied to during fitting, which changes the fitted
coefficients (and, per the Ada comment, the fit's actual optimality).
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cv_rolling_blocked import prepare, fit_score, blocked_cv, valid_range  # noqa: E402
from lte_forward import (regression_factors, lte_response, filter9point,  # noqa: E402
                          Config, read_resp)
from signal_operators import delayed_difference  # noqa: E402
import json

IDX = "tpi"
ROOT = Path(__file__).resolve().parent


def is_uncompensated(idx: str) -> bool:
    params = json.loads((ROOT / idx / "lt.exe.p").read_text())
    resp = read_resp(ROOT / idx / "lt.exe.resp")
    cfg = Config(resp, {})
    return cfg.get("UNCOMPENSATED", False)


def build_design(t, fv, m, trend_on):
    """Exact column order/formula match to lte_forward.regression_factors,
    but returns the raw matrix (not the fitted coefficients) so it can be
    transformed before fitting."""
    cols = [np.ones_like(fv), fv]
    for k in range(len(m)):
        cols.append(np.sin(2 * math.pi * m[k] * fv))
        cols.append(np.cos(2 * math.pi * m[k] * fv))
    if trend_on:
        cols += [np.cos(2 * math.pi * t), np.sin(2 * math.pi * t),
                 np.cos(4 * math.pi * t), np.sin(4 * math.pi * t),
                 t, (t - t[0]) ** 2.0]
    return np.column_stack(cols)


def one_shot_lag_diff(x: np.ndarray, ir: float) -> np.ndarray:
    """x(t) - ir*x(t-12), first 12 rows unchanged -- matches the reverse-
    order Ada/Python loop (each output uses the ORIGINAL, untransformed
    value 12 back, not a recursively-corrected one)."""
    out = x.copy()
    out[12:] = x[12:] - ir * x[:-12]
    return out


def fit_uncompensated(prep: dict, train_idx: np.ndarray, test_idx: np.ndarray
                       ) -> dict:
    dates, data, forcing, m = prep["dates"], prep["data"], prep["forcing"], prep["m"]
    ir = prep["ir"]
    A = build_design(dates, forcing, m, prep["trend_on"])
    A_L = np.apply_along_axis(one_shot_lag_diff, 0, A, ir)
    tr = train_idx[data[train_idx] != 0.0]
    te = test_idx[data[test_idx] != 0.0]
    beta = np.linalg.solve(A_L[tr].T @ A_L[tr], A_L[tr].T @ data[tr])
    raw_modulation = A @ beta
    model = one_shot_lag_diff(raw_modulation, ir) if ir != 0.0 else raw_modulation
    for _ in range(prep["f9"]):
        model = filter9point(model)
    train_r = float(np.corrcoef(model[tr], data[tr])[0, 1])
    val_r = float(np.corrcoef(model[te], data[te])[0, 1])
    return dict(train_r=train_r, val_r=val_r, model=model, beta=beta)


def fit_standard(prep, train_idx, test_idx):
    """The OLD (default/non-uncompensated) regression, applied anyway --
    what you'd get if you naively reused nino4's exact recipe without
    checking the resp file first."""
    dates, data, forcing, m = prep["dates"], prep["data"], prep["forcing"], prep["m"]
    level, k0, amp, phase, annual, trend, accel = regression_factors(
        dates[train_idx], forcing[train_idx], data[train_idx], m, prep["trend_on"])
    model = lte_response(forcing, dates, m, amp, phase, level, k0, trend,
                          accel, prep["nonlin"], annual)
    if prep["ir"] != 0.0:
        model = model.copy()
        for i in range(len(model) - 1, 11, -1):
            model[i] -= prep["ir"] * model[i - 12]
    for _ in range(prep["f9"]):
        model = filter9point(model)
    tr = train_idx[data[train_idx] != 0.0]
    te = test_idx[data[test_idx] != 0.0]
    return dict(train_r=float(np.corrcoef(model[tr], data[tr])[0, 1]),
                val_r=float(np.corrcoef(model[te], data[te])[0, 1]),
                model=model)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--unfiltered", action="store_true")
    args = ap.parse_args()

    uncomp = is_uncompensated(IDX)
    prep = prepare(IDX)
    if args.unfiltered:
        prep = dict(prep)
        prep["data"] = prep["data_raw"]
        prep["f9"] = 0
    dates, data, forcing, m = prep["dates"], prep["data"], prep["forcing"], prep["m"]
    print(f"[{IDX}] UNCOMPENSATED (from resp)={uncomp}")
    print(f"[{IDX}] {len(m)} winding terms: {np.round(m, 4).tolist()}")
    print(f"[{IDX}] nonlin={prep['nonlin']}  trend_on={prep['trend_on']}  "
          f"f9={prep['f9']}  ir={prep['ir']:.4f}")
    t_min, t_max = valid_range(prep)
    print(f"[{IDX}] valid date range: {t_min:.2f} - {t_max:.2f}")

    test_idx = np.nonzero((dates >= 2000.0) & (dates < 2010.0))[0]
    train_idx = np.nonzero((dates < 2000.0) | (dates >= 2010.0))[0]

    std = fit_standard(prep, train_idx, test_idx)
    unc = fit_uncompensated(prep, train_idx, test_idx)
    print(f"[{IDX}] 2000-2010 holdout, STANDARD-mode regression "
          f"(ignores UNCOMPENSATED flag):")
    print(f"    train_r={std['train_r']:+.4f}   val_r={std['val_r']:+.4f}")
    print(f"[{IDX}] 2000-2010 holdout, UNCOMPENSATED-mode regression "
          f"(matches production's actual resp setting):")
    print(f"    train_r={unc['train_r']:+.4f}   val_r={unc['val_r']:+.4f}")

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

    # full-record robustness, using whichever mode matches production
    best_fit = fit_uncompensated if uncomp else fit_standard
    folds = []
    t_min, t_max = valid_range(prep)
    edges = np.linspace(t_min, t_max, 7)
    for k in range(6):
        b0, b1 = edges[k], edges[k + 1]
        te = np.nonzero((dates >= b0) & (dates < b1))[0]
        tr = np.nonzero((dates < b0 - 0.5) | (dates >= b1 + 0.5))[0]
        r = best_fit(prep, tr, te)
        folds.append((b0, b1, r["val_r"]))
    print(f"[{IDX}] 6-block CV across full record "
          f"({'uncompensated' if uncomp else 'standard'} mode):")
    for b0, b1, cc in folds:
        print(f"    test={b0:.1f}-{b1:.1f}  cc={cc:+.4f}")
    ccs = np.array([f[2] for f in folds])
    print(f"    -> mean={ccs.mean():+.4f}  median={np.median(ccs):+.4f}  "
          f"min={ccs.min():+.4f}  max={ccs.max():+.4f}")

    # 12mo delayed-difference OF FORCING (distinct check from the
    # IR/UNCOMPENSATED basis-folding above -- same diagnostic as nino4)
    forcing_diff = np.concatenate([np.zeros(12), delayed_difference(forcing, 12)])
    prep_diff = dict(prep)
    prep_diff["forcing"] = forcing_diff
    diff_fit = best_fit(prep_diff, train_idx, test_idx)
    print(f"[{IDX}] 12mo delayed-difference-of-FORCING variant "
          f"({'uncompensated' if uncomp else 'standard'} mode, same term "
          f"count): train_r={diff_fit['train_r']:+.4f}  "
          f"val_r={diff_fit['val_r']:+.4f}  (plain-forcing val_r="
          f"{(unc if uncomp else std)['val_r']:+.4f} for reference)")

    # plot using the mode that matches production
    chosen = unc if uncomp else std
    model_full = chosen["model"]
    fig, axes = plt.subplots(2, 1, figsize=(11, 7.5))
    nz = data != 0.0
    label = "raw (unfiltered)" if args.unfiltered else "F9-filtered"
    axes[0].plot(dates[nz], data[nz], color="0.25", linewidth=0.8,
                 label=f"real TPI ({label})")
    axes[0].plot(dates[nz], model_full[nz], color="teal", linewidth=0.9,
                 alpha=0.85,
                 label=f"closed-form winding-transfer model "
                       f"({'uncompensated' if uncomp else 'standard'} mode, "
                       f"trained outside 2000-2010)")
    axes[0].axvspan(2000, 2010, color="gold", alpha=0.15, label="held-out window")
    axes[0].set_title(f"TPI ({label}): closed-form non-autonomous "
                       f"winding-transfer fit, full record  "
                       f"(holdout val_r={chosen['val_r']:+.3f})",
                       fontsize=10, fontweight="bold")
    axes[0].legend(fontsize=8, loc="upper left")
    axes[0].set_xlim(1880, 2025)

    zmask = (dates >= 1998) & (dates < 2012) & nz
    axes[1].plot(dates[zmask], data[zmask], color="0.25", linewidth=1.3,
                 marker="o", markersize=2, label=f"real TPI ({label})")
    axes[1].plot(dates[zmask], model_full[zmask], color="teal", linewidth=1.3,
                 alpha=0.85, label="model (held out here)")
    axes[1].axvspan(2000, 2010, color="gold", alpha=0.15)
    axes[1].set_title("zoom: 1998-2012 (2000-2010 = never seen by the fit)",
                       fontsize=10)
    axes[1].legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    out_name = ("tpi_winding_regression_unfiltered.png" if args.unfiltered
                else "tpi_winding_regression.png")
    out_path = ROOT / IDX / out_name
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"[{IDX}] saved {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
