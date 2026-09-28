#!/usr/bin/env python3
"""brestexcl_full_cycle_overlay.py — STEP 1 of a redo, per direct
feedback that the previous script's raw-PSMSL refit only weakly tracked
even the data it was fit on. This script does the obvious thing first:
fit the closed-form winding model against brestexcl's OWN (filtered)
data -- exactly the same recipe already validated for nino4/tpi/iode/
tna -- and visually confirm it actually tracks brestexcl.dat well over
the FULL 1880-2020 record BEFORE extending anything backward or making
any pre-1880 claim.
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
from lte_forward import regression_factors, lte_response, filter9point  # noqa: E402
from brestexcl_resolved_shallow_water import build_forcing_continuous  # noqa: E402
from lte_forward import year_length, Config, read_resp  # noqa: E402
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
    prep = prepare(IDX)
    dates, data, forcing, m = prep["dates"], prep["data"], prep["forcing"], prep["m"]
    nz = data != 0.0

    # --- fit against brestexcl's OWN data (its actual, filtered target),
    # using ALL available rows -- the same recipe already validated for
    # nino4/tpi/iode/tna, not the raw-PSMSL refit from the last script --
    level, k0, amp, phase, annual, trend, accel = regression_factors(
        dates[nz], forcing[nz], data[nz], m, prep["trend_on"])
    model_raw = lte_response(forcing, dates, m, amp, phase, level, k0,
                              trend, accel, prep["nonlin"], annual)
    model = decompensate_f9(model_raw, prep["ir"], prep["f9"])

    train_r = float(np.corrcoef(model[nz], data[nz])[0, 1])
    print(f"[{IDX}] fit against brestexcl's OWN data, full 1880-2020 "
          f"in-sample correlation: r={train_r:+.4f}")

    # also the honest 2000-2010-style holdout, matching every other index
    # checked this session, for a fair apples-to-apples number
    test_idx = np.nonzero((dates >= 2000.0) & (dates < 2010.0) & nz)[0]
    train_idx = np.nonzero(((dates < 2000.0) | (dates >= 2010.0)) & nz)[0]
    lvl2, k02, amp2, ph2, ann2, tr2, ac2 = regression_factors(
        dates[train_idx], forcing[train_idx], data[train_idx], m, prep["trend_on"])
    model2_raw = lte_response(forcing, dates, m, amp2, ph2, lvl2, k02, tr2,
                               ac2, prep["nonlin"], ann2)
    model2 = decompensate_f9(model2_raw, prep["ir"], prep["f9"])
    val_r = float(np.corrcoef(model2[test_idx], data[test_idx])[0, 1])
    print(f"[{IDX}] 2000-2010 GENUINE holdout (fit excludes it): "
          f"val_r={val_r:+.4f}")

    # --- visual check #1, the thing asked for: does this model actually
    # track brestexcl.dat well across the WHOLE record it was fit to? ----
    fig, axes = plt.subplots(2, 1, figsize=(12, 7))
    axes[0].plot(dates[nz], data[nz], color="0.25", linewidth=0.8,
                 label="real brestexcl.dat")
    axes[0].plot(dates[nz], model[nz], color="teal", linewidth=0.9, alpha=0.85,
                 label="closed-form winding model (fit on ALL 1880-2020 data)")
    axes[0].set_title(f"brestexcl: full-record in-sample fit  (r={train_r:+.3f}) "
                       f"-- does this actually track, unlike the raw-PSMSL refit?",
                       fontsize=10, fontweight="bold")
    axes[0].legend(fontsize=8, loc="upper left")

    zmask = (dates >= 1998) & (dates < 2012) & nz
    axes[1].plot(dates[zmask], data[zmask], color="0.25", linewidth=1.3,
                 marker="o", markersize=2, label="real brestexcl.dat")
    axes[1].plot(dates[zmask], model2[zmask], color="crimson", linewidth=1.3,
                 alpha=0.85, label=f"model, 2000-2010 GENUINELY held out (val_r={val_r:+.3f})")
    axes[1].axvspan(2000, 2010, color="gold", alpha=0.15)
    axes[1].set_title("zoom 1998-2012, with a real (not in-sample) holdout window",
                       fontsize=10)
    axes[1].legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    out1 = ROOT / IDX / "brestexcl_full_cycle_step1_insample_check.png"
    fig.savefig(out1, dpi=140)
    plt.close(fig)
    print(f"[saved] {out1}")

    if train_r < 0.4:
        print("[STOP] in-sample fit is weak -- NOT proceeding to the "
              "pre-1880 extension until this is understood; reporting "
              "only the visual check as requested.")
        return 0

    # --- only now, having confirmed the model tracks brestexcl.dat well,
    # extend the SAME (all-data) fit backward through the PSMSL pre-1880
    # era and plot the WHOLE 1807-2020 cycle continuously. Raw PSMSL is
    # shown alongside for qualitative visual reference ONLY (already
    # established: it is not directly comparable in value to
    # brestexcl.dat, r=0.26) -- no correlation claim is made on it here.
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

    params = json.loads((ROOT / IDX / "lt.exe.p").read_text())
    resp = read_resp(ROOT / IDX / "lt.exe.resp")
    cfg = Config(resp, {})
    yl = year_length(cfg.get("YEAR", 0.0), float(params.get("year", 0.0)))
    forcing_ext = build_forcing_continuous(IDX, yl, dates_psmsl)

    # SAME all-data fit (level,k0,amp,phase,trend,accel from above),
    # evaluated over the extended date range -- this is a pure
    # extrapolation of the ALREADY-FIT brestexcl model, nothing refit.
    model_ext_raw = lte_response(forcing_ext, dates_psmsl, m, amp, phase,
                                  level, k0, trend, accel, prep["nonlin"], annual)
    model_ext = decompensate_f9(model_ext_raw, prep["ir"], prep["f9"])

    fig2, ax = plt.subplots(figsize=(13, 5.5))
    ax.plot(dates[nz], data[nz], color="0.25", linewidth=0.8,
            label="real brestexcl.dat (1880-2020, what the model was fit to)")
    ax.plot(dates_psmsl, model_ext, color="teal", linewidth=0.8, alpha=0.85,
            label="SAME model, extrapolated across the entire 1807-2020 cycle")
    ax_r = ax.twinx()
    show_psmsl = ~missing
    ax_r.plot(dates_psmsl[show_psmsl], vals_psmsl[show_psmsl], color="0.75",
              linewidth=0.5, alpha=0.6, zorder=0,
              label="raw PSMSL (right axis, different units/processing -- "
                    "qualitative reference only, NOT a claimed match)")
    ax_r.set_ylabel("raw PSMSL (mm, arbitrary datum)", color="0.5")
    ax.axvspan(1807, 1880, color="crimson", alpha=0.06)
    ax.set_title(f"brestexcl: the SAME all-data-fit model, whole 1807-2020 cycle "
                 f"(1880-2020 in-sample r={train_r:+.3f}; pre-1880 is a pure "
                 f"extrapolation, not yet scored against anything)",
                 fontsize=10, fontweight="bold")
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax_r.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, fontsize=7, loc="upper left")
    fig2.tight_layout()
    out2 = ROOT / IDX / "brestexcl_full_cycle_overlay.png"
    fig2.savefig(out2, dpi=140)
    plt.close(fig2)
    print(f"[saved] {out2}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
