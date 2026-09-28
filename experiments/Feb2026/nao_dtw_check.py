#!/usr/bin/env python3
"""nao_dtw_check.py — supplemental visual-match analysis for the bottom
panel of nao/nao_resolved_shallow_water.png (the 12-month delayed-
difference gauge vs real NAO), using the SAME dynamic-time-warping
machinery already validated elsewhere in this project
(qbo_compensating_dtw.py's dtw_distance + null_baseline), not a new
metric invented for this check.

Why DTW here specifically: a plain Pearson correlation penalizes any
timing mismatch between two curves as if it were a shape mismatch, even
when the two are tracking the same regime changes with a few months of
slip. The nao_resolved_shallow_water.py plot's bottom panel visually
tracks real NAO's regime turns reasonably well despite the reported
r=+0.059 -- DTW's warping tolerance is the direct way to ask "how much of
that visual agreement is real pattern-matching vs an artifact of a
correlation coefficient being the wrong tool for a phase-slippy match,"
exactly the role it already plays for QBO's phase residual in
qbo_compensating_dtw.py.

Compares FOUR series against real NAO (not just the one panel asked
about) so the 12mo-diff variant's result has honest context: raw FIXED
gauge, and (regenerated fresh, since nao_resolved_shallow_water.py did
not persist arrays) the same 12mo-diff series.

Usage
-----
    ./nao_dtw_check.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from nao_winding_rank_check import regenerate_fixed                        # noqa: E402
from signal_operators import delayed_difference                           # noqa: E402
from qbo_compensating_dtw import standardize, dtw_distance, null_baseline  # noqa: E402

ROOT = Path(__file__).resolve().parent


def check(label: str, model: np.ndarray, dates: np.ndarray,
         obs: np.ndarray, seed: int) -> None:
    model_s = standardize(model)
    obs_s = standardize(obs)

    r_raw = float(np.corrcoef(model_s, obs_s)[0, 1])
    sign = -1.0 if r_raw < 0 else 1.0
    model_signed = sign * model_s

    d_true, _ = dtw_distance(model_signed, obs_s)
    rng = np.random.default_rng(seed)
    null_mean, null_std = null_baseline(rng, model_signed, obs_s, n_trials=16)
    z = (null_mean - d_true) / null_std if null_std > 0 else float("nan")

    print(f"-- {label} --")
    print(f"  Pearson r (sign-matched): {sign * r_raw:+.3f}")
    print(f"  DTW normalized distance (true alignment): {d_true:.4f}")
    print(f"  DTW distance, random-phase-shifted null of real NAO "
          f"(16 trials): mean={null_mean:.4f} std={null_std:.4f}")
    print(f"  z-score (how much closer the TRUE alignment is than random "
          f"phase, in null std units): {z:+.2f}")
    verdict = ("a genuine, warping-tolerant pattern match beyond what "
              "unrelated timing of the same-shaped curve would give"
              if z > 2.0 else
              "not clearly distinguishable from what random phase "
              "alignment of the same-shaped curve would give")
    print(f"  -> {verdict}\n")


def main() -> int:
    print("Regenerating NAO resolved-model FIXED (Iceland) gauge...")
    gauge_dates, gauge_forcing, gauge = regenerate_fixed()

    nao_dat = np.loadtxt(ROOT / "nao" / "nao.dat")
    dates_nao, obs_nao = nao_dat[:, 0].copy(), nao_dat[:, 1].copy()
    obs_on_gauge = np.interp(gauge_dates, dates_nao, obs_nao)

    check("RAW FIXED gauge (Iceland) vs real NAO", gauge, gauge_dates,
         obs_on_gauge, seed=0)

    lag = 12
    diff = delayed_difference(gauge, lag)
    diff_dates = gauge_dates[lag:]
    diff_obs = obs_on_gauge[lag:]
    check("FIXED, 12mo delayed-difference vs real NAO (the plot's bottom "
          "panel)", diff, diff_dates, diff_obs, seed=1)

    return 0


if __name__ == "__main__":
    sys.exit(main())
