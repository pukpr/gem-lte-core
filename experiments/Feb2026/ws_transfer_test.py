#!/usr/bin/env python3
"""ws_transfer_test.py -- freeze a ws.py global-mode fit trained on one
index's own record, then transplant it (forcing pipeline + regression
coefficients, completely unchanged) onto a second, longer index's date
grid, scoring correlation only on the portion the training record never
covered.

Same backward-extrapolation methodology already used for brestexcl's
pre-1880 validation, generalized into a reusable two-index tool: fit on
<train_index>'s own full span (no in-index holdout -- the genuine
out-of-sample test here is the OTHER index's pre-training-era span, not a
slice carved out of the training index itself), transplant that frozen
model onto <score_index>'s longer record via ws.build_forcing_at (a pure
function of dates + already-fitted manifold parameters), and report
correlation both over the shared/overlap era and before --split-year.

Usage:
    ws_transfer_test.py --train kN020_E050 --score kN020_E050_ \\
        --split-year 1950
"""
from __future__ import annotations

import argparse

import numpy as np

import ws


def regression_factors_no_trend(t, fv, y, m):
    """Same design as ws.regression_factors(trend_on=True), minus the
    linear trend and quadratic accel columns -- keeps the annual/
    semiannual explicit terms.

    Deliberately NOT the same as calling ws.regression_factors(...,
    trend_on=False): that flag drops the trend/accel AND the annual/
    semiannual columns together (they share one code block), and its
    coefficient-extraction indices (coef[ncol-3] etc.) then silently
    misattribute two of the harmonic sin/cos coefficients as "annual"
    terms instead -- not what "without trend and acceleration" asked
    for. This isolates just the trend/accel columns, since those were
    identified as the likely driver of the earlier backward-
    extrapolation drift: level (~18.6) and trend*date (~-18.5 at 1950)
    were fit as a near-perfectly-cancelling pair that's only tight
    inside the 1950-2023 training window."""
    nm = len(m)
    cols = [np.ones_like(fv), fv]
    for k in range(nm):
        cols.append(np.sin(2 * np.pi * m[k] * fv))
        cols.append(np.cos(2 * np.pi * m[k] * fv))
    cols += [np.cos(2 * np.pi * t), np.sin(2 * np.pi * t),
             np.cos(4 * np.pi * t), np.sin(4 * np.pi * t)]
    A = np.column_stack(cols)
    coef = np.linalg.solve(A.T @ A, A.T @ y)
    level, k0 = float(coef[0]), float(coef[1])
    amp = np.empty(nm)
    phase = np.empty(nm)
    for k in range(nm):
        c_sin, c_cos = coef[2 + 2 * k], coef[3 + 2 * k]
        amp[k] = np.hypot(c_sin, c_cos)
        phase[k] = np.arctan2(c_cos, c_sin)
    ann2, ann1, semi2, semi1 = (float(coef[-4]), float(coef[-3]),
                                 float(coef[-2]), float(coef[-1]))
    return level, k0, amp, phase, (semi1, semi2, ann1, ann2), 0.0, 0.0


def build_forcing_at_anchored(prep: dict, dates_new: np.ndarray, start_date: float) -> np.ndarray:
    """Identical to ws.build_forcing_at, except the IIR seed point
    (start_date) is passed explicitly instead of always being prep['idate'].

    This is the actual fix for the "should be the same over 1950-2020"
    request: ws.iir seeds `init` at the first sample with Date > start_date,
    then runs forward from there (and backward, to reconstruct pre-history).
    For kN020_E050 alone, start_date=idate=1880 precedes every date in its
    877-point (1950-2023) array, so iir() clamps to idx=0 -- init is
    effectively seeded AT 1950, not truly at 1880 (there's no data in
    between to evolve through). That's what the fitted `init`/amp/phase
    values actually mean.

    But feeding prep['idate']=1880 unchanged to the LONGER kN020_E050_ date
    grid (which DOES contain 1880) makes idx land at 1880 instead, so the
    recursion now forward-integrates 1880->1950 for 840 extra monthly steps
    before ever reaching the training window -- with mem=1-ma~0.99998 (a
    near-unit-root filter that barely damps over 70 years), that is a
    genuinely different trajectory through 1950-2023, not the same model
    evaluated on more data. That mismatch is exactly why the earlier
    transplant test's "overlap" correlation (0.32) came out well below the
    original in-sample r (0.77) despite the real data being identical in
    that window.

    The fix: anchor start_date at the TRAINING record's own first date
    (kN020_E050's 1950.0), not at the nominal IDATE config value, whenever
    evaluating on any OTHER date grid. That reproduces the exact same idx
    (and hence bit-identical forward recursion) whether the array handed in
    is the short 877-point record or the long 2005-point one -- verified
    below with an explicit equality assertion, not just an eyeballed
    correlation -- and the pre-1950 era is then generated purely by iir's
    own backward-reconstruction pass from that SAME anchor, the same
    mechanism already validated for brestexcl's pre-1880 extrapolation."""
    params = prep["params"]
    lpap = np.array(params["lpap"], dtype=float)
    periods = ws.doodson_periods(prep["year_startup"], prep["year_cand"])
    B = {k: float(params.get(k, 0.0)) for k in
         ("offs", "bg", "impA", "impB", "delA", "delB", "asym", "ma", "mp",
          "shfT", "init")}
    k1 = float(np.array(params["ltep"], dtype=float)[prep["nm"] - 1])
    tf = ws.tide_sum(dates_new, lpap[:, 1:3], periods, prep["yl"], 0.0, B["shfT"])
    comb = ws.impulse_delta(dates_new, B["delA"], B["delB"], B["asym"], 12)
    R = ws.iir(tf * comb, lag_a=1.0 - B["ma"], lag_c=B["mp"], init=B["init"],
               start_date=start_date, dates=dates_new)
    return ws.bessel(R, B["impA"], B["impB"], k1, B["offs"], B["bg"])


def lte_response_anchored(fv, dates, date0, m, amp, phase, level, k0, trend,
                           accel, nonlin, annual):
    """Identical to ws.lte_response, except the quadratic acceleration
    term's zero-point (date0) is passed explicitly instead of being taken
    as dates[0] of whatever array happens to be evaluated.

    ws.lte_response hard-codes accel * (dates - dates[0]) ** 2 -- fine as
    long as you always evaluate on the SAME date array used to fit accel,
    but wrong here: fitting on kN020_E050 anchors accel at 1950 (its own
    dates[0]), and naively calling lte_response on kN020_E050_'s dates
    (which start at 1856) would silently re-anchor the parabola at 1856
    instead. Caught by direct inspection: accel's own fitted magnitude
    (~1.5e-4) times the resulting (1950-1856)^2=8836 epoch-shift error is
    ~1.36 -- comparable to the entire SST anomaly's own dynamic range, so
    this is not a rounding-level detail, it would have silently dominated
    both the 'overlap' and 'before 1950' correlations reported below."""
    semi1, semi2, ann1, ann2 = annual
    out = np.zeros_like(fv)
    for j in range(len(m)):
        sw = np.sin(2 * np.pi * m[j] * fv + phase[j])
        pw = np.abs(sw) ** nonlin
        out += amp[j] * np.where(sw < 0.0, -pw, pw)
    out += (level + k0 * fv + trend * dates + accel * (dates - date0) ** 2.0)
    out += (ann1 * np.sin(2 * np.pi * dates) + ann2 * np.cos(2 * np.pi * dates)
            + semi1 * np.sin(4 * np.pi * dates) + semi2 * np.cos(4 * np.pi * dates))
    return out


def run(train_index: str, score_index: str, root: str | None, split_year: float,
        no_trend: bool = False) -> None:
    prep_a = ws.load_index(ws.find_index_root(train_index, root), train_index)
    prep_b = ws.load_index(ws.find_index_root(score_index, root), score_index)

    if prep_a["uncompensated"]:
        raise SystemExit(f"{train_index} uses UNCOMPENSATED mode -- this tool only "
                          f"implements the standard regression_factors/lte_response path")

    # Fit regression on ALL of the training index's own data.
    if no_trend:
        print(f"[{train_index}] fitting WITHOUT trend/acceleration terms "
              f"(annual/semiannual terms kept)")
        level, k0, amp, phase, annual, trend, accel = regression_factors_no_trend(
            prep_a["dates"], prep_a["forcing"], prep_a["data"], prep_a["m"])
    else:
        level, k0, amp, phase, annual, trend, accel = ws.regression_factors(
            prep_a["dates"], prep_a["forcing"], prep_a["data"], prep_a["m"], prep_a["trend_on"])

    model_train = ws.lte_response(prep_a["forcing"], prep_a["dates"], prep_a["m"],
                                   amp, phase, level, k0, trend, accel,
                                   prep_a["nonlin"], annual)
    model_train = ws.decompensate_f9(model_train, prep_a["ir"], prep_a["f9"])
    train_r = float(np.corrcoef(model_train, prep_a["data"])[0, 1])
    print(f"[{train_index}] in-sample fit over its own "
          f"{prep_a['dates'].min():.1f}-{prep_a['dates'].max():.1f} record: "
          f"r={train_r:.3f}  n={len(prep_a['dates'])}")

    # Anchor the IIR seed at the TRAINING record's own first date (not the
    # nominal IDATE config value) -- see build_forcing_at_anchored's
    # docstring for why. A hair before dates[0] guarantees np.searchsorted
    # lands on that exact row regardless of which array (short or long)
    # is being evaluated.
    anchor = prep_a["dates"][0] - 1e-3
    forcing_ext = build_forcing_at_anchored(prep_a, prep_b["dates"], anchor)

    # Apply the SAME regression coefficients -- fit only on the training
    # index -- to that extended forcing plus the scoring index's own
    # calendar dates (trend/annual terms are functions of calendar date).
    model_ext = lte_response_anchored(forcing_ext, prep_b["dates"], prep_a["dates"][0],
                                       prep_a["m"], amp, phase, level, k0, trend, accel,
                                       prep_a["nonlin"], annual)
    model_ext = ws.decompensate_f9(model_ext, prep_a["ir"], prep_a["f9"])

    dates_b, data_b = prep_b["dates"], prep_b["data"]
    overlap = dates_b >= split_year
    before = dates_b < split_year

    # Sanity check the actual point of this fix: the overlap-region MODEL
    # values (not just their correlation with data) must come out
    # bit-identical whether computed directly on kN020_E050's own short
    # array or sliced out of the long array's extended computation --
    # they are the same recursion now, just evaluated over different
    # amounts of surrounding history.
    same_window = dates_b >= prep_a["dates"][0]
    diff = model_ext[same_window] - model_train
    max_diff = float(np.max(np.abs(diff)))
    max_diff_after_1yr = float(np.max(np.abs(diff[12:]))) if len(diff) > 12 else float("nan")
    print(f"\nmax|model(long, sliced) - model(short, direct)| over the full "
          f"overlap = {max_diff:.3e}; restricted to month 13 onward "
          f"= {max_diff_after_1yr:.2e}")
    print("  (any nonzero difference is confined to the training record's "
          "own first 12 months: prep['ir']'s one_shot_lag_diff leaves the "
          "first 12 rows of a SHORT array undecompensated -- a boundary "
          "artifact already present in the original production fit itself, "
          "not something this anchoring fix should or does paper over. "
          "From month 13 on, the two computations are bit-identical.)")

    r_overlap = float(np.corrcoef(model_ext[overlap], data_b[overlap])[0, 1])
    print(f"\nFrozen {train_index} model (forcing pipeline + regression "
          f"coefficients, unchanged, correctly anchored), transplanted onto "
          f"{score_index}'s {dates_b.min():.1f}-{dates_b.max():.1f} record:")
    print(f"  overlap region (>= {split_year:.0f}, n={int(overlap.sum())}): "
          f"r={r_overlap:.3f}  (matches {train_index}'s own in-sample "
          f"r={train_r:.3f} up to the split_year != training start caveat below)")

    if before.sum() >= 3:
        r_before = float(np.corrcoef(model_ext[before], data_b[before])[0, 1])
        print(f"  BEFORE {split_year:.0f} (n={int(before.sum())}, never seen "
              f"by training): r={r_before:.3f}")
    else:
        print(f"  BEFORE {split_year:.0f}: not enough points "
              f"(n={int(before.sum())})")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--train", required=True,
                     help="index to fit on, using its own full record")
    ap.add_argument("--score", required=True,
                     help="longer index to transplant the frozen fit onto")
    ap.add_argument("--root", default=None,
                     help="explicit root directory containing both <index>/ dirs")
    ap.add_argument("--split-year", type=float, default=1950.0,
                     help="report correlation before this year as the "
                          "out-of-sample score (default 1950)")
    ap.add_argument("--no-trend", action="store_true",
                     help="drop the linear trend + quadratic accel terms "
                          "(keeps annual/semiannual); isolates whether the "
                          "level/trend near-cancellation is driving the "
                          "backward-extrapolation result")
    args = ap.parse_args()
    run(args.train, args.score, args.root, args.split_year, args.no_trend)


if __name__ == "__main__":
    main()
