#!/usr/bin/env python3
"""pysr_optimizer.py -- standalone PySR symbolic-regression fit for an
lte_results.csv file.

Usage:
    python3 pysr_optimizer.py path/to/lte_results.csv

lte_results.csv has no header, columns (1-indexed, as the file is usually
described):
    1: time
    2: model   (the existing GEM-LTE regression's own reconstruction)
    3: data    (the real observed series)
    4: manifold (the shared non-autonomous forcing term)
    5+: not used here

This script fits Column 3 (data) using ONLY Column 4 (manifold) and
Column 1 (time) as input variates -- it never shows PySR the existing
model's own output (Column 2) or any column past 4. Column 2 is used only
afterwards, as a final, uninvolved comparison.

WHY TIME AND MANIFOLD ARE NOT BOTH HANDED TO PySR THE SAME WAY:
Handing PySR raw (t, manifold) as two interchangeable inputs to search
over freely lets it fall into a cheap AUTONOMOUS local optimum --
something like cos(t*omega) -- because a bare-time sinusoid is a "free"
low-complexity feature from a pure curve-fitting standpoint, even though
it has nothing to do with the real lunisolar forcing. Once evolution
finds that local optimum it rarely escapes it, since any mutation toward
an actual multi-term manifold-based winding sum looks worse mid-mutation
(a fitness-valley problem), regardless of how long the search runs.

This mirrors the real GEM-LTE model structure itself: in the Ada
regression, calendar time enters ONLY linearly/quadratically (the
"trend"/"accel" terms), and every oscillatory ("winding") term is
sin(M * manifold + phase) -- time never appears directly inside a trig
function. So this script enforces that same separation instead of
leaving it to chance:
  1. Fit a polynomial in time (degree --trend-degree, default 1) to the
     data by ordinary least squares -- this uses Column 1 exactly the
     way the real model uses it, as a trend/accel term outside any
     oscillation.
  2. Subtract that trend to get a residual.
  3. Run PySR on the residual using Column 4 (manifold) as the ONLY
     input variable, so it cannot fall back on bare time and is forced
     to express the oscillatory structure as a function of the shared
     manifold -- the actual non-autonomous winding search.
  4. Recombine trend + PySR(manifold) for the final comparison against
     Column 2 and Column 3.
"""
from __future__ import annotations

import argparse

import numpy as np


def load_columns(path: str):
    raw = np.loadtxt(path, delimiter=",")
    time = raw[:, 0]
    model = raw[:, 1]
    data = raw[:, 2]
    manifold = raw[:, 3]
    return time, model, data, manifold


def corr(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.corrcoef(a, b)[0, 1])


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv_path", help="path to an lte_results.csv file")
    ap.add_argument("--trend-degree", type=int, default=1,
                     help="degree of the polynomial-in-time trend removed "
                          "before the manifold search (default: 1, i.e. "
                          "linear trend; use 2 to also include an accel term)")
    ap.add_argument("--niterations", type=int, default=500,
                     help="PySR search iterations (default: 500 -- this "
                          "search needs to run long to find multi-term "
                          "winding sums, not just one sinusoid)")
    ap.add_argument("--populations", type=int, default=60,
                     help="PySR populations (default: 60)")
    ap.add_argument("--population-size", type=int, default=50,
                     help="individuals per population (default: 50)")
    ap.add_argument("--maxsize", type=int, default=45,
                     help="max equation complexity (default: 45 -- large "
                          "enough to hold several Amp*sin(M*manifold+phase) "
                          "terms added together)")
    ap.add_argument("--procs", type=int, default=4,
                     help="worker processes for PySR (default: 4)")
    args = ap.parse_args()

    time, model, data, manifold = load_columns(args.csv_path)

    # --- Step 1-2: remove the time trend by ordinary least squares -----
    t0 = time[0]
    t_c = time - t0
    trend_powers = np.column_stack([t_c ** k for k in range(args.trend_degree + 1)])
    trend_coeffs, *_ = np.linalg.lstsq(trend_powers, data, rcond=None)
    trend_fit = trend_powers @ trend_coeffs
    residual = data - trend_fit

    print(f"[pysr_optimizer] {len(data)} rows from {args.csv_path}")
    print(f"[pysr_optimizer] degree-{args.trend_degree} time trend coeffs "
          f"(constant, linear, ...): {trend_coeffs}")
    print(f"[pysr_optimizer] trend explains R={corr(trend_fit, data):+.4f} "
          f"of Column 3 (data) on its own; residual std "
          f"{residual.std():.4f} (data std {data.std():.4f})")

    # --- Step 3: symbolic-regress the residual using manifold ONLY -----
    X = manifold.reshape(-1, 1)
    y = residual

    from pysr import PySRRegressor

    sr = PySRRegressor(
        niterations=args.niterations,
        populations=args.populations,
        population_size=args.population_size,
        procs=args.procs,
        binary_operators=["+", "-", "*", "/"],
        unary_operators=["sin", "cos", "square"],
        nested_constraints={
            "sin": {"sin": 0, "cos": 0},
            "cos": {"sin": 0, "cos": 0},
        },
        maxsize=args.maxsize,
        model_selection="accuracy",  # lowest-loss candidate, not PySR's own
                                      # complexity-penalized "best" pick --
                                      # a genuine multi-term winding sum is
                                      # legitimately more complex than a
                                      # single sinusoid, so let the data
                                      # (via the Pareto front below) decide,
                                      # not PySR's parsimony heuristic
        progress=True,
    )

    print("\n[pysr_optimizer] fitting the detrended residual from Column 4 "
          "(manifold) only -- this is the non-autonomous winding search")
    sr.fit(X, y, variable_names=["manifold"])

    # --- Step 4: report the WHOLE Pareto front, not just one pick -------
    # PySR's internal "best" trades loss against complexity using a
    # heuristic that has no notion of "does this look like a physical
    # winding sum" -- it happily prefers a single bare cos(k*manifold)
    # over a multi-term sum that fits far better. Report every candidate's
    # actual R against the real data/model instead, so the winding
    # structure (if present) is visible directly rather than hidden behind
    # that heuristic.
    front = sr.equations_.copy()
    print("\n[pysr_optimizer] full Pareto front (all candidates found), "
          "recombined with the time trend and scored against real data:")
    print(f"  {'complexity':>10s}  {'loss':>10s}  {'R vs data':>10s}  "
          f"{'R vs model':>10s}  equation")
    best_row = None
    best_r = -2.0
    for i, row in front.iterrows():
        wp = sr.predict(X, index=i)
        full_pred = trend_fit + wp
        r_data = corr(full_pred, data)
        r_model = corr(full_pred, model)
        print(f"  {row['complexity']:>10d}  {row['loss']:>10.4g}  "
              f"{r_data:>+10.4f}  {r_model:>+10.4f}  {row['equation']}")
        if r_data > best_r:
            best_r = r_data
            best_row = i

    print(f"\n[pysr_optimizer] best winding equation by ACTUAL R vs data "
          f"(index {best_row} in the Pareto front):")
    print(sr.sympy(best_row))

    winding_pred = sr.predict(X, index=best_row)
    pred = trend_fit + winding_pred

    r_fit = corr(pred, data)
    r_windings_only = corr(winding_pred, residual)
    r_model_data = corr(model, data)
    r_pred_model = corr(pred, model)

    print("\n[pysr_optimizer] ---- final comparison (Column 2 not used in the fit) ----")
    print(f"  windings-only vs detrended residual : R = {r_windings_only:+.4f}")
    print(f"  trend + SR windings vs Column 3 (data)  : R = {r_fit:+.4f}")
    print(f"  Column 2  vs Column 3 (data)  : R = {r_model_data:+.4f}  (existing LTE model, for reference)")
    print(f"  trend + SR windings vs Column 2 (model) : R = {r_pred_model:+.4f}")


if __name__ == "__main__":
    main()
