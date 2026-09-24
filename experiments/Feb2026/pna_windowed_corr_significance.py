#!/usr/bin/env python3
"""pna_windowed_corr_significance.py -- two independent, rigorous checks on
the pna EXCLUDE=1970-1985 fit's "solid 50-month windowed correlation"
(pnasite1970-1985.png, bottom-left panel):

(1) AR(1)-surrogate significance test on the windowed-correlation trace
    itself: is the observed level of the real Model-vs-Data rolling
    correlation (outside the excluded 1970-1985 gap) explainable by
    chance autocorrelated noise, or is it a genuine, non-chance level of
    shared structure? Same discipline as winding_rank.py's own AR(1)-
    floor significance test, applied to a windowed-correlation trace
    instead of a ridge scan.

(2) Nested-model AIC/BIC comparison: does the fitted 6-winding LTE
    manifold (ltep + harm, i.e. 2 free base windings + 4 integer-harmonic
    multiples of them) explain pna's data better than a plain
    trend+annual+semiannual baseline with ZERO windings, once the
    parameter-count penalty is paid honestly -- including an
    effective-sample-size correction for residual autocorrelation
    (raw AIC/BIC on N=917 monthly points massively overstates
    confidence for a highly autocorrelated climate series).

Both models are refit via the SAME OLS machinery (ws.regression_factors)
on the SAME EXCLUDE=TRUE training subset (both outer flanks, 1970-1985
excluded) that the real Ada optimizer search actually used -- this is a
clean, honest apples-to-apples nested comparison, not a re-use of the
Ada-fitted amp/phase (which came from a stochastic search, not a plain
OLS solve restricted to this exact training set).
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

import ws

HERE = Path(__file__).resolve().parent
INDEX = "pna"
EXCLUDE_LO, EXCLUDE_HI = 1970.0, 1985.0
WINDOW = 50  # months, matches plot.py's own rolling-correlation window
N_SURROGATES = 2000
RNG_SEED = 0


def ar1_surrogate(x: np.ndarray, phi: float, sigma: float, rng: np.random.Generator) -> np.ndarray:
    """Same construction as build_k_sst.ar1_surrogate: a zero-mean AR(1)
    process matched to x's own lag-1 autocorrelation and standard
    deviation."""
    n = len(x)
    s = np.empty(n)
    s[0] = rng.normal(0.0, x.std())
    innov_sigma = sigma * math.sqrt(max(1e-12, 1.0 - phi ** 2))
    for i in range(1, n):
        s[i] = phi * s[i - 1] + rng.normal(0.0, innov_sigma)
    return s


def rolling_corr(a: np.ndarray, b: np.ndarray, window: int) -> np.ndarray:
    n = len(a) - window + 1
    out = np.empty(n)
    for i in range(n):
        out[i] = np.corrcoef(a[i:i + window], b[i:i + window])[0, 1]
    return out


def aic_bic(rss: float, n: int, k: int) -> tuple[float, float]:
    # classic Gaussian-likelihood OLS AIC/BIC (the additive -n/2*ln(2*pi)-n/2
    # constant is dropped -- it's identical for both models being compared,
    # so it cancels in every delta reported below)
    ll_term = n * math.log(rss / n)
    return ll_term + 2 * k, ll_term + k * math.log(n)


def main() -> None:
    root = ws.find_index_root(INDEX, None)
    prep = ws.load_index(root, INDEX)
    dates, data, forcing, m = prep["dates"], prep["data"], prep["forcing"], prep["m"]
    trend_on = True  # this run's own config (TREND=TRUE)

    excl_mask = (dates >= EXCLUDE_LO) & (dates < EXCLUDE_HI)
    train_mask = ~excl_mask

    print(f"=== {INDEX}: N={len(dates)} months, {dates.min():.2f}-{dates.max():.2f}, "
          f"training excludes [{EXCLUDE_LO},{EXCLUDE_HI}) -> "
          f"{train_mask.sum()} training / {excl_mask.sum()} held-out months ===\n")

    # -----------------------------------------------------------------
    # Part 1: AR(1)-surrogate significance test on the windowed
    # correlation trace, against the REAL fitted Ada Model (not a
    # refit) -- this is testing the actual saved lt.exe.p, exactly what
    # pnasite1970-1985.png plotted.
    # -----------------------------------------------------------------
    import pandas as pd
    df = pd.read_csv(root / INDEX / "lte_results.csv", header=None)
    model_col = df.iloc[:, 1].to_numpy(dtype=float)
    data_col = df.iloc[:, 2].to_numpy(dtype=float)
    time_col = df.iloc[:, 0].to_numpy(dtype=float)

    real_corrs = rolling_corr(model_col, data_col, WINDOW)
    corr_time = time_col[WINDOW - 1:]
    eval_mask = (corr_time < EXCLUDE_LO) | (corr_time >= EXCLUDE_HI)
    real_mean = float(np.mean(real_corrs[eval_mask]))
    real_min = float(np.min(real_corrs[eval_mask]))

    resid = data_col - data_col.mean()
    phi = float(np.corrcoef(resid[:-1], resid[1:])[0, 1])
    sigma = float(data_col.std())
    rng = np.random.default_rng(RNG_SEED)
    surrogate_means = np.empty(N_SURROGATES)
    surrogate_mins = np.empty(N_SURROGATES)
    for i in range(N_SURROGATES):
        surr = ar1_surrogate(data_col, phi, sigma, rng)
        sc = rolling_corr(model_col, surr, WINDOW)
        surrogate_means[i] = np.mean(sc[eval_mask])
        surrogate_mins[i] = np.min(sc[eval_mask])

    z = (real_mean - surrogate_means.mean()) / surrogate_means.std()
    percentile = float(np.mean(surrogate_means < real_mean)) * 100.0
    z_min = (real_min - surrogate_mins.mean()) / surrogate_mins.std()
    percentile_min = float(np.mean(surrogate_mins < real_min)) * 100.0

    print("--- Part 1: windowed-correlation significance (AR(1) null) ---")
    print(f"  real Data lag-1 autocorrelation (phi) used for surrogates: {phi:.4f}")
    print(f"  MEAN windowed correlation, real Model vs real Data "
          f"(outside excluded gap): {real_mean:.4f}")
    print(f"  {N_SURROGATES} AR(1) surrogates, mean-of-window stat: "
          f"mean={surrogate_means.mean():.4f} std={surrogate_means.std():.4f}")
    print(f"  z-score (mean stat): {z:.2f}   percentile: {percentile:.2f}%")
    print(f"  WORST-CASE (minimum) windowed correlation outside the gap: {real_min:.4f}")
    print(f"  {N_SURROGATES} AR(1) surrogates, min-of-window stat: "
          f"mean={surrogate_mins.mean():.4f} std={surrogate_mins.std():.4f}")
    print(f"  z-score (min stat): {z_min:.2f}   percentile: {percentile_min:.2f}%\n")

    # -----------------------------------------------------------------
    # Part 2: nested AIC/BIC, full (6-winding) vs baseline (0-winding),
    # both refit via plain OLS on the SAME EXCLUDE=1970-1985 training
    # subset.
    # -----------------------------------------------------------------
    t_tr, fv_tr, y_tr = dates[train_mask], forcing[train_mask], data[train_mask]

    A_full = ws.build_design(t_tr, fv_tr, m, trend_on)
    coef_full, res_full, *_ = np.linalg.lstsq(A_full, y_tr, rcond=None)
    rss_full = float(np.sum((y_tr - A_full @ coef_full) ** 2))
    k_full_linear = A_full.shape[1]          # level,k0 + 12 sin/cos + 6 trend/seasonal
    k_full_nonlinear = 2                     # ltep[0], ltep[1] (the 2 FREE base windings)
    k_full = k_full_linear + k_full_nonlinear

    A_base = ws.build_design(t_tr, fv_tr, np.array([]), trend_on)
    coef_base, *_ = np.linalg.lstsq(A_base, y_tr, rcond=None)
    rss_base = float(np.sum((y_tr - A_base @ coef_base) ** 2))
    k_base = A_base.shape[1]                 # level,k0 + 6 trend/seasonal, zero windings

    n_raw = len(y_tr)
    resid_base = y_tr - A_base @ coef_base
    phi_r = float(np.corrcoef(resid_base[:-1], resid_base[1:])[0, 1])
    phi_r = max(-0.98, min(0.98, phi_r))
    n_eff = max(k_full + 2, int(round(n_raw * (1.0 - phi_r) / (1.0 + phi_r))))

    aic_full_raw, bic_full_raw = aic_bic(rss_full, n_raw, k_full)
    aic_base_raw, bic_base_raw = aic_bic(rss_base, n_raw, k_base)
    aic_full_eff, bic_full_eff = aic_bic(rss_full, n_eff, k_full)
    aic_base_eff, bic_base_eff = aic_bic(rss_base, n_eff, k_base)

    print("--- Part 2: nested AIC/BIC, full (6-winding) vs baseline (trend+seasonal only) ---")
    print(f"  training subset (both outer flanks, gap excluded): n_raw={n_raw} months")
    print(f"  baseline residual lag-1 autocorrelation: {phi_r:.4f}  ->  n_eff={n_eff}")
    print(f"  full model:     k={k_full} ({k_full_linear} OLS + {k_full_nonlinear} nonlinear winding) "
          f"  RSS={rss_full:.4f}")
    print(f"  baseline model: k={k_base}   RSS={rss_base:.4f}   "
          f"(RSS ratio full/baseline = {rss_full / rss_base:.4f})\n")
    print(f"  {'':14s}{'raw N='+str(n_raw):>16s}{'eff N='+str(n_eff):>16s}")
    print(f"  {'delta AIC':14s}{aic_base_raw - aic_full_raw:16.2f}{aic_base_eff - aic_full_eff:16.2f}"
          f"   (positive = full model preferred; >10 = decisive)")
    print(f"  {'delta BIC':14s}{bic_base_raw - bic_full_raw:16.2f}{bic_base_eff - bic_full_eff:16.2f}"
          f"   (positive = full model preferred; >10 = decisive)")
    print(f"\n  NOTE: 4 of the full model's 6 windings are integer-harmonic multiples "
          f"of ltep[1] (harm={prep['params'].get('harm')}), a discrete/combinatorial "
          f"choice not captured by k_full's continuous-parameter count -- flagged, "
          f"not silently absorbed into k.")


if __name__ == "__main__":
    main()
