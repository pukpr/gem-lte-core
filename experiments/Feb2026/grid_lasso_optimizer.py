#!/usr/bin/env python3
"""grid_lasso_optimizer.py -- standalone grid + LASSO winding-number search
for an lte_results.csv file, validated with MANY chronological
train/test folds (rolling walk-forward and/or blocked contiguous k-fold),
not a single train/validation/test split.

Usage:
    python3 grid_lasso_optimizer.py path/to/lte_results.csv
    python3 grid_lasso_optimizer.py path/to/lte_results.csv --mode blocked --n-blocks 8
    python3 grid_lasso_optimizer.py path/to/lte_results.csv --plot

Same column convention and same input-variate restriction throughout this
project's other optimizer scripts (columns 1-indexed, no header):
    1: time
    2: model    (existing GEM-LTE regression's own reconstruction)
    3: data     (the real observed series -- what we fit)
    4: manifold (the shared non-autonomous forcing term)
    5+: not used

WHY GRID + LASSO INSTEAD OF FREE-FORM SYMBOLIC REGRESSION:
pysr_optimizer.py showed that letting a genetic SR search invent the
winding multiplier M via continuous constant optimization gets stuck
composing several small-M sinusoids algebraically rather than proposing
one clean high-M term directly. The fix here is the one already
validated elsewhere in this project (sindy_discovery.py,
cv_ridge_transfer.py): precompute a dense grid of
sin(M*manifold)/cos(M*manifold) and let a LINEAR sparse solver (LASSO)
pick which survive -- "select columns from a fixed dictionary" instead
of "invent a frequency by nonlinear constant search".

WHY MANY FOLDS, NOT ONE SPLIT (this is the point of this version):
A single chronological train/validation/test split answers "did this one
boundary happen to overfit" but says nothing about whether that pattern
holds generally, nor whether a different boundary would look better or
worse purely by chance. This project's own cv_rolling_blocked.py makes
exactly this argument for the resolved-PDE forward model; the same
argument applies here unchanged, so this script reuses its two disciplined
fold schemes directly:

  rolling (walk-forward): strictly causal. Train only on the past (an
    expanding or fixed-length sliding window), score on the immediately
    -following block, advance, repeat. This is what "does this model
    forecast" actually means for a time series.

  blocked (contiguous k-fold): split the record into K contiguous
    chunks; for each, train on all the OTHER chunks (past and future)
    and score on the held-out one, with an embargo buffer dropped from
    training on either side of the test block. Not a forecast test, but
    a fairer generalization test than one hand-picked split, and the
    standard alternative to shuffled k-fold for autocorrelated series.

Every fold repeats the FULL pipeline independently -- trend removal,
coarse LASSO dictionary selection, refine, final Amplitude/Phase refit --
using ONLY that fold's training rows, then scores strictly on that fold's
held-out rows. Fold correlations are then summarized (mean/median/sd/
min/max), the same convention cv_rolling_blocked.py already reports in.
A model that only "works" on one lucky split will show up here as a wide
spread across folds instead of one deceptively good number.
"""
from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent


def load_columns(path: str):
    raw = np.loadtxt(path, delimiter=",")
    return raw[:, 0], raw[:, 1], raw[:, 2], raw[:, 3]  # time, model, data, manifold


def corr(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 2 or a.std() == 0.0 or b.std() == 0.0:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


def fit_at_M(manifold: np.ndarray, target: np.ndarray, M: float) -> float:
    """Single-mode OLS R: target ~ c0 + c1*sin(M*manifold) + c2*cos(M*manifold)."""
    design = np.column_stack([np.ones_like(manifold),
                               np.sin(M * manifold), np.cos(M * manifold)])
    coef, *_ = np.linalg.lstsq(design, target, rcond=None)
    return corr(design @ coef, target)


def refine(manifold: np.ndarray, target: np.ndarray, m0: float, dm: float,
           window_mult: float, refine_dm: float) -> float:
    """Local fine-grid scan around a coarse survivor m0 -- confined to a
    small window so it never reintroduces the collinearity a global fine
    grid would have."""
    lo, hi = max(refine_dm, m0 - window_mult * dm), m0 + window_mult * dm
    grid = np.arange(lo, hi + refine_dm / 2, refine_dm)
    scores = [fit_at_M(manifold, target, m) ** 2 for m in grid]
    return float(grid[int(np.argmax(scores))])


def winding_design(man: np.ndarray, refined_m: list[float]) -> np.ndarray:
    cols = [np.ones_like(man)]
    for m in refined_m:
        cols.append(np.sin(m * man)); cols.append(np.cos(m * man))
    return np.column_stack(cols)


# ---------------------------------------------------------------------------
# Harmonic-family scan: a Doodson-style periodogram over candidate
# FUNDAMENTAL frequencies M0, jointly fitting a whole family of harmonics
# h*M0 (h=1..max_harmonic) at once, rather than one M at a time.
#
# WHY THIS IS NEEDED (found by testing on real Baltic data): Baltic's own
# real windings K3..K7 are EXACT integer multiples (1x,6x,9x,10x,11x) of
# one backbone M0~0.2076 -- a genuine harmonic series, not independent
# tones. Verified directly: each harmonic's OWN marginal correlation with
# the detrended residual is only R~0.03-0.09, but fitting all of them
# TOGETHER reaches R~0.30 (matching the real model). The per-M LASSO
# dictionary in fit_windings() evaluates each grid column's marginal
# contribution independently, so a family like this -- individually weak,
# only strong jointly, phase-locked at simple integer ratios -- is exactly
# the structure L1 selection is known to miss. Scanning candidate
# fundamentals and scoring the JOINT fit of their whole harmonic family
# sidesteps that: it turns "is M any good alone" into "is this WHOLE
# non-sinusoidal periodic family any good", which is the right question
# for a Doodson-style response.
# ---------------------------------------------------------------------------

def family_design(man: np.ndarray, m0: float, max_harmonic: int) -> np.ndarray:
    cols = [np.ones_like(man)]
    for h in range(1, max_harmonic + 1):
        cols.append(np.sin(h * m0 * man)); cols.append(np.cos(h * m0 * man))
    return np.column_stack(cols)


def family_r2(man: np.ndarray, target: np.ndarray, m0: float, max_harmonic: int) -> float:
    design = family_design(man, m0, max_harmonic)
    coef, *_ = np.linalg.lstsq(design, target, rcond=None)
    return corr(design @ coef, target) ** 2


def harmonic_family_scan(man: np.ndarray, residual: np.ndarray, m_grid: np.ndarray,
                         max_harmonic: int, refine_dm: float) -> float | None:
    """Coarse scan over m_grid for the best-supported fundamental M0 (fixed
    harmonic count at every candidate, so R^2 is comparable across M0),
    then a local fine-grid refine around the winner scoring the SAME joint
    family fit (not a single-mode score)."""
    scores = [family_r2(man, residual, m0, max_harmonic) for m0 in m_grid]
    m0_coarse = float(m_grid[int(np.argmax(scores))])
    dm = float(m_grid[1] - m_grid[0]) if len(m_grid) > 1 else refine_dm
    lo, hi = max(refine_dm, m0_coarse - dm), m0_coarse + dm
    fine_grid = np.arange(lo, hi + refine_dm / 2, refine_dm)
    fine_scores = [family_r2(man, residual, m0, max_harmonic) for m0 in fine_grid]
    return float(fine_grid[int(np.argmax(fine_scores))])


def sparsify_family(man: np.ndarray, residual: np.ndarray, m0: float,
                    max_harmonic: int, keep_frac: float) -> list[float]:
    """Fit the full H-harmonic family once, then keep only the harmonics
    whose amplitude is at least keep_frac of the family's largest -- a
    magnitude-based parsimony pass playing LASSO's role WITHIN the family,
    since LASSO across the whole grid can't select this family at all."""
    design = family_design(man, m0, max_harmonic)
    coef, *_ = np.linalg.lstsq(design, residual, rcond=None)
    amps = [float(np.hypot(coef[1 + 2 * h], coef[2 + 2 * h])) for h in range(max_harmonic)]
    if max(amps) <= 0.0:
        return []
    thresh = keep_frac * max(amps)
    return [(h + 1) * m0 for h, a in enumerate(amps) if a >= thresh]


# ---------------------------------------------------------------------------
# Fit trend + harmonic-family scan + LASSO dictionary + refine + final
# refit on a given set of rows -- the one routine every fold, and the
# final full-data fit, both call. Returns everything needed to EVALUATE
# the fitted model on any (t, manifold) pair, not just the rows it was
# fit on.
# ---------------------------------------------------------------------------

def fit_windings(t_tr: np.ndarray, d_tr: np.ndarray, man_tr: np.ndarray,
                 args) -> dict:
    t0 = t_tr[0]
    trend_powers_tr = np.column_stack([(t_tr - t0) ** k for k in range(args.trend_degree + 1)])
    trend_coeffs, *_ = np.linalg.lstsq(trend_powers_tr, d_tr, rcond=None)
    residual_tr = d_tr - trend_powers_tr @ trend_coeffs

    m_grid = np.arange(args.dm, args.m_max + args.dm / 2, args.dm)

    # --- Stage 1: harmonic-family scan (catches phase-locked families a
    # per-M LASSO dictionary structurally can't select) ---------------------
    family_m: list[float] = []
    if args.max_harmonic > 0:
        m0_star = harmonic_family_scan(man_tr, residual_tr, m_grid,
                                       args.max_harmonic, args.refine_dm)
        family_m = sparsify_family(man_tr, residual_tr, m0_star,
                                   args.max_harmonic, args.harmonic_keep_frac)

    # --- Stage 2: per-M LASSO dictionary, for independent tones outside
    # the family (e.g. Baltic's K1/K2, which are not clean multiples of
    # its own backbone) -------------------------------------------------
    cols, m_of_col = [], []
    for m in m_grid:
        cols.append(np.sin(m * man_tr)); m_of_col.append(m)
        cols.append(np.cos(m * man_tr)); m_of_col.append(m)
    Theta_tr = np.column_stack(cols)

    from sklearn.linear_model import LassoCV
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    inner_cv = max(2, min(args.inner_cv, len(d_tr) // 50))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pipe = make_pipeline(
            StandardScaler(),
            LassoCV(cv=inner_cv, max_iter=200_000, n_alphas=100),
        )
        pipe.fit(Theta_tr, residual_tr)
    nonzero = np.flatnonzero(pipe.named_steps["lassocv"].coef_)
    survivor_m = sorted({round(m_of_col[i], 6) for i in nonzero})
    lasso_refined_m = sorted(refine(man_tr, residual_tr, m0, args.dm,
                                    args.refine_window, args.refine_dm)
                             for m0 in survivor_m)

    # --- Combine both stages, refit ONCE jointly on the original residual --
    # Merge tolerance is deliberately NOT tied to --dm: shrinking it along
    # with a fine scan grid (as `dm/2` used to) lets near-duplicate,
    # highly collinear columns back into the final design matrix, which
    # blows up the unregularized refit's coefficients (verified directly:
    # dm=0.02/0.01 on a clean synthetic Baltic target gave "amplitudes" in
    # the millions/billions with near-identical adjacent M's, despite R
    # staying ~0.999 -- classic multicollinearity, not a real fit).
    # --min-separation is the actual "these are the same winding" distance,
    # independent of how finely M was scanned to find it.
    raw_combined = sorted(family_m + lasso_refined_m)
    refined_m: list[float] = []
    for m in raw_combined:
        if refined_m and abs(m - refined_m[-1]) < args.min_separation:
            continue
        refined_m.append(m)

    winding_coef = np.zeros(1)
    if refined_m:
        Design_tr = winding_design(man_tr, refined_m)
        winding_coef, *_ = np.linalg.lstsq(Design_tr, residual_tr, rcond=None)

    return dict(t0=t0, trend_degree=args.trend_degree, trend_coeffs=trend_coeffs,
                refined_m=refined_m, winding_coef=winding_coef,
                n_train=len(d_tr), train_span=(float(t_tr.min()), float(t_tr.max())))


def eval_fit(fit: dict, t: np.ndarray, man: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Evaluate a fit_windings() result on arbitrary (t, manifold) rows.
    Returns (full prediction, trend-only, windings-only)."""
    powers = np.column_stack([(t - fit["t0"]) ** k for k in range(fit["trend_degree"] + 1)])
    trend = powers @ fit["trend_coeffs"]
    winding = (winding_design(man, fit["refined_m"]) @ fit["winding_coef"]
               if fit["refined_m"] else np.zeros_like(t))
    return trend + winding, trend, winding


def winding_table(fit: dict) -> list[tuple[float, float, float]]:
    """[(M, Amplitude, Phase), ...] in the same convention as the real Ada
    model / k_amp_phase.json."""
    rows = []
    for i, m in enumerate(fit["refined_m"]):
        c_sin, c_cos = fit["winding_coef"][1 + 2 * i], fit["winding_coef"][2 + 2 * i]
        rows.append((m, float(np.hypot(c_sin, c_cos)), float(np.arctan2(c_sin, c_cos))))
    return rows


# ---------------------------------------------------------------------------
# One fold: fit on TRAIN rows only (via fit_windings), score strictly on
# TEST rows.
# ---------------------------------------------------------------------------

def fit_fold(t_tr, mo_tr, d_tr, man_tr, t_te, mo_te, d_te, man_te,
             args) -> dict | None:
    m_grid = np.arange(args.dm, args.m_max + args.dm / 2, args.dm)
    n_dict_cols = 2 * len(m_grid)
    min_train = max(50, 4 * n_dict_cols // 2)  # a few points per free coef
    if len(d_tr) < min_train or len(d_te) < 3:
        return None

    fit = fit_windings(t_tr, d_tr, man_tr, args)
    pred_te, _, _ = eval_fit(fit, t_te, man_te)
    return dict(
        n_train=len(d_tr), n_test=len(d_te),
        train_span=fit["train_span"], test_span=(float(t_te.min()), float(t_te.max())),
        n_survivors=len(fit["refined_m"]), refined_m=fit["refined_m"],
        cc_data=corr(pred_te, d_te),
        cc_model=corr(pred_te, mo_te),
        cc_model_data=corr(mo_te, d_te),
    )


# ---------------------------------------------------------------------------
# Outer fold schemes (reused verbatim in spirit from cv_rolling_blocked.py)
# ---------------------------------------------------------------------------

def rolling_cv(time, model, data, manifold, args) -> list[dict]:
    t_min, t_max = float(time.min()), float(time.max())
    folds = []
    test_start = t_min + args.initial_years
    while test_start < t_max:
        test_end = min(test_start + args.horizon_years, t_max)
        if test_end - test_start < args.horizon_years * 0.5:
            break
        if args.window == "sliding":
            train_start = max(t_min, test_start - args.window_years)
        else:
            train_start = t_min
        train_idx = np.nonzero((time >= train_start) & (time < test_start))[0]
        test_idx = np.nonzero((time >= test_start) & (time < test_end))[0]
        result = fit_fold(time[train_idx], model[train_idx], data[train_idx], manifold[train_idx],
                          time[test_idx], model[test_idx], data[test_idx], manifold[test_idx], args)
        if result is not None:
            folds.append(result)
        test_start += args.step_years
    return folds


def blocked_cv(time, model, data, manifold, args) -> list[dict]:
    t_min, t_max = float(time.min()), float(time.max())
    edges = np.linspace(t_min, t_max, args.n_blocks + 1)
    folds = []
    for k in range(args.n_blocks):
        b_start, b_end = edges[k], edges[k + 1]
        test_idx = np.nonzero((time >= b_start) & (time < b_end))[0]
        train_idx = np.nonzero((time < b_start - args.embargo_years) |
                               (time >= b_end + args.embargo_years))[0]
        result = fit_fold(time[train_idx], model[train_idx], data[train_idx], manifold[train_idx],
                          time[test_idx], model[test_idx], data[test_idx], manifold[test_idx], args)
        if result is not None:
            result["block"] = k
            folds.append(result)
    return folds


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def summarize(vals: list[float]) -> str:
    if not vals:
        return "no scoreable folds"
    arr = np.array(vals)
    return (f"mean={arr.mean():+.3f}  median={np.median(arr):+.3f}  "
            f"sd={arr.std():.3f}  min={arr.min():+.3f}  max={arr.max():+.3f}  n={len(arr)}")


def print_folds(label: str, folds: list[dict]) -> None:
    print(f"  {label}:")
    if not folds:
        print("    (no folds -- record too short, or --dm/--m-max too rich, for these settings)")
        return
    for i, f in enumerate(folds):
        tr = f"{f['train_span'][0]:.1f}-{f['train_span'][1]:.1f}"
        te = f"{f['test_span'][0]:.1f}-{f['test_span'][1]:.1f}"
        print(f"    fold {i:2d}  train={tr:>13} (n={f['n_train']:4d})  "
              f"test={te:>13} (n={f['n_test']:3d})  "
              f"survivors={f['n_survivors']:2d}  "
              f"cc(vs data)={f['cc_data']:+.3f}  cc(vs model)={f['cc_model']:+.3f}  "
              f"[model-vs-data ref={f['cc_model_data']:+.3f}]")
    print(f"    -> vs data : {summarize([f['cc_data'] for f in folds])}")
    print(f"    -> vs model: {summarize([f['cc_model'] for f in folds])}")


def plot_folds(csv_path: str, rolling: list[dict], blocked: list[dict], out_path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    panels = [(f, t, c) for f, t, c in
              [(rolling, "rolling (walk-forward)", "tab:blue"),
               (blocked, "blocked (contiguous k-fold)", "tab:orange")] if f]
    if not panels:
        return
    fig, axes = plt.subplots(1, len(panels), figsize=(6.5 * len(panels), 4.2), squeeze=False)
    for ax, (folds, title, color) in zip(axes[0], panels):
        ccs = [f["cc_data"] for f in folds]
        ref = [f["cc_model_data"] for f in folds]
        centers = [0.5 * (f["test_span"][0] + f["test_span"][1]) for f in folds]
        width = (centers[1] - centers[0]) * 0.8 if len(centers) > 1 else 5.0
        ax.axhline(0.0, color="0.6", linewidth=0.8)
        ax.bar(centers, ccs, width=width, color=color, alpha=0.85, label="SR windings vs data")
        ax.plot(centers, ref, "k.--", linewidth=1, label="existing model vs data (ref)")
        ax.axhline(np.mean(ccs), color="black", linestyle=":", linewidth=1,
                   label=f"mean={np.mean(ccs):+.3f}")
        ax.set_title(title, fontweight="bold", fontsize=10)
        ax.set_xlabel("test-window center (year)")
        ax.set_ylabel("held-out correlation vs Column 3 (data)")
        ax.set_ylim(-1.0, 1.0)
        ax.legend(fontsize=8)
    fig.suptitle(f"{Path(csv_path).parent.name}: grid+LASSO winding fit, "
                 f"held-out correlation per fold", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"  saved {out_path}")


def fold_stability(all_folds: list[dict], dm: float) -> list[tuple[float, int, int]]:
    """Bin every fold's surviving M's to the dictionary's own grid spacing
    and count how many DISTINCT folds each bin appears in -- an M that
    only ever shows up in the final full-data fit and never recurs across
    independent folds is fold-specific noise, not a robust winding."""
    from collections import Counter
    counts: Counter[float] = Counter()
    for f in all_folds:
        seen_bins = {round(m / dm) * dm for m in f["refined_m"]}
        for b in seen_bins:
            counts[b] += 1
    return sorted(((m, n, len(all_folds)) for m, n in counts.items()),
                  key=lambda x: (-x[1], x[0]))


def report_final_fit(csv_path: str, all_folds: list[dict], args,
                     target_label: str, reference_label: str) -> None:
    time, model, data, manifold = load_columns(csv_path)
    if args.target == "model":
        data, model = model, data  # fit Column 2, keep Column 3 as reference-only
    fit = fit_windings(time, data, manifold, args)
    pred, trend, winding = eval_fit(fit, time, manifold)

    print("\n[grid_lasso_optimizer] ==== FINAL FIT (all data, this is the "
          "deployable model) ====")
    print("  fit on the FULL record -- NOT held out, so its own R below is "
          "optimistic; trust the fold statistics above for the honest "
          "generalization estimate, this is only the equation itself.")
    print(f"  degree-{args.trend_degree} time trend coeffs "
          f"(constant, linear, ...): {fit['trend_coeffs']}")

    stability = fold_stability(all_folds, args.min_separation) if all_folds else []
    stable_bins = {m for m, n, total in stability if n >= max(2, total // 3)}

    print("\n  winding table (M, Amplitude, Phase) -- 'seen in N/T folds' "
          "says how often an independent fold's OWN search also surfaced "
          "a survivor near this M; low counts mean this term is likely "
          "specific to fitting the full record, not a robust winding:")
    print(f"  {'M':>10s}  {'Amplitude':>12s}  {'Phase(rad)':>12s}  {'seen in folds':>14s}")
    for m, amp, phase in winding_table(fit):
        b = round(m / args.min_separation) * args.min_separation
        n_seen = next((n for mb, n, total in stability if abs(mb - b) < 1e-9), 0)
        total = len(all_folds)
        flag = "" if not all_folds else (" <-- STABLE" if b in stable_bins else " <-- fold-specific")
        print(f"  {m:>10.5f}  {amp:>12.6f}  {phase:>12.6f}  "
              f"{f'{n_seen}/{total}':>14s}{flag}")
    if not fit["refined_m"]:
        print("  (no winding terms survived on the full record either)")

    r_windings_only = corr(winding, data - trend)
    r_fit = corr(pred, data)
    r_model_data = corr(model, data)
    r_pred_model = corr(pred, model)
    print(f"\n  in-sample (NOT held out): windings-only vs detrended residual R={r_windings_only:+.4f}")
    print(f"  in-sample (NOT held out): full fit vs {target_label:<22s} R={r_fit:+.4f}  <- fit target")
    print(f"  {reference_label} vs {target_label}                              R={r_model_data:+.4f}  (reference only, not fit)")
    print(f"  in-sample (NOT held out): full fit vs {reference_label:<22s} R={r_pred_model:+.4f}  (reference only)")


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv_path", help="path to an lte_results.csv file")
    ap.add_argument("--mode", choices=["rolling", "blocked", "both"], default="both")
    ap.add_argument("--initial-years", type=float, default=40.0,
                     help="rolling: length of the first training window (yr)")
    ap.add_argument("--step-years", type=float, default=10.0,
                     help="rolling: how far the origin advances per fold (yr)")
    ap.add_argument("--horizon-years", type=float, default=10.0,
                     help="rolling: length of each held-out test block (yr)")
    ap.add_argument("--window", choices=["expanding", "sliding"], default="expanding",
                     help="rolling: grow the training window (default) or keep "
                          "it a fixed --window-years length")
    ap.add_argument("--window-years", type=float, default=None,
                     help="rolling: training-window length when --window sliding "
                          "(default: --initial-years)")
    ap.add_argument("--n-blocks", type=int, default=6,
                     help="blocked: number of contiguous blocks")
    ap.add_argument("--embargo-years", type=float, default=1.0,
                     help="blocked: buffer dropped from training data on each "
                          "side of the held-out block")
    ap.add_argument("--trend-degree", type=int, default=1,
                     help="degree of the polynomial-in-time trend removed "
                          "before the manifold search, per fold (default: 1)")
    ap.add_argument("--m-max", type=float, default=5.0,
                     help="max winding number scanned (default 5.0)")
    ap.add_argument("--dm", type=float, default=0.03,
                     help="COARSE grid step for the LASSO dictionary and the "
                          "harmonic-family scan (default 0.03 -- verified on "
                          "a clean synthetic Baltic target: 0.1 fragments "
                          "the real K6=1.2454 winding into two nearby bins, "
                          "0.03 (with --min-separation decoupled from --dm) "
                          "recovers it almost exactly, diff 0.0006, at "
                          "roughly half the runtime cost of going to 0.02)")
    ap.add_argument("--refine-dm", type=float, default=0.001,
                     help="fine grid step for the local refine step (default 0.001)")
    ap.add_argument("--refine-window", type=float, default=2.0,
                     help="refine window half-width, in multiples of --dm (default 2.0)")
    ap.add_argument("--min-separation", type=float, default=0.05,
                     help="minimum M spacing enforced in the FINAL combined "
                          "winding table, independent of --dm (default 0.05). "
                          "This is deliberately NOT tied to --dm: scanning "
                          "with a fine --dm to localize a peak precisely, "
                          "then merging survivors at only that same fine "
                          "spacing, lets near-duplicate/collinear columns "
                          "into the final unregularized refit -- verified "
                          "directly to blow amplitudes up into the millions "
                          "while R stays deceptively ~0.999. Keep this at a "
                          "physically-sensible minimum winding separation "
                          "regardless of how fine --dm is set.")
    ap.add_argument("--inner-cv", type=int, default=5,
                     help="KFold splits used INSIDE each fold's training rows "
                          "to pick the LASSO alpha (default 5; auto-reduced "
                          "for small folds)")
    ap.add_argument("--max-harmonic", type=int, default=12,
                     help="harmonics jointly fit per candidate fundamental in "
                          "the harmonic-family scan (default 12; set to 0 to "
                          "disable the family scan and use only the per-M "
                          "LASSO dictionary)")
    ap.add_argument("--harmonic-keep-frac", type=float, default=0.15,
                     help="within the winning harmonic family, keep only "
                          "harmonics with amplitude >= this fraction of the "
                          "family's largest (default 0.15)")
    ap.add_argument("--target", choices=["data", "model"], default="data",
                     help="which column to fit: 'data' (Column 3, the real "
                          "default) or 'model' (Column 2, the existing "
                          "model's own deterministic output -- a "
                          "calibration check: since Column 2 is generated "
                          "from known M's, the search SHOULD cleanly "
                          "recover them here; if it can't, that is a bug "
                          "in the search, not a weak-signal issue)")
    ap.add_argument("--plot", action="store_true", help="save a per-fold correlation PNG")
    ap.add_argument("--outdir", type=Path, default=None)
    args = ap.parse_args()

    if args.window_years is None:
        args.window_years = args.initial_years

    time, model, data, manifold = load_columns(args.csv_path)
    if args.target == "model":
        fit_target, fit_reference = model, data
        target_label, reference_label = "Column 2 (model)", "Column 3 (data)"
    else:
        fit_target, fit_reference = data, model
        target_label, reference_label = "Column 3 (data)", "Column 2 (model)"
    print(f"[grid_lasso_optimizer] {len(data)} rows from {args.csv_path}, "
          f"record {time.min():.1f}-{time.max():.1f}, fitting {target_label} "
          f"({reference_label} kept as reference only)")

    rolling_folds, blocked_folds = [], []
    if args.mode in ("rolling", "both"):
        rolling_folds = rolling_cv(time, fit_reference, fit_target, manifold, args)
        print_folds("rolling (walk-forward)", rolling_folds)
    if args.mode in ("blocked", "both"):
        blocked_folds = blocked_cv(time, fit_reference, fit_target, manifold, args)
        print_folds("blocked (contiguous k-fold)", blocked_folds)

    print("\n[grid_lasso_optimizer] a model that only 'works' on one lucky "
          "split shows up here as a wide spread (large sd, min far below "
          "mean) across folds, not as one deceptively good number.")

    report_final_fit(args.csv_path, rolling_folds + blocked_folds, args,
                     target_label, reference_label)

    if args.plot and (rolling_folds or blocked_folds):
        out_dir = args.outdir if args.outdir is not None else Path(args.csv_path).resolve().parent
        plot_folds(args.csv_path, rolling_folds, blocked_folds,
                  out_dir / "grid_lasso_optimizer.png")


if __name__ == "__main__":
    main()
