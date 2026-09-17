#!/usr/bin/env python3
"""sindy_discovery.py — sparse structure discovery (Brunton et al.'s SINDy,
core algorithm only: sequential thresholded least squares, hand-implemented
rather than adding the pysindy dependency, matching this project's own
"implement it from scratch" convention, e.g. wavelet_scalogram.py's own
Morlet CWT) applied to a question the rest of this session's tools all
presuppose the answer to: is a response of the form
sin(2*pi*M*Forcing(t)+phase) actually the right SPARSE family to explain
the data, or does the data need something else -- in particular, genuine
autonomous memory (the system referencing its own past), rather than pure
external (non-autonomous) forcing?

Framing (the user's point, worth stating precisely rather than assuming):
Doodson-style harmonic tidal analysis has always been non-autonomous --
sea level is represented as a sum of sinusoids at EXTERNALLY PRESCRIBED
astronomical frequencies, not derived from an autonomous free oscillator.
That's 100+ years of established, unglamorous practice. sin(k*M(t)) is a
nonlinear generalization of exactly that same non-autonomous idea (the
sinusoid's argument is still an externally-prescribed function of time,
just no longer a bare linear phase) -- it is NOT a claim that the ocean is
some exotic self-sustained nonlinear oscillator. The nonlinear MATH makes
it look more exotic than it is; the underlying claim is parsimonious and
lands squarely in the same non-autonomous tradition.

This script makes that distinction the actual variable under test, rather
than an assumption: build ONE combined candidate library with both kinds
of term side by side --

  non-autonomous (Doodson/LTE-style): sin(2*pi*M*Forcing), cos(2*pi*M*Forcing)
    over a wavenumber grid, plus Forcing, Forcing^2

  autonomous (genuine memory): the series' own lagged values,
    data(t - 1mo), data(t - 3mo), data(t - 6mo), data(t - 12mo)

-- and let sparse regression (no functional form assumed, no side
preferred) pick whichever actually survives thresholding. If the survivors
are almost entirely from the non-autonomous side, that's the sharpest
evidence available that this really is "Doodson's approach, nonlinearly
generalized" rather than something requiring exotic self-sustained
dynamics.

Validated on the SAME synthetic ground truth as synthetic_ridge_recovery.py
first (known k_true, real Forcing(t), realistic AR(1) noise) before ever
touching real data -- if the method can't cleanly prefer the correct
non-autonomous terms when we already know that's the truth, it has no
business being trusted on data where we don't know the answer.

Usage
-----
    ./sindy_discovery.py validate                 # synthetic ground truth
    ./sindy_discovery.py validate --k-true 1.25 0.45 --snr 3
    ./sindy_discovery.py real baltic nino4         # apply to real indices
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cv_rolling_blocked import prepare, DEFAULT_INDICES        # noqa: E402
from synthetic_ridge_recovery import make_synthetic, ar1_noise  # noqa: E402

ROOT = Path(__file__).resolve().parent
AR_LAGS_MONTHS = [1, 3, 6, 12]


# ---------------------------------------------------------------------------
# Library construction
# ---------------------------------------------------------------------------

def build_library(forcing: np.ndarray, data: np.ndarray, m_grid: np.ndarray,
                  ar_lags: list[int]) -> tuple[np.ndarray, list[str], list[str],
                                               np.ndarray]:
    """Combined [non-autonomous | autonomous] candidate matrix, trimmed to
    align with the longest AR lag (lagged terms lose their first
    max(ar_lags) rows). Returns (Theta, names, kind_per_column, y) where
    kind is "auto" or "nonauto" per column, and y is data trimmed the same way."""
    max_lag = max(ar_lags)
    n = len(data)

    cols, names, kinds = [], [], []
    # non-autonomous: plain Forcing terms
    cols.append(np.ones(n)); names.append("const"); kinds.append("nonauto")
    cols.append(forcing.copy()); names.append("Forcing"); kinds.append("nonauto")
    cols.append(forcing ** 2); names.append("Forcing^2"); kinds.append("nonauto")
    # non-autonomous: Doodson/LTE-style winding library
    for m in m_grid:
        cols.append(np.sin(2 * np.pi * m * forcing))
        names.append(f"sin(2pi*{m:.3f}*F)")
        kinds.append("nonauto")
        cols.append(np.cos(2 * np.pi * m * forcing))
        names.append(f"cos(2pi*{m:.3f}*F)")
        kinds.append("nonauto")
    # autonomous: the series' own lagged past (genuine memory candidates)
    for lag in ar_lags:
        cols.append(np.concatenate([np.full(lag, np.nan), data[:-lag]]))
        names.append(f"data(t-{lag}mo)")
        kinds.append("auto")

    Theta = np.column_stack(cols)[max_lag:]
    y = data[max_lag:].copy()
    return Theta, names, kinds, y


def winding_m(name: str) -> float | None:
    """Extract M from a "sin(2pi*{M}*F)" or "cos(2pi*{M}*F)" term name;
    None for any other (non-winding) library term."""
    if "(2pi*" not in name:
        return None
    return float(name.split("*")[1])


# ---------------------------------------------------------------------------
# STLSQ: SINDy's core algorithm, standardized columns, ridge-stabilized
# ---------------------------------------------------------------------------

def stlsq(Theta: np.ndarray, y: np.ndarray, threshold: float,
         n_iters: int = 15, ridge: float = 1e-6) -> np.ndarray:
    """Sequential thresholded least squares on STANDARDIZED columns (so one
    threshold is comparable across terms of very different raw scale --
    sin/cos are O(1), Forcing and its lags are not). Returns coefficients
    in the STANDARDIZED basis; caller maps back to named terms via the
    same scaling."""
    n_col = Theta.shape[1]
    reg = ridge * np.eye(n_col)
    coef = np.linalg.solve(Theta.T @ Theta + reg, Theta.T @ y)
    active = np.ones(n_col, dtype=bool)
    for _ in range(n_iters):
        small = active & (np.abs(coef) < threshold)
        if not np.any(small):
            break
        active[small] = False
        coef[~active] = 0.0
        if not np.any(active):
            break
        sub = Theta[:, active]
        coef[active] = np.linalg.solve(
            sub.T @ sub + ridge * np.eye(sub.shape[1]), sub.T @ y)
    return coef


def sparse_fit(Theta: np.ndarray, y: np.ndarray, names: list[str],
               kinds: list[str], threshold: float) -> dict:
    mu, sigma = Theta.mean(axis=0), Theta.std(axis=0)
    sigma[sigma == 0.0] = 1.0
    Theta_s = (Theta - mu) / sigma
    y_mu, y_sigma = y.mean(), y.std()
    y_s = (y - y_mu) / y_sigma

    coef_s = stlsq(Theta_s, y_s, threshold)
    fitted_s = Theta_s @ coef_s
    r2 = 1.0 - np.sum((y_s - fitted_s) ** 2) / np.sum((y_s - y_s.mean()) ** 2)

    survivors = [(names[i], kinds[i], float(coef_s[i]))
                for i in np.argsort(-np.abs(coef_s)) if coef_s[i] != 0.0]
    auto_var = float(np.sum(coef_s[np.array(kinds) == "auto"] ** 2))
    nonauto_var = float(np.sum(coef_s[np.array(kinds) == "nonauto"] ** 2))
    total = auto_var + nonauto_var if (auto_var + nonauto_var) > 0 else 1.0
    return dict(r2=float(r2), survivors=survivors,
                auto_share=auto_var / total, nonauto_share=nonauto_var / total,
                n_survivors=len(survivors))


# ---------------------------------------------------------------------------
# Properly-nulled autonomy test
#
# The raw "autonomous variance share" above is misleading on its own: at
# the AR(1) persistence this project's own data actually has (rho~0.97,
# per cv_ridge_transfer.py's own floor), a lag-1 term explains rho^2 = 0.94
# of variance from pure, mechanism-free noise persistence ALONE -- no real
# memory required. The SNR sweep below reproduces exactly this: as noise
# grows, "autonomous share" balloons even though the synthetic ground
# truth is 100% non-autonomous by construction. Comparing against ordinary
# noise, not against zero, is required -- the same AR(1)-floor discipline
# cv_ridge_transfer.py and winding_rank.py already apply to ridge
# significance, now applied to the autonomy question instead.
# ---------------------------------------------------------------------------

def fit_nonauto_only(forcing: np.ndarray, data: np.ndarray, m_grid: np.ndarray,
                     threshold: float) -> np.ndarray:
    """STLSQ restricted to just the non-autonomous (Forcing-only) columns;
    returns the RAW-scale residual data - fitted (no standardization
    needed for a residual used only for its own further analysis)."""
    cols = [np.ones_like(forcing), forcing, forcing ** 2]
    for m in m_grid:
        cols.append(np.sin(2 * np.pi * m * forcing))
        cols.append(np.cos(2 * np.pi * m * forcing))
    Theta = np.column_stack(cols)
    mu, sigma = Theta.mean(axis=0), Theta.std(axis=0)
    sigma[sigma == 0.0] = 1.0
    Theta_s = (Theta - mu) / sigma
    y_mu, y_sigma = data.mean(), data.std()
    y_s = (data - y_mu) / y_sigma
    coef_s = stlsq(Theta_s, y_s, threshold)
    fitted = (Theta_s @ coef_s) * y_sigma + y_mu
    return data - fitted


def lag_design(x: np.ndarray, ar_lags: list[int]) -> tuple[np.ndarray, np.ndarray]:
    max_lag = max(ar_lags)
    cols = [x[max_lag - lag: len(x) - lag] for lag in ar_lags]
    return np.column_stack(cols), x[max_lag:]


def ols_r2(X: np.ndarray, y: np.ndarray) -> float:
    coef, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    fitted = X @ coef
    ss_res = np.sum((y - fitted) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0


def autonomy_significance(residual: np.ndarray, ar_lags: list[int],
                          n_surrogates: int, seed: int) -> dict:
    """How much of the residual (what's left after the best non-autonomous
    fit) do its own lags explain -- and is that more than an AR(1) process
    of the SAME measured persistence would explain by chance alone?"""
    resid = residual - residual.mean()
    rho = float(np.clip(np.corrcoef(resid[:-1], resid[1:])[0, 1], 0.0, 0.999))
    X_real, y_real = lag_design(resid, ar_lags)
    r2_real = ols_r2(X_real, y_real)

    rng = np.random.default_rng(seed)
    n = len(resid)
    sigma_eps = resid.std() * np.sqrt(max(1.0 - rho ** 2, 1e-9))
    null_r2 = np.empty(n_surrogates)
    for i in range(n_surrogates):
        eps = rng.normal(0.0, sigma_eps, n)
        surrogate = np.empty(n)
        surrogate[0] = resid[0]
        for t in range(1, n):
            surrogate[t] = rho * surrogate[t - 1] + eps[t]
        X_s, y_s = lag_design(surrogate, ar_lags)
        null_r2[i] = ols_r2(X_s, y_s)

    percentile = float(np.mean(null_r2 < r2_real) * 100.0)
    return dict(rho=rho, r2_real=r2_real, null_mean=float(null_r2.mean()),
                null_p95=float(np.percentile(null_r2, 95)),
                percentile=percentile,
                significant=bool(r2_real > np.percentile(null_r2, 95)))


# ---------------------------------------------------------------------------
# Validation on synthetic ground truth
# ---------------------------------------------------------------------------

def validate(idx: str, k_true: list[float], amp: list[float], phase: list[float],
            snr: float, ar1_rho: float, m_max: float, dm: float,
            threshold: float, seed: int) -> None:
    prep = prepare(idx)
    dates, forcing = prep["dates"], prep["forcing"]
    signal = make_synthetic(forcing, k_true, amp, phase)
    rng = np.random.default_rng(seed)
    noise_std = 0.0 if np.isinf(snr) else float(np.std(signal)) / np.sqrt(snr)
    data = signal + ar1_noise(len(dates), ar1_rho, noise_std, rng)

    m_grid = np.arange(dm, m_max + dm / 2, dm)
    Theta, names, kinds, y = build_library(forcing, data, m_grid, AR_LAGS_MONTHS)
    result = sparse_fit(Theta, y, names, kinds, threshold)

    snr_label = "inf (noiseless)" if np.isinf(snr) else f"{snr:.2g}"
    print(f"-- SINDy validation on synthetic ground truth ({idx}'s Forcing, "
          f"k_true={k_true}, SNR={snr_label}) --")
    print(f"  library: {len(names)} candidates "
          f"({sum(1 for k in kinds if k == 'nonauto')} non-autonomous / "
          f"{sum(1 for k in kinds if k == 'auto')} autonomous), "
          f"threshold={threshold}")
    print(f"  R^2 of sparse fit: {result['r2']:.3f}   "
          f"survivors: {result['n_survivors']}")
    print(f"  raw explained-variance share (MISLEADING alone, see below): "
          f"non-autonomous={result['nonauto_share']:.1%}  "
          f"autonomous={result['auto_share']:.1%}")
    print("  surviving terms (standardized coefficient, largest first):")
    for name, kind, c in result["survivors"][:15]:
        print(f"    [{kind:7s}] {name:24s} {c:+.3f}")
    recovered = sorted({winding_m(nm) for nm, k, c in result["survivors"]
                        if k == "nonauto" and winding_m(nm) is not None})
    if recovered:
        nearest = [min(k_true, key=lambda t: abs(t - r)) for r in recovered]
        errs = [abs(r - t) for r, t in zip(recovered, nearest)]
        print(f"  recovered winding numbers (from surviving sin/cos terms): "
              f"{[f'{r:.3f}' for r in recovered]}")
        print(f"  best match to true k_true={k_true}: "
              f"max error {max(errs):.3f} ({max(errs)/dm:.1f} grid cells)")
    else:
        print("  ** no non-autonomous winding term survived thresholding -- "
              "recovery FAILED at this SNR/threshold **")

    resid = fit_nonauto_only(forcing, data, m_grid, threshold)
    sig = autonomy_significance(resid, AR_LAGS_MONTHS, n_surrogates=200, seed=seed)
    print(f"  AR(1)-nulled autonomy test on the non-autonomous fit's residual "
          f"(rho={sig['rho']:.2f}):")
    print(f"    real lag-R^2={sig['r2_real']:.3f}  vs  matched-AR(1)-noise "
          f"null: mean={sig['null_mean']:.3f}, p95={sig['null_p95']:.3f}  "
          f"(percentile {sig['percentile']:.0f})")
    if sig["significant"]:
        print("    -> residual has MORE lag-structure than ordinary noise "
              "persistence explains -- ground truth is 100% non-autonomous "
              "here, so this would be a FALSE POSITIVE at this threshold/SNR.")
    else:
        print("    -> residual's lag-structure is fully consistent with "
              "ordinary AR(1) noise persistence -- correctly finds no real "
              "autonomous memory, matching the known ground truth.")


# ---------------------------------------------------------------------------
# Application to real data
# ---------------------------------------------------------------------------

def run_real(idx: str, m_max: float, dm: float, threshold: float,
            outdir: Path | None) -> None:
    prep = prepare(idx)
    forcing, data = prep["forcing"], prep["data_raw"]
    m_grid = np.arange(dm, m_max + dm / 2, dm)
    Theta, names, kinds, y = build_library(forcing, data, m_grid, AR_LAGS_MONTHS)
    result = sparse_fit(Theta, y, names, kinds, threshold)

    print(f"-- SINDy on real {idx} data --")
    print(f"  library: {len(names)} candidates, threshold={threshold}")
    print(f"  R^2 of sparse fit: {result['r2']:.3f}   "
          f"survivors: {result['n_survivors']}")
    print(f"  raw explained-variance share (MISLEADING alone, see below): "
          f"non-autonomous={result['nonauto_share']:.1%}  "
          f"autonomous={result['auto_share']:.1%}")
    print("  surviving terms (standardized coefficient, largest first):")
    for name, kind, c in result["survivors"][:15]:
        print(f"    [{kind:7s}] {name:24s} {c:+.3f}")

    resid = fit_nonauto_only(forcing, data, m_grid, threshold)
    r2_nonauto_only = 1.0 - np.var(resid) / np.var(data)
    print(f"  non-autonomous-only fit (Forcing library alone, no lag terms "
          f"available at all): R^2={r2_nonauto_only:.3f}")
    sig = autonomy_significance(resid, AR_LAGS_MONTHS, n_surrogates=200, seed=0)
    print(f"  AR(1)-nulled autonomy test on the non-autonomous fit's residual "
          f"(rho={sig['rho']:.2f}):")
    print(f"    real lag-R^2={sig['r2_real']:.3f}  vs  matched-AR(1)-noise "
          f"null: mean={sig['null_mean']:.3f}, p95={sig['null_p95']:.3f}  "
          f"(percentile {sig['percentile']:.0f})")
    if sig["significant"]:
        print("    -> residual has MORE lag-structure than ordinary AR(1) "
              "noise persistence explains: there IS evidence of genuine "
              "autonomous memory beyond the forced (Doodson/LTE-style) "
              "response -- worth scrutinizing before trusting a purely "
              "forced-response story here.")
    else:
        print("    -> residual's lag-structure is fully consistent with "
              "ordinary AR(1) noise persistence at this rho: no evidence "
              "of genuine autonomous memory beyond the external forcing -- "
              "consistent with the parsimonious, Doodson-tradition, "
              "non-autonomous explanation.")

    plot_survivors(idx, result, m_grid,
                  (outdir if outdir is not None else ROOT / idx) /
                  "sindy_discovery.png")


def plot_survivors(idx: str, result: dict, m_grid: np.ndarray,
                   out_path: Path) -> None:
    winding = [(winding_m(nm), c) for nm, k, c in result["survivors"]
              if k == "nonauto" and winding_m(nm) is not None]
    fig, ax = plt.subplots(figsize=(9, 4))
    if winding:
        ms, cs = zip(*winding)
        ax.stem(ms, np.abs(cs), basefmt=" ")
    ax.set_xlim(0, m_grid[-1])
    ax.set_xlabel("surviving winding number M")
    ax.set_ylabel("|standardized coefficient|")
    ax.set_title(f"{idx}: SINDy-surviving non-autonomous winding terms "
                f"(auto share={result['auto_share']:.0%}, "
                f"R^2={result['r2']:.2f})", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"\n  saved {out_path}")


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="mode", required=True)

    v = sub.add_parser("validate", help="run on synthetic ground truth")
    v.add_argument("--idx", default="baltic")
    v.add_argument("--k-true", type=float, nargs="+", default=[1.250])
    v.add_argument("--amp", type=float, nargs="+", default=None)
    v.add_argument("--phase", type=float, nargs="+", default=None)
    v.add_argument("--snr", type=float, default=3.0,
                   help="signal-std/noise-std (use inf for noiseless)")
    v.add_argument("--ar1-rho", type=float, default=0.97)
    v.add_argument("--seed", type=int, default=0)

    r = sub.add_parser("real", help="apply to real index data")
    r.add_argument("indices", nargs="*", default=DEFAULT_INDICES)
    r.add_argument("--outdir", type=Path, default=None)

    for p in (v, r):
        p.add_argument("--m-max", type=float, default=5.0)
        p.add_argument("--dm", type=float, default=0.05,
                       help="wavenumber grid step for the library (default "
                            "0.05 -- coarser than cv_ridge_transfer.py's "
                            "0.01 default; a finer grid makes adjacent "
                            "sin/cos(M*F) columns nearly collinear, which "
                            "destabilizes which single one STLSQ keeps)")
        p.add_argument("--threshold", type=float, default=0.08)

    args = ap.parse_args()
    if args.mode == "validate":
        n = len(args.k_true)
        amp = args.amp if args.amp is not None else [1.0] * n
        phase = args.phase if args.phase is not None else [0.0] * n
        validate(args.idx, args.k_true, amp, phase, args.snr, args.ar1_rho,
                args.m_max, args.dm, args.threshold, args.seed)
    else:
        for idx in args.indices:
            run_real(idx, args.m_max, args.dm, args.threshold, args.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
