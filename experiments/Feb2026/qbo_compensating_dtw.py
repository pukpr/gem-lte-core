#!/usr/bin/env python3
"""Does the second-largest-amplitude Doodson term (9.1207d, annually
aliased to ~22yr) explain the shape of the slow multi-decadal wander in
qbo30's winding-phase residual (qbo_phase_lock.py)?

Reconstruction: from qbo30/lt.exe.p's "lpap" array ([period_days,
amplitude, phase] per Doodson term, matching Tide_Sum's own
Amplitude*Cos(2*pi*(Year_Len/Period)*Time + Phase) convention), evaluate
each fast candidate cosine directly at the record's actual monthly sample
times -- no manual aliasing needed, undersampling a ~9-day wave at a
~30-day cadence produces the slow apparent beat on its own (the
"stroboscopic sampling" mechanism already described in this directory's
own QBO_TIDAL_THEORY_VERIFICATION.md). Also built directly: the explicitly
demodulated slow cosine (frequency taken mod 1 cycle/year), which is what
Tide_Sum computes internally when ALIAS=true (see Aliased_Period,
gem-lte-primitives.adb:41) -- same idea, no sampling required.

Comparison tool: classic Dynamic Time Warping via an explicit O(N*M)
dynamic-programming table (no library) -- lets the candidate's slow wander
be compared against the phase residual's shape even if the two aren't
aligned in absolute time or locally stretched, which is exactly the
"fix wander" property DTW is for. A random-phase-shifted null version of
each candidate gives the baseline DTW cost an unrelated same-shaped signal
would achieve, since DTW's warping freedom will find *some* alignment
between almost any two curves -- the raw distance alone isn't meaningful
without that comparison.
"""
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from wavelet_scalogram import load_columns
from qbo_phase_lock import locked_M, phase_residual

ROOT = Path(__file__).resolve().parent
YEAR_LEN = 365.242484


def standardize(x):
    x = np.asarray(x, dtype=float)
    return (x - x.mean()) / x.std()


def dtw_distance(x, y):
    """Classic DTW via explicit dynamic programming. Returns
    (normalized_distance, path) where path is a list of (i, j) index pairs
    from (0,0) to (n-1, m-1)."""
    n, m = len(x), len(y)
    INF = np.inf
    D = np.full((n + 1, m + 1), INF)
    D[0, 0] = 0.0
    for i in range(1, n + 1):
        xi = x[i - 1]
        for j in range(1, m + 1):
            cost = abs(xi - y[j - 1])
            D[i, j] = cost + min(D[i - 1, j], D[i, j - 1], D[i - 1, j - 1])

    # backtrack
    i, j = n, m
    path = []
    while (i, j) != (0, 0):
        path.append((i - 1, j - 1))
        choices = [(D[i - 1, j], (i - 1, j)), (D[i, j - 1], (i, j - 1)),
                   (D[i - 1, j - 1], (i - 1, j - 1))]
        _, (i, j) = min(choices, key=lambda c: c[0])
    path.append((0, 0))
    path.reverse()
    return D[n, m] / (n + m), path


def lpap_term(idx: str, target_period_days: float, tol=0.01):
    p = json.loads((ROOT / idx / "lt.exe.p").read_text())
    best = min(p["lpap"], key=lambda e: abs(e[0] - target_period_days))
    if abs(best[0] - target_period_days) > tol:
        raise ValueError(f"no lpap term within {tol}d of "
                          f"{target_period_days} days (closest: {best[0]})")
    return tuple(best)


def raw_signal(year, terms):
    s = np.zeros_like(year)
    for period_d, amp, phase in terms:
        freq_cyc_per_yr = YEAR_LEN / period_d
        s += amp * np.cos(2.0 * np.pi * freq_cyc_per_yr * year + phase)
    return s


def aliased_signal(year, terms):
    s = np.zeros_like(year)
    for period_d, amp, phase in terms:
        freq = YEAR_LEN / period_d
        freq_aliased = freq - round(freq)  # matches Long_Float'Remainder(freq,1.0)
        s += amp * np.cos(2.0 * np.pi * freq_aliased * year + phase)
    return s


def null_baseline(rng, residual_s, candidate, n_trials=8):
    """DTW cost of the REAL residual against random-phase-shifted
    (circularly rolled) copies of the candidate -- same amplitude/frequency
    content as the real candidate, unrelated timing. This is the honest
    comparison point: DTW's warping freedom can partially align almost any
    two same-shaped curves, so the raw distance for the true (unshifted)
    candidate only means something relative to this."""
    costs = []
    for _ in range(n_trials):
        shift = rng.integers(len(candidate) // 4, 3 * len(candidate) // 4)
        y = standardize(np.roll(candidate, shift))
        d, _ = dtw_distance(residual_s, y)
        costs.append(d)
    return np.mean(costs), np.std(costs)


def main():
    idx = "qbo30"
    year, dt, model, obs, forcing = load_columns(idx)
    M, decay, mP = locked_M(idx)
    residual, slope = phase_residual(year, forcing, M)
    residual_s = standardize(residual)

    t9a = lpap_term(idx, 9.1207)
    t9b = lpap_term(idx, 9.1085)
    t27 = lpap_term(idx, 27.2122)
    print(f"[{idx}] lpap terms used: 9.1207d={t9a}  9.1085d={t9b}  27.2122d={t27}")

    candidates = {
        "raw 9.1207d alone (monthly-sampled)": raw_signal(year, [t9a]),
        "raw 9.1207d+9.1085d (monthly-sampled)": raw_signal(year, [t9a, t9b]),
        "aliased 9.1207d+9.1085d (demodulated, ~22yr+~10yr)":
            aliased_signal(year, [t9a, t9b]),
        "aliased 27.2122d (demodulated, ~2.37yr) -- expected NOT to match":
            aliased_signal(year, [t27]),
    }

    rng = np.random.default_rng(0)
    print()
    print(f"{'candidate':55s} {'DTW dist':>10s} {'null mean':>10s} "
          f"{'null std':>9s} {'z-score':>8s}")
    results = {}
    for name, sig in candidates.items():
        sig_s = standardize(sig)
        d, path = dtw_distance(residual_s, sig_s)
        null_mean, null_std = null_baseline(rng, residual_s, sig)
        z = (null_mean - d) / null_std if null_std > 0 else float("nan")
        results[name] = (d, path, sig_s, null_mean, null_std, z)
        print(f"{name:55s} {d:10.4f} {null_mean:10.4f} {null_std:9.4f} "
              f"{z:8.2f}")

    # best (lowest-distance, i.e. most improved over its own null) candidate
    best_name = min(results, key=lambda k: results[k][0])
    d, path, sig_s, null_mean, null_std, z = results[best_name]

    fig, axes = plt.subplots(3, 1, figsize=(11, 11))

    axes[0].plot(year, residual_s, color="tab:blue", label="phase residual (standardized)")
    axes[0].plot(year, sig_s, color="tab:red", alpha=0.8,
                 label=f"best candidate: {best_name}")
    axes[0].set_title(f"{idx}: phase residual vs. best-matching candidate "
                       f"(direct overlay, no warping)", fontweight="bold")
    axes[0].legend(loc="upper right", fontsize=8)
    axes[0].set_ylabel("standardized")

    path_arr = np.array(path)
    axes[1].plot(year[path_arr[:, 0]], year[path_arr[:, 1]], color="black",
                 linewidth=0.7)
    axes[1].plot([year.min(), year.max()], [year.min(), year.max()],
                 color="gray", linestyle="--", linewidth=0.7,
                 label="no-warp diagonal")
    axes[1].set_xlabel("residual time (yr)")
    axes[1].set_ylabel("candidate time (yr)")
    axes[1].set_title("DTW warping path -- deviation from the diagonal is "
                       "how much timing wander is being absorbed",
                       fontweight="bold")
    axes[1].legend(loc="upper left", fontsize=8)

    warped_sig = sig_s[path_arr[:, 1]]
    warped_year = year[path_arr[:, 0]]
    axes[2].plot(warped_year, residual_s[path_arr[:, 0]], color="tab:blue",
                 label="phase residual")
    axes[2].plot(warped_year, warped_sig, color="tab:red", alpha=0.8,
                 label=f"{best_name} (DTW-warped onto residual)")
    axes[2].set_title("after DTW alignment", fontweight="bold")
    axes[2].set_xlabel("year")
    axes[2].legend(loc="upper right", fontsize=8)

    fig.tight_layout()
    out = ROOT / idx / "qbo_compensating_dtw.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"\n  saved {out}")


if __name__ == "__main__":
    main()
