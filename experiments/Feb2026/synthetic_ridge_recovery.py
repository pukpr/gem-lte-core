#!/usr/bin/env python3
"""synthetic_ridge_recovery.py — a deterministic guide before trusting the
toolchain on real (ambiguous, red-noise-confounded) climate data: embed a
KNOWN winding number into synthetic data built on the REAL Baltic forcing
manifold, then run the exact same ridge-transfer machinery
(cv_ridge_transfer.py) against it and check whether it recovers the
answer we already know is there.

Why this matters right now: the spectral-slope check (spectral_slope_check.py)
came back genuinely ambiguous, confounded with ordinary AR(1) red noise --
a real signal and a stochastic look-alike can produce similar-looking
spectra. This script sidesteps that entirely by using a signal with a
KNOWN, EXACTLY-EMBEDDED mechanism, so "does the method recover the right
answer" has a ground-truth to check against instead of an open question.
It's also the right prerequisite for anything downstream (an Arnold-tongue
sweep of a nonlinear parametrically-forced network, or SINDy structure
discovery): both are only as trustworthy as the diagnostic used to read
their output, so that diagnostic needs its own recovery test first.

Construction: reuse the REAL Forcing(t) manifold from cv_rolling_blocked.
prepare() (so the frequency content is exactly as complicated/bursty as
the real thing, not a toy single sinusoid), and build

    data_synth(t) = sum_i amp_i * sin(2*pi * k_true_i * Forcing(t) + phase_i)
                    + noise(t)

matching the exact sin(2*pi*M*Forcing+phase) convention already used by
gem-lte-primitives.adb's LTE(), mode_evidence.fit_at_M, and
cv_ridge_transfer.ridge_scan -- so recovery is apples-to-apples with how
real indices are actually scored, not a different convention that happens
to look similar.

noise(t) is AR(1)-persistent (not white) at a chosen correlation rho and
target SNR, because the real record's own persistence (rho ~ 0.97-0.98,
measured in cv_ridge_transfer.py's own AR(1) floor) is what actually
threatens recovery -- a white-noise-only test would understate the real
difficulty, the same anti-conservative mistake WINDING_SCALOGRAM_FEASIBILITY.md
already flagged for white-noise normalization in general.

Two things get characterized as a function of SNR:
  1. single-shot recovery: does a whole-record ridge scan find k_true
     (within one grid cell, matching this project's own established
     synthetic-recovery standard)?
  2. transfer recovery: does cv_ridge_transfer's 3-interval split still
     converge all three independently-optimized thirds on k_true with
     near-zero regret, the same test that gave the striking real-Baltic
     result -- and at what SNR does that agreement start to break down?

Usage
-----
    ./synthetic_ridge_recovery.py
    ./synthetic_ridge_recovery.py --k-true 1.25 0.45 --amp 1.0 0.4
    ./synthetic_ridge_recovery.py --snr-list inf 10 3 1 0.3 0.1
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
from cv_rolling_blocked import prepare, valid_range  # noqa: E402
from cv_ridge_transfer import make_intervals, ridge_scan, ar1_floor, find_ridge  # noqa: E402
from mode_evidence import fit_at_M  # noqa: E402

ROOT = Path(__file__).resolve().parent


def make_synthetic(forcing: np.ndarray, k_true: list[float], amp: list[float],
                   phase: list[float]) -> np.ndarray:
    """sum_i amp_i * sin(2*pi*k_true_i*Forcing + phase_i) -- the exact
    convention used everywhere else in this toolchain."""
    y = np.zeros_like(forcing)
    for k, a, ph in zip(k_true, amp, phase):
        y += a * np.sin(2 * np.pi * k * forcing + ph)
    return y


def ar1_noise(n: int, rho: float, target_std: float, rng: np.random.Generator
             ) -> np.ndarray:
    if target_std == 0.0:
        return np.zeros(n)
    innovation_std = target_std * np.sqrt(max(1.0 - rho ** 2, 1e-9))
    eps = rng.normal(0.0, innovation_std, n)
    x = np.empty(n)
    x[0] = rng.normal(0.0, target_std)
    for i in range(1, n):
        x[i] = rho * x[i - 1] + eps[i]
    return x


def whole_record_recovery(forcing: np.ndarray, data: np.ndarray,
                          k_true: list[float], m_max: float, dm: float
                          ) -> dict:
    m_grid = np.arange(dm, m_max + dm / 2, dm)
    curve = ridge_scan(forcing, data, m_grid)
    i_star = int(np.argmax(curve))
    m_hat = float(m_grid[i_star])
    nearest_true = min(k_true, key=lambda k: abs(k - m_hat))
    return dict(m_hat=m_hat, cc=float(curve[i_star]), nearest_true=nearest_true,
                error=abs(m_hat - nearest_true), grid_cells_off=abs(m_hat - nearest_true) / dm)


def transfer_recovery(dates: np.ndarray, forcing: np.ndarray, data: np.ndarray,
                      k_true: list[float], m_max: float, dm: float,
                      n_intervals: int, n_reps: int, seed: int) -> dict:
    prep_synth = dict(dates=dates, forcing=forcing, data=data)
    intervals = make_intervals(prep_synth, n_intervals)
    m_grid = np.arange(dm, m_max + dm / 2, dm)
    ridges = []
    for iv in intervals:
        curve = ridge_scan(iv["forcing"], iv["data"], m_grid)
        floor, _rho = ar1_floor(iv["forcing"], iv["data"], m_grid, n_reps, seed + iv["i"])
        ridges.append(find_ridge(m_grid, curve, floor))
    m_hats = [r["m"] for r in ridges]
    pairwise = [abs(m_hats[i] - m_hats[j])
               for i in range(n_intervals) for j in range(i + 1, n_intervals)]
    regrets = []
    for i in range(n_intervals):
        for j in range(n_intervals):
            if i == j:
                continue
            transfer_cc = fit_at_M(intervals[j]["forcing"], intervals[j]["data"],
                                   ridges[i]["m"])
            regrets.append(ridges[j]["cc"] - transfer_cc)
    nearest_true = [min(k_true, key=lambda k: abs(k - m)) for m in m_hats]
    errors = [abs(m - t) for m, t in zip(m_hats, nearest_true)]
    return dict(m_hats=m_hats, mean_pairwise_disagreement=float(np.mean(pairwise)) if pairwise else 0.0,
                mean_regret=float(np.mean(regrets)) if regrets else float("nan"),
                mean_error=float(np.mean(errors)), all_significant=all(r["significant"] for r in ridges))


def parse_snr(s: str) -> float:
    return float("inf") if s.lower() == "inf" else float(s)


def run(idx: str, k_true: list[float], amp: list[float], phase: list[float],
        snr_list: list[float], ar1_rho: float, m_max: float, dm: float,
        n_intervals: int, n_reps: int, seed: int, outdir: Path | None) -> None:
    prep = prepare(idx)
    dates, forcing = prep["dates"], prep["forcing"]
    signal = make_synthetic(forcing, k_true, amp, phase)
    signal_std = float(np.std(signal))
    rng = np.random.default_rng(seed)

    print(f"-- synthetic recovery on {idx}'s real Forcing(t), embedded "
          f"k_true={k_true}, amp={amp} --")
    print(f"  signal std: {signal_std:.4f}   noise: AR(1) rho={ar1_rho:.2f}")

    rows = []
    for snr in snr_list:
        noise_std = 0.0 if np.isinf(snr) else signal_std / np.sqrt(snr)
        noise = ar1_noise(len(dates), ar1_rho, noise_std, rng)
        data_synth = signal + noise

        whole = whole_record_recovery(forcing, data_synth, k_true, m_max, dm)
        transfer = transfer_recovery(dates, forcing, data_synth, k_true, m_max,
                                     dm, n_intervals, n_reps, seed)
        rows.append(dict(snr=snr, whole=whole, transfer=transfer))

        snr_label = "inf (noiseless)" if np.isinf(snr) else f"{snr:.2g}"
        print(f"\n  SNR={snr_label}:")
        print(f"    whole-record: M_hat={whole['m_hat']:.3f}  "
              f"(true {whole['nearest_true']:.3f}, off by "
              f"{whole['grid_cells_off']:.1f} grid cells)  cc={whole['cc']:+.3f}")
        print(f"    {n_intervals}-way transfer: M_hats={[f'{m:.3f}' for m in transfer['m_hats']]}  "
              f"mean pairwise disagreement={transfer['mean_pairwise_disagreement']:.3f}  "
              f"mean error vs true={transfer['mean_error']:.3f}  "
              f"mean regret={transfer['mean_regret']:+.3f}  "
              f"all AR1-significant={transfer['all_significant']}")

    plot_recovery(idx, rows, k_true, dm,
                 (outdir if outdir is not None else ROOT / idx) /
                 "synthetic_ridge_recovery.png")


def plot_recovery(idx: str, rows: list[dict], k_true: list[float], dm: float,
                  out_path: Path) -> None:
    finite_snrs = [r["snr"] for r in rows if not np.isinf(r["snr"])]
    plot_snr = [1e3 if np.isinf(r["snr"]) else r["snr"] for r in rows]

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    ax = axes[0]
    ax.plot(plot_snr, [r["whole"]["grid_cells_off"] for r in rows], "o-",
           label="whole-record")
    ax.plot(plot_snr, [r["transfer"]["mean_error"] / dm for r in rows], "s-",
           label=f"transfer mean (per interval)")
    ax.set_xscale("log")
    ax.axhline(1.0, color="0.6", linestyle=":", linewidth=1,
              label="1 grid cell (this project's own recovery standard)")
    ax.set_xlabel("SNR (signal std / noise std, log; rightmost = noiseless)")
    ax.set_ylabel("recovery error (grid cells)")
    ax.set_title(f"{idx}: recovery error vs noise", fontsize=10)
    ax.legend(fontsize=7)

    ax2 = axes[1]
    ax2.plot(plot_snr, [r["transfer"]["mean_pairwise_disagreement"] for r in rows],
             "o-", color="tab:purple", label="mean |M*_i - M*_j|")
    ax2.plot(plot_snr, [r["transfer"]["mean_regret"] for r in rows], "s-",
            color="tab:green", label="mean transfer regret")
    ax2.axhline(0.0, color="0.6", linewidth=0.8)
    ax2.set_xscale("log")
    ax2.set_xlabel("SNR (log; rightmost = noiseless)")
    ax2.set_ylabel("disagreement / regret")
    ax2.set_title(f"{idx}: transfer stability vs noise", fontsize=10)
    ax2.legend(fontsize=7)

    fig.suptitle(f"Synthetic ground truth: k_true={k_true} embedded in "
                f"{idx}'s real Forcing(t)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"\n  saved {out_path}")


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--idx", default="baltic",
                     help="whose real Forcing(t) to build synthetic data on "
                          "(default baltic)")
    ap.add_argument("--k-true", type=float, nargs="+", default=[1.250],
                     help="winding number(s) to embed (default: baltic's own "
                          "cv_ridge_transfer.py finding, 1.250)")
    ap.add_argument("--amp", type=float, nargs="+", default=None,
                     help="amplitude per k-true (default: 1.0 for each)")
    ap.add_argument("--phase", type=float, nargs="+", default=None,
                     help="phase per k-true, radians (default: 0.0 for each)")
    ap.add_argument("--snr-list", type=str, nargs="+",
                     default=["inf", "10", "3", "1", "0.3", "0.1"],
                     help="signal-std/noise-std ratios to sweep, 'inf' = "
                          "noiseless (default: inf 10 3 1 0.3 0.1)")
    ap.add_argument("--ar1-rho", type=float, default=0.97,
                     help="noise persistence (default 0.97, matching this "
                          "project's own measured Baltic/AR(1) floor)")
    ap.add_argument("--m-max", type=float, default=5.0)
    ap.add_argument("--dm", type=float, default=0.01)
    ap.add_argument("--n-intervals", type=int, default=3)
    ap.add_argument("--n-reps", type=int, default=40)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--outdir", type=Path, default=None)
    args = ap.parse_args()

    n = len(args.k_true)
    amp = args.amp if args.amp is not None else [1.0] * n
    phase = args.phase if args.phase is not None else [0.0] * n
    if len(amp) != n or len(phase) != n:
        raise SystemExit("--amp and --phase must match --k-true in length")
    snr_list = [parse_snr(s) for s in args.snr_list]

    run(args.idx, args.k_true, amp, phase, snr_list, args.ar1_rho, args.m_max,
       args.dm, args.n_intervals, args.n_reps, args.seed, args.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
