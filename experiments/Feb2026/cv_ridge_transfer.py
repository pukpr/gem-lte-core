#!/usr/bin/env python3
"""cv_ridge_transfer.py — a structural (not predictive) cross-validation:
optimize the forcing manifold's resonant wavenumber on one interval, then
ask whether that same ridge location -- not the fitted values, the
LOCATION in M -- also shows up as the best-supported wavenumber on a
disjoint interval it never saw.

Rolling/blocked CV (cv_rolling_blocked.py) asks "does the fitted model
predict held-out VALUES". This asks a different, complementary question:
if the shared tidal manifold is real (WINDING_NARRATIVE.md's central
claim), its resonant wavenumber is a property of the external forcing, not
of whichever stretch of noisy climate data happened to be fit -- so
independently optimizing on two non-overlapping intervals should converge
on close to the same M, and using interval A's ridge on interval B's data
should cost little relative to letting B pick its own ridge. If the
manifold were instead an artifact of overfitting one interval's noise, the
ridge location would drift around with no reason to agree, and A's M
would transfer no better than a randomly-chosen one.

Method, per interval:
  1. Sweep candidate wavenumbers M over a grid and, at each M, fit
     offset + k0*Forcing + amp*sin(2*pi*M*Forcing) + amp*cos(2*pi*M*Forcing)
     by plain OLS against that interval's own data (mode_evidence.py's
     fit_at_M/sweep_M, reused verbatim -- this is the same single-mode
     "ridge scan" primitive that script uses to sanity-check one fitted
     mode, applied here over the whole grid on two separate intervals).
  2. The interval's ridge M* is the sweep's argmax.
  3. Significance: compare against matched-length AR(1) surrogates (rho
     estimated from that interval's own detrended data, matching
     winding_rank.py's red-noise-floor philosophy) run through the same
     sweep -- a ridge that doesn't clear the 95th-percentile AR(1) floor
     is flagged, not trusted.

Then, for every ordered pair of intervals (i predicting j):
  - transfer cc  = fit_at_M(forcing_j, data_j, M*_i)   -- freeze i's ridge,
    refit only amplitude/phase on j, score on j
  - own cc       = fit_at_M(forcing_j, data_j, M*_j)   -- j's own best
  - regret       = own cc - transfer cc                 (0 = perfect transfer)
  - chance level = median of j's own sweep curve         (what a blindly
    chosen M would typically score on j -- the baseline "no information"
    transfer would need to beat)

This is exploratory: a small number of intervals on a noisy climate index,
not a controlled experiment. Read the regret/chance numbers as evidence,
not proof, of manifold stability.

Usage
-----
    ./cv_ridge_transfer.py baltic nino4
    ./cv_ridge_transfer.py nino4 --n-intervals 3 --m-max 3 --dm 0.005
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
from cv_rolling_blocked import prepare, valid_range, DEFAULT_INDICES  # noqa: E402
from mode_evidence import fit_at_M                                    # noqa: E402

ROOT = Path(__file__).resolve().parent


def make_intervals(prep: dict, n_intervals: int) -> list[dict]:
    """n_intervals equal-length contiguous, non-overlapping chunks spanning
    the valid (non-placeholder) date range."""
    dates = prep["dates"]
    t_min, t_max = valid_range(prep)
    edges = np.linspace(t_min, t_max, n_intervals + 1)
    intervals = []
    for i in range(n_intervals):
        lo, hi = edges[i], edges[i + 1]
        idx = np.nonzero((dates >= lo) & (dates < hi + (1e-9 if i == n_intervals - 1 else 0))
                         & (prep["data"] != 0.0))[0]
        intervals.append(dict(i=i, lo=lo, hi=hi, idx=idx,
                              forcing=prep["forcing"][idx], data=prep["data"][idx]))
    return intervals


def ridge_scan(forcing: np.ndarray, data: np.ndarray, m_grid: np.ndarray
               ) -> np.ndarray:
    return np.array([fit_at_M(forcing, data, m) for m in m_grid])


def ar1_floor(forcing: np.ndarray, data: np.ndarray, m_grid: np.ndarray,
              n_reps: int, seed: int) -> tuple[np.ndarray, float]:
    """95th-percentile sweep curve from AR(1) surrogates matched to
    `data`'s own (detrended) lag-1 autocorrelation and variance -- the
    red-noise floor a real ridge needs to clear to be trusted, following
    winding_rank.py's rationale but estimating rho per-interval instead of
    assuming a fixed 0.97."""
    n = len(data)
    t = np.arange(n, dtype=float)
    trend = np.polyfit(t, data, 1)
    resid = data - np.polyval(trend, t)
    if resid.std() == 0.0:
        return np.zeros_like(m_grid), 0.0
    rho = float(np.clip(np.corrcoef(resid[:-1], resid[1:])[0, 1], 0.0, 0.99))
    sigma_eps = resid.std() * math_sqrt(1.0 - rho ** 2)
    rng = np.random.default_rng(seed)
    curves = np.empty((n_reps, len(m_grid)))
    for r in range(n_reps):
        eps = rng.normal(0.0, sigma_eps, n)
        surrogate = np.empty(n)
        surrogate[0] = resid[0]
        for i in range(1, n):
            surrogate[i] = rho * surrogate[i - 1] + eps[i]
        curves[r] = ridge_scan(forcing, surrogate, m_grid)
    return np.percentile(curves, 95, axis=0), rho


def math_sqrt(x: float) -> float:
    return x ** 0.5 if x > 0.0 else 0.0


def find_ridge(m_grid: np.ndarray, curve: np.ndarray, floor: np.ndarray
               ) -> dict:
    i_star = int(np.argmax(curve))
    return dict(m=float(m_grid[i_star]), cc=float(curve[i_star]),
                significant=bool(curve[i_star] > floor[i_star]),
                floor_at_ridge=float(floor[i_star]))


def run_index(idx: str, m_max: float, dm: float, n_intervals: int,
              n_reps: int, seed: int, outdir: Path | None) -> None:
    idx_dir = ROOT / idx
    if not (idx_dir / "lt.exe.p").exists():
        print(f"  (skipping {idx}: no lt.exe.p)")
        return
    prep = prepare(idx)
    intervals = make_intervals(prep, n_intervals)
    m_grid = np.arange(dm, m_max + dm / 2, dm)  # skip M=0 (degenerate: sin/cos(0)=const)

    print(f"-- {idx} --  ({n_intervals} interval(s), M in [0,{m_max}] step {dm})")
    curves, floors, ridges = [], [], []
    for iv in intervals:
        curve = ridge_scan(iv["forcing"], iv["data"], m_grid)
        floor, rho = ar1_floor(iv["forcing"], iv["data"], m_grid, n_reps, seed + iv["i"])
        ridge = find_ridge(m_grid, curve, floor)
        curves.append(curve); floors.append(floor); ridges.append(ridge)
        flag = "significant" if ridge["significant"] else "NOT sig. vs AR1 floor"
        print(f"  interval {iv['i']}  {iv['lo']:.1f}-{iv['hi']:.1f}  "
              f"(n={len(iv['idx']):4d}, AR1 rho={rho:.2f})  "
              f"ridge M*={ridge['m']:.3f}  cc={ridge['cc']:+.3f}  ({flag})")

    print("  transfer (train ridge on row, score+chance on column):")
    header = "         " + "".join(f"{'->'+str(j):>12}" for j in range(n_intervals))
    print(header)
    for i in range(n_intervals):
        cells = []
        for j in range(n_intervals):
            if i == j:
                cells.append(f"{'own':>12}")
                continue
            transfer_cc = fit_at_M(intervals[j]["forcing"], intervals[j]["data"],
                                   ridges[i]["m"])
            own_cc = ridges[j]["cc"]
            regret = own_cc - transfer_cc
            cells.append(f"{regret:+.3f}({transfer_cc:+.2f})".rjust(12))
        print(f"  from {i:2d}: " + "".join(cells))
    for j in range(n_intervals):
        chance = float(np.median(curves[j]))
        print(f"  interval {j} chance level (median over M-grid): {chance:+.3f}  "
              f"vs own ridge cc {ridges[j]['cc']:+.3f}")
    dm_agree = [abs(ridges[i]["m"] - ridges[j]["m"])
                for i in range(n_intervals) for j in range(i + 1, n_intervals)]
    if dm_agree:
        print(f"  |M*_i - M*_j| across all interval pairs: "
              f"mean={np.mean(dm_agree):.3f}  max={np.max(dm_agree):.3f}  "
              f"(grid step dm={dm:g})")

    out_dir = outdir if outdir is not None else idx_dir
    plot_ridges(idx, m_grid, curves, floors, ridges, intervals,
                out_dir / "cv_ridge_transfer.png")


def plot_ridges(idx: str, m_grid: np.ndarray, curves: list[np.ndarray],
                floors: list[np.ndarray], ridges: list[dict],
                intervals: list[dict], out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(11, 5.5))
    colors = plt.cm.tab10(np.linspace(0, 1, len(intervals)))
    for iv, curve, floor, ridge, color in zip(intervals, curves, floors,
                                              ridges, colors):
        label = f"{iv['lo']:.0f}-{iv['hi']:.0f}  (M*={ridge['m']:.3f})"
        ax.plot(m_grid, curve, color=color, linewidth=1.2, label=label)
        ax.plot(m_grid, floor, color=color, linewidth=0.7, linestyle=":",
                alpha=0.6)
        marker = "o" if ridge["significant"] else "x"
        ax.plot(ridge["m"], ridge["cc"], marker=marker, color=color,
                markersize=8, markeredgecolor="black", markeredgewidth=0.6)
    ax.axhline(0.0, color="0.7", linewidth=0.8)
    ax.set_xlabel("candidate winding number  M")
    ax.set_ylabel("single-mode OLS correlation")
    ax.set_title(
        f"{idx}: ridge-location scan per interval "
        f"(dotted = matched-AR(1) 95th-pct floor; o = clears floor, x = doesn't)",
        fontsize=10)
    ax.legend(fontsize=8, title="interval (fit range)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"  saved {out_path}")


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("indices", nargs="*", default=DEFAULT_INDICES)
    ap.add_argument("--m-max", type=float, default=5.0,
                     help="max winding number scanned (default 5.0, matching "
                          "winding_scalogram.py's own default)")
    ap.add_argument("--dm", type=float, default=0.01,
                     help="winding-number grid resolution (default 0.01)")
    ap.add_argument("--n-intervals", type=int, default=2,
                     help="number of equal contiguous intervals to compare "
                          "(default 2: first half vs second half)")
    ap.add_argument("--n-reps", type=int, default=40,
                     help="AR(1) surrogate repetitions for the noise floor")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--outdir", type=Path, default=None)
    args = ap.parse_args()

    for idx in args.indices:
        run_index(idx, args.m_max, args.dm, args.n_intervals, args.n_reps,
                  args.seed, args.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
