#!/usr/bin/env python3
"""cv_rolling_blocked.py — rolling-origin (walk-forward) and blocked
time-series cross-validation for the GEM-LTE forward model.

Every other CV-flavored knob in this project (TRAIN_START/TRAIN_END/EXCLUDE
in lt.exe.resp, --cv in lte_forward.py) does ONE train/test split, chosen by
hand, and reports one correlation. That's fine for headline numbers but
says nothing about whether the fit is stable across the record or just
lucky on the interval someone happened to pick. This script automates many
splits, walking through the record in two different disciplined ways:

  rolling (walk-forward): strictly causal -- train only on the past, score
    on the immediately-following block, then advance. This is what "does
    this model forecast" actually means for a time series: test data is
    always chronologically after train data, exactly once each.

  blocked (contiguous k-fold): split the record into K contiguous chunks;
    for each chunk, train on ALL the others (past AND future) and score on
    the held-out chunk. This is NOT a forecast test -- future data is
    allowed to inform the fit -- but it's a fairer test of "does the fitted
    manifold generalize to an unseen interval" than a single hand-picked
    split, and it's a standard alternative to plain (shuffled) k-fold for
    autocorrelated series, since blocks stay temporally contiguous.

Both reuse lte_forward.py's ported pipeline (Doodson periods -> jerk ->
tide sum -> impulse comb -> IIR -> Bessel -> regression -> LTE response),
so a fold's "fit" is byte-for-byte the same computation the Ada optimizer
performs for one candidate parameter set -- only the train/test row
selection changes. The forcing series itself never depends on the fold: it
is a fixed function of the astronomical parameters and the date grid, not
of which points get regressed. Only step 7 onward (the OLS fit and the
correlation) is fold-dependent.

Caveat worth keeping in mind when reading the numbers: `data` is smoothed
with a 9-point (1,2,1)/4 boxcar (F9) BEFORE the fold split, which mixes
each point with roughly its four neighbors on each side (~0.33yr for
monthly data). That's a small amount of information leaking across a
train/test boundary that sits inside the smoothed run. --embargo (blocked
mode only; default 0.5yr) drops a buffer of points adjacent to each test
block from its training set to cover this. Rolling folds don't get an
embargo since a walk-forward split is already one-directional -- the only
leakage there would be from a handful of points just past the boundary,
which the F9 window already reaches barely at all.

Usage
-----
    ./cv_rolling_blocked.py baltic nino4
    ./cv_rolling_blocked.py nino4 --mode rolling --initial-years 30 \\
        --step-years 10 --horizon-years 10
    ./cv_rolling_blocked.py baltic --mode blocked --n-blocks 8 --embargo-years 1.0
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lte_forward import (                                    # noqa: E402
    Config, read_resp, doodson_periods, year_length, JERK_REFERENCE_PERIOD,
    tide_sum, impulse_delta, iir, bessel, regression_factors, lte_response,
    filter9point,
)

ROOT = Path(__file__).resolve().parent
DEFAULT_INDICES = ["baltic", "nino4"]


# ---------------------------------------------------------------------------
# Shared prep: build the forcing manifold once per index, independent of
# any train/test split (steps 1-6 of lte_forward.forward()).
# ---------------------------------------------------------------------------

def prepare(idx: str) -> dict:
    """Reproduces lte_forward.forward()'s steps 1-6 with no CV-window
    applied yet, so the same (dates, data, forcing, m) can be reused across
    every fold. `data` is the F9-filtered series, matching what forward()
    regresses against; `data_raw` is the .dat file's own series before F9.
    `tide_raw` (raw Doodson tide sum, pre impulse-comb/IIR/Bessel) and
    `forcing_pre_bessel` (post-IIR, pre-Bessel) are also returned so a
    caller can tell how much of Forcing's own spectral shape comes from
    the pipeline's literal linear integrator (IIR) versus anything
    downstream."""
    idx_dir = ROOT / idx
    params = json.loads((idx_dir / "lt.exe.p").read_text())
    resp = read_resp(idx_dir / "lt.exe.resp")
    resp["CLIMATE_INDEX"] = str(idx_dir / f"{idx}.dat")
    overrides = {}
    if "harm" in params:
        # lt.exe.resp's own NH line is stale/hand-edited (often left at
        # "1 1 1 1", which repeats the SAME base frequency instead of
        # harmonics of it and makes the regression design singular) --
        # lt.exe.p's "harm" field is what the optimizer actually converged
        # on. Same override param_survey.compute_winding applies.
        overrides["NH"] = " ".join(str(int(float(h))) for h in params["harm"])
    cfg = Config(resp, overrides)

    dat = Path(cfg.get("CLIMATE_INDEX", f"{idx}.dat"))
    raw = np.loadtxt(dat)
    dates, data = raw[:, 0].copy(), raw[:, 1].copy()
    data_raw = data.copy()  # pre-F9, as stored in the .dat file

    lpap = np.array(params["lpap"], dtype=float)
    lt = np.array(params["ltep"], dtype=float)
    year_cand = float(params.get("year", 0.0))
    B = {k: float(params.get(k, 0.0)) for k in
         ("offs", "bg", "impA", "impB", "impC", "delA", "delB", "asym",
          "ann1", "ann2", "IR", "ma", "mp", "shfT", "init")}

    year_startup = cfg.get("YEAR", 0.0)
    jerk = cfg.get("JERK", 0.0)
    sampling = int(cfg.get("SAMPLING", 12.0))
    f9 = cfg.get("F9", 1)
    idate = cfg.get("IDATE", 0.0)
    nonlin = cfg.get("NONLIN", 1.0)
    nm = cfg.get("NM", len(lt))
    nh_str = cfg.get("NH", "")
    harms = [int(float(h)) for h in str(nh_str).split()]
    trend_on = cfg.get("TREND", True)

    if f9 > 0:
        for _ in range(f9):
            data = filter9point(data)

    periods = doodson_periods(year_startup, year_cand)
    ap = lpap[:, 1:3].copy()
    if jerk != 0.0:
        w_norm, w_der = 1.0 - jerk, jerk * JERK_REFERENCE_PERIOD / periods
        ap[:, 0] *= np.sqrt(w_norm ** 2 + w_der ** 2)
        ap[:, 1] += np.arctan2(w_der, w_norm)

    tf = tide_sum(dates, ap, periods, year_length(year_startup, year_cand),
                  scaling=0.0, integ=B["shfT"])
    forcing = iir(tf * impulse_delta(dates, B["delA"], B["delB"], B["asym"],
                                     sampling),
                  lag_a=1.0 - B["ma"], lag_c=B["mp"], init=B["init"],
                  start_date=idate, dates=dates)
    forcing_pre_bessel = forcing.copy()  # after IIR, before the Bessel step
    m = np.concatenate([lt[:nm], [lt[nm - 1] * h for h in harms]])
    forcing = bessel(forcing, B["impA"], B["impB"], m[nm - 1],
                     B["offs"], B["bg"])

    return dict(idx=idx, dates=dates, data=data, forcing=forcing, m=m,
                nonlin=nonlin, trend_on=trend_on, f9=f9, ir=B["IR"],
                data_raw=data_raw, tide_raw=tf,
                forcing_pre_bessel=forcing_pre_bessel)


def valid_range(prep: dict) -> tuple[float, float]:
    """(t_min, t_max) over non-placeholder (data != 0) points -- mirrors
    lte_forward.cc()'s own leading/trailing zero trim, so fold boundaries
    aren't quietly padded out into a run of sentinel zeros."""
    nz = np.nonzero(prep["data"] != 0.0)[0]
    return float(prep["dates"][nz[0]]), float(prep["dates"][nz[-1]])


# ---------------------------------------------------------------------------
# Fit on one index subset, score on another
# ---------------------------------------------------------------------------

def fit_score(prep: dict, train_idx: np.ndarray, test_idx: np.ndarray
              ) -> dict | None:
    """Regress on train_idx only, build the model over the whole record,
    then score strictly on test_idx -- a genuine held-out-interval
    correlation, unlike forward()'s own cc() (which always scores the same
    rows it fit). Returns None if either side has too few usable points."""
    dates, data, forcing, m = (prep["dates"], prep["data"], prep["forcing"],
                                prep["m"])
    train_idx = train_idx[data[train_idx] != 0.0]
    test_idx = test_idx[data[test_idx] != 0.0]
    if len(train_idx) < 2 * (2 * len(m) + 7) or len(test_idx) < 3:
        return None  # not enough rows to identify the regression / to score

    level, k0, amp, phase, annual, trend, accel = regression_factors(
        dates[train_idx], forcing[train_idx], data[train_idx], m,
        prep["trend_on"])
    model = lte_response(forcing, dates, m, amp, phase, level, k0, trend,
                         accel, prep["nonlin"], annual)
    if prep["ir"] != 0.0:
        model = model.copy()
        for i in range(len(model) - 1, 11, -1):
            model[i] -= prep["ir"] * model[i - 12]
    for _ in range(prep["f9"]):
        model = filter9point(model)

    x, y = model[test_idx], data[test_idx]
    if np.std(x) == 0.0 or np.std(y) == 0.0:
        return None
    r = float(np.corrcoef(x, y)[0, 1])
    return dict(cc=r, n_train=len(train_idx), n_test=len(test_idx),
                train_span=(float(dates[train_idx].min()),
                            float(dates[train_idx].max())),
                test_span=(float(dates[test_idx].min()),
                           float(dates[test_idx].max())))


# ---------------------------------------------------------------------------
# Rolling-origin (walk-forward) CV
# ---------------------------------------------------------------------------

def rolling_cv(prep: dict, initial_years: float, step_years: float,
               horizon_years: float, window: str = "expanding",
               window_years: float | None = None) -> list[dict]:
    dates = prep["dates"]
    t_min, t_max = valid_range(prep)
    folds = []
    test_start = t_min + initial_years
    while test_start < t_max:
        test_end = min(test_start + horizon_years, t_max)
        if test_end - test_start < horizon_years * 0.5:
            break  # trailing sliver too short to score meaningfully
        if window == "sliding" and window_years is not None:
            train_start = max(t_min, test_start - window_years)
        else:
            train_start = t_min
        train_idx = np.nonzero((dates >= train_start) & (dates < test_start))[0]
        test_idx = np.nonzero((dates >= test_start) & (dates < test_end))[0]
        result = fit_score(prep, train_idx, test_idx)
        if result is not None:
            folds.append(result)
        test_start += step_years
    return folds


# ---------------------------------------------------------------------------
# Blocked (contiguous k-fold) CV
# ---------------------------------------------------------------------------

def blocked_cv(prep: dict, n_blocks: int, embargo_years: float = 0.5
               ) -> list[dict]:
    dates = prep["dates"]
    t_min, t_max = valid_range(prep)
    edges = np.linspace(t_min, t_max, n_blocks + 1)
    folds = []
    for k in range(n_blocks):
        b_start, b_end = edges[k], edges[k + 1]
        test_idx = np.nonzero((dates >= b_start) & (dates < b_end))[0]
        train_idx = np.nonzero((dates < b_start - embargo_years) |
                               (dates >= b_end + embargo_years))[0]
        result = fit_score(prep, train_idx, test_idx)
        if result is not None:
            result["block"] = k
            folds.append(result)
    return folds


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def summarize(ccs: list[float]) -> str:
    if not ccs:
        return "no scoreable folds"
    arr = np.array(ccs)
    return (f"mean={arr.mean():+.3f}  median={np.median(arr):+.3f}  "
            f"sd={arr.std():.3f}  min={arr.min():+.3f}  max={arr.max():+.3f}"
            f"  n={len(arr)}")


def print_folds(label: str, folds: list[dict]) -> None:
    print(f"  {label}:")
    if not folds:
        print("    (no folds -- record too short for these settings)")
        return
    for i, f in enumerate(folds):
        tr = f"{f['train_span'][0]:.1f}-{f['train_span'][1]:.1f}"
        te = f"{f['test_span'][0]:.1f}-{f['test_span'][1]:.1f}"
        print(f"    fold {i:2d}  train={tr:>13}  (n={f['n_train']:4d})   "
              f"test={te:>13}  (n={f['n_test']:3d})   cc={f['cc']:+.3f}")
    print(f"    -> {summarize([f['cc'] for f in folds])}")


def plot_folds(idx: str, rolling: list[dict], blocked: list[dict],
               out_path: Path) -> None:
    panels = [(rolling, "rolling (walk-forward)", "tab:blue"),
              (blocked, "blocked (contiguous k-fold)", "tab:orange")]
    panels = [(f, t, c) for f, t, c in panels if f]
    if not panels:
        return
    fig, axes = plt.subplots(1, len(panels), figsize=(6 * len(panels), 4.2),
                              squeeze=False)
    for ax, (folds, title, color) in zip(axes[0], panels):
        ccs = [f["cc"] for f in folds]
        centers = [0.5 * (f["test_span"][0] + f["test_span"][1])
                   for f in folds]
        ax.axhline(0.0, color="0.6", linewidth=0.8)
        ax.bar(centers, ccs, width=(centers[1] - centers[0]) * 0.8
               if len(centers) > 1 else 5.0, color=color, alpha=0.85)
        ax.axhline(np.mean(ccs), color="black", linestyle="--", linewidth=1,
                   label=f"mean={np.mean(ccs):+.3f}")
        ax.set_title(f"{idx}: {title}", fontweight="bold", fontsize=10)
        ax.set_xlabel("test-window center (year)")
        ax.set_ylabel("held-out correlation")
        ax.set_ylim(-1.0, 1.0)
        ax.legend(fontsize=8)
    fig.suptitle(f"{idx}: held-out (never-fitted) correlation per fold",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"  saved {out_path}")


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("indices", nargs="*", default=DEFAULT_INDICES)
    ap.add_argument("--mode", choices=["rolling", "blocked", "both"],
                     default="both")
    ap.add_argument("--initial-years", type=float, default=30.0,
                     help="rolling: length of the first training window (yr)")
    ap.add_argument("--step-years", type=float, default=10.0,
                     help="rolling: how far the origin advances per fold (yr)")
    ap.add_argument("--horizon-years", type=float, default=10.0,
                     help="rolling: length of each held-out test block (yr)")
    ap.add_argument("--window", choices=["expanding", "sliding"],
                     default="expanding",
                     help="rolling: grow the training window (default) or "
                          "keep it a fixed --window-years length")
    ap.add_argument("--window-years", type=float, default=None,
                     help="rolling: training-window length when --window "
                          "sliding (default: --initial-years)")
    ap.add_argument("--n-blocks", type=int, default=6,
                     help="blocked: number of contiguous blocks")
    ap.add_argument("--embargo-years", type=float, default=0.5,
                     help="blocked: buffer dropped from training data on "
                          "each side of the held-out block")
    ap.add_argument("--plot", action="store_true",
                     help="save a per-index PNG of fold correlations")
    ap.add_argument("--outdir", type=Path, default=None)
    args = ap.parse_args()

    if args.window_years is None:
        args.window_years = args.initial_years

    for idx in args.indices:
        print(f"-- {idx} --")
        idx_dir = ROOT / idx
        if not (idx_dir / "lt.exe.p").exists():
            print(f"  (skipping {idx}: no lt.exe.p)")
            continue
        prep = prepare(idx)
        t_min, t_max = valid_range(prep)
        print(f"  record: {t_min:.2f}-{t_max:.2f}  "
              f"({sum(prep['data'] != 0.0)} nonzero of {len(prep['data'])} pts)")

        rolling_folds, blocked_folds = [], []
        if args.mode in ("rolling", "both"):
            rolling_folds = rolling_cv(
                prep, args.initial_years, args.step_years,
                args.horizon_years, window=args.window,
                window_years=args.window_years)
            print_folds("rolling (walk-forward)", rolling_folds)
        if args.mode in ("blocked", "both"):
            blocked_folds = blocked_cv(prep, args.n_blocks,
                                       args.embargo_years)
            print_folds("blocked (contiguous k-fold)", blocked_folds)

        if args.plot and (rolling_folds or blocked_folds):
            out_dir = args.outdir if args.outdir is not None else idx_dir
            plot_folds(idx, rolling_folds, blocked_folds,
                       out_dir / "cv_rolling_blocked.png")

    return 0


if __name__ == "__main__":
    sys.exit(main())
