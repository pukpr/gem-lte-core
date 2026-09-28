#!/usr/bin/env python3
"""withheld_interval_report.py -- post-run characterization of a fitted
model's withheld/test interval, meant to be read as a short blurb after a
run rather than eyeballed off a time-series plot.

Motivation (verbatim from the session this was built in): a single whole-
interval Pearson CC over a long (e.g. 15-year) withheld window conflates
two very different questions -- "how long does the fit stay good before
decorrelating" and "do the excursions have the right SHAPE even where the
exact values drift" -- the same conflation operational weather/climate
forecasts avoid by reporting skill-vs-lead-time curves instead of one
number. It also can't by itself catch a known failure mode: an overfit
model whose excursions have blown up in amplitude can still show a
deceptively good Pearson R by sign-coincidence alone.

This script reports four complementary numbers for one withheld
interval, each targeting one of those:

0. PEAK SHORT-WINDOW BURST -- caught live on pna's 2001-2016 gap: the
   50-month window used for the dead-reckoning horizon (below) has a blind spot
   covering the first ~4 years of ANY interval (it needs that many
   points before its first value even exists), and even where it can
   report, a 50-month average smears out a short (1-2 year) burst of
   genuinely strong local agreement. A separate, much shorter rolling
   window (default 18 months) is scanned across the WHOLE interval
   (including the 50-month window's own blind spot) to catch exactly
   that. Its significance is judged against the surrogate ensemble's own
   MAXIMUM correlation anywhere in the interval (not a single point) --
   the fair, multiple-comparisons-aware null for "is there a real burst
   somewhere," not just "is this one point high."

1. DEAD-RECKONING HORIZON (years) -- NOT "coherence" (that implies a
   stable phase-lock the model would need to actively maintain). What's
   really being measured is closer to ship's dead reckoning: right at a
   training-adjacent edge the model still has a good fix, and extrapolating
   from it stays accurate for a while on the strength of that anchor alone,
   before drift accumulates and it decorrelates -- "there is always SOME
   correlation for a little bit" near a real anchor point, not a claim
   that the two signals are staying locked in phase throughout (renamed
   from an earlier "coherence time" after exactly this correction).
   How far INTO the withheld interval, scanned from a training-adjacent
   edge, the rolling correlation stays above an AR(1)-surrogate
   significance floor. Two numbers per edge: a threshold-crossing horizon
   (where it first drops out, sustained) and a fitted exponential
   e-folding tau (smoother, less sensitive to one noisy dip). A one-sided
   interval (a forward tail, e.g. VALIDATE's own [Last, D'Last] window)
   has one edge (the start). A two-sided interval (an EXCLUDE=TRUE gap,
   bounded by training data on BOTH sides) has two, and each is scanned
   INDEPENDENTLY -- forward from the start, backward from the end --
   never merged onto one "distance from nearest edge" axis. An earlier
   version of this script did merge them, which lets a bad patch near one
   edge get reported as the horizon even though the other edge is still
   fine; confirmed live and fixed.

2. DTW "PHONEME" MATCH -- Model vs Data's own dtw_distance (this
   project's normalized, CC-like DTW score: 1.0 = perfect alignment, from
   lte_forward.py, a faithful port of the Ada METRIC=DTW primitive),
   reported as a z-score against an AR(1)-surrogate null -- exactly the
   same discipline winding_rank.py and this project's resolved-PDE
   checks already use, so "the excursions look right" becomes a real,
   falsifiable statistic instead of an eyeball call. DTW tolerates the
   timing wobble that kills Pearson CC first, so it's expected to survive
   longer into the withheld interval than the CC-based dead-reckoning horizon.

3. OVERFIT / BLOWN-UP-EXCURSION CHECK -- Model's own standard deviation
   over the withheld interval, as a ratio to Data's. Pearson CC is scale-
   invariant and can look good purely by sign-coincidence even when the
   model's excursions are wildly over- or under-scaled; DTW's own
   distance is NOT scale-invariant, so a real blowup should show up as a
   poor DTW z-score even when CC looks fine. The report flags this
   combination explicitly when it occurs.

Usage:
    withheld_interval_report.py <index_dir> --start YYYY --end YYYY
        [--window 50] [--dtw-window 3] [--n-surr 500] [--out FILE.png]

Example (this session's own case):
    withheld_interval_report.py sam --start 2000 --end 2015
    withheld_interval_report.py sam --start 2015 --end 2026.5   # VALIDATE's own tail
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from lte_forward import dtw_distance

EDGE_EPS_YEARS = 1.0  # within this of the record's own last date -> one-sided


def load_series(index_dir: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    raw = np.loadtxt(index_dir / "lte_results.csv", delimiter=",")
    return raw[:, 0], raw[:, 1], raw[:, 2]  # dates, model, data


def ar1_surrogate(n: int, phi: float, sigma: float, rng: np.random.Generator) -> np.ndarray:
    s = np.empty(n)
    s[0] = rng.normal(0.0, sigma)
    innov = sigma * math.sqrt(max(1e-12, 1.0 - phi * phi))
    for i in range(1, n):
        s[i] = phi * s[i - 1] + rng.normal(0.0, innov)
    return s


def rolling_corr(a: np.ndarray, b: np.ndarray, window: int) -> np.ndarray:
    n = len(a) - window + 1
    out = np.empty(n)
    for i in range(n):
        out[i] = np.corrcoef(a[i:i + window], b[i:i + window])[0, 1]
    return out


def find_local_peaks(t: np.ndarray, corr: np.ndarray, min_value: float,
                      min_separation_years: float) -> list[tuple[float, float]]:
    """Local maxima of `corr` above `min_value`, at least
    `min_separation_years` apart (keeping the higher of any two closer
    than that -- avoids listing the same burst 4-5 times just because
    it's a few points wide). Simple/dependency-free: no scipy."""
    is_peak = np.zeros(len(corr), dtype=bool)
    for i in range(1, len(corr) - 1):
        if corr[i] >= corr[i - 1] and corr[i] >= corr[i + 1] and corr[i] >= min_value:
            is_peak[i] = True
    candidates = sorted(zip(t[is_peak], corr[is_peak]), key=lambda p: -p[1])
    kept: list[tuple[float, float]] = []
    for tp, cp in candidates:
        if all(abs(tp - kt) >= min_separation_years for kt, _ in kept):
            kept.append((tp, cp))
    return sorted(kept)


def fit_exp_decay(distance: np.ndarray, corr: np.ndarray) -> float | None:
    """Fit corr ~ A*exp(-distance/tau); returns tau (years) or None if the
    fit is degenerate (e.g. corr never really decays, or goes negative
    immediately). Uses abs(corr) so a sign flip doesn't break the log-fit;
    only the ENVELOPE/magnitude decay is being characterized here."""
    y = np.abs(corr)
    y = np.clip(y, 1e-3, None)  # keep the log finite
    try:
        coef = np.polyfit(distance, np.log(y), 1)
    except Exception:
        return None
    slope = coef[0]
    if slope >= -1e-6:  # not decaying (or growing) -- no meaningful tau
        return None
    return -1.0 / slope


def dead_reckoning_horizon(distance: np.ndarray, corr: np.ndarray, floor: float,
                            sustain: int = 3) -> float | None:
    """First distance (scanning outward, ascending) at which `sustain`
    consecutive points fall below `floor` -- avoids reporting a single
    noisy dip as the horizon. None if it never happens (stays above floor
    throughout). Returns 0.0, not the window's own minimum measurable
    distance, if it is ALREADY below floor at the closest point available
    -- otherwise the horizon would be misreported as e.g. "4.08 years"
    (the rolling window's own warm-up requirement) when the honest
    finding is "never tracks, even at the edge," which reads as the
    opposite of what it means. BUG FIX: this distinction was documented
    in this docstring but not actually implemented -- caught live on
    pna, where the very first measurable point already failed and was
    reported as "coherent for 4.08 years."""
    order = np.argsort(distance)
    d, c = distance[order], corr[order]
    below = c < floor
    if below[:sustain].all():
        return 0.0
    for i in range(len(below) - sustain + 1):
        if below[i:i + sustain].all():
            return float(d[i])
    return None


def run(index_dir: Path, start: float, end: float, window: int,
        dtw_window: int, n_surr: int, out: Path | None, seed: int,
        short_window: int = 18) -> None:
    dates, model, data = load_series(index_dir)
    name = index_dir.name
    record_end = dates.max()
    one_sided = (record_end - end) < EDGE_EPS_YEARS

    mask = (dates >= start) & (dates < end)
    n = int(mask.sum())
    if n < window:
        raise SystemExit(f"withheld interval [{start},{end}) has only {n} points, "
                          f"< rolling window ({window}) -- widen the interval or "
                          f"shrink --window")
    m_win, d_win = model[mask], data[mask]
    t_win = dates[mask]

    print(f"=== {name}: withheld interval [{start:.2f}, {end:.2f})  "
          f"n={n} months  ({'one-sided tail' if one_sided else 'two-sided gap'}) ===")

    whole_cc = float(np.corrcoef(m_win, d_win)[0, 1])
    print(f"\nwhole-interval Pearson CC: {whole_cc:+.4f}  (the single number that "
          f"conflates the two questions below)")

    # ------------------------------------------------------------------
    # Shared AR(1) surrogate ensemble, fit to the withheld Data itself.
    # ------------------------------------------------------------------
    resid = d_win - d_win.mean()
    phi = float(np.clip(np.corrcoef(resid[:-1], resid[1:])[0, 1], -0.98, 0.98))
    sigma = float(d_win.std())
    rng = np.random.default_rng(seed)

    surrogates = [ar1_surrogate(n, phi, sigma, rng) + d_win.mean() for _ in range(n_surr)]

    # ------------------------------------------------------------------
    # 0. Peak short-window burst -- catches short (1-2 year) local
    # tracking the 50-month window (below) both structurally cannot see
    # for its own first ~4 years, and would smear out even where it can.
    # ------------------------------------------------------------------
    peak_roll = None
    peaks: list[tuple[float, float]] = []
    if n >= short_window:
        peak_roll = rolling_corr(m_win, d_win, short_window)
        peak_t = t_win[short_window - 1:]
        surr_peak_roll = np.array([rolling_corr(m_win, s, short_window) for s in surrogates])
        # Fair null for "is there a real burst ANYWHERE": each surrogate's
        # own single best moment, not a per-point comparison -- otherwise
        # a long, heavily-overlapping window scan will manufacture
        # "significant" peaks out of red noise alone by pure multiple
        # comparisons.
        surr_burst_null = surr_peak_roll.max(axis=1)
        real_burst = float(peak_roll.max())
        burst_z = (real_burst - surr_burst_null.mean()) / surr_burst_null.std()
        burst_pct = float(np.mean(surr_burst_null < real_burst)) * 100.0

        peaks = find_local_peaks(peak_t, peak_roll, min_value=0.5,
                                  min_separation_years=short_window / 12.0)

        print(f"\n--- 0. Peak short-window burst ({short_window}-month rolling CC) ---")
        print(f"  strongest burst anywhere in interval: r={real_burst:+.3f} at "
              f"t={peak_t[int(np.argmax(peak_roll))]:.2f}")
        print(f"  {n_surr} AR(1) surrogates' OWN best-anywhere burst: "
              f"mean={surr_burst_null.mean():.3f} std={surr_burst_null.std():.3f}")
        print(f"  z-score: {burst_z:+.2f}   percentile: {burst_pct:.1f}%  "
              f"(fair null: best burst ANYWHERE, not a single point)")
        if peaks:
            print(f"  distinct local peaks (r>=0.5, >={short_window/12.0:.1f}y apart):")
            for tp, cp in peaks:
                print(f"    t={tp:6.2f}   r={cp:+.3f}")
        else:
            print(f"  no local peaks reached r>=0.5")

    # ------------------------------------------------------------------
    # 1. Dead-reckoning horizon: rolling CC, scanned INWARD from each
    # training-adjacent edge independently. BUG FIX: an earlier version
    # used distance = min(t-start, end-t) and sorted ALL points (both
    # edges) together before scanning for the first sustained drop -- for
    # a two-sided gap this silently interleaves two UNRELATED decay
    # processes onto one axis, so a bad patch near one edge can be
    # reported as "the" horizon even though the other edge is fine
    # (confirmed live: sam's 2000-2015 gap reported 0.08y because the
    # 2015 edge is poor, even though the 2000 edge stays good for ~4
    # years on its own). Each edge is now scanned independently, in its
    # own natural direction (forward from start / backward from end),
    # and both are reported separately for a two-sided interval.
    # ------------------------------------------------------------------
    real_roll = rolling_corr(m_win, d_win, window)
    roll_t = t_win[window - 1:]

    surr_roll_all = np.array([rolling_corr(m_win, s, window) for s in surrogates])
    # Pooled significance floor: real correlation must clear this to still
    # count as "tracking" -- 2 sigma above the surrogate ensemble's own
    # mean, pooled across all lags/surrogates (matches this project's
    # standing AR(1)-floor convention elsewhere, e.g. winding_rank.py).
    floor = float(surr_roll_all.mean() + 2.0 * surr_roll_all.std())

    print(f"\n--- 1. Dead-reckoning horizon (each training-adjacent edge, scanned independently) ---")
    print(f"  AR(1)-surrogate significance floor (rolling CC, {window}-month window): {floor:+.4f}")

    def report_edge(horizon_edge: float | None, tau_edge: float | None, label: str) -> None:
        if horizon_edge is None:
            print(f"  [{label}] threshold-crossing: never drops below floor "
                  f"(tracks through all {end - start:.1f} years shown)")
        elif horizon_edge == 0.0:
            print(f"  [{label}] threshold-crossing: never reaches significance -- "
                  f"already below floor at the closest measurable point")
        else:
            print(f"  [{label}] dead-reckoning horizon: {horizon_edge:.2f} years")
        print(f"  [{label}] fitted e-folding tau: "
              f"{'no clear decay found' if tau_edge is None else f'{tau_edge:.2f} years'}")

    left_dist = roll_t - start
    horizon_left = dead_reckoning_horizon(left_dist, real_roll, floor)
    tau_left = fit_exp_decay(left_dist, real_roll)
    label_left = "since start" if one_sided else "from the START edge"
    report_edge(horizon_left, tau_left, label_left)

    horizon_right = tau_right = None
    if not one_sided:
        right_dist = end - roll_t
        horizon_right = dead_reckoning_horizon(right_dist, real_roll, floor)
        tau_right = fit_exp_decay(right_dist, real_roll)
        report_edge(horizon_right, tau_right, "from the END edge")

    # ------------------------------------------------------------------
    # 2. DTW "phoneme" pattern match, whole-interval + windowed.
    # ------------------------------------------------------------------
    real_dtw = dtw_distance(m_win, d_win, dtw_window)
    surr_dtw = np.array([dtw_distance(m_win, s, dtw_window) for s in surrogates])
    dtw_z = (real_dtw - surr_dtw.mean()) / surr_dtw.std()
    dtw_pct = float(np.mean(surr_dtw < real_dtw)) * 100.0

    print(f"\n--- 2. DTW pattern match (shape, not exact value) ---")
    print(f"  DTW score (this project's normalized scale, 1.0=perfect): {real_dtw:+.4f}")
    print(f"  {n_surr} AR(1) surrogates: mean={surr_dtw.mean():+.4f} std={surr_dtw.std():.4f}")
    print(f"  z-score: {dtw_z:+.2f}   percentile: {dtw_pct:.1f}%")

    windowed_dtw_z = None
    if n >= window * 2:
        n_roll = n - window + 1
        real_dtw_roll = np.array([dtw_distance(m_win[i:i + window], d_win[i:i + window], dtw_window)
                                   for i in range(n_roll)])
        surr_dtw_roll = np.array([[dtw_distance(m_win[i:i + window], s[i:i + window], dtw_window)
                                    for i in range(n_roll)] for s in surrogates])
        windowed_dtw_z = (real_dtw_roll - surr_dtw_roll.mean(axis=0)) / surr_dtw_roll.std(axis=0)

    # ------------------------------------------------------------------
    # 3. Overfit / blown-up-excursion check.
    # ------------------------------------------------------------------
    var_ratio = float(m_win.std() / d_win.std())
    print(f"\n--- 3. Excursion-scale check (catches spurious-R overfitting) ---")
    print(f"  std(Model)/std(Data) over withheld interval: {var_ratio:.3f}  "
          f"(1.0 = matched; far from 1.0 = excursions blown up or deflated)")
    blown_up = var_ratio > 2.0 or var_ratio < 0.5
    spurious_r = whole_cc > 0.4 and (dtw_z < 0.0 or blown_up)
    if blown_up:
        print(f"  FLAG: excursion amplitude is {'>2x' if var_ratio > 2.0 else '<0.5x'} "
              f"Data's own -- a real scale mismatch regardless of R.")
    if spurious_r:
        print(f"  WARNING: whole-interval CC ({whole_cc:+.4f}) looks good, but DTW's "
              f"scale-sensitive z-score ({dtw_z:+.2f}) and/or the excursion-scale ratio "
              f"above suggest this may be a SPURIOUS/overfit correlation -- large, "
              f"badly-scaled excursions coincidentally sign-matching real data -- "
              f"rather than genuine pattern skill. Trust the DTW z-score and the "
              f"dead-reckoning horizon over the raw R here.")
    else:
        print(f"  no spurious-R flag: DTW and CC broadly agree on fit quality.")

    # ------------------------------------------------------------------
    # Plot: one panel per edge (rolling CC + rolling DTW-z, sorted by
    # that edge's own distance -- each edge plotted independently, never
    # merged), plus the raw overlay.
    # ------------------------------------------------------------------
    if out is not None:
        def plot_edge(ax, dist: np.ndarray, horizon_edge: float | None, title: str) -> None:
            order = np.argsort(dist)
            ax.plot(dist[order], real_roll[order], color="tab:blue", lw=1.2,
                    label="rolling Pearson CC")
            ax.axhline(floor, color="0.4", ls="--", lw=0.8, label="AR(1) floor")
            if horizon_edge is not None:
                ax.axvline(horizon_edge, color="tab:red", ls=":", lw=1.2,
                           label=f"dead-reckoning horizon ({horizon_edge:.2f}y)")
            if windowed_dtw_z is not None:
                axb = ax.twinx()
                axb.plot(dist[order], windowed_dtw_z[order], color="tab:green",
                         lw=1.0, alpha=0.8, label="rolling DTW z-score")
                axb.set_ylabel("DTW z-score", color="tab:green")
                axb.tick_params(axis="y", labelcolor="tab:green")
                axb.axhline(0.0, color="tab:green", lw=0.5, alpha=0.3)
            ax.set_xlabel("years")
            ax.set_ylabel("rolling Pearson CC", color="tab:blue")
            ax.tick_params(axis="y", labelcolor="tab:blue")
            ax.set_title(title)
            ax.legend(loc="upper right", fontsize=7)

        n_rows = 2 if one_sided else 3
        fig, axes = plt.subplots(n_rows, 1, figsize=(10, 3.3 * n_rows))
        if one_sided:
            ax_left, ax_overlay = axes
        else:
            ax_left, ax_right, ax_overlay = axes
            plot_edge(ax_right, end - roll_t, horizon_right,
                      f"{name}: dead-reckoning horizon (CC) / pattern-match (DTW-z), "
                      f"scanned backward from the END edge ({end:g})")

        plot_edge(ax_left, roll_t - start, horizon_left,
                  f"{name}: dead-reckoning horizon (CC) / pattern-match (DTW-z), "
                  f"scanned forward from the START edge ({start:g})"
                  if not one_sided else
                  f"{name}: dead-reckoning horizon (CC) and pattern-match (DTW-z) "
                  f"since the start of the withheld interval")

        ax_overlay.plot(t_win, d_win, color="tab:blue", lw=0.8, alpha=0.7, label="Data")
        ax_overlay.plot(t_win, m_win, color="tab:red", lw=0.8, label="Model")
        if peak_roll is not None and peaks:
            for tp, cp in peaks:
                ax_overlay.axvspan(tp - short_window / 24.0, tp + short_window / 24.0,
                                    color="gold", alpha=0.25)
            ax_overlay.plot([], [], color="gold", alpha=0.5, lw=8,
                             label=f"local burst (r>=0.5, {short_window}mo)")
        ax_overlay.set_title(f"withheld interval overlay  (CC={whole_cc:+.3f}  DTW-z={dtw_z:+.2f}  "
                              f"std-ratio={var_ratio:.2f})")
        ax_overlay.set_xlabel("Year")
        ax_overlay.legend(fontsize=8)

        fig.tight_layout()
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=140)
        print(f"\nsaved {out}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("index_dir", type=Path, help="directory containing lte_results.csv")
    ap.add_argument("--start", type=float, required=True, help="withheld interval start (year)")
    ap.add_argument("--end", type=float, required=True, help="withheld interval end (year)")
    ap.add_argument("--window", type=int, default=50, help="rolling-correlation window, months (default 50)")
    ap.add_argument("--dtw-window", type=int, default=3, help="DTW Sakoe-Chiba band width (default 3)")
    ap.add_argument("--n-surr", type=int, default=500, help="AR(1) surrogates (default 500)")
    ap.add_argument("--out", type=Path, default=None, help="save a diagnostic plot here")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--short-window", type=int, default=18,
                     help="rolling-correlation window, months, for the peak local-burst "
                          "check (default 18) -- independent of --window, which is "
                          "sized for the sustained dead-reckoning horizon, not short bursts")
    args = ap.parse_args()
    run(args.index_dir, args.start, args.end, args.window, args.dtw_window,
        args.n_surr, args.out, args.seed, args.short_window)


if __name__ == "__main__":
    main()
