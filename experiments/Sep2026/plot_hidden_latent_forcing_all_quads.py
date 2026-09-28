#!/usr/bin/env python3
"""Variation of ../Feb2026/plot_hidden_latent_forcing.py, scaled from the
7 hardcoded flagship series to all 89 k_sst grid quads. Same underlying
signal (lte_results.csv column 4 -- the post-Bessel "hidden latent
forcing" / autonomous manifold, per this project's convention that it
should be invariant across INIT_DATE/IDATE) and the same core math
(zero-mean/unit-RMS normalization, Pearson correlation), but the
visualization changes shape at this scale:

  - 89 overlapping raw-colored lines is unreadable (categorical color
    only supports ~8 distinguishable hues) -- so individual quads are
    drawn thin/low-alpha in a single hue, with the cross-quad ensemble
    mean drawn bold on top (spaghetti + mean, the standard many-series
    pattern).
  - The correlation inset becomes the main event: an 89x89 diverging
    heatmap (matplotlib coolwarm, matching the original script's own
    choice), geographically ordered (latitude then longitude) so any
    spatial coherence structure is visible, with no per-cell numeric
    text (89x89 = 7921 cells -- illegible at any font size; magnitude by
    color is the right encoding here, not a label on every cell).

For each cell, reads lte_results.csv from whichever of Feb2026/<cell>/
or Sep2026/<cell>/ has the newer mtime (this session's fixes landed in
different places depending on how each cell got solved -- TEST_ONLY
verifications in Feb2026, donor-seeded automated attempts in Sep2026 --
so "freshest" is the correct proxy for "reflects the final accepted
lt.exe.p", not a fixed preference for one directory root).
"""
import csv
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, "/home/paul/eval/gem-lte-core/experiments/Sep2026")
import sweep

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

NAME_RE = re.compile(r"^k([NS])(\d{3})_([EW])(\d{3})$")


def decode(name: str):
    m = NAME_RE.match(name)
    ns, lat_s, ew, lon_s = m.groups()
    return float(lat_s) * (1 if ns == "N" else -1), float(lon_s) * (1 if ew == "E" else -1)


def freshest_results_csv(cell: str) -> Path | None:
    feb = sweep.FEB / cell / "lte_results.csv"
    sep = sweep.HERE / cell / "lte_results.csv"
    candidates = [p for p in (feb, sep) if p.exists()]
    if not candidates:
        return None
    return max(candidates, key=lambda p: p.stat().st_mtime)


def read_time_and_forcing(path: Path) -> tuple[list[float], list[float]]:
    times, forcings = [], []
    with path.open(newline="", encoding="utf-8") as f:
        for row in csv.reader(f, skipinitialspace=True):
            if not row or len(row) < 4:
                continue
            times.append(float(row[0]))
            forcings.append(float(row[3]))
    return times, forcings


def normalize_excursions(values: list[float]) -> list[float]:
    mean = sum(values) / len(values)
    excursions = [v - mean for v in values]
    rms = math.sqrt(sum(v * v for v in excursions) / len(excursions))
    if rms == 0.0:
        return excursions
    return [v / rms for v in excursions]


def pearson_correlation(left: list[float], right: list[float]) -> float:
    n = min(len(left), len(right))
    lv, rv = left[:n], right[:n]
    lm, rm = sum(lv) / n, sum(rv) / n
    cov = sum((a - lm) * (b - rm) for a, b in zip(lv, rv))
    lvar = sum((a - lm) ** 2 for a in lv)
    rvar = sum((b - rm) ** 2 for b in rv)
    denom = math.sqrt(lvar * rvar)
    return cov / denom if denom else 0.0


def detrend_and_renormalize(times: list[float], forcing: list[float]) -> list[float]:
    """Remove a per-series linear trend, THEN rescale to unit RMS -- per
    this project's own documented lesson (see memory:
    feedback_decompose_nuisance_before_reporting.md) that a shared trend
    must be decomposed out before crediting a raw correlation as
    evidence of the interesting (non-trend) shared signal."""
    t = np.asarray(times)
    y = np.asarray(forcing)
    coeffs = np.polyfit(t, y, 1)
    resid = y - np.polyval(coeffs, t)
    rms = np.sqrt(np.mean(resid ** 2))
    return list(resid / rms) if rms > 0 else list(resid)


def main() -> None:
    cells = sorted((p.name for p in sweep.FEB.iterdir()
                     if p.is_dir() and NAME_RE.match(p.name)),
                    key=lambda c: decode(c))  # geographic order: lat then lon

    series = []  # (cell, times, normalized_forcing)
    skipped = []
    for cell in cells:
        path = freshest_results_csv(cell)
        if path is None:
            skipped.append(cell)
            continue
        try:
            times, forcing = read_time_and_forcing(path)
            if len(times) < 2:
                skipped.append(cell)
                continue
            series.append((cell, times, normalize_excursions(forcing)))
        except Exception as exc:
            skipped.append(f"{cell} ({exc})")

    print(f"loaded {len(series)} quads, skipped {len(skipped)}: {skipped}")

    # ---- Panel 1: spaghetti + ensemble mean --------------------------
    fig, (ax_spaghetti, ax_heat, ax_heat_dt) = plt.subplots(
        3, 1, figsize=(20, 38), gridspec_kw={"height_ratios": [1, 2.2, 2.2]})

    for cell, times, forcing in series:
        ax_spaghetti.plot(times, forcing, color="#2a78d6", alpha=0.10, linewidth=0.6)

    start = max(t[0] for _, t, _ in series)
    stop = min(t[-1] for _, t, _ in series)
    n_months = int((stop - start) * 12) + 1
    grid_t = [start + i / 12 for i in range(n_months)]
    stack = np.full((len(series), n_months), np.nan)
    for row, (cell, times, forcing) in enumerate(series):
        stack[row, :] = np.interp(grid_t, times, forcing, left=np.nan, right=np.nan)
    ensemble_mean = np.nanmean(stack, axis=0)
    ax_spaghetti.plot(grid_t, ensemble_mean, color="#0d366b", linewidth=2.2,
                       label=f"ensemble mean (n={len(series)})")
    ax_spaghetti.set_title(
        f"Hidden latent forcing (manifold, lte_results.csv col 4) — "
        f"all {len(series)} k_sst quads, normalized to unit RMS",
        fontsize=16)
    ax_spaghetti.set_xlabel("Decimal time")
    ax_spaghetti.set_ylabel("Normalized manifold (unit RMS)")
    ax_spaghetti.grid(True, alpha=0.3)
    ax_spaghetti.legend(loc="upper right", fontsize=11)

    # ---- Panel 2: 89x89 geographically-ordered correlation heatmap ---
    labels = [c for c, _, _ in series]
    n = len(series)
    corr = np.eye(n)
    for i in range(n):
        for j in range(i + 1, n):
            r = pearson_correlation(series[i][2], series[j][2])
            corr[i, j] = corr[j, i] = r

    im = ax_heat.imshow(corr, vmin=-1, vmax=1, cmap="coolwarm")
    ax_heat.set_xticks(range(n), labels, rotation=90, fontsize=4.5)
    ax_heat.set_yticks(range(n), labels, fontsize=4.5)
    ax_heat.set_title(
        "RAW pairwise Pearson correlation of normalized manifolds "
        "(geographic order: latitude then longitude)", fontsize=14)
    cbar = fig.colorbar(im, ax=ax_heat, fraction=0.03, pad=0.01)
    cbar.set_label("Pearson r", fontsize=11)

    mean_offdiag = (corr.sum() - n) / (n * n - n)

    # ---- Panel 3: DETRENDED correlation -- per this project's own
    # documented lesson (decompose a shared nuisance trend before
    # crediting a raw correlation as evidence of real shared structure)
    detrended_series = [detrend_and_renormalize(t, f) for _, t, f in series]
    corr_dt = np.eye(n)
    for i in range(n):
        for j in range(i + 1, n):
            r = pearson_correlation(detrended_series[i], detrended_series[j])
            corr_dt[i, j] = corr_dt[j, i] = r
    mean_offdiag_dt = (corr_dt.sum() - n) / (n * n - n)

    im2 = ax_heat_dt.imshow(corr_dt, vmin=-1, vmax=1, cmap="coolwarm")
    ax_heat_dt.set_xticks(range(n), labels, rotation=90, fontsize=4.5)
    ax_heat_dt.set_yticks(range(n), labels, fontsize=4.5)
    ax_heat_dt.set_title(
        "DETRENDED pairwise Pearson correlation (each series' own linear "
        "trend removed before correlating -- isolates non-trend coherence)",
        fontsize=14)
    cbar2 = fig.colorbar(im2, ax=ax_heat_dt, fraction=0.03, pad=0.01)
    cbar2.set_label("Pearson r", fontsize=11)

    # Outlier cells: lowest mean correlation to the rest of the ensemble,
    # by the detrended metric (the more honest one) -- highlight on both
    # heatmaps' tick labels so they're easy to spot.
    mean_per_cell_dt = (corr_dt.sum(axis=1) - 1) / (n - 1)
    outlier_idx = set(np.argsort(mean_per_cell_dt)[:6])
    for ax in (ax_heat, ax_heat_dt):
        for i, lbl in enumerate(ax.get_yticklabels()):
            if i in outlier_idx:
                lbl.set_color("#d03b3b")
                lbl.set_fontweight("bold")
        for i, lbl in enumerate(ax.get_xticklabels()):
            if i in outlier_idx:
                lbl.set_color("#d03b3b")
                lbl.set_fontweight("bold")

    outlier_names = [labels[i] for i in sorted(outlier_idx, key=lambda i: mean_per_cell_dt[i])]
    fig.text(0.5, 0.005,
              f"Mean off-diagonal correlation: RAW={mean_offdiag:.3f}  "
              f"DETRENDED={mean_offdiag_dt:.3f}  (raw is inflated by a shared "
              f"secular trend common to all quads -- detrended isolates the "
              f"non-trend coherence). Lowest-coherence outliers (red labels): "
              f"{', '.join(outlier_names)}",
              ha="center", fontsize=12)

    plt.tight_layout(rect=(0, 0.015, 1, 1))
    out = sweep.HERE / "manifold_comparison_all_quads.png"
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")
    print(f"mean off-diagonal correlation: RAW={mean_offdiag:.4f}  DETRENDED={mean_offdiag_dt:.4f}")
    print(f"outliers (detrended, lowest-6): {outlier_names}")


if __name__ == "__main__":
    main()
