#!/usr/bin/env python3
"""power_spectrum.py -- standard calendar-time power spectrum (periodogram)
of an index's lte_results.csv, Data (col 3) vs Model (col 2) against the
col 1 time base, with an AR(1) red-noise significance floor and a
band-power breakdown.

This is deliberately a DIFFERENT axis from winding_scalogram.py: that tool
transforms against the model's own Forcing variable (M winds around
Forcing(t), not calendar time), which is the right tool for asking "does
this index's own fitted winding numbers show up as persistent structure."
This script instead asks the plainer question "what does a standard
FFT-based power spectrum in cycles/year look like, and how much of the
Data's real power does the Model reproduce, band by band" -- a good
complement, not a replacement (confirmed directly this session: a rich,
"continuous winding ridge" forcing-domain scalogram does NOT imply uniform
fidelity across calendar-time frequency bands -- kN020_E050's model
reproduced ~65-68% of the real power at interannual/annual timescales but
only ~5-6% at the fastest, sub-5-month band, despite dense winding-domain
structure throughout).

AR(1) red noise (not white noise) is used for the significance floor,
matching winding_rank.ar1_floor's and winding_scalogram.noise_floor's own
convention -- monthly climate series are strongly autocorrelated, and a
white-noise floor is anti-conservative for this kind of data (verified
directly in winding_scalogram.py's own docstring: pure AR1 red noise with
no real signal produces a false ~4-bit "bright ridge" against a
white-noise floor).

Usage:
    power_spectrum.py <index_dir> [--fs 12] [--n-surr 200] [--outdir .]

<index_dir> must contain lte_results.csv (col 1 = date, col 2 = model,
col 3 = data), e.g. as produced by any lt.exe run in this project.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BANDS = [(0.05, 0.30), (0.30, 0.70), (0.70, 1.30), (1.30, 2.50),
         (2.50, 4.00), (4.00, 6.00)]


def detrend_window(x: np.ndarray) -> np.ndarray:
    """Linear detrend + Hann window, power-normalized so the window itself
    doesn't bias the overall power level."""
    t = np.arange(len(x))
    coef = np.polyfit(t, x, 1)
    xd = x - np.polyval(coef, t)
    win = np.hanning(len(x))
    return xd * win / np.sqrt(np.mean(win ** 2))


def ar1_coef(x: np.ndarray) -> float:
    return float(np.corrcoef(x[:-1], x[1:])[0, 1])


def periodogram(x: np.ndarray, fs: float) -> tuple[np.ndarray, np.ndarray]:
    n = len(x)
    freqs = np.fft.rfftfreq(n, d=1.0 / fs)
    power = np.abs(np.fft.rfft(detrend_window(x))) ** 2 / n
    return freqs, power


def ar1_floor(x: np.ndarray, fs: float, n_surr: int = 200,
              seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """AR(1) red-noise surrogate floor, matched to x's own lag-1
    autocorrelation. Returns (mean, p99) periodograms, scaled to match
    x's own in-band total power so the floor is a fair reference curve
    rather than an arbitrary-amplitude null."""
    n = len(x)
    t = np.arange(n)
    xd = x - np.polyval(np.polyfit(t, x, 1), t)
    rho = ar1_coef(xd)
    s = float(np.sqrt(max(1e-9, 1.0 - rho * rho)))
    rng = np.random.default_rng(seed)
    freqs = np.fft.rfftfreq(n, d=1.0 / fs)
    surrogates = np.empty((n_surr, len(freqs)))
    for k in range(n_surr):
        e = rng.standard_normal(n)
        z = np.zeros(n)
        for i in range(1, n):
            z[i] = rho * z[i - 1] + s * e[i]
        surrogates[k] = np.abs(np.fft.rfft(detrend_window(z))) ** 2 / n
    mean = surrogates.mean(axis=0)
    p99 = np.percentile(surrogates, 99, axis=0)
    _, p_real = periodogram(x, fs)
    scale = np.sum(p_real[1:]) / max(1e-12, np.sum(mean[1:]))
    return mean * scale, p99 * scale, rho


def find_significant_peaks(freqs: np.ndarray, power: np.ndarray,
                            floor_p99: np.ndarray, f_min: float = 0.05
                            ) -> list[tuple[float, float, float]]:
    """Local maxima exceeding the AR1 p99 floor. Returns
    (freq, power, power/floor_mean_at_that_bin) sorted by power desc."""
    peaks = []
    for i in range(2, len(freqs) - 2):
        if freqs[i] < f_min:
            continue
        if power[i] > power[i - 1] and power[i] > power[i + 1] and power[i] > floor_p99[i]:
            peaks.append((freqs[i], power[i]))
    peaks.sort(key=lambda p: -p[1])
    return peaks


def band_power_table(freqs: np.ndarray, p_data: np.ndarray,
                      p_model: np.ndarray) -> list[dict]:
    total_data = p_data[1:].sum()
    total_model = p_model[1:].sum()
    rows = []
    for lo, hi in BANDS:
        m = (freqs >= lo) & (freqs < hi)
        pd = p_data[m].sum()
        pm = p_model[m].sum()
        rows.append(dict(lo=lo, hi=hi, data_frac=pd / total_data if total_data else 0.0,
                          model_frac=pm / total_model if total_model else 0.0,
                          ratio=pm / pd if pd > 0 else float("nan")))
    return rows


def run(index_dir: Path, fs: float, n_surr: int, outdir: Path,
        linear: bool = False) -> None:
    csv_path = index_dir / "lte_results.csv"
    raw = np.loadtxt(csv_path, delimiter=",")
    dates, model, data = raw[:, 0], raw[:, 1], raw[:, 2]
    n = len(dates)
    name = index_dir.name

    dt = float(np.median(np.diff(dates)))
    print(f"[{name}] n={n}  date range {dates.min():.1f}-{dates.max():.1f}  "
          f"sampling dt={dt:.5f}yr (~{round(1/dt)} samples/yr)")

    freqs, p_data = periodogram(data, fs)
    _, p_model = periodogram(model, fs)
    floor_mean, floor_p99, rho = ar1_floor(data, fs, n_surr=n_surr)
    print(f"[{name}] data AR(1) rho = {rho:.3f}  ({n_surr} surrogates)")

    peaks = find_significant_peaks(freqs, p_data, floor_p99)
    total_power = p_data[1:].sum()
    print(f"\n[{name}] {len(peaks)} AR1-significant Data peaks "
          f"(>p99 red-noise floor):")
    print(f"  {'freq(cyc/yr)':>13s} {'period':>12s} {'power':>10s} "
          f"{'frac of total':>14s} {'model/data':>11s}")
    for f, p in peaks[:20]:
        i = int(np.argmin(np.abs(freqs - f)))
        period_str = f"{12.0/f:.2f} mo" if f > 0 else "inf"
        ratio = p_model[i] / p if p > 0 else float("nan")
        print(f"  {f:13.4f} {period_str:>12s} {p:10.4g} {p/total_power:14.3f} {ratio:11.3f}")

    print(f"\n[{name}] band-power breakdown (data's share of its own total "
          f"variance vs model/data power ratio):")
    print(f"  {'band (cyc/yr)':>16s} {'data frac':>10s} {'model frac':>11s} {'model/data':>11s}")
    for row in band_power_table(freqs, p_data, p_model):
        print(f"  {row['lo']:6.2f}-{row['hi']:<6.2f}   {row['data_frac']:10.3f} "
              f"{row['model_frac']:11.3f} {row['ratio']:11.3f}")

    to_plot = (lambda p: p) if linear else (lambda p: np.log2(np.maximum(p, 1e-12)))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    for ax, power, label, color in [(ax1, p_data, "Data", "tab:blue"),
                                     (ax1, p_model, "Model", "tab:orange")]:
        ax.plot(freqs, to_plot(power), color=color, lw=1.0, label=label)
    ax1.plot(freqs, to_plot(floor_mean), color="0.4",
              lw=1.0, linestyle="--", label="AR1 red-noise floor (mean)")
    ax1.plot(freqs, to_plot(floor_p99), color="0.4",
              lw=0.8, linestyle=":", label="AR1 floor (p99)")
    for f, _ in peaks[:10]:
        ax1.axvline(f, color="red", lw=0.5, alpha=0.4)
    ax1.set_ylabel("power" if linear else "log2 power")
    ax1.set_title(f"{name}: calendar-time power spectrum, Data vs Model "
                  f"(col1 time base, col3 vs col2 of lte_results.csv)")
    ax1.legend(fontsize=8)

    ratio = np.where(p_data > 1e-12, p_model / np.maximum(p_data, 1e-12), np.nan)
    ax2.plot(freqs, ratio, color="tab:green", lw=1.0)
    ax2.axhline(1.0, color="0.6", lw=0.8)
    ax2.set_ylim(0, 2)
    ax2.set_xlabel("frequency (cycles/year)")
    ax2.set_ylabel("model / data power ratio")
    ax2.set_title("model fidelity by frequency (1.0 = perfect match)")

    fig.tight_layout()
    outdir.mkdir(parents=True, exist_ok=True)
    suffix = "_power_spectrum_linear.png" if linear else "_power_spectrum.png"
    out_path = outdir / f"{name}{suffix}"
    fig.savefig(out_path, dpi=140)
    print(f"\n[{name}] saved {out_path}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("index_dir", type=Path,
                     help="directory containing lte_results.csv")
    ap.add_argument("--fs", type=float, default=12.0,
                     help="samples per year (default 12, monthly)")
    ap.add_argument("--n-surr", type=int, default=200,
                     help="number of AR(1) surrogates for the noise floor")
    ap.add_argument("--outdir", type=Path, default=None,
                     help="where to save the plot (default: index_dir)")
    ap.add_argument("--linear", action="store_true",
                     help="plot the top panel's power on a linear scale "
                          "instead of the default log2 scale")
    args = ap.parse_args()
    outdir = args.outdir if args.outdir is not None else args.index_dir
    run(args.index_dir, args.fs, args.n_surr, outdir, linear=args.linear)


if __name__ == "__main__":
    main()
