#!/usr/bin/env python3
"""spectral_slope_check.py — does the "unexplained long-period" residual
get systematically redder (more power at low frequency) than what's
already built into the forcing manifold, consistent with the 1/omega^2
material-coordinate compliance derived in the two-layer LTE review
(k_i' = k_i/omega_i^2), on top of everything the pipeline's own
construction of Forcing already does?

Two things have to be handled carefully, both raised directly by the user
rather than assumed:

1. The index .dat files (e.g. baltic.dat) are NOT raw sea level. They have
   already had the strong diurnal/semidiurnal astronomical tide and the
   annual/semi-annual seasonal cycle removed upstream, before this
   pipeline ever sees them -- what's left is specifically the long-period
   residual the LTE manifold is trying to explain. Diurnal/semidiurnal
   periods (<1 day) are moot for this check regardless, since these series
   are sampled monthly (well above their Nyquist already excludes them);
   what matters here is that the annual (1 cycle/yr) and semi-annual
   (2 cycles/yr) bands are a NOTCH, not a genuine absence of variability --
   fitting a spectral slope through them (or judging "how red" the
   spectrum is by comparing power at 1/yr against low frequencies) would
   be comparing against an artificially-suppressed reference, not a fair
   one. Both bands are excluded (with a guard width) from every slope fit
   and shown hatched, not silently dropped, on the plot.

2. `Forcing` (this project's own "manifold" variable) is not raw
   astronomical forcing -- by the time cv_rolling_blocked.prepare()
   returns it, it has already passed through impulse_delta (a monthly
   Dirac-comb gate -- broadband spectral spreading) and, crucially, `iir`,
   which IS a literal linear (leaky) integrator: y[i] = x[i] + Mem*y[i-1]
   - ..., Mem = 1-ma close to 1. An integrator has its own low-frequency
   amplification built in, independent of any real ocean physics. So a
   red `data` spectrum only supports the stratification/compliance
   argument if `data` is REDDER than `Forcing` already is -- if data's
   slope just matches Forcing's, the pipeline's own IIR step could explain
   all of it with no additional ocean-side compliance required. This
   script therefore compares THREE spectra, in pipeline order:

     tide_raw            (raw Doodson harmonic sum -- no impulse/IIR/Bessel)
       -> forcing_pre_bessel  (+ impulse comb + IIR integrator)
         -> forcing          (+ Bessel nonlinearity -- what modes actually see)
           -> data           (the real, already-notched residual)

   A slope that gets progressively steeper at each stage, with data's
   final slope clearly steeper than forcing's, is the signature the
   1/omega^2 argument actually predicts. A slope that's already fully
   explained by tide_raw -> forcing is not evidence of anything beyond
   the pipeline's own construction.

Method: a single-taper (Tukey) periodogram of each (detrended) series,
log-frequency-binned to tame variance, with a power-law (log-log OLS)
slope fit over a chosen band, excluding annual/semi-annual guard bands.

Usage
-----
    ./spectral_slope_check.py baltic
    ./spectral_slope_check.py baltic nino4 --long-period-band 3 40
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.signal.windows import tukey
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cv_rolling_blocked import prepare, DEFAULT_INDICES  # noqa: E402

ROOT = Path(__file__).resolve().parent
ANNUAL_NOTCH = [(1.0, 0.15), (2.0, 0.15)]  # (cycles/yr, half-width) -- annual + semi-annual


def periodogram(dates: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Single-taper (Tukey) periodogram of a detrended series, uniformly
    sampled at dt = median diff(dates) years. Returns (freq cycles/yr,
    power), positive frequencies only, DC dropped."""
    n = len(x)
    dt = float(np.median(np.diff(dates)))
    x = x - np.polyval(np.polyfit(dates, x, 1), dates)  # remove linear trend
    win = tukey(n, alpha=0.1)
    xw = x * win
    # normalize so window tapering doesn't bias absolute power (only
    # relative/slope comparisons are used, but keep it honest anyway)
    xw *= np.sqrt(n / np.sum(win ** 2))
    spec = np.fft.rfft(xw)
    freq = np.fft.rfftfreq(n, d=dt)  # cycles/yr
    power = (np.abs(spec) ** 2) * dt / n
    return freq[1:], power[1:]  # drop DC


def log_bin(freq: np.ndarray, power: np.ndarray, n_bins: int = 40
            ) -> tuple[np.ndarray, np.ndarray]:
    """Geometric-spaced binning of the periodogram (mean power per bin) --
    standard variance reduction for fitting a power-law slope to a noisy
    periodogram."""
    edges = np.geomspace(freq[0], freq[-1], n_bins + 1)
    f_out, p_out = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (freq >= lo) & (freq < hi)
        if np.any(sel):
            f_out.append(np.exp(np.mean(np.log(freq[sel]))))
            p_out.append(np.mean(power[sel]))
    return np.array(f_out), np.array(p_out)


def notch_mask(freq: np.ndarray, notches: list[tuple[float, float]]
              ) -> np.ndarray:
    """True where freq is CLEAR of every (center, half-width) notch band."""
    clear = np.ones_like(freq, dtype=bool)
    for center, half in notches:
        clear &= np.abs(freq - center) > half
    return clear


def fit_slope(freq: np.ndarray, power: np.ndarray, band: tuple[float, float]
              ) -> dict:
    """OLS slope of log(power) vs log(freq) within [1/band[1], 1/band[0]]
    cycles/yr (band given as a PERIOD range in years, e.g. (3, 40) =
    3-40 year periods), excluding the annual/semi-annual notches."""
    f_lo, f_hi = 1.0 / band[1], 1.0 / band[0]
    sel = (freq >= f_lo) & (freq <= f_hi) & notch_mask(freq, ANNUAL_NOTCH)
    if sel.sum() < 4:
        return dict(slope=float("nan"), intercept=float("nan"), n=int(sel.sum()))
    x, y = np.log(freq[sel]), np.log(power[sel])
    slope, intercept = np.polyfit(x, y, 1)
    resid = y - (slope * x + intercept)
    r2 = 1.0 - np.sum(resid ** 2) / np.sum((y - y.mean()) ** 2)
    return dict(slope=float(slope), intercept=float(intercept), n=int(sel.sum()),
                r2=float(r2))


def run(idx: str, long_band: tuple[float, float], short_band: tuple[float, float],
        outdir: Path | None) -> None:
    prep = prepare(idx)
    dates = prep["dates"]
    series = {
        "tide_raw (pre impulse/IIR/Bessel)": prep["tide_raw"],
        "forcing_pre_bessel (post-IIR)": prep["forcing_pre_bessel"],
        "forcing (post-Bessel, final manifold)": prep["forcing"],
        "data (real, already annual/semi-annual-filtered residual)": prep["data_raw"],
    }

    print(f"-- {idx}: spectral slope by pipeline stage --")
    print(f"  long-period band: {long_band[0]:.0f}-{long_band[1]:.0f} yr periods "
          f"(annual/semi-annual notches excluded, half-width "
          f"{ANNUAL_NOTCH[0][1]:.2f} cyc/yr)")
    print(f"  short-period band: {short_band[0]:.0f}-{short_band[1]:.0f} yr periods")

    results = {}
    binned = {}
    for name, x in series.items():
        freq, power = periodogram(dates, x)
        fb, pb = log_bin(freq, power)
        binned[name] = (fb, pb)
        long_fit = fit_slope(fb, pb, long_band)
        short_fit = fit_slope(fb, pb, short_band)
        results[name] = (long_fit, short_fit)
        print(f"  {name}:")
        print(f"    long-period slope:  {long_fit['slope']:+.2f}  "
              f"(R^2={long_fit.get('r2', float('nan')):.2f}, n={long_fit['n']} bins)")
        print(f"    short-period slope: {short_fit['slope']:+.2f}  "
              f"(R^2={short_fit.get('r2', float('nan')):.2f}, n={short_fit['n']} bins)")

    names = list(series)
    data_long = results[names[-1]][0]["slope"]
    forcing_long = results[names[2]][0]["slope"]
    tide_long = results[names[0]][0]["slope"]
    print(f"\n  verdict (long-period band, more negative = redder):")
    print(f"    tide_raw slope:        {tide_long:+.2f}")
    print(f"    forcing (final) slope: {forcing_long:+.2f}  "
          f"(pipeline's own IIR/Bessel steepening vs tide_raw: "
          f"{forcing_long - tide_long:+.2f})")
    print(f"    data slope:            {data_long:+.2f}  "
          f"(additional steepening beyond forcing: {data_long - forcing_long:+.2f})")
    if data_long < forcing_long - 0.3:
        print("    -> data is meaningfully redder than the forcing manifold "
              "itself: consistent with genuine additional (ocean-side) "
              "low-frequency compliance, not just inherited pipeline shape.")
    elif data_long > forcing_long + 0.3:
        print("    -> data is actually LESS red than the forcing manifold: "
              "no evidence here of extra low-frequency compliance beyond "
              "what the pipeline's own construction already imposes.")
    else:
        print("    -> data's slope is close to forcing's: the red character, "
              "such as it is, looks largely inherited from the forcing "
              "construction (notably the IIR integrator), not demonstrated "
              "as an additional ocean-side effect.")

    plot_spectra(idx, binned, long_band, short_band,
                out_path=(outdir if outdir is not None else ROOT / idx) /
                "spectral_slope_check.png")


def plot_spectra(idx: str, binned: dict, long_band, short_band,
                 out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 6))
    colors = ["0.6", "tab:blue", "tab:orange", "tab:red"]
    for (name, (f, p)), color in zip(binned.items(), colors):
        ax.plot(f, p, "o-", color=color, markersize=3, linewidth=1,
               label=name)
    ax.set_xscale("log")
    ax.set_yscale("log")
    for center, half in ANNUAL_NOTCH:
        ax.axvspan(center - half, center + half, color="0.85", zorder=0)
    ax.axvspan(1 / long_band[1], 1 / long_band[0], color="tab:green",
              alpha=0.08, zorder=0, label="long-period fit band")
    ax.set_xlabel("frequency (cycles/year, log)")
    ax.set_ylabel("power (log)")
    ax.set_title(f"{idx}: spectrum by pipeline stage "
                f"(grey bands = excluded annual/semi-annual notches)",
                fontsize=10)
    ax.legend(fontsize=7, loc="lower left")
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"\n  saved {out_path}")


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("indices", nargs="*", default=DEFAULT_INDICES)
    ap.add_argument("--long-period-band", type=float, nargs=2, default=(3.0, 40.0),
                     metavar=("MIN_YR", "MAX_YR"),
                     help="period range (years) treated as 'unexplained "
                          "long-period' for the primary slope fit (default 3-40)")
    ap.add_argument("--short-period-band", type=float, nargs=2, default=(0.6, 2.5),
                     metavar=("MIN_YR", "MAX_YR"),
                     help="comparison period range (default 0.6-2.5 yr, "
                          "excluding the annual/semi-annual notches within it)")
    ap.add_argument("--outdir", type=Path, default=None)
    args = ap.parse_args()

    for idx in args.indices:
        run(idx, tuple(args.long_period_band), tuple(args.short_period_band),
           args.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
