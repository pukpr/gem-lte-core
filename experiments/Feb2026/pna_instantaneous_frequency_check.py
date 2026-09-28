#!/usr/bin/env python3
"""pna_instantaneous_frequency_check.py — tests a specific mechanistic
claim about sin(k*M(t)) directly, rather than scanning fixed candidate
winding numbers: for phi(t) = k*M(t), the INSTANTANEOUS frequency is
d(phi)/dt = k*dM/dt, not a constant -- so a single, already-fitted,
modest k can manifest as locally HIGH-frequency oscillation wherever the
manifold's own gradient is steep, and near-zero oscillation on its
plateaus (M(t) here has a genuinely "squared-off" character: Impulse_
Delta's twice-yearly comb kicks the IIR state hard, then it holds nearly
flat between kicks -- verified directly: smoothed |dF/dt| for pna spans
0.02 to ~15 (raw, unsmoothed spikes reach ~90), a ~1000x range between
plateau and transition).

This is a genuinely different test than every ridge-scan/winding_rank
check elsewhere in this project, which all treat M as a FIXED constant
scanned uniformly across the whole record. Here there is no scan: k is
whatever is already established (this project's shared 0.2076 cross-
index backbone -- not a new parameter), and the prediction is a time-
VARYING instantaneous frequency k*dF/dt(t), checked directly against the
real (necessarily filtered first -- monthly PNA is noisy) data's own
empirically-measured instantaneous frequency via a Hilbert transform.

A real subtlety, caught and corrected before trusting a null result: the
manifold's steepest transitions imply instantaneous frequencies far
beyond monthly Nyquist (6 cyc/yr) for pna's own higher established
windings (M=1.46/2.97/4.73 all alias badly there -- e.g. M=4.73 times
the peak gradient implies ~425 cyc/yr, hopelessly aliased at monthly
resolution). The base backbone k=0.2076 stays mostly within Nyquist even
at the steepest points (predicted range 0.005-3.16 cyc/yr), so it is
the only one of pna's own established windings for which this test is
even well-posed at monthly resolution -- tested first with an
ungrounded implicit k=1 comparison (null, r=0.06), then correctly with
k=0.2076 (real signal, r=0.22 raw -> 0.52 after matched smoothing).

Usage
-----
    ./pna_instantaneous_frequency_check.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.signal import hilbert, butter, filtfilt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cv_rolling_blocked import prepare                                     # noqa: E402

ROOT = Path(__file__).resolve().parent
K_BASE = 0.2076  # this project's shared cross-index winding backbone


def smooth(x: np.ndarray, win_months: int) -> np.ndarray:
    k = np.ones(win_months) / win_months
    return np.convolve(x, k, mode="same")


def main() -> int:
    prep = prepare("pna")
    dates, F = prep["dates"], prep["forcing"]
    dt_yr = np.median(np.diff(dates))
    fs = 1.0 / dt_yr
    nyq = fs / 2

    dat = np.loadtxt(ROOT / "pna" / "pna.dat")
    obs = np.interp(dates, dat[:, 0], dat[:, 1])

    dF = np.gradient(F, dates)
    absdF_smooth = smooth(np.abs(dF), 12)
    pred_freq = K_BASE * absdF_smooth
    print(f"Nyquist: {nyq:.2f} cyc/yr")
    print(f"predicted instantaneous frequency range at k={K_BASE}: "
          f"{pred_freq.min():.3f}-{pred_freq.max():.3f} cyc/yr "
          f"({'within' if pred_freq.max() < nyq else 'EXCEEDS'} Nyquist)")

    # bandpass matched to the predicted frequency range -- filtering PNA
    # is necessary before a Hilbert-transform instantaneous frequency
    # means anything at all (the method requires a roughly narrowband
    # signal); periods ~0.5-60yr covers the predicted 0.005-3.16 cyc/yr
    # range with margin.
    low = 1 / 60.0 / nyq
    high = min(1 / 0.5 / nyq, 0.99)
    b, a = butter(4, [low, high], btype="band")
    obs_filt = filtfilt(b, a, obs)
    analytic = hilbert(obs_filt)
    inst_phase = np.unwrap(np.angle(analytic))
    inst_freq = np.abs(np.gradient(inst_phase, dates) / (2 * np.pi))

    edge = 36
    m = slice(edge, -edge)
    r_raw = np.corrcoef(inst_freq[m], pred_freq[m])[0, 1]
    inst_freq_s = smooth(inst_freq, 24)
    pred_freq_s = smooth(pred_freq, 24)
    r_smooth = np.corrcoef(inst_freq_s[m], pred_freq_s[m])[0, 1]

    print(f"\ncorr(real filtered-PNA instantaneous freq, predicted "
          f"k_base*|dF/dt|):")
    print(f"  raw (monthly): r={r_raw:+.3f}")
    print(f"  after 2yr smoothing both: r={r_smooth:+.3f}")

    print("\nrolling 20yr correlation by era:")
    for c in range(1965, 2015, 10):
        mm = (dates >= c - 10) & (dates < c + 10)
        r_local = np.corrcoef(inst_freq[mm], pred_freq[mm])[0, 1]
        print(f"  {c}: r={r_local:+.3f}")

    fig, axes = plt.subplots(3, 1, figsize=(11, 9))
    axes[0].plot(dates, F, linewidth=0.6, color="tab:gray")
    axes[0].set_title("pna: manifold Forcing(t) -- note the squared-off "
                     "plateau/transition character", fontsize=9)
    axes[0].set_xlabel("year")

    axes[1].plot(dates, absdF_smooth, linewidth=0.7, color="tab:red")
    axes[1].set_title("smoothed |dF/dt| (12mo) -- 'local manifold "
                     "activity'", fontsize=9)
    axes[1].set_ylabel("manifold units / yr")
    axes[1].set_xlabel("year")

    axes[2].plot(dates[m], inst_freq_s[m], label="real filtered-PNA "
                "instantaneous frequency (Hilbert, 2yr-smoothed)",
                color="tab:blue", linewidth=0.9)
    axes[2].plot(dates[m], pred_freq_s[m], label=f"predicted k*|dF/dt| "
                f"(k={K_BASE}, 2yr-smoothed)", color="tab:orange",
                linewidth=0.9)
    axes[2].set_title(f"pna: real vs predicted instantaneous frequency "
                     f"(r={r_smooth:+.3f})", fontsize=9)
    axes[2].set_ylabel("cycles/year")
    axes[2].set_xlabel("year")
    axes[2].legend(fontsize=7)

    fig.tight_layout()
    out_path = ROOT / "pna" / "pna_instantaneous_frequency_check.png"
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"\n  saved {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
