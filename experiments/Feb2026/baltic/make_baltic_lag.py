#!/usr/bin/env python3
"""Time-shifted Baltic MSL with the most negative correlation to NAO, for 2-leg fits.

    python3 make_baltic_lag.py [--lag MONTHS]      default: scan -3..+3 months, pick the best

The shift is band-limited (Fourier phase shift of the quadratic-detrended series,
reflect-padded, trend shifted with it), so a fractional lag does not smooth the
series the way linear interpolation does, which would otherwise inflate |r| at
half-month lags. Positive lag: Baltic lags NAO, and the value written at date t
is Baltic at t + lag. Writes baltic_lag<lag>.dat (e.g. baltic_lag0p50.dat) on
baltic.dat's own dates.
"""
import argparse
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lag", type=float, default=None)
    a = ap.parse_args()

    b = np.loadtxt(HERE / "baltic.dat")
    n = np.loadtxt(HERE.parent / "nao" / "nao.dat")
    kb = np.rint(b[:, 0] * 12).astype(int)
    assert np.all(np.diff(kb) == 1), "baltic.dat is not contiguous monthly"
    coef = np.polyfit(kb, b[:, 1], 2)
    x0 = b[:, 1] - np.polyval(coef, kb)
    N = len(x0)
    freq = np.fft.rfftfreq(2 * N)
    spec = np.fft.rfft(np.r_[x0, x0[::-1]])

    def shifted(lag):
        return np.fft.irfft(spec * np.exp(2j * np.pi * freq * lag), 2 * N)[:N]

    nd = dict(zip(np.rint(n[:, 0] * 12).astype(int), n[:, 1]))
    ok = np.isin(kb, list(nd))
    ny = np.array([nd[i] for i in kb[ok]])

    def cc(lag):
        return np.corrcoef(shifted(lag)[ok], ny)[0, 1]

    lag = a.lag
    if lag is None:
        lags = np.round(np.arange(-3.0, 3.001, 0.05), 2)
        lag = float(lags[np.argmin([cc(l) for l in lags])])
    y = shifted(lag) + np.polyval(coef, kb + lag)
    name = f"baltic_lag{lag:.2f}".replace("-", "m").replace(".", "p") + ".dat"
    with open(HERE / name, "w") as f:
        for t, v in zip(b[:, 0], y):
            f.write(f"{t:.4f}  {v:.4f}\n")
    print(f"lag {lag:+.2f} mo: r with NAO {cc(lag):+.4f} (unshifted {cc(0.0):+.4f}) -> {name}")


if __name__ == "__main__":
    main()
