#!/usr/bin/env python3
"""Multivariate IAAFT surrogates for the surrogate parsimony test (see DECLARATION.md).
Common span 1880.0-2022.083 (1706 months). The same random phase rotation is applied to all three
series in a draw (keeps cross-spectra); then each series alternates 200 times between imposing its
original Fourier amplitudes (keeping current phases) and rank-matching its original values.
Writes surr/<index>_s00.dat (real, truncated) ... <index>_s19.dat, plus planted-control series."""
import numpy as np
from pathlib import Path
H = Path(__file__).resolve().parent; S = H / "surr"; S.mkdir(exist_ok=True)
IDX = ("amo", "pdo", "nino4"); N = 1706; NDRAW = 19; SEED = 20261006
def load(i):
    a = np.loadtxt(H / "seeds" / i / f"{i}.dat")[:N]; return a[:, 0], a[:, 1]
def iaaft_multi(X, rng, iters=200):
    n = X.shape[1]; F = np.fft.rfft(X, axis=1); amp = np.abs(F); srt = np.sort(X, axis=1)
    rot = rng.uniform(0, 2 * np.pi, F.shape[1]); rot[0] = 0
    if n % 2 == 0: rot[-1] = 0
    Y = np.fft.irfft(amp * np.exp(1j * (np.angle(F) + rot)), n=n, axis=1)
    for _ in range(iters):
        for k in range(X.shape[0]):                       # rank-match values
            Y[k, np.argsort(Y[k])] = srt[k]
        G = np.fft.rfft(Y, axis=1); Y = np.fft.irfft(amp * np.exp(1j * np.angle(G)), n=n, axis=1)
    for k in range(X.shape[0]): Y[k, np.argsort(Y[k])] = srt[k]
    return Y
def write(path, t, x): np.savetxt(path, np.column_stack([t, x]), fmt=["%.6f", "%.6f"], delimiter="\t")
if __name__ == "__main__":
    T = {}; X = []
    for i in IDX:
        t, x = load(i); T[i] = t; X.append(x)
    assert all(np.allclose(T[i], T["amo"]) for i in IDX), "date columns differ"
    X = np.array(X); rng = np.random.default_rng(SEED)
    for k, i in enumerate(IDX): write(S / f"{i}_s00.dat", T[i], X[k])
    for d in range(1, NDRAW + 1):
        Y = iaaft_multi(X, rng)
        for k, i in enumerate(IDX): write(S / f"{i}_s{d:02d}.dat", T[i], Y[k])
    # positive control: real fit's model output + IAAFT of its residual (phase-locked by construction)
    for i in IDX:
        a = np.loadtxt(H / "seeds" / i / "lte_results.csv", delimiter=",")[:N]; m = a[:, 1]; res = X[IDX.index(i)] - m
        r = iaaft_multi(res[None, :], np.random.default_rng(SEED + 1))[0]; write(S / f"{i}_planted.dat", T[i], m + r)
    # quality report
    for k, i in enumerate(IDX):
        Y = np.loadtxt(S / f"{i}_s01.dat")[:, 1]
        pX, pY = np.abs(np.fft.rfft(X[k]))**2, np.abs(np.fft.rfft(Y))**2
        print(f"{i}: spectrum CC real vs s01 {np.corrcoef(np.log(pX[1:]), np.log(pY[1:]))[0,1]:.3f}; values identical {np.allclose(np.sort(X[k]), np.sort(Y))}; time-domain CC {np.corrcoef(X[k], Y)[0,1]:+.3f}")
    Y = np.array([np.loadtxt(S / f"{i}_s01.dat")[:, 1] for i in IDX])
    print("cross-correlation real:", np.round(np.corrcoef(X)[np.triu_indices(3, 1)], 3), " s01:", np.round(np.corrcoef(Y)[np.triu_indices(3, 1)], 3))
