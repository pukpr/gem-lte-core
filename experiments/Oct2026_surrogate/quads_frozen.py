#!/usr/bin/env python3
"""DECLARATION_QUADS.md: frozen-manifold low-winding test across the 89 quads (no fitting)."""
import numpy as np, sys
from pathlib import Path
from math import comb
import surrogates as SG
H = Path(__file__).resolve().parent; E = H.parent; L = E / "Oct2026_layer"
FG = np.round(np.arange(0.5, 40.0001, 0.05), 2); LOW = FG <= 6.0
SHIFTS = (5, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47, 50, 53, 56, 59)
a = np.loadtxt(H / "seeds" / "amo" / "lte_results.csv", delimiter=","); tA, FA = a[:, 0], a[:, 3]
m = (tA >= 1950 - 1e-6); tF, F2 = tA[m], FA[m]; n = len(tF)
E_ = np.exp(1j * 2 * np.pi * np.outer(FG, F2))              # winding basis on the frozen manifold
def prep(y): x = np.arange(len(y)); y = y - np.polyval(np.polyfit(x, y, 1), x); return y / y.std()
def stats(y, Fbasis=E_):
    P = np.abs(Fbasis @ y)**2; L1, L2 = P.sum(), np.sqrt((P**2).sum()); k = len(P)
    return P[LOW].sum() / P.sum(), (np.sqrt(k) - L1 / L2) / (np.sqrt(k) - 1), FG[np.argmax(P)]
def shifted_basis(s): return np.exp(1j * 2 * np.pi * np.outer(FG, np.roll(F2, 12 * s)))
SB = [shifted_basis(s) for s in SHIFTS]
def binom_p(h, N, p=0.05): return sum(comb(N, k) * p**k * (1 - p)**(N - k) for k in range(h, N + 1))
def quad_series(q):
    d = np.loadtxt(L / q / f"{q}.dat"); idx = [np.argmin(np.abs(d[:, 0] - x)) for x in tF]
    assert np.allclose(d[idx, 0], tF, atol=0.02); return prep(d[idx, 1])
if __name__ == "__main__":
    rng = np.random.default_rng(7)
    # positive control first
    hits = [0, 0]
    for k in range(50):
        e = np.zeros(n); w = rng.normal(size=n)
        for t in range(1, n): e[t] = 0.9 * e[t - 1] + w[t]
        y = prep(0.3 * np.sin(2 * np.pi * 4.0 * F2) / np.sin(2 * np.pi * 4.0 * F2).std() + e / e.std())
        s0 = stats(y)[0]
        n1 = [stats(prep(SG.iaaft_multi(y[None, :], np.random.default_rng(1000 + k * 19 + j), iters=100)[0]))[0] for j in range(19)]
        n2 = [stats(y, B)[0] for B in SB]
        hits[0] += s0 > max(n1); hits[1] += s0 > max(n2)
    print(f"POSITIVE CONTROL (planted k=4.0, SNR 0.3, AR1 0.9 noise): hit rate N1 {hits[0]}/50, N2 {hits[1]}/50 (need >= 40)")
    quads = sorted(p.name for p in L.iterdir() if p.is_dir() and p.name.startswith("k") and (p / f"{p.name}.dat").exists())
    rec = []
    for q in quads:
        y = quad_series(q); r = stats(y)
        sur = [prep(SG.iaaft_multi(y[None, :], np.random.default_rng(5000 + 19 * quads.index(q) + j), iters=200)[0]) for j in range(19)]
        n1 = np.array([stats(s) for s in sur]); n2 = np.array([stats(y, B) for B in SB])
        neg = n1[0, 0] > n1[1:, 0].max()
        rec.append((q, r, n1, n2, neg))
    N = len(rec)
    for lab, j in (("S_low (primary)", 0), ("S_H Hoyer", 1)):
        h1 = sum(x[1][j] > x[2][:, j].max() for x in rec); h2 = sum(x[1][j] > x[3][:, j].max() for x in rec)
        p1 = np.median([100 * (x[2][:, j] < x[1][j]).mean() for x in rec]); p2 = np.median([100 * (x[3][:, j] < x[1][j]).mean() for x in rec])
        print(f"{lab}: hits vs N1 (data IAAFT) {h1}/{N} (binomial p {binom_p(h1, N):.1e}), median percentile {p1:.0f} | hits vs N2 (shifted clock) {h2}/{N} (p {binom_p(h2, N):.1e}), median percentile {p2:.0f}")
    hn = sum(x[4] for x in rec); print(f"NEGATIVE CONTROL (surrogate as real, N1, S_low): hits {hn}/{N} (binomial p {binom_p(hn, N):.2f}; expected ~5%)")
    pk = np.array([x[1][2] for x in rec]); print(f"peak winding of real quads: median {np.median(pk):.2f}; share with peak in 3.5-4.5: {np.mean((pk>=3.5)&(pk<=4.5)):.2f}; N2 shifted-clock peaks in 3.5-4.5: {np.mean([(x[3][:,2]>=3.5)&(x[3][:,2]<=4.5) for x in rec]):.2f}")
    print(f"real S_low median {np.median([x[1][0] for x in rec]):.3f} vs N1 {np.median([x[2][:,0].mean() for x in rec]):.3f}, N2 {np.median([x[3][:,0].mean() for x in rec]):.3f}")
    import json; json.dump([dict(quad=x[0], real=list(map(float, x[1])), n1_max=float(x[2][:, 0].max()), n2_max=float(x[3][:, 0].max())) for x in rec], open(H / "quads_frozen.json", "w"), indent=1)
