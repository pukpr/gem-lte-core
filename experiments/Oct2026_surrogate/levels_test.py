#!/usr/bin/env python3
"""DECLARATION_LEVELS.md: level-structure (LTE Modulation) surrogate test, low windings only."""
import numpy as np
from pathlib import Path
import surrogates as SG
H = Path(__file__).resolve().parent; S = H / "surr"; N = SG.N; B = 12
CROSS = {"amo": "pdo", "pdo": "amo", "nino4": "amo"}; KS = np.round(np.arange(0.25, 10.0001, 0.01), 2)
def det(y): x = np.arange(len(y)); return y - np.polyval(np.polyfit(x, y, 1), x)
def eta2(y, F):
    y = det(y); edges = np.quantile(F, np.linspace(0, 1, B + 1)); b = np.clip(np.searchsorted(edges, F, side="right") - 1, 0, B - 1)
    mu = np.array([y[b == j].mean() for j in range(B)]); e = ((mu[b] - y.mean())**2).sum() / ((y - y.mean())**2).sum()
    n = len(y); return 1 - (1 - e) * (n - 1) / (n - B)
def dbic(y, F, t):
    y = det(y); n = len(y); X0 = np.column_stack([np.ones(n), F, t - t.mean()])
    r0 = y - X0 @ np.linalg.lstsq(X0, y, rcond=None)[0]; best = (np.inf, None)
    for k in KS:
        X = np.column_stack([X0, np.sin(2*np.pi*k*F), np.cos(2*np.pi*k*F)]); r = y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
        d = n * np.log((r @ r) / (r0 @ r0)) + 3 * np.log(n)
        if d < best[0]: best = (d, k)
    return best
def rank(x0, xs): xs = np.array(xs); return 1 + int((xs >= x0).sum()), (x0 - xs.max()) / xs.std(), xs.mean(), xs.std()
def F2_of(i): a = np.loadtxt(H / "seeds" / i / "lte_results.csv", delimiter=",")[:N]; return a[:, 0], a[:, 3]
if __name__ == "__main__":
    for i in SG.IDX:
        t, Fo = F2_of(i); _, Fc = F2_of(CROSS[i]); Y = [np.loadtxt(S / f"{i}_s{d:02d}.dat")[:, 1] for d in range(20)]
        pl = np.loadtxt(S / f"{i}_planted.dat")[:, 1]
        pls = [SG.iaaft_multi(pl[None, :], np.random.default_rng(300 + d), iters=100)[0] for d in range(19)]
        r = rank(eta2(pl, Fc), [eta2(s, Fc) for s in pls]); print(f"[{i}] POSITIVE CONTROL  S1 cross: planted rank {r[0]}/20, margin {r[1]:+.2f} sd")
        r = rank(eta2(Y[1], Fc), [eta2(s, Fc) for s in Y[2:]]); print(f"[{i}] NEGATIVE CONTROL  S1 cross: s01 rank {r[0]}/19, margin {r[1]:+.2f} sd")
        for lab, F in (("S1 cross manifold (" + CROSS[i] + ")", Fc), ("S1 own manifold", Fo)):
            v = [eta2(y, F) for y in Y]; r = rank(v[0], v[1:])
            print(f"[{i}] {lab:26s}: real eta2_adj {v[0]:.3f}, rank {r[0]}/20, surrogates {r[2]:.3f}+-{r[3]:.3f}, margin over best {r[1]:+.2f} sd")
        for lab, F in (("cross", Fc), ("own", Fo)):
            bb = [dbic(y, F, t) for y in Y]; r = rank(-bb[0][0], [-b[0] for b in bb[1:]])
            print(f"[{i}] S2 {lab:5s} penalized winding: real dBIC {bb[0][0]:+.1f} at k={bb[0][1]}, rank {r[0]}/20; surrogates dBIC {np.mean([b[0] for b in bb[1:]]):+.1f}, k chosen {sorted(b[1] for b in bb[1:])}")
