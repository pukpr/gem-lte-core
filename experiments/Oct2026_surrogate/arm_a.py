#!/usr/bin/env python3
"""Arm A (DECLARATION.md): fixed layered manifold (column 4 of each index's current fit), regression
[1, F2, t, cos/sin 2 pi t, cos/sin 4 pi t, sin/cos 2 pi k F2] with k on a 0.25..10 grid (step 0.01)
chosen by training CC; F9 [1/4,1/2,1/4] applied to the data and every column; train = outside 2000-05.
Runs the positive control first (planted series ranked against IAAFT surrogates of itself), then the
real ranking, then the negative control (surrogate s01 ranked against s02..s19)."""
import numpy as np
from pathlib import Path
import surrogates as SG
H = Path(__file__).resolve().parent; S = H / "surr"; N = SG.N; KS = np.round(np.arange(0.25, 10.0001, 0.01), 2)
def f9(y): return np.convolve(y, [.25, .5, .25], "same")
def fit(t, F2, y):
    tr = ~((t >= 2000) & (t < 2005)); te = ~tr; base = [np.ones_like(t), F2, t - 1950, np.cos(2*np.pi*t), np.sin(2*np.pi*t), np.cos(4*np.pi*t), np.sin(4*np.pi*t)]
    Y = f9(y); best = None
    for k in KS:
        X = np.column_stack([f9(c) for c in base + [np.sin(2*np.pi*k*F2), np.cos(2*np.pi*k*F2)]])
        c, *_ = np.linalg.lstsq(X[tr], Y[tr], rcond=None); m = X @ c
        cc = np.corrcoef(m[tr], Y[tr])[0, 1]
        if best is None or cc > best[0]: best = (cc, np.corrcoef(m[te], Y[te])[0, 1], k)
    return best
def rank_report(lab, x0, xs):
    xs = np.array(xs); r = 1 + int((xs >= x0).sum())
    return f"{lab}: {x0:.4f} rank {r}/{len(xs)+1}, surrogates {xs.mean():.4f}+-{xs.std():.4f} (max {xs.max():.4f}), margin over best {(x0-xs.max())/xs.std():+.2f} sd"
if __name__ == "__main__":
    for i in SG.IDX:
        a = np.loadtxt(H / "seeds" / i / "lte_results.csv", delimiter=",")[:N]; t, F2 = a[:, 0], a[:, 3]
        # positive control
        pl = np.loadtxt(S / f"{i}_planted.dat")[:, 1]
        plS = SG.iaaft_multi(pl[None, :], np.random.default_rng(99), iters=100)[0:1]
        pls = [fit(t, F2, SG.iaaft_multi(pl[None, :], np.random.default_rng(100 + d), iters=100)[0])[0] for d in range(19)]
        p0 = fit(t, F2, pl)[0]
        print(f"[{i}] POSITIVE CONTROL " + rank_report("planted train CC", p0, pls))
        R = [fit(t, F2, np.loadtxt(S / f"{i}_s{d:02d}.dat")[:, 1]) for d in range(20)]
        print(f"[{i}] REAL            " + rank_report("train CC", R[0][0], [x[0] for x in R[1:]]) + f"  | best k real {R[0][2]}, test CC real {R[0][1]:+.3f} vs surrogates {np.mean([x[1] for x in R[1:]]):+.3f}")
        print(f"[{i}] REAL (test)     " + rank_report("test CC", R[0][1], [x[1] for x in R[1:]]))
        print(f"[{i}] NEGATIVE CONTROL " + rank_report("s01 train CC", R[1][0], [x[0] for x in R[2:]]))
