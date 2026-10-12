#!/usr/bin/env python3
"""How tightly does dLOD constrain each tidal line, and how does that compare with what the AMO+PDO two-leg search did
to it? (2026-10-12)
1. Least squares of the measured dLOD (Feb2026/dlod3.dat, 1962-2019) on all 42 table lines at the model year: each
   line's amplitude with its standard error (residual-based, inflated for the residual's autocorrelation).
2. For each line: seed and final amplitude in the search, the change in units of dLOD's standard error, and the loss
   in the overall dLOD correlation if that line ALONE is changed by the same factor (from the fitted table).
usage: dlod_police.py seed.p final.p"""
import json, sys, numpy as np
E = "/home/paul/eval/gem-lte-core/experiments/"; a = json.load(open(sys.argv[1])); b = json.load(open(sys.argv[2])); D = np.loadtxt(E + "Feb2026/dlod3.dat"); t, y = D[:, 0], D[:, 1] - D[:, 1].mean()
YL = 365.2422484 + 0.00405305 + a["year"]; P = np.array([x[0] for x in a["lpap"]]); w = 2 * np.pi * YL / P
X = np.column_stack([np.ones_like(t), t - t.mean()] + [f(wi * t) for wi in w for f in (np.sin, np.cos)]); beta, *_ = np.linalg.lstsq(X, y, rcond=None); fit = X @ beta; r = y - fit; cc = lambda u, v: float(np.corrcoef(u, v)[0, 1])
rho = cc(r[:-1], r[1:]); neff = len(t) * (1 - rho) / (1 + rho); s2 = r.var() * len(t) / neff; cov = s2 * np.linalg.inv(X.T @ X) * 1.0; cc0 = cc(fit, y)
print(f"dLOD: {len(t)} samples {t[0]:.1f}-{t[-1]:.1f}; 42-line fit CC {cc0:.5f}; residual lag-1 autocorrelation {rho:.2f} (effective samples {neff:.0f})\n")
print(f"{'period (d)':>10s} {'dLOD amplitude':>14s} {'+- se':>8s} {'seed':>8s} {'final':>8s} {'change':>7s} {'change / se':>11s} {'dLOD CC loss if this line alone':>32s}")
rows = []
for j in range(len(P)):
    s, c = beta[2 + 2 * j], beta[3 + 2 * j]; A = np.hypot(s, c); se = np.sqrt(0.5 * (cov[2 + 2 * j, 2 + 2 * j] + cov[3 + 2 * j, 3 + 2 * j])); A0, A1 = a["lpap"][j][1], b["lpap"][j][1]; f = A1 / A0
    b2 = beta.copy(); b2[2 + 2 * j] *= f; b2[3 + 2 * j] *= f; rows.append((abs(f - 1), P[j], A, se, A0, A1, f, abs(A1 - A0) / se, cc0 - cc(X @ b2, y)))
for _, p, A, se, A0, A1, f, zc, loss in sorted(rows, reverse=True)[:14]: print(f"{p:10.4g} {A:14.5f} {se:8.5f} {A0:8.5f} {A1:8.5f} {100 * (f - 1):+6.0f}% {zc:11.1f} {loss:32.6f}")
tot = beta.copy()
for j in range(len(P)): f = b["lpap"][j][1] / a["lpap"][j][1]; tot[2 + 2 * j] *= f; tot[3 + 2 * j] *= f
print(f"\nall amplitude changes together: dLOD CC {cc0:.5f} -> {cc(X @ tot, y):.5f}")
big = [r_ for r_ in rows if r_[7] > 3]; print(f"lines moved by more than 3 standard errors of dLOD's own estimate: {len(big)} of {len(rows)}; by more than 10: {sum(r_[7] > 10 for r_ in rows)}")
strong = [r_ for r_ in rows if r_[2] > 0.02]; print(f"the {len(strong)} strong lines (dLOD amplitude > 0.02) changed by at most {100 * max(r_[0] for r_ in strong):.1f}%")
