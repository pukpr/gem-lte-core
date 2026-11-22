#!/usr/bin/env python3
"""Which tidal lines and manifold parameters separate two fits' dLOD agreement?

    python3 compare_tidal_factors.py RUN_A RUN_B [--top 15]

dLOD agreement (lt.exe's dLOD: line) is CC(R2, M2) over dlod_compare.csv, where R2 is the
observed dLOD regressed on the fit's own tidal periods and M2 the fit's tidal sum
(gem-lte-primitives-solution.adb). Both are sums of the same lines, so each run's
columns are regressed jointly on sin/cos of the lpap periods to recover per-line
amplitude and phase; M2 is scaled to R2 (CC is scale-free) and the misfit
M2 - s R2 is split by line. 1 - CC^2 is the misfit power share, so a line's share of
it says how much that line costs the run's dLOD agreement.
Also prints the manifold scalars and windings side by side, and compares the two
runs' forcing manifolds F(t) (lte_results.1.csv column 4).
"""
import argparse
import json
from pathlib import Path

import numpy as np


def line_fit(t_yr, y, periods_d, year_days):
    w = 2 * np.pi * year_days / np.asarray(periods_d)       # rad / yr
    X = np.column_stack([f(np.outer(t_yr - 1990.0, w)) for f in (np.sin, np.cos)])
    coef, *_ = np.linalg.lstsq(X, y - y.mean(), rcond=None)
    s, c = coef[:len(w)], coef[len(w):]
    return s + 1j * c, X @ coef


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--top", type=int, default=15)
    a = ap.parse_args()
    runs = [Path(a.a), Path(a.b)]
    P = [json.load(open(r / "lt.exe.p")) for r in runs]
    per = [np.array([l[0] for l in p["lpap"]]) for p in P]
    assert len(per[0]) == len(per[1]) and np.allclose(per[0], per[1], rtol=1e-6), "different line sets"
    periods = per[0]
    year_days = [365.2422 * (1 + p.get("year", 0.0)) for p in P]

    lines, misfit = [], []
    print("dLOD agreement decomposition (per-line complex amplitudes from dlod_compare.csv):")
    for r, yd in zip(runs, year_days):
        c = np.loadtxt(r / "dlod_compare.csv", delimiter=",")
        t, R2, M2 = c[:, 0], c[:, 1], c[:, 2]
        gR, fR = line_fit(t, R2, periods, yd)
        gM, fM = line_fit(t, M2, periods, yd)
        cc = np.corrcoef(R2, M2)[0, 1]
        s = np.dot(M2 - M2.mean(), R2 - R2.mean()) / np.dot(R2 - R2.mean(), R2 - R2.mean())
        d = gM - s * gR                                       # misfit per line
        pw = np.abs(d) ** 2 / np.sum(np.abs(d) ** 2)
        print(f"  {r.name}: CC {cc:.5f} (1-CC^2 = {1 - cc * cc:.5f}); line fit explains R2 {1 - np.var(R2 - R2.mean() - fR) / np.var(R2):.6f}, "
              f"M2 {1 - np.var(M2 - M2.mean() - fM) / np.var(M2):.6f}; scale s {s:.3f}")
        lines.append((gR, gM, s))
        misfit.append((d, pw))

    # rank lines by the difference in misfit power (absolute, in units of R2 variance)
    varR = [np.sum(np.abs(l[0] * l[2]) ** 2) for l in lines]
    absmis = [np.abs(m[0]) ** 2 / v for m, v in zip(misfit, varR)]
    order = np.argsort(-(absmis[1] - absmis[0]))
    names = [r.name[:22] for r in runs]
    print(f"\nlines ranked by extra misfit in {names[1]} over {names[0]} "
          "(misfit power / tidal-sum power; ratio and phase are each run's line relative to the dLOD reference):")
    print(f"  {'period d':>10s} {'ref amp':>8s} | {'A ratio':>7s} {'A dphi':>6s} {'A mis':>7s} | {'B ratio':>7s} {'B dphi':>6s} {'B mis':>7s} | {'B-A':>7s}")
    for i in order[:a.top]:
        row = []
        for (gR, gM, s), am in zip(lines, absmis):
            ratio = abs(gM[i]) / abs(s * gR[i]) if abs(gR[i]) > 0 else np.nan
            dphi = np.degrees(np.angle(gM[i] / gR[i])) if abs(gR[i]) > 0 else np.nan
            row.append(f"{ratio:7.3f} {dphi:+6.0f} {am[i]:7.4f}")
        print(f"  {periods[i]:10.4f} {abs(lines[0][0][i]):8.4f} | " + " | ".join(row) + f" | {absmis[1][i] - absmis[0][i]:+7.4f}")
    print(f"  total misfit share: {names[0]} {absmis[0].sum():.4f}, {names[1]} {absmis[1].sum():.4f}")

    print("\nmanifold scalars:")
    keys = [k for k in P[0] if not isinstance(P[0][k], list)]
    for k in keys:
        va, vb = P[0][k], P[1].get(k)
        flag = "" if vb is None or np.isclose(va, vb, rtol=0.05, atol=1e-12) else "  <-- differs"
        print(f"  {k:6s} {va:+14.6g} {vb if vb is None else f'{vb:+14.6g}'}{flag}")
    for k in ("ltep", "harm"):
        print(f"  {k:6s} {P[0].get(k)}  |  {P[1].get(k)}")

    F = [np.loadtxt(r / "lte_results.1.csv", delimiter=",")[:, [0, 3]] for r in runs]
    ka = {int(round(t * 12)): f for t, f in F[0]}
    kb = {int(round(t * 12)): f for t, f in F[1]}
    ks = sorted(set(ka) & set(kb))
    fa, fb = np.array([ka[k] for k in ks]), np.array([kb[k] for k in ks])
    t = np.array(ks) / 12
    da, db = np.diff(fa), np.diff(fb)
    print(f"\nforcing manifold F(t), monthly {t[0]:.0f}-{t[-1]:.0f}: range A {np.ptp(fa):.2f}, B {np.ptp(fb):.2f}; "
          f"CC(F_A, F_B) {np.corrcoef(fa, fb)[0, 1]:+.3f}; CC of monthly increments dF {np.corrcoef(da, db)[0, 1]:+.3f}")
    fr = np.fft.rfftfreq(len(da), 1 / 12)
    for name, x in (("A", da), ("B", db)):
        X = np.abs(np.fft.rfft((x - x.mean()) * np.hanning(len(x)))) ** 2
        top = np.argsort(-X[1:])[:6] + 1
        print(f"  dF spectrum peaks {name}: " + ", ".join(f"{1 / fr[i]:.2f} yr ({X[i] / X[1:].sum():.2f})" for i in sorted(top, key=lambda i: -X[i])))
    Xa = np.fft.rfft((da - da.mean()) * np.hanning(len(da)))
    Xb = np.fft.rfft((db - db.mean()) * np.hanning(len(db)))
    for lo, hi in ((0, 1 / 10), (1 / 10, 1 / 2), (1 / 2, 1.0), (1.0, 6.0)):
        s_ = (fr >= lo) & (fr < hi)
        coh = abs(np.sum(Xa[s_] * np.conj(Xb[s_]))) / np.sqrt(np.sum(abs(Xa[s_]) ** 2) * np.sum(abs(Xb[s_]) ** 2))
        print(f"  band {('%.1f' % (1 / hi)) if hi else 'inf'}-{('%.1f' % (1 / lo)) if lo else 'inf'} yr: power A {np.sum(abs(Xa[s_])**2)/np.sum(abs(Xa)**2):.2f}, "
              f"B {np.sum(abs(Xb[s_])**2)/np.sum(abs(Xb)**2):.2f}, coherence {coh:.2f}")


if __name__ == "__main__":
    main()
