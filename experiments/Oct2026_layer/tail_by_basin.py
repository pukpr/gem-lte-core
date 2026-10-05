#!/usr/bin/env python3
"""Warm/cold tail capture of the LAYER=1 quad composites by basin and for the uniform
positive-gain runs (gain 1.0 = _backcast, 1.5/2.0 = _alpha_g*). Linearly detrended monthly,
capture = model mean / data mean over the data's top / bottom decile; CI from resampling years."""
import sys, numpy as np
from pathlib import Path
sys.path.insert(0, "/home/paul/eval/gem-lte-core/experiments/Sep2026"); import global_sst_from_quads as g
L = Path(__file__).resolve().parent
names = sorted(p.name for p in L.iterdir() if p.is_dir() and g.NAME_RE.match(p.name))
W, _ = g.ocean_area_weights(names, g.load_kaplan(g.DEFAULT_NC))
def basin(n):
    lat, lon = g.quad_centre(n)
    if abs(lat) >= 60: return "polar"
    if lon >= 130 or lon <= -75: return "Pacific"
    if -75 < lon <= 15: return "Atlantic"
    return "Indian"
RUNS = {"gain 1.0": L / "_backcast", "gain 1.5": L / "_alpha_g1.5", "gain 2.0": L / "_alpha_g2.0"}
A = {k: {n: np.loadtxt(v / n / "lte_results.csv", delimiter=",") for n in names} for k, v in RUNS.items()}
t = A["gain 1.0"][names[0]][:, 0]; rng = np.random.default_rng(2)
g.MIN_YEAR = 1880.0; keys, sst, lat, lon = g.load_kaplan(g.DEFAULT_NC); km = dict(zip(keys.tolist(), g.global_mean(sst, lat)))
Dk = np.array([km.get(g.month_key(x), np.nan) for x in t])
def det(y, m): x = t[m] - t[m].mean(); return y[m] - np.polyval(np.polyfit(x, y[m], 1), x)
def stats(D, M, tt):
    hi, lo = D >= np.quantile(D, .9), D <= np.quantile(D, .1)
    w, c = M[hi].mean() / D[hi].mean(), M[lo].mean() / D[lo].mean(); yrs = np.floor(tt).astype(int); U = np.unique(yrs); diffs = []
    for _ in range(1000):
        p = np.concatenate([np.where(yrs == y)[0] for y in rng.choice(U, len(U))]); Dp, Mp = D[p], M[p]
        h, l = Dp >= np.quantile(Dp, .9), Dp <= np.quantile(Dp, .1); diffs.append(Mp[h].mean() / Dp[h].mean() - Mp[l].mean() / Dp[l].mean())
    a, b = np.percentile(diffs, [2.5, 97.5])
    return f"CC {np.corrcoef(D, M)[0,1]:+.2f} amp {M.std()/D.std():.2f} | warm {w:.2f} cold {c:.2f} diff {w-c:+.2f} [{a:+.2f},{b:+.2f}] | top1% {np.percentile(M,99)/np.percentile(D,99):.2f} bot1% {np.percentile(M,1)/np.percentile(D,1):.2f} | near-max months {int((M > M.max()-0.1*np.ptp(M)).sum())} near-min {int((M < M.min()+0.1*np.ptp(M)).sum())}"
for per, (lo_, hi_) in (("1950-2023 fit", (1950, 2024)), ("1880-1949 back-cast", (1880, 1950))):
    print(f"\n##### {per}")
    for grp in ("GLOBAL (Kaplan)", "Pacific", "Atlantic", "Indian", "polar"):
        sel = names if grp.startswith("GLOBAL") else [n for n in names if basin(n) == grp]; ws = sum(W[n] for n in sel)
        for k in RUNS:
            Dq = sum(W[n] / ws * A[k][n][:, 2] for n in sel) if not grp.startswith("GLOBAL") else Dk
            M = sum(W[n] / ws * A[k][n][:, 1] for n in sel)
            m = (t >= lo_) & (t < hi_) & np.isfinite(Dq)
            if per.startswith("1880") and not grp.startswith("GLOBAL"):
                m &= np.abs(Dq) > 0   # quad data before 1950 = the 1880-2023 Kaplan box series
            print(f"  {grp:16s} {k}: {stats(det(Dq, m), det(M, m), t[m])}")
