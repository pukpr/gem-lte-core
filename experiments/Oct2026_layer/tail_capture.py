#!/usr/bin/env python3
"""Tail-capture asymmetry of the global-SST quad composite (2026-10-03).
For linearly detrended monthly global SST D and composite model M in a period:
  capture(tail) = mean(M) / mean(D) over the months where D is in that tail
  (top / bottom 5, 10, 20 % of D); the overall amplitude ratio is std(M)/std(D).
Also the tail regression slopes (D on M restricted to tail months is biased, so we use
capture) and a symmetric check conditioned on the MODEL's tails, plus a block bootstrap
(resampling whole years) for the warm-minus-cold capture difference.
Sets: LAYER (Oct2026_layer/_backcast), original (Feb2026_sweep/_backcast)."""
import sys, numpy as np
from pathlib import Path
sys.path.insert(0, "/home/paul/eval/gem-lte-core/experiments/Sep2026"); import global_sst_from_quads as g
E = Path("/home/paul/eval/gem-lte-core/experiments")
def composite(root):
    names = sorted(p.name for p in root.iterdir() if p.is_dir() and g.NAME_RE.match(p.name))
    W, _ = g.ocean_area_weights(names, g.load_kaplan(g.DEFAULT_NC)); ws = sum(W.values())
    bc = {n: np.loadtxt(root / "_backcast" / n / "lte_results.csv", delimiter=",") for n in names}
    return bc[names[0]][:, 0], sum(W[n] / ws * bc[n][:, 1] for n in names)
g.MIN_YEAR = 1880.0; keys, sst, lat, lon = g.load_kaplan(g.DEFAULT_NC); km = dict(zip(keys.tolist(), g.global_mean(sst, lat)))
rng = np.random.default_rng(1)
def capture(D, M, q, cond="data"):
    X = D if cond == "data" else M
    hi, lo = X >= np.quantile(X, 1 - q), X <= np.quantile(X, q)
    if cond == "data": return M[hi].mean() / D[hi].mean(), M[lo].mean() / D[lo].mean()
    return D[hi].mean() / M[hi].mean(), D[lo].mean() / M[lo].mean()
for lab, root in (("LAYER=1 (Oct2026_layer)", E / "Oct2026_layer"), ("original (Feb2026_sweep)", E / "Feb2026_sweep")):
    t, M0 = composite(root); D0 = np.array([km.get(g.month_key(x), np.nan) for x in t])
    print(f"\n=== {lab}")
    for lo_, hi_ in ((1950, 2024), (1880, 1950)):
        m = (t >= lo_) & (t < hi_) & np.isfinite(D0); tt = t[m]; x = tt - tt.mean()
        D = D0[m] - np.polyval(np.polyfit(x, D0[m], 1), x); M = M0[m] - np.polyval(np.polyfit(x, M0[m], 1), x)
        sk = lambda y: ((y - y.mean())**3).mean() / y.std()**3
        print(f"  {lo_}-{hi_-1}: amplitude ratio std(M)/std(D) {M.std()/D.std():.2f}; skew data {sk(D):+.2f} model {sk(M):+.2f}; CC {np.corrcoef(D, M)[0,1]:+.2f}")
        for q in (0.05, 0.10, 0.20):
            w, c = capture(D, M, q); wm, cm = capture(D, M, q, "model")
            yrs = np.floor(tt).astype(int); U = np.unique(yrs); diffs = []
            for _ in range(2000):
                pick = np.concatenate([np.where(yrs == y)[0] for y in rng.choice(U, len(U))])
                a, b = capture(D[pick], M[pick], q); diffs.append(a - b)
            lo95, hi95 = np.percentile(diffs, [2.5, 97.5])
            print(f"    tails {int(q*100):2d}%: data-conditioned capture  warm {w:.2f}  cold {c:.2f}  warm-cold {w-c:+.2f} [95% CI {lo95:+.2f}, {hi95:+.2f}]"
                  f" | mean excursion data warm {D[D>=np.quantile(D,1-q)].mean():+.3f} cold {D[D<=np.quantile(D,q)].mean():+.3f}")
        # largest positive data events: model fraction at each
        if lo_ == 1950:
            d3 = np.convolve(D, np.ones(3) / 3, "same"); m3 = np.convolve(M, np.ones(3) / 3, "same"); ev = []
            for y in range(1951, 2023):
                w_ = (tt >= y) & (tt < y + 1)
                i = np.where(w_)[0][np.argmax(d3[w_])]; j = np.where(w_)[0][np.argmin(d3[w_])]; ev.append((y, d3[i], m3[i], d3[j], m3[j]))
            ev = np.array(ev); top = ev[np.argsort(-ev[:, 1])[:8]]; bot = ev[np.argsort(ev[:, 3])[:8]]
            print("    8 warmest years (3-mo peak, data/model): " + ", ".join(f"{int(r[0])} {r[1]:+.2f}/{r[2]:+.2f}" for r in top) + f"  -> mean ratio {np.mean(top[:,2]/top[:,1]):.2f}")
            print("    8 coldest years (3-mo trough, data/model): " + ", ".join(f"{int(r[0])} {r[3]:+.2f}/{r[4]:+.2f}" for r in bot) + f"  -> mean ratio {np.mean(bot[:,4]/bot[:,3]):.2f}")
