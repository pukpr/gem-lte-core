"""Held-out 2015-2021 pattern check per level: r, model/obs size, peak matching (+-60 d), best lag, DTW (+-90 d band)
against a randomly time-shifted model null. Compares several runs of the same five ERA legs (2026-10-27)."""
import sys, json, numpy as np
from scipy.signal import find_peaks
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
# usage: heldout_pattern.py OUT tag=dir [tag=dir ...]   (first tag is drawn with the observations)
OUT = sys.argv[1]; RUNS = [a.split("=", 1) for a in sys.argv[2:]]; LV = [10, 20, 30, 50, 70]
def load(d, i):
    a = np.loadtxt(f"{d}/lte_results.{i}.csv", delimiter=","); m = a[:, 2] != 0; return a[m, 0], a[m, 1], a[m, 2]
def sm(x, n=3): return np.convolve(x, np.ones(n) / n, "same")
def peaks(x, ref):  # positive peaks, prominence 0.5 sd of obs
    p, _ = find_peaks(sm(x), prominence=0.5 * ref.std(), distance=12); return p
def dtw(x, y, band):
    n = len(x); D = np.full((n + 1, n + 1), np.inf); D[0, 0] = 0
    for i in range(1, n + 1):
        for j in range(max(1, i - band), min(n, i + band) + 1):
            D[i, j] = abs(x[i-1] - y[j-1]) + min(D[i-1, j], D[i, j-1], D[i-1, j-1])
    return D[n, n] / n
res = {}; fig, ax = plt.subplots(5, 1, figsize=(13, 13))
for k, L in enumerate(LV, 1):
    r = {}
    for n_, (tag, d) in enumerate(RUNS):
        t, m, o = load(d, k); h = (t >= 2015) & (t <= 2021); tr = ~h
        oh, mh = o[h], m[h]; po, pm = peaks(oh, oh), peaks(mh, oh)
        match = np.mean([np.min(np.abs(t[h][pm] - t[h][q])) <= 60 / 365.25 for q in po]) if len(pm) else 0
        false = np.mean([np.min(np.abs(t[h][po] - t[h][q])) > 60 / 365.25 for q in pm]) if len(pm) else 0
        lags = np.arange(-36, 37); cc = [np.corrcoef(np.roll(mh, l)[40:-40], oh[40:-40])[0, 1] for l in lags]
        z = lambda v: (v - v.mean()) / v.std(); rnd = np.random.default_rng(0)
        dd = dtw(z(oh), z(mh), 18)   # +-90 d warping
        null = [dtw(z(oh), z(np.roll(mh, rnd.integers(60, len(mh) - 60))), 18) for _ in range(30)]
        r[tag] = {"r": float(np.corrcoef(mh, oh)[0, 1]), "sd_ratio_held": float(mh.std() / oh.std()), "sd_ratio_train": float(m[tr].std() / o[tr].std()),
                  "max_ratio_held": float(np.abs(mh).max() / np.abs(oh).max()), "obs_peaks": len(po), "model_peaks": len(pm),
                  "obs_peaks_matched_60d": float(match), "model_peaks_unmatched": float(false),
                  "best_lag_d": int(lags[int(np.argmax(cc))] * 5), "best_lag_r": float(max(cc)), "dtw90": float(dd), "dtw90_null": [float(np.mean(null)), float(np.std(null))]}
        if n_ == 0: ax[k-1].plot(t[h], oh, "k", lw=1.2, label="ERA obs"); ax[k-1].plot(t[h][po], oh[po], "kv")
        ax[k-1].plot(t[h], mh, lw=1, label=f"{tag} r {r[tag]['r']:+.2f}", color=["#2a78d6", "#d94a4a", "#1d9a6c", "#a35bd6"][n_])
    ax[k-1].set_title(f"{L} hPa, held out 2015-2021", loc="left", fontsize=10); ax[k-1].legend(fontsize=8, ncol=3, frameon=False)
    res[L] = r
fig.tight_layout(); fig.savefig(f"{OUT}.png", dpi=95); json.dump(res, open(f"{OUT}.json", "w"), indent=1)
for L, r in res.items():
    for tag, v in r.items():
        print(f"{L:3d} {tag:9s} r {v['r']:+.2f} | sd ratio held {v['sd_ratio_held']:.2f} train {v['sd_ratio_train']:.2f} | max ratio {v['max_ratio_held']:.2f} | peaks obs {v['obs_peaks']} model {v['model_peaks']} matched {v['obs_peaks_matched_60d']:.2f} spurious {v['model_peaks_unmatched']:.2f} | best lag {v['best_lag_d']:+d} d r {v['best_lag_r']:+.2f} | DTW90 {v['dtw90']:.3f} null {v['dtw90_null'][0]:.3f}±{v['dtw90_null'][1]:.3f}")
