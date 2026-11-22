#!/usr/bin/env python3
"""How much of a leg's model-data mismatch is +-1..3 month timing of the excursions?

    python3 dtw_alignment.py RUN_DIR [--legs 1 2] [--windows 1 2 3] [--held 1960 1970] [--null 40]

For each leg's lte_results.<k>.csv (col 2 model, col 3 data):
  CC          plain correlation
  DTW score   lt.exe's METRIC=DTW, (D(-data, data) - D(model, data)) / D(-data, data), |x - y| cost,
              Sakoe-Chiba band of +-w months (gem-lte-primitives.adb DTW_Distance)
  warped CC   correlation of model and data after the optimal +-w month warp (z-scored, squared cost)
  null        the same warped CC for the model cut into 5-yr blocks and shuffled: warping inflates
              any correlation, so only the excess over the null credits the model's shape
  offsets     share of months warped by 0, 1, 2+ months, and mean |offset| on the 20% largest
              |data| excursions vs the rest
by whole record (data present), held-out window and 36-yr era.
"""
import argparse
from pathlib import Path

import numpy as np


def dtw(x, y, w, cost):
    """Banded DTW; returns total cost and the warping offsets j - i along the path, per i."""
    n = len(x)
    INF = np.inf
    D = np.full((n + 1, 2 * w + 1), INF)   # D[i, j - i + w], i, j 1-based
    D[0, w] = 0.0
    for i in range(1, n + 1):
        lo, hi = max(1, i - w), min(n, i + w)
        for j in range(lo, hi + 1):
            c = cost(x[i - 1], y[j - 1])
            k = j - i + w
            best = INF
            if k - 1 >= 0:
                best = min(best, D[i, k - 1])          # (i, j-1)
            best = min(best, D[i - 1, k])              # (i-1, j-1)
            if k + 1 <= 2 * w:
                best = min(best, D[i - 1, k + 1])      # (i-1, j)
            D[i, k] = c + best
    # backtrack
    i, k = n, w
    pairs = []
    while i > 0:
        j = i + k - w
        pairs.append((i - 1, j - 1))
        cands = []
        if k - 1 >= 0 and j - 1 >= 1:
            cands.append((D[i, k - 1], i, k - 1))
        if i - 1 >= 1 or (i - 1 == 0 and k == w):
            cands.append((D[i - 1, k], i - 1, k))
        if k + 1 <= 2 * w and i - 1 >= 1:
            cands.append((D[i - 1, k + 1], i - 1, k + 1))
        if not cands:
            break
        _, i, k = min(cands)
    return D[n, w], pairs[::-1]


def dtw_score(m, d, w):
    a = lambda u, v: abs(u - v)
    return (dtw(-d, d, w, a)[0] - dtw(m, d, w, a)[0]) / dtw(-d, d, w, a)[0]


def warped_cc(m, d, w):
    zm, zd = (m - m.mean()) / m.std(), (d - d.mean()) / d.std()
    _, pairs = dtw(zm, zd, w, lambda u, v: (u - v) ** 2)
    i, j = np.array(pairs).T
    return np.corrcoef(zm[i], zd[j])[0, 1], i, j


def block_shuffle(m, block, rng):
    nb = int(np.ceil(len(m) / block))
    order = rng.permutation(nb)
    return np.concatenate([m[b * block:(b + 1) * block] for b in order])[:len(m)]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("--legs", nargs="+", type=int, default=[1, 2])
    ap.add_argument("--windows", nargs="+", type=int, default=[1, 2, 3])
    ap.add_argument("--held", nargs=2, type=float, default=(1960, 1970))
    ap.add_argument("--null", type=int, default=40)
    ap.add_argument("--results", default=None, help="directory holding lte_results.<k>.csv (default RUN_DIR)")
    a = ap.parse_args()
    run = Path(a.run)
    res = Path(a.results) if a.results else run
    rng = np.random.default_rng(3)
    print(f"{run.name}{'' if res == run else ' (' + str(res) + ')'}")
    for leg in a.legs:
        r = np.loadtxt(res / f"lte_results.{leg}.csv", delimiter=",")
        ok = r[:, 2] != 0
        i0, i1 = np.where(ok)[0][[0, -1]]
        t, m, d = r[i0:i1 + 1, 0], r[i0:i1 + 1, 1], r[i0:i1 + 1, 2]
        y0 = t[0]
        segs = [("all", np.ones(len(t), bool)), (f"held {a.held[0]:g}-{a.held[1]:g}", (t >= a.held[0]) & (t <= a.held[1]))]
        segs += [(f"{y0 + 36 * e:.0f}-{y0 + 36 * (e + 1):.0f}", (t >= y0 + 36 * e) & (t < y0 + 36 * (e + 1))) for e in range(4)]
        print(f"\n leg {leg}: {t[0]:.1f}-{t[-1]:.1f}")
        print(f"  {'segment':14s} {'CC':>6s} " + " ".join(f"{'DTW' + str(w):>6s}" for w in a.windows) +
              " " + " ".join(f"{'wCC' + str(w):>6s} {'gain':>6s} {'ampN':>6s} {'z':>5s}" for w in a.windows))
        for name, s in segs:
            if s.sum() < 24:
                continue
            ms, ds = m[s], d[s]
            cc = np.corrcoef(ms, ds)[0, 1]
            sc = [dtw_score(ms, ds, w) for w in a.windows]
            row = []
            for w in a.windows:
                wc = warped_cc(ms, ds, w)[0]
                nn = a.null if name == "all" else max(10, a.null // 4)
                nulls = []
                for _ in range(nn):
                    sh = block_shuffle(ms, 60, rng)
                    nulls.append(warped_cc(sh, ds, w)[0] - np.corrcoef(sh, ds)[0, 1])
                # amplitude-error null: same model, residual replaced by a phase-randomized copy
                # (same spectrum and plain CC, no timing error); the real data's warped gain above
                # this null is what +-w month timing misalignment explains
                zm, zd = (ms - ms.mean()) / ms.std(), (ds - ds.mean()) / ds.std()
                beta = np.dot(zm, zd) / np.dot(zm, zm)
                res_ = zd - beta * zm
                F_ = np.fft.rfft(res_)
                amp = []
                for _ in range(nn):
                    ph = rng.uniform(0, 2 * np.pi, len(F_)); ph[0] = 0
                    syn = beta * zm + np.fft.irfft(np.abs(F_) * np.exp(1j * (np.angle(F_) + ph)), len(res_))
                    amp.append(warped_cc(ms, syn, w)[0] - np.corrcoef(ms, syn)[0, 1])
                gain = wc - cc
                z_ = (gain - np.mean(amp)) / np.std(amp)
                row.append(f"{wc:+.3f} {gain:+.3f} {np.mean(amp):+.3f} {z_:+5.1f}")
            print(f"  {name:14s} {cc:+.3f} " + " ".join(f"{v:+.3f}" for v in sc) + "  " + "  ".join(row))
        # where does the warp act?
        w = max(a.windows)
        _, i, j = warped_cc(m, d, w)
        off = np.zeros(len(t))
        for ii, jj in zip(i, j):
            off[ii] = jj - ii
        big = np.abs(d - d.mean()) >= np.percentile(np.abs(d - d.mean()), 80)
        share = [np.mean(np.abs(off) == 0), np.mean(np.abs(off) == 1), np.mean(np.abs(off) >= 2)]
        print(f"  warp +-{w}: months offset 0 / 1 / 2+: {share[0]:.2f} / {share[1]:.2f} / {share[2]:.2f}; "
              f"mean |offset| on largest-20% |data| excursions {np.abs(off[big]).mean():.2f} vs rest {np.abs(off[~big]).mean():.2f}; "
              f"mean offset (data lags model) {off.mean():+.2f}")
    print("\n  gain = wCC - CC for the real data; ampN = the same gain when the misfit is pure amplitude error "
          "(phase-randomized residual, same spectrum and CC); z = (gain - ampN) / sd: > 2 means real +-w month "
          "timing errors in the misfit")


if __name__ == "__main__":
    main()
