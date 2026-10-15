#!/usr/bin/env python3
"""The SOI as the model's own Tahiti-minus-Darwin differential (two-leg Darwin/Tahiti fit, deseasonalized, full records),
compared with the official SOI series (2026-10-15):
   CRU SOI   https://crudata.uea.ac.uk/cru/data/soi/soi.dat (Ropelewski-Jones, same two stations, from 1866)
   project SOI  Feb2026/soi/soi.dat (the monthly SOI the SOI fits use, from 1880)
Model differential = model Tahiti - model Darwin (each leg's model, in hPa); also the data differential from the two
station files, as a check of how well the station anomalies reproduce each official SOI. The official series are given
the same 1-2-1 smoothing as the fitted data (F9 = 1); scores without it are shown too. Held out: 2000-2005, which
neither leg trained on."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
cc = lambda a, b: float(np.corrcoef(a, b)[0, 1]); f121 = lambda x: np.r_[x[0], .25 * x[:-2] + .5 * x[1:-1] + .25 * x[2:], x[-1]]
D = np.loadtxt("lte_results.1.csv", delimiter=","); T = np.loadtxt("lte_results.2.csv", delimiter=","); kd = {round(t * 12): (m, y) for t, m, y in D[:, :3]}
rows = [(t, mt - kd[round(t * 12)][0], yt - kd[round(t * 12)][1]) for t, mt, yt in T[:, :3] if round(t * 12) in kd]; A = np.array(rows); t, md, dd = A[:, 0], A[:, 1], A[:, 2]
cru = {}
for r in (l.split() for l in open("cru_soi.dat") if l.split()):
    for m in range(12):
        v = float(r[1 + m])
        if v > -99: cru[int(r[0]) * 12 + m] = v
P = np.loadtxt("../../../Feb2026/soi/soi.dat"); proj = {int(round(x * 12)): v for x, v in P}
def series(src):
    k = np.array([int(round(x * 12)) in src for x in t]); v = np.array([src[int(round(x * 12))] for x in t[k]]); return k, v
g = (t >= 2000) & (t < 2005); print(f"model differential (Tahiti - Darwin), {len(t)} months {t[0]:.2f}-{t[-1]:.2f}")
fig, ax = plt.subplots(2, 1, figsize=(17, 8.6), facecolor="#fcfcfb", gridspec_kw=dict(hspace=0.3)); z = lambda v: (v - v.mean()) / v.std()
for i, (nm, src) in enumerate((("CRU SOI", cru), ("project monthly SOI", proj))):
    k, v = series(src); vs = f121(v); tt, m_, d_, gg = t[k], md[k], dd[k], g[k]; era = lambda f, a0, a1: cc(m_[(tt >= a0) & (tt < a1) & ~gg], f[(tt >= a0) & (tt < a1) & ~gg])
    print(f"\n{nm}: {int(k.sum())} common months {tt[0]:.2f}-{tt[-1]:.2f}")
    print(f"   model differential vs {nm} (smoothed like the data): train {cc(m_[~gg], vs[~gg]):+.3f}, held out 2000-05 {cc(m_[gg], vs[gg]):+.3f} | unsmoothed: train {cc(m_[~gg], v[~gg]):+.3f}, held out {cc(m_[gg], v[gg]):+.3f}")
    print(f"   by era (smoothed, outside 2000-05): " + "  ".join(f"{a0}-{a1}: {era(vs, a0, a1):+.2f}" for a0, a1 in ((1860, 1880), (1880, 1920), (1920, 1960), (1960, 2000), (2005, 2027)) if ((tt >= a0) & (tt < a1)).sum() > 60))
    print(f"   check, station-data differential vs {nm}: CC {cc(d_, vs):+.3f} (the ceiling the model differential could reach)")
    A_ = ax[i]; A_.set_facecolor("#fcfcfb"); A_.axvspan(2000, 2005, color="#e4e3df", alpha=0.7, lw=0); A_.plot(tt, z(vs), color="#2a78d6", lw=1.0, label=f"{nm} (1-2-1 smoothed, standardized)"); A_.plot(tt, z(m_), color="#eb6834", lw=1.0, label="model: Tahiti minus Darwin (standardized)")
    A_.set_title(f"{nm} against the model differential: train {cc(m_[~gg], vs[~gg]):+.3f}, held out 2000-05 {cc(m_[gg], vs[gg]):+.3f}", loc="left", fontsize=11); A_.legend(frameon=False, fontsize=9, ncol=2, loc="upper left"); A_.set_xlim(1950, tt[-1])
    for s_ in ("top", "right"): A_.spines[s_].set_visible(False)
fig.savefig("soi_differential.png", dpi=100, bbox_inches="tight", facecolor="#fcfcfb"); print("\nplot: soi_differential.png (from 1950)")
