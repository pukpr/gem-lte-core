#!/usr/bin/env python3
"""global_sst_mlr.py plus a volcanic aerosol predictor: stratospheric AOD (ext_data/aero_saod.dat) at lag L years.
Negative L (aerosol placed before the eruption) is the placebo. Same target, trend, 2000-2005 exclusion and modeled /
observed index sets as global_sst_mlr.py. Prints the lag scan, the residual around each eruption with and without the
aerosol term, and plots the best positive lag.  usage: global_sst_mlr_aero.py DIR NAME1 NAME2 ..."""
import sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pathlib import Path
E = Path("/home/paul/eval/gem-lte-core/experiments"); sys.path.insert(0, str(E / "Sep2026")); import global_sst_from_quads as G
d, names = sys.argv[1], sys.argv[2:]; cc = lambda x, y: float(np.corrcoef(x, y)[0, 1]); f121 = lambda x: np.r_[x[0], .25 * x[:-2] + .5 * x[1:-1] + .25 * x[2:], x[-1]]
G.MIN_YEAR = 1880.0; keys, sst, lat, lon = G.load_kaplan(E / "Feb2026/k_sst/sst.mean.anom.nc"); gm = np.asarray(G.global_mean(sst, lat), float); gk = {int(k): v for k, v in zip(keys, gm) if np.isfinite(v)}
L = [np.loadtxt(f"{d}/lte_results.{i + 1}.csv", delimiter=",") for i in range(len(names))]; mk = lambda t: int(G.month_key(float(t))); common = sorted(set(gk).intersection(*[set(mk(t) for t in a[:, 0]) for a in L]))
pick = lambda a, col: np.array([dict((mk(t), v) for t, v in a[:, [0, col]])[k] for k in common]); M = np.column_stack([pick(a, 1) for a in L]); O = np.column_stack([pick(a, 2) for a in L]); t = np.array([dict((mk(x), x) for x in L[0][:, 0])[k] for k in common]); y = f121(np.array([gk[k] for k in common]))
g = (t >= 2000) & (t < 2005); tr = ~g; T = np.column_stack([np.ones_like(t), t - t.mean(), (t - t.mean()) ** 2]); dq = lambda v: v - T @ np.linalg.lstsq(T[tr], v[tr], rcond=None)[0]; yd = dq(y)
A = np.loadtxt(E / "ext_data/aero_saod.dat"); aod = lambda lag: np.interp(t - lag, A[:, 0], A[:, 1])
def fit(P):
    X = np.column_stack([T] + [p for p in P if p is not None]); c = np.linalg.lstsq(X[tr], y[tr], rcond=None)[0]; return X @ c, c
base_m, _ = fit([M]); base_o, _ = fit([O])
print(f"{'aerosol lag (yr)':>17s} {'modeled: CC, trend removed':>27s} {'held out':>9s} {'coef (C/AOD)':>13s} {'observed: trend removed':>24s} {'coef':>8s}")
print(f"{'none':>17s} {cc(dq(base_m)[tr], yd[tr]):+27.3f} {cc(dq(base_m)[g], yd[g]):+9.3f} {'':>13s} {cc(dq(base_o)[tr], yd[tr]):+24.3f}")
best = None
for lag in np.arange(-2.0, 2.01, 0.25):
    m, c = fit([M, aod(lag)]); o, co = fit([O, aod(lag)]); r = cc(dq(m)[tr], yd[tr])
    print(f"{lag:+17.2f} {r:+27.3f} {cc(dq(m)[g], yd[g]):+9.3f} {c[-1]:+13.2f} {cc(dq(o)[tr], yd[tr]):+24.3f} {co[-1]:+8.2f}")
    if lag >= 0 and (best is None or r > best[1]): best = (lag, r)
lag = best[0]; m, c = fit([M, aod(lag)]); o, co = fit([O, aod(lag)])
print(f"\nbest physical lag {lag:+.2f} yr; mean residual (data - fit, C), modeled without | with aerosol; observed without | with:")
for nm, a0, a1 in (("Krakatoa 1883", 1883.7, 1886.0), ("Santa Maria 1902", 1903.0, 1905.0), ("Ksudach 1907", 1907.3, 1909.5), ("Katmai 1912", 1912.5, 1914.0),
                   ("Agung 1963", 1963.3, 1966.0), ("~1981 (pre El Chichon)", 1980.0, 1982.2), ("El Chichon 1982", 1982.3, 1985.0), ("Pinatubo 1991", 1991.5, 1994.0)):
    w = (t >= a0) & (t < a1); f = lambda v: f"{np.mean((y - v)[w]):+.3f}"
    print(f"   {nm:24s} {a0:.1f}-{a1:.1f}:  {f(base_m)} | {f(m)}     {f(base_o)} | {f(o)}")
INK, MUTED, SURF, GRID, BL, OR, GR = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df", "#2a78d6", "#eb6834", "#1baf7a"
fig, ax = plt.subplots(2, 1, figsize=(17, 8.5), facecolor=SURF, sharex=True, gridspec_kw=dict(hspace=0.22, height_ratios=[3, 1]))
for A_ in ax:
    A_.set_facecolor(SURF); A_.axvspan(2000, 2005, color=GRID, alpha=0.7, lw=0); A_.tick_params(colors=MUTED); A_.axhline(0, color=GRID, lw=0.8)
    for s_ in ("top", "right"): A_.spines[s_].set_visible(False)
ax[0].plot(t, yd, color=BL, lw=1.1, label="global mean SST, quadratic trend removed")
ax[0].plot(t, dq(base_m), color=MUTED, lw=0.9, alpha=0.8, label=f"modeled indices ({cc(dq(base_m)[tr], yd[tr]):+.3f})")
ax[0].plot(t, dq(m), color=OR, lw=1.2, label=f"modeled indices + aerosol at {lag:+.2f} yr ({cc(dq(m)[tr], yd[tr]):+.3f})")
ax[0].set_ylabel("deg C", color=MUTED); ax[0].legend(frameon=False, fontsize=9.5, ncol=3, loc="upper left")
ax[0].set_title(f"Global mean SST, trend removed: {', '.join(names)} modeled indices with and without a volcanic aerosol term", loc="left", fontsize=11.5, color=INK)
ax[1].fill_between(A[:, 0], A[:, 1], color=MUTED, alpha=0.5, lw=0); ax[1].set_ylabel("stratospheric AOD", color=MUTED); ax[1].set_xlim(t[0], t[-1])
for yr, nm in ((1883.6, "Krakatoa"), (1902.8, "Santa María"), (1907.3, "Ksudach"), (1912.4, "Katmai"), (1963.2, "Agung"), (1982.3, "El Chichón"), (1991.5, "Pinatubo")):
    ax[1].text(yr, ax[1].get_ylim()[1] * (0.72 if nm == "Ksudach" else 0.9), nm, fontsize=8, color=MUTED, ha="center")
out = f"{d}/global_sst_mlr_aero.png"; fig.savefig(out, dpi=100, bbox_inches="tight", facecolor=SURF); print("plot:", out)
