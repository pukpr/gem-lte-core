#!/usr/bin/env python3
"""Which windings make the 5-day model's fast cycling, and does that cycling match the data's? (2026-10-15, per user)
The saved 5-day model (lte_results.2.csv) is decomposed exactly: it is regressed on the winding basis sin/cos(2 pi k F)
of the low terms (slow, base, harmonics up to 30) and of the high terms (above 30), plus level, F and the annual terms,
each column passed through the model's own one-year lag and 1-2-1 filter. The two fitted parts are then compared with
the data band by band, train and held out 2000-2005. Also the monthly leg, and a zoom figure."""
import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
cc = lambda a, b: float(np.corrcoef(a, b)[0, 1]); f121 = lambda x: np.r_[x[0], .25 * x[:-2] + .5 * x[1:-1] + .25 * x[2:], x[-1]]; p = json.load(open("lt.exe.p")); kb = abs(p["ltep"][-1]); ks = p["ltep"][0]; ir = p["IR"]
H = sorted(int(h) for h in p["harm"]); LO = [h for h in H if h <= 30]; HI = [h for h in H if h > 30]; sm = lambda v, k: np.convolve(v, np.ones(k) / k, mode="same")
def mask(x, h):
    c = np.r_[0, np.cumsum(x)]; i = np.arange(len(x)); A = np.maximum(0, i - h); B = np.minimum(len(x) - 1, i + h); return x > (c[B + 1] - c[A]) / (B - A + 1)
def split(i, N):
    a = np.loadtxt(f"lte_results.{i}.csv", delimiter=","); t, m, y, F = a[:, 0], a[:, 1], a[:, 2], a[:, 3]; sc = lambda k: [np.sin(2 * np.pi * k * F), np.cos(2 * np.pi * k * F)]
    base = [np.ones_like(t), F, t - t.mean(), np.cos(2 * np.pi * t), np.sin(2 * np.pi * t), np.cos(4 * np.pi * t), np.sin(4 * np.pi * t)]; low = sc(ks) + sc(kb) + [c for h in LO for c in sc(h * kb)]; high = [c for h in HI for c in sc(h * kb)]
    mk = mask(m, int(10 * N / 2)); tr = lambda x: f121(np.r_[x[:N], x[N:] - ir * x[:-N] * mk[:-N]]); X = np.column_stack([tr(c) for c in base + low + high]); c = np.linalg.lstsq(X, m, rcond=None)[0]
    nb, nl = len(base), len(low); P_low = X[:, :nb + nl] @ c[:nb + nl]; P_high = X[:, nb + nl:] @ c[nb + nl:]; return t, m, y, P_low, P_high, cc(P_low + P_high, m)
print(f"windings: slow {ks:.4f}, base {kb:.4f}; low harmonics {LO}; high harmonics {HI}\n")
for i, nm, N in ((2, "SOI 5-day", 73), (1, "CRU SOI monthly", 12)):
    t, m, y, L, Hh, rec = split(i, N); g = (t >= 2000) & (t < 2005); print(f"{nm}: decomposition reproduces the saved model at CC {rec:.4f}; spread of model {m.std():.2f}, low part {L.std():.2f}, high part {Hh.std():.2f}, data {y.std():.2f}")
    bands = [("whole series", lambda v: v)] + ([("faster than 6 months", lambda v: v - sm(v, 37)), ("faster than 2 months", lambda v: v - sm(v, 13)), ("2 to 6 months", lambda v: sm(v, 13) - sm(v, 37))] if N == 73 else [("faster than 1 year", lambda v: v - sm(v, 13)), ("faster than 6 months", lambda v: v - sm(v, 7))])
    for bn, f in bands:
        fy, fm, fl, fh = f(y), f(m), f(L), f(Hh); print(f"   {bn:22s} model train {cc(fm[~g], fy[~g]):+.3f} held out {cc(fm[g], fy[g]):+.3f} | low part {cc(fl[~g], fy[~g]):+.3f} / {cc(fl[g], fy[g]):+.3f} | high part {cc(fh[~g], fy[~g]):+.3f} / {cc(fh[g], fy[g]):+.3f} | high part's share of the model's variance in this band {100 * fh.var() / fm.var():.0f}%")
    print()
t, m, y, L, Hh, _ = split(2, 73); z = lambda v: (v - v.mean()) / y.std()
INK, MUTED, SURF, GRID, BL, OR, GR = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df", "#2a78d6", "#eb6834", "#1baf7a"; fig, ax = plt.subplots(2, 1, figsize=(17, 8.5), facecolor=SURF, gridspec_kw=dict(hspace=0.3))
for A_, (a0, a1) in zip(ax, ((2000, 2006), (2009, 2016))):
    k = (t >= a0) & (t < a1); A_.set_facecolor(SURF); A_.axvspan(2000, 2005, color=GRID, alpha=0.6, lw=0) if a0 < 2005 else None
    A_.plot(t[k], (y[k] - y.mean()) / y.std(), color=BL, lw=1.0, label="SOI 5-day"); A_.plot(t[k], (L[k] - L.mean() + Hh[k] - Hh.mean()) / y.std() * y.std() / m.std(), color=OR, lw=1.2, label="model (low + high)"); A_.plot(t[k], (L[k] - L.mean()) / m.std(), color=GR, lw=1.2, label="model, low harmonics only")
    A_.set_xlim(a0, a1); A_.legend(frameon=False, fontsize=9, ncol=3, loc="upper left"); A_.set_title(f"{a0}-{a1}" + ("  (2000-2005 held out)" if a0 < 2005 else ""), loc="left", fontsize=11, color=INK)
    for s in ("top", "right"): A_.spines[s].set_visible(False)
fig.savefig("fast_split.png", dpi=100, bbox_inches="tight", facecolor=SURF); print("plot: fast_split.png")
