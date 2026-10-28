#!/usr/bin/env python3
"""Can fast windings riding on the manifold plateaus carry the observed stripes and their tilt? (2026-10-27)
LTE picture: a fast winding sin(2 pi k F(t) + phi_L) on a plateau of slope s runs at k*s cycles/yr, so a per-level
phase phi_L shifts it by only (d phi)/(2 pi k s) -- days -- while slow windings, stepping between plateaus, turn the
same phase change into months. The nine-leg QBO fit has |k| <= 2.65 (k*s <= 1.1 cpy), so its stripes come from the
manifold steps instead and are vertical (winding_lags.py). Here the fast windings are added directly:
  target  daily ERA5 u (0N 0E 00 UTC) at 20, 30, 50, 70 hPa, high-passed (minus 120-day running mean)
  basis   sin, cos (2 pi k F(t)), k = 15..300 (k*s spans about 4-60 day periods), F = the nine-leg manifold
          (5-day ERA70 leg grid, linear inside plateaus, step at the midpoint of its 5-day gap)
  fit     ridge regression per level (per-level amplitudes and phases, like the legs), trained outside 2015-2021,
          lambda chosen on 2009-2015 inside the training years; scored on 2015-2021
  lag     adjacent-level lag of the fitted stripe series (whole record), compared with the observed 1.5-2.4 d
Controls (same k grid, same ridge): RAMP = straight F with the median plateau slope (a plain Fourier basis, no
manifold); SHIFT = the real F delayed by one year (same staircase statistics, wrong mapping to time).
    python3 fast_windings.py   -> fast_windings.png, fast_windings.json"""
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
RUN = HERE.parent / "legs" / "qbo9_kit_era70full_ex2015_2021"
INK, MUTED, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
LV = [20, 30, 50, 70]
a = np.loadtxt(RUN / "lte_results.9.csv", delimiter=","); t5, F5 = a[:, 0], a[:, 3]
td = np.loadtxt(HERE / "era_qbo_u30_point0E_daily.dat")[:, 0]
td = td[(td >= t5[0]) & (td <= t5[-1])]


def daily_F(t5, F5, td):
    i = np.clip(np.searchsorted(t5, td) - 1, 0, len(t5) - 2); x0, x1 = t5[i], t5[i + 1]; y0, y1 = F5[i], F5[i + 1]
    step = np.abs(y1 - y0) / (x1 - x0) > 5
    lin = y0 + (y1 - y0) * (td - x0) / (x1 - x0)
    return np.where(step, np.where(td < 0.5 * (x0 + x1), y0, y1), lin)


F = daily_F(t5, F5, td)
s_med = float(np.median(np.abs(np.diff(F5) / np.diff(t5))[np.abs(np.diff(F5) / np.diff(t5)) < 5]))
MANIFOLDS = {"manifold F": F, "RAMP (no manifold)": s_med * (td - td[0]), "SHIFT (F delayed 1 yr)": np.interp(td - 1.0, td, F, left=F[0])}
K = np.arange(15, 301, 1.0)


def hp120(x):
    w = np.ones(120) / 120; return x - np.convolve(x, w, "same") / np.convolve(np.ones_like(x), w, "same")


Y = {L: hp120(np.interp(td, *np.loadtxt(HERE / f"era_qbo_u{L}_point0E_daily.dat").T)) for L in LV}
test = (td >= 2015) & (td <= 2021); val = (td >= 2009) & (td < 2015); tr_inner = ~test & ~val; train = ~test


def basis(FF):
    ph = 2 * np.pi * np.outer(FF, K); return np.hstack([np.sin(ph), np.cos(ph)])


def ridge(X, y, lam):
    A = X.T @ X; return np.linalg.solve(A + lam * np.trace(A) / A.shape[0] * np.eye(A.shape[0]), X.T @ y)


def ccf_lag(x, y, M=10):
    x = (x - x.mean()) / x.std(); y = (y - y.mean()) / y.std(); n = len(x); Ls = np.arange(-M, M + 1)
    c = np.array([np.mean(x[max(0, -l):n - max(0, l)] * y[max(0, l):n - max(0, -l)]) for l in Ls]); i = int(np.argmax(c))
    off = 0.0
    if 0 < i < len(c) - 1:
        d = c[i - 1] - 2 * c[i] + c[i + 1]; off = 0.5 * (c[i - 1] - c[i + 1]) / d if d < 0 else 0.0
    return float(Ls[i] + off), float(c[i])


out = {}; fits = {}
for name, FF in MANIFOLDS.items():
    X = basis(FF); X = X - X[train].mean(0); res = {"levels": {}}; fits[name] = {}
    for L in LV:
        y = Y[L] - Y[L][train].mean(); best = None
        for lam in (1e-3, 1e-2, 1e-1, 1, 10):
            b = ridge(X[tr_inner], y[tr_inner], lam); r = float(np.corrcoef(X[val] @ b, y[val])[0, 1])
            if best is None or r > best[1]: best = (lam, r)
        b = ridge(X[train], y[train], best[0]); fit = X @ b; fits[name][L] = fit
        res["levels"][L] = {"lambda": best[0], "val_r": best[1], "train_r": float(np.corrcoef(fit[train], y[train])[0, 1]),
                            "test_r_2015_2021": float(np.corrcoef(fit[test], y[test])[0, 1])}
    res["fitted_stripe_lags_d"] = {f"{u}->{l}": ccf_lag(fits[name][u], fits[name][l]) for u, l in zip(LV[:-1], LV[1:])}
    out[name] = res
out["observed_stripe_lags_d"] = {f"{u}->{l}": ccf_lag(Y[u], Y[l]) for u, l in zip(LV[:-1], LV[1:])}
json.dump({"k_range": [15, 300], "plateau_slope": s_med, **out}, open(HERE / "fast_windings.json", "w"), indent=1)

print(f"plateau slope {s_med:.3f}/yr; k {K[0]:.0f}..{K[-1]:.0f} -> periods {365.25 / (K[-1] * 0.42):.1f}..{365.25 / (K[0] * 0.27):.0f} d")
for name in MANIFOLDS:
    r = out[name]
    print(f"== {name}")
    for L, v in r["levels"].items():
        print(f"  {L} hPa  lambda {v['lambda']:g}  train r {v['train_r']:+.2f}  val r {v['val_r']:+.2f}  TEST 2015-21 r {v['test_r_2015_2021']:+.2f}")
    print("  fitted stripe lags (lower behind upper): " + "  ".join(f"{p} {l:+.1f} d (r {c:.2f})" for p, (l, c) in r["fitted_stripe_lags_d"].items()))
print("observed stripe lags: " + "  ".join(f"{p} {l:+.1f} d (r {c:.2f})" for p, (l, c) in out["observed_stripe_lags_d"].items()))

fig, ax = plt.subplots(2, 1, figsize=(15, 8), facecolor=SURF)
m = (td >= 2016.0) & (td <= 2017.0)
ax[0].plot(td[m], Y[30][m], color=INK, lw=1, label="observed 30 hPa, high-passed")
for (name, col) in (("manifold F", "#d94a4a"), ("RAMP (no manifold)", "#2a78d6"), ("SHIFT (F delayed 1 yr)", "#1d9a6c")):
    ax[0].plot(td[m], fits[name][30][m], color=col, lw=1.1, label=f"{name}: test r {out[name]['levels'][30]['test_r_2015_2021']:+.2f}")
ax[0].set_title("(a) fast-winding fits, held-out year 2016, 30 hPa", loc="left", fontsize=10.5, color=INK); ax[0].legend(fontsize=8.5, frameon=False, ncol=2)
x = np.arange(len(LV)); wdt = 0.27
for i, (name, col) in enumerate((("manifold F", "#d94a4a"), ("RAMP (no manifold)", "#2a78d6"), ("SHIFT (F delayed 1 yr)", "#1d9a6c"))):
    ax[1].bar(x + (i - 1) * wdt, [out[name]["levels"][L]["test_r_2015_2021"] for L in LV], wdt, color=col, label=name)
ax[1].axhline(0, color=MUTED, lw=0.6); ax[1].set_xticks(x); ax[1].set_xticklabels([f"{L} hPa" for L in LV])
ax[1].set_title("(b) held-out 2015–2021 correlation with the observed stripes", loc="left", fontsize=10.5, color=INK); ax[1].legend(fontsize=8.5, frameon=False)
for a_ in ax: a_.set_facecolor(SURF); a_.tick_params(colors=MUTED)
fig.tight_layout(); fig.savefig(HERE / "fast_windings.png", dpi=105, facecolor=SURF)
print("wrote fast_windings.png")
