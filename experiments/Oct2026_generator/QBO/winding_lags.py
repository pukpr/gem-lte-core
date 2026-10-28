#!/usr/bin/env python3
"""Which part of the LTE solution makes the slow descent and which the fast stripes? Per-winding lags between levels
in the nine-leg QBO fit (2026-10-27).
model_L(t) = ... + sum_k amp_{L,k} sin(2 pi k F(t) + phase_{L,k}); the windings k and the manifold F are shared by all
legs, amp and phase are per leg (lt.exe.windings.<leg>.json). F is a folded staircase: plateaus with a slight slope s,
steps at the impulses. A phase difference d(phi) between two levels becomes a time lag d(phi) / (2 pi f), with f the
frequency the winding actually runs at:
  plateau frequency  f_plat = k * s                          (inside a plateau)
  stepped frequency  f_eff  = dominant frequency of sin(2 pi k F(t)) over the record (plateau drift + steps, aliased)
Both predicted lags are listed, and the lag is also MEASURED: each winding's component is built per level on the
5-day F grid (ERA70 leg), and adjacent levels are cross-correlated. The whole model is then split as the data were
(low-pass 200 d = QBO phase; high-pass 120 d = stripes) and its lags are measured the same way.
Levels: the six KIT legs (70, 50, 40, 30, 20, 15 hPa).
Reconstruction (checked exact, max error ~1e-12, against lte_results column 2 on each leg's own grid):
  m(t) = level + k0 F + trend t + accel (t - accel_ref)^2 + sum_j A_j sin(2 pi k_j F + phi_j) + per-leg annual/semiannual,
  model(t) = F9[ m(t) - IR m(t - 1 yr) ]  (IR = 0.79, one 1-2-1 pass).
The per-leg annual/semiannual coefficients are not saved (lt.exe.p holds leg 1's only); they are recovered here by
least squares against column 2, which is exact because everything else is known. On the 5-day grid the 1-yr lag is
73 samples; the F9 pass (symmetric) is left out there, which does not change lags.
    python3 winding_lags.py   -> winding_lags.png, winding_lags.json"""
import json
from pathlib import Path

import matplotlib
import numpy as np
from scipy.signal import butter, filtfilt

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
RUN = HERE.parent / "legs" / "qbo9_kit_era70full_ex2015_2021"
INK, MUTED, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
LV = [70, 50, 40, 30, 20, 15]
a = np.loadtxt(RUN / "lte_results.9.csv", delimiter=","); t, F = a[:, 0], a[:, 3]       # 5-day manifold grid
dt = np.median(np.diff(t))
W = {L: np.array(json.load(open(RUN / f"lt.exe.windings.{i}.json"))["k_amp_phase"]) for i, L in enumerate(LV, 1)}
K = W[70][:, 0]; assert all(np.allclose(W[L][:, 0], K) for L in LV)
s = np.diff(F) / dt; plat = np.abs(s) < 5; s_plat = float(np.median(np.abs(s[plat])))


def ccf_lag(x, y, maxlag):
    """lag (yr) of y behind x at max correlation, parabola-refined."""
    x = (x - x.mean()) / x.std(); y = (y - y.mean()) / y.std(); n = len(x); M = int(round(maxlag / dt)); L = np.arange(-M, M + 1)
    c = np.array([np.mean(x[max(0, -l):n - max(0, l)] * y[max(0, l):n - max(0, -l)]) for l in L]); i = int(np.argmax(c))
    off = 0.0
    if 0 < i < len(c) - 1:
        d = c[i - 1] - 2 * c[i] + c[i + 1]; off = 0.5 * (c[i - 1] - c[i + 1]) / d if d < 0 else 0.0
    return (L[i] + off) * dt, float(c[i])


fr = np.linspace(0.05, 6, 3000)
rows = []
for j, k in enumerate(K):
    comp = np.sin(2 * np.pi * k * F); w = np.hanning(len(t))
    spec = np.array([abs(np.sum(w * (comp - comp.mean()) * np.exp(-2j * np.pi * f * t))) for f in fr]); f_eff = float(fr[np.argmax(spec)])
    f_plat = abs(k) * s_plat
    pairs = []
    for La, Lb in zip(LV[:-1], LV[1:]):                  # La lower level (higher p), Lb the level above
        dphi = float(np.angle(np.exp(1j * (W[Lb][j, 2] - W[La][j, 2]))))
        xa = W[La][j, 1] * np.sin(2 * np.pi * k * F + W[La][j, 2]); xb = W[Lb][j, 1] * np.sin(2 * np.pi * k * F + W[Lb][j, 2])
        lag, r = ccf_lag(xb, xa, 0.5 / f_eff)            # lag of the LOWER level behind the upper one; + = descent
        pairs.append({"pair": f"{Lb}->{La}", "dphi_rad": dphi, "measured_lag_d": lag * 365.25, "r": r})
    amp_share = float(np.mean([abs(W[L][j, 1]) for L in LV]))
    rows.append({"k": float(k), "f_plateau_cpy": f_plat, "f_eff_cpy": f_eff, "period_eff_months": 12 / f_eff, "mean_abs_amp": amp_share,
                 "pairs": pairs, "lag_15_to_70_measured_d": float(sum(p["measured_lag_d"] for p in pairs))})

# whole model, exact reconstruction, evaluated on the 5-day manifold grid
f121 = lambda x: np.r_[x[0], 0.25 * x[:-2] + 0.5 * x[1:-1] + 0.25 * x[2:], x[-1]]


def ir(x, n, IR):
    y = x.copy(); y[n:] = x[n:] - IR * x[:-n]; return y


mod = {}; check = {}
for i, L in enumerate(LV, 1):
    w = json.load(open(RUN / f"lt.exe.windings.{i}.json")); sc = w["secular_context"]
    b = np.loadtxt(RUN / f"lte_results.{i}.csv", delimiter=","); tm, Mm, Fm = b[:, 0], b[:, 1], b[:, 3]
    base = lambda tt, FF: (w["level"] + w["k0"] * FF + w["trend"] * tt + w["accel"] * (tt - sc["accel_ref"]) ** 2
                           + sum(A * np.sin(2 * np.pi * k * FF + ph) for k, A, ph in w["k_amp_phase"]))
    harm = lambda tt: [f(2 * np.pi * h * tt) for h in (1, 2) for f in (np.sin, np.cos)]
    resp = lambda x: f121(ir(x, 12, w["IR"]))
    B = np.array([resp(c) for c in harm(tm)]).T; sl = slice(13, -1)
    c = np.linalg.lstsq(B[sl], (Mm - resp(base(tm, Fm)))[sl], rcond=None)[0]
    check[L] = float(np.max(np.abs(Mm - resp(base(tm, Fm)) - B @ c)[sl]))
    mod[L] = ir(base(t, F) + sum(ci * hi for ci, hi in zip(c, harm(t))), 73, w["IR"])
print("reconstruction max error per level:", {L: f"{v:.1e}" for L, v in check.items()})
fs = 1 / (dt * 365.25)                                                        # samples per day
bl, al = butter(3, 1 / 200, fs=fs); bh, ah = butter(3, 1 / 120, btype="high", fs=fs)
whole = {"QBO phase (low-pass 200 d)": [], "stripes (high-pass 120 d)": []}
mod = {L: v[73:] for L, v in mod.items()}; t, F = t[73:], F[73:]           # drop the first year (no IR lag term)
for La, Lb in zip(LV[:-1], LV[1:]):
    whole["QBO phase (low-pass 200 d)"].append(ccf_lag(filtfilt(bl, al, mod[Lb]), filtfilt(bl, al, mod[La]), 1.2))
    whole["stripes (high-pass 120 d)"].append(ccf_lag(filtfilt(bh, ah, mod[Lb]), filtfilt(bh, ah, mod[La]), 0.15))
# where does the model's high-pass content sit in time? share of high-pass variance within +/- 10 d of a step
steps = a[1:, 0][~plat]; hp = filtfilt(bh, ah, mod[30])
near = np.min(np.abs(t[:, None] - steps[None, :]), axis=1) <= 10 / 365.25
hp_step_share = float(np.sum(hp[near] ** 2) / np.sum(hp ** 2)); near_frac = float(near.mean())

out = {"plateau_slope_per_yr": s_plat, "windings": rows,
       "whole_model": {k: [{"pair": f"{b}->{a_}", "lag_d": l * 365.25, "r": r} for (a_, b), (l, r) in zip(zip(LV[:-1], LV[1:]), v)] for k, v in whole.items()},
       "highpass_variance_within_10d_of_steps": hp_step_share, "time_fraction_within_10d_of_steps": near_frac}
json.dump(out, open(HERE / "winding_lags.json", "w"), indent=1)

print(f"plateau slope |dF/dt| median {s_plat:.3f} per yr; steps per yr {len(steps) / (t[-1] - t[0]):.2f}")
print(f"{'k':>6s} {'f_plat':>7s} {'f_eff':>6s} {'T_eff':>7s} {'|amp|':>6s}  measured lag per pair, lower behind upper (d)            15->70")
for r in rows:
    print(f"{r['k']:+6.2f} {r['f_plateau_cpy']:7.2f} {r['f_eff_cpy']:6.2f} {r['period_eff_months']:6.1f}m {r['mean_abs_amp']:6.2f}  "
          + " ".join(f"{p['measured_lag_d']:+7.1f}" for p in r["pairs"]) + f"   {r['lag_15_to_70_measured_d']:+7.1f}")
for kname, v in out["whole_model"].items():
    print(f"whole model, {kname}: " + "  ".join(f"{p['pair']} {p['lag_d']:+.1f} d (r {p['r']:.2f})" for p in v))
print(f"high-pass (stripe) variance within 10 d of a manifold step: {hp_step_share:.0%} (those windows are {near_frac:.0%} of the time)")

fig, ax = plt.subplots(1, 3, figsize=(17, 5.6), facecolor=SURF, gridspec_kw=dict(width_ratios=[1.15, 1, 1.2], wspace=0.3))
m = (t >= 2000) & (t <= 2004)
ax[0].plot(t[m], F[m], color=INK, lw=1.4)
ax[0].set_title(f"(a) manifold F(t): plateaus slope ~{s_plat:.2f}/yr, steps at the impulses", loc="left", fontsize=10.5, color=INK)
ax[0].set_xlabel("year", color=MUTED); ax[0].set_ylabel("F", color=MUTED)
x = np.array([r["f_eff_cpy"] for r in rows]); y = np.array([r["lag_15_to_70_measured_d"] for r in rows]) / 30.44
sz = 40 + 300 * np.array([r["mean_abs_amp"] for r in rows]) / max(r["mean_abs_amp"] for r in rows)
ax[1].scatter(x, y, s=sz, color="#2a78d6", edgecolor=INK, lw=0.5, zorder=3)
for r, xi, yi in zip(rows, x, y): ax[1].annotate(f"k={r['k']:+.2f}", (xi, yi), xytext=(5, 4), textcoords="offset points", fontsize=8, color=MUTED)
ax[1].axhline(12, color="#d94a4a", ls="--", lw=1); ax[1].text(ax[1].get_xlim()[1] if False else 0.02, 12.4, "observed QBO, 15→70 hPa: ~12 months", fontsize=8, color="#d94a4a")
ax[1].axhline(0, color=MUTED, lw=0.6)
ax[1].set_xlabel("frequency the winding runs at, f_eff (cycles/yr)", color=MUTED); ax[1].set_ylabel("lag 15→70 hPa (months), + = descends", color=MUTED)
ax[1].set_title("(b) per-winding descent (size ∝ amplitude)", loc="left", fontsize=10.5, color=INK)
hpm = filtfilt(bh, ah, mod[30]); m2 = (t >= 2000) & (t <= 2004)
ax[2].plot(t[m2], hpm[m2], color="#d94a4a", lw=1.1, label="model 30 hPa, high-passed 120 d")
for st in steps[(steps >= 2000) & (steps <= 2004)]: ax[2].axvline(st, color=MUTED, lw=0.7, ls=":")
ax[2].set_title(f"(c) model stripes sit on the steps ({hp_step_share:.0%} of variance within ±10 d of a step)", loc="left", fontsize=10.5, color=INK)
ax[2].set_xlabel("year (dotted: manifold steps)", color=MUTED); ax[2].legend(fontsize=8.5, frameon=False)
for a_ in ax: a_.set_facecolor(SURF); a_.tick_params(colors=MUTED)
fig.savefig(HERE / "winding_lags.png", dpi=105, facecolor=SURF, bbox_inches="tight")
print("wrote winding_lags.png")
