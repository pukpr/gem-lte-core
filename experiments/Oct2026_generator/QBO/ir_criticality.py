#!/usr/bin/env python3
"""Is the IR term a real annual feedback or just a filter that trades off against the windings? (2026-10-27)
model(t) = F9[ m(t) - IR m(t - 1 yr) ], m = level + k0 F + trend t + accel (t - ref)^2 + sum_k A_k sin(2 pi k F + phi_k)
+ annual/semiannual. Every term of m except those in F is closed under the 1-yr lag (a constant, t, t^2 and the
annual harmonics map back into their own span), so IR can only matter through F(t - 1 yr):
  1. degeneracy  per winding k, R^2 of sin(2 pi k F(t - 1 yr)) on {sin, cos}(2 pi k F(t)). R^2 ~ 1 means the lagged
                 copy is just a re-phased, re-scaled version of the same winding, so IR is absorbed by the amplitudes
                 and phases (a filter, conflated with the windings); R^2 << 1 means IR adds structure of its own.
  2. profile     with the manifold F and the k set fixed, the linear coefficients are refit by least squares for each
                 IR on a grid (per ERA level, 2015-2021 held out) -> train r, held-out r, held-out DTW (z-scored,
                 +-90 d warping band). A sharp shared optimum = criticality; a flat curve = IR is a nuisance filter.
  3. free lag    IR = 0 with the lagged windings sin/cos(2 pi k F(t - 1 yr)) added as free regressors: the
                 unconstrained version of what a single IR imposes.
Manifolds: the five-level ERA fits searched with CC (IR 0.944) and DTW (IR 0.891), and the nine-leg fit (IR 0.789), each with its own k set,
on the shared ERA 5-day grid (one year = 73 samples; F9 = one 1-2-1 pass, as lt.exe).
    python3 ir_criticality.py   -> ir_criticality.png, ir_criticality.json"""
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
LEGS = HERE.parent / "legs"
INK, MUTED, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
LV = [10, 20, 30, 50, 70]
LAG = 73
IRS = np.array([0, 0.2, 0.4, 0.5, 0.6, 0.7, 0.79, 0.85, 0.9, 0.944, 0.97, 0.99])
MANIFOLDS = {"5-level CC fit (IR 0.944)": ("qbo5_era_ex2015_2021", 5), "5-level DTW fit (IR 0.891)": ("qbo5_era_dtw_ex2015_2021", 5),
             "9-leg fit (IR 0.789)": ("qbo9_kit_era70full_ex2015_2021", 9)}

f121 = lambda x: np.r_[x[0], 0.25 * x[:-2] + 0.5 * x[1:-1] + 0.25 * x[2:], x[-1]]


def lag1(x): return np.r_[np.full(LAG, np.nan), x[:-LAG]]


def dtw(x, y, band=18):
    n = len(x); D = np.full((n + 1, n + 1), np.inf); D[0, 0] = 0
    for i in range(1, n + 1):
        for j in range(max(1, i - band), min(n, i + band) + 1):
            D[i, j] = abs(x[i - 1] - y[j - 1]) + min(D[i - 1, j], D[i, j - 1], D[i - 1, j - 1])
    return D[n, n] / n


def z(v): return (v - v.mean()) / v.std()


data = {}
for i, L in enumerate(LV, 1):
    b = np.loadtxt(LEGS / "qbo5_era_ex2015_2021" / f"lte_results.{i}.csv", delimiter=","); data[L] = (b[:, 0], b[:, 2])
out = {}
fig, ax = plt.subplots(2, 4, figsize=(21, 8.5), facecolor=SURF)
for col, (name, (run, legF)) in enumerate(MANIFOLDS.items()):
    g = np.loadtxt(LEGS / run / f"lte_results.{legF}.csv", delimiter=","); tg, Fg = g[:, 0], g[:, 3]
    w1 = json.load(open(LEGS / run / "lt.exe.windings.1.json")); K = np.array([r[0] for r in w1["k_amp_phase"]]); ref = w1["secular_context"]["accel_ref"]
    # 1. degeneracy, on the full grid
    dF = Fg[LAG:] - Fg[:-LAG]; deg = []
    for k in K:
        y = np.sin(2 * np.pi * k * Fg[:-LAG]); X = np.c_[np.sin(2 * np.pi * k * Fg[LAG:]), np.cos(2 * np.pi * k * Fg[LAG:])]
        c = np.linalg.lstsq(X, y, rcond=None)[0]; deg.append(float(1 - np.var(y - X @ c) / np.var(y)))
    yF = Fg[:-LAG]; XF = np.c_[np.ones(len(dF)), Fg[LAG:], tg[LAG:]]
    degF = float(1 - np.var(yF - XF @ np.linalg.lstsq(XF, yF, rcond=None)[0]) / np.var(yF))
    res = {"k": K.tolist(), "IR_fit": w1["IR"], "dF_1yr_mean": float(dF.mean()), "dF_1yr_sd": float(dF.std()),
           "R2_lagged_winding_on_own_span": deg, "R2_lagged_F_on_F_t_1": degF, "levels": {}}
    for L in LV:
        t, d = data[L]; n = min(len(t), len(tg)); t, d, F = t[:n], d[:n], Fg[:n]
        m_terms = [np.ones(n), F, t, (t - ref) ** 2] + [f(2 * np.pi * k * F) for k in K for f in (np.sin, np.cos)] \
            + [f(2 * np.pi * h * t) for h in (1, 2) for f in (np.sin, np.cos)]
        ok = np.arange(n) >= LAG + 1; held = (t >= 2015) & (t <= 2021); tr = ok & ~held; ho = ok & held
        r_tr, r_ho, d_ho = [], [], []
        for IR in IRS:
            X = np.array([f121(np.nan_to_num(x - IR * lag1(x))) for x in m_terms]).T
            c = np.linalg.lstsq(X[tr], d[tr], rcond=None)[0]; mod = X @ c
            r_tr.append(float(np.corrcoef(mod[tr], d[tr])[0, 1])); r_ho.append(float(np.corrcoef(mod[ho], d[ho])[0, 1]))
            d_ho.append(dtw(z(d[ho]), z(mod[ho])))
        lagged = [f(2 * np.pi * k * lag1(F)) for k in K for f in (np.sin, np.cos)] + [lag1(F)]
        X = np.array([f121(np.nan_to_num(x)) for x in m_terms + lagged]).T
        c = np.linalg.lstsq(X[tr], d[tr], rcond=None)[0]; mod = X @ c
        free = {"train_r": float(np.corrcoef(mod[tr], d[tr])[0, 1]), "held_r": float(np.corrcoef(mod[ho], d[ho])[0, 1]), "held_dtw": dtw(z(d[ho]), z(mod[ho]))}
        res["levels"][L] = {"train_r": r_tr, "held_r": r_ho, "held_dtw": d_ho, "free_lag": free,
                            "best_IR_train": float(IRS[int(np.argmax(r_tr))]), "best_IR_held": float(IRS[int(np.argmax(r_ho))])}
        cl = plt.cm.viridis(LV.index(L) / 4)
        ax[0, col].plot(IRS, r_tr, color=cl, ls=":", lw=1); ax[0, col].plot(IRS, r_ho, color=cl, lw=1.6, marker="o", ms=3, label=f"{L} hPa")
        ax[1, col].plot(IRS, d_ho, color=cl, lw=1.6, marker="o", ms=3, label=f"{L} hPa")
        ax[1, col].axhline(free["held_dtw"], color=cl, lw=0.6, ls="--")
    ax[0, col].axvline(w1["IR"], color=MUTED, lw=0.8); ax[1, col].axvline(w1["IR"], color=MUTED, lw=0.8)
    ax[0, col].set_title(f"{name}\nheld-out r (solid), train r (dotted) vs IR", loc="left", fontsize=10, color=INK)
    ax[1, col].set_title("held-out DTW (lower = better); dashed = free lagged windings", loc="left", fontsize=10, color=INK)
    ax[0, col].legend(fontsize=8, frameon=False); ax[1, col].set_xlabel("IR", color=MUTED)
    out[name] = res
    ax[0, 3].bar(np.arange(len(K)) + (col - 1) * 0.28, deg, 0.28, label=name, color=["#d94a4a", "#1d9a6c", "#2a78d6"][col])
ax[0, 3].set_title("R^2 of the 1-yr-lagged winding on its own span\n(1 = IR fully absorbed by amplitude/phase)", loc="left", fontsize=10, color=INK)
ax[0, 3].set_xlabel("winding index (k sorted as in the fit)", color=MUTED); ax[0, 3].legend(fontsize=8, frameon=False)
fr = np.linspace(0, 2, 400)
for IR in (0.5, 0.79, 0.944):
    ax[1, 3].plot(fr, np.abs(1 - IR * np.exp(-2j * np.pi * fr)), lw=1.4, label=f"IR {IR}")
ax[1, 3].axvspan(0.38, 0.48, color="#e4e3df"); ax[1, 3].set_xlabel("frequency (cpy); grey = QBO band", color=MUTED)
ax[1, 3].set_title("gain of 1 - IR z^(-1 yr)", loc="left", fontsize=10, color=INK); ax[1, 3].legend(fontsize=8, frameon=False)
for a_ in ax.flat: a_.set_facecolor(SURF); a_.tick_params(colors=MUTED)
fig.tight_layout(); fig.savefig(HERE / "ir_criticality.png", dpi=100, facecolor=SURF)
json.dump(out, open(HERE / "ir_criticality.json", "w"), indent=1)
for name, res in out.items():
    print(f"== {name}: 1-yr manifold advance dF {res['dF_1yr_mean']:+.2f} +- {res['dF_1yr_sd']:.2f}; R^2 lagged F on (1, F, t) {res['R2_lagged_F_on_F_t_1']:.3f}")
    print("   R^2 lagged winding on own span: " + " ".join(f"k{k:+.2f}:{r:.2f}" for k, r in zip(res["k"], res["R2_lagged_winding_on_own_span"])))
    print("   IR grid: " + " ".join(f"{x:5.2f}" for x in IRS))
    for L, v in res["levels"].items():
        print(f"   {L:3d} held r " + " ".join(f"{x:+5.2f}" for x in v["held_r"]) + f" | best train {v['best_IR_train']:.2f} held {v['best_IR_held']:.2f}"
              f" | free lag: train {v['free_lag']['train_r']:+.2f} held {v['free_lag']['held_r']:+.2f} DTW {v['free_lag']['held_dtw']:.2f}")
        print(f"       train r " + " ".join(f"{x:+5.2f}" for x in v["train_r"]))
        print(f"       held DTW" + " ".join(f"{x:6.2f}" for x in v["held_dtw"]))
print("gain |1 - IR e^(-i 2 pi f)| at f = 0, 0.43 (QBO), 1 cpy: " + "; ".join(
    f"IR {IR}: " + ", ".join(f"{abs(1 - IR * np.exp(-2j * np.pi * f)):.2f}" for f in (0, 0.43, 1)) for IR in (0.5, 0.79, 0.944)))
