#!/usr/bin/env python3
"""Phase-free pattern test for fast windings: do the stripe spectra stretch with the plateau slope? (2026-10-27)
A fast winding sin(2 pi k F(t) + phi) cannot keep its phase across a manifold step (2 pi k dF is huge for k ~ 100),
so held-out phase skill is the wrong test. What survives is the PATTERN: inside plateau p, F ~ F0 + s_p t, so a fixed
set of windings k shows up at frequencies f = k |s_p|. With |s_p| varying 0.13-0.42 /yr from plateau to plateau, the
stripe spectra should line up across plateaus when plotted against k = f / |s_p| (better than against f, and better
than with the slopes assigned to the wrong plateaus).
Per plateau (the nine-leg manifold's, 5-day ERA70 grid; |s| > 0.1, >= 120 d long) and level (daily ERA5 0N 0E,
20/30/50/70 hPa, plus 10 hPa if present; u minus 120-day running mean): Hann periodogram over 3-90 cpy, log power
whitened by a straight line in log f, then mapped to a common k grid. Consistency C = mean pairwise correlation of
the plateau spectra. C_true (own slopes) is compared with C under permutation of slopes among plateaus (3000 draws),
stratified by season (plateau start in Apr-Sep vs Oct-Mar) x QBO state (sign of the plateau-mean u at that level), so
that neither the seasonal cycle of wave activity nor Doppler shifting by the QBO wind can produce the alignment.
Reference: C in plain frequency (no scaling).
Power check: PLANT=<fraction> PLANT_LEVEL=<hPa> adds sum_k sin(2 pi k F(t)) for k = 60, 110, 170 (equal amplitudes,
F daily: linear inside plateaus, steps at their 5-day gap midpoints), scaled to that fraction of the stripe variance,
to that level's stripes. The phase decorrelates at every step, as for a real fast winding. PLANT_JITTER=<sd> makes
the planted signal's plateau slopes differ from the model's by a random fraction (sd), i.e. an imperfect manifold.
    python3 plateau_scaling.py   -> plateau_scaling.png, plateau_scaling.json"""
import json
import os
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
RUN = HERE.parent / "legs" / "qbo9_kit_era70full_ex2015_2021"
INK, MUTED, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
rng = np.random.default_rng(7)
PLANT = float(os.environ.get("PLANT", 0)); JIT = float(os.environ.get("PLANT_JITTER", 0)); PLANT_LEVEL = int(os.environ.get("PLANT_LEVEL", 30))
LV = [L for L in (10, 20, 30, 50, 70) if (HERE / f"era_qbo_u{L}_point0E_daily.dat").exists()]
LV = [PLANT_LEVEL] if PLANT else LV
a = np.loadtxt(RUN / "lte_results.9.csv", delimiter=","); t5, F5 = a[:, 0], a[:, 3]
d = np.diff(F5) / np.diff(t5); st = np.where(np.abs(d) > 5)[0]; edges = np.r_[0, st + 1, len(t5)]
PL = []
for i0, i1 in zip(edges[:-1], edges[1:]):
    if i1 - i0 < 20: continue
    s = np.polyfit(t5[i0:i1], F5[i0:i1], 1)[0]
    if abs(s) > 0.1 and (t5[i1 - 1] - t5[i0]) * 365.25 >= 120: PL.append((t5[i0] + 5 / 365.25, t5[i1 - 1] - 5 / 365.25, abs(s)))
PL = np.array(PL); SEAS = (((PL[:, 0] % 1) * 12 >= 3) & ((PL[:, 0] % 1) * 12 < 9)).astype(int)
FR = np.linspace(3, 90, 600)                                   # cpy
KG = np.linspace(3 / PL[:, 2].min(), 90 / PL[:, 2].max(), 400)  # common k range covered by every plateau
assert KG[0] < KG[-1], "no common k range"


def hp120(x):
    w = np.ones(120) / 120; return x - np.convolve(x, w, "same") / np.convolve(np.ones_like(x), w, "same")


def C(S):
    Z = (S - S.mean(1, keepdims=True)) / S.std(1, keepdims=True); n = len(Z); v = Z.sum(0)
    return float((v @ v - n * Z.shape[1]) / (n * (n - 1) * Z.shape[1]))


out = {"n_plateaus": int(len(PL)), "abs_slope_range": [float(PL[:, 2].min()), float(PL[:, 2].max())], "k_grid": [float(KG[0]), float(KG[-1])], "levels": {}}
fig, ax = plt.subplots(2, len(LV), figsize=(4.2 * len(LV), 8), facecolor=SURF, squeeze=False)
for j, L in enumerate(LV):
    dd = np.loadtxt(HERE / f"era_qbo_u{L}_point0E_daily.dat"); td, u = dd[:, 0], dd[:, 1]; y = hp120(u)
    if PLANT and L == PLANT_LEVEL:
        i = np.clip(np.searchsorted(t5, td) - 1, 0, len(t5) - 2); x0, x1, y0, y1 = t5[i], t5[i + 1], F5[i], F5[i + 1]
        Fd = np.where(np.abs(y1 - y0) / (x1 - x0) > 5, np.where(td < 0.5 * (x0 + x1), y0, y1), y0 + (y1 - y0) * (td - x0) / (x1 - x0))
        for t0, t1, sl in PL:                            # PLANT_JITTER: true plateau slopes differ from the model's
            m = (td >= t0 - 5 / 365.25) & (td <= t1 + 5 / 365.25); Fd[m] += JIT * rng.normal() * sl * np.sign(np.interp(t0, t5[1:], d)) * (td[m] - t0)
        sig = sum(np.sin(2 * np.pi * k * Fd) for k in (60, 110, 170)); y = y + sig * np.sqrt(PLANT * y.var() / sig.var())
    spec, qbo = [], []
    for t0, t1, s in PL:
        m = (td >= t0) & (td <= t1); x = y[m] - y[m].mean(); w = np.hanning(m.sum()); tt = td[m]
        P = np.abs(np.exp(-2j * np.pi * np.outer(FR, tt)) @ (w * x)) ** 2
        lp = np.log(P); lp = lp - np.polyval(np.polyfit(np.log(FR), lp, 1), np.log(FR)); spec.append(lp); qbo.append(np.sign(u[m].mean()))
    spec = np.array(spec); strata = SEAS * 2 + (np.array(qbo) > 0)
    to_k = lambda slopes: np.array([np.interp(KG * sl, FR, sp) for sp, sl in zip(spec, slopes)])
    c_true = C(to_k(PL[:, 2])); c_freq = C(np.array([np.interp(np.linspace(FR[0], FR[-1], 400), FR, sp) for sp in spec]))
    perm = []
    for _ in range(3000):
        sl = PL[:, 2].copy()
        for g in np.unique(strata):
            idx = np.where(strata == g)[0]; sl[idx] = sl[rng.permutation(idx)]
        perm.append(C(to_k(sl)))
    perm = np.array(perm); p = float(np.mean(perm >= c_true)); z = float((c_true - perm.mean()) / perm.std())
    out["levels"][L] = {"C_k_true": c_true, "C_k_perm_mean": float(perm.mean()), "C_k_perm_sd": float(perm.std()), "z": z, "p": p, "C_freq": c_freq,
                        "strata_sizes": [int((strata == g).sum()) for g in range(4)]}
    ak = ax[0, j]; mk = to_k(PL[:, 2]).mean(0); ak.plot(KG, mk, color="#d94a4a", lw=1.5, label="own slopes")
    sl = PL[:, 2].copy()
    for g in np.unique(strata):
        idx = np.where(strata == g)[0]; sl[idx] = sl[rng.permutation(idx)]
    ak.plot(KG, to_k(sl).mean(0), color=MUTED, lw=1, label="slopes permuted (one draw)")
    ak.set_title(f"{L} hPa: mean whitened spectrum vs k", loc="left", fontsize=10, color=INK); ak.set_xlabel("k = f / |s_p|", color=MUTED)
    if j == 0: ak.legend(fontsize=8, frameon=False); ak.set_ylabel("log power (whitened)", color=MUTED)
    ah = ax[1, j]; ah.hist(perm, bins=40, color="#c9c8c3"); ah.axvline(c_true, color="#d94a4a", lw=2)
    ah.set_title(f"C own {c_true:.4f}, perm {perm.mean():.4f}±{perm.std():.4f}\nz {z:+.2f}, p {p:.3f}", loc="left", fontsize=9.5, color=INK)
    ah.set_xlabel("cross-plateau consistency C (k coordinate)", color=MUTED)
    for a_ in (ak, ah): a_.set_facecolor(SURF); a_.tick_params(colors=MUTED)
tag = f"_plant{PLANT:g}" if PLANT else ""
if PLANT: LV = [PLANT_LEVEL]
fig.tight_layout(); fig.savefig(HERE / f"plateau_scaling{tag}.png", dpi=105, facecolor=SURF)
json.dump(out, open(HERE / f"plateau_scaling{tag}.json", "w"), indent=1)
print(f"{len(PL)} plateaus, |s| {PL[:, 2].min():.2f}-{PL[:, 2].max():.2f} /yr, common k {KG[0]:.0f}-{KG[-1]:.0f}")
for L, v in out["levels"].items():
    print(f"  {L:3d} hPa  C(k, own slopes) {v['C_k_true']:.4f}  vs permuted {v['C_k_perm_mean']:.4f} ± {v['C_k_perm_sd']:.4f}  z {v['z']:+.2f}  p {v['p']:.3f}   "
          f"C(plain f) {v['C_freq']:.4f}  strata {v['strata_sizes']}")
print("wrote plateau_scaling.png")
