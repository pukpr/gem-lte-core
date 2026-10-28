#!/usr/bin/env python3
"""Multi-altitude co-inversion of the plateau stretch: does the data itself define a shared plateau profile?
(2026-10-27)
Fast windings ride the plateau slopes; the model's manifold may have the wrong slopes (too sensitive to trust), but
all altitudes share one manifold, so a real fast-winding pattern would stretch each plateau's stripe spectrum by the
SAME factor at every level. Estimate each plateau's stretch sigma_p blind (no use of the model's slopes), jointly over
a group of levels (sigma_p maximises the summed correlation of each level's spectrum with that level's template from
the other plateaus; 12 sweeps from the median), separately for two disjoint level groups, and compare.
  agreement   Spearman(sigma_A, sigma_B) between groups; null = 3000 permutations within strata of
              season x plateau-length tercile (both affect spectral shape and resolution)
  manifold    Spearman of each group's sigma with the model's |s_p| (does the model have the right profile?)
Groups: far apart {10, 20} vs {50, 70} hPa (stripes barely correlate between them: adjacent-level high-passed lag-0
r is 0.01-0.16), and interleaved {10, 30, 70} vs {20, 50}. Spectra as in plateau_scaling.py (u minus 120-day
running mean, Hann, 3-90 cpy, whitened in log f).
LOCAL=1 (default): refine the model's manifold locally -- sigma_p = |s_p| (1 + eps_p), eps_p in [-0.3, 0.3], starting
at eps = 0; agreement is then tested on the corrections eps_p (manifold-independent: two disjoint level groups must
find the same correction). LOCAL=0: blind search over sigma 0.10-0.50 from the median.
SHIFT_DAYS=<n>: circularly shift every level's daily series by n days (the same for all levels) before cutting the
plateau windows -- the surrogate null: keeps shared weather between levels, window lengths, slopes and any method bias,
and breaks only the link between the data and the manifold's timing. Run joint_plateau_null.sh for the distribution.
Note: with LOCAL=1 sigma correlates with the model's |s| by construction; only the corrections' agreement is informative.
Spectra: sine multitaper, 5 tapers (much lower variance than one Hann taper on a ~170-day plateau).
PLANT=<fraction> [PLANT_JITTER=<sd>] plants k = 60, 110, 170 at every level on a manifold whose plateau slopes are
the model's times (1 + sd N(0,1)) per plateau -- the jittered slopes are the truth the inversion should find.
    python3 joint_plateau_inversion.py   -> joint_plateau_inversion[_plant...].json"""
import json
import os
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
RUN = HERE.parent / "legs" / "qbo9_kit_era70full_ex2015_2021"
rng = np.random.default_rng(23)
SHIFT = int(os.environ.get("SHIFT_DAYS", 0)); PLANT = float(os.environ.get("PLANT", 0)); JIT = float(os.environ.get("PLANT_JITTER", 0))
LV = [10, 20, 30, 50, 70]
GROUPS = {"far {10,20} vs {50,70}": ([10, 20], [50, 70]), "interleaved {10,30,70} vs {20,50}": ([10, 30, 70], [20, 50])}
a = np.loadtxt(RUN / "lte_results.9.csv", delimiter=","); t5, F5 = a[:, 0], a[:, 3]
d = np.diff(F5) / np.diff(t5); st = np.where(np.abs(d) > 5)[0]; edges = np.r_[0, st + 1, len(t5)]
PL = []
for i0, i1 in zip(edges[:-1], edges[1:]):
    if i1 - i0 < 20: continue
    s = np.polyfit(t5[i0:i1], F5[i0:i1], 1)[0]
    if abs(s) > 0.1 and (t5[i1 - 1] - t5[i0]) * 365.25 >= 120: PL.append((t5[i0] + 5 / 365.25, t5[i1 - 1] - 5 / 365.25, s))
PL = np.array(PL); ABS = np.abs(PL[:, 2]); n = len(PL)
SEAS = (((PL[:, 0] % 1) * 12 >= 3) & ((PL[:, 0] % 1) * 12 < 9)).astype(int)
LEN = np.digitize(PL[:, 1] - PL[:, 0], np.quantile(PL[:, 1] - PL[:, 0], [1 / 3, 2 / 3])); STRATA = SEAS * 3 + LEN
JF = 1 + JIT * rng.normal(size=n); TRUE = ABS * JF                    # planted truth (= model slopes if no jitter)
FR = np.linspace(3, 90, 600); KG = np.linspace(30, 180, 300)
LOCAL = os.environ.get("LOCAL", "1") == "1"
EPS = np.arange(-0.30, 0.3001, 0.01)
SIGP = np.array([ABS[p] * (1 + EPS) for p in range(n)]) if LOCAL else np.tile(np.arange(0.10, 0.5001, 0.005), (n, 1))  # per-plateau grid


def hp120(x):
    w = np.ones(120) / 120; return x - np.convolve(x, w, "same") / np.convolve(np.ones_like(x), w, "same")


def z(v): return (v - v.mean()) / v.std()


SPEC = {}
for L in LV:
    dd = np.loadtxt(HERE / f"era_qbo_u{L}_point0E_daily.dat"); td, u = dd[:, 0], dd[:, 1]; y = hp120(u)
    if SHIFT: y = np.roll(y, SHIFT)
    if PLANT:
        i = np.clip(np.searchsorted(t5, td) - 1, 0, len(t5) - 2); x0, x1, y0, y1 = t5[i], t5[i + 1], F5[i], F5[i + 1]
        Fd = np.where(np.abs(y1 - y0) / (x1 - x0) > 5, np.where(td < 0.5 * (x0 + x1), y0, y1), y0 + (y1 - y0) * (td - x0) / (x1 - x0))
        for (t0, t1, sl), jf in zip(PL, JF):
            m = (td >= t0 - 5 / 365.25) & (td <= t1 + 5 / 365.25); Fd[m] += (jf - 1) * sl * (td[m] - t0)
        sg = sum(np.sin(2 * np.pi * k * Fd + rng.uniform(0, 2 * np.pi)) for k in (60, 110, 170)); y = y + sg * np.sqrt(PLANT * y.var() / sg.var())
    sp = []
    for t0, t1, _ in PL:
        m = (td >= t0) & (td <= t1); x = y[m] - y[m].mean(); N = m.sum(); E = np.exp(-2j * np.pi * np.outer(FR, td[m]))
        tap = [np.sqrt(2 / (N + 1)) * np.sin(np.pi * (j + 1) * np.arange(1, N + 1) / (N + 1)) for j in range(5)]
        lp = np.log(np.mean([np.abs(E @ (w * x)) ** 2 for w in tap], 0))
        sp.append(lp - np.polyval(np.polyfit(np.log(FR), lp, 1), np.log(FR)))
    SPEC[L] = np.array(sp)


def invert(levels):
    G = SIGP.shape[1]; idx = np.full(n, G // 2) if LOCAL else np.full(n, int(np.argmin(np.abs(SIGP[0] - np.median(ABS)))))
    MAP = {L: np.array([[z(np.interp(KG * s, FR, sp)) for s in SIGP[p]] for p, sp in enumerate(SPEC[L])]) for L in levels}   # (plateau, grid, k)
    for _ in range(12):
        for p in range(n):
            score = np.zeros(G)
            for L in levels:
                cur = MAP[L][np.arange(n), idx]; T = z(cur.sum(0) - cur[p]); score += MAP[L][p] @ T / len(KG)
            idx[p] = int(np.argmax(score))
    return SIGP[np.arange(n), idx]


def perm_p(x, y, nperm=3000):
    r0 = spearmanr(x, y).correlation; c = 0
    for _ in range(nperm):
        yy = y.copy()
        for g in np.unique(STRATA):
            ii = np.where(STRATA == g)[0]; yy[ii] = yy[rng.permutation(ii)]
        c += spearmanr(x, yy).correlation >= r0
    return float(r0), float(c / nperm)


out = {"n_plateaus": n, "plant": PLANT, "jitter": JIT, "groups": {}}
print(f"{n} plateaus; PLANT {PLANT:g} JITTER {JIT:g}")
for name, (A, B) in GROUPS.items():
    sA, sB = invert(A), invert(B)
    agree = perm_p(sA / ABS, sB / ABS) if LOCAL else perm_p(sA, sB); mA = perm_p(sA, ABS); mB = perm_p(sB, ABS)
    res = {"agreement_rho_p": agree, "A_vs_model_rho_p": mA, "B_vs_model_rho_p": mB}
    if PLANT:
        res["A_vs_truth_rho"] = float(spearmanr(sA, TRUE).correlation); res["B_vs_truth_rho"] = float(spearmanr(sB, TRUE).correlation)
        if LOCAL and JIT: res["A_corr_vs_true_corr"] = float(spearmanr(sA / ABS, JF).correlation); res["B_corr_vs_true_corr"] = float(spearmanr(sB / ABS, JF).correlation)
    out["groups"][name] = res
    print(f"  {name}: agreement rho {agree[0]:+.2f} (p {agree[1]:.3f}) | vs model |s|: A {mA[0]:+.2f} (p {mA[1]:.3f}), B {mB[0]:+.2f} (p {mB[1]:.3f})"
          + (f" | vs planted truth: A {res['A_vs_truth_rho']:+.2f}, B {res['B_vs_truth_rho']:+.2f}" if PLANT else "")
          + (f" | corrections vs true corrections: A {res['A_corr_vs_true_corr']:+.2f}, B {res['B_corr_vs_true_corr']:+.2f}" if PLANT and LOCAL and JIT else ""), flush=True)
sAll = invert(LV); out["all_levels_vs_model_rho_p"] = perm_p(sAll, ABS); out["local"] = LOCAL
out["all_levels_correction_sd"] = float(np.std(sAll / ABS - 1)); out["all_levels_corrections"] = (sAll / ABS - 1).round(3).tolist()
if PLANT: out["all_levels_vs_truth_rho"] = float(spearmanr(sAll, TRUE).correlation)
if PLANT and LOCAL and JIT: out["all_levels_corr_vs_true_corr"] = float(spearmanr(sAll / ABS, JF).correlation)
print(f"  all five levels jointly vs model |s|: rho {out['all_levels_vs_model_rho_p'][0]:+.2f} (p {out['all_levels_vs_model_rho_p'][1]:.3f})"
      + (f" | vs planted truth {out['all_levels_vs_truth_rho']:+.2f}" if PLANT else "")
      + (f" | corrections vs true {out['all_levels_corr_vs_true_corr']:+.2f}" if PLANT and LOCAL and JIT else ""))
tag = (f"_shift{SHIFT}" if SHIFT else "") + ("" if LOCAL else "_blind") + f"_plant{PLANT:g}" + (f"_jit{JIT:g}" if JIT else "") if PLANT else (f"_shift{SHIFT}" if SHIFT else "") + ("" if LOCAL else "_blind")
json.dump(out, open(HERE / f"joint_plateau_inversion{tag}.json", "w"), indent=1)
