#!/usr/bin/env python3
"""Frequency-domain warping: let each plateau choose its own spectral stretch, then ask what the stretches track
(2026-10-27). Companion to plateau_scaling.py, which imposes the manifold's plateau slopes and so fails if they are
off by more than ~15 %.
Each plateau p gets a stretch sigma_p (grid 0.10-0.50, the role of |s_p| in f = k sigma_p). Iterate: template =
mean of the other plateaus' whitened spectra in the k coordinate; sigma_p = argmax correlation of plateau p's
spectrum (mapped with sigma) with the template. 12 sweeps, start at the median slope.
Then compare the inferred sigma_p with
  - the manifold's |s_p| (Spearman; p from 3000 stratified permutations as in plateau_scaling.py), and
  - the plateau-mean u at that level (Doppler shifting by the QBO wind would stretch spectra too).
Final alignment quality C_fit (mean pairwise correlation) is compared with the same inversion run on each
plateau's spectrum circularly shifted by a random amount in log f (destroys any common pattern, keeps the freedom).
PLANT=<fraction> [PLANT_JITTER=<sd>] plants k = 60, 110, 170 at 30 hPa as in plateau_scaling.py.
    python3 plateau_slope_inversion.py [level ...]   -> plateau_slope_inversion.json (+ _plant tag)"""
import json
import os
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
RUN = HERE.parent / "legs" / "qbo9_kit_era70full_ex2015_2021"
rng = np.random.default_rng(11)
PLANT = float(os.environ.get("PLANT", 0)); JIT = float(os.environ.get("PLANT_JITTER", 0))
LV = [int(x) for x in sys.argv[1:]] or ([30] if PLANT else [L for L in (10, 20, 30, 50, 70) if (HERE / f"era_qbo_u{L}_point0E_daily.dat").exists()])
a = np.loadtxt(RUN / "lte_results.9.csv", delimiter=","); t5, F5 = a[:, 0], a[:, 3]
d = np.diff(F5) / np.diff(t5); st = np.where(np.abs(d) > 5)[0]; edges = np.r_[0, st + 1, len(t5)]
PL = []
for i0, i1 in zip(edges[:-1], edges[1:]):
    if i1 - i0 < 20: continue
    s = np.polyfit(t5[i0:i1], F5[i0:i1], 1)[0]
    if abs(s) > 0.1 and (t5[i1 - 1] - t5[i0]) * 365.25 >= 120: PL.append((t5[i0] + 5 / 365.25, t5[i1 - 1] - 5 / 365.25, s))
PL = np.array(PL); ABS = np.abs(PL[:, 2]); SEAS = (((PL[:, 0] % 1) * 12 >= 3) & ((PL[:, 0] % 1) * 12 < 9)).astype(int)
FR = np.linspace(3, 90, 600); SIG = np.arange(0.10, 0.5001, 0.005); KG = np.linspace(3 / 0.10, 90 / 0.50, 300)
KG = np.linspace(30, 180, 300)                                        # k range visible for every sigma in the grid


def hp120(x):
    w = np.ones(120) / 120; return x - np.convolve(x, w, "same") / np.convolve(np.ones_like(x), w, "same")


def z(v): return (v - v.mean()) / v.std()


def invert(spec):
    n = len(spec); sig = np.full(n, np.median(ABS))
    for _ in range(12):
        Km = np.array([z(np.interp(KG * s, FR, sp, left=np.nan, right=np.nan)) for sp, s in zip(spec, sig)])
        tot = np.nansum(Km, 0)
        for p in range(n):
            T = z(tot - np.nan_to_num(Km[p]))
            cand = [np.corrcoef(np.interp(KG * s, FR, spec[p]), T)[0, 1] for s in SIG]
            sig[p] = SIG[int(np.argmax(cand))]
            Km[p] = z(np.interp(KG * sig[p], FR, spec[p])); tot = np.nansum(Km, 0)
    Z = np.array([z(np.interp(KG * s, FR, sp)) for sp, s in zip(spec, sig)]); v = Z.sum(0)
    return sig, float((v @ v - n * Z.shape[1]) / (n * (n - 1) * Z.shape[1]))


out = {}
for L in LV:
    dd = np.loadtxt(HERE / f"era_qbo_u{L}_point0E_daily.dat"); td, u = dd[:, 0], dd[:, 1]; y = hp120(u)
    if PLANT:
        i = np.clip(np.searchsorted(t5, td) - 1, 0, len(t5) - 2); x0, x1, y0, y1 = t5[i], t5[i + 1], F5[i], F5[i + 1]
        Fd = np.where(np.abs(y1 - y0) / (x1 - x0) > 5, np.where(td < 0.5 * (x0 + x1), y0, y1), y0 + (y1 - y0) * (td - x0) / (x1 - x0))
        for t0, t1, sl in PL:
            m = (td >= t0 - 5 / 365.25) & (td <= t1 + 5 / 365.25); Fd[m] += JIT * rng.normal() * sl * (td[m] - t0)
        sg = sum(np.sin(2 * np.pi * k * Fd) for k in (60, 110, 170)); y = y + sg * np.sqrt(PLANT * y.var() / sg.var())
    spec, um = [], []
    for t0, t1, _ in PL:
        m = (td >= t0) & (td <= t1); x = y[m] - y[m].mean(); w = np.hanning(m.sum())
        lp = np.log(np.abs(np.exp(-2j * np.pi * np.outer(FR, td[m])) @ (w * x)) ** 2)
        spec.append(lp - np.polyval(np.polyfit(np.log(FR), lp, 1), np.log(FR))); um.append(u[m].mean())
    spec = np.array(spec); um = np.array(um); strata = SEAS * 2 + (um > 0)
    sig, cfit = invert(spec)
    rho = spearmanr(sig, ABS).correlation; rho_u = spearmanr(sig, um).correlation
    perm = []
    for _ in range(3000):
        sl = ABS.copy()
        for g in np.unique(strata):
            idx = np.where(strata == g)[0]; sl[idx] = sl[rng.permutation(idx)]
        perm.append(spearmanr(sig, sl).correlation)
    perm = np.array(perm)
    null_c = []
    for _ in range(8):                                                # inversion on pattern-destroyed spectra
        sh = np.array([np.roll(sp, rng.integers(40, 560)) for sp in spec]); null_c.append(invert(sh)[1])
    out[L] = {"spearman_sigma_vs_abs_slope": float(rho), "p_perm": float(np.mean(perm >= rho)), "spearman_sigma_vs_u": float(rho_u),
              "C_fit": cfit, "C_fit_null_mean": float(np.mean(null_c)), "C_fit_null_max": float(np.max(null_c)),
              "sigma_ratio_sd": float(np.std(sig / ABS) / np.mean(sig / ABS))}
    v = out[L]
    print(f"{L:3d} hPa  rho(sigma, |s|) {v['spearman_sigma_vs_abs_slope']:+.2f} (p {v['p_perm']:.3f})  rho(sigma, u) {v['spearman_sigma_vs_u']:+.2f}  "
          f"C_fit {v['C_fit']:.4f} vs null {v['C_fit_null_mean']:.4f} (max {v['C_fit_null_max']:.4f})", flush=True)
tag = f"_plant{PLANT:g}" + (f"_jit{JIT:g}" if JIT else "") if PLANT else ""
json.dump(out, open(HERE / f"plateau_slope_inversion{tag}.json", "w"), indent=1)
