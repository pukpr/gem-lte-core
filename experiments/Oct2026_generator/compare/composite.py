#!/usr/bin/env python3
"""Tidal-forcing composites before the comb/integrator/layer (2026-10-06):
  TABLE: sum_j A_j cos(2 pi (YL/P_j) t + phi_j) - shfT f_j sin(...)  [lt.exe Tide_Sum, Cos_Phase=True, as used
         in Calc_Forcing], YL = 365.2422484 + YEAR + p.year, from the 0.833 TABLE fit (table.p);
         also its sin form (the convention the dLOD calibration/gate uses).
  ZONAL: kappa * dV/dt(JD_Of(t) - tau), the Meeus generator as in GEM.Zonal v3 (zonal.p's YEAR clock).
Daily grid over the dLOD record (1962-2019) and over 1880-2023."""
import json, sys, numpy as np
sys.path.insert(0, "/home/paul/eval/gem-lte-core/experiments/Oct2026_generator")
import zonal_full as z
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
YEAR = 0.00405305
pT = json.load(open("table.p")); L = np.array(pT["lpap"]); P, A, PH = L[:, 0], L[:, 1], L[:, 2]
YL_T = 365.2422484 + YEAR + pT["year"]; shf = pT["shfT"]
def table(t, cos=True):
    f = YL_T / P; arg = 2 * np.pi * np.outer(t, f) + PH
    return ((A * np.cos(arg) - shf * f * np.sin(arg)) if cos else (A * np.sin(arg) + shf * f * np.cos(arg))).sum(axis=1)
def zonal(t): return -5.866344 * z.rate(z.jd_of(t, YEAR) - 0.92)
d = np.loadtxt("/home/paul/eval/gem-lte-core/experiments/Feb2026/dlod3.dat"); tL, yL = d[:, 0], d[:, 1]
def lagcc(a, b, L):
    n = len(a); return np.corrcoef(a[L:], b[:n - L])[0, 1] if L >= 0 else np.corrcoef(a[:n + L], b[-L:])[0, 1]
for span, t in (("dLOD record 1962-2019 (daily)", tL), ("1880-2023 (daily)", 1880 + np.arange(int(143.5 * 365.2422)) / 365.2422)):
    Z = zonal(t); Tc = table(t, True); Ts = table(t, False)
    print(f"\n== {span}")
    for lab, T in (("TABLE cos form (climate forcing)", Tc), ("TABLE sin form (dLOD convention)", Ts)):
        lags = np.arange(-10, 11); c = [lagcc(Z, T, L) for L in lags]; k = int(np.argmax(np.abs(c)))
        X = np.column_stack([T, np.ones_like(T)]); cf, *_ = np.linalg.lstsq(X, Z, rcond=None); res = Z - X @ cf
        print(f"  ZONAL vs {lab}: CC lag0 {c[10]:+.4f}; best lag {lags[k]:+d} d (CC {c[k]:+.4f}); amplitude ratio std(ZONAL)/std(TABLE) {Z.std()/T.std():.3f}; residual after scale {100*res.var()/Z.var():.1f}% of ZONAL")
    if span.startswith("dLOD"):
        print(f"  vs observed dLOD: ZONAL {np.corrcoef(Z, yL)[0,1]:+.4f} | TABLE sin form {np.corrcoef(Ts, yL)[0,1]:+.4f} | TABLE cos form {np.corrcoef(Tc, yL)[0,1]:+.4f}")
        # per-constituent comparison in lt.exe's basis (sin convention), ZONAL vs TABLE
        f = YL_T / P; X = np.column_stack([np.sin(2*np.pi*np.outer(t, f)), np.cos(2*np.pi*np.outer(t, f)), np.ones_like(t)])
        c, *_ = np.linalg.lstsq(X, Z, rcond=None); zA = np.hypot(c[:42], c[42:84]); zP = np.arctan2(c[42:84], c[:42])
        s = (zA * np.exp(1j*zP)) @ np.conj(A*np.exp(1j*PH)) / ((A*np.exp(1j*PH)) @ np.conj(A*np.exp(1j*PH)))
        print(f"  per-line, in the sin convention (one complex scale ZONAL -> TABLE: |s| {abs(s):.3f}, angle {np.angle(s):+.2f} rad):")
        wrap = lambda x: (x + np.pi) % (2*np.pi) - np.pi
        for j in np.argsort(-A)[:12]:
            zz = zA[j]*np.exp(1j*zP[j]) / s
            print(f"    {P[j]:9.3f} d: TABLE amp {A[j]:.4f} | ZONAL amp {abs(zz):.4f} ratio {abs(zz)/A[j]:.2f} | phase diff {wrap(np.angle(zz)-PH[j]):+.2f} rad")
        Tw = Ts[:400]; Zw = Z[:400]; tw = t[:400]
        fig, ax = plt.subplots(2, 1, figsize=(12, 6.5), gridspec_kw=dict(hspace=0.35))
        ax[0].plot(tw, (Zw - Zw.mean())/Zw.std(), color="#a05fd0", lw=1.0, label="ZONAL (Meeus, kappa dV/dt)")
        ax[0].plot(tw, (Tw - Tw.mean())/Tw.std(), color="#eb6834", lw=1.0, label="TABLE fit, sin form")
        ax[0].plot(tw, (yL[:400] - yL[:400].mean())/yL[:400].std(), color="#2a78d6", lw=0.8, alpha=0.7, label="observed dLOD")
        ax[0].set_title("first 400 days of the dLOD record (standardized)", loc="left", fontsize=10); ax[0].legend(frameon=False, fontsize=8, ncol=3); ax[0].grid(alpha=.3)
        X = np.column_stack([Ts, np.ones_like(Ts)]); cf, *_ = np.linalg.lstsq(X, Z, rcond=None); res = Z - X @ cf
        fr = np.fft.rfftfreq(len(res), 1.0); Pr = np.abs(np.fft.rfft(res - res.mean()))**2; Pz = np.abs(np.fft.rfft(Z - Z.mean()))**2
        ax[1].semilogy(fr[1:], Pz[1:], color="#a05fd0", lw=0.6, label="ZONAL"); ax[1].semilogy(fr[1:], Pr[1:], color="#555", lw=0.6, label="ZONAL - scaled TABLE (sin form)")
        ax[1].set_xlim(0, 0.2); ax[1].set_xlabel("frequency (cycles/day)"); ax[1].legend(frameon=False, fontsize=8); ax[1].grid(alpha=.3)
        ax[1].set_title("power spectrum, 1962-2019", loc="left", fontsize=10)
        fig.savefig("composite_zonal_vs_table.png", dpi=110, bbox_inches="tight")
print("\nplot: compare/composite_zonal_vs_table.png")
