#!/usr/bin/env python3
"""Where does zonal wavenumber > 0 emerge in the equatorial winds, and which manifold drives it? (2026-10-22)

Data: ERA5 monthly zonal wind, 5S-5N mean at each longitude (2.5 deg), 1-200 hPa, 1979-2025
(SAO/era5_u_eq_monthly.nc and era5_u_eq_monthly_low.nc). Anomalies: monthly climatology removed per longitude.
1. Zonal Fourier decomposition: share of anomaly variance in n = 0, 1, 2, 3 by level.
2. For each level and n = 0, 1, 2 (n >= 1: cos and sin parts), regression on
     QBO windings: the QBO two-leg fit's base 0.1108 with its harmonics 1, 4, 7, 8, 11, 14, 17, 24, 26 (legs/qbo30_point0E,
                   monthly leg manifold)
     ENSO windings: the common manifold's base 0.2017 with harmonics 1..9 (legs/amo_pdo_nino4_iode)
   (9 windings = 18 regressors each). Scores: training r (1979-2010) as a percentile among 300 phase-randomized
   surrogates of the target, and held-out r (2011-2022).
    python3 qbo_low_stellarator.py  -> qbo_low_stellarator.png, stdout
"""
import json
from pathlib import Path

import matplotlib
import netCDF4 as nc
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
SAO = HERE.parent / "experiments/Oct2026_generator/SAO"; LEGS = HERE.parent / "experiments/Oct2026_generator/legs"
INK, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
rng = np.random.default_rng(9)


def load(f):
    ds = nc.Dataset(f); tv = [k for k in ds.variables if "time" in k][0]; pv = [k for k in ds.variables if "pressure" in k][0]
    d = nc.num2date(ds.variables[tv][:], ds.variables[tv].units); T = np.array([x.year + (x.month - 1) / 12 for x in d])
    U = np.asarray(ds.variables["u"][:]).mean(axis=2)                       # time, level, lon
    o = np.argsort(T); return T[o], np.asarray(ds.variables[pv][:]), U[o], np.asarray(ds.variables["longitude"][:])


T, P1, U1, lon = load(SAO / "era5_u_eq_monthly.nc"); _, P2, U2, _ = load(SAO / "era5_u_eq_monthly_low.nc")
P = np.concatenate([P1, P2]); U = np.concatenate([U1, U2], axis=1); o = np.argsort(P); P, U = P[o], U[:, o]
mon = np.round((T % 1) * 12).astype(int) % 12
A = U.copy()
for m in range(12):
    A[mon == m] -= U[mon == m].mean(0)
Z = np.fft.rfft(A, axis=2) / A.shape[2]                                    # time, level, n

print("1. share of anomaly variance by zonal wavenumber")
print(f"{'hPa':>5s}  {'n=0':>5s} {'n=1':>5s} {'n=2':>5s} {'n=3':>5s} {'n>3':>5s}")
share = []
for j, p in enumerate(P):
    pw = np.abs(Z[:, j, :]) ** 2; pw[:, 1:] *= 2; tot = pw.sum(1).mean()
    s = [pw[:, n].mean() / tot for n in range(4)]; share.append(s + [1 - sum(s)])
    print(f"{p:5.0f}  " + " ".join(f"{v:5.2f}" for v in share[-1]))
share = np.array(share)


def manifold(path, col=3):
    a = np.loadtxt(path, delimiter=","); return a[:, 0], a[:, col]


qp = json.load(open(LEGS / "qbo30_point0E/lt.exe.p")); kq = abs(qp["ltep"][0]); hq = [1] + sorted(int(h) for h in qp["harm"])
ep = json.load(open(LEGS / "amo_pdo_nino4_iode/lt.exe.p")); ke = abs(ep["ltep"][-1]); he = list(range(1, 10))
tq, Fq = manifold(LEGS / "qbo30_point0E/lte_results.1.csv"); te_, Fe = manifold(LEGS / "amo_pdo_nino4_iode/lte_results.1.csv")
span = (T >= 1979) & (T < 2022.5); tt = T[span]
SETS = {"QBO windings": (np.interp(tt, tq, Fq), kq, hq), "ENSO windings": (np.interp(tt, te_, Fe), ke, he)}
X = {k: np.column_stack([f(2 * np.pi * h * kk * F) for h in hs for f in (np.cos, np.sin)] + [np.ones_like(tt)])
     for k, (F, kk, hs) in SETS.items()}
tr = tt < 2011.0; ho = ~tr


def score(y, Xm):
    b = np.linalg.lstsq(Xm[tr], y[tr], rcond=None)[0]; m = Xm @ b
    return np.corrcoef(m[tr], y[tr])[0, 1], np.corrcoef(m[ho], y[ho])[0, 1]


def surr(y, n):
    Zf = np.fft.rfft(y - y.mean()); ph = rng.uniform(0, 2 * np.pi, (n, len(Zf))); ph[:, 0] = 0
    return np.fft.irfft(np.abs(Zf) * np.exp(1j * ph), len(y), axis=1) + y.mean()


print("\n2. regression of each wavenumber component on QBO vs ENSO windings (train 1979-2010: r, surrogate pct; held-out 2011-2022: r)")
RES = {}
for j, p in enumerate(P):
    line = f"{p:5.0f} hPa"
    for n in range(3):
        comps = [Z[span, j, 0].real] if n == 0 else [Z[span, j, n].real, Z[span, j, n].imag]
        for k in SETS:
            rt, rh, pc = [], [], []
            for y in comps:
                a, b = score(y, X[k]); S = surr(y, 300)
                sr = np.array([score(s, X[k])[0] for s in S]); rt.append(a); rh.append(b); pc.append(np.mean(sr < a) * 100)
            RES[(p, n, k)] = (np.mean(rt), np.mean(pc), np.mean(rh))
        q, e = RES[(p, n, "QBO windings")], RES[(p, n, "ENSO windings")]
        line += f" | n={n}: QBO {q[0]:.2f} ({q[1]:3.0f}) {q[2]:+.2f}  ENSO {e[0]:.2f} ({e[1]:3.0f}) {e[2]:+.2f}"
    print(line)

fig, axs = plt.subplots(1, 3, figsize=(16, 5.2), facecolor=SURF)
ax = axs[0]
for n, c in zip(range(4), ("#0b0b0b", "#2a78d6", "#eb6834", "#1baf7a")):
    ax.plot(share[:, n], P, marker="o", color=c, lw=1.5, label=f"n = {n}")
ax.plot(share[:, 4], P, marker="o", color="#e87ba4", lw=1.2, label="n > 3")
ax.set_title("share of wind anomaly variance by zonal wavenumber", loc="left", fontsize=10.5, color=INK)
for col, n in ((1, 0), (2, 1)):
    ax = axs[col]
    for k, c in (("QBO windings", "#2a78d6"), ("ENSO windings", "#eb6834")):
        ax.plot([RES[(p, n, k)][2] for p in P], P, marker="o", color=c, lw=1.6, label=f"{k}, held-out r")
        ax.plot([RES[(p, n, k)][1] / 100 for p in P], P, marker="s", ls=":", color=c, lw=1.0, label=f"{k}, surrogate pct / 100")
    ax.axvline(0, color=MUTED, lw=0.6); ax.set_xlim(-0.6, 1.02)
    ax.set_title(f"zonal wavenumber n = {n}: which manifold?", loc="left", fontsize=10.5, color=INK)
axs[1].legend(frameon=False, fontsize=8)
axs[0].legend(frameon=False, fontsize=8.5)
for a in axs:
    a.set_yscale("log"); a.invert_yaxis(); a.set_yticks(P); a.set_yticklabels([f"{p:.0f}" for p in P])
    a.set_facecolor(SURF); a.tick_params(colors=MUTED); a.grid(True, color=GRID, lw=0.5)
    for sp in ("top", "right"):
        a.spines[sp].set_visible(False)
axs[0].set_ylabel("hPa", color=MUTED)
fig.tight_layout(); fig.savefig(HERE / "qbo_low_stellarator.png", dpi=95, facecolor=SURF)
print("plot: qbo_low_stellarator.png")
