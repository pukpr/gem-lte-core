#!/usr/bin/env python3
"""Subsurface (thermocline) response of the tropical Pacific to the tidally forced ENSO model (2026-10-21).

Data: ORAS5 depth of the 20 C isotherm, monthly 1958-2025, 1 x 1 degree boxes 15S-15N, 120E-285E
(experiments/Oct2026_generator/ORAS5/d20_<year>.npz from fetch_d20.py). Anomalies: monthly climatology and a quadratic
trend removed per box; a 3-month running mean. Sign: D20 deeper = positive.
References:
  FORCED   the fitted NINO3.4 model output of the four-region Nino fit (legs/nino4_34_3_12, lte_results.2.csv column 2):
           tidal forcing only, no wind or ocean data enter it.
  OBSERVED Nino3.4 SST (the same run's data column), for comparison.
For each box: correlation with each reference at lags -24..+24 months (positive lag = D20 follows the reference);
the lag of the largest positive correlation (and of the largest negative). Along the equator (2S-2N) and along
off-equatorial bands, the slope of lag against longitude gives the propagation direction and speed:
  eastward at the equator = Kelvin-like; westward off the equator = Rossby-like; together a closed circuit.
Span 1958.0-2022.5 (the manifold's end).
    python3 enso_subsurface_d20.py   -> enso_subsurface_d20.png, digest on stdout
"""
import glob
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
EX = HERE.parent / "experiments"
ORAS = EX / "Oct2026_generator/ORAS5"
NINO = EX / "Oct2026_generator/legs/nino4_34_3_12"
INK, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
LAGS = np.arange(-24, 25)
T0, T1 = 1958.0, 2022.5

files = sorted(glob.glob(str(ORAS / "d20_*.npz")))
T = np.concatenate([np.load(f)["time"] for f in files]); D = np.concatenate([np.load(f)["d20"] for f in files])
LAT, LON = np.load(files[0])["lat"], np.load(files[0])["lon"]
o = np.argsort(T); T, D = T[o], D[o]
keep = (T >= T0 - 1e-6) & (T < T1); T, D = T[keep], D[keep]
mon = np.round((T % 1) * 12).astype(int) % 12
A = D.copy()
for m in range(12):
    A[mon == m] -= np.nanmean(D[mon == m], axis=0)
tt = T - T.mean(); X = np.column_stack([np.ones_like(tt), tt, tt ** 2])
flat = A.reshape(len(T), -1); good = np.all(np.isfinite(flat), axis=0)
flat[:, good] -= X @ np.linalg.lstsq(X, flat[:, good], rcond=None)[0]
k3 = np.ones(3) / 3
flat[:, good] = np.apply_along_axis(lambda v: np.convolve(v, k3, "same"), 0, flat[:, good])
A = flat.reshape(A.shape)

nino = np.loadtxt(NINO / "lte_results.2.csv", delimiter=",")
key = {round(x, 3): i for i, x in enumerate(nino[:, 0])}
idx = np.array([key.get(round(x, 3), -1) for x in T]); assert (idx >= 0).all(), "manifold dates do not cover the D20 span"
REF = {"FORCED (tidal model, Nino3.4)": nino[idx, 1], "OBSERVED Nino3.4": nino[idx, 2]}
for k in REF:
    v = REF[k] - np.polyval(np.polyfit(tt, REF[k], 2), tt); REF[k] = np.convolve(v, k3, "same")


def lagcorr(ref):
    z = (ref - ref.mean()) / ref.std(); n = len(z)
    Y = flat[:, good]; Y = (Y - Y.mean(0)) / Y.std(0)
    C = np.full((len(LAGS), flat.shape[1]), np.nan)
    for i, L in enumerate(LAGS):
        a, b = (z[:n - L], Y[L:]) if L >= 0 else (z[-L:], Y[:n + L])
        C[i, good] = (a[:, None] * b).mean(0)
    return C.reshape(len(LAGS), len(LAT), len(LON))


def band_lag(C, lat0, lat1, sign=+1):
    """per longitude: lag of the extreme correlation (sign) averaged over the latitude band, and that correlation"""
    m = (LAT >= lat0) & (LAT <= lat1); c = sign * np.nanmean(C[:, m, :], axis=1)        # (lags, lon)
    i = np.nanargmax(np.where(np.isfinite(c), c, -9), axis=0)
    return LAGS[i].astype(float), sign * c[i, np.arange(len(LON))]


BANDS = [("equator 2S-2N", -2, 2), ("5-8N", 5, 8), ("5-8S", -8, -5)]
res = {}
print(f"D20 span {T[0]:.2f}-{T[-1]:.2f} ({len(T)} months); boxes with complete data {good.sum()} of {good.size}")
for rk, ref in REF.items():
    C = lagcorr(ref); res[rk] = C
    print(f"\n== reference: {rk}")
    for lab, a, b in BANDS:
        for sgn, nm in ((+1, "deepening"), (-1, "shoaling")):
            lag, r = band_lag(C, a, b, sgn); sel = np.isfinite(r) & (np.abs(r) >= 0.25) & (LON < 280)
            if sel.sum() >= 10:
                s, c0 = np.polyfit(LON[sel], lag[sel], 1)
                spd = (111e3 * np.cos(np.radians((a + b) / 2)) / (s * 30.4 * 86400)) if s != 0 else np.inf
                print(f"   {lab:14s} {nm:10s}: |r| up to {np.nanmax(np.abs(r[sel])):.2f} over {sel.sum()} lon; lag from {lag[sel][0]:+.0f} mo (west) to {lag[sel][-1]:+.0f} mo (east); "
                      f"slope {s:+.3f} mo/deg -> {'eastward' if s > 0 else 'westward'} {abs(spd):.2f} m/s")
            else:
                print(f"   {lab:14s} {nm:10s}: too few longitudes with |r| >= 0.25 ({sel.sum()})")

fig, axs = plt.subplots(2, 2, figsize=(16, 9), facecolor=SURF)
for col, (rk, C) in enumerate(res.items()):
    i0 = list(LAGS).index(0)
    ax = axs[0, col]; im = ax.pcolormesh(LON, LAT, C[i0], cmap="RdBu_r", vmin=-0.8, vmax=0.8, shading="auto")
    ax.set_title(f"{rk}: correlation with D20 at lag 0", loc="left", fontsize=10, color=INK); ax.set_ylabel("latitude", color=MUTED)
    fig.colorbar(im, ax=ax, shrink=0.8)
    ax = axs[1, col]; m = (LAT >= -2) & (LAT <= 2); im = ax.pcolormesh(LON, LAGS, np.nanmean(C[:, m, :], axis=1), cmap="RdBu_r", vmin=-0.8, vmax=0.8, shading="auto")
    ax.set_title(f"{rk}: equatorial (2S-2N) D20 correlation by lag", loc="left", fontsize=10, color=INK)
    ax.set_ylabel("lag (months; + = D20 follows the reference)", color=MUTED); ax.axhline(0, color=INK, lw=0.6)
    fig.colorbar(im, ax=ax, shrink=0.8)
for ax in axs.ravel():
    ax.set_facecolor(SURF); ax.tick_params(colors=MUTED)
    ax.set_xticks([140, 160, 180, 200, 220, 240, 260, 280]); ax.set_xticklabels(["140E", "160E", "180", "160W", "140W", "120W", "100W", "80W"])
fig.suptitle("ORAS5 thermocline depth (20 °C isotherm) vs the tidally forced ENSO model and vs observed Niño3.4", x=0.01, ha="left", fontsize=12, color=INK)
fig.tight_layout(); fig.savefig(HERE / "enso_subsurface_d20.png", dpi=100, facecolor=SURF)
print("\nplot: enso_subsurface_d20.png")
