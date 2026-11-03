#!/usr/bin/env python3
"""Which way does it move: QBO phase vs wave phase vs wave envelope, lag against altitude (2026-10-27).
The QBO regimes (zero-wind lines) descend; that is phase propagation. Vertically propagating equatorial waves (Kelvin,
gravity) carry energy and momentum upward while their phase descends, so the upward part should show in the wave
ENVELOPE (amplitude), not in any phase. Three quantities, as the lag of each level behind the level above
(positive = the lower level lags = downward; negative = upward), accumulated from 70 hPa:
  QBO phase       u low-passed at 200 d
  wave phase      u band-passed 4-20 d (Kelvin band)
  wave envelope   |Hilbert| of the 4-20 d band, high-passed (120 d daily, 60 d hourly) so that the QBO's own
                  modulation of wave amplitude (which descends with the shear zones) is removed
Daily data: ERA5 0N 0E 00 UTC (get_era5_point.py), levels 10 (if present) 20 30 50 70 hPa; lags from adjacent pairs,
over the whole record and over six ~8-yr blocks (spread). Peak lags refined by a parabola through the maximum.
Hourly data: ERA5 hourly u at 10, 30, 50, 70 hPa (adjacent pairs; panel c shows 30-50), 5S-5N at all 72 longitudes (../../../docs/era5/data, half-year files,
1999-2016 as they exist); cross-correlation averaged over grid points, per half-year and pooled. Group transit between
adjacent levels is expected to be under a day, which only the hourly data can resolve. The hourly data also allow a
Kelvin-wave space-time filter (eastward zonal wavenumbers 1-3, periods 4-20 d); its envelope is the modulus of the
one-sided (analytic) filtered field.
    python3 updown_lags.py   -> updown_lags.png, updown_lags.json"""
import json
from pathlib import Path

import matplotlib
import netCDF4 as nc
import numpy as np
from scipy.signal import butter, filtfilt, hilbert

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
HOURLY = HERE.parents[2] / "docs" / "era5" / "data"
INK, MUTED, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
C_QBO, C_PH, C_ENV = "#2a78d6", "#d94a4a", "#1d9a6c"
LV = [L for L in (10, 20, 30, 50, 70) if (HERE / f"era_qbo_u{L}_point0E_daily.dat").exists()]
TG = np.arange(1979.0, 2026.76, 1 / 365.25)
U = {L: np.interp(TG, *np.loadtxt(HERE / f"era_qbo_u{L}_point0E_daily.dat").T) for L in LV}


def ccf(x, y, M):
    """corr(x(t), y(t+l)) for l = -M..M; x, y are (n,) or (n, points), averaged over points."""
    x = (x - x.mean(0)) / x.std(0); y = (y - y.mean(0)) / y.std(0); n = len(x); L = np.arange(-M, M + 1)
    return L, np.array([np.mean(x[max(0, -l):n - max(0, l)] * y[max(0, l):n - max(0, -l)]) for l in L])


def peak(L, c):
    i = int(np.argmax(c))
    if 0 < i < len(c) - 1:
        d = c[i - 1] - 2 * c[i] + c[i + 1]
        return L[i] + (0.5 * (c[i - 1] - c[i + 1]) / d if d < 0 else 0.0), c[i]
    return float(L[i]), c[i]


def filt(x, kind, fs, *p):
    b, a = butter(3, [1 / q for q in p] if kind == "band" else 1 / p[0], btype=kind, fs=fs)
    return filtfilt(b, a, x - x.mean(0), axis=0)


def series(x, fs, env_hp):
    q = filt(x, "low", fs, 200) if fs == 1 else None
    w = filt(x, "band", fs, 20, 4)
    e = filt(np.abs(hilbert(w, axis=0)), "high", fs, env_hp)
    return q, w, e


# ---- daily, adjacent pairs, whole record and six blocks
NB = 6; edges = np.linspace(0, len(TG), NB + 1).astype(int)
daily = {"QBO phase": {}, "wave phase": {}, "wave envelope": {}}
for a, b in zip(LV[:-1], LV[1:]):
    sa, sb = series(U[a], 1, 120), series(U[b], 1, 120)
    for (name, M), xa, xb in zip([("QBO phase", 400), ("wave phase", 15), ("wave envelope", 40)], sa, sb):
        whole = peak(*ccf(xa, xb, M))
        blocks = [peak(*ccf(xa[s:e], xb[s:e], M))[0] for s, e in zip(edges[:-1], edges[1:])]
        daily[name][(a, b)] = {"lag_d": float(whole[0]), "r": float(whole[1]), "block_lags_d": [float(v) for v in blocks]}

# ---- hourly, adjacent pairs 10 <-> 30, 30 <-> 50, 50 <-> 70 hPa (panel c shows 30 <-> 50)
KEYS = ("wave phase", "wave envelope", "Kelvin phase", "Kelvin envelope")


def kelvin(x, nlon=72, kmax=3, p0=20, p1=4):
    """Space-time filter: eastward zonal wavenumbers 1..kmax, periods p1..p0 days (hourly x of shape (nt, nlat*nlon)).
    One-sided in (frequency, wavenumber), so the result is the analytic wave: real part = wave, modulus = envelope.
    Eastward cos(k lon - w t) sits at positive k and negative numpy frequency."""
    nt = len(x); X = np.fft.fft2(x.reshape(nt, -1, nlon) - x.mean(0).reshape(1, -1, nlon), axes=(0, 2))
    f = np.fft.fftfreq(nt, 1 / 24.0)[:, None, None]; k = np.fft.fftfreq(nlon, 1.0 / nlon)[None, None, :]
    keep = (k >= 1) & (k <= kmax) & (f <= -1.0 / p0) & (f >= -1.0 / p1)
    z = 2 * np.fft.ifft2(np.where(keep, X, 0), axes=(0, 2)).reshape(nt, -1)
    return z.real, filt(np.abs(z), "high", 24, 60)


def hourly_pair(up, lo):
    """Per-half-year and pooled cross-correlation of the hourly 5S-5N field at level lo behind level up."""
    hr = {k: [] for k in KEYS}; po = {k: [] for k in KEYS}; us = []
    for fa in sorted(HOURLY.glob(f"era5_u{up}_hourly_*.nc")):
        fb = HOURLY / fa.name.replace(f"u{up}_", f"u{lo}_")
        if not fb.exists():
            continue
        xs = []
        for f in (fa, fb):
            d = nc.Dataset(f); la = np.asarray(d.variables["latitude"][:]); j = np.abs(la) <= 5.01
            xs.append(np.asarray(d.variables["u"][:, 0, j, :], float).reshape(len(d.variables["valid_time"]), -1))
        if len(xs[0]) != len(xs[1]) or len(xs[0]) < 24 * 150:
            continue
        us.append(fa.name.split("_")[-1][:-3]); cut = slice(24 * 10, -24 * 10)
        _, wa, ea = series(xs[0], 24, 60); _, wb, eb = series(xs[1], 24, 60)
        ka, kea = kelvin(xs[0]); kb, keb = kelvin(xs[1])
        for name, x, y in (("wave phase", wa, wb), ("wave envelope", ea, eb), ("Kelvin phase", ka, kb), ("Kelvin envelope", kea, keb)):
            L, c = ccf(x[cut], y[cut], 96); hr[name].append(float(peak(L, c)[0])); po[name].append(c)
    return hr, {k: (np.arange(-96, 97), np.mean(v, 0)) for k, v in po.items() if v}, us


PAIRS = [(a, b) for a, b in ((10, 30), (30, 50), (50, 70))
         if any(HOURLY.glob(f"era5_u{a}_hourly_*.nc")) and any(HOURLY.glob(f"era5_u{b}_hourly_*.nc"))]
HP = {p: hourly_pair(*p) for p in PAIRS}
hourly, hourly_pooled, used = HP[(30, 50)]

# ---- figure
fig, ax = plt.subplots(1, 3, figsize=(17, 6.2), facecolor=SURF, gridspec_kw=dict(width_ratios=[1, 1, 1.15], wspace=0.28))


def cumulative(name, key="lag_d"):
    """lag of each level relative to 70 hPa, built up from the adjacent pairs (70 hPa = 0)."""
    t = {LV[-1]: 0.0}
    for a, b in reversed(list(zip(LV[:-1], LV[1:]))):
        v = daily[name][(a, b)][key]; t[a] = t[b] - v
    return np.array([t[L] for L in LV])


def draw(axx, names, title, xlab):
    for name, col in names:
        t = cumulative(name); axx.plot(t, LV, "-o", color=col, lw=2, ms=6, label=name)
        for k in range(NB):
            tb = {LV[-1]: 0.0}
            for a, b in reversed(list(zip(LV[:-1], LV[1:]))):
                tb[a] = tb[b] - daily[name][(a, b)]["block_lags_d"][k]
            axx.plot([tb[L] for L in LV], LV, "-", color=col, lw=0.7, alpha=0.35)
    axx.axvline(0, color=MUTED, lw=0.8); axx.set_yscale("log"); axx.invert_yaxis(); axx.set_yticks(LV); axx.set_yticklabels([str(L) for L in LV])
    axx.minorticks_off(); axx.set_ylabel("pressure (hPa)", color=MUTED); axx.set_xlabel(xlab, color=MUTED)
    axx.set_title(title, loc="left", fontsize=10.5, color=INK); axx.legend(fontsize=8.5, loc="lower left", frameon=False)
    axx.set_facecolor(SURF); axx.tick_params(colors=MUTED)


draw(ax[0], [("QBO phase", C_QBO)], "(a) QBO phase (u low-passed 200 d), daily ERA5",
     "time of arrival relative to 70 hPa (days)")
draw(ax[1], [("wave phase", C_PH), ("wave envelope", C_ENV)], "(b) 4–20 d waves: phase vs envelope, daily ERA5",
     "time of arrival relative to 70 hPa (days)")
for a in ax[:2]:
    a.text(0.98, 0.98, "left of 0 = arrives before 70 hPa\nline tilting up-left = descends\n(thin: six ~8-yr blocks)",
           transform=a.transAxes, ha="right", va="top", fontsize=8, color=MUTED)
a = ax[2]
for i, (name, col, ls) in enumerate((("wave phase", C_PH, ":"), ("wave envelope", C_ENV, ":"), ("Kelvin phase", C_PH, "-"), ("Kelvin envelope", C_ENV, "-"))):
    if name in hourly_pooled:
        L, c = hourly_pooled[name]; pk, r = peak(L, c)
        a.plot(L, c, color=col, lw=2 if ls == "-" else 1.2, ls=ls,
               label=f"{name}: peak {pk:+.0f} h, r {r:.2f} (half-year median {np.median(hourly[name]):+.0f} h, {sum(v < 0 for v in hourly[name])}/{len(hourly[name])} upward)")
        a.plot(hourly[name], np.full(len(hourly[name]), -0.12 - 0.05 * i), "|", color=col, ms=8)
a.axvline(0, color=MUTED, lw=0.8); a.set_xlabel("lag of 50 hPa behind 30 hPa (hours; + = downward, − = upward)", color=MUTED)
a.set_ylabel("cross-correlation (5S–5N, all longitudes)", color=MUTED); a.legend(fontsize=7.5, loc="upper left", frameon=False); a.set_ylim(-0.35, 1.05)
a.set_title(f"(c) hourly ERA5, 30 ↔ 50 hPa, 4–20 d (dotted: all; solid: Kelvin, east k=1–3), {len(used)} half-years", loc="left", fontsize=10.5, color=INK)
a.text(0.98, 0.04, "ticks: per-half-year peak lags", transform=a.transAxes, ha="right", fontsize=8, color=MUTED)
a.set_facecolor(SURF); a.tick_params(colors=MUTED)
fig.savefig(HERE / "updown_lags.png", dpi=110, facecolor=SURF, bbox_inches="tight")

out = {"daily_levels": LV, "daily": {n: {f"{a}->{b}": v for (a, b), v in d.items()} for n, d in daily.items()},
       **{f"hourly_{a}_{b}": {"half_years": us, "peak_lags_h": hr,
                               "pooled_peak_h": {k: float(peak(*v)[0]) for k, v in po.items()},
                               "pooled_peak_r": {k: float(peak(*v)[1]) for k, v in po.items()}}
          for (a, b), (hr, po, us) in HP.items()}}
json.dump(out, open(HERE / "updown_lags.json", "w"), indent=1)
print("lag of lower level behind upper (+ = downward, - = upward)")
for n, d in daily.items():
    print(f"  {n:14s} " + "  ".join(f"{a}->{b}: {v['lag_d']:+6.1f} d (r {v['r']:.2f}, blocks {min(v['block_lags_d']):+.1f}..{max(v['block_lags_d']):+.1f})" for (a, b), v in d.items()))
for (a, b), (hr, po, us) in HP.items():
    for k, v in hr.items():
        o = out[f"hourly_{a}_{b}"]
        print(f"  hourly {a}->{b} {k:14s} pooled {o['pooled_peak_h'][k]:+.1f} h (r {o['pooled_peak_r'][k]:.2f}), per half-year median {np.median(v):+.1f} h, "
              f"{sum(x < 0 for x in v)}/{len(v)} negative ({len(us)} half-years)")
print("wrote updown_lags.png")
