#!/usr/bin/env python3
"""Observed equatorial QBO time-height section at daily resolution from ERA5 alone (2026-10-27).
ERA5 has no 40 or 15 hPa level, so the section uses the ERA5 levels that exist (10, 20, 30, 50, 70 hPa; one grid point,
0N 0E, 00 UTC, from get_era5_point.py) and fills 40 and 15 hPa by linear interpolation in log-pressure. If the 10 hPa
file is not there yet, the top row (15 hPa) falls back to KIT monthly (PCHIP in time), and the plot says so.
Panels: (a) full record, (b) 2013-2023 zoom, (c) the same zoom high-passed (minus a 120-day running mean) to show the
sub-120-day "stripes". Printed: lag-0 correlation of the high-passed series between adjacent ERA levels, and the lag
of maximum correlation (positive = lower level lags), a measure of how vertically coherent the stripes are.
    python3 era_section.py   -> era_section.png"""
from pathlib import Path

import matplotlib
import numpy as np
from scipy.interpolate import PchipInterpolator

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
INK, MUTED, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
ERA = [L for L in (10, 20, 30, 50, 70) if (HERE / f"era_qbo_u{L}_point0E_daily.dat").exists()]
TG = np.arange(1979.0, 2026.76, 1 / 365.25)
rows = {}
for L in ERA:
    a = np.loadtxt(HERE / f"era_qbo_u{L}_point0E_daily.dat"); rows[L] = np.interp(TG, a[:, 0], a[:, 1], left=np.nan, right=np.nan)
src = {L: "ERA5" for L in ERA}
if 10 not in ERA:
    k = np.loadtxt(HERE / "kit_qbo_u15.dat"); rows[15] = PchipInterpolator(k[:, 0] + 1 / 24, k[:, 1], extrapolate=False)(TG); src[15] = "KIT monthly"
LV = np.array(sorted(rows)); U = np.array([rows[L] for L in LV])
P = np.geomspace(15, 70, 90)
Ug = np.array([np.interp(np.log(P), np.log(LV), U[:, i]) for i in range(len(TG))]).T


def highpass(x, days=120):
    w = np.ones(days) / days; good = np.isfinite(x); y = np.where(good, x, 0.0)
    return x - np.convolve(y, w, "same") / np.maximum(np.convolve(good.astype(float), w, "same"), 1e-9)


H = np.array([highpass(r) for r in Ug])
fig = plt.figure(figsize=(17, 13), facecolor=SURF)
gs = fig.add_gridspec(3, 1, height_ratios=[1, 1, 1], hspace=0.32)


def panel(ax, Z, t0, t1, title, vmax):
    m = (TG >= t0) & (TG <= t1)
    im = ax.pcolormesh(TG[m], P, Z[:, m], cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="gouraud")
    ax.set_yscale("log"); ax.invert_yaxis(); ax.set_yticks([70, 50, 40, 30, 20, 15]); ax.set_yticklabels(["70", "50", "40", "30", "20", "15"])
    ax.minorticks_off(); ax.set_xlim(t0, t1)
    for x in (2016.15, 2020.0):
        if t0 < x < t1: ax.axvline(x, color="#eda100", ls=":", lw=1.4)
    ax.set_ylabel("pressure (hPa)", color=MUTED); ax.set_title(title, loc="left", fontsize=10.5, color=INK)
    ax.set_facecolor(SURF); ax.tick_params(colors=MUTED); fig.colorbar(im, ax=ax, pad=0.02, label="u (m/s)")


lev = ", ".join(f"{L} {src[L]}" for L in LV)
panel(fig.add_subplot(gs[0]), Ug, 1979, 2026.76, f"(a) OBSERVED equatorial u, daily, 0N 0E 00 UTC (levels: {lev}; others interpolated in log-p)", 35)
panel(fig.add_subplot(gs[1]), Ug, 2013, 2023.5, "(b) 2013–2023 (orange: the 2016 and 2020 disruptions)", 35)
ax = fig.add_subplot(gs[2]); panel(ax, H, 2013, 2023.5, "(c) 2013–2023, high-passed (minus 120-day running mean): the sub-120-day stripes", 12)
ax.set_xlabel("year", color=MUTED)
fig.savefig(HERE / "era_section.png", dpi=110, facecolor=SURF, bbox_inches="tight")

print("levels:", {int(L): src[L] for L in LV})
hp = {L: highpass(rows[L]) for L in ERA}
for a, b in zip(ERA[:-1], ERA[1:]):
    x, y = hp[a], hp[b]; ok = np.isfinite(x) & np.isfinite(y); x, y = x[ok] - x[ok].mean(), y[ok] - y[ok].mean()
    lags = np.arange(-30, 31); c = [np.corrcoef(x[max(0, -l):len(x) - max(0, l)], y[max(0, l):len(y) - max(0, -l)])[0, 1] for l in lags]
    print(f"high-passed {a:2d} vs {b:2d} hPa: lag-0 r {c[30]:+.2f}, max r {max(c):+.2f} at lag {lags[int(np.argmax(c))]:+d} d; "
          f"sd {np.std(x):.1f} / {np.std(y):.1f} m/s")
print("wrote era_section.png")
