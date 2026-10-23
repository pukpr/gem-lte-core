#!/usr/bin/env python3
"""Emulated equatorial zonal wind below 70 hPa (150 hPa, 180 deg, central Pacific), 1880 to the manifold's end
(2026-10-22). The stellarator layer: the ENSO windings (base 0.2017, harmonics 1-9) on the common manifold
(legs/amo_pdo_nino4_iode, 1880-2023), amplitudes fitted to ERA5 (5S-5N mean at 180 deg, monthly anomaly,
3-month running mean) over 1979-2010; 2011-2022 is held out. The QBO is absent at this level (qbo_low_stellarator.py).
    python3 low_wind_emulation.py  -> low_wind_emulation.png
"""
import json
from pathlib import Path

import matplotlib
import netCDF4 as nc
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
EX = HERE.parent / "experiments/Oct2026_generator"
LEVEL, LON = 150, 180.0
INK, MUTED, GRID, SURF, BLUE, ORANGE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb", "#2a78d6", "#eb6834"

ds = nc.Dataset(EX / "SAO/era5_u_eq_monthly_low.nc")
tv = [k for k in ds.variables if "time" in k][0]; pv = [k for k in ds.variables if "pressure" in k][0]
d = nc.num2date(ds.variables[tv][:], ds.variables[tv].units); T = np.array([x.year + (x.month - 1) / 12 + 1 / 24 for x in d])
lon = np.asarray(ds.variables["longitude"][:]) % 360; j = int(np.argmin(abs(lon - LON))); k = list(np.asarray(ds.variables[pv][:])).index(LEVEL)
u = np.asarray(ds.variables["u"][:, k, :, j]).mean(1); o = np.argsort(T); T, u = T[o], u[o]
mon = np.round((T % 1) * 12 - 0.5).astype(int) % 12
ua = u - np.array([u[mon == m].mean() for m in range(12)])[mon]; ua = np.convolve(ua, np.ones(3) / 3, "same")

p = json.load(open(EX / "legs/amo_pdo_nino4_iode/lt.exe.p")); kb = abs(p["ltep"][-1])
man = np.loadtxt(EX / "legs/amo_pdo_nino4_iode/lte_results.1.csv", delimiter=","); tm, Fm = man[:, 0], man[:, 3]


def basis(t):
    F = np.interp(t, tm, Fm)
    return np.column_stack([f(2 * np.pi * h * kb * F) for h in range(1, 10) for f in (np.cos, np.sin)] + [np.ones_like(t)])


m = (T < 2022.5); tr = m & (T < 2011); ho = m & (T >= 2011)
b = np.linalg.lstsq(basis(T[tr]), ua[tr], rcond=None)[0]
te = np.arange(tm[0], tm[-1], 1 / 12); em = np.convolve(basis(te) @ b, np.ones(3) / 3, "same")
r_tr = np.corrcoef(np.interp(T[tr], te, em), ua[tr])[0, 1]; r_ho = np.corrcoef(np.interp(T[ho], te, em), ua[ho])[0, 1]
print(f"{LEVEL} hPa at {LON:.0f} deg: training r {r_tr:+.2f} (1979-2010), held-out r {r_ho:+.2f} (2011-2022)")

fig, axs = plt.subplots(2, 1, figsize=(16, 8), facecolor=SURF, gridspec_kw=dict(height_ratios=[1, 1]))
ax = axs[0]
ax.plot(te, em, color=BLUE, lw=0.9, label="emulated (ENSO windings on the common manifold)")
ax.plot(T, ua, color=INK, lw=0.9, alpha=0.75, label="ERA5 (1979-2025)")
ax.set_xlim(1880, 2026); ax.set_title(f"Equatorial zonal wind anomaly at {LEVEL} hPa, 180° (central Pacific): the wavenumber > 0 (stellarator) layer, 1880 to present",
                                      loc="left", fontsize=11, color=INK)
ax.legend(frameon=False, fontsize=9, loc="upper left", ncol=2)
ax = axs[1]; s = te >= 1975
ax.plot(te[s], em[s], color=BLUE, lw=1.4, label="emulated")
ax.plot(T, ua, color=INK, lw=1.1, alpha=0.8, label="ERA5")
ax.axvspan(2011, 2022.5, color=GRID, alpha=0.7, lw=0); ax.text(2011.3, ax.get_ylim()[1] * 0.85, f"held out, r = {r_ho:+.2f}", fontsize=9, color=INK)
ax.text(1976, ax.get_ylim()[1] * 0.85, f"training 1979-2010, r = {r_tr:+.2f}", fontsize=9, color=INK)
ax.set_xlim(1975, 2026); ax.set_title("1975 to present, with ERA5", loc="left", fontsize=10.5, color=INK); ax.legend(frameon=False, fontsize=9, loc="lower left")
for a in axs:
    a.axhline(0, color=MUTED, lw=0.6); a.set_ylabel("u anomaly (m/s), 3-month mean", color=MUTED)
    a.set_facecolor(SURF); a.tick_params(colors=MUTED); a.grid(True, color=GRID, lw=0.5)
    for sp in ("top", "right"):
        a.spines[sp].set_visible(False)
fig.tight_layout(); fig.savefig(HERE / "low_wind_emulation.png", dpi=95, facecolor=SURF)
print("plot: low_wind_emulation.png")
