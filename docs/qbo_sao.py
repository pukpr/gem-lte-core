#!/usr/bin/env python3
"""SAO above, QBO below: the calendar clock and where the 2.37-yr beat takes hold (2026-10-22).

Data: ERA5 monthly zonal wind, zonal mean over 5S-5N, 1-70 hPa, 1979-2025
(experiments/Oct2026_generator/SAO/era5_u_eq_monthly.nc, fetch_sao.py).
1. SAO: per level, share of variance in the mean annual cycle and in its semiannual harmonic; year-to-year phase
   lock of the semiannual harmonic (Rayleigh R over years); inequality of the two cycles in the climatology.
2. QBO: per level, the wind minus its climatology, 5-month running mean; reversals with +-3 m/s hysteresis; their
   phase in the 2.370-yr beat (Rayleigh R, mean phase) as in qbo_phase_lock.py; interval-shuffle null for R.
   The mean phase against height gives the descent and the level where the reversal phase equals the drive's beat
   maximum for either sign of the drive (qbo_abs_phase.py: 30 hPa lags 21.8 or 7.6 months).
3. Month of QBO westerly onsets at 10/7/5 hPa against the SAO's westerly months.
    python3 qbo_sao.py  -> qbo_sao.png, stdout
"""
from pathlib import Path

import matplotlib
import netCDF4 as nc
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
F = HERE.parent / "experiments/Oct2026_generator/SAO/era5_u_eq_monthly.nc"
F0 = (365.2422 / 27.212221) % 1.0; P0 = 1 / F0
INK, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
rng = np.random.default_rng(5)

ds = nc.Dataset(F)
tv = [k for k in ds.variables if "time" in k][0]; pv = [k for k in ds.variables if "pressure" in k or k == "level"][0]
times = nc.num2date(ds.variables[tv][:], ds.variables[tv].units)
T = np.array([d.year + (d.month - 1) / 12 for d in times]); P = np.asarray(ds.variables[pv][:])
U = np.asarray(ds.variables["u"][:]).mean(axis=(2, 3))          # (time, level): zonal and 5S-5N mean
o = np.argsort(T); T, U = T[o], U[o]; lv = np.argsort(P); P, U = P[lv], U[:, lv]
mon = np.round((T % 1) * 12).astype(int) % 12; yrs = np.floor(T).astype(int)
clim = np.array([[U[mon == m, j].mean() for m in range(12)] for j in range(len(P))])      # level x month

print("1. SAO / annual cycle by level")
print(f"{'hPa':>5s} {'clim share':>10s} {'semiann share':>13s} {'semiann R over yrs':>18s} {'SAO W peaks (m/s, months)':>28s} {'E troughs':>20s}")
for j, p in enumerate(P):
    u = U[:, j]; cl = clim[j][mon]
    h2 = np.fft.rfft(clim[j])[2]; sh2 = 2 * abs(h2 / 12) ** 2 / 2 / u.var()
    ph = []
    for y in np.unique(yrs):
        m = yrs == y
        if m.sum() == 12: ph.append(np.angle(np.sum(u[m] * np.exp(-2j * np.pi * 2 * np.arange(12) / 12))))
    R = abs(np.mean(np.exp(1j * np.array(ph))))
    c = clim[j]; w1, w2 = c[1:6].argmax() + 1, c[7:12].argmax() + 7; e1, e2 = np.r_[c[10:], c[:4]].argmin(), c[4:10].argmin() + 4
    e1 = (e1 + 10) % 12
    print(f"{p:5.0f} {cl.var() / u.var():10.2f} {sh2:13.2f} {R:18.2f}   {c[w1]:+5.1f} ({w1 + 1:2d}) / {c[w2]:+5.1f} ({w2 + 1:2d})   {c[e1]:+5.1f} ({e1 + 1:2d}) / {c[e2]:+5.1f} ({e2 + 1:2d})")


def reversals(t, u, thr=3.0):
    out = []; state = np.sign(u[0]) or 1.0; last0 = t[0]
    for i in range(1, len(u)):
        if np.sign(u[i]) != np.sign(u[i - 1]): last0 = t[i - 1] + (t[i] - t[i - 1]) * u[i - 1] / (u[i - 1] - u[i])
        if state < 0 and u[i] > thr: out.append((last0, 1)); state = 1.0
        elif state > 0 and u[i] < -thr: out.append((last0, -1)); state = -1.0
    return np.array([a for a, _ in out]), np.array([b for _, b in out])


print("\n2. QBO (deseasonalized, 5-month mean): lock to the 2.370-yr beat by level")
rows = []
for j, p in enumerate(P):
    a = U[:, j] - clim[j][mon]; a = np.convolve(a, np.ones(5) / 5, "same")
    t, s = reversals(T[2:-2], a[2:-2])
    z = np.mean(np.exp(2j * np.pi * (F0 * (t - 1953.0) + np.where(s > 0, 0, 0.5))))
    d = np.diff(t); nul = []
    for _ in range(2000):
        tt = t[0] + np.concatenate([[0.0], np.cumsum(rng.permutation(d))]); ss = s[0] * (-1.0) ** np.arange(len(tt))
        nul.append(abs(np.mean(np.exp(2j * np.pi * (F0 * (tt - 1953.0) + np.where(ss > 0, 0, 0.5))))))
    ph = (np.angle(z) / (2 * np.pi)) % 1
    rows.append((p, len(t), abs(z), np.mean(np.array(nul) >= abs(z)), ph, a.std()))
    print(f"{p:5.0f} hPa: {len(t):3d} reversals, QBO sd {a.std():5.1f} m/s, R {abs(z):.2f} (shuffle p {rows[-1][3]:.3f}), mean phase {ph:.2f} cycles")
ph30 = [r[4] for r in rows if r[0] == 30][0]
print("\n   phase relative to 30 hPa, in months (negative = earlier, i.e. above):")
print("   " + "  ".join(f"{r[0]:.0f}:{(((r[4] - ph30 + 0.5) % 1) - 0.5) * P0 * 12:+.1f}" for r in rows))
for lab, lag30 in (("drive sign as written", 21.8), ("drive sign reversed", 7.6)):
    print(f"   trigger for '{lab}': {lag30:.1f} months above 30 hPa in time; levels reach that lead at: ",
          ", ".join(f"{r[0]:.0f} hPa {(((r[4] - ph30 + 0.5) % 1) - 0.5) * P0 * 12:+.1f}" for r in rows if r[0] <= 10))

print("\n3. Month of QBO westerly onsets (J..D counts) vs the SAO westerly months")
for p in (10, 7, 5):
    j = list(P).index(p); a = U[:, j] - clim[j][mon]; a = np.convolve(a, np.ones(5) / 5, "same")
    t, s = reversals(T[2:-2], a[2:-2]); mw = np.round((t[s > 0] % 1) * 12).astype(int) % 12; me = np.round((t[s < 0] % 1) * 12).astype(int) % 12
    print(f"   {p:2d} hPa westerly onsets: {' '.join(str(c) for c in np.bincount(mw, minlength=12))}   easterly onsets: {' '.join(str(c) for c in np.bincount(me, minlength=12))}")

fig, axs = plt.subplots(1, 3, figsize=(16, 5.2), facecolor=SURF, gridspec_kw=dict(width_ratios=[1.6, 1, 1]))
ax = axs[0]; m = T >= 1995
im = ax.pcolormesh(T[m], P, U[m].T, cmap="RdBu_r", vmin=-40, vmax=40, shading="auto"); ax.set_yscale("log"); ax.invert_yaxis()
ax.set_yticks(P); ax.set_yticklabels([f"{p:.0f}" for p in P]); ax.set_ylabel("hPa", color=MUTED)
ax.set_title("ERA5 equatorial zonal wind: SAO on top, QBO below", loc="left", fontsize=10.5, color=INK); fig.colorbar(im, ax=ax, pad=0.01, label="m/s")
ax = axs[1]; ax.plot([r[2] for r in rows], P, marker="o", color="#2a78d6", label="R, reversals vs 2.37-yr beat")
ax.plot([np.sqrt(max(0, 1 - 0)) * 0 + (clim[j].var() / U[:, j].var()) for j in range(len(P))], P, marker="s", color="#eb6834", label="calendar (climatology) share")
ax.set_yscale("log"); ax.invert_yaxis(); ax.set_yticks(P); ax.set_yticklabels([f"{p:.0f}" for p in P]); ax.set_xlim(0, 1)
ax.set_title("clock above, beat below", loc="left", fontsize=10.5, color=INK); ax.legend(frameon=False, fontsize=8.5)
ph_up = np.unwrap(2 * np.pi * np.array([r[4] for r in rows])[::-1])[::-1] / (2 * np.pi)    # continuous from 70 hPa upward
rel = (ph_up - ph_up[list(P).index(30)]) * P0 * 12
print("   unwrapped timing relative to 30 hPa (months): " + "  ".join(f"{p:.0f}:{v:+.1f}" for p, v in zip(P, rel)))
ax = axs[2]; lock = np.array([r[2] for r in rows]) > 0.2
ax.plot(rel[lock], P[lock], marker="o", color=INK); ax.scatter(rel[~lock], P[~lock], facecolor="none", edgecolor=INK, s=40, label="no lock (R < 0.2)")
for lag, lab, c in ((-21.8, "trigger, sign as written", "#e87ba4"), (-7.6, "trigger, sign reversed", "#1baf7a")):
    ax.axvline(lag, color=c, lw=1.4, ls="--", label=lab)
ax.set_yscale("log"); ax.invert_yaxis(); ax.set_yticks(P); ax.set_yticklabels([f"{p:.0f}" for p in P])
ax.set_xlabel("reversal timing relative to 30 hPa (months)", color=MUTED); ax.set_title("descent of the reversals", loc="left", fontsize=10.5, color=INK)
ax.legend(frameon=False, fontsize=8.5)
for a_ in axs:
    a_.set_facecolor(SURF); a_.tick_params(colors=MUTED); a_.grid(True, color=GRID, lw=0.5)
fig.tight_layout(); fig.savefig(HERE / "qbo_sao.png", dpi=95, facecolor=SURF)
print("plot: qbo_sao.png")
