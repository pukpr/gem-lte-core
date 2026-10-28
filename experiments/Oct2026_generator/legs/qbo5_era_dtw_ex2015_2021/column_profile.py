#!/usr/bin/env python3
"""All-ERA5 QBO time-height section, 10 hPa (top) to 70 hPa (bottom), observed vs the five-leg shared-manifold model
(2026-10-27). Every row is an ERA5 5-day leg (0N 0E, 00 UTC): 10, 20, 30, 50, 70 hPa; 40 and 15 hPa are not ERA5
levels and are filled by linear interpolation in log-pressure. Style as ../qbo9_kit_era_ex2015_2021/column_profile.py.
Also prints descent lags measured as for the data elsewhere (u low-passed 200 d, adjacent levels, lower behind upper,
parabola-refined) and the same for the high-passed (< 120 d) stripes (5-day sampling: lags below ~5 d unresolved).
    python3 column_profile.py   -> column_profile.png, column_profile.json"""
import json
import re
from pathlib import Path

import matplotlib
import numpy as np
from scipy.signal import butter, filtfilt

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
INK, MUTED, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
legs = {}
for k in range(1, 40):
    f = HERE / f"lt.exe.resp.{k}"
    if not f.exists(): break
    L = int(re.search(r"u(\d+)", f.read_text().split()[1]).group(1)); a = np.loadtxt(HERE / f"lte_results.{k}.csv", delimiter=",")
    legs[L] = a
LV = sorted(legs, reverse=True)                                     # 70 ... 10
key = {L: np.round(legs[L][:, 0] * 73).astype(int) for L in LV}            # align the legs on their shared 5-day bins
common = sorted(set.intersection(*[set(v[legs[L][:, 2] != 0]) for L, v in key.items()]))
rows = {L: np.searchsorted(key[L], common) for L in LV}; t = legs[LV[0]][rows[LV[0]], 0]
OBS = np.array([legs[L][rows[L], 2] for L in LV]); MOD = np.array([legs[L][rows[L], 1] for L in LV])
gain = OBS.std(1) / MOD.std(1); MODS = (MOD - MOD.mean(1, keepdims=True)) * gain[:, None] + MOD.mean(1, keepdims=True)
P = np.geomspace(10, 70, 90); lv = np.log(np.array(LV, float))
grid = lambda U: np.array([np.interp(np.log(P), lv[::-1], U[::-1, i]) for i in range(U.shape[1])]).T
fig = plt.figure(figsize=(17, 17), facecolor=SURF); gs = fig.add_gridspec(4, 2, height_ratios=[1, 1, 1, 1.05], hspace=0.34, wspace=0.08)


def panel(ax, U, t0, t1, title):
    m = (t >= t0) & (t <= t1)
    im = ax.pcolormesh(t[m], P, grid(U)[:, m], cmap="RdBu_r", vmin=-35, vmax=35, shading="gouraud")
    ax.set_yscale("log"); ax.invert_yaxis(); ax.set_yticks([70, 50, 40, 30, 20, 15, 10]); ax.set_yticklabels(["70", "50", "40", "30", "20", "15", "10"]); ax.minorticks_off()
    ax.set_xlim(t0, t1)
    for x in (2015, 2021): ax.axvline(x, color=INK, ls="--", lw=1)
    for x in (2016.15, 2020.0):
        if t0 < x < t1: ax.axvline(x, color="#eda100", ls=":", lw=1.4)
    ax.set_ylabel("pressure (hPa)", color=MUTED); ax.set_title(title, loc="left", fontsize=10.5, color=INK); ax.set_facecolor(SURF); ax.tick_params(colors=MUTED)
    return im


for row, (U, title) in enumerate(((OBS, "(a) OBSERVED equatorial u, ERA5 5-day at 10/20/30/50/70 hPa (0N 0E, 00 UTC)"),
                                  (MOD, "(b) MODEL: five ERA legs on one shared manifold (2015–2021 excluded from fitting)"),
                                  (MODS, "(c) MODEL, each level rescaled to the observed s.d. (gain " + ", ".join(f"{L}: {g:.2f}" for L, g in zip(LV, gain)) + ")"))):
    ax = fig.add_subplot(gs[row, :]); im = panel(ax, U, 1979, t[-1], title); fig.colorbar(im, ax=ax, pad=0.035, label="u (m/s)")
ax = fig.add_subplot(gs[3, 0]); panel(ax, OBS, 2013, 2023.5, "(d) observed, 2013–2023"); ax.set_xlabel("year", color=MUTED)
ax = fig.add_subplot(gs[3, 1]); panel(ax, MODS, 2013, 2023.5, "(e) model (rescaled), 2013–2023"); ax.set_ylabel(""); ax.set_yticklabels([]); ax.set_xlabel("year", color=MUTED)
fig.savefig(HERE / "column_profile.png", dpi=110, facecolor=SURF, bbox_inches="tight")


def lag(x, y, M):
    x = (x - x.mean()) / x.std(); y = (y - y.mean()) / y.std(); n = len(x); Ls = np.arange(-M, M + 1)
    c = np.array([np.mean(x[max(0, -l):n - max(0, l)] * y[max(0, l):n - max(0, -l)]) for l in Ls]); i = int(np.argmax(c))
    d = c[i - 1] - 2 * c[i] + c[i + 1] if 0 < i < len(c) - 1 else 0
    return float((Ls[i] + (0.5 * (c[i - 1] - c[i + 1]) / d if d < 0 else 0)) * 365.25 / 73), float(c[i])


bl, al = butter(3, 1 / 200, fs=73 / 365.25); bh, ah = butter(3, 1 / 120, btype="high", fs=73 / 365.25)
out = {"gain": dict(zip(map(str, LV), gain.round(3).tolist())), "lags": {}}
for nm, f in (("QBO phase (low-pass 200 d)", lambda x: filtfilt(bl, al, x - x.mean())), ("stripes (high-pass 120 d)", lambda x: filtfilt(bh, ah, x - x.mean()))):
    M = 80 if "QBO" in nm else 6
    for src, U in (("obs", OBS), ("model", MOD)):
        pr = [lag(f(U[i + 1]), f(U[i]), M) for i in range(len(LV) - 1)]           # upper (i+1) -> lower (i)
        out["lags"][f"{nm} {src}"] = {f"{LV[i + 1]}->{LV[i]}": v for i, v in enumerate(pr)}
        print(f"{nm:28s} {src:5s} " + "  ".join(f"{LV[i + 1]}->{LV[i]} {l:+6.1f} d (r {c:.2f})" for i, (l, c) in enumerate(pr))
              + f"   cumulative {LV[-1]}->{LV[0]} {sum(l for l, _ in pr):+.0f} d")
json.dump(out, open(HERE / "column_profile.json", "w"), indent=1)
print("gain:", out["gain"]); print("wrote column_profile.png")
