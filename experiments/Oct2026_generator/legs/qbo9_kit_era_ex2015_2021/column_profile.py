#!/usr/bin/env python3
"""Classic QBO time-height section, 15 hPa (top) to 70 hPa (bottom): observed vs the nine-leg shared-manifold model
(2026-10-26). Style as docs/model_column.png (RdBu_r, +/-35 m/s, log-pressure axis).

Finest resolution available: a 5-day grid (73 per year), 1979-2026.
  30 / 50 / 70 hPa  ERA5 5-day legs (data and their fitted model); ERA70 ends 2020.45, after which the KIT 70 hPa leg
                    is used (interpolated)
  40 / 20 / 15 hPa  KIT monthly legs, interpolated in time to the 5-day grid (PCHIP, mid-month centred)
Vertically: linear in log-pressure between the six levels onto 80 rows. No smoothing.
Panel (c) rescales each model level to the observed standard deviation (amplitude only; phases and timing untouched):
the regression fit carries 64-86 % of the observed amplitude, which washes out the raw model on a +/-35 m/s scale.
Reads lte_results.<k>.csv (decimal year, model, data) for the legs named in lt.exe.resp.<k>.
    python3 column_profile.py   -> column_profile.png"""
import re
from pathlib import Path

import matplotlib
import numpy as np
from scipy.interpolate import PchipInterpolator

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
INK, MUTED, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
LEVELS = [70, 50, 40, 30, 20, 15]
TG = np.arange(1979.0, 2026.70, 1 / 73.0) + 1 / 146.0           # 5-day bin centres
legs = {}
for k in range(1, 40):
    f = HERE / f"lt.exe.resp.{k}"
    if not f.exists():
        break
    name = f.read_text().split()[1]; L = int(re.search(r"u(\d+)", name).group(1)); src = "ERA" if name.startswith("era") else "KIT"
    a = np.loadtxt(HERE / f"lte_results.{k}.csv", delimiter=","); ok = a[:, 2] != 0
    shift = 1 / 146.0 if src == "ERA" else 1 / 24.0                # bin / month centres
    legs[(src, L)] = (a[ok, 0] + shift, a[ok, 1], a[ok, 2])


def on_grid(src, L, col):
    d, m, y = legs[(src, L)]; v = m if col == "model" else y
    out = np.full(len(TG), np.nan); inside = (TG >= d[0]) & (TG <= d[-1])
    out[inside] = np.interp(TG[inside], d, v) if src == "ERA" else PchipInterpolator(d, v)(TG[inside])
    return out


def column(col):
    rows, source = [], []
    for L in LEVELS:
        if ("ERA", L) in legs:
            v = on_grid("ERA", L, col); kit = on_grid("KIT", L, col); gap = np.isnan(v); v[gap] = kit[gap]
            source.append("E" if not gap.any() else f"E to {legs[('ERA', L)][0][-1]:.1f}, then K")
        else:
            v = on_grid("KIT", L, col); source.append("K")
        rows.append(v)
    return np.array(rows), source


obs, src = column("data"); mod, _ = column("model")
gain = np.array([np.nanstd(o) / np.nanstd(m) for o, m in zip(obs, mod)])
mod_s = (mod - np.nanmean(mod, axis=1, keepdims=True)) * gain[:, None] + np.nanmean(mod, axis=1, keepdims=True)
P = np.geomspace(15, 70, 80); lp = np.log(P); lv = np.log(np.array(LEVELS, float))
order = np.argsort(lv)
interp = lambda U: np.array([np.interp(lp, lv[order], U[order, i]) for i in range(U.shape[1])]).T
Uo, Um, Us = interp(obs), interp(mod), interp(mod_s)

fig = plt.figure(figsize=(17, 17), facecolor=SURF)
gs = fig.add_gridspec(4, 2, height_ratios=[1, 1, 1, 1.05], hspace=0.34, wspace=0.08)


def panel(ax, U, t0, t1, title, cbar=False):
    m = (TG >= t0) & (TG <= t1)
    im = ax.pcolormesh(TG[m], P, U[:, m], cmap="RdBu_r", vmin=-35, vmax=35, shading="gouraud")
    ax.set_yscale("log"); ax.invert_yaxis(); ax.set_yticks(LEVELS); ax.set_yticklabels([str(L) for L in LEVELS])
    ax.minorticks_off(); ax.set_xlim(t0, t1)
    for x in (2015, 2021):
        ax.axvline(x, color=INK, ls="--", lw=1)
    if t0 < 2015 < t1:
        ax.text(2018, 16.2, "excluded from fitting", ha="center", va="top", fontsize=8, color=INK,
                bbox=dict(facecolor="white", alpha=0.7, lw=0, pad=1))
    for x in (2016.15, 2020.0):
        if t0 < x < t1: ax.axvline(x, color="#eda100", ls=":", lw=1.4)
    ax.set_ylabel("pressure (hPa)", color=MUTED); ax.set_title(title, loc="left", fontsize=10.5, color=INK)
    ax.set_facecolor(SURF); ax.tick_params(colors=MUTED)
    return im


ax = fig.add_subplot(gs[0, :]); im = panel(ax, Uo, 1979, 2026.7, "(a) OBSERVED equatorial zonal wind: ERA5 5-day at 30/50/70 hPa, KIT monthly at 40/20/15 hPa (interpolated)")
for L, s in zip(LEVELS, src):
    ax.text(1.002, (np.log(70) - np.log(L)) / (np.log(70) - np.log(15)), s, transform=ax.transAxes, fontsize=7, color=MUTED, va="center")
fig.colorbar(im, ax=ax, pad=0.035, label="u (m/s)")
ax = fig.add_subplot(gs[1, :]); im = panel(ax, Um, 1979, 2026.7, "(b) MODEL: nine legs on one shared manifold (2015–2021 excluded from fitting), same grid and sources")
fig.colorbar(im, ax=ax, pad=0.035, label="u (m/s)")
ax = fig.add_subplot(gs[2, :]); im = panel(ax, Us, 1979, 2026.7, "(c) MODEL, each level rescaled to the observed s.d. (gain " +
                                         ", ".join(f"{L}: {g:.2f}" for L, g in zip(LEVELS, gain)) + "): timing and descent only")
fig.colorbar(im, ax=ax, pad=0.035, label="u (m/s)")
ax = fig.add_subplot(gs[3, 0]); panel(ax, Uo, 2013, 2023.5, "(d) observed, 2013–2023 (orange: the 2016 and 2020 disruptions)")
ax.set_xlabel("year", color=MUTED)
ax = fig.add_subplot(gs[3, 1]); panel(ax, Us, 2013, 2023.5, "(e) model (rescaled as in c), 2013–2023"); ax.set_ylabel(""); ax.set_yticklabels([])
ax.set_xlabel("year", color=MUTED)
fig.savefig(HERE / "column_profile.png", dpi=110, facecolor=SURF, bbox_inches="tight")
print("sources by level:", dict(zip(LEVELS, src))); print("model gain to observed s.d.:", dict(zip(LEVELS, np.round(gain, 2))))
print("wrote column_profile.png")
