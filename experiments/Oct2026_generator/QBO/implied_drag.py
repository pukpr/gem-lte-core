#!/usr/bin/env python3
"""Test (a): the momentum budget implied by a phase-only (barber-pole) solution, observed vs the nine-leg LTE model
(2026-10-27). The LTE fit gives u(z,t) as a phase pattern; the zonal momentum equation turns any u(z,t) into the wave
drag that must be deposited, G = du/dt + u/tau, and the upward flux that must arrive from below,
F(z) = F(top) + integral_z^top rho G dz' (F(top) = 0: whatever leaves above 15 hPa is not counted).
Data: the six KIT monthly legs (70, 50, 40, 30, 20, 15 hPa) of ../legs/qbo9_kit_era70full_ex2015_2021, observed and
model columns on the same months. z = -H ln(p/1000), H = 7 km; rho = p/(R T) with a tropical T profile.
A pure barber pole gives G = c du/dz (drag tracks shear), so that correlation alone only restates the descent.
The physical content is in what a translating pattern does NOT fix:
  - the magnitude of the required flux (mPa) and whether the model carries it,
  - the westerly/easterly asymmetry: observed westerly shear zones descend faster and more regularly (narrow Kelvin
    absorption) than easterly ones,
  - G of the model vs G observed, level by level, training years vs the held-out 2015-2021.
tau = infinity (mean-flow damping neglected) and tau = 1 yr are both reported.
    python3 implied_drag.py   -> implied_drag.png, implied_drag.json"""
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
RUN = HERE.parent / "legs" / "qbo9_kit_era70full_ex2015_2021"
INK, MUTED, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
LV = np.array([70, 50, 40, 30, 20, 15.0]); TK = np.array([197, 205, 210, 216, 221, 224.0])
H, R, DAY = 7000.0, 287.0, 86400.0
Z = -H * np.log(LV / 1000.0); RHO = LV * 100 / (R * TK)

cols = []
for k in range(1, 7):
    a = np.loadtxt(RUN / f"lte_results.{k}.csv", delimiter=","); cols.append(a)
t = cols[0][:, 0]; assert all(np.allclose(c[:, 0], t) for c in cols)
ok = np.all([c[:, 2] != 0 for c in cols], axis=0)
OBS = np.array([c[:, 2] for c in cols]); MOD = np.array([c[:, 1] for c in cols])
t, OBS, MOD = t[ok], OBS[:, ok], MOD[:, ok]
dt = np.median(np.diff(t)) * 365.25 * DAY
HELD = (t >= 2015) & (t <= 2021)


def budget(U, tau_yr=None):
    G = np.gradient(U, axis=1) / dt                                   # m/s^2
    if tau_yr: G = G + U / (tau_yr * 365.25 * DAY)
    S = np.gradient(U, Z, axis=0)                                      # du/dz, 1/s
    F = np.zeros_like(U)                                               # Pa (upward flux needed below each level)
    for i in range(len(LV) - 2, -1, -1):
        F[i] = F[i + 1] + 0.5 * (RHO[i] * G[i] + RHO[i + 1] * G[i + 1]) * (Z[i + 1] - Z[i])
    return G, S, F


res = {}
for tau in (None, 1.0):
    Go, So, Fo = budget(OBS, tau); Gm, Sm, Fm = budget(MOD, tau); key = "tau_inf" if tau is None else "tau_1yr"
    r = {}
    for i, L in enumerate(LV):
        cc = lambda a, b, m: float(np.corrcoef(a[m], b[m])[0, 1])
        r[int(L)] = {"G_model_vs_obs_train": cc(Gm[i], Go[i], ~HELD), "G_model_vs_obs_held": cc(Gm[i], Go[i], HELD),
                     "G_vs_shear_obs": cc(Go[i], So[i], np.ones_like(HELD)), "G_vs_shear_model": cc(Gm[i], Sm[i], np.ones_like(HELD)),
                     "rms_G_obs_m_s_day": float(np.std(Go[i]) * DAY), "rms_G_model_m_s_day": float(np.std(Gm[i]) * DAY)}
    # westerly / easterly asymmetry at interior levels: accelerations in westerly vs easterly shear
    asym = {}
    for nm, G, S in (("obs", Go, So), ("model", Gm, Sm)):
        g, s = G[1:-1].ravel() * DAY, S[1:-1].ravel()
        asym[nm] = {"mean_G_in_westerly_shear": float(g[s > 0].mean()), "mean_G_in_easterly_shear": float(g[s < 0].mean()),
                    "p95_westerly_accel": float(np.percentile(g, 95)), "p05_easterly_accel": float(np.percentile(g, 5))}
        asym[nm]["ratio_west_to_east"] = -asym[nm]["p95_westerly_accel"] / asym[nm]["p05_easterly_accel"]
    res[key] = {"levels": r, "asymmetry": asym,
                "F70_obs_rms_mPa": float(np.std(Fo[0]) * 1e3), "F70_model_rms_mPa": float(np.std(Fm[0]) * 1e3),
                "F70_obs_p95_mPa": float(np.percentile(Fo[0], 95) * 1e3), "F70_obs_p05_mPa": float(np.percentile(Fo[0], 5) * 1e3),
                "F70_model_p95_mPa": float(np.percentile(Fm[0], 95) * 1e3), "F70_model_p05_mPa": float(np.percentile(Fm[0], 5) * 1e3),
                "F70_model_vs_obs_train": float(np.corrcoef(Fm[0][~HELD], Fo[0][~HELD])[0, 1]),
                "F70_model_vs_obs_held": float(np.corrcoef(Fm[0][HELD], Fo[0][HELD])[0, 1])}
    if tau is None:
        keep = (Go, Gm, Fo, Fm)
json.dump(res, open(HERE / "implied_drag.json", "w"), indent=1)

Go, Gm, Fo, Fm = keep
fig = plt.figure(figsize=(17, 13), facecolor=SURF); gs = fig.add_gridspec(3, 2, hspace=0.38, wspace=0.12, width_ratios=[1.6, 1])
P = np.geomspace(15, 70, 60)
grid = lambda A: np.array([np.interp(np.log(P), np.log(LV[::-1]), A[::-1, j]) for j in range(A.shape[1])]).T
m = (t >= 2000) & (t <= 2024)
for row, (A, title) in enumerate(((Go, "(a) OBSERVED implied drag G = du/dt (KIT monthly), m/s per day"),
                                  (Gm, "(b) MODEL implied drag (nine-leg LTE, KIT legs), m/s per day"))):
    ax = fig.add_subplot(gs[row, 0]); v = 0.6
    im = ax.pcolormesh(t[m], P, grid(A * DAY)[:, m], cmap="PuOr_r", vmin=-v, vmax=v, shading="gouraud")
    ax.contour(t[m], P, grid(OBS if row == 0 else MOD)[:, m], levels=[0], colors=INK, linewidths=0.8)
    ax.set_yscale("log"); ax.invert_yaxis(); ax.set_yticks(LV); ax.set_yticklabels([f"{L:g}" for L in LV]); ax.minorticks_off()
    for x in (2015, 2021): ax.axvline(x, color=INK, ls="--", lw=1)
    ax.set_title(title + " (black: u = 0)", loc="left", fontsize=10.5, color=INK); ax.set_ylabel("pressure (hPa)", color=MUTED)
    ax.set_facecolor(SURF); ax.tick_params(colors=MUTED); fig.colorbar(im, ax=ax, pad=0.015)
ax = fig.add_subplot(gs[2, 0])
ax.plot(t[m], Fo[0][m] * 1e3, color=INK, lw=1.2, label="observed"); ax.plot(t[m], Fm[0][m] * 1e3, color="#d94a4a", lw=1.2, label="model")
for x in (2015, 2021): ax.axvline(x, color=INK, ls="--", lw=1)
ax.axhline(0, color=MUTED, lw=0.6)
rr = res["tau_inf"]
ax.set_title(f"(c) net upward flux required at 70 hPa to drive 70–15 hPa (mPa); r model vs obs: train {rr['F70_model_vs_obs_train']:+.2f}, held-out {rr['F70_model_vs_obs_held']:+.2f}",
             loc="left", fontsize=10.5, color=INK)
ax.set_ylabel("mPa", color=MUTED); ax.legend(fontsize=9, frameon=False, ncol=2); ax.set_facecolor(SURF); ax.tick_params(colors=MUTED)
ax = fig.add_subplot(gs[0, 1])
x = np.arange(len(LV))
for off, key, col, lab in ((-0.2, "G_model_vs_obs_train", "#2a78d6", "train"), (0.2, "G_model_vs_obs_held", "#eda100", "held-out 2015–21")):
    ax.bar(x + off, [rr["levels"][int(L)][key] for L in LV], 0.38, color=col, label=lab)
ax.set_xticks(x); ax.set_xticklabels([f"{L:g}" for L in LV]); ax.set_xlabel("hPa", color=MUTED); ax.set_ylim(0, 1)
ax.set_title("(d) r(G model, G observed) by level", loc="left", fontsize=10.5, color=INK); ax.legend(fontsize=9, frameon=False)
ax.set_facecolor(SURF); ax.tick_params(colors=MUTED)
ax = fig.add_subplot(gs[1, 1])
ax.bar(x - 0.2, [rr["levels"][int(L)]["rms_G_obs_m_s_day"] for L in LV], 0.38, color=INK, label="observed")
ax.bar(x + 0.2, [rr["levels"][int(L)]["rms_G_model_m_s_day"] for L in LV], 0.38, color="#d94a4a", label="model")
ax.set_xticks(x); ax.set_xticklabels([f"{L:g}" for L in LV]); ax.set_xlabel("hPa", color=MUTED)
ax.set_title("(e) rms implied drag (m/s per day)", loc="left", fontsize=10.5, color=INK); ax.legend(fontsize=9, frameon=False)
ax.set_facecolor(SURF); ax.tick_params(colors=MUTED)
ax = fig.add_subplot(gs[2, 1])
for nm, col, off in (("obs", INK, -0.2), ("model", "#d94a4a", 0.2)):
    A = rr["asymmetry"][nm]; ax.bar(np.array([0, 1]) + off, [A["p95_westerly_accel"], -A["p05_easterly_accel"]], 0.38, color=col, label=f"{nm} (W/E {A['ratio_west_to_east']:.2f})")
ax.set_xticks([0, 1]); ax.set_xticklabels(["westerly accel. (95th pct)", "easterly accel. (−5th pct)"])
ax.set_title("(f) asymmetry, 50–20 hPa (m/s per day)", loc="left", fontsize=10.5, color=INK); ax.legend(fontsize=9, frameon=False)
ax.set_facecolor(SURF); ax.tick_params(colors=MUTED)
fig.savefig(HERE / "implied_drag.png", dpi=105, facecolor=SURF, bbox_inches="tight")

for key, rr in res.items():
    print(f"== {key}")
    for L, v in rr["levels"].items():
        print(f"  {L:3d} hPa  r(Gm,Go) train {v['G_model_vs_obs_train']:+.2f} held {v['G_model_vs_obs_held']:+.2f} | r(G,shear) obs {v['G_vs_shear_obs']:+.2f} model {v['G_vs_shear_model']:+.2f} | rms G obs {v['rms_G_obs_m_s_day']:.3f} model {v['rms_G_model_m_s_day']:.3f} m/s/d")
    for nm, A in rr["asymmetry"].items():
        print(f"  {nm:5s} westerly p95 {A['p95_westerly_accel']:+.3f}, easterly p05 {A['p05_easterly_accel']:+.3f} m/s/d, W/E ratio {A['ratio_west_to_east']:.2f}")
    print(f"  F(70 hPa) rms obs {rr['F70_obs_rms_mPa']:.2f} mPa (p05 {rr['F70_obs_p05_mPa']:+.2f}, p95 {rr['F70_obs_p95_mPa']:+.2f}); model {rr['F70_model_rms_mPa']:.2f} "
          f"(p05 {rr['F70_model_p05_mPa']:+.2f}, p95 {rr['F70_model_p95_mPa']:+.2f}); r train {rr['F70_model_vs_obs_train']:+.2f}, held {rr['F70_model_vs_obs_held']:+.2f}")
print("wrote implied_drag.png")
