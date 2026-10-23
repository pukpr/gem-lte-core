#!/usr/bin/env python3
"""MODEL time-height field of the equatorial zonal wind, 1 to 200 hPa, at 180 deg (2026-10-22).
Structure only, not a fit to observations. Three stacked components, each with a vertical envelope:
  SAO   (1-5 hPa)    calendar-locked semiannual wind with unequal cycles (the 1 hPa ERA5 climatology: Jan easterlies
                     -39, Jul -20, Apr westerlies +32, Oct +23 m/s), decaying below the stratopause
  QBO   (3-70 hPa)   the switched device: comb pulses (solar strobe x lunar ecliptic latitude, qbo_shoulders.py) into a
                     saturating edge just below the SAO, then strobe-gated inward penetration (shoulders), descent
                     ~1 km/month; amplitude envelope peaking near 20 hPa
  n > 0 (100-200 hPa) the ENSO-window standing response at 180 deg (low_wind_emulation.py: ENSO windings on the common
                     manifold), fading above 100 hPa
Heights: z = 7 ln(1000/p) km.
    python3 model_column.py  -> model_column.png
"""
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
_g = {"__file__": str(HERE / "qbo_shoulders.py")}
exec(open(HERE / "qbo_shoulders.py").read().split("B = {")[0], _g)          # pulses, core, warped, TG
pulses, core, warped, TG = _g["pulses"], _g["core"], _g["warped"], _g["TG"]
x = np.array(json.load(open(HERE / "qbo_shoulders.json"))["SATWARP"]["x"])
d0, rho, kappa, _, _, _, _, tau_i, m1, p1, m2, p2 = x
INK, MUTED, SURF = "#0b0b0b", "#52514e", "#fcfcfb"

P = np.geomspace(1, 200, 90); Z = 7 * np.log(1000 / P)
T0, T1 = 1990.0, 2016.0; tg = np.arange(T0, T1, 1 / 24)

# --- SAO: unequal semiannual cycle (amplitude at 1 hPa), decaying downward with an e-folding of ~4 km below 45 km
c = np.array([-39, -30, 5, 32, 15, -10, -20, -10, 10, 23, 5, -25], float)        # monthly climatology at 1 hPa (approx.)
mon = (tg % 1) * 12
sao_t = np.interp(mon, np.arange(13) - 0.5 + 1, np.r_[c, c[0]], period=12)
env_sao = np.exp(-np.clip(48.4 - Z, 0, None) / 4.0)

# --- QBO: device edge just below the SAO (3 hPa, z = 40.7 km), gated penetration downward
tn, a = pulses(d0, rho); ue = core(tn, a / np.abs(a).max(), kappa, tau_i)
z_edge = 7 * np.log(1000 / 3.0); W, ETA, TAU_D = 12.0, 3.0, 12.0                    # km/yr, km^2/yr, yr
qbo = np.zeros((len(P), len(tg)))
for i, z in enumerate(Z):
    depth = z_edge - z
    if depth <= 0.2:
        v = np.interp(tg, TG, ue)
    else:
        v = np.interp(tg, TG, warped(ue, depth, W, ETA, TAU_D, m1, p1, m2, p2))
    qbo[i] = v / (np.std(v) + 1e-9)
env_qbo = 20 * np.exp(-0.5 * ((Z - 26.0) / 6.5) ** 2) * (Z < z_edge + 2)
env_qbo *= np.clip((z_edge + 2 - Z) / 3, 0, 1)
env_qbo *= np.clip((Z - 16.1) / (18.6 - 16.1), 0, 1) ** 0.7                         # fades out between 70 and 100 hPa

# --- n > 0 layer: ENSO windings at 180 deg (fitted to ERA5 150 hPa, low_wind_emulation.py)
_l = {"__file__": str(HERE / "low_wind_emulation.py")}
exec(open(HERE / "low_wind_emulation.py").read().split("fig, axs")[0], _l)
enso = np.interp(tg, _l["te"], _l["em"])
env_n1 = 2.4 * np.clip((Z - 18.4) / (13.2 - 18.4), 0, 1) ** 1.5                    # 0 at 70 hPa (18.6 km) -> 1.6 at 150 hPa
env_n1[Z < 12.0] = 2.4

U = env_sao[:, None] * sao_t[None, :] + env_qbo[:, None] * qbo + env_n1[:, None] * enso[None, :]

fig = plt.figure(figsize=(17, 8.5), facecolor=SURF); gs = fig.add_gridspec(1, 2, width_ratios=[5, 1], wspace=0.04)
ax = fig.add_subplot(gs[0])
im = ax.pcolormesh(tg, P, U, cmap="RdBu_r", vmin=-35, vmax=35, shading="gouraud")
ax.set_yscale("log"); ax.invert_yaxis(); ax.set_yticks([1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 150, 200])
ax.set_yticklabels(["1", "2", "3", "5", "7", "10", "20", "30", "50", "70", "100", "150", "200"])
ax.set_ylabel("pressure (hPa)", color=MUTED); ax.set_xlabel("year", color=MUTED)
for p_, lab in ((2.0, "SAO: the clock (calendar-locked, unequal half-cycles)"), (25, "QBO: switched reversals, gated descent (RFP-like, n = 0)"),
                (150, "below the tropopause: n > 0 standing response (stellarator-like, ENSO windings)")):
    ax.text(T0 + 0.3, p_, lab, fontsize=10, color=INK, va="center", bbox=dict(fc=SURF, ec="none", alpha=0.75, pad=2))
for p_ in (4.0, 85.0):
    ax.axhline(p_, color=INK, lw=0.7, ls=":")
ax.set_title("MODEL equatorial zonal wind at 180°: SAO clock above, switched QBO descending through the stratosphere, n > 0 layer below",
             loc="left", fontsize=11.5, color=INK)
fig.colorbar(im, ax=ax, pad=0.01, label="u (m/s)")
ax2 = fig.add_subplot(gs[1], sharey=ax)
ax2.plot(30 * env_sao, P, color="#eb6834", lw=2, label="SAO")
ax2.plot(env_qbo, P, color="#2a78d6", lw=2, label="QBO")
ax2.plot(env_n1 * np.std(enso), P, color="#1baf7a", lw=2, label="n > 0 (ENSO)")
ax2.set_xlabel("amplitude (m/s)", color=MUTED); ax2.legend(frameon=False, fontsize=9, loc="lower right")
ax2.set_title("component envelopes", loc="left", fontsize=10, color=INK); ax2.tick_params(labelleft=False, colors=MUTED)
for a_ in (ax, ax2):
    a_.set_facecolor(SURF); a_.tick_params(colors=MUTED)
fig.savefig(HERE / "model_column.png", dpi=95, facecolor=SURF, bbox_inches="tight")
print("plot: model_column.png")
