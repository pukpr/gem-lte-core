#!/usr/bin/env python3
"""QBO device with a saturating flux core: staircase reversals with semiannual shoulders (2026-10-22).

Extends qbo_tokamak.py (same astronomy, Green's-function penetration, data). The edge is no longer a +-1 switch:
  core flux   Phi_n = clip(lambda * Phi_{n-1} + kappa * w_n * a_n, -1, +1), held between pulses (volt-second steps
              on a saturable transformer core), lambda = exp(-dt / tau_i)
  pulses      a_n = sin(beta_moon) * sin(dec_sun) * (solar distance)^3 at two triggers a year (day d0 and d0 + 182.6);
              the second ("weak") trigger is scaled by rho, w_n = 1 or rho
  satellites  optional: each strobe also fires two satellite pulses at +-delta days with relative weight sigma
              (sub-banding of the stroboscopic impulse)
Models compared (each refit with the same objective):
  THRESH      the earlier comb + hysteresis switch (no shoulders possible)
  SAT         saturating core, unequal triggers
  SAT_SAT     SAT + satellite pulses
Objective (1953-2005): mean over 30/50/70 hPa of [corr(u) + corr(3-month change of u)] / 2, so the timing and size of
the steps count, not only the sign. Held out 2006-2026. Also: correlation of the model's calendar-month profile of
|du/dt| with the observed one, and the share of "shoulder" months.
    python3 qbo_shoulders.py  -> qbo_shoulders.png, qbo_shoulders.json
"""
import json
from pathlib import Path

import matplotlib
import numpy as np
from scipy.optimize import differential_evolution
from scipy.signal import fftconvolve

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
__file__ = str(HERE / "qbo_tokamak.py")
exec(open(HERE / "qbo_tokamak.py").read().split("BOUNDS = ")[0])   # astro, comb_drive, relax, penetrate, DATA, TG, DT, DZ
INK, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
TR0, TR1, HO0, HO1 = 1953.0, 2006.0, 2006.0, 2026.8
LV = (30, 50, 70)


def pulses(d0, rho, sigma=0.0, delta=0.0):
    base = np.concatenate([np.floor(TG[0]) + k + d0 / 365.2422 + np.array([0.0, 0.5]) for k in range(int(TG[-1] - TG[0]) + 1)])
    wt = np.tile([1.0, rho], len(base) // 2)
    tn, wn = [base], [wt]
    if sigma > 0:
        for s in (-1, 1):
            tn.append(base + s * delta / 365.2422); wn.append(wt * sigma)
    tn, wn = np.concatenate(tn), np.concatenate(wn)
    k = (tn >= TG[0]) & (tn <= TG[-1]); tn, wn = tn[k], wn[k]
    _, sds, dist, sb = astro(tn)
    o = np.argsort(tn)
    return tn[o], (sb * sds * dist * wn)[o]


def core(tn, a, kappa, tau_i):
    idx = np.minimum(((tn - TG[0]) / DT).astype(int), len(TG) - 1)
    lam = np.exp(-DT / tau_i); phi = 0.0; out = np.empty_like(TG); j = 0
    for i in range(len(TG)):
        phi *= lam
        while j < len(idx) and idx[j] == i:
            phi = min(1.0, max(-1.0, phi + kappa * a[j])); j += 1
        out[i] = phi
    return out


def warped(ue, z, w, eta, tau_d, m1, p1, m2, p2):
    """penetration whose rate (pinch and diffusion together) is gated through the year:
    rate(t) = exp(m1 cos 2pi(t - p1) + m2 cos 4pi(t - p2)), normalized to mean 1; solved in the warped time
    tau = integral of rate, where the transport is uniform (damping stays in real time, approximately)"""
    r = np.exp(m1 * np.cos(2 * np.pi * (TG - p1)) + m2 * np.cos(4 * np.pi * (TG - p2))); r /= r.mean()
    tau = np.concatenate([[0.0], np.cumsum(r[:-1]) * DT]) + TG[0]
    grid = np.arange(tau[0], tau[-1], DT)
    u_tau = np.interp(grid, tau, ue)
    s = np.arange(1, int(15 / DT)) * DT
    g = z / np.sqrt(4 * np.pi * eta * s ** 3) * np.exp(-(z - w * s) ** 2 / (4 * eta * s) - s / tau_d) * DT
    v = fftconvolve(u_tau, np.concatenate([[0.0], g]))[: len(grid)]
    return np.interp(tau, grid, v)


def run(name, x):
    if name == "THRESH":
        d0, theta, w, eta, dd, tau_d, tau_i = x
        ue = relax(comb_drive(d0, tau_i), theta, "threshold")
    elif name in ("WARP", "SATWARP"):
        if name == "WARP":
            d0, theta, w, eta, dd, tau_d, tau_i, m1, p1, m2, p2 = x; ue = relax(comb_drive(d0, tau_i), theta, "threshold")
        else:
            d0, rho, kappa, w, eta, dd, tau_d, tau_i, m1, p1, m2, p2 = x
            tn, a = pulses(d0, rho); ue = core(tn, a / np.abs(a).max(), kappa, tau_i)
        return ue, {L: warped(ue, dd + DZ[L], w, eta, tau_d, m1, p1, m2, p2) for L in LV}
    else:
        d0, rho, kappa, w, eta, dd, tau_d, tau_i = x[:8]
        sig, dl = (x[8], x[9]) if name == "SAT_SAT" else (0.0, 0.0)
        tn, a = pulses(d0, rho, sig, dl); ue = core(tn, a / np.abs(a).max(), kappa, tau_i)
    return ue, {L: penetrate(ue, dd + DZ[L], w, eta, tau_d) for L in LV}


def metrics(out, a, b):
    r, rd = [], []
    for L in LV:
        t, y = DATA[L]; m = (t >= a) & (t < b); mv = np.interp(t[m], TG, out[L]); yy = y[m]
        r.append(np.corrcoef(mv, yy)[0, 1]); rd.append(np.corrcoef(mv[3:] - mv[:-3], yy[3:] - yy[:-3])[0, 1])
    return np.array(r), np.array(rd)


def season(t, y):
    d = np.abs(np.diff(y)); mon = np.round((t[1:] % 1) * 12).astype(int) % 12
    return np.array([d[mon == k].mean() for k in range(12)]) / d.mean()


def shoulder_share(y):
    z = (y - np.median(y)) / np.percentile(np.abs(y - np.median(y)), 90); d = np.abs(np.diff(z))
    return np.mean((np.abs(z[1:]) < 0.4) & (d < np.percentile(d, 40)))


B = {"THRESH": [(0, 365.24), (0, 1.2), (1, 40), (0.05, 60), (0, 12), (0.3, 20), (0.1, 4)],
     "SAT": [(0, 365.24), (0, 1), (0.05, 3), (1, 40), (0.05, 60), (0, 12), (0.3, 20), (0.2, 20)]}
B["SAT_SAT"] = B["SAT"] + [(0, 1), (5, 60)]
GATE = [(0, 3), (0, 1), (0, 3), (0, 0.5)]
B["WARP"] = B["THRESH"] + GATE
B["SATWARP"] = B["SAT"] + GATE
RES = {}
for name in ("THRESH", "SAT", "SAT_SAT", "WARP", "SATWARP"):
    def f(x):
        r, rd = metrics(run(name, x)[1], TR0, TR1); return -0.5 * (r.mean() + rd.mean())
    opt = differential_evolution(f, B[name], seed=5, maxiter=70, popsize=12, tol=1e-7, polish=False)
    ue, out = run(name, opt.x)
    r, rd = metrics(out, TR0, TR1); h, hd = metrics(out, HO0, HO1)
    tg = np.arange(1953.0, 2024.2, 1 / 12)
    sc = [np.corrcoef(season(tg, np.interp(tg, TG, out[L])), season(tg, np.interp(tg, *DATA[L])))[0, 1] for L in LV]
    sh = [shoulder_share(np.interp(tg, TG, out[L])) for L in LV]
    RES[name] = dict(x=opt.x.tolist(), train=r.tolist(), train_d=rd.tolist(), held=h.tolist(), held_d=hd.tolist(), season_r=sc, shoulders=sh)
    print(f"{name:8s} train r {np.round(r, 3)} d3 {np.round(rd, 3)} | held r {np.round(h, 3)} d3 {np.round(hd, 3)} | "
          f"seasonal |du/dt| profile r {np.round(sc, 2)} | shoulder months {np.round(np.array(sh) * 100, 1)}%")
    print(f"         params {np.round(opt.x, 3)}")
obs_sh = [shoulder_share(np.interp(np.arange(1953.0, 2024.2, 1 / 12), *DATA[L])) for L in LV]
print(f"observed shoulder months {np.round(np.array(obs_sh) * 100, 1)}%")
json.dump(RES, open(HERE / "qbo_shoulders.json", "w"), indent=1)

# figure: 1985-2005 at 30 and 50 hPa, observed vs THRESH vs SAT(_SAT); plus seasonal profiles
best = max(("SAT", "SAT_SAT", "WARP", "SATWARP"), key=lambda k: np.mean(RES[k]["train"]) + np.mean(RES[k]["train_d"]))
fig = plt.figure(figsize=(16, 10), facecolor=SURF); gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 0.9])
for row, L in enumerate((30, 50)):
    ax = fig.add_subplot(gs[row, :]); t, y = DATA[L]; m = (t >= 1985) & (t < 2005)
    z = lambda v: (v - v.mean()) / v.std()  # noqa: E731
    ax.plot(t[m], z(y[m]), color=INK, lw=1.4, label=f"observed {L} hPa")
    for name, col in (("THRESH", "#eda100"), (best, "#2a78d6")):
        out = run(name, np.array(RES[name]["x"]))[1][L]; s = (TG >= 1985) & (TG < 2005)
        ax.plot(TG[s], z(out[s]), color=col, lw=1.4, alpha=0.9, label=f"{name.lower()} (train r {RES[name]['train'][LV.index(L)]:.2f})")
    xb = RES[best]["x"]
    if best in ("WARP", "SATWARP"):
        m1, p1, m2, p2 = xb[-4:]; tt = np.arange(1985, 2005, 1 / 73)
        rr = np.exp(m1 * np.cos(2 * np.pi * (tt - p1)) + m2 * np.cos(4 * np.pi * (tt - p2))); rr /= rr.mean()
        ax.fill_between(tt, -2.6, -2.6 + 0.35 * rr, color="#1baf7a", alpha=0.35, lw=0, label="penetration rate (gate)")
    ax.set_xlim(1985, 2005); ax.legend(frameon=False, fontsize=8.5, ncol=3, loc="upper left"); ax.set_ylabel("sd units", color=MUTED)
    ax.set_title(f"{L} hPa, 1985-2005: reversals with shoulders (green band at the bottom: gated penetration rate)", loc="left", fontsize=10.5, color=INK)
tg = np.arange(1953.0, 2024.2, 1 / 12)
for k, L in enumerate(LV):
    ax = fig.add_subplot(gs[2, k]); mo = np.arange(1, 13)
    ax.plot(mo, season(tg, np.interp(tg, *DATA[L])), color=INK, marker="o", lw=1.5, label="observed")
    for name, col in (("THRESH", "#eda100"), (best, "#2a78d6")):
        ax.plot(mo, season(tg, np.interp(tg, TG, run(name, np.array(RES[name]["x"]))[1][L])), color=col, marker="o", lw=1.3, label=name.lower())
    ax.set_xticks(mo); ax.set_xticklabels("JFMAMJJASOND"); ax.set_title(f"{L} hPa: |du/dt| by calendar month", loc="left", fontsize=10, color=INK)
    if k == 0: ax.legend(frameon=False, fontsize=8)
for a in fig.axes:
    a.set_facecolor(SURF); a.tick_params(colors=MUTED); a.grid(True, color=GRID, lw=0.5)
    for sp in ("top", "right"):
        a.spines[sp].set_visible(False)
fig.tight_layout(); fig.savefig(HERE / "qbo_shoulders.png", dpi=95, facecolor=SURF)
print("plot: qbo_shoulders.png  (best saturating model:", best, ")")
