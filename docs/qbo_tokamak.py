#!/usr/bin/env python3
"""QBO as an axisymmetric toroidal device with a reversing poloidal drive (2026-10-22). Derivation: qbo_tokamak.md.

Structure (all cases):
  drive D(t)      the poloidal drive (zonal, m = 0 triggers)
  edge relaxation the toroidal flow at the top edge flips sign when D crosses +-theta (hysteresis, the
                  reversed-field-pinch threshold); theta = 0 is plain sign(D); variant "linear" uses D itself
  penetration     the edge state penetrates downward by advection-diffusion with damping (tokamak current
                  penetration with an inward pinch):  du/dt = -w du/dz + eta d2u/dz2 - u/tau_d,
                  solved with the first-passage Green's function from the edge
  levels          30, 50, 70 hPa at depths d, d + 3.5 km, d + 5.9 km below the edge (log-pressure heights)
Drives:
  SIMPLE  D = cos(2 pi (t - t_phi) / 2.370 yr): the draconic alias, period fixed from the Doodson numbers
  COMB    impulses twice a year (trigger day of year fitted, second half a year later), weight
          sin(dec_moon) * sin(dec_sun) * (solar distance factor), then a leaky integrator (tau_i).
          The lunar declination is computed from Meeus' mean elements, so its phase is astronomical and its
          amplitude carries the 18.6-yr nodal modulation; the sign of the solar declination alternates between
          the two triggers, which moves the alias from 1.73 yr to 2.37 yr (see qbo_tokamak.md).
Fit: differential evolution on the mean correlation over the three levels, 1953-2005; held out 2006-2026.
    python3 qbo_tokamak.py   -> qbo_tokamak.png, qbo_tokamak.json, table on stdout
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
LEGS = HERE.parent / "experiments/Oct2026_generator/legs"
FILES = {30: LEGS / "qbo30_point0E/kit_qbo_u30.dat", 50: LEGS / "qbo50/qbo50.dat", 70: LEGS / "qbo70_point0E/kit_qbo_u70.dat"}
DZ = {30: 0.0, 50: 3.5, 70: 5.9}
DT = 1 / 73.0
TG = np.arange(1900.0, 2027.0, DT)
TR0, TR1, HO0, HO1 = 1953.0, 2006.0, 2006.0, 2026.8
P_DRAC = 365.2422 / 27.212221        # draconic months per year
ALIAS = P_DRAC - round(P_DRAC)        # signed: alias with an annual comb
rad = np.radians


def astro(t):
    T = (t - 2000.0) * 365.25 / 36525.0
    Lp, D, M, Mp, F = (rad(v) for v in (218.3164477 + 481267.88123421 * T, 297.8501921 + 445267.1114034 * T,
                                         357.5291092 + 35999.0502909 * T, 134.9633964 + 477198.8675055 * T,
                                         93.2720950 + 483202.0175233 * T))
    lam = Lp + rad(6.289) * np.sin(Mp) + rad(1.274) * np.sin(2 * D - Mp) + rad(0.658) * np.sin(2 * D) + rad(0.214) * np.sin(2 * Mp) - rad(0.186) * np.sin(M)
    beta = rad(5.128) * np.sin(F) + rad(0.281) * np.sin(Mp + F) + rad(0.278) * np.sin(Mp - F) + rad(0.173) * np.sin(2 * D - F)
    eps = rad(23.439)
    sdm = np.sin(beta) * np.cos(eps) + np.cos(beta) * np.sin(eps) * np.sin(lam)
    lams = rad(280.466 + 36000.770 * T) + rad(1.915) * np.sin(M) + rad(0.020) * np.sin(2 * M)
    return sdm, np.sin(eps) * np.sin(lams), (1 + 0.0167 * np.cos(M)) ** 3, np.sin(beta)


def comb_drive(doy, tau_i, alternate=True, lunar="beta"):
    tn = np.concatenate([np.floor(TG[0]) + k + doy / 365.2422 + np.array([0.0, 0.5]) for k in range(int(TG[-1] - TG[0]) + 1)])
    tn = tn[(tn >= TG[0]) & (tn <= TG[-1])]
    sdm, sds, dist, sb = astro(tn)
    a = (sb if lunar == "beta" else sdm) * (sds if alternate else np.abs(sds)) * dist
    imp = np.zeros_like(TG); np.add.at(imp, np.minimum(((tn - TG[0]) / DT).astype(int), len(TG) - 1), a)
    r = np.exp(-DT / tau_i); D = np.zeros_like(TG); acc = 0.0
    for i, x in enumerate(imp):
        acc = r * acc + x; D[i] = acc
    return D


def relax(D, theta, mode):
    z = (D - D.mean()) / D.std()
    if mode == "linear":
        return z
    u = np.empty_like(z); s = 1.0 if z[0] >= 0 else -1.0
    for i, v in enumerate(z):
        if v > theta: s = 1.0
        elif v < -theta: s = -1.0
        u[i] = s
    return u


def penetrate(ue, z, w, eta, tau_d):
    s = np.arange(1, int(15 / DT)) * DT                     # 15-yr kernel
    g = z / np.sqrt(4 * np.pi * eta * s ** 3) * np.exp(-(z - w * s) ** 2 / (4 * eta * s) - s / tau_d) * DT
    return fftconvolve(ue, np.concatenate([[0.0], g]))[: len(TG)]


DATA = {}
for L, f in FILES.items():
    a = np.loadtxt(f); a = a[np.isfinite(a[:, 1])]; DATA[L] = (a[:, 0], a[:, 1])


def model(case, x, mode):
    if case == "SIMPLE":
        tphi, theta, w, eta, d, tau_d = x; D = np.cos(2 * np.pi * (TG - tphi) / (1 / abs(ALIAS)))
    else:
        doy, theta, w, eta, d, tau_d, tau_i = x; D = comb_drive(doy, tau_i, True, "decl" if case == "COMB_DECL" else "beta")
    ue = relax(D, theta, mode)
    return D, ue, {L: penetrate(ue, d + DZ[L], w, eta, tau_d) for L in FILES}


def scores(out, a, b):
    r = []
    for L, (t, y) in DATA.items():
        m = (t >= a) & (t < b); r.append(np.corrcoef(np.interp(t[m], TG, out[L]), y[m])[0, 1])
    return np.array(r)


BOUNDS = {"SIMPLE": [(0, 1 / abs(ALIAS)), (0, 1.2), (1, 40), (0.1, 60), (0, 12), (0.3, 20)],
          "COMB_DECL": [(0, 365.24), (0, 1.2), (1, 40), (0.1, 60), (0, 12), (0.3, 20), (0.1, 4)],
          "COMB": [(0, 365.24), (0, 1.2), (1, 40), (0.1, 60), (0, 12), (0.3, 20), (0.1, 4)]}
RES = {}
for case in ("SIMPLE", "COMB", "COMB_DECL"):
    for mode in (("threshold", "linear") if case != "COMB_DECL" else ("threshold",)):
        bounds = BOUNDS[case] if mode == "threshold" else [b if i != 1 else (0, 0) for i, b in enumerate(BOUNDS[case])]
        f = lambda x: -scores(model(case, x, mode)[2], TR0, TR1).mean()  # noqa: E731
        opt = differential_evolution(f, bounds, seed=3, maxiter=60, popsize=12, tol=1e-6, polish=False)
        D, ue, out = model(case, opt.x, mode)
        RES[(case, mode)] = dict(x=opt.x.tolist(), train=scores(out, TR0, TR1).tolist(), held=scores(out, HO0, HO1).tolist())
        print(f"{case:7s} {mode:9s} train {np.round(RES[(case, mode)]['train'], 3)} held {np.round(RES[(case, mode)]['held'], 3)}  params {np.round(opt.x, 3)}")


# reversal statistics at 30 hPa: mean period from sign changes; descent lag 30 -> 70 from cross-correlation
def period(t, y):
    s = np.sign(y - 0); c = np.where(np.diff(s) != 0)[0]; return 2 * np.mean(np.diff(t[c])) if len(c) > 2 else np.nan


def lag(t, a, b):
    ls = np.arange(-24, 25); r = [np.corrcoef(a[max(0, -k):len(a) - max(0, k)], b[max(0, k):len(b) - max(0, -k)])[0, 1] for k in ls]
    return ls[int(np.argmax(r))]


tg = np.arange(1953.0, 2024.2, 1 / 12)
obs = {L: np.interp(tg, *DATA[L]) for L in DATA}
print(f"\npredicted draconic alias period {1 / abs(ALIAS):.3f} yr; observed mean period at 30 hPa {period(tg, obs[30]):.3f} yr, 70 hPa {period(tg, obs[70]):.3f} yr")
print(f"observed descent lag 30 -> 70 hPa {lag(tg, obs[30], obs[70])} months")
for k in RES:
    out = model(k[0], np.array(RES[k]["x"]), k[1])[2]; mm = {L: np.interp(tg, TG, out[L]) for L in out}
    RES[k].update(period30=period(tg, mm[30]), lag=int(lag(tg, mm[30], mm[70])))
    print(f"   {k[0]:7s} {k[1]:9s}: model mean period 30 hPa {RES[k]['period30']:.3f} yr, descent lag 30->70 {RES[k]['lag']} months")
if RES[("COMB", "threshold")]:
    doy = RES[("COMB", "threshold")]["x"][0]
    print(f"   COMB trigger days of year: {doy:.0f} and {(doy + 182.62) % 365.24:.0f}  (solstices ~172 / ~355, equinoxes ~80 / ~266)")
json.dump({f"{a}|{b}": v for (a, b), v in RES.items()}, open(HERE / "qbo_tokamak.json", "w"), indent=1)

# spectra of the drives: annual comb with and without the alternating solar sign
x = RES[("COMB", "threshold")]["x"]
INK, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
fig = plt.figure(figsize=(15, 11), facecolor=SURF); gs = fig.add_gridspec(4, 2, height_ratios=[1, 1, 1, 1.1])
ax = fig.add_subplot(gs[0, :]); m = TG >= 1950
fr = np.fft.rfftfreq(m.sum(), DT); w_ = np.hanning(m.sum())
for lab, alt, lun, col in (("lunar ecliptic latitude (draconic), alternating solar sign", True, "beta", "#2a78d6"), ("lunar ecliptic latitude, same sign at both triggers", False, "beta", "#eb6834"),
                          ("lunar declination (tropical month), alternating solar sign", True, "decl", "#1baf7a")):
    Dd = comb_drive(x[0], 0.1, alt, lun)[m]; P = np.abs(np.fft.rfft((Dd - Dd.mean()) * w_)) ** 2
    ax.plot(fr, P / P.max(), color=col, lw=1.4, label=lab)
for f0, lab in ((abs(ALIAS), "2.37 yr"), (1 - abs(ALIAS), "1.73 yr"), (abs(ALIAS) - 1 / 18.613, "−nodal"), (abs(ALIAS) + 1 / 18.613, "+nodal"), (365.2422 / 27.321582 - 13, "2.72 yr (tropical)")):
    ax.axvline(f0, color=MUTED, lw=0.7, ls=":"); ax.text(f0, 1.02, lab, ha="center", fontsize=8, color=MUTED)
ax.set_xlim(0, 1.2); ax.set_xlabel("cycles per year", color=MUTED); ax.set_ylabel("power (norm.)", color=MUTED)
ax.set_title("Semiannual comb x lunar drive: draconic latitude + alternating solar sign -> 2.37 yr; same sign -> 1.73 yr; declination -> 2.72 yr", loc="left", fontsize=11, color=INK)
ax.legend(frameon=False, fontsize=9)
for row, L in enumerate((30, 50, 70), start=1):
    ax = fig.add_subplot(gs[row, :]); t, y = DATA[L]; mm = (t >= 1953)
    ax.plot(t[mm], y[mm] / np.std(y[mm]), color=INK, lw=1.0, label=f"observed {L} hPa")
    for k, col in ((("SIMPLE", "threshold"), "#eda100"), (("COMB", "threshold"), "#2a78d6")):
        out = model(k[0], np.array(RES[k]["x"]), k[1])[2][L]; s = TG >= 1953
        v = out[s]; ax.plot(TG[s], (v - v.mean()) / v.std(), color=col, lw=1.3, alpha=0.85,
                             label=f"{k[0].lower()} (held-out r {RES[k]['held'][[30, 50, 70].index(L)]:+.2f})")
    ax.axvspan(HO0, HO1, color=GRID, alpha=0.6, lw=0); ax.set_xlim(1953, 2027); ax.legend(frameon=False, fontsize=8, ncol=3, loc="upper left")
    ax.set_ylabel(f"{L} hPa (sd)", color=MUTED)
for a in fig.axes:
    a.set_facecolor(SURF); a.tick_params(colors=MUTED); a.grid(True, color=GRID, lw=0.5)
    for sp in ("top", "right"):
        a.spines[sp].set_visible(False)
fig.tight_layout(); fig.savefig(HERE / "qbo_tokamak.png", dpi=95, facecolor=SURF)
print("plot: qbo_tokamak.png")
