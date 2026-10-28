"""QBO from an equinox-impulsed draconic forcing, fitted ONLY to the amplitude spectrum; time domain shown, not fitted.
(2026-10-26)  Pipeline (standalone version of Impulse_Delta -> IIR -> sin(k L) used in QBO_FORMULATION.md):
  beta(t)  Moon's ecliptic-latitude factor sin F (Meeus mean argument of latitude; real draconic phase)
  comb     impulses at the actual equinox instants (Meeus), weight 1 at March (S->N solar crossing), rho at September
  L(t)     leaky integrator of the impulses beta(t_n) x weight, time constant Ti (Ti >> 1 yr = sample-and-hold)
  u(t)     A sin(k L + d), then a one-pole response low-pass (tau), sampled at mid-month
Fit      parameters (rho, Ti, k, d, tau) by differential evolution on the SHAPE of the smoothed log10 amplitude spectrum
         (0.1-3 cpy, +/-0.015 cpy running mean, Hann window) of the deseasonalised observed series; A then from the mean
         log offset. Nothing from the time domain is used: u and -u have the same spectrum, so both are drawn.
Data     KIT/FU Berlin radiosonde monthly zonal wind (Canton/Gan/Singapore), default 30 hPa, 1953-2026.
Output   qbo_spectral_fit_<level>_<crit>.png / .json
         Criterion "log" (default): shape of log10 amplitude, every bin equal; "lin": linear amplitude (energy-weighted,
         dominated by the QBO peak and its sidebands), scale-free (model rescaled by LS before comparing).
Response "sat" (5th argument): y = tanh(G (sin(kL+d) + b)) / tanh(G) before the low-pass, adding a saturation gain G
         (G -> 0: the linear model; large G: square wave with plateaus) and a bias b (unequal westerly/easterly phases).
usage: python3 qbo_spectral_fit.py [level hPa, default 30] [log|lin] [k max, default 60] [seed, default 7] [lin|sat]"""
import datetime as dt, json, os, re, sys
import numpy as np
from scipy.optimize import differential_evolution
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
LEVEL = int(sys.argv[1]) if len(sys.argv) > 1 else 30
CRIT = sys.argv[2] if len(sys.argv) > 2 else "log"
KMAX = float(sys.argv[3]) if len(sys.argv) > 3 else 60.0
SEED = int(sys.argv[4]) if len(sys.argv) > 4 else 7
RESP = sys.argv[5] if len(sys.argv) > 5 else "lin"
TAG = f"{LEVEL}_{CRIT}" + ("" if KMAX == 60 else f"_k{KMAX:g}") + ("" if SEED == 7 else f"_s{SEED}") + ("" if RESP == "lin" else f"_{RESP}")
SRC = os.path.join(HERE, "../../experiments/Oct2026_generator/QBO/kit_qbo.dat")
LEVELS = [70, 50, 40, 30, 20, 15, 10]; ENDS = [16, 23, 30, 37, 44, 51, 58]

# ---- observed series (deseasonalised: each calendar month's mean removed) -----------------------------------------
rows = {}
for ln in open(SRC).read().splitlines()[9:]:
    if len(ln) < 20 or not ln[:5].strip().isdigit():
        continue
    yy, mo = int(ln[6:8]), int(ln[8:10]); yr = (1900 if yy >= 50 else 2000) + yy
    v = {m.end() + 10: int(m.group()) for m in re.finditer(r"-?\d+", ln[10:])}; rows[(yr, mo)] = v.get(ENDS[LEVELS.index(LEVEL)], None)
keys = [k for k in sorted(rows) if rows[k] is not None]
t = np.array([y + (m - 0.5) / 12 for y, m in keys]); month = np.array([m for _, m in keys]); u_obs = np.array([rows[k] / 10 for k in keys], float)
assert np.allclose(np.diff(t), 1 / 12, atol=1e-6), "gap in the monthly record"
clim = np.array([u_obs[month == m].mean() for m in range(1, 13)]); u_obs = u_obs - clim[month - 1]
print(f"KIT {LEVEL} hPa: {len(t)} months {t[0]:.2f}-{t[-1]:.2f}, sd {u_obs.std():.1f} m/s")
JDm = np.array([dt.datetime(y, m, 15).toordinal() + 1721424.5 for y, m in keys])           # mid-month Julian day


def equinoxes(y):
    Y = (y - 2000) / 1000
    return (2451623.80984 + 365242.37404 * Y + 0.05169 * Y ** 2 - 0.00411 * Y ** 3 - 0.00057 * Y ** 4,
            2451810.21715 + 365242.01767 * Y - 0.11575 * Y ** 2 + 0.00337 * Y ** 3 + 0.00078 * Y ** 4)


EQ = np.array([(e, k) for y in range(1900, 2028) for k, e in enumerate(equinoxes(y))])     # long spin-up before 1953
BETA = np.sin(np.radians(93.2720950 + 483202.0175233 * (EQ[:, 0] - 2451545.0) / 36525))


def model(p):
    rho, Ti, k, d, tau = p[:5]
    w = np.where(EQ[:, 1] == 0, 1.0, rho) * BETA
    L_n = np.empty(len(EQ)); acc = 0.0
    for i in range(len(EQ)):                                     # leaky integrator, updated at each impulse
        if i: acc *= np.exp(-(EQ[i, 0] - EQ[i - 1, 0]) / (Ti * 365.25))
        acc += w[i]; L_n[i] = acc
    idx = np.searchsorted(EQ[:, 0], JDm, side="right") - 1
    L = L_n[idx] * np.exp(-(JDm - EQ[idx, 0]) / (Ti * 365.25))
    x = np.sin(k * L + d)
    if RESP == "sat":
        G, b = p[5], p[6]; x = np.tanh(G * (x + b)) / np.tanh(G)
    if tau > 0.05:                                               # one-pole response, tau in months
        a = np.exp(-1 / tau); y = np.empty_like(x); y[0] = x[0]
        for i in range(1, len(x)): y[i] = a * y[i - 1] + (1 - a) * x[i]
        x = y
    return x - x.mean(), L


N = len(t); NF = 8 * N; f = np.fft.rfftfreq(NF, 1 / 12); W = np.hanning(N)
band = (f >= 0.1) & (f <= 3.0); sm = max(1, int(round(0.015 / (f[1] - f[0]))))


def logspec(x):
    a = np.abs(np.fft.rfft((x - x.mean()) * W, NF)); a = np.convolve(a, np.ones(2 * sm + 1) / (2 * sm + 1), mode="same")
    return np.log10(a[band] + 1e-12)


S_obs = logspec(u_obs)


A_obs = 10 ** S_obs


def cost(p):
    if CRIT == "lin":
        m_ = 10 ** logspec(model(p)[0]); s_ = np.dot(m_, A_obs) / np.dot(m_, m_)
        return float(np.sqrt(np.mean((s_ * m_ - A_obs) ** 2)) / np.sqrt(np.mean(A_obs ** 2)))   # relative rms
    d_ = logspec(model(p)[0]) - S_obs; return float(np.sqrt(np.mean((d_ - d_.mean()) ** 2)))


bounds = [(-1, 1), (0.3, 20), (0.2, KMAX), (0, 2 * np.pi), (0, 6)] + ([(0.05, 20), (-0.8, 0.8)] if RESP == "sat" else [])
fit = differential_evolution(cost, bounds, seed=SEED, maxiter=300, popsize=25, tol=1e-7, polish=True)
p = fit.x; xm, L = model(p)
mA = 10 ** logspec(xm); A = float(np.dot(mA, A_obs) / np.dot(mA, mA)) if CRIT == "lin" else 10 ** float(np.mean(S_obs - logspec(xm))); um = A * xm
rho, Ti, k, d, tau = p[:5]
satx = f"  G {p[5]:.2f}  b {p[6]:+.2f}" if RESP == "sat" else ""
print(f"spectral fit{satx}: rho {rho:+.2f}  Ti {Ti:.2f} yr  k {k:.2f}  d {np.degrees(d):.0f} deg  tau {tau:.2f} mo  A {A:.1f} m/s;"
      f" criterion {CRIT}: misfit {fit.fun:.3f}")
q = (365.242199 / 27.212221) % 1
def fpk(x):
    m = (f[band] > 0.3) & (f[band] < 0.6); return float(f[band][m][np.argmax(logspec(x)[m])])
print(f"main peak: observed {fpk(u_obs):.3f} cpy, model {fpk(um):.3f} cpy (draconic alias {q:.4f})")
json.dump({"level_hPa": LEVEL, "rho": rho, "Ti_yr": Ti, "k": k, "d_deg": float(np.degrees(d)), "tau_months": tau, "A_ms": A, "response": RESP, **({"G": float(p[5]), "b": float(p[6])} if RESP == "sat" else {}),
           "criterion": CRIT, "misfit": fit.fun, "note": "fitted to amplitude spectrum only; time domain not used"},
          open(os.path.join(HERE, f"qbo_spectral_fit_{TAG}.json"), "w"), indent=1)

fig = plt.figure(figsize=(15, 13)); gs = fig.add_gridspec(4, 1, height_ratios=[1.1, 1, 1, 0.7])
a0 = fig.add_subplot(gs[0]); fb = f[band]
a0.plot(fb, 10 ** S_obs, color="k", lw=1.2, label=f"observed KIT {LEVEL} hPa (deseasonalised)")
a0.plot(fb, A * 10 ** logspec(xm), color="C3", lw=1, label="model (spectral fit)")
for x0, lab in ((q, "q"), (1 - q, "1-q"), (1 + q, "1+q"), (2 - q, "2-q"), (2 + q, "2+q"), (2 * q, "2q"), (3 * q, "3q")):
    a0.axvline(x0, color="0.6", lw=0.6, ls=":"); a0.text(x0, 1.0, lab, transform=a0.get_xaxis_transform(), fontsize=7, ha="center", va="bottom", color="0.4")
a0.set_yscale("log"); a0.set_xlabel("frequency (cycles/yr)"); a0.set_ylabel("amplitude (smoothed)"); a0.legend(fontsize=8)
a0.set_title(f"FIT STEP (only, criterion {CRIT}): misfit {fit.fun:.3f}.  rho {rho:+.2f}, Ti {Ti:.1f} yr, k {k:.2f}, "
             f"d {np.degrees(d):.0f} deg, tau {tau:.1f} mo" + (f", saturation G {p[5]:.1f}, b {p[6]:+.2f}" if RESP == "sat" else ""), fontsize=10)
half = len(t) // 2
for ax_, sl, ttl in ((fig.add_subplot(gs[1]), slice(0, half), "first half"), (fig.add_subplot(gs[2]), slice(half, None), "second half")):
    ax_.plot(t[sl], u_obs[sl], color="k", lw=1.3, label="observed")
    ax_.plot(t[sl], um[sl], color="C3", lw=1.1, label="model (no time-domain fitting)")
    ax_.plot(t[sl], -um[sl], color="C0", lw=0.8, ls="--", alpha=0.6, label="-model (same spectrum)")
    ax_.axhline(0, color="0.7", lw=0.5); ax_.set_ylabel("m/s"); ax_.set_title(f"VISUAL COMPARISON ({ttl}): generated from the spectral fit", fontsize=10)
    ax_.legend(fontsize=8, ncol=3, loc="upper right"); ax_.set_xlim(t[sl][0], t[sl][-1])
a3 = fig.add_subplot(gs[3]); a3.plot(t, L, color="C2", lw=0.8)
a3.set_ylabel("L(t)"); a3.set_xlabel("year"); a3.set_title("forcing manifold L(t): equinox-impulsed draconic factor, leaky-integrated", fontsize=10); a3.set_xlim(t[0], t[-1])
fig.tight_layout(); fig.savefig(os.path.join(HERE, f"qbo_spectral_fit_{TAG}.png"), dpi=110); print(f"wrote qbo_spectral_fit_{TAG}.png")
