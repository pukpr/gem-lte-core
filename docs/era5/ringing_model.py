"""Corrected "atmospheric ringing" model of the ERA5 eastward k=1 tropical surface-pressure signal (2026-10-26).
Recreates the two-panel figure of geoenergymath.com/the-just-so-story-narrative/ (forcing in the time domain,
model vs data power spectrum), correcting the forcing by what wn1_spectrum.py / sideband_test.py / seam_test.py found:

  blog:      daily impulse train  x  (1 + m cos 2pi t / 13.66 d)   + ringing modes
  corrected: D(t)  solar daily waveform, harmonics S_n (n = 1..11) with +/-h, +/-2h seasonal neighbours  (LS fit)
           + L(t)  lunar M2 (eastward and westward) multiplied by the ERA5 12-h analysis-window function m(UTC hour),
                   reset at the 09|10 and 21|22 UTC 4D-Var seams (per-hour fit), plus O1 (= K1 - Mf, the real 13.66 d line)
           + E(t)  stochastic daily patterns: two families, each ONE fixed waveform with a slowly wandering random
                   amplitude s(t) (complex Lorentzian / AR(1)): an even family (12-h periodic, harmonics 2..10) and an odd
                   family (24-h, half-wave antisymmetric, harmonics 3..11). Waveforms = leading eigenvector of the
                   demodulated envelopes across harmonics (coherence 0.73-0.91 even, 0.56-0.81 odd, <=0.16 between);
                   they produce the broadened pedestals. The measured amplitude s(t) is phase-randomised for the simulation. (A multiplicative 12-h modulation of the slow field was tested
                   and is ~1e4 too weak; independent per-window ramps give broadband noise, not pedestals.)
           + R(t)  ringing free modes: complex damped oscillators driven by white noise, f0 and half-width d from
                   Lorentzian fits (S&H 2020 Eq. 5) to the eastward spectrum near 33 h, 9.4 h and 7.3 h
           + B(t)  red background, random phase, smooth fit to the eastward spectrum away from lines and modes
Everything is estimated from wn1_coef.npy (the record wn1_spectrum.py last wrote; contiguous hourly).
Spectra by the paper's method (365-day yearly DFTs, detrended, averaged).
Output  ringing_model.png, ringing_model.json"""
import json, os
import numpy as np
from scipy.optimize import least_squares
from scipy.signal import lfilter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(20261026)
SR, HR, PR = 1 / 27.321582, 1 / 365.242199, 1 / 3232.6054     # moon, sun, lunar perigee (8.85 yr) rates, cpd
MF, MSF = 2 * SR, 2 * SR - 2 * HR
M2, O1 = 2 - 2 * SR + 2 * HR, 1 - 2 * SR + HR
D_ = np.load(os.path.join(HERE, "wn1_coef.npy")); t = D_[:, 0]; c = D_[:, 1] + 1j * D_[:, 2]
assert (np.diff(t) == 1).all() and t[0] % 24 == 0
N = len(t); td = (t - t[0]) / 24.0; hr = t.astype(int) % 24; c = c - c.mean()
Y0 = 1970 + t[0] / 8766; print(f"{N} hours, {N / 8766:.1f} yr from {Y0:.1f}")


def yearly(x):
    """paper method: first 365 d of each block of 8760 h (record starts 1 Jan), detrended, 2|q|^2/df averaged."""
    n = 8760; S = []
    for i0 in range(0, len(x) - n + 1, 8784 if False else n):
        s = x[i0:i0 + n]; A = np.vstack([np.arange(n), np.ones(n)]).T; s = s - A @ np.linalg.lstsq(A, s, rcond=None)[0]
        S.append(2 * 365.0 * np.abs(np.fft.fft(s) / n) ** 2)
    f = -np.fft.fftfreq(n, 1 / 24.0); o = np.argsort(f); return f[o], np.mean(S, axis=0)[o]   # f > 0 = eastward


fE, PE = yearly(c)

# ---- D: solar daily waveform (eastward lines) ---------------------------------------------------------------------
sol = [n + k * HR for n in range(1, 12) for k in (-2, -1, 0, 1, 2)]
X = np.exp(-2j * np.pi * np.outer(td, sol)); aS = np.linalg.lstsq(X, c, rcond=None)[0]; D = X @ aS

# ---- L: lunar M2 x analysis-window function (per UTC hour), plus O1 ----------------------------------------------
def per_hour(sign):
    a = np.empty(24, complex)
    for h in range(24):
        m = hr == h; tt = td[m]
        Xh = np.column_stack([np.exp(sign * -2j * np.pi * f * tt) for f in (M2, M2 - 2 * HR, 2 + 2 * HR, 2 - 2 * HR)] + [np.ones(m.sum())])
        a[h] = np.linalg.lstsq(Xh, c[m], rcond=None)[0][0]
    return a


aE, aW = per_hour(1), per_hour(-1)
aO1 = np.vdot(np.exp(-2j * np.pi * O1 * td), c - D) / N
L = aE[hr] * np.exp(-2j * np.pi * M2 * td) + aW[hr] * np.exp(2j * np.pi * M2 * td) + aO1 * np.exp(-2j * np.pi * O1 * td)

# ---- seam statistics from data (anomaly from the mean solar day) --------------------------------------------------
def seam_steps(x):
    day = np.array([x[hr == h].mean() for h in range(24)]); a = x - day[hr]; d = np.diff(a)
    exc = d[1:-1] - 0.5 * (d[:-2] + d[2:]); into = hr[2:-1]
    return exc, into


east = lambda x: np.fft.ifft(np.fft.fft(x) * (np.fft.fftfreq(len(x)) <= 0))    # keep eastward (numpy f <= 0) side
exc, into = seam_steps(east(c - L))
Jm = lambda e, i: np.array([np.abs(e[i == h]).mean() for h in range(24)])
Jd = Jm(exc, into); bgJ = np.median(Jd)
sk = exc[(into == 10) | (into == 22)]                      # successive seams, alternating 10 and 22 UTC
phi = np.vdot(sk[:-1], sk[1:]) / np.vdot(sk, sk); phi = float(abs(phi)) * np.exp(1j * np.angle(phi))
print(f"seam step excess, eastward part (data): 10 UTC {Jd[10] / bgJ:.1f}x, 22 UTC {Jd[22] / bgJ:.1f}x; seam-to-seam correlation {abs(phi):.2f}")

# ---- R: Lorentzian fits to the eastward spectrum (paper Eq. 5) ----------------------------------------------------
def lorentz_fit(f0, w):
    m = (fE > f0 - w) & (fE < f0 + w) & (np.abs(fE - np.round(fE)) > 0.05)
    ff, pp = fE[m], PE[m]
    mod = lambda q: q[0] * q[1] / ((ff - q[2]) ** 2 + q[1] ** 2) + q[3] * ff + q[4]
    q0 = [pp.max() * 0.05 * 0.3, 0.03, f0, 0.0, np.median(pp)]
    r = least_squares(lambda q: np.log(np.maximum(mod(q), 1e-30)) - np.log(pp), q0,
                      bounds=([0, 0.01, f0 - w / 2, -np.inf, 0], [np.inf, 0.3, f0 + w / 2, np.inf, np.inf]))
    A, d, fo = r.x[:3]; return {"f0_cpd": float(fo), "period_h": float(24 / fo), "halfwidth_cpd": float(d),
                                "peak_psd": float(A / d), "decay_days": float(1 / (2 * np.pi * d))}


modes = {"33h Kelvin": lorentz_fit(0.727, 0.15), "9.4h": lorentz_fit(2.55, 0.2), "7.3h": lorentz_fit(3.29, 0.2)}
for k, v in modes.items():
    print(f"mode {k:10s} f0 {v['f0_cpd']:.3f} cpd ({v['period_h']:.1f} h), half-width {v['halfwidth_cpd']:.3f} cpd, e-fold {v['decay_days']:.1f} d")
fs = [v["f0_cpd"] for v in modes.values()]; print(f"triad check f(33h) + f(9.4h) = {fs[0] + fs[1]:.3f} vs f(7.3h) = {fs[2]:.3f} cpd")


def osc(f0, d):
    """complex damped oscillator, eastward rotation, white-noise driven, unit innovation"""
    rho = np.exp((-2 * np.pi * d - 2j * np.pi * f0) / 24.0); e = (rng.standard_normal(N) + 1j * rng.standard_normal(N)) / np.sqrt(2)
    return lfilter([1.0], [1.0, -rho], e)


# ---- B: background = running median of the data spectrum, away from lines, pedestals (+/-0.3 cpd) and modes ---------
excl = (np.abs(fE - np.round(fE)) < 0.3) & (np.round(fE) >= 2) | (np.abs(fE - np.round(fE)) < 0.06) & (np.round(fE) >= 1)
for v in modes.values(): excl |= np.abs(fE - v["f0_cpd"]) < 4 * v["halfwidth_cpd"]
mb = (fE > 0) & ~excl; fb_, pb_ = fE[mb], np.log(PE[mb]); half = 40
med = np.array([np.median(pb_[max(0, i - half):i + half + 1]) for i in range(len(pb_))])
Gb = lambda f: np.exp(np.interp(f, fb_, med))


def colored(P):
    """random-phase series with eastward-only spectrum shape P(f) (f in cpd, eastward > 0)"""
    f = -np.fft.fftfreq(N, 1 / 24.0); amp = np.where(f > 0, np.sqrt(P(np.abs(f) + 1e-6)), 0.0)
    return np.fft.ifft(amp * (rng.standard_normal(N) + 1j * rng.standard_normal(N)))


B = colored(Gb)

# ---- E: stochastic daily-pattern families from demodulated envelopes ------------------------------------------------
MSF_lines = [n + sg * MSF for n in range(1, 12) for sg in (-1, 1)]
Xl = np.exp(-2j * np.pi * np.outer(td, sol + MSF_lines)); resid = c - Xl @ np.linalg.lstsq(Xl, c, rcond=None)[0]
fF = np.fft.fftfreq(N, 1 / 24.0); FR = np.fft.fft(resid); BW = 0.3
demod = lambda q: np.fft.ifft(np.roll(FR, int(round(q * N / 24))) * (np.abs(fF) < BW))   # eastward q cpd -> 0


def family(qs):
    Eq = np.array([demod(q) for q in qs]); C = Eq @ Eq.conj().T / N
    lam, V = np.linalg.eigh(C); g = V[:, -1]; sq = g.conj() @ Eq                 # leading pattern and its amplitude s(t)
    # Lorentzian (centred) + floor fitted to the periodogram of s within the band
    S = np.abs(np.fft.fft(sq)) ** 2 / N / 24.0; m = (np.abs(fF) < 0.25) & (fF != 0); ff, pp = fF[m], S[m]
    nb = 400; o = np.argsort(ff); ff, pp = ff[o], pp[o]; k = len(ff) // nb
    fb = ff[:k * nb].reshape(nb, k).mean(1); pb = pp[:k * nb].reshape(nb, k).mean(1)
    r = least_squares(lambda q: np.log(q[0] * q[1] / (fb ** 2 + q[1] ** 2) + q[2]) - np.log(pb), [pb.max() * 0.03, 0.03, pb.min()],
                      bounds=([0, 0.003, 0], [np.inf, 0.3, np.inf]))
    A, d, fl = r.x
    return {"q": list(qs), "g": g, "share": float(lam[-1] / lam.sum()), "halfwidth_cpd": float(d), "var": float(np.mean(np.abs(sq) ** 2)),
            "efold_days": float(1 / (2 * np.pi * d)), "s": sq}


fams = {"even (12-h)": family([2, 4, 6, 8, 10, 12]), "odd (24-h)": family([3, 5, 7, 9, 11])}
E = np.zeros(N, complex); Es = {}
for k, F_ in fams.items():
    sim = np.fft.ifft(np.fft.fft(F_["s"]) * np.exp(2j * np.pi * rng.random(N)))     # phase-randomised measured amplitude
    Es[k] = sum(gq * sim * np.exp(-2j * np.pi * q * td) for q, gq in zip(F_["q"], F_["g"])); E += Es[k]
    print(f"family {k}: leading-pattern share {F_['share']:.2f}, envelope half-width {F_['halfwidth_cpd']:.3f} cpd "
          f"(Lorentzian e-fold {F_['efold_days']:.1f} d, descriptive only), rms {np.sqrt(F_['var']):.4f} hPa")
hh = np.linspace(0, 24, 289)
wave = {k: sum(gq * np.exp(-2j * np.pi * q * hh / 24) for q, gq in zip(F_["q"], F_["g"])) for k, F_ in fams.items()}


def scale_to(x, target_ratio_P, f0, f1):
    """scale a component so its yearly spectrum median in [f0,f1] matches target"""
    fx, Px = yearly(x); m = (fx > f0) & (fx < f1); return np.sqrt(np.median(target_ratio_P[m]) / np.median(Px[m]))


# calibrate amplitudes: background to the data background; modes to the Lorentzian peak excess; seam ramps to the
# measured seam-step excess (two passes, the step statistic is ~linear in the ramp amplitude)
B *= scale_to(B, Gb(fE), 5.2, 5.8)
R = np.zeros(N, complex); Rs = {}
for k, v in modes.items():
    z = osc(v["f0_cpd"], v["halfwidth_cpd"]); fz, Pz = yearly(z)
    pk = np.median(Pz[np.abs(fz - v["f0_cpd"]) < 0.3 * v["halfwidth_cpd"] + 0.003]); z *= np.sqrt(v["peak_psd"] / pk); Rs[k] = z; R += z
model = D + L + E + R + B
fM, PM = yearly(model)
eM, iM = seam_steps(east(model - L)); JM = Jm(eM, iM)
print(f"model seam step excess (not fitted): 10 UTC {JM[10] / np.median(JM):.1f}x, 22 UTC {JM[22] / np.median(JM):.1f}x (data {Jd[10] / bgJ:.1f}x, {Jd[22] / bgJ:.1f}x)")

# blog-style alternative: whole daily waveform x (1 + m cos(2 pi Mf t)), m from the observed even sideband ratio
mBlog = 2 * 0.22; blog = D * (1 + mBlog * np.cos(2 * np.pi * MF * td)) + R + B
fB, PB = yearly(blog)

# line check on the model vs data, same joint-LS approach as sideband_test.py
def sidebands(x):
    fr = [n + s for n in range(1, 11) for s in (0, -MSF, MSF, -MF, MF)]
    Xl = np.exp(-2j * np.pi * np.outer(td, fr)); a = np.linalg.lstsq(Xl, x, rcond=None)[0]
    return 2 * np.abs(a).reshape(10, 5)


SD, SM, SB = sidebands(c), sidebands(model), sidebands(blog)
print(" n  data: -MSf/car +MSf/car -Mf/car | corrected model | blog model (13.66 d)")
for n in range(10):
    print(f"{n + 1:2d}  " + " ".join(f"{x:7.3f}" for x in SD[n, 1:4] / SD[n, 0]) + " | " + " ".join(f"{x:7.3f}" for x in SM[n, 1:4] / SM[n, 0])
          + " | " + " ".join(f"{x:7.3f}" for x in SB[n, 1:4] / SB[n, 0]))
def sb_err(S):
    """rms log10 error of the 40 sideband/carrier ratios (n = 1..10; -MSf, +MSf, -Mf, +Mf) against the data"""
    return float(np.sqrt(np.mean((np.log10(S[:, 1:] / S[:, :1]) - np.log10(SD[:, 1:] / SD[:, :1])) ** 2)))


def ped(fx, Px):
    """pedestal contrast at n = 2..11: median PSD within +/-0.15 cpd (central +/-0.01 excluded) over flanks at 0.3-0.45"""
    out = []
    for n in range(2, 12):
        d = np.abs(fx - n); band = (d > 0.01) & (d < 0.15); flank = (d > 0.3) & (d < 0.45)
        for sg in (-1, 1): band &= np.abs(fx - (n + sg * MSF)) > 0.004
        out.append(float(np.median(Px[band]) / np.median(Px[flank])))
    return out


errM, errB = sb_err(SM), sb_err(SB); pD, pM, pB = ped(fE, PE), ped(fM, PM), ped(fB, PB)
print(f"sideband-ratio rms log10 error vs data: corrected {errM:.2f}, blog-style {errB:.2f}")
print("pedestal contrast n=2..11  data:     ", " ".join(f"{x:.2f}" for x in pD))
print("                           corrected:", " ".join(f"{x:.2f}" for x in pM))
print("                           blog:     ", " ".join(f"{x:.2f}" for x in pB))
json.dump({"record": [float(Y0), float(N / 8766)], "modes": modes, "families": {k: {kk: (vv if kk != "g" else [[float(z.real), float(z.imag)] for z in vv]) for kk, vv in F_.items() if kk != "s"} for k, F_ in fams.items()},
           "seam_step_excess_data": [float(Jd[10] / bgJ), float(Jd[22] / bgJ)],
           "seam_step_excess_model": [float(JM[10] / np.median(JM)), float(JM[22] / np.median(JM))],
           "blog_m": mBlog, "sideband_rms_log10_err": {"corrected": errM, "blog": errB}, "pedestal_contrast": {"data": pD, "corrected": pM, "blog": pB},
           "sideband_ratios": {"data": SD.tolist(), "corrected": SM.tolist(), "blog": SB.tolist()},
           "window_function_M2_east": (2 * np.abs(aE)).tolist(), "window_function_M2_west": (2 * np.abs(aW)).tolist()},
          open(os.path.join(HERE, "ringing_model.json"), "w"), indent=1)

# ---- figure ---------------------------------------------------------------------------------------------------------
fig = plt.figure(figsize=(14, 15)); gs = fig.add_gridspec(4, 3, height_ratios=[1.15, 0.8, 1.35, 1.0])
d40 = np.arange(40 * 24); x40 = td[d40]
a = fig.add_subplot(gs[0, :])
p = lambda z: 2 * z.real                                      # eastward k=1 pressure at 0E
a.plot(x40, p(D[d40] + L[d40]), color="C3", lw=0.8, label="deterministic forcing D + L at 0E (solar daily waveform + lunar)")
a.plot(x40, p(L[d40]) * 3 - 0.10, color="C0", lw=0.7, label="lunar part L x3 (M2 through the 12-h analysis windows), offset")
env = 2 * np.abs(aW[hr[d40]]) * 3; a.plot(x40, env - 0.10, color="C0", lw=1.2, alpha=0.6, label="window envelope |m(UTC hour)| x3: 12-h sawtooth, reset 10 & 22 UTC")
sn = 2 * np.abs(aS[sol.index(2)] + aE[hr[d40]] * np.exp(-2j * np.pi * (M2 - 2) * x40)) * 3
a.plot(x40, sn - 0.19, color="k", lw=1, label="semidiurnal |S2 + M2| x3: spring-neap beat, MSf = 14.77 d")
a.plot(x40, 0.012 * np.cos(2 * np.pi * MF * x40 + np.pi) - 0.19 + np.mean(sn), color="g", ls="--", lw=1, label="blog: 13.66 d (Mf) modulation (phase arbitrary)")
a.set_xlim(0, 40); a.set_xticks(range(0, 41, 2)); a.set_xlabel(f"days from 1 Jan {Y0:.0f}"); a.set_ylabel("hPa")
a.set_title("Corrected forcing, eastward k=1, 20S-20N (ERA5 1999-2016 fits)", fontsize=11); a.legend(fontsize=7.5, ncol=2, loc="upper right")
a = fig.add_subplot(gs[1, :])
for k, col in zip(Rs, ("C4", "C2", "C1")):
    a.plot(x40, p(Rs[k][d40]), lw=0.7, color=col, label=f"ringing mode {k}: {modes[k]['period_h']:.1f} h, e-fold {modes[k]['decay_days']:.1f} d")
a.plot(x40, p(E[d40]), lw=0.6, color="0.3", label="stochastic daily patterns E (even + odd families)")
a.set_xlim(0, 40); a.set_xticks(range(0, 41, 2)); a.set_ylabel("hPa"); a.legend(fontsize=7.5, ncol=4, loc="upper right")
a.set_title("Noise-driven components (one realisation): free ringing modes and the two wandering daily patterns", fontsize=10)
ai = a.inset_axes([0.005, 0.56, 0.17, 0.42])
for k, col in zip(wave, ("C0", "C3")):
    w_ = wave[k] * np.sqrt(fams[k]["var"]); ai.plot(hh, 2 * w_.real, color=col, lw=1, label=k)
for x in (9.5, 21.5): ai.axvline(x, color="k", lw=0.6, ls=":")
ai.set_xticks([0, 6, 9.5, 12, 18, 21.5, 24]); ai.set_xticklabels(["0", "6", "9|10", "12", "18", "21|22", "24"], fontsize=6)
ai.tick_params(labelsize=6); ai.legend(fontsize=5.5, loc="lower left"); ai.set_title("family waveforms vs UTC hour", fontsize=6.5)
a = fig.add_subplot(gs[2, :])
m12 = (fE > 0) & (fE <= 12)
a.semilogy(fE[m12], PE[m12], color="0.6", lw=0.6, label="ERA5 data")
a.semilogy(fM[m12], PM[m12], color="C3", lw=0.5, alpha=0.85, label="corrected model")
for k, v in modes.items(): a.annotate(k, (v["f0_cpd"], v["peak_psd"] * 2), fontsize=8, ha="center")
a.set_xlim(0, 12); a.set_xticks(range(13)); a.set_xlabel("eastward frequency (cpd)"); a.set_ylabel("PSD (hPa$^2$/cpd)")
a.set_title(f"Power spectrum, paper method. Sideband-ratio rms log10 error vs data: corrected {errM:.2f}, blog-style (13.66 d) {errB:.2f}", fontsize=10); a.legend(fontsize=8)
for j, (ctr, ttl) in enumerate(((1.0, "about 1 cpd"), (3.0, "about 3 cpd (odd)"), (4.0, "about 4 cpd (even)"))):
    a = fig.add_subplot(gs[3, j]); m = (fE > ctr - 0.14) & (fE < ctr + 0.14)
    a.semilogy(fE[m], PE[m], color="0.5", lw=1.4, label="data"); a.semilogy(fM[m], PM[m], color="C3", lw=0.9, label="corrected")
    a.semilogy(fB[m], PB[m], color="g", lw=0.8, ls="--", label="blog: x(1+m cos Mf)")
    for s, col in ((-MSF, "C0"), (MSF, "C0"), (-MF, "g"), (MF, "g")): a.axvline(ctr + s, color=col, lw=0.6, ls=":")
    a.set_title(f"{ttl}: blue = MSf, green = Mf", fontsize=9.5); a.set_xlabel("eastward frequency (cpd)")
    if j == 0: a.legend(fontsize=7.5)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "ringing_model.png"), dpi=110); print("wrote ringing_model.png")
