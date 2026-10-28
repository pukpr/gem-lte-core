"""Test of the "asymmetric (sawtooth) daily carrier" reading of the even-harmonic sidebands (2026-10-26).
Claim (geoenergymath.com/the-just-so-story-narrative/): double sidebands at the 13.66 d fortnight appear only about
the EVEN harmonics of the daily tide, the daily waves having "an asymmetric sawtooth impulse character", perhaps
a "1/2 day (semi-diurnal) modulation".

If one daily waveform g(t) (harmonics S_n) is multiplied by a fortnightly factor (1 + m cos W t), every harmonic gets
the SAME pair of sidebands, each of relative amplitude m/2, odd or even, whatever the shape of g (sawtooth or not).
So the test is the sideband/carrier amplitude ratio as a function of n, and whether the two sidebands are equal.

Method: joint least-squares fit to the eastward k=1 coefficient c(t) (from wn1_spectrum.py, wn1_coef.npy) of complex
lines at S_n (n = 1..11), their annual/semiannual neighbours S_n +/- h, 2h, and the sidebands S_n +/- Mf and
S_n +/- MSf, plus O1/Q1/OO1. Noise per line = rms amplitude of single-frequency fits to the residual at 60 nearby
off-line frequencies.  Eastward wave amplitude = 2|a| (cos(lon - w t) gives c = e^{-i w t}/2).
Output  sideband_test.png, sideband_test.json"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
SR, HR, PR = 1 / 27.321582, 1 / 365.242199, 1 / 3232.6054     # moon, sun, lunar perigee (8.85 yr) rates, cpd
MF, MSF = 2 * SR, 2 * SR - 2 * HR
D = np.load(os.path.join(HERE, "wn1_coef.npy")); t = D[:, 0]; c = D[:, 1] + 1j * D[:, 2]
td = (t - t[0]) / 24.0; c = c - c.mean(); T = td[-1] - td[0]
print(f"{len(t)} hours, {T / 365.25:.1f} yr from {1970 + t[0] / 8766:.1f}")

lines = {}
for n in range(1, 12):
    lines[f"S{n}"] = n
    for k in (1, 2):
        lines[f"S{n}-{k}h"] = n - k * HR; lines[f"S{n}+{k}h"] = n + k * HR
    lines[f"S{n}-Mf"] = n - MF; lines[f"S{n}+Mf"] = n + MF; lines[f"S{n}-MSf"] = n - MSF; lines[f"S{n}+MSf"] = n + MSF
lines.update({"O1": 1 - 2 * SR + HR, "Q1": 1 - 3 * SR + HR + PR, "OO1": 1 + 2 * SR + HR})
names = list(lines); fr = np.array([lines[k] for k in names])
X = np.exp(-2j * np.pi * np.outer(td, fr))                      # eastward lines (negative numpy frequency)
a, *_ = np.linalg.lstsq(X, c, rcond=None); r = c - X @ a


def noise(f0):
    d = np.concatenate([np.linspace(-0.03, -0.008, 30), np.linspace(0.008, 0.03, 30)]); f = f0 + d
    f = f[np.min(np.abs(f[:, None] - fr[None, :]), axis=1) > 5 / T]
    E = np.exp(-2j * np.pi * np.outer(f, td)); return float(np.sqrt(np.mean(np.abs(E @ r / len(td)) ** 2)))


out = {}
for k, f0, ai in zip(names, fr, a):
    out[k] = {"cpd": float(f0), "amp_hPa": float(2 * abs(ai)), "phase_deg": float(np.degrees(np.angle(ai))),
              "noise_hPa": 2 * noise(f0)}
print(f"{'n':>2s} {'carrier hPa':>12s} | {'-MSf/car':>9s} {'+MSf/car':>9s} {'-Mf/car':>9s} {'+Mf/car':>9s} | sideband amp/noise  -MSf +MSf -Mf +Mf")
tab = []
for n in range(1, 12):
    C = out[f"S{n}"]["amp_hPa"]; row = [n, C]
    sb = [out[f"S{n}{s}"] for s in ("-MSf", "+MSf", "-Mf", "+Mf")]
    row += [s["amp_hPa"] / C for s in sb] + [s["amp_hPa"] / s["noise_hPa"] for s in sb] + [s["noise_hPa"] / C for s in sb]
    tab.append(row)
    print(f"{n:2d} {C:12.5f} | " + " ".join(f"{x:9.4f}" for x in row[2:6]) + " | " + " ".join(f"{x:5.1f}" for x in row[6:10]))
tab = np.array(tab)
for k in ("O1", "Q1", "OO1"):
    print(f"{k:4s} amp {out[k]['amp_hPa']:.5f} hPa, {out[k]['amp_hPa'] / out[k]['noise_hPa']:.1f} x noise")

# what a single fortnightly modulation of the whole daily waveform would predict: take the ratio measured at the
# best even harmonic and apply it to every carrier
odd = tab[0::2]; even = tab[1::2]
mE = np.median(even[:, 2:4]); pred = {int(n): float(mE * C) for n, C in tab[:, :2]}
summary = {"years": float(T / 365.25), "median_even_MSf_ratio": float(mE),
           "odd_MSf_ratio_median": float(np.median(odd[:, 2:4])),
           "odd_MSf_sidebands_predicted_by_uniform_modulation_hPa": {k: v for k, v in pred.items() if k % 2},
           "lower_over_upper_MSf_even_n": [float(x) for x in even[:, 2] / even[:, 3]],
           "even_power_fraction_of_carriers": float((tab[1::2, 1] ** 2).sum() / (tab[:, 1] ** 2).sum()),
           "lines": out}
json.dump(summary, open(os.path.join(HERE, "sideband_test.json"), "w"), indent=1)
print(f"median even-n MSf sideband/carrier = {mE:.3f}; odd-n = {summary['odd_MSf_ratio_median']:.4f}")
print("lower/upper MSf sideband ratio at n = 2,4,6,8,10:", np.round(even[:, 2] / even[:, 3], 2))

fig, ax = plt.subplots(2, 2, figsize=(12.5, 9))
n = tab[:, 0]; ev = n % 2 == 0
a0 = ax[0, 0]; a0.semilogy(n[ev], tab[ev, 1], "o", color="C0", label="even n"); a0.semilogy(n[~ev], tab[~ev, 1], "s", color="C3", label="odd n")
nn = np.arange(1, 12); a0.semilogy(nn, tab[0, 1] / nn, "k:", lw=1, label="sawtooth 1/n (from n=1)")
a0.semilogy(nn, tab[0, 1] / nn ** 2, "k--", lw=0.8, label="1/n$^2$ (triangle-like)")
a0.set_xlabel("harmonic n of the solar day"); a0.set_ylabel("eastward k=1 amplitude (hPa)"); a0.legend(fontsize=8)
a0.set_title("Daily carrier harmonics S$_n$ (joint LS fit)", fontsize=10)
a1 = ax[0, 1]
for j, (lab, mk, col) in enumerate((("$-$MSf (14.77 d)", "v", "C0"), ("+MSf", "^", "C0"), ("$-$Mf (13.66 d)", "v", "C1"), ("+Mf", "^", "C1"))):
    a1.semilogy(n + (j - 1.5) * 0.08, tab[:, 2 + j], mk, color=col, ms=6, mfc=col if j < 2 else "none", label=lab)
a1.semilogy(n, np.median(tab[:, 10:14], axis=1), "k_", ms=14, label="noise / carrier")
a1.axhline(mE, color="C0", lw=0.6, ls=":"); a1.text(11.2, mE, "uniform-\nmodulation\nprediction", fontsize=7, va="center", color="C0")
a1.set_xticks(nn); a1.set_xlabel("harmonic n"); a1.set_ylabel("sideband / carrier amplitude"); a1.legend(fontsize=8, ncol=2)
a1.set_title("A fortnightly modulation of one daily waveform gives the SAME ratio at every n", fontsize=10)
# reconstructed daily waveform at lon 0 from the carriers: p(t) = sum 2|a_n| cos(2 pi n t - phase)
tt = np.linspace(0, 2, 961); W = lambda sel: sum(out[f"S{k}"]["amp_hPa"] * np.cos(2 * np.pi * k * tt + np.radians(out[f"S{k}"]["phase_deg"])) for k in sel)
a2 = ax[1, 0]; a2.plot(tt * 24, W(range(1, 12)), "k", lw=1.2, label="all n = 1..11")
a2.plot(tt * 24, W(range(2, 12, 2)), "C0", lw=1, label="even n only (12-h periodic)"); a2.plot(tt * 24, W(range(1, 12, 2)), "C3", lw=1, label="odd n only (half-wave antisymmetric)")
a2.set_xlabel("hour (UTC), two days"); a2.set_ylabel("hPa"); a2.legend(fontsize=8); a2.set_title("Eastward k=1 daily waveform at 0E from the S$_n$ fit", fontsize=10)
a3 = ax[1, 1]; a3.plot(tt * 24, W(range(2, 12)), "k", lw=1, label="n >= 2"); a3.plot(tt * 24, W(range(2, 12, 2)), "C0", lw=1, label="even n >= 2")
a3.plot(tt * 24, W(range(3, 12, 2)), "C3", lw=1, label="odd n >= 3"); a3.legend(fontsize=8); a3.set_xlabel("hour (UTC)")
a3.set_title("Same, fundamental removed: the impulse/asymmetric shape", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "sideband_test.png"), dpi=120); print("wrote sideband_test.png")
