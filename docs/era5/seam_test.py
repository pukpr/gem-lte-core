"""Are the even-harmonic lunar sidebands an ERA5 assimilation-window artifact?  (2026-10-26)
ERA5 runs 4D-Var in two 12-h windows, 09-21 and 21-09 UTC, and its forecast model has NO lunar gravitational forcing
(Sakazaki & Hamilton 2020, sec. 2a), so the lunar M2 tide enters only through the analysis increments.  Anything that
makes the analysed field depend on the position within the 12-h window multiplies the lunar signal by a 12-h periodic
function m(t): an eastward line A e^{-i w t} becomes sum_q B_q e^{-i (w + q) t}, q = 0, +/-2, +/-4 ... cpd, i.e. lines
at M2 + 2j = S_{2k} - MSf, and the westward M2 line becomes eastward lines at q - M2 = S_{2k} + MSf.  Only even
harmonics, and with a 14.77 d (MSf) spacing, not 13.66 d (Mf).

Tests, all on the eastward k=1 coefficient c(t) from wn1_spectrum.py (wn1_coef.npy, contiguous hourly):
  1. step excess by UTC hour: |d_h - (d_{h-1} + d_{h+1})/2| of the anomaly from the mean solar day; seams -> spikes
  2. M2 amplitude fitted separately at each UTC hour, eastward and westward: m(h)
  3. consistency check only: the Fourier coefficients B_q of the per-hour M2 amplitude reproduce the joint LS sideband
     amplitudes of sideband_test.py. This is an IDENTITY (splitting by hour of day and transforming over the 24 hours is
     the same as fitting the lines M2 + q), NOT an independent prediction. What is informative is the TIME-domain shape
     of m(h) in 2.: a 12-h sawtooth reset exactly at the 09|10 and 21|22 UTC seams is what puts sidebands out to n = 10.
Output  seam_test.png, seam_test.json"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
SR, HR = 1 / 27.321582, 1 / 365.242199
M2, MSF = 2 - 2 * SR + 2 * HR, 2 * SR - 2 * HR
D = np.load(os.path.join(HERE, "wn1_coef.npy")); t = D[:, 0]; c = D[:, 1] + 1j * D[:, 2]
assert (np.diff(t) == 1).all() and t[0] % 24 == 0, "needs a contiguous hourly record starting at 00 UTC"
hr = t.astype(int) % 24; td = (t - t[0]) / 24.0; c = c - c.mean()

# 1. step excess by UTC hour
day = np.array([c[hr == h].mean() for h in range(24)]); a = c - day[hr]; d = np.diff(a)
J = np.abs(d[1:-1] - 0.5 * (d[:-2] + d[2:])); into = hr[2:-1]
step = np.array([J[into == h].mean() for h in range(24)])              # step from hour h-1 into hour h

# 2. per-UTC-hour M2 amplitude (eastward: e^{-i w t}; westward: e^{+i w t}), both fitted together with the
#    neighbouring solar lines so the strong S2 does not leak in
def per_hour(sign):
    out = np.empty(24, complex)
    for h in range(24):
        m = hr == h; tt = td[m]
        X = np.column_stack([np.exp(sign * -2j * np.pi * M2 * tt), np.exp(sign * -2j * np.pi * (M2 - 2 * HR) * tt),
                             np.exp(sign * -2j * np.pi * (2 + 2 * HR) * tt), np.exp(sign * -2j * np.pi * (2 - 2 * HR) * tt), np.ones(m.sum())])
        out[h] = np.linalg.lstsq(X, c[m], rcond=None)[0][0]
    return out


aE, aW = per_hour(1), per_hour(-1)
hh = np.arange(24)
B = lambda a, q: np.mean(a * np.exp(2j * np.pi * q * hh / 24))       # a(h) = sum_q B_q e^{-i 2 pi q h/24}

S = json.load(open(os.path.join(HERE, "sideband_test.json")))["lines"]
rows = []
for k in range(1, 6):                                   # carriers S2..S10
    n = 2 * k
    lo_pred = 2 * abs(B(aE, n - 2))                      # eastward M2 shifted by q = n-2  ->  S_n - MSf
    up_pred = 2 * abs(B(aW, n + 2))                      # westward M2 at q - M2, q = n+2 ->  S_n + MSf
    lo, up = S[f"S{n}-MSf"], S[f"S{n}+MSf"]
    rows.append({"n": n, "lower_obs": lo["amp_hPa"], "lower_from_hourly": lo_pred, "lower_noise": lo["noise_hPa"],
                 "upper_obs": up["amp_hPa"], "upper_from_hourly": up_pred, "upper_noise": up["noise_hPa"]})
odd = [{"n": n, "lower_obs": S[f"S{n}-MSf"]["amp_hPa"], "lower_from_hourly": 2 * abs(B(aE, n - 2)),
        "upper_obs": S[f"S{n}+MSf"]["amp_hPa"], "upper_from_hourly": 2 * abs(B(aW, n + 2))} for n in (1, 3, 5, 7)]
print("step excess by UTC hour (x background):", " ".join(f"{h:02d}:{step[h] / np.median(step):.1f}" for h in range(24)))
print("eastward M2 |a(h)| hPa:", np.round(2 * np.abs(aE), 4)); print("westward M2 |a(h)| hPa:", np.round(2 * np.abs(aW), 4))
print("identity check (not a prediction): joint-LS sideband vs per-hour-fit Fourier coefficient")
print(" n   lower LS  hourly  (noise)  |  upper LS  hourly  (noise)")
for r in rows:
    print(f"{r['n']:2d}  {r['lower_obs']:.5f} {r['lower_from_hourly']:.5f} ({r['lower_noise']:.5f}) | {r['upper_obs']:.5f} {r['upper_from_hourly']:.5f} ({r['upper_noise']:.5f})")
print("odd n:")
for r in odd:
    print(f"{r['n']:2d}  {r['lower_obs']:.5f} {r['lower_from_hourly']:.5f} | {r['upper_obs']:.5f} {r['upper_from_hourly']:.5f}")
obs = np.array([[r["lower_obs"], r["upper_obs"]] for r in rows]).ravel(); pr = np.array([[r["lower_from_hourly"], r["upper_from_hourly"]] for r in rows]).ravel()
rlog = float(np.max(np.abs(obs / pr - 1))); print(f"identity holds to max relative difference {rlog:.1e}")
json.dump({"step_excess_by_hour": step.tolist(), "M2_east_amp_by_hour": (2 * np.abs(aE)).tolist(),
           "M2_west_amp_by_hour": (2 * np.abs(aW)).tolist(), "even": rows, "odd": odd, "identity_max_rel_diff": rlog},
          open(os.path.join(HERE, "seam_test.json"), "w"), indent=1)

fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
a0 = ax[0]; a0.bar(hh - 0.5, step / np.median(step), width=0.9, color=["C3" if h in (10, 22) else "0.6" for h in hh])
a0.set_xticks(range(0, 25, 3)); a0.set_xlabel("UTC hour (step from h-1 to h)"); a0.set_ylabel("step excess / median")
a0.set_title("Hour-to-hour jumps: ERA5 4D-Var seams at 09|10 and 21|22 UTC", fontsize=10)
a1 = ax[1]; a1.plot(hh, 2 * np.abs(aE), "o-", label="eastward M2"); a1.plot(hh, 2 * np.abs(aW), "s-", label="westward M2")
for x in (9.5, 21.5): a1.axvline(x, color="C3", lw=0.8, ls="--")
a1.set_xticks(range(0, 24, 3)); a1.set_xlabel("UTC hour"); a1.set_ylabel("k=1 M2 amplitude (hPa)"); a1.legend(fontsize=8)
a1.set_title("Lunar M2 amplitude fitted separately by UTC hour", fontsize=10)
a2 = ax[2]; nn = [r["n"] for r in rows]
a2.semilogy(nn, [r["lower_obs"] for r in rows], "v-", color="C0", label="S$_n$ $-$ MSf (= M2 + S$_{n-2}$)")
a2.semilogy(nn, [r["upper_obs"] for r in rows], "^-", color="C1", label="S$_n$ + MSf (= S$_{n+2}$ $-$ M2)")
a2.semilogy([r["n"] for r in odd], [r["lower_obs"] for r in odd], "x", color="0.4", label="odd n, either side")
a2.semilogy([r["n"] for r in odd], [r["upper_obs"] for r in odd], "x", color="0.4")
a2.semilogy(nn, [r["lower_noise"] for r in rows], "k_", ms=14, label="noise")
a2.set_xlabel("harmonic n"); a2.set_ylabel("eastward k=1 line amplitude (hPa)"); a2.legend(fontsize=8)
a2.set_title("MSf sidebands: flat to n = 10 at even n, absent at odd n\n(the frequency-domain image of the 12-h reset)", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "seam_test.png"), dpi=120); print("wrote seam_test.png")
