"""Eastward zonal-wavenumber-1 spectrum of ERA5 hourly tropical surface pressure (20S-20N), 2026-10-25.
Reproduces the analysis discussed at geoenergymath.com/the-just-so-story-narrative/ (after Sakazaki & Hamilton 2020, JAS):
look for ~13.66 day (Mf) satellite sidebands about the daily carrier.

Steps: cos(lat)-weighted mean over 20S-20N -> p(t, lon); wavenumber-1 coefficient c(t) = <p e^{-i lon}>_lon;
the FFT of the complex c(t) separates propagation: with this sign convention an EASTWARD wave cos(lon - w t)
appears at NEGATIVE frequency -w, a westward one (the sun-following migrating tide) at +w (checked below on a synthetic).
Whole-record periodogram (Hann taper), so line frequencies are resolved to ~1/T: with T of decades the
S1 +/- Mf sidebands (0.92680 / 1.07320 cpd) are separable from the lunar K1 +/- Mf = O1 / OO1 (0.92954 / 1.07594 cpd).
Input  data/era5_sp_hourly_*.nc (fetch_sp.py).  Output  wn1_spectrum.png, wn1_lines.json, wn1_coef.npy
usage: python3 wn1_spectrum.py [first year, default 1979] [last year, default 2016 as in the paper]"""
import glob, json, os, sys
import netCDF4 as nc
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
YA = int(sys.argv[1]) if len(sys.argv) > 1 else 1979
YB = int(sys.argv[2]) if len(sys.argv) > 2 else 2016
# mean-longitude rates (cpd): moon s, sun h, lunar perigee p; tidal lines as Doodson combinations with the solar day
SR, HR, PR = 1 / 27.321582, 1 / 365.242199, 1 / 3232.6054     # moon, sun, lunar perigee (8.85 yr) rates, cpd
MF, MSF, MM = 2 * SR, 2 * SR - 2 * HR, SR - PR                       # 13.661 d, 14.765 d, 27.555 d
S1, K1, P1, O1, OO1, Q1 = 1.0, 1 + HR, 1 - HR, 1 - 2 * SR + HR, 1 + 2 * SR + HR, 1 - 3 * SR + HR + PR
S2, K2, M2, N2, L2 = 2.0, 2 + 2 * HR, 2 - 2 * SR + 2 * HR, 2 - 3 * SR + 2 * HR + PR, 2 - SR + 2 * HR - PR
# identities: O1 = K1 - Mf, OO1 = K1 + Mf, M2 = S2 - MSf, S2 + Mf = K2 + MSf; S1 -/+ Mf = O1/OO1 - h


def load():
    T, C = [], []
    for fn in sorted(glob.glob(os.path.join(HERE, "data", "era5_sp_hourly_*.nc"))):
        if not YA <= int(os.path.basename(fn)[15:19]) <= YB:
            continue
        ds = nc.Dataset(fn)
        tv = ds.variables["valid_time" if "valid_time" in ds.variables else "time"]
        hours = nc.date2num(nc.num2date(tv[:], tv.units), "hours since 1970-01-01")
        lat = np.asarray(ds.variables["latitude"][:]); lon = np.asarray(ds.variables["longitude"][:])
        sp = np.asarray(ds.variables["sp"][:], float).squeeze() / 100.0           # Pa -> hPa, (time, lat, lon)
        w = np.cos(np.deg2rad(lat)); p = (sp * w[None, :, None]).sum(axis=1) / w.sum()
        C.append((p * np.exp(-1j * np.deg2rad(lon))[None, :]).mean(axis=1)); T.append(np.asarray(hours))
        ds.close()
    t = np.concatenate(T); c = np.concatenate(C); o = np.argsort(t); t, c = t[o], c[o]
    keep = np.concatenate([[True], np.diff(t) > 0]); t, c = t[keep], c[keep]
    gaps = np.diff(t) != 1
    if gaps.any():
        print(f"note: {gaps.sum()} gaps in the hourly record (years still downloading?); the direct DFT uses actual times")
    return t, c


def check_sign():
    lon = np.deg2rad(np.arange(0, 360, 5.0)); t = np.arange(240) / 24.0
    east = np.cos(lon[None, :] - 2 * np.pi * 1.0 * t[:, None])
    c = (east * np.exp(-1j * lon)[None, :]).mean(axis=1); f = np.fft.fftfreq(len(t), 1 / 24.0)
    assert np.isclose(f[np.argmax(np.abs(np.fft.fft(c)))], -1.0), "sign convention"


def dft(t, c, freqs):
    """Direct DFT at the actual sample times (hours), Hann taper over the full span; freqs in cpd, signed
    (numpy sign: eastward lines at negative f). Works across gaps without interpolation.  Returns hPa^2 per cpd."""
    td = (t - t[0]) / 24.0; w = np.hanning(len(t)) if not (np.diff(t) != 1).any() else 0.5 - 0.5 * np.cos(2 * np.pi * td / td[-1])
    x = (c - c.mean()) * w; P = np.empty(len(freqs))
    for i0 in range(0, len(freqs), 64):
        fr = freqs[i0:i0 + 64]; P[i0:i0 + 64] = np.abs(np.exp(-2j * np.pi * np.outer(fr, td)) @ x) ** 2
    return P / (w ** 2).sum() / 24.0


def yearly_spectrum(t, c):
    """Sakazaki & Hamilton (2020) section 2b: first 365 days of each calendar year, linear trend removed,
    untapered DFT (df = 1/365 cpd), P = 2|q|^2/df, averaged over the years. Same e^{-i lon} projection as above."""
    import datetime as dt
    base = dt.datetime(1970, 1, 1); n = 365 * 24; S = []; yrs = []
    for y in range(int(1970 + t[0] / 8766) - 1, int(1970 + t[-1] / 8766) + 2):
        h0 = (dt.datetime(y, 1, 1) - base).total_seconds() / 3600; i0 = np.searchsorted(t, h0)
        if i0 + n > len(t) or t[i0] != h0 or t[i0 + n - 1] != h0 + n - 1:
            continue
        x = c[i0:i0 + n]; tt = np.arange(n); A = np.vstack([tt, np.ones(n)]).T
        x = x - A @ np.linalg.lstsq(A, x, rcond=None)[0]
        q = np.fft.fft(x) / n; S.append(2 * 365.0 * np.abs(q) ** 2); yrs.append(y)
    return np.fft.fftfreq(n, 1 / 24.0), (np.mean(S, axis=0) if S else None), yrs


def peak(f, P, f0, half):
    m = (f > f0 - half) & (f < f0 + half); i = np.argmax(P[m]); return float(f[m][i]), float(P[m][i])


def bg(f, P, f0, lo, hi):
    m = (np.abs(f - f0) > lo) & (np.abs(f - f0) < hi); return float(np.median(P[m]))


check_sign()
t, c = load(); np.save(os.path.join(HERE, "wn1_coef.npy"), np.column_stack([t, c.real, c.imag]))
years = len(t) / 24 / 365.25; df = 1 / (len(t) / 24.0)
y0 = 1970 + t[0] / 24 / 365.25; print(f"{len(t)} hours, {years:.2f} yr from {y0:.2f}; Rayleigh resolution {df:.2e} cpd")
span = (t[-1] - t[0]) / 24.0; df = 1 / span; print(f"span {span / 365.25:.2f} yr; resolution 1/span = {df:.2e} cpd")
res_step = df / 4

# candidate lines: every harmonic n = 1..6 with Mf and MSf sidebands, plus the named lunar/lunisolar tides
cand = {}
for n in range(1, 7):
    cand.update({f"S{n}": float(n), f"S{n}-Mf": n - MF, f"S{n}+Mf": n + MF, f"S{n}-MSf": n - MSF, f"S{n}+MSf": n + MSF})
cand.update({"K1": K1, "P1": P1, "O1": O1, "OO1": OO1, "Q1": Q1, "K2": K2, "M2": M2, "N2": N2, "L2": L2,
             "M3": 1.5 * M2, "M4": 2 * M2, "MS4": M2 + 2, "S1-Mm": 1 - MM, "S1+Mm": 1 + MM})
half = 2.5 * df; res = {"hours": int(len(t)), "years": years, "start": y0, "span_yr": span / 365.25, "rayleigh_cpd": df,
                        "snr_note": "peak power over the median of 82 frequencies 0.01-0.03 cpd either side; "
                                    "for a chi^2(2) background SNR 10 is p~1e-3 per independent bin", "lines": {}}
bgf = np.concatenate([np.linspace(-0.03, -0.01, 41), np.linspace(0.01, 0.03, 41)])
for direction, sgn in (("east", -1), ("west", 1)):
    for k, f0 in cand.items():
        ff = f0 + np.arange(-half, half + res_step, res_step); pk = dft(t, c, sgn * ff); i = np.argmax(pk)
        fp, pp = float(ff[i]), float(pk[i]); b = float(np.median(dft(t, c, sgn * (f0 + bgf))))
        res["lines"][f"{direction}:{k}"] = {"nominal_cpd": f0, "peak_cpd": fp, "offset_in_bins": (fp - f0) / df,
                                            "power": pp, "amp_hPa": float(np.sqrt(pp * df)), "snr": pp / b}
print(f"whole record: line search +/-{half / df:.1f} Rayleigh bins; offset in bins; SNR = peak / local median")
for direction in ("east", "west"):
    print(f"  {direction:5s}" + "".join(f"{k:>16s}" for k in ("carrier", "-Mf", "+Mf", "-MSf", "+MSf")))
    for n in range(1, 7):
        row = [res["lines"][f"{direction}:S{n}{x}"] for x in ("", "-Mf", "+Mf", "-MSf", "+MSf")]
        print(f"  n={n}  " + "".join(f"{r['snr']:9.1f} ({r['offset_in_bins']:+4.1f})" for r in row))
    print("  lunar " + "  ".join(f"{k} {res['lines'][f'{direction}:{k}']['snr']:.1f}" for k in
                                 ("K1", "P1", "O1", "OO1", "Q1", "M2", "N2", "L2", "K2", "M3", "M4", "MS4")))

# paper method (Fig. 4 replica) + AR(1) background fitted to the log spectrum, harmonics +/-0.05 cpd excluded (their Eq. 4)
fy, Py, yrs = yearly_spectrum(t, c)
fe_y = -fy; m12 = (fe_y > 0) & (fe_y <= 12); fx, Px = fe_y[m12], Py[m12]; o = np.argsort(fx); fx, Px = fx[o], Px[o]
keep = np.abs(fx - np.round(fx)) > 0.05; keep &= np.round(fx) > 0
from scipy.optimize import least_squares
G = lambda q, f: q[0] ** 2 / (1 + q[1] ** 2 - 2 * q[1] * np.cos(2 * np.pi * f / 24.0))
fit = least_squares(lambda q: np.log(G(q, fx[keep])) - np.log(Px[keep]), [0.05, 0.9], bounds=([1e-6, 0], [10, 0.99999]))
res["yearly"] = {"years": yrs, "ar1_sigma": float(fit.x[0]), "ar1_r": float(fit.x[1]), "lines": {}}
for k, f0 in cand.items():
    i = np.argmin(np.abs(fx - f0)); j = slice(max(i - 1, 0), i + 2); jj = j.start + int(np.argmax(Px[j]))
    b = np.median(Px[(np.abs(fx - f0) > 0.01) & (np.abs(fx - f0) < 0.03)])
    res["yearly"]["lines"][f"east:{k}"] = {"bin_cpd": float(fx[jj]), "bin_offset": int(jj - i), "snr": float(Px[jj] / b),
                                           "over_ar1": float(Px[jj] / G(fit.x, fx[jj]))}
print(f"paper method: {len(yrs)} yearly spectra {yrs[0]}-{yrs[-1]}, AR(1) r = {fit.x[1]:.5f}")
for n in range(1, 5):
    print(f"  n={n} " + "  ".join(f"{x or 'carrier'} {res['yearly']['lines'][f'east:S{n}{x}']['snr']:.1f}"
                                  for x in ("", "-Mf", "+Mf", "-MSf", "+MSf")))
json.dump(res, open(os.path.join(HERE, "wn1_lines.json"), "w"), indent=1)


def mark(a, lines):
    for name, fr, col in lines:
        a.axvline(fr, color=col, lw=0.8, alpha=0.6, ls="--")
        a.text(fr, 1.0, name, rotation=90, va="top", ha="right", fontsize=7, color=col, transform=a.get_xaxis_transform())


def panel(a, f, P, lo, hi, title, lines):
    m = (f > lo) & (f < hi); a.semilogy(f[m], P[m], lw=0.6, color="k"); a.set_xlim(lo, hi)
    mark(a, [l for l in lines if lo < l[1] < hi]); a.set_title(title, fontsize=9.5); a.set_xlabel("eastward frequency (cpd)")


LN = [("S1", S1, "r"), ("K1", K1, "b"), ("P1", P1, "b"), ("O1", O1, "b"), ("OO1", OO1, "b"), ("Q1", Q1, "b"),
      ("S1-Mf", S1 - MF, "orange"), ("S1+Mf", S1 + MF, "orange"), ("S1-MSf", S1 - MSF, "g"), ("S1+MSf", S1 + MSF, "g"),
      ("S2", S2, "r"), ("K2", K2, "b"), ("M2", M2, "b"), ("N2", N2, "b"), ("L2", L2, "b"),
      ("S2-Mf", S2 - MF, "orange"), ("S2+Mf=K2+MSf", S2 + MF, "orange"), ("S2+MSf", S2 + MSF, "g")]
fig, axs = plt.subplots(4, 2, figsize=(13, 16))
a = axs[0, 0]; a.semilogy(fx, Px, lw=0.5, color="k"); a.semilogy(fx, G(fit.x, fx), color="r", lw=1)
a.set_xlim(0, 12); a.set_xlabel("eastward frequency (cpd)"); a.set_ylabel("PSD (hPa$^2$/cpd)")
a.set_title(f"Replica of S&H 2020 Fig. 4: eastward k=1 surface pressure, 20S-20N\n"
            f"ERA5 hourly, {len(yrs)} yearly spectra averaged {yrs[0]}-{yrs[-1]}; red = AR(1) fit", fontsize=9.5)
panel(axs[0, 1], fx, Px, 0.84, 1.16, "Paper method about 1 cpd (bins of 1/365 cpd)", LN)
panel(axs[1, 0], fx, Px, 1.84, 2.16, "Paper method about 2 cpd", LN)
W = 0.0065
for a, ctr, ttl in ((axs[1, 1], 0.9282, "lower 1 cpd sideband: S1-Mf vs O1 = K1-Mf"),
                    (axs[2, 0], 1.0745, "upper 1 cpd sideband: S1+Mf vs OO1 = K1+Mf"),
                    (axs[2, 1], 1.0014, "carrier: S1 vs K1 (and P1)"),
                    (axs[3, 0], 1.9295, "lower 2 cpd sideband: S2-Mf vs M2 = S2-MSf"),
                    (axs[3, 1], 2.0705, "upper 2 cpd sideband: S2+Mf (= K2+MSf) vs S2+MSf")):
    fz = np.arange(ctr - W, ctr + W, df / 5)
    panel(a, fz, dft(t, c, -fz), ctr - W, ctr + W, f"Whole record ({span / 365.25:.0f} yr, 1/T = {df:.1e} cpd): {ttl}", LN)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "wn1_spectrum.png"), dpi=120); print("wrote wn1_spectrum.png")
