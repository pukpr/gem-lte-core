"""Wavenumber-0 (zonal-mean) long-period lunar lines in ERA5 tropical surface pressure (2026-10-26).
At k > 0 a zonally symmetric long-period tide appears only as an offset of a sidereal carrier (O1 = K1 - Mf,
M2 = K2 - Mf = S2 - MSf). At k = 0, the wavenumber of the QBO, the long-period species can appear directly, or as the
k - k = 0 rectification of k > 0 tidal pairs (K1 x O1* -> Mf, K1 x P1* -> 2h). The draconic month enters the degree-2
potential only through 18.6-yr NODAL SATELLITES (Mf' = 2s + n, 13.633 d; 2s + 2n, 13.606 d), which a record longer than
18.6 yr is needed to separate from Mf.

Input   data/era5_sp_hourly_*.nc (fetch_sp.py), years Y0..Y1 (default 1979-2016; gaps are fine for least squares)
Series  cos(lat)-weighted 20S-20N, zonal-mean surface pressure, daily means (hPa)
Fit     joint least squares: quadratic trend, annual harmonics 1..6, and the lines below; noise per line = rms of
        single-frequency fits to the residual at 60 nearby off-line frequencies (0.0015-0.008 cpd away)
Lines   potential lines: Mm (s - p), Mf (2s), Mf' (2s + n), Mf'' (2s + 2n), MSf (2s - 2h), Mtm (3s - p), MSm (s - 2h + p)
        controls (NOT degree-2 potential lines): draconic month (s + n, 27.2122 d), tropical month (s, 27.3216 d)
Output  k0_longperiod_<Y0>_<Y1>.json
usage: python3 k0_longperiod.py [Y0] [Y1]"""
import glob, json, os, sys
import netCDF4 as nc
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
Y0 = int(sys.argv[1]) if len(sys.argv) > 1 else 1979
Y1 = int(sys.argv[2]) if len(sys.argv) > 2 else 2016
S, H, P, NN = 1 / 27.321582, 1 / 365.242199, 1 / 3232.6054, 1 / 6798.383   # moon, sun, perigee, node (regression) rates
LINES = {"Mm": S - P, "Mf": 2 * S, "Mf'(2s+n)": 2 * S + NN, "Mf''(2s+2n)": 2 * S + 2 * NN, "MSf": 2 * S - 2 * H,
         "Mtm": 3 * S - P, "MSm": S - 2 * H + P, "draconic(s+n)": S + NN, "tropical(s)": S}

T, Z = [], []
for fn in sorted(glob.glob(os.path.join(HERE, "data", "era5_sp_hourly_*.nc"))):
    if not Y0 <= int(os.path.basename(fn)[15:19]) <= Y1:
        continue
    ds = nc.Dataset(fn); tv = ds.variables["valid_time" if "valid_time" in ds.variables else "time"]
    hours = np.asarray(nc.date2num(nc.num2date(tv[:], tv.units), "hours since 1970-01-01"))
    lat = np.asarray(ds.variables["latitude"][:]); w = np.cos(np.deg2rad(lat))
    sp = np.asarray(ds.variables["sp"][:], float).squeeze() / 100.0
    Z.append((sp.mean(axis=2) * w[None, :]).sum(axis=1) / w.sum()); T.append(hours); ds.close()
t = np.concatenate(T); z = np.concatenate(Z); o = np.argsort(t); t, z = t[o], z[o]
day = np.floor(t / 24).astype(int); ud, inv, cnt = np.unique(day, return_inverse=True, return_counts=True)
zd = np.bincount(inv, z) / cnt; keep = cnt == 24; td = ud[keep].astype(float); zd = zd[keep]
span = td[-1] - td[0]
print(f"{len(td)} complete days, {1970 + td[0] / 365.2425:.2f}-{1970 + td[-1] / 365.2425:.2f} (span {span / 365.25:.1f} yr), "
      f"1/T = {1 / span:.2e} cpd; nodal spacing n = {NN:.2e} cpd = {NN * span:.1f} Rayleigh bins")

x = (td - td.mean()) / span
base = [np.ones_like(x), x, x ** 2]
for k in range(1, 7):
    base += [np.cos(2 * np.pi * k * H * td), np.sin(2 * np.pi * k * H * td)]
names = list(LINES)
cols = []
for k in names:
    cols += [np.cos(2 * np.pi * LINES[k] * td), np.sin(2 * np.pi * LINES[k] * td)]
X = np.column_stack(base + cols); coef = np.linalg.lstsq(X, zd, rcond=None)[0]; r = zd - X @ coef
nb = len(base); fr = np.array([LINES[k] for k in names])


def noise(f0):
    d = np.concatenate([np.linspace(-0.008, -0.0015, 30), np.linspace(0.0015, 0.008, 30)]); f = f0 + d
    f = f[np.min(np.abs(f[:, None] - fr[None, :]), axis=1) > 3 / span]
    amps = [2 * abs(np.mean(r * np.exp(-2j * np.pi * fi * td))) for fi in f]
    return float(np.sqrt(np.mean(np.square(amps))))


out = {"years": [Y0, Y1], "days": int(len(td)), "span_yr": float(span / 365.25), "lines": {}}
print(f"{'line':14s} {'period d':>9s} {'amp hPa':>9s} {'noise':>8s} {'amp/noise':>9s}")
for i, k in enumerate(names):
    a, b = coef[nb + 2 * i], coef[nb + 2 * i + 1]; amp = float(np.hypot(a, b)); nz = noise(LINES[k])
    out["lines"][k] = {"cpd": LINES[k], "period_d": 1 / LINES[k], "amp_hPa": amp, "noise_hPa": nz, "snr": amp / nz,
                       "phase_deg": float(np.degrees(np.arctan2(-b, a)))}
    print(f"{k:14s} {1 / LINES[k]:9.4f} {amp:9.5f} {nz:8.5f} {amp / nz:9.2f}")
# equilibrium-tide reference amplitude ratios (Cartwright-Tayler): Mf' / Mf = 0.414, Mf'' / Mf = 0.039, MSf / Mf = 0.0861, Mm / Mf = 0.5257
ref = {"Mf'(2s+n)": 0.414, "Mf''(2s+2n)": 0.039, "MSf": 0.0861, "Mm": 0.5257}
L = out["lines"]
print("ratio to Mf (data vs equilibrium potential):", ", ".join(f"{k} {L[k]['amp_hPa'] / L['Mf']['amp_hPa']:.2f} vs {v}" for k, v in ref.items()))
json.dump(out, open(os.path.join(HERE, f"k0_longperiod_{Y0}_{Y1}.json"), "w"), indent=1)
