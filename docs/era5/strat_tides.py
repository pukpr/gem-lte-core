"""Tidal lines from the surface into the QBO levels: ERA5 hourly, 20S-20N, same analysis for every field (2026-10-26).
Question: how do the lunisolar tidal factors change in importance between surface pressure and the stratospheric zonal
(and meridional) wind, in particular the Mf / MSf satellite lines and the k = 0 (zonal-mean) long-period lines that a
QBO-type forcing would need?

Fields   sp (surface pressure, hPa; fetch_sp.py) and u, v at 70, 50, 30, 10 hPa (m/s; fetch_uv.py), each 5 deg, hourly
Indices  sym  = cos(lat)-weighted 20S-20N mean;  anti = mean(0-20N) - mean(20S-0)
Waves    k = 0 (zonal mean), and k = 1, 2 split into eastward / westward (c_k = <x e^{-ik lon}>; eastward at numpy f < 0)
Lines    joint least squares on the hourly series (trend + annual harmonics 1..4 included for k = 0):
           k = 0 long period: Mm, Mf, Mf' (2s+n), MSf, Mtm, MSm; controls draconic (s+n), tropical (s)
           k = 0 and k >= 1 tidal: S1, P1, K1, K1', O1, O1', Q1, S2, K2, T2, R2, M2, M2', N2, MS4, and S_n +/- Mf, +/- MSf (n=1..4)
         noise = rms single-frequency fit to the residual at 60 nearby off-line frequencies (window scales with period)
Importance (unit-free, comparable across fields): amplitude / noise, and lunar / solar ratios (O1/K1, M2/S2, M2/S2-family,
         sideband / carrier). Seam check per series: step excess at 09|10 and 21|22 UTC (eastward part not separated).
Output   strat_tides_<Y0>_<Y1>.json (+ printed tables)
usage: python3 strat_tides.py [Y0] [Y1] [fields, e.g. sp,u30,u50]   (missing fields are skipped)"""
import glob, json, os, sys
import netCDF4 as nc
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
Y0 = int(sys.argv[1]) if len(sys.argv) > 1 else 1999
Y1 = int(sys.argv[2]) if len(sys.argv) > 2 else 2016
FIELDS = sys.argv[3].split(",") if len(sys.argv) > 3 else ["sp", "u70", "u50", "u30", "u10", "v70", "v50", "v30", "v10"]
s, h, p, n = 1 / 27.321582, 1 / 365.242199, 1 / 3232.6054, 1 / 6798.383       # moon, sun, perigee, node rates (cpd)
MF, MSF = 2 * s, 2 * s - 2 * h
TIDAL = {"S1": 1, "P1": 1 - h, "K1": 1 + h, "K1'": 1 + h + n, "O1": 1 + h - 2 * s, "O1'": 1 + h - 2 * s - n, "Q1": 1 + h - 3 * s + p,
         "S2": 2, "K2": 2 + 2 * h, "T2": 2 - h, "R2": 2 + h, "M2": 2 + 2 * h - 2 * s, "M2'": 2 + 2 * h - 2 * s - n,
         "N2": 2 + 2 * h - 3 * s + p, "MS4": 4 + 2 * h - 2 * s}
for k in range(1, 5):
    TIDAL.update({f"S{k}-Mf": k - MF, f"S{k}+Mf": k + MF, f"S{k}-MSf": k - MSF, f"S{k}+MSf": k + MSF})
for k in (3, 4):
    TIDAL[f"S{k}"] = float(k)
ALIAS = {}                                               # exact duplicates: keep the named tidal line (M2 = S2-MSf, MS4 = S4-MSf)
for q in list(TIDAL):
    for r in list(TIDAL):
        if r == q:
            break
        if abs(TIDAL[q] - TIDAL[r]) < 1e-9:
            ALIAS[q] = r; del TIDAL[q]; break
LONG = {"Mm": s - p, "Mf": MF, "Mf'": MF + n, "MSf": MSF, "Mtm": 3 * s - p, "MSm": s - 2 * h + p,
        "draconic(ctl)": s + n, "tropical(ctl)": s}


def load(field):
    var, lev = (field, None) if field == "sp" else (field[0], field[1:])
    pat = f"era5_sp_hourly_*.nc" if field == "sp" else f"era5_{var}{lev}_hourly_*.nc"
    T, X = [], []
    for fn in sorted(glob.glob(os.path.join(HERE, "data", pat))):
        yr = int(os.path.basename(fn).split("_hourly_")[1][:4])
        if not Y0 <= yr <= Y1:
            continue
        ds = nc.Dataset(fn); tv = ds.variables["valid_time" if "valid_time" in ds.variables else "time"]
        T.append(np.asarray(nc.date2num(nc.num2date(tv[:], tv.units), "hours since 1970-01-01")))
        x = np.asarray(ds.variables[var][:], float).squeeze()
        X.append(x / 100.0 if field == "sp" else x); lat = np.asarray(ds.variables["latitude"][:]); lon = np.asarray(ds.variables["longitude"][:]); ds.close()
    if not T:
        return None
    t = np.concatenate(T); x = np.concatenate(X); o = np.argsort(t); t, x = t[o], x[o]
    keep = np.r_[True, np.diff(t) > 0]; t, x = t[keep], x[keep]
    w = np.cos(np.deg2rad(lat)); nh, sh = lat > 0, lat < 0
    sym = (x * w[None, :, None]).sum(1) / w.sum()
    anti = (x[:, nh] * w[nh][None, :, None]).sum(1) / w[nh].sum() - (x[:, sh] * w[sh][None, :, None]).sum(1) / w[sh].sum()
    return t, np.deg2rad(lon), {"sym": sym, "anti": anti}


def fit(td, y, lines, sign, extra_daily=False):
    """joint LS of complex lines e^{sign * -i 2 pi f t} (sign=+1 eastward / k=0 complex, -1 westward); returns amps, noise"""
    names = list(lines); fr = np.array([lines[k] for k in names])
    cols = [np.exp(-sign * 2j * np.pi * f * td) for f in fr]
    base = [np.ones_like(td), (td - td.mean()) / np.ptp(td)]
    if extra_daily:                                      # k = 0 real series: seasonal harmonics as nuisance
        base += [np.exp(-sign * 2j * np.pi * k * h * td) for k in range(-4, 5) if k]
    X = np.column_stack(base + cols); a = np.linalg.lstsq(X, y, rcond=None)[0]; r = y - X @ a; a = a[len(base):]
    span = np.ptp(td); out = {}
    for k, f0, ak in zip(names, fr, a):
        wdt = max(0.01, 0.15 * f0) if f0 < 0.2 else 0.01
        f = f0 + np.concatenate([np.linspace(-wdt, -wdt / 5, 30), np.linspace(wdt / 5, wdt, 30)])
        f = f[np.min(np.abs(f[:, None] - fr[None, :]), axis=1) > 3 / span]
        nz = np.sqrt(np.mean([abs(np.mean(r * np.exp(sign * 2j * np.pi * fi * td))) ** 2 for fi in f]))
        mult = 2.0 if extra_daily else 2.0                # amplitude of the real wave = 2|a| in both conventions here
        out[k] = {"amp": float(mult * abs(ak)), "noise": float(mult * nz), "snr": float(abs(ak) / nz)}
    return out


def seam(t, y):
    hr = t.astype(int) % 24; ok = np.r_[np.diff(t) == 1, False]
    day = np.array([y[hr == q].mean() for q in range(24)]); a = y - day[hr]; d = np.diff(a)
    good = ok[:-1] & np.r_[ok[1:-1], False]
    e = np.abs(d[1:-1] - 0.5 * (d[:-2] + d[2:])); into = hr[2:-1]; g = good[1:]
    J = np.array([e[(into == q) & g[:len(e)]].mean() for q in range(24)]); return float(J[10] / np.median(J)), float(J[22] / np.median(J))


res = {"years": [Y0, Y1], "aliases": ALIAS, "fields": {}}
print("exact aliases dropped:", ALIAS)
for field in FIELDS:
    L = load(field)
    if L is None:
        print(f"{field}: no files for {Y0}-{Y1}, skipped"); continue
    t, lon, idx = L; td = (t - t[0]) / 24.0; R = res["fields"][field] = {"hours": int(len(t)), "span_yr": float(np.ptp(td) / 365.25)}
    print(f"===== {field}: {len(t)} hours, {np.ptp(td) / 365.25:.1f} yr")
    for ix, x in idx.items():
        for k in (0, 1, 2):
            c = (x * np.exp(-1j * k * lon)[None, :]).mean(1); c = c - c.mean()
            if k == 0:
                y = c.real; daily = y.reshape(-1)          # zonal mean: real series, lines as complex pairs below
                rows = fit(td, y.astype(complex), {**LONG, **{q: v for q, v in TIDAL.items() if v < 2.5}}, +1, extra_daily=True)
                R[f"{ix}_k0"] = rows; R[f"{ix}_k0_seam"] = seam(t, y)
            else:
                for direc, sg in (("east", +1), ("west", -1)):
                    R[f"{ix}_k{k}_{direc}"] = fit(td, c, TIDAL, sg)
        # compact print
        r0 = R[f"{ix}_k0"]
        print(f"  {ix} k=0 long-period SNR: " + "  ".join(f"{q} {r0[q]['snr']:.1f}" for q in LONG) + f"   seam {R[f'{ix}_k0_seam'][0]:.1f}x/{R[f'{ix}_k0_seam'][1]:.1f}x")
        print(f"  {ix} k=0 tidal SNR:       " + "  ".join(f"{q} {r0[q]['snr']:.1f}" for q in ("S1", "K1", "O1", "O1'", "S2", "M2", "S1-Mf", "S1-MSf", "S2+MSf")))
        for k in (1, 2):
            for direc in ("east", "west"):
                r = R[f"{ix}_k{k}_{direc}"]
                rat = lambda a, b: r[a]["amp"] / max(r[b]["amp"], 1e-30)
                print(f"  {ix} k={k} {direc}: SNR " + " ".join(f"{q} {r[q]['snr']:.1f}" for q in ("S1", "K1", "O1", "O1'", "S2", "M2", "N2", "MS4", "S2+MSf", "S4+MSf", "S2-Mf", "S2+Mf"))
                      + f" | O1/K1 {rat('O1', 'K1'):.3f} M2/S2 {rat('M2', 'S2'):.3f}")
json.dump(res, open(os.path.join(HERE, f"strat_tides_{Y0}_{Y1}.json"), "w"), indent=1)
print(f"wrote strat_tides_{Y0}_{Y1}.json")
