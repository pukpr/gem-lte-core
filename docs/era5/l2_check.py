"""Migrating lunar semidiurnal tide L2 (westward k=2 M2) in ERA5 hourly surface pressure, 1999-2016 (2026-10-26).
ERA5's forecast model has no lunar gravitational forcing (Sakazaki & Hamilton 2020, sec. 2a), so any L2 can only come
from assimilated observations. Checks: amplitude vs the classical barometric value (~0.07 hPa at the equator),
N2/M2 (perigee / anomalistic-month fingerprint, equilibrium 0.191), and the amplitude fitted separately for each UTC
hour (window dependence; seams at 09|10 and 21|22 UTC).
usage: python3 l2_check.py   (data from data/era5_sp_hourly_*.nc)"""
import glob, os, numpy as np, netCDF4 as nc
T, C2, CL = [], [], []
for fn in sorted(glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "era5_sp_hourly_*.nc"))):
    if not 1999 <= int(os.path.basename(fn)[15:19]) <= 2016: continue
    ds = nc.Dataset(fn); tv = ds.variables["valid_time"]; hrs = np.asarray(nc.date2num(nc.num2date(tv[:], tv.units), "hours since 1970-01-01"))
    lat = np.asarray(ds.variables["latitude"][:]); lon = np.deg2rad(np.asarray(ds.variables["longitude"][:])); sp = np.asarray(ds.variables["sp"][:], float).squeeze() / 100
    w = np.cos(np.deg2rad(lat)); p = (sp * w[None, :, None]).sum(1) / w.sum()
    C2.append((p * np.exp(-2j * lon)[None, :]).mean(1)); CL.append((sp[:, lat == 0, :].squeeze() * np.exp(-2j * lon)[None, :]).mean(1)); T.append(hrs); ds.close()
t = np.concatenate(T); o = np.argsort(t); t = t[o]; td = (t - t[0]) / 24; hr = t.astype(int) % 24
s, h = 1 / 27.321582, 1 / 365.242199; M2 = 2 - 2 * s + 2 * h
for name, c in (("20S-20N mean", np.concatenate(C2)[o]), ("equator 0N", np.concatenate(CL)[o])):
    c = c - c.mean()
    fr = [M2, 2, 2 + 2 * h, 2 - 2 * h, 2 + h, 2 - h, 2 - 3 * s + 2 * h + 1 / 3232.6054]
    X = np.exp(2j * np.pi * np.outer(td, fr)); a = np.linalg.lstsq(X, c, rcond=None)[0]   # westward k=2: e^{+i w t}
    amp = 2 * abs(a[0])                                                                   # wave amplitude in hPa
    ph = []
    per = []
    for hh in range(24):
        m = hr == hh; Xh = np.exp(2j * np.pi * np.outer(td[m], fr)); ah = np.linalg.lstsq(Xh, c[m], rcond=None)[0][0]; per.append(2 * abs(ah)); ph.append(np.degrees(np.angle(ah / a[0])))
    per = np.array(per)
    print(f"{name}: westward k=2 M2 (L2) amplitude {amp:.4f} hPa; S2 {2*abs(a[1]):.3f} hPa; N2 {2*abs(a[6]):.4f} (N2/M2 {abs(a[6]/a[0]):.3f}, eq 0.191)")
    print(f"   per-UTC-hour amplitude: min {per.min():.4f} max {per.max():.4f} (max/min {per.max()/per.min():.2f}); phase spread {np.ptp(ph):.0f} deg")
    print("   ", " ".join(f"{x:.3f}" for x in per))
