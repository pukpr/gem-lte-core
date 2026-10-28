"""Physical nonlinear mixing vs analysis-window multiplication: where do the S2 x M2 products sit in wavenumber? (2026-10-26)
A quadratic interaction of the migrating tides S2 and M2 (both westward k=2) puts the sum MS4 at westward k=4 and the
difference MSf at k=0; a UTC-locked analysis-window multiplier (the same at all longitudes, k=0) keeps the wavenumber,
putting MS4 at k=2. Fits S2, M2, K2, MS4, S4, M4, 2SM2, MSf, Mf jointly at westward k = 0, 2, 4 of the symmetric
20S-20N surface pressure, 1999-2016, with noise from nearby off-line frequencies.
usage: python3 kmix_test.py"""
import glob, os, numpy as np, netCDF4 as nc
T, X = [], []
for fn in sorted(glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "era5_sp_hourly_*.nc"))):
    if not 1999 <= int(os.path.basename(fn)[15:19]) <= 2016: continue
    ds = nc.Dataset(fn); tv = ds.variables["valid_time"]; T.append(np.asarray(nc.date2num(nc.num2date(tv[:], tv.units), "hours since 1970-01-01")))
    lat = np.asarray(ds.variables["latitude"][:]); lon = np.deg2rad(np.asarray(ds.variables["longitude"][:])); X.append(np.asarray(ds.variables["sp"][:], float).squeeze() / 100); ds.close()
t = np.concatenate(T); o = np.argsort(t); t = t[o]; x = np.concatenate(X)[o]; w = np.cos(np.deg2rad(lat)); p = (x * w[None, :, None]).sum(1) / w.sum()
td = (t - t[0]) / 24; s, h = 1 / 27.321582, 1 / 365.242199; M2 = 2 + 2*h - 2*s; MSF = 2*s - 2*h
L = {"S2": 2, "M2": M2, "K2": 2 + 2*h, "MS4 (S2+M2)": 2 + M2, "S4": 4, "M4 (2M2)": 2 * M2, "2SM2 (2S2-M2)": 4 - M2, "MSf (S2-M2)": MSF, "Mf": 2*s}
fr = np.array(list(L.values()))
print("westward amplitude (hPa) and amp/noise, symmetric 20S-20N, 1999-2016")
for k in (0, 2, 4):
    c = (p * np.exp(-1j * k * lon)[None, :]).mean(1); c -= c.mean()
    sg = -1                                                      # westward: e^{+i w t}
    X_ = np.column_stack([np.exp(-sg * 2j * np.pi * f * td) for f in fr] + [np.ones_like(td)])
    if k == 0:
        X_ = np.column_stack([X_] + [np.exp(2j * np.pi * j * h * td) for j in range(-4, 5) if j])
    a = np.linalg.lstsq(X_, c, rcond=None)[0][:len(fr)]; r = c - np.column_stack([np.exp(-sg * 2j * np.pi * f * td) for f in fr]) @ a
    out = []
    for name, f0, ak in zip(L, fr, a):
        wd = 0.01 if f0 > 0.2 else 0.006
        ff = f0 + np.concatenate([np.linspace(-wd, -wd / 5, 30), np.linspace(wd / 5, wd, 30)]); ff = ff[np.min(np.abs(ff[:, None] - fr[None, :]), 1) > 3 / td[-1]]
        nz = np.sqrt(np.mean([abs(np.mean(r * np.exp(sg * 2j * np.pi * fi * td))) ** 2 for fi in ff]))
        out.append(f"{name} {2*abs(ak):.5f}({abs(ak)/nz:.1f})")
    print(f" k={k}: " + "  ".join(out))
