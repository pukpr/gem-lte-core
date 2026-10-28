"""Lunar declination (draconic) content in the ANTISYMMETRIC k=1 diurnal band, ERA5 1979-2016 (2026-10-26).
Degree-2 potential: the k=0 term is ~ delta^2 and symmetric about the equator; the k=1 diurnal term is ~ sin(2 delta) and
ANTISYMMETRIC, so the declination envelope delta(t), which contains the draconic month through the 5.15 deg inclination
(18.6-yr nodal satellites, O1' = O1 - n with equilibrium ratio 0.189 ~ i/eps), shows up linearly only there.
Series  k=1 coefficient of p(0-20N) - p(20S-0), cos(lat)-weighted, hourly 1979-2016 (gaps fine for least squares)
        plus the k=0 antisymmetric daily mean (degree-3 draconic-fundamental control)
Fit     joint LS of the listed lines, east and west separately; noise = rms single-frequency fit to the residual at
        60 nearby off-line frequencies.  Controls S1-Mf, P1-Mf, OO1, M1 are not (or barely) potential lines.
usage: python3 antisym_lunar.py   (prints; data from data/era5_sp_hourly_*.nc)"""
import glob, os, numpy as np, netCDF4 as nc
T, CA, CS, ZA = [], [], [], []
for fn in sorted(glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "era5_sp_hourly_*.nc"))):
    if not 1979 <= int(os.path.basename(fn)[15:19]) <= 2016: continue
    ds = nc.Dataset(fn); tv = ds.variables["valid_time"]; hrs = np.asarray(nc.date2num(nc.num2date(tv[:], tv.units), "hours since 1970-01-01"))
    lat = np.asarray(ds.variables["latitude"][:]); lon = np.deg2rad(np.asarray(ds.variables["longitude"][:])); sp = np.asarray(ds.variables["sp"][:], float).squeeze() / 100
    w = np.cos(np.deg2rad(lat)); n, s_ = lat > 0, lat < 0
    pn = (sp[:, n] * w[n][None, :, None]).sum(1) / w[n].sum(); ps = (sp[:, s_] * w[s_][None, :, None]).sum(1) / w[s_].sum()
    a = pn - ps; e = np.exp(-1j * lon)[None, :]
    CA.append((a * e).mean(1)); ZA.append(a.mean(1)); T.append(hrs); ds.close()
t = np.concatenate(T); o = np.argsort(t); t = t[o]; ca = np.concatenate(CA)[o]; za = np.concatenate(ZA)[o]
td = (t - t[0]) / 24; span = td[-1]; ca -= ca.mean()
s, h, p, n = 1 / 27.321582, 1 / 365.242199, 1 / 3232.6054, 1 / 6798.383
L = {"S1": 1, "P1": 1 - h, "K1": 1 + h, "S1-2h": 1 - 2*h, "S1+2h": 1 + 2*h, "O1": 1 + h - 2*s, "O1'(-n)": 1 + h - 2*s - n,
     "K1'(+n)": 1 + h + n, "Q1": 1 + h - 3*s + p, "OO1": 1 + h + 2*s, "M1(s+n?)": 1 + h - s, "S1-Mf": 1 - 2*s, "P1-Mf": 1 - h - 2*s,
     "S2": 2, "M2": 2 + 2*h - 2*s, "K2": 2 + 2*h}
names = list(L); fr = np.array([L[k] for k in names])
print(f"antisymmetric (0-20N minus 20S-0) k=1, {span/365.25:.1f} yr")
for direction, sg in (("west", -1), ("east", 1)):
    X = np.exp(-sg * 2j * np.pi * np.outer(td, fr)); a = np.linalg.lstsq(X, ca, rcond=None)[0]; r = ca - X @ a
    def noise(f0):
        f = f0 + np.concatenate([np.linspace(-0.01, -0.002, 30), np.linspace(0.002, 0.01, 30)]); f = f[np.min(np.abs(f[:, None] - fr[None, :]), axis=1) > 3 / span]
        return 2 * np.sqrt(np.mean([abs(np.mean(r * np.exp(sg * 2j * np.pi * fi * td))) ** 2 for fi in f]))
    A = {k: 2 * abs(x) for k, x in zip(names, a)}
    print(f"--- {direction}: " + "  ".join(f"{k} {A[k]:.4f}({A[k] / noise(L[k]):.1f})" for k in names))
    print(f"    O1'/O1 = {A["O1'(-n)"] / A['O1']:.3f} (eq 0.189); Q1/O1 = {A['Q1'] / A['O1']:.3f} (eq 0.194); K1/O1 = {A['K1'] / A['O1']:.2f} (lunar-only eq 1.41)")
# k=0 antisymmetric daily: degree-3 draconic fundamental control
day = np.floor(t / 24).astype(int); ud, inv, cnt = np.unique(day, return_inverse=True, return_counts=True); zd = (np.bincount(inv, za) / cnt)[cnt == 24]; dd = ud[cnt == 24].astype(float)
x = (dd - dd.mean()) / (dd[-1] - dd[0]); B = [np.ones_like(x), x, x**2] + [f(2*np.pi*k*h*dd) for k in range(1, 7) for f in (np.cos, np.sin)]
LL = {"draconic(s+n)": s + n, "tropical(s)": s, "Mm": s - p, "Mf": 2*s}
C = [f(2*np.pi*v*dd) for v in LL.values() for f in (np.cos, np.sin)]; X = np.column_stack(B + C); cf = np.linalg.lstsq(X, zd, rcond=None)[0]; r = zd - X @ cf
for i, (k, v) in enumerate(LL.items()):
    amp = np.hypot(*cf[len(B) + 2*i: len(B) + 2*i + 2]); f = v + np.concatenate([np.linspace(-0.008, -0.0015, 30), np.linspace(0.0015, 0.008, 30)])
    nz = np.sqrt(np.mean([(2 * abs(np.mean(r * np.exp(-2j * np.pi * fi * dd)))) ** 2 for fi in f]))
    print(f"k=0 antisym daily {k:14s} amp {amp:.4f} hPa  {amp / nz:.2f} x noise")
