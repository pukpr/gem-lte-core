"""Observed IERS C04 LOD (tides included): long-period tidal line ratios vs dlod3.dat, Meeus, NINO4 fit,
and the interannual LOD peaks."""
import json, sys, numpy as np
from scipy.ndimage import gaussian_filter1d
sys.path.insert(0, '../../../docs/periodogram'); import peak_estimate as pe
a = np.loadtxt((sys.argv[1] if len(sys.argv) > 1 else 'eopc04.1962-now'), comments='#', usecols=(4, 12))
mjd, lod = a[:, 0], a[:, 1] * 1000.0
t = 1858.87 + (mjd - 0) / 365.2422 - 0.0     # MJD 0 = 1858-11-17
slow = gaussian_filter1d(lod, 60.0); hp = lod - slow
YD = 365.2422
lines = [13.66079, 13.63341, 13.60611, 13.77727, 27.55455, 27.66677, 27.44324, 27.09268, 9.13295, 9.12069,
         9.55689, 7.09581, 6.85940, 14.76530, 31.81209, 5.64271, 27.32158, 27.21222, 13.71877]
def fit(m):
    B = np.column_stack([np.ones(m.sum())] + [g(2*np.pi*YD/P*t[m]) for P in lines for g in (np.sin, np.cos)])
    c = np.linalg.lstsq(B, hp[m], rcond=None)[0]; r = hp[m] - B @ c
    A = {P: np.hypot(c[1+2*i], c[2+2*i]) for i, P in enumerate(lines)}
    # 1-sigma amplitude noise from residual rms
    s = r.std() * np.sqrt(2 / m.sum())
    return A, s
G = '../../../experiments/Oct2026_generator/legs/nino4_slow_then_harm'
d = np.loadtxt('../../../experiments/Feb2026/dlod3.dat')
Bd = np.column_stack([np.ones_like(d[:, 0]), d[:, 0]] + [g(2*np.pi*YD/P*d[:, 0]) for P in lines for g in (np.sin, np.cos)])
cd = np.linalg.lstsq(Bd, d[:, 1], rcond=None)[0]; Ad = {P: np.hypot(cd[2+2*i], cd[3+2*i]) for i, P in enumerate(lines)}
m = json.load(open('meeus_lines.json')); R = {float(k): v[0] for k, v in m['rate'].items()}
p = json.load(open('nino4_hires_fit.p')); F = {P: abs(x) for P, x, ph in p['lpap'] if P > 0}
near = lambda D, P: D[min(D, key=lambda q: abs(q-P))] if min(abs(q-P) for q in D) < 2e-3 else np.nan
all_ = np.ones_like(t, bool); h1 = t < 1994; h2 = ~h1
A, s = fit(all_); A1, _ = fit(h1); A2, _ = fit(h2)
Mf = 13.66079
print(f'C04 1962-{t[-1]:.1f}: Mf amplitude {A[Mf]:.4f} ms, line-amplitude noise ~{s:.4f} ms ({s/A[Mf]:.4f} of Mf)')
print(f"{'P (d)':>9s} {'C04/Mf':>7s} {'1962-93':>7s} {'1994-26':>7s} {'dlod3':>7s} {'Meeus':>7s} {'NINO4':>7s}")
for P in lines:
    print(f'{P:9.4f} {A[P]/A[Mf]:7.4f} {A1[P]/A1[Mf]:7.4f} {A2[P]/A2[Mf]:7.4f} {Ad[P]/Ad[Mf]:7.4f} {near(R, P)/near(R, Mf):7.4f} {near(F, P)/near(F, Mf):7.4f}')
# interannual LOD peaks: monthly means of the slow part, quadratic detrend, 1.5-10 yr
mo = np.floor((t - t[0]) * 12).astype(int); x = np.bincount(mo, slow) / np.bincount(mo)
x = x[:-1]
f, Pk = pe.multitaper(pe.prep(x)); P_ = Pk.mean(0)
pe.BAND = (1/10, 1/1.5)
print('\ninterannual LOD multitaper peaks (yr):', ', '.join(f'{1/v:.2f}' for v in pe.band_peaks(f, P_, 8)))
for o in (60, 120):
    fb, Pb = pe.burg_psd(pe.prep(x), o); print(f'Burg order {o}:', ', '.join(f'{1/v:.2f}' for v in pe.band_peaks(fb, Pb, 8)))
