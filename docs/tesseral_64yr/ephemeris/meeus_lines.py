"""Line content of the Meeus zonal potential at the periods of the NINO4 fit, in lt.exe's convention
A sin(2 pi (YD/P) t + phi), t = decimal model year."""
import json, sys, numpy as np
sys.path.insert(0, '../../../experiments/Oct2026_generator'); import zonal_full as z
G = '../../../experiments/Oct2026_generator/legs/nino4_slow_then_harm'
p = json.load(open('nino4_hires_fit.p'))
YD = 365.2422484 + 0.00405305 + p['year']
PREC, TP0 = 1.3969713, (2448074.761995 - 2451545.0) / 36525.0
def pot(jd, w_sun=0.46):
    jd = np.atleast_1d(jd); T, L, M, Mm, D, F, A1, A2, A3, eps, E = z.args(jd)
    V = np.column_stack([D, M, Mm, F, A1, A2, A3, L]) * z.D2R
    lam = (L + z.series(z.TL, V, E) - PREC * (T - TP0)) * z.D2R   # sidereal, as the table / ZONAL default
    bet = z.series(z.TB, V, E) * z.D2R; r = (385000.56 + z.series(z.TR, V, E)) / 385000.56
    sd = np.sin(bet) * np.cos(eps * z.D2R) + np.cos(bet) * np.sin(eps * z.D2R) * np.sin(lam)
    rs, ss = z.sun(jd)
    return r ** -3 * (3 * sd * sd - 1) / 2 + w_sun * rs ** -3 * (3 * ss * ss - 1) / 2
t = np.arange(1880, 2025, 0.5 / YD)          # half-day steps in model years
jd = z.jd_of(t, YD - 365.2422484)
Vp = pot(jd); R = np.gradient(Vp, t) / YD    # per day
Ps = sorted(P for P, _, _ in p['lpap'] if P > 0)
out = {}
for name, y in (('potential', Vp), ('rate', R)):
    B = [np.ones_like(t)]
    for P in Ps:
        w = 2 * np.pi * YD / P * t; B += [np.sin(w), np.cos(w)]
    B = np.column_stack(B); c = np.linalg.lstsq(B, y - y.mean(), rcond=None)[0][1:].reshape(-1, 2)
    out[name] = {P: (float(np.hypot(*ab)), float(np.arctan2(ab[1], ab[0]))) for P, ab in zip(Ps, c)}  # a sin + b cos = A sin(w + phi)
    out[name + '_r2'] = float(1 - ((y - y.mean() - B[:, 1:] @ c.ravel()) ** 2).mean() / y.var())
json.dump({k: ({str(P): v for P, v in d.items()} if isinstance(d, dict) else d) for k, d in out.items()}, open('meeus_lines.json', 'w'), indent=1)
print('fit periods explain', {k: round(v, 4) for k, v in out.items() if k.endswith('_r2')})
