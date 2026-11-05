"""Hemispherically odd forcing from the same Meeus ephemeris: does it supply the fit's non-zonal lines?
  S1 = (a/r)^3 sin(dec)          (odd response to the lunar potential)
  P3 = (a/r)^4 (5 sin^3 - 3 sin)/2 (degree-3 zonal potential)"""
import json, numpy as np
exec(open('meeus_lines.py').read().split("t = np.arange")[0])
def moon(jd):
    jd = np.atleast_1d(jd); T, L, M, Mm, D, F, A1, A2, A3, eps, E = z.args(jd)
    V = np.column_stack([D, M, Mm, F, A1, A2, A3, L]) * z.D2R
    lam = (L + z.series(z.TL, V, E) - PREC * (T - TP0)) * z.D2R
    bet = z.series(z.TB, V, E) * z.D2R; r = (385000.56 + z.series(z.TR, V, E)) / 385000.56
    return 1/r, np.sin(bet) * np.cos(eps * z.D2R) + np.cos(bet) * np.sin(eps * z.D2R) * np.sin(lam)
t = np.arange(1880, 2025, 0.5 / YD); jd = z.jd_of(t, YD - 365.2422484)
ir, sd = moon(jd)
sers = {'S1': ir**3 * sd, 'P3': ir**4 * (5*sd**3 - 3*sd)/2}
fit = {P: (abs(a), ph + (np.pi if a < 0 else 0)) for P, a, ph in p['lpap'] if P > 0}
Ps = sorted(fit)
m = json.load(open('meeus_lines.json')); Rz = {float(k): v[0] for k, v in m['rate'].items()}
B = np.column_stack([np.ones_like(t)] + [g(2*np.pi*YD/P*t) for P in Ps for g in (np.sin, np.cos)])
Tr = min(Ps, key=lambda q: abs(q-27.3217))
res = {}
for nm, y in sers.items():
    yr = np.gradient(y, t) / YD                      # rate, as for the zonal forcing
    c = np.linalg.lstsq(B, yr - yr.mean(), rcond=None)[0][1:].reshape(-1, 2)
    res[nm] = {P: (np.hypot(*ab), np.arctan2(ab[1], ab[0])) for P, ab in zip(Ps, c)}
    r2 = 1 - ((yr - yr.mean() - B[:, 1:] @ c.ravel())**2).mean() / yr.var()
    print(f'{nm}: the fit periods explain {100*r2:.2f}% of its rate')
aT = fit[Tr][0]
print(f"\n{'P (d)':>9s} {'fit/Trop':>8s} {'S1/Trop':>8s} {'dphi S1':>7s} {'P3/Trop':>8s}   zonal(Meeus) present?")
for P in sorted(Ps, key=lambda q: -res['S1'][q][0]):
    a, ph = fit[P]; s, sp = res['S1'][P]; q = res['P3'][P][0]
    sT = res['S1'][Tr][0]; qT = res['P3'][Tr][0]
    d = np.degrees((ph - sp - (fit[Tr][1] - res['S1'][Tr][1]) + np.pi) % (2*np.pi) - np.pi)
    if s/sT > 0.003 or a/aT > 0.1:
        print(f'{P:9.4f} {a/aT:8.3f} {s/sT:8.4f} {d:+7.0f} {q/qT:8.4f}   {"yes" if Rz[P]/max(Rz.values())>1e-3 else "no"}')
