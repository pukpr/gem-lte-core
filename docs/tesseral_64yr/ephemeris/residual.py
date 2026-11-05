"""What Meeus has that the NINO4 line table does not: residual lines, their aliases, and their share
of the forcing once sampled once a year."""
import json, itertools, numpy as np
exec(open('meeus_lines.py').read().split("Ps = sorted")[0])     # pot(), t, jd, Vp, R, YD, p
Ps = sorted(P for P, _, _ in p['lpap'] if P > 0)
def design(tt, Ps):
    return np.column_stack([np.ones_like(tt)] + [g(2*np.pi*YD/P*tt) for P in Ps for g in (np.sin, np.cos)])
B = design(t, Ps); c = np.linalg.lstsq(B, R, rcond=None)[0]; res = R - B @ c
print(f'rate residual after the fit periods: {100*res.var()/R.var():.3f}% of variance')
# greedy line extraction on the residual (half-day series)
fs = []; r = res.copy(); T = t[-1]-t[0]
grid = np.concatenate([np.arange(0.2, 60, 1/(4*T))])            # cycles per model year, periods >= ~6 d
Y_ = YD
for _ in range(20):
    # coarse search via FFT on zero-padded series
    n = len(r); N = 1 << int(np.ceil(np.log2(8*n)))
    F = np.abs(np.fft.rfft(r*np.hanning(n), N)); fq = np.fft.rfftfreq(N, t[1]-t[0])
    m = (fq > 0.15) & (fq < 75)   # periods 4.9 d .. 6.7 yr
    f0 = fq[m][np.argmax(F[m])]
    fine = np.linspace(f0-1/(4*T), f0+1/(4*T), 201)
    pw = [abs(np.exp(-2j*np.pi*f*t) @ r) for f in fine]; f0 = fine[int(np.argmax(pw))]
    fs.append(f0)
    Bx = np.column_stack([np.ones_like(t)] + [g(2*np.pi*f*t) for f in fs for g in (np.sin, np.cos)])
    cx = np.linalg.lstsq(Bx, res, rcond=None)[0]; r = res - Bx @ cx
amps = np.hypot(cx[1::2], cx[2::2])
iMf = Ps.index(min(Ps, key=lambda q: abs(q-13.6608))); aMf = np.hypot(c[1+2*iMf], c[2+2*iMf])
# Doodson-style label with T(trop), D(drac), A(anom), Y(year)
Tr, Dr, An = 27.321582, 27.212221, 27.554550
def lab(P):
    best = None
    for i, j, k, mm in itertools.product(range(-5, 6), range(-5, 6), range(-5, 6), range(-4, 5)):
        g = i/Tr + j/Dr + k/An + mm/YD
        if g <= 0: continue
        e = abs(1/g - P)
        cst = abs(i)+abs(j)+abs(k)+abs(mm)
        if e < 0.02*P*P/(T*YD) and (best is None or (cst, e) < best[:2]): best = (cst, e, i, j, k, mm)
    return '' if best is None else f'{best[2]:+d}T{best[3]:+d}D{best[4]:+d}A{best[5]:+d}Y'
print(f"\nresidual lines (amplitude relative to Mf in the Meeus rate; gain-weighted = x (13.66/P)):")
print(f"{'P (d)':>9s} {'amp/Mf':>8s} {'x f gain':>8s}  {'combination':20s} alias(yr)")
for f, a in sorted(zip(fs, amps), key=lambda x: -x[1]):
    P = YD/f; al = abs(f - round(f))
    print(f'{P:9.4f} {a/aMf:8.4f} {a/aMf*13.6608/P:8.4f}  {lab(P):20s} {1/al:8.2f}')
# annual sampling: share of residual in the once-a-year forcing, for each calendar month of the impulse
print('\nonce-a-year sampling (impulse month 1..12): residual / full forcing variance')
yrs = np.arange(1880, 2025)
out = []
for mo in range(12):
    ts = yrs + (mo + 0.5)/12
    Rs = np.interp(ts, t, R); rs = np.interp(ts, t, res)
    out.append(rs.var()/Rs.var())
print('  ' + ' '.join(f'{100*v:5.1f}%' for v in out))
