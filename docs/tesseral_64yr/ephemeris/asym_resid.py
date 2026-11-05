import itertools, numpy as np
exec(open('asym.py').read().split("res = {}")[0])
y = np.gradient(sers['S1'], t) / YD; y = y - y.mean()
c = np.linalg.lstsq(B, y, rcond=None)[0]; base = B @ c; res = y - base
iT = 1 + 2*Ps.index(Tr); aT = np.hypot(c[iT], c[iT+1])
print(f'S1 rate residual after the fit periods: {100*res.var()/y.var():.2f}%')
Ttot = t[-1]-t[0]; fs = []; r = res.copy()
for _ in range(14):
    n = len(r); N = 1 << int(np.ceil(np.log2(8*n)))
    F = np.abs(np.fft.rfft(r*np.hanning(n), N)); fq = np.fft.rfftfreq(N, t[1]-t[0])
    m = (fq > 0.05) & (fq < 75); f0 = fq[m][np.argmax(F[m])]
    fine = np.linspace(f0-1/(4*Ttot), f0+1/(4*Ttot), 201)
    f0 = fine[int(np.argmax([abs(np.exp(-2j*np.pi*f*t) @ r) for f in fine]))]; fs.append(f0)
    Bx = np.column_stack([np.ones_like(t)] + [g(2*np.pi*f*t) for f in fs for g in (np.sin, np.cos)])
    cx = np.linalg.lstsq(Bx, res, rcond=None)[0]; r = res - Bx @ cx
amps = np.hypot(cx[1::2], cx[2::2])
Tr_, Dr, An = 27.321582, 27.212221, 27.554550
def lab(P):
    best = None
    for i, j, k, mm in itertools.product(range(-4, 5), range(-4, 5), range(-4, 5), range(-4, 5)):
        g = i/Tr_ + j/Dr + k/An + mm/YD
        if g <= 0: continue
        e = abs(1/g - P); cst = abs(i)+abs(j)+abs(k)+abs(mm)
        if e < 0.05*P*P/(Ttot*YD) and (best is None or (cst, e) < best[:2]): best = (cst, e, i, j, k, mm)
    return '' if best is None else f'{best[2]:+d}T{best[3]:+d}D{best[4]:+d}A{best[5]:+d}Y'
print(f"{'P (d)':>9s} {'amp/Trop':>8s}  {'combination':18s} alias(yr)")
for f, a in sorted(zip(fs, amps), key=lambda x: -x[1]):
    P = YD/f; al = abs(f - round(f))
    print(f'{P:9.4f} {a/aT:8.4f}  {lab(P):18s} {1/al:8.2f}')
