import json, numpy as np
P=[p for p in json.load(open('map_points.json')) if p['err_a']<25 and p['err_p']<25]
lat=np.radians([p['lat'] for p in P]); lon=np.array([p['lon'] for p in P]); A=np.array([p['amp'] for p in P]); ph=np.array([p['phase'] for p in P])
print(len(P),'cells; lat range',np.degrees(lat).min().round(1),np.degrees(lat).max().round(1),'; southern',(lat<0).sum())
shapes={'tesseral |sin 2phi| (P21)':abs(np.sin(2*lat)),'zonal |3sin^2-1|/2 (P20)':abs(3*np.sin(lat)**2-1)/2,
        'sectorial cos^2 (P22)':np.cos(lat)**2,'odd degree-1 |sin phi| (P10)':abs(np.sin(lat)),'constant':np.ones_like(lat)}
print('\namplitude vs latitude, one-parameter fits A*shape:')
for k,s in shapes.items():
    a=(s@A)/(s@s); r=A-a*s; print(f'  {k:30s} A={a:5.1f} mm  rms {np.sqrt((r**2).mean()):5.2f}  R2 {1-r.var()/A.var() if k!="constant" else 0:+.2f}  r {np.corrcoef(s,A)[0,1] if s.std()>0 else 0:+.2f}')
wrap=lambda x:(x+180)%360-180
print('\nphase vs longitude (circular fit, phase = m*lon + c [+180 south]):')
for south_flip in (False,True):
    for m in (-2,-1,0,1,2):
        best=None
        for c in np.arange(-180,180,0.5):
            r=wrap(ph-(m*lon+c+180*(lat<0)*south_flip)); s=np.sqrt((r**2).mean())
            if best is None or s<best[0]: best=(s,c)
        print(f'  m={m:+d} south+180={south_flip!s:5s}  c={best[1]:+6.1f}  circular rms {best[0]:5.1f} deg')
# free slope
best=None
for m in np.arange(-3,3,0.01):
    for c in np.arange(-180,180,2):
        r=wrap(ph-(m*lon+c)); s=np.sqrt((r**2).mean())
        if best is None or s<best[0]: best=(s,m,c)
print(f'  free slope: m={best[1]:+.2f} c={best[2]:+.0f} rms {best[0]:.1f}')
print('\nrandom-phase rms would be ~104 deg')
