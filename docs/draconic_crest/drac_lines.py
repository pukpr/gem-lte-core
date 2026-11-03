import numpy as np
d2r=np.pi/180; eps=23.44*d2r; inc=5.145*d2r
t=np.arange(0,186*365.2422,0.25)              # days, 10 nodal cycles
lam=2*np.pi*t/27.321582                       # tropical (ecliptic longitude)
F=2*np.pi*t/27.212221                         # draconic (argument of latitude)
beta=np.arcsin(np.sin(inc)*np.sin(F))
sind=np.sin(beta)*np.cos(eps)+np.cos(beta)*np.sin(eps)*np.sin(lam)
def amp(x,P):                                 # line amplitude by projection
    w=2*np.pi*t/P; x=x-x.mean()
    return 2*np.hypot((x*np.cos(w)).mean(),(x*np.sin(w)).mean())
Ps={'trop 27.3216':27.321582,'drac 27.2122':27.212221,'Mf 13.6608':13.660791,'Mf-sat 13.6334':1/(1/27.321582+1/27.212221),'18.6 yr':6798.38}
z2=sind**2                                    # degree-2 zonal (k=0) potential ~ P2(sin d)
for name,x in [('sin^2 d  (deg-2 zonal)',z2),('sin d    (linear, odd)',sind),('beta',np.sin(beta)),
               ('Mf-band x sin d (product)',(z2-z2.mean())*sind)]:
    print(f"{name:28s}"+"  ".join(f"{k}:{amp(x,P):.4f}" for k,P in Ps.items()))
# solar equinox-signed semiannual gate: sign of solar declination (alternates at each equinox)
sun=np.sign(np.sin(2*np.pi*t/365.2422))
for name,x in [('sign(sun) x sin d',sun*sind),('sign(sun) x beta',sun*np.sin(beta))]:
    print(f"{name:28s} 2.370yr:{amp(x,365.2422/0.42199):.4f}  2.715yr:{amp(x,365.2422/0.36834):.4f}")
