#!/usr/bin/env python3
"""Closed-form zonal (long-period) tidal forcing from a few fundamental lunar/solar arguments.

V(t) = sum over Moon, Sun of  w_b * (abar_b / r_b)^3 * (3 sin^2(delta_b) - 1) / 2
with w_sun / w_moon = (M_sun / a_sun^3) / (M_moon / a_moon^3) = 0.4600 (degree-2 tidal ratio).
Moon: main periodic terms of the truncated lunar theory (Meeus, Astronomical Algorithms ch. 47)
in the fundamental arguments D, M, M', F (and the node, through F and the mean longitude).
Sun: Keplerian orbit with the equation of the centre.
Every tidal constituent (Mf, Mm, Mt, their 18.6-yr nodal sidebands, cross terms) then follows from
this closed form; nothing per-constituent is free.
"""
import numpy as np
D2R = np.pi / 180
def fundamentals(jd):
    T = (jd - 2451545.0) / 36525.0
    Lp = 218.3164477 + 481267.88123421 * T            # Moon mean longitude
    D = 297.8501921 + 445267.1114034 * T               # mean elongation
    M = 357.5291092 + 35999.0502909 * T                # Sun mean anomaly
    Mp = 134.9633964 + 477198.8675055 * T              # Moon mean anomaly
    F = 93.2720950 + 483202.0175233 * T                # Moon argument of latitude
    eps = 23.439291 - 0.0130042 * T                    # obliquity
    return T, Lp, D, M, Mp, F, eps
# (D, M, M', F) multipliers with longitude coeff (deg), distance coeff (km)  -- Meeus table 47.A, main terms
LR = [(0,0,1,0,6.288774,-20905.355),(2,0,-1,0,1.274027,-3699.111),(2,0,0,0,0.658314,-2955.968),(0,0,2,0,0.213618,-569.925),
      (0,1,0,0,-0.185116,48.888),(0,0,0,2,-0.114332,-3.149),(2,0,-2,0,0.058793,246.158),(2,-1,-1,0,0.057066,-152.138),
      (2,0,1,0,0.053322,-170.733),(2,-1,0,0,0.045758,-204.586),(0,1,-1,0,-0.040923,-129.620),(1,0,0,0,-0.034720,108.743),
      (0,1,1,0,-0.030383,104.755),(2,0,0,-2,0.015327,10.321),(0,0,1,2,-0.012528,0.0),(0,0,1,-2,0.010980,79.661),
      (4,0,-1,0,0.010675,-34.782),(0,0,3,0,0.010034,-23.210),(4,0,-2,0,0.008548,-21.636),(2,1,-1,0,-0.007888,24.208),
      (2,1,0,0,-0.006766,30.824),(1,0,-1,0,-0.005163,-8.379),(1,1,0,0,0.004987,-16.675),(2,-1,1,0,0.004036,-12.831),
      (2,0,2,0,0.003994,-10.445),(4,0,0,0,0.003861,-11.650),(2,0,-3,0,0.003665,14.403),(0,1,-2,0,-0.002689,-7.003),
      (2,0,-1,2,-0.002602,0.0),(2,-1,-2,0,0.002390,10.056),(1,0,1,0,-0.002348,6.322),(2,-2,0,0,0.002236,-9.884)]
# latitude terms (D, M, M', F, coeff deg) -- Meeus table 47.B, main terms
BT = [(0,0,0,1,5.128122),(0,0,1,1,0.280602),(0,0,1,-1,0.277693),(2,0,0,-1,0.173237),(2,0,-1,1,0.055413),
      (2,0,-1,-1,0.046271),(2,0,0,1,0.032573),(0,0,2,1,0.017198),(2,0,1,-1,0.009266),(0,0,2,-1,0.008822),
      (2,-1,0,-1,0.008216),(2,0,-2,-1,0.004324),(2,0,1,1,0.004200),(2,1,0,-1,-0.003359),(2,-1,-1,1,0.002463)]
def moon(jd, nterms=None):
    T, Lp, D, M, Mp, F, eps = fundamentals(jd); E = 1 - 0.002516 * T
    lam = Lp.copy(); r = np.full_like(jd, 385000.56); beta = np.zeros_like(jd)
    for d, m, mp, f, cl, cr in LR[:nterms]:
        arg = (d * D + m * M + mp * Mp + f * F) * D2R; e = E ** abs(m)
        lam += e * cl * np.sin(arg); r += e * cr * np.cos(arg)
    for d, m, mp, f, cb in BT[:nterms]:
        beta += (E ** abs(m)) * cb * np.sin((d * D + m * M + mp * Mp + f * F) * D2R)
    lam, beta, eps = lam * D2R, beta * D2R, eps * D2R
    sdec = np.sin(beta) * np.cos(eps) + np.cos(beta) * np.sin(eps) * np.sin(lam)
    return r / 385000.56, sdec
def sun(jd):
    T, Lp, D, M, Mp, F, eps = fundamentals(jd)
    L0 = 280.46646 + 36000.76983 * T; Mr = M * D2R
    C = (1.914602 - 0.004817 * T) * np.sin(Mr) + 0.019993 * np.sin(2 * Mr) + 0.000289 * np.sin(3 * Mr)
    e = 0.016708634 - 0.000042037 * T; nu = Mr + C * D2R
    r = (1.000001018 * (1 - e * e)) / (1 + e * np.cos(nu))
    sdec = np.sin(eps * D2R) * np.sin((L0 + C) * D2R)
    return r, sdec
W_SUN = 0.4600
def zonal_potential(jd, w_sun=W_SUN, nterms=None):
    rm, sm = moon(jd, nterms); rs, ss = sun(jd)
    return rm ** -3 * (3 * sm ** 2 - 1) / 2 + w_sun * rs ** -3 * (3 * ss ** 2 - 1) / 2

# --- exact conventions shared with the Ada GEM.Zonal package ---------------------------------
T0_DLOD, JD0_DLOD, STEP_DLOD = 1962.00547580006, 2437667.5, 0.0027379   # dlod3.dat grid: 1962-01-03 00:00 UTC, daily
def jd_of(t): return JD0_DLOD + (np.asarray(t) - T0_DLOD) / STEP_DLOD
def rate(jd, h=0.5, **kw): return (zonal_potential(jd + h, **kw) - zonal_potential(jd - h, **kw)) / (2 * h)   # per day
