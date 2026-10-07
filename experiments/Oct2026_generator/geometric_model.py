#!/usr/bin/env python3
"""Concise geometric zonal-tide model (2026-10-06): the Moon as a Keplerian ellipse (eccentricity e)
inclined at i to the ecliptic, plus (optionally) the three main solar perturbations of the orbit
(evection, variation, annual equation); the Sun as a Keplerian ellipse. Cross-harmonics arise only from
(a/r)^3 (3 sin^2 dec - 1)/2. Rates: the fundamental arguments (linear). Constants:
  e = 0.0549, i = 5.145 deg, obliquity 23.44 deg, Sun/Moon weight 0.46, Sun e = 0.0167,
  evection 1.274 deg / -0.0096 a (= 1.274 deg x the longitude-to-distance ratio of the equation of centre),
  variation 0.658 deg / -0.0077 a, annual equation -0.186 deg."""
import numpy as np
D2R = np.pi / 180
def fund(jd):
    T = (jd - 2451545.0) / 36525.0
    return (218.3164477 + 481267.88123421 * T, 297.8501921 + 445267.1114034 * T, 357.5291092 + 35999.0502909 * T,
            134.9633964 + 477198.8675055 * T, 93.2720950 + 483202.0175233 * T, 280.46646 + 36000.76983 * T)
def potential(jd, e=0.0549, inc=5.145, ev=1.274, var=0.658, ann=-0.186, ev_r=-0.0096, var_r=-0.0077, eps=23.44, w_sun=0.46, es=0.0167):
    L, D, M, Mp, F, Ls = [x * D2R for x in fund(jd)]
    lam = L + 2 * e * np.sin(Mp) + 1.25 * e * e * np.sin(2 * Mp) + (ev * np.sin(2 * D - Mp) + var * np.sin(2 * D) + ann * np.sin(M)) * D2R
    r = 1 - e * np.cos(Mp) + 0.5 * e * e * (1 - np.cos(2 * Mp)) + ev_r * np.cos(2 * D - Mp) + var_r * np.cos(2 * D)
    bet = inc * D2R * np.sin(F + 2 * e * np.sin(Mp))           # latitude follows the true argument of latitude
    sdec = np.sin(bet) * np.cos(eps * D2R) + np.cos(bet) * np.sin(eps * D2R) * np.sin(lam)
    lam_s = Ls + 2 * es * np.sin(M); rs = 1 - es * np.cos(M); sds = np.sin(eps * D2R) * np.sin(lam_s)
    return r ** -3 * (3 * sdec ** 2 - 1) / 2 + w_sun * rs ** -3 * (3 * sds ** 2 - 1) / 2
