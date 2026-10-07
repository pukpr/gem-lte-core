#!/usr/bin/env python3
"""Zonal tidal generator v2 (2026-10-06): full Meeus ch. 47 tables taken from lteMod's
gem-ephemeris.adb (62 longitude, 66 latitude, 46 distance terms, with the A1/A2/A3/L corrections),
with fixes: (1) a duplicated latitude argument corrected to 2D-M-M'+F, (2) L converted to radians
inside the L-F, L, L+-M' terms, (3) obliquity T^2 coefficient -1.638889e-7 deg (was the linear
coefficient), (4) eccentricity factor E (E^|m|) applied to terms containing M. Polynomial fundamental
arguments as in gem-ephemeris.adb. Sun: Keplerian. Clock: tropical year from 1880-01-01 00:00 UTC."""
import json, numpy as np
from pathlib import Path
D2R = np.pi / 180; TROP = 365.24219; JD1880 = 2407715.5
T_ = json.load(open(Path(__file__).with_name("meeus_full.json")))
def tab(name):
    return (np.array([t["c"] for t in T_[name]]), np.array([t["mult"] for t in T_[name]], float),
            np.array([t["fn"] == "sin" for t in T_[name]]))
TL, TB, TR = tab("lambda"), tab("B"), tab("R")
T_ANCHOR, JD_ANCHOR = 1990.5, 2448074.761995
def jd_of(t, year=0.0):   # model year = 365.2422484 + YEAR days, anchored at 1990.5 (same as GEM.Zonal)
    return JD_ANCHOR + (np.asarray(t, float) - T_ANCHOR) * (365.2422484 + year)
def args(jd):
    T = (np.asarray(jd, float) - 2451545.0) / 36525.0
    L = 218.3164477 + (481267.88123421 + (-0.0015786 + (1.85584e-6 - 1.53388e-8 * T) * T) * T) * T
    M = 357.5291092 + (35999.0502909 + (-0.0001536 - 4.0833e-8 * T) * T) * T
    Mm = 134.9633964 + (477198.8675055 + (0.0087414 + (1.43741e-6 - 6.7972e-8 * T) * T) * T) * T
    D = 297.8501921 + (445267.1114034 + (-1.8819e-3 + (1.831945e-6 - 8.84447e-9 * T) * T) * T) * T
    F = 93.272095 + (483202.0175233 + (-0.0036539 + (-2.8361e-7 + 1.15833e-9 * T) * T) * T) * T
    A1 = 119.75 + 131.849 * T; A2 = 53.09 + 479264.290 * T; A3 = 313.45 + 481266.484 * T
    eps = 23.43929111 + (-0.0130041667 + (-1.638889e-7 + 5.036111e-7 * T) * T) * T
    E = 1 - (0.002516 + 0.0000074 * T) * T
    return T, L, M, Mm, D, F, A1, A2, A3, eps, E
def series(tb, V, E):          # V: (n, 8) arguments in radians in the order D, M, Mm, F, A1, A2, A3, L
    c, mult, is_sin = tb
    ang = V @ mult.T; ef = E[:, None] ** np.abs(mult[:, 1])[None, :]
    return (ef * c * np.where(is_sin, np.sin(ang), np.cos(ang))).sum(axis=1)
def moon(jd):
    jd = np.atleast_1d(np.asarray(jd, float)); T, L, M, Mm, D, F, A1, A2, A3, eps, E = args(jd)
    V = np.column_stack([D, M, Mm, F, A1, A2, A3, L]) * D2R
    lam = (L + series(TL, V, E)) * D2R; bet = series(TB, V, E) * D2R; r = 385000.56 + series(TR, V, E)
    sdec = np.sin(bet) * np.cos(eps * D2R) + np.cos(bet) * np.sin(eps * D2R) * np.sin(lam)
    return r / 385000.56, sdec
def sun(jd):
    jd = np.atleast_1d(np.asarray(jd, float)); T, L, M, Mm, D, F, A1, A2, A3, eps, E = args(jd)
    L0 = 280.46646 + 36000.76983 * T; Mr = M * D2R
    C = (1.914602 - 0.004817 * T) * np.sin(Mr) + 0.019993 * np.sin(2 * Mr) + 0.000289 * np.sin(3 * Mr)
    e = 0.016708634 - 0.000042037 * T; nu = Mr + C * D2R
    r = (1.000001018 * (1 - e * e)) / (1 + e * np.cos(nu)); sdec = np.sin(eps * D2R) * np.sin((L0 + C) * D2R)
    return r, sdec
def potential(jd, w_sun=0.4600):
    rm, sm = moon(jd); rs, ss = sun(jd)
    return rm ** -3 * (3 * sm ** 2 - 1) / 2 + w_sun * rs ** -3 * (3 * ss ** 2 - 1) / 2
def rate(jd, h=0.5): return (potential(np.asarray(jd) + h) - potential(np.asarray(jd) - h)) / (2 * h)
