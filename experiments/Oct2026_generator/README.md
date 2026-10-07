# Closed-form zonal tidal generator (2026-10-06)

zonal_generator.py: V(t) = sum over Moon, Sun of (abar/r)^3 (3 sin^2 delta - 1)/2, with
Sun/Moon = 0.4600. Moon from the main terms of Meeus ch. 47 (32 longitude/distance and 15 latitude
terms) in D, M, M', F; Sun from the Keplerian orbit. No per-constituent parameters.

Checks against Feb2026/dlod3.dat (daily, 1962-01-03 .. 2019; dLOD is a RATE, so dV/dt is compared):
- CC(dV/dt, dLOD) = -0.969 with only an overall scale, sign and a 1-day alignment; the free
  42-constituent regression (GEM.dLOD) reaches 0.971.
- In lt.exe's 42-constituent basis, ONE complex scale maps the generator onto the LOD-only table:
  major lines agree to 0.98-1.09 in amplitude and within 0.07 rad in phase (Mf, Mf nodal sideband,
  Mt, Mt nodal, Mm, 14.765, 7.096, 9.557, 6.859, 13.777). The nodal-sideband ratios are fixed by
  astronomy: Mf 0.414 (LOD-only 0.409), Mt 0.412 (0.413).
- Exception: 182.623 d (Ssa) is 2.1x the generator in LOD-only, with a -0.92 rad phase: non-tidal
  seasonal LOD (atmosphere), inherited identically by all three climate fits.
Climate fits vs astronomy: the large lines are within about 6-16% and 0.15 rad. The small lines
deviate more, often in the same direction for AMO, PDO and NINO4 (6.859 d about 1.8-2.3x;
13.777 d 1.4-2.2x; 9.557 d 1.1-1.5x; 13.606 d 0.1-0.6x). They share seed ancestry, so the
common direction may be inherited rather than independent.

## v2 (2026-10-06, afternoon): merged with lteMod's GEM.Ephemeris
- `zonal_full.py` (reference) and `src/gem-zonal.adb` (generated from `meeus_full.json`) now use the
  full Meeus ch. 47 tables from `~/eval/lteMod/src/gem-ephemeris.adb`: 62 longitude, 66 latitude and
  46 distance terms, plus the A1/A2/A3/L corrections, with polynomial fundamental arguments.
- Fixes relative to gem-ephemeris.adb: (1) a duplicated latitude argument (2D-M-M'-F appeared twice)
  corrected to 2D-M-M'+F (coefficient 0.002463); (2) L converted to radians inside the L-F, L and
  L+-M' terms (it was added to radian arguments in degrees); (3) obliquity T^2 coefficient corrected
  to -1.638889e-7 deg (the file had the linear coefficient there); (4) eccentricity factor E applied
  to terms containing M. Not used here: lteMod's JULDAT, which uses float instead of integer division.
- Clock: tropical year from 1880-01-01 00:00 UTC, JD = 2407715.5 + (t - 1880) * 365.24219, for dLOD
  and the climate files alike. dLOD prefers it (CC 0.97009 vs 0.96731 for 365.25 and 0.96906 for
  the old fitted YEAR); over 1880-2023 the alternatives drift by 0.18-1.12 days (0.08-0.51 rad of Mf).
- Calibration: ZONAL_KAPPA = -5.866344, ZONAL_TAU = 0.92 d, dLOD CC 0.97014.
- Checks: Ada matches Python to 2e-13 on all 20,829 dLOD dates; TIDES=TABLE reproduces the old fits
  exactly. AMO with nothing refit: 0.606/0.709/0.515 (v1 0.590/0.682/0.410).
