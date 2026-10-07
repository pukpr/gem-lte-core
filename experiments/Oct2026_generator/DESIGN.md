# Reduced-DOF tidal forcing: closed-form zonal generator for dLOD and GEM-LTE

*Design document, 2026-10-06. Status: Python reference implementation done
(`zonal_generator.py`); the Ada integration is specified below and not yet built.*

## 1. Purpose

GEM-LTE currently drives everything with a table of 42 tidal constituents (`lpap` in `lt.exe.p`),
each with a free amplitude and phase: 84 real parameters per fitted series. This document replaces
that table with the **closed-form zonal tidal potential**, computed from a handful of orbital
generating factors. Every constituent, including the cross-harmonic and nodal-sideband terms, then
follows from the expansion instead of being fitted line by line.

Why now (results of 2026-10-06, see `README.md`):

| finding | consequence |
|---|---|
| The generator's time derivative reproduces the daily dLOD record at **CC = -0.969**, with only an overall scale, sign and 1-day alignment. The free 42-line regression reaches 0.971. | Nearly all 84 per-line parameters are redundant for dLOD. |
| One complex scale maps the generator onto the LOD-only constituent table: major lines agree to 0.98-1.09 in amplitude and within 0.07 rad in phase. Nodal-sideband ratios match astronomy (Mf 0.414 vs 0.409; Mt 0.412 vs 0.413). | The empirical LOD table *is* astronomy plus one response factor. |
| Under the annual impulse comb the 42 lines collapse to **28 distinguishable alias frequencies**. Msf is Mf with an annual harmonic (same 3.8-yr alias); Ssa aliases to a constant. | Within an alias group only the combined amplitude and phase can matter. |
| The climate fits depart from astronomy mainly on small lines that sit within about one resolution width of a major alias group (6.859 d next to Mm/Mf; 13.777 d next to the 31.8 d group). AMO/PDO/NINO4 fits are knife-edge sensitive to those departures, and they are not shared across indices. | Per-line tuning is an uncontrolled degree of freedom, likely compensating for neighbouring groups. The generator removes it. |

## 2. Physics

### 2.1 The zonal (long-period) potential
For each body b (Moon, Sun), the degree-2 zonal part of the tide-generating potential is

    V_b(t) = w_b * (abar_b / r_b(t))^3 * (3 sin^2 delta_b(t) - 1) / 2

where r_b is the body's distance, abar_b its mean distance, and delta_b its declination. The
relative weight w_sun / w_moon = (M_sun / a_sun^3) / (M_moon / a_moon^3) = 0.4600. The total is
V(t) = V_moon(t) + V_sun(t).

### 2.2 Generating factors
The Moon's distance and declination come from a short series in the fundamental arguments
(Meeus, *Astronomical Algorithms*, ch. 47; Delaunay arguments):

- D: mean elongation (synodic month)
- M: Sun's mean anomaly (anomalistic year)
- M': Moon's mean anomaly (anomalistic month, 27.555 d)
- F: Moon's argument of latitude (draconic month, 27.212 d). Through F and the mean longitude,
  the node regression (18.61 yr) enters.
- eps: obliquity of the ecliptic

The ecliptic longitude lambda and latitude beta of the Moon, and its distance r, are sums of terms
c * sin or cos(d*D + m*M + m'*M' + f*F). The current implementation uses 32 longitude/distance
terms and 15 latitude terms. Declination follows from
sin(delta) = sin(beta) cos(eps) + cos(beta) sin(eps) sin(lambda).
The Sun is a Keplerian orbit with the equation of the centre.

Expanding (abar/r)^3 (3 sin^2 delta - 1) generates every zonal constituent at once: Mm from the
distance (M'), Mf from the declination squared (2F and 2 x mean longitude), Mt and the higher terms
as products, and the 18.6-yr nodal sidebands from the node's modulation of the declination
amplitude. Their amplitude and phase ratios are therefore fixed by the orbit.

### 2.3 The time derivative (why dLOD is dV/dt)
The solid-Earth response to the zonal tide changes the polar moment of inertia, so the
length-of-day perturbation is proportional to the potential: delta LOD(t) = kappa * V(t). The
calibration series `Feb2026/dlod3.dat` is a **rate** (the code notes "dLOD is already a rate"), so
the model compares with

    dLOD(t) = kappa * dV/dt (t - tau)

where kappa is a scale (with sign) and tau is a timing alignment. Empirically tau = +1 day (dV/dt
leads dLOD by one sample). The derivative is evaluated by a central difference with step
h = 0.5 day. The truncation error is below 1e-4 of the Mf line amplitude, negligible. An analytic
derivative of the series is an option later.

Every constituent's rate amplitude is its potential amplitude times its angular frequency. The
fortnightly and shorter lines therefore dominate the rate. This explains why the 42-line table is
dominated by Mf, its nodal sideband, Mt and Mm, and why the very long-period lines (18.6 yr etc.)
carry negligible amplitude in dLOD.

### 2.4 Angular-momentum sharing (the LOD-wander mechanism)
Per the 2026-10-02 LOD-wander result, the observed LOD is angular momentum shared between the
lunar torque, ocean sloshing and the solid Earth. That gives the design a split by frequency band:

- **Tidal band (days to months):** LOD is dominated by the solid-body response to V. This is the
  band that the generator plus a single scale reproduces (CC 0.969). The generator is calibrated
  **only** here.
- **Decadal band:** the ocean's aliased, integrated response exchanges angular momentum with the
  solid Earth (AMO model vs -LOD, r = +0.70). No fixed astronomical-to-solid-body transfer applies.
  Decadal LOD is therefore a **prediction target** for validation, never an input to the forcing.

The forcing that drives the ocean is the astronomical generator itself, not the LOD-fitted table.
Using LOD-fitted constituents would partly feed the ocean's own momentum signal back in as its
forcing.

## 3. Parameter tiers

The generator's coefficients are empirical (fitted to lunar and solar observations), so they may
be tweaked. They stay *structured*: a tweak changes a generating factor, never a single
constituent, and the expansion carries the change consistently into every line that depends on it.

| tier | parameters | default | freedom | how set |
|---|---|---|---|---|
| 0 Fixed astronomy | fundamental-argument rates (D, M, M', F); obliquity rate | Meeus | none | ephemeris |
| 1 Generating factors (tweakable) | w_sun/moon (0.460); lunar eccentricity term (6.289 deg / 20905 km); evection (1.274 / 3699); variation (0.658 / 2956); annual equation (-0.186); inclination (5.128 deg); obliquity at epoch (23.439 deg) | Meeus | up to 7, bounded (e.g. +-5% of each coefficient; inclination +-0.05 deg) | fit to dLOD in the tidal band only; each tweak logged |
| 1b Clock | timing alignment tau (days); epoch mapping of decimal year to JD | tau = 1.0 d | 1 | fit to dLOD |
| 2 Response before the comb | gain kappa (sign, scale); optional frequency slope s of the gain | from dLOD | 1-2 | kappa from dLOD; s fitted jointly on climate, shared by all series |
| 3 Comb and integration (existing) | impulse timing delB, shape delA/asym; leaky integrator mA, mP; init | current | as now, **shared across all series** | joint climate fit |
| 4 Climate outer (existing) | layer k_f, A, phi_f; Bessel k_b, B, phi_b; regressed windings k_j + harmonics; linear terms | current | per LAYER=1 | layer/Bessel shared; windings per series |

### DOF accounting

| | per series now | with the generator |
|---|---|---|
| tidal constituents | 84 (42 amp + phase) | 0 |
| generator tweaks | - | at most 7 + 1 clock, **shared** by all series and fixed by dLOD |
| response | (absorbed in the 84) | 1-2 shared |
| comb, integrator, layer, Bessel | about 12, nominally shared | about 12, shared |
| outer windings + linear | about 15 | about 15 |

The 89 quads and the climate indices then share everything except their outer terms, which makes
the common-manifold question well posed by construction.

## 4. Implementation in lt.exe

### 4.1 New mode
- Environment / resp key `FORCING=ZONAL` (default `TABLE` = current behaviour, byte-for-byte
  unchanged). Add `FORCING` to the `Options` enum in `src/gem.adb`.
- New package `GEM.Zonal` (spec + body) containing: fundamental arguments as functions of JD;
  the lunar series (coefficients in a constant table, with the tier-1 multipliers applied); the solar
  orbit; `V(jd)`; and `dVdt(jd)` = (V(jd + h) - V(jd - h)) / 2h.
- Time mapping: one function `JD_Of (t : decimal year)`, the single place that defines the clock. (v2: tropical year from 1880-01-01 00:00 UTC, see README.)
  It must reproduce dlod3.dat's convention (a uniform daily grid from 1962-01-03 00:00 UTC, step
  0.0027379 yr) for calibration, and the monthly convention of the climate files for forcing. Both
  conventions are documented in code. TT - UT (about 1 minute) is ignored.

### 4.2 Where it plugs in
`Calc_Forcing` currently evaluates the constituent sum
sum_j A_j cos(2 pi (YL/P_j) t + phi_j) at the model's time steps, before the annual impulse comb
and the leaky integrator. In `ZONAL` mode that sum is replaced by
kappa * (1 + s * g(f)) * dVdt(JD_Of(t) - tau), evaluated analytically at exactly the same time steps.
Here g is the optional gain slope; until it is implemented, s = 0. The impulse comb, integration,
Bessel stage, layer and regression are untouched.

`YEAR` keeps its role for the comb's year length (the impulse recurs once per tropical year), but no
longer re-times the tidal lines. Their clock is now the ephemeris, which removes the "YEAR detuning
re-phases every constituent" artefact found in the wrong-clock test.

### 4.3 Calibration against dLOD
- In `ZONAL` mode `GEM.dLOD` fits only kappa and tau (and optionally the tier-1 tweaks) to
  dlod3.dat, in the tidal band (high-pass at 1/200 d^-1, so decadal LOD stays held out).
- Gate: the dLOD CC stays as now, against the ZONAL calibration (expected about 0.97).
- Ssa (182.6 d) is ignored. It aliases to a constant under the annual comb and is absorbed by the
  level and layer phase, and its LOD amplitude is mostly atmospheric.

### 4.4 Parameter storage
Add a `zonal` object to `lt.exe.p`: kappa, tau, s, and the tier-1 multipliers (default 1.0).
`lpap` is retained but ignored in ZONAL mode, so old fits still load. The windings JSON records the
mode, the tweaks and the calibration CC.

### 4.5 Independent reference
`zonal_generator.py` stays the reference engine. A test compares Ada `dVdt` against Python at 1,000
random dates (agreement to 1e-10 relative), and the Ada calibration against the Python
calibration (kappa, tau, CC). Different language, same equations: this checks the implementation,
not the physics.

## 5. Validation plan (per the acceptance-gate protocols)

Declared before any climate fit in ZONAL mode:

1. **Planted truth (done in Python):** the generator reproduces the LOD-only table (major lines
   within 9% / 0.07 rad; nodal ratios).
2. **dLOD calibration:** ZONAL CC >= 0.965 with kappa and tau only. Then report how much the tier-1
   tweaks add. If more than 0.005, report which factor moved and by how much against its
   ephemeris uncertainty.
3. **Climate refit, low freedom:** AMO, PDO and NINO4 with the shared ZONAL forcing. The
   round-robin v2 criteria (noise floor eps, sharing penalty, success / plateau / budget exits) decide
   whether a common manifold is reached.
4. **Frozen-manifold certification:** level-structure and penalized-winding tests against forcing-side
   nulls. The nulls are a fake Moon (lunar rates scaled by +-1-3%, re-calibrated to dLOD, cold
   start, equal budget) and whole-year clock shifts. Power established on plants first.
5. **Held-out geophysics:** the model's ocean angular-momentum proxy must predict decadal LOD wander
   (the band excluded from calibration), and the 1880-1949 climate back-cast.
6. **Residual diagnosis:** if the climate fits cannot recover their structure, report the residual
   per alias group (28 groups, about 10 with material amplitude). That localizes any missing physics
   to specific aliased periods instead of to 42 free lines.

## 6. Risks and open questions

- **Lunar-theory truncation:** the 32/15-term series is good to about 10 arc-seconds in longitude.
  The effect on the zonal potential is far below the dLOD noise. More terms can be added from the
  same table if residuals point to a specific argument.
- **Climate-file clock:** monthly climate dates (e.g. 1880.0833) need an exact JD mapping. Each day
  of offset moves Mf's raw phase by about 0.46 rad (2 pi / 13.66). Under the comb only the sampled
  instants matter, and those are set by delB, so the mapping must be fixed once and documented.
- **Confirm Calc_Forcing's order of operations** (tidal sum, then impulse comb, then leaky
  integration) when implementing. Section 4.2 assumes that order.
- **Groups near one resolution width:** some alias groups sit within 0.008 cyc/yr of each other
  (Mm 3.9 yr vs 6.859 d at 4.0 yr). The record cannot separate them, so a residual there is
  ambiguous between neighbours by nature.
- **Continuous pathway:** the alias-group degeneracy assumes all forcing reaches the ocean through
  the annual comb. If any continuous pathway exists, Msf and Mf act separately there. ZONAL mode
  handles that automatically, because it generates the full continuous V(t).
