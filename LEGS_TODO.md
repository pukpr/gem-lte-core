# Leg groups: to do

Companion to `LEGS.md`. First written 2026-10-12; reviewed 2026-10-13 after the multi-index runs, with
the 5-day against monthly SOI question back in view.

## Where the legs stand

| Done | Evidence |
|---|---|
| N legs (`LEGS N`, 2 to 9), per-leg response files | 2, 3 and 4-leg runs; each leg file bit-identical to its single run |
| Central choice of the saving group; `SELECT_TOP` tie-break on held-out score | every run since AMO+PDO |
| Per-leg settings beyond the index file: `SAMPLING`, intervals, `LEG_IR`, `ALPHA_GAIN`, `UNCOMPENSATED` | SOI 5-day + monthly; AMO+NAO; four Nino indices |
| Test mode exits with an idle pipe on input | fixed |

What the multi-index runs established:

- **One manifold carries AMO, PDO, NINO4 and IODE** (+0.49 to +0.73 with trends removed), and separately
  **NINO4, NINO34, NINO3, NINO12** (+0.56 to +0.70). The ninth harmonic of the 0.2017 base leads in every
  Pacific index.
- **For different indices, separate amplitudes and phases are right**, and comparing them is informative
  (PDO and NINO4 agree at +0.86; the four Nino indices at +0.83 to +0.98; AMO stands apart).
- **AMO alone does not pin the manifold** to better than about 0.05 cycles of the base; NAO and Baltic
  each pair with it on manifolds that far apart.
- **Timestamps matter:** dates rounded to four decimals (0.012 days) move the manifold by 0.06 cycles.
  Index files must carry exact sample dates.
- **dLOD agreement drifts down slowly** with long free searches (0.9991 to 0.9982). Changes so far are
  within dLOD's own per-line errors except one line; `legs/dlod_police.py` checks this.

## The SOI question, restated

Two legs on the same index (5-day SOI calibrated on 1992-2005, monthly SOI on the rest of 1880-2026)
are a different problem from two indices. There the shared windings are not enough: each leg regresses
its own amplitudes and phases, and the 5-day + monthly result agreed between legs at only +0.33
overall and -0.07 for harmonics above 30. So that run does not show that the fast excursions seen in
1992-2005 exist, relative to the manifold, anywhere else.

The aim: **pin down how far the fast excursions are stationary relative to the manifold over the whole
span, given that they are calibrated on 1992-2005.**

### First measurement (2026-10-13, no code change)

`legs/soi_stationarity/`: the 5-day fit is evaluated on the 5-day grid from 1880 (`pad_run.py`), split
into LOW (harmonics up to 30) and HIGH (81, 132, 134, 159, 215, 227), averaged to months and compared
with the monthly SOI era by era. HIGH is tested against the SOI minus LOW; "gain" is the slope of that
target on HIGH (1 = the calibrated amplitude holds, 0 = gone); z is against random high harmonics
calibrated the same way.

| Era | HIGH: CC | Gain | z |
|---|---|---|---|
| 1880-1910 | +0.04 | 0.17 | +0.4 |
| 1910-1940 | -0.04 | -0.16 | -0.5 |
| 1940-1970 | +0.08 | 0.35 | +2.0 |
| 1970-1992 | -0.07 | -0.35 | -0.9 |
| 1992-2005 (calibration) | +0.26 | 1.09 | +2.1 |
| 2005-2026 | +0.09 | 0.36 | +1.2 |
| All outside 1992-2005 | +0.01 | 0.05 | +0.6 |

At 5-day resolution after the window: gain 0.22 in 2005-2012, 0.51 in 2012-2019 (z +2.4), 0.00 in
2019-2027. (`stationarity2.py`, LOW taken from the long monthly record; `stationarity.py`, with LOW also
calibrated on the 13 years, gives the same picture, and there LOW itself fails outside the window.)

Reading: on this manifold the fast response calibrated in 1992-2005 keeps a fraction of its amplitude
for a decade or so afterwards and none in the earlier record. This does not yet separate two causes:

- the fast response is not stationary; or
- it is, but the manifold is not accurate enough away from the window. The top harmonic turns 46 times
  per unit of manifold, so a level error of 0.01 is half a cycle, and the multi-index runs show
  equally good manifolds 0.05 cycles of the base apart.

### Planted-term control (2026-10-13)

`legs/soi_stationarity/planted.py`, `planted_dlod.py`: a known fast term (the six high harmonics with their
1992-2005 amplitudes) applied unchanged over 1880-2026 on the fit's own manifold, plus the slow part and
noise with the real residual's spectrum and size; monthly data = month means. Same measurement.

| Case | Gain outside 1992-2005 (1 = recovered) |
|---|---|
| Planted, exact manifold | 0.86 +- 0.16 (z +7 to +9) |
| Nothing planted | -0.01 +- 0.15 |
| Real data | 0.05 |
| Planted; yearly level error random-walking from 1998 at 0.0005 / 0.002 / 0.005 per year | 0.84 / 0.42 / 0.01 |
| Planted; tide amplitudes moved by 1% / 3% / 10% / 100% of dLOD's per-line error | 0.3 to 0.9 / 0 / 0 / 0 |

- The measurement works: a stationary fast term of this size would be recovered at +7 sd through real
  noise and month averaging.
- It needs the manifold right to about 0.2 cycles of the top harmonic, i.e. yearly levels good to about
  0.0005 manifold units against yearly steps of about 6.
- A table known only to dLOD's accuracy is short of that by about a factor of 100. The manifold's
  working precision comes from the climate fit itself, and for the fast harmonics only the 13
  calibration years supply it.
- In the window the fast part "fits" whatever the manifold (gain 1.2 to 1.5 with a wrong manifold, 1.3
  with nothing planted), so the in-window fit is not evidence.
- So the real data's null outside the window does not show non-stationarity. It cannot be decided by
  transferring a fixed manifold; the search has to be asked for a manifold on which the fast term holds
  everywhere (items 1 to 3 below), and the control for THAT is to plant a term, perturb the start and
  see whether the search finds its way back.

### What is needed to separate them

*Built 2026-10-13:* `LEG_COMMON`, `LEG_SCALE`, `LEG_SCORE FAST` (see `LEGS.md`): the fast terms come from
the calibration leg and are scored on the other legs. The slow terms are still regressed per leg, not
jointly. First runs: `legs/soi_stationarity/common_planted` (recovery control) and `common_real`.

1. **Common response for same-index legs** (`LEG_COMMON`), in two parts that match the measurement
   above:
   - LOW terms: one joint regression over both legs (normal equations summed in the `Exchange`), so
     the slow response is fixed by the long record and is the same at both samplings.
   - HIGH terms: regressed on the calibration leg only and applied unchanged on the other leg, which
     then scores them. The split point (`LEG_COMMON_SPLIT`, a harmonic number) is a setting.
   With this the search itself optimizes the manifold for a fast response that holds everywhere,
   instead of for two independent fits.
2. **Month-mean evaluation for the monthly leg.** A monthly value is a mean over the month; the model
   is evaluated at one date. For harmonics above about 30 the two differ. Needs an oversampled
   evaluation (model on a 5-day grid, averaged) in the monthly leg: `OVERSAMPLE N`.
3. **Stationarity as the composite.** With 1 and 2 in place, the monthly leg's score for the HIGH part
   is the stationarity measure, per era if the monthly leg is split into several legs by interval
   (N legs already allow this: one 5-day calibration leg, several monthly era legs).
4. **Manifold-accuracy control.** Repeat the measurement with HIGH limited to its lowest member (81)
   and to synthetic data built from the fitted model plus noise: if a known stationary fast response is
   not recovered outside the window either, the limit is the manifold, not the physics.
5. **Per-leg scale factor**, so the two series need not be in identical units.
6. **Report the leg-to-leg response agreement at save**, next to `COMPOSITE`.

## Other open items

7. **Composite options:** `min` across legs as an alternative to the weighted mean.
8. **Cross-interval scores in the status line**, for display.
9. **Held-out intervals:** a single six-year window swings widely (AMO+NAO: -0.56 with 2000-2005 held
   out, +0.46 with 1980-1985). Rotate several windows, or use era legs.
10. **dLOD bound:** keep each line within a set number of dLOD standard errors during the search.
11. **Asymmetry for the eastern Nino indices:** harmonics, `ALPHA_GAIN` and a power law do not give the
    skew; it is concentrated in 1976-1999 (see `legs/nino_skew/east_west.py`).
12. **`VALIDATE`, `ALTERNATE`, `SPLIT_TRAINING`** with legs; **`TIDES ZONAL`** cache with two time grids.
13. **Plot scripts and `lte_gui.py`:** accept `lte_results.N.csv`.
14. ~~First test run in a fresh directory printed nothing.~~ Not a program fault: my launch line put the
    scoring run in the background, so its output was never shown.
