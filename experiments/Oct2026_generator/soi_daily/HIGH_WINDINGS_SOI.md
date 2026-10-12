# High-order windings in the SOI: an appraisal

Question: are the rapid excursions in the daily SOI part of the unified LTE/tidal response, carried by
windings too fast to see at monthly sampling, or are they weather?

Verdict on the evidence to date (2026-10-11): **not established.** High windings leave a real,
measurable imprint in the years they are fitted to, and the mechanism by which they would hide at
monthly sampling is confirmed. But every out-of-sample and control test so far fails to distinguish
them from arbitrary fast terms. What is established is the slow and interannual part of the unified
response, at both poles and both time scales.

All numbers are from `soi_daily/` and `soi_universal/` (scripts named at the end).

## Definitions

- **Winding rate.** A winding k cycles k × (plateau slope) times per month between yearly impulses.
- **High winding.** More than half a cycle per month, the monthly Nyquist limit. At monthly sampling
  such a term appears only as an irregular month-to-month pattern, never as a visible oscillation.
- **Unified model.** Tidal table within the dLOD gate, the shared manifold from 1880, one low
  winding plus the ~0.2 backbone and its integer harmonics.

## How much fast variance there is to explain

| Period | Share of daily SOI variance |
|---|---|
| Under 10 days | 18% |
| 10 to 30 days | 24% |
| 30 to 90 days | 20% |
| Over 90 days | 38% |

## Evidence that high windings contribute

1. **In-sample imprint.** In the best 5-day fit of 1992–2005 (CC 0.75), eight windings with periods
   of 0.75 to 1.8 months raise the fit from 0.705 to 0.750 and correlate +0.36 with what the low
   windings leave.
2. **The alias is real and survives monthly averaging.** Averaged to calendar months, that
   high-winding part keeps 44% of its variance, moves to periods of 2 to 12 months, and correlates
   +0.38 with the monthly SOI's residual (+0.29 for the 1992–2015 fit). So fast windings can act at
   monthly resolution without being visible as cycles.
3. **Sloped plateaus generate the right amount of sub-annual variance.** With friction ×5 to ×20 or
   a leak of 0.02, the model's variance at 60 days to 1 year is 24–33% against 21% in the data, and
   the sub-annual correlation after the training window is +0.10 to +0.18, against −0.09 with flat
   plateaus.
4. **A fortnightly line.** In the daily SOI, Mf (13.66 d) stands at 2.2 times the local background
   (96th percentile), coming from Tahiti (98th) and not Darwin. One line of ten tested.

## Evidence that they are not yet distinguishable from noise

1. **The margin over surrogates comes from the low windings.** The 0.75 fit is 3.5 sd above
   phase-randomized surrogates (which reach 0.57). That margin is in place after the four lowest
   windings (0.60 against 0.46). The remaining fourteen raise the real and surrogate scores equally
   (to 0.76 and 0.57).
2. **Fitted high windings do no better than random ones.** On the monthly record outside the fitted
   window, the eight fitted high windings improve the regression by +0.6 sd relative to 200 random
   sets from the same range; evaluated as monthly means, −2.3 sd. The 1992–2015 fit gives −1.0 and
   −0.8 sd.
3. **Nothing faster than six months survives out of sample.**

   | Test | Faster than 6 months, trained | Held out |
   |---|---|---|
   | Free 5-day fit, after its window | +0.47 | −0.02 to +0.04 |
   | Unified model, harmonics to 20, six legs | +0.06 to +0.14 | −0.30 to −0.02 |
   | Unified model, harmonics to 250, six legs | +0.08 to +0.13 | −0.25 to +0.02 |

4. **Below 20 days the model has almost no amplitude.** Spread 0.07 to 0.19 against 0.84 to 1.10 in
   the data; training correlation at most +0.12.
5. **Time-warping does not help.** With ±15 days of DTW slack the fast part aligns at +0.65, and
   surrogates align at +0.71.
6. **The model has no lag between the poles; the data do.** See the next section. At zero lag the
   two stations' fast variance is uncorrelated, and the fitted models are in step (+0.77), so the
   model's fast component does not have the structure the data show.

## The poles are linked in the fast band, with a lag (corrected 2026-10-11)

An earlier version of this note concluded from the zero-lag correlation (+0.04 daily) that each
station's fast variance is local. That inference was wrong: a lag scan shows a significant travelling
relationship. Daily station anomalies, 1991–2026, Darwin against sign-reversed Tahiti
(`pole_lag.py`, `pole_lag.png`):

| Time scale | Lag 0 | Largest | At lag | Most negative | At lag | Against surrogates |
|---|---|---|---|---|---|---|
| 20 to 60 days | −0.14 | +0.32 | Tahiti 15 d earlier | −0.25 | Tahiti 10 d later | +9.9 sd |
| 15-day smoothed, sub-annual | −0.04 | +0.27 | 15 d earlier | −0.10 | 5 d later | +11.3 sd |
| Faster than 20 days | −0.10 | +0.10 | 10 d earlier | −0.10 | 0 | +7.8 sd |
| Slower than 1 year | +0.71 | +0.71 | 10 d earlier | | | +3.3 sd |

- In the 20-to-60-day band the correlation swings from +0.32 to −0.25 within 25 days of lag, so the
  shared variance has a period near 50 days. At that period "reversed Tahiti 15 days earlier" and
  "same-sign Tahiti 10 days later" are the same observation; the phase at other periods decides
  between them (next subsection).
- The relationship is far above chance (z near 10, allowing for the search over lags).
- A delay of two to three months does not show it in these data: at +60 to +90 days the fast-band
  correlations are between −0.06 and +0.02. On the long monthly records the full-series correlation
  is +0.32 at zero lag and +0.35 with Tahiti one month earlier, falling to +0.29 and +0.25 at two and
  three months later.
- **The unified fits do not reproduce the lag.** In the 20-to-60-day band the Darwin and Tahiti
  models correlate +0.77 at zero lag, where the data have −0.13 and peak 15 days off. The model drives
  both poles from one manifold with no phase between them.

So the fast variance is not independent weather at each pole; it contains a wave travelling between
them. What remains open is whether that wave is part of the tidal response.

### One see-saw with a 13-day offset fits every period (2026-10-12)

Phase of raw Tahiti pressure behind raw Darwin pressure, by period (`pole_phase.py`,
`pole_phase_fit.py`; cross-spectrum of the daily anomalies):

| Period (days) | 20 | 28 | 39 | 53 | 73 | 105 | 150 | 211 | 383 | 512 | 683 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Coherence | 0.23 | 0.21 | 0.27 | 0.38 | 0.47 | 0.17 | 0.26 | 0.25 | 0.49 | 0.67 | 0.53 |
| Measured phase (deg) | +3 | +32 | +12 | +59 | +95 | +112 | +113 | +141 | +157 | +146 | −177 |
| See-saw, 13-day offset | −57 | +11 | +60 | +91 | +116 | +135 | +149 | +158 | +168 | +171 | +173 |

| Description | Weighted rms phase error, all periods | Under 100 days only |
|---|---|---|
| Sign-reversed see-saw, reversed Tahiti 13 days ahead of Darwin | 28° | 35° (12.75 d) |
| Fixed spatial phase (a wave pattern between two longitudes) | 55° | 41° |
| Same-sign travelling delay | 121° | 61° |

- The phase runs smoothly from 180° at long periods toward zero at short ones. One sign reversal with
  one time offset of 13 days describes it from 20 days to 2 years. This agrees with the author's
  estimate of about 12 days.
- A same-sign delay, as an eastward wave passing first Darwin and then Tahiti would give, is rejected.
- **Relation to the MJO** (`mjo_lag.py`). Darwin's sub-annual pressure correlates +0.34 with the NOAA
  MJO index at 140E at zero lag. But the index's own phase between 140E and 120W is +63° (zero-lag
  correlation +0.33): it is not reversed over the Darwin–Tahiti span, while the stations are. An
  earlier version of this note called the fast pole-to-pole signal "the MJO's pressure signal"; the
  station phase does not support that. Within 45 to 75 days alone the two descriptions are close
  (measured 59° to 95°; see-saw 91° to 116°; MJO pattern about 55°); across all periods the see-saw
  wins.
- **Consequence for the model.** The unified fits already have the two poles sign-reversed on one
  manifold. What they lack is the 13-day offset, which is invisible at monthly resolution (hence the
  +0.88 agreement there) and decisive at 5 days.

## What is established instead

| Unified SOI fit (nothing free but harmonics) | Monthly, 1880–2026 | 5-day, 1991–2026 |
|---|---|---|
| Train | +0.57 | +0.60 |
| Held out 2000–2005 | +0.61 | +0.43 |
| dLOD agreement | 0.9994 | |
| Manifold correlation with the universal one | 1.0000 | |

Darwin and Tahiti, fitted separately on the same manifold, have the same response with opposite
sign: their winding amplitudes and phases agree at +0.88, although their monthly data correlate only
+0.31. The ninth harmonic of the backbone leads in all six fits.

## Assessment

- The claim that rapid SOI excursions belong to the unified response is **plausible in mechanism**
  (items 1 to 3 above) and **unsupported in validation** (items 1 to 6 below them).
- The honest statement today: the unified model explains the SOI down to about six months. Below
  that the two poles remain one sign-reversed see-saw with a 13-day offset, so the fast variance is
  organised in the way the model assumes; the model's fast terms have not yet been shown to predict it.
- The fast part cannot be reached by higher integer harmonics on the unified manifold, whose
  plateaus are nearly flat. It needs a sloped plateau, which is a departure from strict adoption.

## Tests that would settle it

1. **Sloped plateau as one declared departure.** Hold the table, manifold drift and backbone; free
   only the second pole `mb`. Require that the slow fit and dLOD agreement are unchanged and that the
   faster-than-six-month band becomes positive on held-out years.
2. **Random-winding control, out of sample.** Any high winding claimed must beat random windings of
   the same number on years it was not fitted to.
3. **The offset between the poles.** Fit Darwin and sign-reversed Tahiti with one response and a
   single time offset (13 days expected, reversed Tahiti first). The fast band should then score at
   both poles on held-out years with shared amplitudes.
4. **Monthly values formed as monthly means.** Evaluate the model on a fine grid and average, as the
   monthly SOI is made; point sampling scrambles fast windings (+0.15 against +0.38 in the fitted
   window).
5. **Longer daily records.** The Mf hint rests on 35 years. Daily station pressure before 1991 would
   test it directly.

## Sources

`prep_soi.py` (variance by period, tidal lines), `plateau2.py`, `plateau3.py`, `plateau4.py` (sloped
plateaus), `fastfix.py` (DTW control), `rr/attribution.py` (surrogates, term-by-term), `rr/highwind.py`
(alias and random-winding control), `rr/reconcile.py` (time shifts, filtering, initialization),
`../soi_universal/soi_universal.py`, `cross_scale.py`, `rr/rr_uni.py` (unified model, six legs).
