# 100 hPa wind at 0E (5-day) with the unified SOI seed, 150 hPa at 180E as second leg (2026-10-27)

> **CRITICAL AREA FOR FURTHER ANALYSIS (on hold, 2026-10-27).** This is the finest-resolution ENSO fit available (5-day ERA5
> wind) and the held-out window shows coherence at both scales. Paused by request; resume with the "Next" items below.
> Note for that work: at the 5-day scale the held-out skill with the seasonal cycle removed is +0.25 (raw +0.41, mostly
> seasonal). The +0.73 to +0.79 is on 1-year means.

## Setup

- **Leg 1, U100_0E:** ERA5 u at 100 hPa, 0N 0E, 00 UTC, 5-day means (SAMPLING 73), 1979.0–2026.8
  (`../../QBO/era_qbo_u100_point0E_5day.dat`). Raw wind, seasonal cycle included.
- **Leg 2, U150_180:** ERA5 monthly anomaly at 150 hPa, 180E (`../wind_soi_src/u150_180.dat`), as in `../wind150_soi5d`.
- **Seed:** `../wind_soi_src/seed.p`, the unified SOI fit. Settings are identical to `../wind150_soi5d` (IMPULSE_EXACT,
  INIT_DATE 1880, MAXH 20, CC). 2011–2016 is held out on both legs. The search ran 30 min on 16 threads.
- `deseason.py` (scratch, logic below): for leg 1 the raw score is mostly the seasonal cycle. So each leg is also
  scored after subtracting a training-years climatology per sample slot from both model and data, and on 1-yr
  (13-point) running means of that.

## Scores (train / held out 2011–16)

| run | U100 raw | U100 seasonal-only | U100 deseasoned | U100 deseasoned 1-yr | U150 |
|---|---|---|---|---|---|
| SOI seed, nothing searched | +0.47 / +0.37 | +0.43 / +0.36 | +0.28 / +0.11 | +0.71 / +0.30 | +0.66 / +0.89 |
| QBO 9-leg seed, nothing searched (`../wind100_150_qboseed`) | +0.41 / +0.27 | same | +0.19 / −0.06 | +0.42 / +0.40 | +0.43 / −0.45 |
| searched (`lt.exe.p`) | +0.49 / +0.42 | same | +0.32 / +0.26 | +0.80 / +0.79 | +0.73 / +0.89 |
| searched, slow winding reset to seed (`../wind100_150_slowback`) | +0.49 / +0.41 | same | +0.32 / +0.25 | +0.81 / +0.73 | +0.75 / +0.87 |

## Findings

- **The SOI seed describes 100 hPa better than the QBO seed,** with nothing searched. It also keeps 150 hPa at +0.89
  held out, where the QBO manifold gives −0.45.
- **After the search, the interannual 100 hPa wind is predicted out of sample:** +0.73 to +0.79 on 1-yr means of the
  deseasoned series. That is up from +0.30 for the seed. The held-out window has only about 6 independent yearly
  values, so treat this as suggestive. At the 5-day scale the deseasoned skill is modest (+0.25 held out). The plot
  (`../wind100_150_slowback/legs.png`) shows the model following the slow swings, not the 5-day spikes.
- **The search moved the slow winding off the unified model:** −0.1376 (0.68× base) → −0.2056 (1.02× base). That
  is almost degenerate with the base winding 0.2016, and both legs then carry large, nearly cancelling amplitudes on
  the pair (100 hPa 3.45 / 3.26, 150 hPa 10.5 / 9.0). **It isn't needed:** with the slow winding put back, every score
  is within 0.06 and the 150 hPa training score is higher.
- **Where the gain comes from.** The ann1/ann2/sem1 and IR values in `lt.exe.p` are leg 1's own regression
  coefficients (the 100 hPa seasonal cycle), not forcing terms, so their large change is expected. The shared model
  moved only slightly:
  - impulse shape: impA −3.40 → −3.37, impB −4.05 → −4.03, impC 0.012 → 0.020, offs 0.018 → 0.023;
  - tidal amplitudes by at most 0.0007;
  - dLOD stays at 0.9993.
  The held-out gain therefore comes from small changes to the impulse/manifold within the dLOD constraint, not from
  new windings: the harmonic set is unchanged from the seed.
- **What the 100 hPa signal is.** Deseasoned monthly data, 13-month means:
  - against 150 hPa at 180E, r −0.44 at lag 0;
  - against the QBO at 50 hPa, ≤ +0.37 and only at a 5-month lag.
  So at 0E the interannual 100 hPa wind is mainly the ENSO/Walker pattern with reversed sign (the n = 1 standing
  wave), plus a weaker QBO part. That fits the stellarator-layer picture: below 70 hPa the zonal mean carries only
  0.23 of the variance.
- **The per-winding response is not clearly shared.**
  - With the slow winding on-model, the two legs' amplitudes correlate +0.69 across 12 windings, but amplitude and
    phase together only −0.21.
  - The searched fit's +0.99 / −0.76 is driven by the degenerate pair, so it is not evidence.
  - Contrast `../wind150_soi5d`, where wind and SOI amplitudes agreed +0.77.

## Next

- Rerun the search with the slow winding locked (or the unified-model gate used for SOI) so it can't drift to the
  degenerate pair.
- Add the 5-day SOI as a third leg, to anchor the structure as in `../wind150_soi5d`.
