# Nine-leg QBO fit: KIT 70–15 hPa + ERA5 5-day 30/50/70 hPa, 2015–2021 excluded (2026-10-26)

**Setup.** `make_leg_run.py`.
- **Legs 1–6:** KIT monthly at 70, 50, 40, 30, 20 and 15 hPa.
- **Legs 7–9:** ERA5 5-day means at one equatorial grid point (0°E, 00 UTC; `../../QBO/get_era5_point.py`) at 30, 50
  and 70 hPa, `SAMPLING 73`.
- **Hold-out:** `EXCLUDE=true TRAIN_START=2015 TRAIN_END=2021` on every leg.
- **Seed and settings:** the seed is the six-leg result (`../qbo6_kit_ex2015_2021/result.p`, kept as
  `seed_qbo6_result.p`); template and environment as in that run (`BZONE=TRUE IDATE=1920.9 TIDE_LEGACY=TRUE`).
- **Search:** 3600 s on 16 threads. Result: `result.p`. Scores: `excluded_region.py seed|result` (legs read from
  `lt.exe.resp.<k>`).

**ERA 70 hPa covers only 1979–2020.45** (the existing `../qbo70_point0E` file). It includes the 2016 disruption and the
start of the 2020 one, but not late 2020–2021. A full 1979–2026 series is downloading
(`../../QBO/get_era5_point.py 70 0 0 12 00`); rerun with it when it arrives.

## Scores (held-out = 2015–2021)

| leg | train r | held-out r (seed → result) | rms / sd | 2016 disruption | between 2017–19.5 | 2020 disruption |
|---|---|---|---|---|---|---|
| KIT70 | 0.72 | 0.61 → 0.62 | 0.94 | +0.27 | +0.90 | +0.89 |
| KIT50 | 0.78 | 0.74 → 0.74 | 0.68 | +0.10 | +0.92 | +0.76 |
| **KIT40** | 0.82 | 0.79 → **0.81** | **0.57** | −0.07 | +0.91 | +0.91 |
| **KIT30** | 0.85 | 0.78 → **0.79** | 0.75 | **+0.54** | +0.91 | **+0.95** |
| KIT20 | 0.78 | 0.77 → 0.78 | 0.76 | +0.25 | +0.96 | +0.70 |
| KIT15 | 0.74 | 0.70 → 0.70 | 0.75 | −0.09 | +0.93 | +0.83 |
| ERA30 | 0.84 | 0.69 → 0.71 | 0.77 | +0.28 | +0.82 | +0.54 |
| ERA50 | 0.75 | 0.66 → 0.67 | 0.79 | +0.03 | +0.91 | +0.54 |
| ERA70 (to 2020.45) | 0.65 | 0.55 → 0.56 | 0.95 | +0.22 | +0.82 | +0.22 |

## Findings

- **Adding the ERA legs did not change the altitude ranking for the excluded window:** KIT 40 hPa (r 0.81, lowest
  rms/sd 0.57) and KIT 30 hPa (0.79) fit best, then 20 and 50 hPa; the column edges (15, 70 hPa) fit worst. The
  shared manifold moved little; the KIT held-out scores rose by 0.00–0.02.
- **ERA fits worse than KIT at the same level,** at 30 hPa by 0.08 and at 50 hPa by 0.07 held out, with a train/held-out
  drop of 0.08–0.14. The ERA legs are a single grid point sampled once a day: they carry short-period variability and
  local noise that a monthly station mean averages out (clearest at 70 hPa). At 5-day resolution the model's own
  impulse-and-hold steps also show as abrupt jumps (e.g. ERA30, 2017.8 and 2018.2).
- **The 2016 disruption is still not reproduced.** Best is KIT30 at +0.54, up from +0.41 in the six-leg run; the
  others range from −0.09 to +0.28. ERA30 shows the event as a brief drop toward zero inside a westerly; the model
  instead goes fully easterly and recovers.
- **The 2020 disruption** is carried in KIT (0.70–0.95) but not in ERA (0.22–0.54). Over 2019.5–2020.5 the ERA30
  westerly fades slowly while the model turns easterly. KIT's monthly averaging hides part of that mismatch.
- **2017–2019.5 is fitted everywhere** (r 0.82–0.96).

## Caveats

- The six-leg seed inherits the original `qbo30_point0E` seed's exposure to 2015–21 at 30 hPa. The held-out scores
  are not strictly out-of-sample (see `../qbo6_kit_ex2015_2021/RESULTS.md`).
- The sub-window correlations rest on 15–30 monthly values (about 110 five-day values for ERA); differences of about
  ±0.2 are not meaningful.
- One one-hour search. ERA70 is truncated (above).

## Time-height section (`column_profile.py` → `column_profile.png`)

The classic QBO section, 15 hPa (top) to 70 hPa (bottom), 1979–2026, on the finest common grid (5-day, 73 per year):
- **30 / 50 / 70 hPa:** ERA5 5-day legs (ERA70 to 2020.45, then KIT70).
- **40 / 20 / 15 hPa:** KIT monthly, interpolated in time (PCHIP).
- Linear in log-pressure between levels; no smoothing; style as `docs/model_column.png`.

Panels: (a) observed, (b) model as fitted, (c) model with each level rescaled to the observed s.d. (gains 1.19–1.57;
the regression carries 64–86 % of the observed amplitude), (d, e) 2013–2023 zooms of observed vs rescaled model.

What it shows:
- **Cycle timing and descent are right (corrected 2026-10-27).** An earlier version of this note said the model's bands
  were nearly vertical. That was misread from the colour section, where the lower model amplitude, the stripes and the
  mix of monthly and 5-day rows hide the tilt. Measured lags behind 15 hPa (training years, full-ERA70 run) are:

  | | 20 hPa | 30 hPa | 40 hPa | 50 hPa | 70 hPa |
  |---|---|---|---|---|---|
  | observed (months) | 1 | 4 | 7 | 9 | 12 |
  | model (months) | 1 | 3 | 6 | 8 | 12 |

  The per-leg winding phases reproduce the descent. A sharper measure, summing lags between adjacent levels (KIT monthly,
  low-pass 200 d), gives 10.3 months from 15 to 70 hPa observed and 7.0 months in the model (about 70 % of the
  observed rate; `../../QBO/winding_lags.py`, `../../QBO/UPDOWN.md`).
- **No level-dependent manifold shift needed.** A time-shifted copy of the 30 hPa model, standing in for a manifold
  whose impulse timing moves with altitude, does not beat each level's own phase-adjusted fit on the held-out window
  (70 / 50 / 40 / 20 / 15 hPa: 0.41 / 0.65 / 0.79 / 0.76 / 0.61 vs 0.62 / 0.74 / 0.81 / 0.77 / 0.70). It is better on
  training years only at 15–20 hPa. The common manifold holds.
- **Extra high-frequency structure:** vertical stripes at the 5-day ERA levels (the model's impulse-and-hold steps).
- **2016:** the observed signature is an easterly tongue cutting into the westerly at 40–50 hPa while the westerly
  holds above. The model instead has easterlies at 15–30 hPa through 2015–mid 2016 and no tongue.
- **2017–2019:** the bands match, consistent with the leg scores.
- The remaining mismatches are the 2016 disruption, the amplitude (64–86 % of observed) and the stripe structure: the
  real sub-120-day stripes are only weakly timed by the model, and the model's stripes are too vertically coherent.
