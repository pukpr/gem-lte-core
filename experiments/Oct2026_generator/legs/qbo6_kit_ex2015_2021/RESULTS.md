# Six-leg QBO fit with the 2015–2021 anomalies excluded (2026-10-26)

**Setup.** `make_leg_run.py`. Six KIT radiosonde legs (70, 50, 40, 30, 20, 15 hPa, monthly 1953–2026) on one shared
manifold and set of windings; each leg regresses its own amplitudes.
- **Exclusion:** `EXCLUDE=true TRAIN_START=2015 TRAIN_END=2021`, so the anomaly years (the 2016 and 2020
  disruptions) are held out on every leg.
- **Settings:** template and seed from `qbo30_point0E` (its `lt.exe.resp` and `result.p`), with
  `BZONE=TRUE IDATE=1920.9 TIDE_LEGACY=TRUE` as in that run. The seed is kept as `seed_qbo30_point0E.p`.
- **Search:** 3600 s on 16 threads. Result: `result.p`.
- **Data:** the 40, 20 and 15 hPa files were generated from `../../QBO/kit_qbo.dat` (fixed value columns). The same
  parser reproduces the existing 30/50/70 files exactly.
- **Analysis:** `excluded_region.py [tag]` → `excluded_region_<tag>.png/.json`, for `seed` and `result`.

## Scores

| hPa | train r (seed → result) | held-out 2015–21 r (seed → result) | rms / sd | 2016 disruption, 2015.5–17.0 | between, 2017.0–19.5 | 2020 disruption, 2019.5–20.75 |
|---|---|---|---|---|---|---|
| 70 | 0.50 → 0.72 | 0.29 → **0.61** | 0.94 | +0.27 | +0.88 | +0.89 |
| 50 | 0.62 → 0.79 | 0.49 → **0.74** | 0.69 | +0.06 | +0.93 | +0.75 |
| 40 | 0.76 → 0.82 | 0.74 → **0.79** | 0.60 | −0.23 | +0.93 | +0.88 |
| 30 | 0.86 → 0.85 | 0.83 → **0.78** | 0.76 | **+0.41** | +0.91 | **+0.94** |
| 20 | 0.71 → 0.78 | 0.68 → **0.77** | 0.78 | +0.20 | **+0.96** | +0.70 |
| 15 | 0.58 → 0.74 | 0.45 → **0.70** | 0.76 | −0.09 | +0.94 | +0.83 |

(`leg_scores2.py` gives the same held-out values within ±0.01; it treats the window endpoints slightly differently.)

## Which altitudes fit the excluded region best

- **Over the whole excluded window, 40 and 30 hPa fit best** (r 0.79 and 0.78; 40 hPa has the lowest rms/sd, 0.60),
  then 20 (0.77) and 50 (0.74). The edges are weakest: 15 hPa (0.70) and especially 70 hPa (0.61).
- **2017–2019.5, between the disruptions, is fitted at every altitude** (r 0.88–0.96). Outside the two disruptions,
  the excluded window is ordinary for the model.
- **The 2020 disruption is mostly carried** (r 0.70–0.94; best 30 and 70 hPa). Visually, though, the model misses its
  details: at 30 hPa it dips to an easterly that the data only approach, and at 50–70 hPa the sharp early-2020
  easterly dip is under-reproduced.
- **The 2016 disruption is not reproduced at any altitude** (r −0.23 to +0.41, best 30 hPa). The figure shows two
  failure modes:
  - **15–30 hPa:** the model brings in the scheduled easterly in 2016, while the observed westerly stalls and
    persists into 2017.
  - **40–50 hPa:** the observed short, sharp easterly dip in early 2016 appears in the model only weakly and late.

  In both cases the observation departs from the schedule. This matches the control reading (`docs/qbo_control.md`):
  the 2016 disruption is one out-of-schedule input, and 2020 is a milder departure that the scheduled forcing almost
  absorbs.

## Caveats

- **The seed is not blind to 2015–21:** `qbo30_point0E` was fitted on 1953–2026 minus 2000–05, so it saw 2015–21 at
  30 hPa. The search excluded 2015–21 and moved every leg (training r rose at five levels). Still, the held-out
  values, especially at 30 hPa, are not strictly out-of-sample. A seed fitted with 2015–21 excluded from the start would
  remove this.
- The sub-window correlations rest on 15–30 monthly values each, so differences of ±0.2 between levels are not
  meaningful.
- One search (one hour). `SELECT_TOP 3` picked candidate 1 (held-out 0.705 mean) over candidate 7 (0.732) on
  training score; held-out scores were not used for selection.
