# Nine-leg QBO fit with the full ERA5 70 hPa series (2026-10-26)

Same as `../qbo9_kit_era_ex2015_2021` (KIT 70–15 hPa monthly + ERA5 5-day point 30/50/70 hPa, 2015–2021 excluded),
with these differences:
- **ERA 70 hPa:** the complete 1979–2026.79 series (`../../QBO/era_qbo_u70_point0E_5day.dat`, 00 UTC, made like
  30/50) replaces the truncated file, which ended in 2020.45 and came from a different extraction (r 0.96 with the
  new series).
- **Seed:** `../qbo9_kit_era_ex2015_2021/result.p`, kept as `seed_qbo9_result.p`.
- **Search:** 1 h, 16 threads. Outputs: `result.p`, `excluded_region_result.*`, `column_profile.png`.

The search moved every leg by ≤ 0.01. Held-out 2015–21 r: KIT40 0.81, KIT30 0.79, KIT20 0.77, KIT50 0.74, KIT15 0.70,
KIT70 0.62; ERA30 0.71, ERA50 0.67, ERA70 0.52. ERA70 sub-windows: 2016 disruption +0.16, between +0.75, 2020
disruption +0.58. The altitude ranking and the other findings are unchanged; see the earlier `RESULTS.md`.

**Correction (2026-10-27):** the descent is reproduced (model lags behind 15 hPa of 1 / 3 / 6 / 8 / 12 months against
observed 1 / 4 / 7 / 9 / 12 at 20 / 30 / 40 / 50 / 70 hPa). A delayed common model does not beat the per-leg
phase-adjusted fits held out. See `../qbo9_kit_era_ex2015_2021/RESULTS.md`. Summing lags between adjacent levels
instead (low-pass 200 d) gives 7.0 months from 15 to 70 hPa in the model against 10.3 observed: the model descends at
about 70 % of the observed rate (`../../QBO/winding_lags.py`).
