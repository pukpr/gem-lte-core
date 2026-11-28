# Baltic 3pt CC >= 0.70 Meeus-Free `TIDES=POLY` (`Z52`) Champion Preset

Reproduces the 4-leg Atlantic (`NAO`, `baltic_3pt_na`, `AMO`, `TNA`) Meeus-free
`TIDES=POLY` (`POLY_TIER=Z52`) champion model (`F9=1`, `TREND=TRUE`, `ACCEL=TRUE`,
`REFINE_LAM=0.30`, `POLY_N5_RIDGE=3.0`, `IR_SIDE=0`).

## Quick Verification (`TEST_ONLY`)

```bash
cd experiments/Dec2026/baltic_070_poly
ulimit -s unlimited
F9=1 TREND=TRUE ACCEL=TRUE EXCLUDE=TRUE TRAIN_START=1960.0 TRAIN_END=1970.0 \
  POLY_N5_RIDGE=3.0 ../lte_run.sh --test --tier Z52
```

## Expected Metrics (1960–1970 Held-Out Fold)
- **Baltic 3pt**: Train CC = `+0.703`, Held-Out (1960–1970) CC = `+0.287`, `MaxA = 63.1`
- **ATL4 Held-Out**: `4/4` positive (`NAO +0.057`, `Baltic +0.287`, `AMO +0.079`, `TNA +0.013`, Mean HO = `+0.109`)
- **dLOD Fidelity**: `0.9620` (vs `0.9544` for 52-line `TABLE` R5)
- **5-Fold Rolling 10-Yr OOS CV (1960–2010)**: `5/5` positive Baltic folds (Mean Train `+0.704`, Mean OOS `+0.203`, ATL4 Mean OOS `+0.108`)
