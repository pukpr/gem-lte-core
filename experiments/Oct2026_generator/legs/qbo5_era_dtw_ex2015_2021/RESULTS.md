# Five-level ERA QBO fit, METRIC=DTW, 2015–2021 held out (2026-10-27)

## Setup

- **Legs:** ERA5 u at 0N 0E, 5-day samples (SAMPLING 73), at 10, 20, 30, 50 and 70 hPa. All five share one manifold.
- **Seed:** `qbo9_kit_era70full_ex2015_2021/result.p`.
- **Environment:** BZONE=TRUE, IDATE=1920.9, TIDE_LEGACY=TRUE.
- **Search:** 3600 s on 16 threads, 14:53–15:54. The settings match `../qbo5_era_ex2015_2021` (the same search with METRIC=CC) except for the metric.

## Scores

Pearson r, train / held-out 2015–2021, from `run.out`:

| leg | seed | CC search | DTW search |
|---|---|---|---|
| ERA10 | 0.667 / 0.621 | 0.771 / 0.110 | 0.749 / **0.571** |
| ERA20 | 0.768 / 0.679 | 0.840 / 0.322 | 0.823 / **0.527** |
| ERA30 | 0.844 / 0.706 | 0.874 / 0.274 | 0.875 / **0.399** |
| ERA50 | 0.753 / 0.669 | 0.807 / 0.283 | 0.808 / **0.531** |
| ERA70 | 0.640 / 0.519 | 0.699 / 0.132 | 0.669 / **0.460** |

The DTW search improves training about as much as CC does, and keeps most of the held-out skill. It is still below the seed on held-out at every level.

## Held-out pattern check (`heldout_pattern.py` → `heldout_pattern.png/.json`)

DTW is the mean |Δ| of the z-scored series with a ±90-day warping band. The null shifts the model randomly in time (30 draws).

| hPa | DTW: seed / CC / DTW run | null (DTW run) | sd held, model/obs: seed / CC / DTW run | obs peaks matched within 60 d (DTW run) |
|---|---|---|---|---|
| 10 | 0.36 / 0.79 / 0.44 | 0.89 ± 0.13 | 0.90 / 0.72 / 0.98 | 0.67 |
| 20 | 0.42 / 0.66 / 0.49 | 0.86 ± 0.11 | 1.01 / 0.73 / 0.82 | 1.00 |
| 30 | 0.36 / 0.61 / 0.54 | 0.86 ± 0.10 | 0.83 / 0.68 / 0.68 | 0.75 |
| 50 | 0.30 / 0.61 / 0.43 | 0.85 ± 0.20 | 0.92 / 0.75 / 1.09 | 0.75 |
| 70 | 0.41 / 0.62 / 0.41 | 0.75 ± 0.13 | 0.66 / 0.57 / 0.70 | 0.47 |

- **No overfitting signature.** The held-out amplitude ratios stay at or below the training ratios, apart from 10 and 50 hPa, which are about 1. The largest held-out peaks are 0.6–1.3 times the observed ones.
- **Pattern beats the null at every level.** It does so by 2.3–3.4 sd, against 0–2 sd for the CC run.
- **30 hPa is weakest.** Its best lag is −100 d (r 0.39), so the timing there is off by about 3 months.

## Descent (`column_profile.py` → `column_profile.png/.json`)

QBO-phase lags between adjacent levels, low-pass 200 d, cumulative 10→70 hPa:

| | 50→70 | 30→50 | 20→30 | 10→20 | 10→70 |
|---|---|---|---|---|---|
| observed | 66 d | 111 d | 76 d | 119 d | 372 d |
| seed | 47 d | 78 d | 45 d | 67 d | 237 d (64%) |
| CC search | 49 d | 91 d | 54 d | 81 d | 275 d (74%) |
| DTW search | 48 d | 84 d | 56 d | 74 d | 261 d (70%) |

Both searches make the descent slower, moving it towards the observed rate. Most of the shortfall is at 10→20 and 30→50 hPa.

Stripes (high-pass 120 d) are still essentially vertical in the model: 0.4–0.8 d per level pair (r 0.77–0.83). The exception is 30→50 at 8.3 d, with r of only 0.17. The observed values are 2.7–4.9 d, at r 0.15–0.33.

## What the search changed

- **IR** went from 0.789 to 0.891 (the CC search went to 0.944).
- **Second winding** went from k −0.553 to −1.310 (CC: −2.945). The other eight windings moved by at most 0.04.

The fixed-manifold IR profile (`../../QBO/ir_criticality.py`) refits all linear coefficients for IR from 0 to 0.99. Held-out r at IR 0.79:

| manifold | 10 | 20 | 30 | 50 | 70 hPa |
|---|---|---|---|---|---|
| 9-leg (seed) | 0.60 | 0.67 | 0.70 | 0.67 | 0.52 |
| DTW search | 0.54 | 0.52 | 0.39 | 0.53 | 0.44 |
| CC search | 0.10 | 0.31 | 0.26 | 0.29 | 0.15 |

Each curve is flat for IR ≳ 0.4. So the ranking of the three runs comes from the manifold and the second winding, not from IR.

## Reading

- The metric matters. CC rewarded transitions in the training years that do not carry over to the held-out years. DTW, being phase-tolerant, kept the held-out pattern.
- Neither search beat the seed's manifold on held-out data. The ERA-only legs seem to pull the second winding away from k ≈ −0.55. That winding is what the seed relies on for 2015–21.
- Possible next run: DTW with the second winding locked at the seed's value, or a DTW run seeded from this result with IR locked at 0.79. Either would check whether the held-out loss is the k −0.55 → −1.31 move alone.
