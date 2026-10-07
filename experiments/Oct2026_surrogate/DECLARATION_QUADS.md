# Frozen-manifold low-winding test across all quads: declaration (2026-10-06, before running)

## Idea (user)
With the manifold FROZEN, the quads' data should show clearly delineated low-winding structure
against the latent level (the "LTE Modulation" / "modulation on latent" panels). No model is fitted
to the quads at all.

## Frozen manifold
F2 = column 4 of AMO's current fit (seeds/amo/lte_results.csv, LAYER=1), over the quads' common
dates 1950.0-2022.583. The quads never tuned AMO's fit.

## Data
The 89 quad series Oct2026_layer/<quad>/<quad>.dat (Kaplan boxes, 1950-2023), linearly detrended,
standardized.

## Winding spectrum (as in lt.exe ME_Power_Spectrum / the plot panel)
P(f) = |sum_t y(t) exp(i 2 pi f F2(t))|^2 for f = 0.50 ... 40.00 (step 0.05) in cycles per unit of
F2. Windings below 0.5 are excluded because they overlap the level and trend.
S_low (primary): the fraction of sum P in 0.5 <= f <= 6.0.
S_H (secondary): Hoyer sparsity of P over the full grid, (sqrt(n) - L1/L2)/(sqrt(n) - 1).
Also reported: the f of the largest peak.

## Nulls (19 each)
N1, data side: univariate IAAFT surrogate of each quad (200 iterations), same frozen F2.
N2, forcing side ("wrong clock"): the real data against F2 circularly shifted by a whole number of
years, s in {5, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47, 50, 53, 56, 59}. Whole years
keep the annual alignment and the level structure, so only the tidal timing moves.

## Rule (fixed now)
For each quad and null, the real value counts as a hit if it ranks 1 of 20 on S_low. Under the null
the hit rate is 5%. The claim is supported if, for BOTH N1 and N2, the number of hits among the 89
quads gives a one-sided binomial p < 0.001 against 5%. Also reported: the median percentile of the
real value within its null, and the same counts for S_H.

## Controls
Negative: one N1 surrogate of each quad treated as "real" against the other 18, which must give a
hit rate consistent with 5%. Positive: a planted quad y = sin(2 pi 4.0 F2) + red noise (AR1 0.9,
SNR 0.3), which must be a hit against both nulls in at least 80% of 50 plants.
