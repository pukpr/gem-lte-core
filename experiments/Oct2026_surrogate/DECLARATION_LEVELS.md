# Low-winding / level-structure surrogate test: declaration (2026-10-06, before running)

## Rationale (user direction)
High windings act like a nonlinear feature layer and can fit almost any series, so they are
either excluded or explicitly penalized. The test targets what the LTE Modulation chart shows:
the data organized by the level of the latent layered manifold F2.

## Data, manifolds, surrogates
The same 19 multivariate IAAFT surrogates and real series (draw 0) as DECLARATION.md (common span
1880.0-2022.083). Each series is linearly detrended over the span. F2 is column 4 of a FIXED fit.
- own manifold: the index's own current fit (biased toward the real series, which tuned it);
- cross manifold, the primary case: another index's fit (AMO data on the PDO manifold, PDO data on
  the AMO manifold, NINO4 data on the AMO manifold), which the target's data never tuned.

## Statistics
S1 (primary): adjusted correlation ratio of the data by F2 level,
  eta2 = between-bin variance / total variance, with F2 cut into B = 12 equal-count bins,
  adjusted: eta2_adj = 1 - (1 - eta2)(n - 1)/(n - B).
  No windings, no regression on F2 beyond the bin means, so the degrees of freedom are the 12 bin
  means.
S2 (secondary, penalized windings): BIC improvement from adding ONE winding sin/cos(2 pi k F2),
  k on 0.25..10 (step 0.01), to a baseline of [1, F2, trend], counting 3 extra parameters
  (amplitude, phase, k) with the BIC log(n) penalty. Report dBIC (negative = winding justified)
  and the k selected.

## Rule (fixed now)
S1 with the cross manifold: the real series must rank 1 of 20 for at least 2 of the 3 indices.
S1 with the own manifold and S2 are reported, but by themselves they do not decide the verdict.

## Controls
Positive: the planted series (real fit output plus IAAFT noise), scored on the cross manifold,
must rank 1 of 20 against its own IAAFT surrogates under S1. Negative: s01 ranked against
s02..s19 (S1, cross manifold) must not rank 1.
