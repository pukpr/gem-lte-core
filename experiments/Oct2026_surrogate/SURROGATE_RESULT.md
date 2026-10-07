# Surrogate-data parsimony test: result (2026-10-06)

Declared in DECLARATION.md before any surrogate existed. 19 multivariate IAAFT surrogates per
index, plus the real series (draw 0), common span 1880.0-2022.083.

## Arm A: fixed manifold, one free winding (independent Python engine)
| index | positive control | real train CC: rank, margin over best surrogate | negative control |
|---|---|---|---|
| AMO | rank 1/20 (+2.1 sd) | rank 1/20, 0.761 vs 0.458 +- 0.098 (+1.7 sd) | rank 4/19 (fails correctly) |
| PDO | rank 1/20 (+3.0 sd) | rank 1/20, 0.573 vs 0.248 +- 0.044 (+5.2 sd) | rank 8/19 |
| NINO4 | rank 1/20 (+7.2 sd) | rank 1/20, 0.599 vs 0.240 +- 0.044 (+5.5 sd) | rank 10/19 |
As declared, Arm A favours the real series (the manifold was tuned on them), so this PASS is
weak.
Post-hoc, not certifying: the real series all chose k = 4.06 / 4.02 / 3.99 on the free grid,
while surrogates scatter over 0.25-10 (only 1, 3 and 2 of 19 fall within 3.94-4.12).

## Arm B: full 300 s lt.exe search, equal budget, cross-seeded start (decisive arm)
| index | train (primary): real, rank, surrogates | validate rank | test (2000-05) rank |
|---|---|---|---|
| AMO | 0.696, rank 1/20, 0.656 +- 0.017 (max 0.686; +0.63 sd) | 6/20 | 6/20 |
| PDO | 0.611, rank 13/20, 0.615 +- 0.011 | 9/20 | 3/20 |
| NINO4 | 0.413, rank 3/20, 0.375 +- 0.048 (max 0.519) | 2/20 | 1/20 (0.539 vs 0.009 +- 0.229) |
Rule: the claim is supported only if at least 2 of 3 indices rank first on train. Result: 1 of 3,
with AMO by a thin margin. Arm B FAILS.

## Verdict: claim NOT supported at full search freedom
With the full search (nominally ~150 search entries, but the tidal ones stay locked by dLOD to ~1% in amplitude and ~0.5 deg in phase, so ~7 nonlinear + ~13 linear effective parameters), the formulation fits spectrum-matched,
phase-scrambled series nearly as well as the real indices. In-sample CC at that number of degrees
of freedom is therefore not evidence of tidal phase-locking. Arm B surrogate searches also ended
with regressed k of about 4 (3.93-4.30), as the real series did: the search barely moves k from the
seed, so k = 4 in Arm B reflects the seed, not the data.

Bias noted: Arm B starts were tuned on a different REAL index, which could only favour the real
series. The failure is therefore robust to that bias.

Residual positive signals (secondary or post-hoc only): the NINO4 test window ranks first in Arm B;
Arm A real series rank first with a common k of about 4.

Implication for parsimony: the test needs the restructured estimator, in which only the six shared
inner parameters are searched (jointly over indices) and each index's outer sum is solved by
regression, applied identically to surrogate triplets. Only then are the degrees of freedom small
enough for a fit to be able to fail.

# Level-structure test (DECLARATION_LEVELS.md), low windings only, 2026-10-06
Raw output: levels_result.txt.

S1 (primary): adjusted correlation ratio of the detrended data by 12 F2-level bins, no windings.
| index | cross manifold: real, rank, margin over best | own manifold: real, rank, margin | positive control (cross) | negative control |
|---|---|---|---|---|
| AMO (PDO manifold) | 0.366, 1/20, +0.81 sd (surr 0.158+-0.088) | 0.456, 1/20, +2.08 sd | rank 1 | rank 5/19 |
| PDO (AMO manifold) | 0.113, 3/20 (surr 0.067+-0.032) | 0.263, 1/20, +5.17 sd | rank 2: FAILED | rank 9/19 |
| NINO4 (AMO manifold) | 0.105, 1/20, +0.23 sd (surr 0.058+-0.025) | 0.212, 1/20, +3.71 sd | rank 1 | rank 9/19 |
Declared rule (cross manifold, at least 2 of 3 ranked first): PASSES (AMO and NINO4), with thin
margins. PDO's positive control failed on the AMO manifold, so the cross test has no power for PDO:
its rank 3 is not evidence against.

S2 (penalized single winding, BIC with 3 parameters): own manifold dBIC real -483 / -526 / -560
against surrogates -73 / -33 / -37 (rank 1 for all three), chosen k = 4.06 / 4.01 / 3.99 (real);
surrogates choose k scattered over 0.25-10. Cross manifold: PDO and NINO4 rank 1 (k 3.62, 3.68),
AMO rank 6.

Reading: with high windings excluded, the data are organized by the latent level beyond what any
spectrum-matched surrogate achieves, clearly on each index's own manifold and marginally on a
manifold the index never tuned (AMO, NINO4).

# YEAR detuning, re-search arm Y2 (DECLARATION_YEAR.md), 2026-10-06
P (level eta2_adj) and training CC after 300 s re-searches; dLOD in brackets:
- AMO: delta=0 0.484 / 0.797 (0.988); +-0.01: 0.038, 0.127 (dLOD -0.77); +-0.05: 0.118, 0.154 (-0.16).
- PDO: delta=0 0.273 / 0.754 (0.989); +-0.01: 0.078, 0.046; +-0.05: 0.077, 0.082.
- NINO4: delta=0 0.212 / 0.734 (0.988); +-0.01: 0.054, 0.036; +-0.05: 0.140, 0.129.
Y2 rule passes for all three, so Y1 and Y2 both pass. CAVEAT (raised by the simulated referee): the
delta=0 run starts at its own long-searched optimum, and no detuned run re-pinned its tides to
dLOD (dLOD stays negative), so Y2 favours the incumbent. A cold-start "fake Moon" null (lunar
frequencies scaled by 1-3%) is needed for lunar specificity.
