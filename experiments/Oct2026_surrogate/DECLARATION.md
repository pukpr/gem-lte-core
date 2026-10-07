# Surrogate-data parsimony test: declaration (2026-10-06 01:20, before any surrogate exists)

## Claim under test
The cascaded LAYER=1 formulation
  sum_j a_j sin(2 pi k_j A sin(2 pi k_f (F + B sin(2 pi k_b F + phi_b)) + phi_f) + phi_j)
fits AMO, PDO and NINO4 because they are phase-locked to the tidal forcing F, not because the
formulation is flexible enough to fit any series with the same spectrum.

## Series and span
Feb2026/{amo, pdo, nino4}/<index>.dat, current fits (lt.exe.p and lt.exe.resp as of the snapshot
copied into this directory). Common span 1880.0-2022.083 (1706 months, the PDO length). The real
series are truncated to that span and pass through exactly the same procedure as the surrogates
("draw 0").

## Surrogates
Multivariate IAAFT (iterated amplitude-adjusted Fourier transform, 200 iterations): each series
keeps its power spectrum and its value distribution. The SAME random phases are used for all three
series in a draw, so their cross-correlations are kept. Phase-locking to F is what is destroyed.
19 draws (seed 20261006). With the real series that makes 20 per index, so a real series ranked
first among 20 corresponds to p = 0.05.

## Two arms

Arm A, fixed manifold, outer degrees of freedom only (independent Python engine).
Each index's layered manifold (column 4 of its current fit) is held fixed. For the real series
and each surrogate the same regression is solved: level, k0*F2, linear trend,
annual/semiannual cos/sin, and ONE winding sin/cos(2 pi k F2), with k chosen by grid
(0.25 to 10.00, step 0.01) for best training CC. F9 (3-point [1/4, 1/2, 1/4]) is applied to
the data and to every column. Training = all months outside 2000-2005; test = 2000-2005.
Caveat stated in advance: the manifold was tuned on the real series, so Arm A favours the
real series; a FAIL in Arm A is therefore strong evidence against the claim, while a PASS is
weak.

Arm B, full lt.exe search, equal budget, cross-seeded start.
Real and surrogates start from the SAME parameters, taken from a different index's fit so the
start is not tuned to the target: AMO from the PDO fit, PDO from the AMO fit, NINO4 from the AMO
fit (each index keeps its own lt.exe.resp, apart from NM/NH, which come with the seed's .p). One
300 s CC search each (4 threads, the frozen Oct2026_layer_ir/lt.exe binary, settings from the
index's lte_run.sh), then a TEST_ONLY re-score. Scores: train, validate, test.

## Statistic and rejection rule (fixed now)
Primary: training CC. Secondary: validate and test (2000-05).
For each index and arm, the claim PASSES only if the real series' primary score is higher than
all 19 surrogates (rank 1 of 20, p <= 0.05). It is reported as the rank plus the margin over the
best surrogate, in units of the surrogate standard deviation.
Overall: the claim is supported only if Arm B passes for at least 2 of the 3 indices. Arm A is
reported, but by itself it cannot support the claim.

## Controls
Positive (run first): a planted series y = the real fit's own model output (column 2 of
lte_results.csv) plus IAAFT noise matching the real residual, so it is phase-locked by
construction. Arm A must rank it first against IAAFT surrogates of itself. If it does not,
the harness lacks power and Arm A is void.
Negative: a single surrogate treated as if it were the "real" series, ranked against the other
18 using the same code. It must NOT pass (expected rank roughly uniform). Run for Arm A for all
three indices.
