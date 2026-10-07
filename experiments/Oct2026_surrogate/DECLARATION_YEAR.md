# Forcing-side "wrong clock" test via YEAR detuning: declaration (2026-10-06, before running)

## Rationale
YEAR (resp) is added to the year length (365.2422484 days + YEAR), and every tidal period is
converted to calendar years through it. Offsetting YEAR by delta detunes the whole tidal forcing
against the annual impulse comb and the calendar, while keeping the data, the model structure and
its degrees of freedom. An error of delta days drifts Mf's phase by about
2 pi x 26.7 x delta / 365 rad per year: about 0.6 rad over 140 yr at 0.01 day, a full scramble near
0.05 day.

## Setup
AMO, PDO, NINO4: the current fits (snapshot in seeds/), each index's own resp and lte_run.sh, and
the frozen binary ../Oct2026_layer_ir/lt.exe. Baseline YEAR = 0.00405305 (all three). Offsets
delta in {0, +-0.0005, +-0.001, +-0.002, +-0.005, +-0.01, +-0.02, +-0.05} days.

Arm Y1, fixed parameters (TEST_ONLY; regression re-solved): the wrong clock with the same model.
Arm Y2, equal-budget re-search: 300 s search at delta in {0, +-0.01, +-0.05} from the same start, so
the model can re-adapt to the wrong clock. Runs with dLOD < 0.99 are reported, but flagged.

## Statistics (low-winding, as in DECLARATION_LEVELS.md, on each run's own F2 = column 4)
P (primary): adjusted correlation ratio of the detrended data by 12 F2-level bins.
Secondary: training CC, validate and test from lt.exe, dLOD, and the BIC-penalized single-winding
dBIC.

## Rule (fixed now)
Y1: delta = 0 must be the best value of P among all offsets, for at least 2 of the 3 indices.
Y2: delta = 0 must beat both |delta| = 0.05 runs on P, for at least 2 of the 3 indices.
The tidal-timing claim is supported only if both Y1 and Y2 pass.

## Control
dLOD must be best at, or within 0.001 of, delta = 0 in Y1. That shows the detuning is felt by the
LOD calibration, which is the physical anchor.
