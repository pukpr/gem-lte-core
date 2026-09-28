# Refactor Recovery Plan

## Purpose

The current checkout, `~/eval/gem-lte-core`, is based on local commits
`807541a` and `87a8efb`.  It diverged from the prior remote/refactor history
at `2c43664`.

The displaced history is preserved in this checkout as:

```text
recovery/origin-master-before-force-push
```

The corresponding working copy is:

```text
~/refactor/gem-lte-core
```

`~/refactor/GEM-LTE` is a separate repository and is not the source for this
plan.  Treat `~/refactor/gem-lte-core` and the recovery branch as read-only
references until individual changes have been selected and validated.

## Current Capabilities to Preserve

Do not replace the current implementation wholesale.  It contains integrated
model and fitting behavior not present in the refactor copy:

- Dynamic dLOD/year calibration, including year adjustment before dLOD phase
  calibration and `CAL_LOD`.
- Forcing controls: `IMPULSE`, `FRICTION`, and bounded `JERK`.
- Ridge-regularized regression via `RIDGE`.
- Uncompensated regression/spectral scoring via `UNCOMPENSATED`.
- Alternate validation geometries: `COVERAGE`, `ENCLOSING`, and `ENCLOSED`.
- Modulation and search safeguards: `MIN_LT` validation and
  `Force_Harmonic` duplicate-period recovery.
- The expanded February 2026 research scripts for reconstruction, cross
  validation, PySR, QBO, multi-index, wavelet, winding, and parameter work.

These features should be the behavioral baseline for every refactor recovery
decision.

## Recovery Candidates

### 1. Recover the CLI and Regression Gates

**Priority: high**

Review and selectively recover the following from `~/refactor/gem-lte-core`:

```text
experiments/Feb2026/lte_run_cli.py
experiments/Feb2026/lte_config.py
experiments/Feb2026/lte_physics.py
experiments/Feb2026/lte_primitives.py
experiments/Feb2026/lte_solver.py
experiments/Feb2026/lte_types.py
experiments/Feb2026/check_lte_results_csv.py
experiments/Feb2026/test_amo.sh
experiments/Feb2026/expect_test_amo.exp
experiments/Feb2026/expect_test_amo_csv.exp
```

First adapt their input/output contracts to the current Ada model and
experiment layout.  The regression checks should establish a repeatable
baseline for later native-code refactoring rather than dictate the current
model behavior.

### 2. Recover Port Tracing as a Diagnostic Facility

**Priority: high**

The refactor copy's `PORT_TRACE`/`PORT_TRACE_N` support in
`src/gem-lte-primitives-solution.adb` writes bounded CSV snapshots for:

- input data;
- raw tidal sum;
- impulses;
- IIR forcing;
- post-Bessel forcing; and
- generated model output.

This is useful for comparing the current forcing pipeline against historical
results.  Reintroduce it behind its environment flag, ensure all output files
are closed on failure, and avoid changing results when tracing is disabled.

### 3. Evaluate `MF_SAME` Independently

**Priority: medium**

The refactor copy implements a distinct Bessel/forcing branch selected by
`MF_SAME`.  Recover it only as an explicit, opt-in experiment and compare it
against the current default forcing path using the recovered regression gates.

Do not make `MF_SAME` the default unless it demonstrably improves the selected
metrics without regressing the current baseline.

### 4. Reconcile the `Param_B` Data Model Before Refactoring

**Priority: high**

The refactor copy removed legacy fields:

```text
bg, ImpC, Ann1, Ann2, Sem1, Sem2, IR, Year
```

and added:

```text
dd, fbS, fbC
```

The current checkout still depends on the legacy model and adds behavior
around year correction, impulse processing, and uncompensated scoring.
Changing this record changes persisted parameter layout and can invalidate
saved `.p` files.

Before any field is removed or added:

1. Document every current read/write site and experiment that depends on it.
2. Define a versioned on-disk representation or an explicit migration path.
3. Add a compatibility test for representative saved parameter files.
4. Migrate one feature at a time, retaining old behavior as the default until
   results are compared.

`dd`, `fbS`, and `fbC` should not be restored merely because they existed in
the refactor copy; recover each only with a documented equation, call site,
and regression result.

### 5. Do Not Recover Unwired Random-Descent Experiments

**Priority: do not recover as-is**

The refactor copy declares and implements:

- `Phase_Lock` and `Adaptive_Phase_Lock`;
- `Enhanced_Nonlinear_Structure` and its adaptive variant; and
- combined and array variants of both.

There are no call sites outside those declarations and definitions.  They
therefore had no runtime effect.  Do not merge them until a concrete model
integration point, expected physical behavior, bounds, and regression cases
are defined.

### 6. Preserve Historical Results and Documentation Separately

**Priority: medium**

The refactor copy contains valuable but mostly non-executable artifacts:

- `docs/correlation_matrix.html` and PDF;
- LTE formulation documentation;
- site datasets and output files; and
- hundreds of `docs/cross/*.png` correlation images.

Recover these in a documentation/results change separate from native source
changes.  Do not mix generated images, result files, and code refactors in
the same commit.

## Suggested Execution Order

1. Keep `recovery/origin-master-before-force-push` and
   `~/refactor/gem-lte-core` available as comparison sources.
2. Establish the CLI/regression-gate baseline against the current checkout.
3. Add `PORT_TRACE` without changing default calculations.
4. Evaluate `MF_SAME` behind an opt-in flag.
5. Decide whether the `Param_B` layout needs a migration; implement that
   migration only after compatibility coverage exists.
6. Recover documentation and generated results independently.
7. Leave phase-locking and enhanced-nonlinearity helpers out unless a
   subsequent experiment supplies real call sites and measurable acceptance
   criteria.

## Refactor Acceptance Criteria

For every recovered capability:

- The default current behavior remains unchanged unless an intentional default
  change is documented.
- Environment flags have an explicit default and validation for invalid
  values.
- Saved parameter files remain readable or have a documented migration.
- Regression gates compare the current baseline, the recovered option, and
  any changed output artifacts.
- Code, experiment scripts, generated results, and documentation land in
  separate, reviewable commits.
