# The QBO Formulation, So Far

This documents the Quasi-Biennial Oscillation (QBO) work already present in
this project — mostly built before today's session (Baltic/AMO/ENSO), and
not touched by it. Two distinct provenance threads exist and are kept
separate below: a theoretical derivation authored by a Gemini CLI session
(`qbo30/QBO_TIDAL_THEORY_VERIFICATION.md`, `GEMINI_VERIFICATION.md`), and a
set of independent verification/diagnostic scripts in this project's own
Python layer (`qbo_phase_lock.py`, `qbo_isolated_iir.py`,
`qbo_compensating_dtw.py`, `qbo_phase_sensitivity.py`). Every numeric result
below marked "verified today" was obtained by actually re-running the
corresponding script during this write-up, not copied from an old comment.

## 1. What QBO is, and why it fits this project's framework

QBO is the ~28-month (mean), quasi-periodic reversal of equatorial
stratospheric zonal winds — famous for a mean period that is neither fixed
nor obviously tied to any single forcing, and for a documented anomaly (the
2015–2018 disruption) where the regular alternation briefly broke down. The
dataset used here is Singapore 30 hPa zonal wind (`qbo30/qbo30.dat`), the
standard long equatorial QBO record.

Structurally, QBO is treated with the exact same architecture as the
AMO/ENSO/Baltic work: a calibrated tidal forcing manifold `L(t)` (column 4
of `lte_results.csv`, called `Forcing` elsewhere in this project) built from
the same Doodson-argument tidal periods, fed through the same
`Impulse_Delta` (twice-yearly Dirac comb) → `IIR` (leaky-integrator-with-
drag) → `Bessel`/sinusoidal-modulation pipeline as `gem-lte-primitives.adb`
already implements for every other index.

## 2. The astronomical basis for QBO's mean period (~2.37 years)

Authored by Gemini CLI (`qbo30/QBO_TIDAL_THEORY_VERIFICATION.md`), and
structurally the same mechanism as this session's own comb-alias study for
AMO (a fast tidal line aliased against a coarser sampling comb producing an
emergent slow beat) — arrived at independently for a different index.

Two separate derivations, converging on the same number:

- **Draconic month vs. the annual cycle.** The draconic month (27.2122d)
  fits into a year 365.2425/27.2122 ≈ 13.422 times. Stroboscopic (annual)
  sampling subtracts the integer part, leaving a beat frequency of 0.422
  cycles/year → beat period **≈ 2.3699 years**.
- **Synodic-draconic beat vs. the calendar month.** The ~29.4-day
  synodic-draconic beat, itself beat against the 30.4369-day calendar
  month, gives a beat frequency of 1.1587×10⁻³ /day → **≈ 2.363 years**.

Both land within 0.3% of each other and close to QBO's real documented mean
period (~28 months = 2.33 years) — a legitimate, checkable astronomical
coincidence-or-mechanism, in the same family as this session's own
Mt-vs-calendar-year (40:1) finding for AMO.

## 3. The regression form and the "calendar-month winding" mechanism

From `GEMINI_VERIFICATION.md`, the fitted model (`lte_results.csv` column
2) has the structure:

```
Model(t) = a·L(t) + b·sin(2π·c·L(t) + d)
```

— a direct linear pass-through of the forcing (`a·L(t)`, baseline
background flow) plus a winding/sinusoidal term, exactly the
`sin(k·M(t))` family this whole project is built on, with `L(t)` playing
the role of the manifold.

The distinctive QBO-specific piece is how `c` (the winding number) is set.
Rather than a single fixed constant fit once (as for AMO/Baltic/nino4),
`c` is derived from **the reciprocal of `L(t)`'s own local plateau slope**:

```
c = 1 / ΔL
```

so that the phase advances by exactly 2π (one full winding) per calendar
month along each plateau — a genuine "rotating reference frame" /
coordinate-demodulation trick: it doesn't claim the atmosphere cares about
the human calendar, it uses the fixed monthly sampling cadence to cancel
the high-frequency tidal carrier out of the argument, leaving only the slow
QBO envelope to resolve. Deviations of `c·ΔL` from exactly 1.0 (rather than
the mean value) are what carry the slow secular structure (e.g. 18.6-year
nodal-scale drift).

Three versions of `c` were compared, with a **1988–1998 holdout** withheld
from fitting for all three (from `GEMINI_VERIFICATION.md`):

| variant | `c` | train r | holdout r |
|---|---|---|---|
| fixed global (median plateau slope) | ≈18.64 | 0.763 | 0.762 |
| fixed global (optimized) | ≈18.68 | 0.830 | 0.791 |
| **dynamic, instantaneous** `c(t)` | tracks local `ΔL(t)` | **0.924** | **0.822** |

The dynamic version — where the winding number adapts to the forcing's own
local slope rather than being pinned to one global constant — gives both
the best fit and (with the fixed-global case as the honest baseline) real,
non-collapsing out-of-sample skill.

## 4. QBO's variable period: phase jumps, not a wobbling beat frequency

Also from the Gemini derivation: real QBO cycles vary from 20–36 months
(mean 28), while the astronomical beat above is a rigid 2.37 years. The
proposed reconciliation is not that the beat frequency itself wobbles, but
that **discrete phase-slips** occur at the sudden jumps in `L(t)` between
plateaus (`ΔΦ = 2π·c·ΔL_jump mod 2π`, empirically ≈185.6° — close to a full
phase reversal). A forward slip shortens the local cycle toward 20 months;
a backward/destructive slip lengthens it toward 36. This is presented as a
hypothesis with a clear mechanism, not independently re-verified today.

## 5. Independent verification layer (this project's own scripts)

Separate from the Gemini-authored theory above, four diagnostic scripts
exist that check specific, falsifiable pieces of it against the live Ada
model output. Re-run today for this write-up:

**`qbo_phase_lock.py` — is the oscillator actually phase-locked?**
The production Ada model (`LOCKF` mode) sets the base winding number
directly from the IIR filter's own drag parameter: `M(NM) = 1/(Decay·mP)`.
Current calibration gives `M(NM) = -26.609` (this has clearly been
re-optimized since the docstring's stale comment of 18.92 — the live value
is what's reported here). A standard PLL/circle-map diagnostic (unwrap the
winding phase, remove the linear trend, check whether the residual stays
bounded or grows) finds the residual **bounded** (std 11.86 → 10.70 rad,
first half vs second half — flat, not growing): the signature of a locked
oscillator, not a slipping one. The `mP`-only free-running rate is
**−2.22 years**, close to the real lunar-calendar beat of **2.37 years**
derived above. Realized full-cycle periods: model median 2.00 years vs.
real data median 2.33 years — same ballpark, model runs somewhat fast.

This is a direct, concrete link back to today's AMO work: `mP` is the same
Coulomb-friction-style IIR drag parameter (`copysign(lag_c, ...)`) that
turned out to be a sensitively-tuned amplifier for AMO's own comb-alias
resonance. Here it's doing double duty — QBO's production configuration
locks the winding number *directly* to that same friction parameter, one
concrete step further than AMO's case where friction only amplified an
already-existing resonance.

**`qbo_compensating_dtw.py` — does the second-largest tidal line explain
the slow multi-decadal wander? Result: no, not shown here.** Tests whether
the 9.1207-day Doodson term (aliased to ~22 years under monthly sampling)
shape-matches the phase residual's slow wander, via DTW against a
random-phase-shifted null. Re-run today: every candidate (raw 9.1207d
alone, the 9.1207d+9.1085d pair, their aliased/demodulated versions, and a
27.2122d negative control) scored **|z| < 1** against its own null — no
candidate beat chance. This is an honest open gap: the slow wander in the
phase residual is not yet explained by this specific candidate mechanism.

**`qbo_phase_sensitivity.py` — the 2015–2018 QBO disruption.** Tests how
much the model's reversal timing near the real, documented 2015–2018
disruption moves under tiny perturbations of the two dominant lunar
periods. Re-run today: a **0.20-day** shift in the 27.2122-day draconic
period, or a **1.3-day** phase shift, moves reversal timing by ~6 months —
extreme sensitivity, at perturbation scales far smaller than any real
ephemeris uncertainty. The script's own stated caveat, worth preserving
exactly: this cuts both ways. It's consistent with "the 2015–2018 anomaly
is unremarkable phase noise in a marginally-resolved system" (a monthly
sample of a ~29.4-day mechanism is barely coarser than the driving period
itself), but equally consistent with "this construction sits near an
ill-conditioned/bifurcation-like regime," which would itself be a
reliability concern independent of whether real lunar forcing is involved.
The sensitivity measurement doesn't distinguish between these readings.

**`qbo_isolated_iir.py` — methodological note, not a result.** Establishes
that pushing an isolated subset of Doodson terms through the real, *nonlinear*
`Impulse_Delta`→`IIR` pipeline is not the same as decomposing the real
42-term `Forcing` column (IIR is nonlinear in its drag term, so
`IIR(sum) ≠ sum(IIR)`). The isolated 2-term port validates at r=0.409
against the real Forcing column — informative on its own terms, explicitly
not a certified decomposition.

## 6. Daily-resolution extension (Singapore QVBO)

A separate, smaller thread (`fit_singapore_qvbo_by_lte_hidden.py`,
`evaluate_singapore_qvbo_multiscale.py`) fits the same
`u = A + B·x + C·sin(D·x+E)` form directly to *daily* Singapore 30 hPa wind,
using the monthly `L(t)` hidden forcing linearly interpolated to each day's
midpoint, then tests whether adding intra-month phase features (deviations
from the monthly anchor, zero at month boundaries by construction) improves
daily skill beyond the baseline monthly-resolution model. Not re-run today;
included here for completeness of what infrastructure exists.

## 7. Summary verdict

| Claim | Status |
|---|---|
| QBO's ~2.37yr mean period = draconic/annual tidal beat | Two independent derivations agree to 0.3%; same aliasing family as this session's AMO comb-alias finding |
| `sin(k·L(t))` regression with dynamic, slope-tracking `c(t)` | Real, holdout-validated skill (r=0.822 out-of-sample) |
| Production model's winding is a genuinely locked oscillator | Verified today — bounded phase residual, not slipping |
| Winding number tied directly to the IIR friction parameter (`mP`) | Confirmed today from live `lt.exe.p` — same friction mechanism family as AMO's Coulomb-friction amplifier |
| Variable QBO period explained by discrete phase-slips at plateau jumps | Plausible, stated mechanism — not independently re-verified today |
| 9.1207d term explains the slow multi-decadal phase-residual wander | **Not supported** — re-run today, no candidate beat its null (all \|z\|<1) |
| 2015–2018 disruption timing sensitivity | Extreme (re-verified today: ~6 months of timing shift from a 0.2-day period perturbation) — genuinely ambiguous between "ordinary phase noise" and "ill-conditioned construction," by the diagnostic's own design |

The QBO thread is, if anything, further along than the AMO/ENSO work in one
specific respect: it already has a real out-of-sample holdout validation
(1988–1998) with a materially better-than-chance result, and an explicit,
already-built phase-lock diagnostic distinguishing "genuinely resonant" from
"merely correlated." What it does *not* yet have is the moving-gauge /
resolved-fluid-dynamics treatment built today for Baltic/AMO/ENSO — that
extension hasn't been attempted for QBO.
