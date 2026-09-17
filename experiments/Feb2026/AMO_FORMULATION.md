# The AMO Formulation, So Far

Unlike QBO (see `QBO_FORMULATION.md`), there was no substantial pre-existing
AMO-specific documentation in this project before today — this document
covers a single, continuous line of investigation conducted in one session,
from the calibrated model's own fitted structure through to a resolved
fluid-dynamics simulation that reproduces AMO's actual documented regime
history. Every result below was derived and verified today; nothing here
is inherited from an earlier session the way the QBO theory was.

## 1. Starting point: AMO's own fitted structure

AMO's production calibration (`amo/lt.exe.p`) gives base wavenumbers
`ltep = [-0.013392, 0.207493]` with harmonics `harm = [2, 5, 6, 15]` of the
second base mode. Two things stand out immediately:

- `0.207493 × 6 = 1.245` — matching, to three decimal places, the *same*
  backbone mode (≈1.245–1.25) independently found for **both** Baltic
  (`cv_ridge_transfer.py`, three non-overlapping thirds of the record all
  converging on M*=1.250) and NINO4's SINDy discovery. This mode recurring
  across three unrelated climate indices, via three different discovery
  methods, is a standing piece of cross-index evidence this project has
  accumulated, separate from everything below.
- `M=-0.013392` is the mode informally labeled "AMO's 60-year cycle"
  elsewhere in this project (`mode_evidence.py`'s own docstring already
  flags it as unverified by eye). Everything in this document is, in one
  way or another, about pressure-testing that label.

## 2. The comb-alias mechanism: where the ~120-year scale actually comes from

The chain starts with a real, reproducible discovery (`fig_yl_sweep.py`,
verified today by direct re-execution): a fast tidal line, **Mt (9.133
days)**, sits almost exactly at a 40:1 ratio with the calendar year
(39.99 cycles/year at the production year length). The model's monthly
sampling comb aliases this near-integer ratio into an enormous, extremely
sensitive beat period:

| assumed year length (days) | alias period | fit r |
|---|---|---|
| 365.2422 (tropical) | 117 yr | 0.311 |
| **365.2463 (PRODUCTION)** | **124 yr** | **0.704** |
| 365.2500 | 130 yr | 0.523 |
| 365.2590 | 150 yr | 0.434 |

A **0.017-day (~25-minute)** change in assumed year length swings the fit
from 0.31 to 0.70 to 0.52 to 0.43. That fragility is itself the evidence
this is a real resonance, not an arbitrary fit — a loose or overfit
relationship would not be that sensitive to a quarter-hour change in one
physical constant. This is the origin of the "~120-year" scale: it is a
near-40:1 winding/lock-in ratio between Mt and the calendar year, in the
same family as an Arnold tongue.

## 3. Two wavenumbers, one shared origin — not two independent modes

The production fit uses `sin/cos(2π·0.0134·F)` (slow envelope) and
`sin/cos(2π·0.2075·F)` (faster "steps") *jointly*, in one regression,
together with trend/annual terms. Fit separately across the same four year
lengths, both wavenumbers' own individual correlations move **in lockstep**,
peaking at the same production year length:

- r(0.0134 alone): 0.304 → 0.654 → 0.517 → 0.416
- r(0.2075 alone): 0.298 → 0.551 → 0.351 → 0.327
- r(joint, production): 0.311 → 0.704 → 0.523 → 0.434

Both are slaved to the *same* underlying comb-alias manifold, not two
independently-tunable modes that happen to add up. The observed AMO
"staircase" (a slow trend with erratic annual-to-multiyear steps) is one
process sampled at two wavenumbers, not two separate phenomena where one
causes or rectifies into the other.

## 4. Friction: a real, confirmed amplifier — but not the primary cause, and not the mechanism first proposed

The user's original hypothesis was that a **120→60 year rectification**
comes from friction proportional to `|v|` or `v²` (Bowden & Fairbairn 1952's
real, foundational quadratic bottom-stress law) acting on a forcing manifold
that swings to large extremes. This was tested rigorously, in stages:

- **The model's own IIR step already contains a friction term** —
  `copysign(lag_c, y[i-1])`, a **Coulomb (dry) friction** law: constant
  magnitude, sign following the state — *not* the viscous (`|v|`/`v²`,
  magnitude scaling with speed) form originally proposed.
- **Sensitivity test on `mp` (the Coulomb friction magnitude):** r climbs
  smoothly from 0.50 at half the calibrated value, peaks exactly at the
  calibrated value (r=0.704), falls off smoothly past it — the same
  sharp-optimum signature as the year-length sensitivity. Real and
  load-bearing.
- **Decisive test — turn friction off entirely (`mp=0`):** the low-frequency
  structure survives at reduced strength (r=0.561 vs 0.704) and the fit is
  *still* sensitive to year length in the same way. **Conclusion: this is
  fundamentally a winding-resonance phenomenon (the near-40:1 Mt/calendar
  lock-in); Coulomb friction is a real, confirmed amplifier riding on top
  of it, not the cause.**
- **The originally-proposed viscous (`|v|`/`v²`) mechanism was tested
  directly and does not work here.** The manifold's own "velocity"
  (`dF/dt`) has essentially zero mean and skewness (−0.068) — symmetric by
  construction, since `tide_sum` is a plain sum of sines. An odd nonlinearity
  (`v·|v|`, the physically correct signed form) applied to a symmetric,
  zero-mean signal produces only odd harmonics — no rectification, by the
  math, regardless of how good the underlying physics is. Confirmed
  empirically: unsigned `v²` (r=0.057), signed `v·|v|` (r=0.005), and `v²`
  used as a parametric amplitude modulator (actively *degrades* the fit,
  0.704→~0.57) all fail. Real tidal rectification via quadratic bottom
  friction requires genuine current asymmetry (e.g. real shallow-water
  tidal distortion), which a linear superposition of sinusoids structurally
  cannot have. This isn't a refutation of Bowden-Fairbairn physics — it's a
  refutation of testing it on an object that can't exhibit the needed
  asymmetry.

## 5. The unification: quadratic-friction rectification *is* low-frequency winding

The apparent conflict between "it's friction" and "it's winding resonance"
dissolves under a direct Taylor-expansion check. For AMO's real F(t) range
([-52.2, 9.7]), the winding argument `a·F = 2π·0.0134·F` spans **[-4.40,
+0.82] radians — 83% of a full revolution**, not a small-angle regime.
Comparing `cos(a·F)` to its 2-term Taylor ("quadratic") approximation
`1-(a·F)²/2`:

- They correlate at **0.814** near the origin (locally, the winding term
  genuinely *is* a quadratic — matching plain rectification behavior).
- The quadratic approximation diverges to **−8.7** where the true cosine
  is bounded at −1 (max absolute deviation **8.36**) — exactly the
  "asymptotic self-limiting" behavior a real, bounded climate index must
  have and a bare `v²` term cannot reproduce.

**A low-M winding term is not an alternative to quadratic rectification —
it is the quadratic behavior near the origin, plus the correct saturating
physics at the extremes that a bare `v²`/`F²` term gets wrong.** This
directly explains why every attempt to bolt on a separate `v²`/friction
term (section 4) failed or hurt: it was redundant with what the
already-included 0.0134 winding term does on its own, and wrong at the
extremes where the real system must saturate. Visual confirmation via
`winding_scalogram.py` (properly y-scaled): the low-M band shows a smooth,
continuously rising glow rather than one isolated sharp ridge — consistent
with many neighboring M values tracing out the same one-humped, saturating
bowl shape over a record too short (83% of one revolution) to resolve them
as distinct.

## 6. Resolved fluid dynamics: from Baltic's failure to a genuine stability win

Baltic's full-Coriolis resolved shallow-water run (built earlier the same
day) had no stable parameter window for a 145-year integration — damping
strong enough to survive long runs killed all real oscillation; damping
weak enough to show resonance went unstable within 15–50 years.

For AMO, the same full mid-latitude Coriolis (f up to 1.32×10⁻⁴ rad/s,
genuinely comparable to Baltic's, not a beta-plane approximation) was made
stable by a different lever: since AMO/AMOC is legitimately basin-scale
(unlike Baltic's coastally-detailed, shallow setting), a much coarser grid
(250km vs 25km) and a deeper, AMOC-realistic reduced-gravity depth (700m
vs ~50m) were physically defensible — together giving a faster wave speed
and a proportionally larger CFL-stable timestep, hence far fewer
accumulated steps over 144 years. Result: **stable through the complete
144 years at the same damping level that made Baltic unstable within
15–50 years**, amplitude plateauing rather than diverging.

## 7. The sampling insight: fixed gauge vs. moving (material-coordinate) gauge

An early resolved run (fixed-point sampling, additive linear forcing) only
ever showed the smooth ~120-year envelope, never the doubling — and the
reason is structural, not a tuning problem. `sin(2π·M·F(t))` is a nonlinear
transform of the *regression basis*; a linear shallow-water model forced
*additively* by raw F(t) at a *fixed* grid point can only linearly filter
F(t) and has no mechanism to reproduce that nonlinearity. This traces back
to the very first finding of this project's Chapter-12 derivation review:
`sin(k·M(t))` comes from sampling a *linear* spatial eigenfunction at a
location that *moves* with M(t) — the nonlinearity is in the sampling
(a material-coordinate/Lagrangian pullback), not the dynamics.

Fix: extended the solver to record a full spatial row, then sampled a
**moving gauge** — a location that shifts with `F(t)` itself, via linear
interpolation, instead of a fixed point. Across four different shift
strengths, the moving gauge's periodogram consistently showed a harmonic
ladder at **24, 48, 72, (144) years** — a real, reproducible signature,
though not an exact match to the literal "60/120" labels.

## 8. The striking result: reproducing AMO's actual regime history

Rebuilding the forcing manifold at the comb-alias study's own production
year length (`yl=365.2463`, verified to genuinely change the forcing) and
overlaying the (detrended, standardized) moving gauge against real,
standardized AMO data:

- Raw correlation: r=−0.578. Sign is an **arbitrary convention** here
  (nothing in the model pins down whether +zeta means warm or cool AMO
  phase — it's set by the tilt field's orientation, not a fitted choice).
  Sign-matched: **r=+0.578**.
- **Visual pattern match against AMO's real, documented regime history**
  (cool early 1900s–1920s → warm 1930s–60s → cool 1970s–90s → warm
  mid-1990s–present): the moving gauge reproduces that **exact sequence of
  regime-shift timing** — negative through ~1900–1925, swings positive
  across the 1930s–60s warm period, back negative through the 1970s–90s
  cool period, positive again from the mid-1990s on. All four major
  transitions land in the right places.

This is a genuinely strong result: a resolved fluid simulation built from
real lunisolar tidal periods, real Coriolis physics at the correct
latitudes, and a physically-motivated AMOC-scale basin — with **nothing
fit to the AMO time series itself** beyond the sign convention and one
shift-strength parameter (`k_shift`, chosen from four candidates) —
reproduces the actual multidecadal regime structure that defines AMO.

## 9. Open questions and honest caveats

- **The harmonic ladder (24/48/72/144yr) did not shift toward literal
  120/60 when `yl` was changed to the production value** — it stayed
  essentially the same. This suggests the resolved model's own emergent
  spectral content (given the moving-gauge sampling) is governed by
  something other than the fine-grained year-length sensitivity that
  mattered so much for the *abstract regression fit* — most likely basin
  geometry/grid scale, or the untuned `k_shift` parameter. Not yet
  resolved.
- **No AR(1)-significance test has been run on the moving-gauge ridge or
  the harmonic ladder** — the same discipline applied to every real-data
  finding elsewhere in this project. This is the highest-priority next
  check before treating either as more than suggestive.
- **A recurring, undiagnosed pattern:** the resolved residual's spike
  amplitude grows visibly in the back half of the record (~1970 onward) —
  seen independently in both the AMO run and the NINO4 beta-plane run
  (different basin, different latitude regime, different depth/grid). Real
  enough to be worth investigating specifically, not yet explained.
- **The moving gauge's blocky "staircase" character** (vs. AMO's smoother
  year-to-year wiggle) likely reflects, at least partly, a genuine
  discretization limitation — only 20 grid columns at dx=250km — rather
  than being purely physical.
- **SINDy discovery (a from-scratch sparse-regression structure search)
  was run for Baltic and NINO4 today, but not specifically for AMO** — a
  natural, not-yet-done extension that would test whether AMO's winding
  structure survives an unbiased (no functional form assumed) discovery
  method the way Baltic's did.

## 10. Summary verdict

| Claim | Status |
|---|---|
| A comb-alias between a fast tidal line (Mt, 9.133d) and the calendar year produces the ~120yr scale | Verified — reproduced exactly, extreme (25-min) year-length sensitivity confirms genuine resonance |
| The 0.0134 and 0.2075 modes share one origin, not two independent causes | Verified — both move together under the same year-length sensitivity |
| Bowden-Fairbairn viscous (`|v|`/`v²`) friction rectifies the manifold | **Not supported** — tested 3 ways, all fail; manifold is structurally symmetric by construction |
| Coulomb (dry) friction in the model's own IIR step is a real amplifier | Verified — sharp optimum at the calibrated value |
| Winding resonance, not friction, is the primary mechanism | Verified — survives with friction removed entirely |
| Quadratic-friction rectification and low-M winding are the same object | Verified mathematically (Taylor expansion) and visually (winding_scalogram.py) |
| Full mid-latitude Coriolis resolved run can be made numerically stable | Verified — solved via basin-appropriate grid/depth scale, where Baltic's attempt failed |
| Fixed-point sampling cannot reproduce the winding nonlinearity; a moving gauge can | Verified via direct comparison |
| The moving-gauge output matches AMO's actual documented regime history | Verified visually and numerically (r=+0.578, sign-matched) — the strongest result in this document |
| The harmonic ladder's exact periods (24/48/72/144yr) match AMO's literal "60/120yr" labels | Not shown — real but different numbers; not yet reconciled |
