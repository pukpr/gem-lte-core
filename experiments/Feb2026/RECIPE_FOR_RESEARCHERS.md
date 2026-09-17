# A Modest-Compute Recipe: Calibrated Tidal Manifold → Climate-Index Winding Numbers

**Who this is for**: anyone with a laptop, an LLM coding assistant, and an
afternoon. This is not a pitch for a grant, a cluster, or a bespoke model.
It's a recipe — replicate it first, in a sitting, then decide for yourself
whether it's worth pursuing further.

## The claim, stated plainly

Seven climate indices spanning four ocean basins (NINO4, PDO, AMO, TNA,
NAO, IODE, and a Baltic Sea mean-sea-level aggregate) — each individually a
well-documented, decades-long modeling headache — can be fit with a
*low-degree-of-freedom* regression (a handful of base winding numbers, a
few integer harmonics of them, plus independent annual/trend terms) once
they're regressed against a **single shared forcing manifold**. That
manifold isn't fit to the climate data at all. It's the Laplace Tidal
Equation forcing basis, and it's calibrated once, up front, against a
completely different, independently measurable signal: Earth's delta-LOD
(length-of-day) record. See `WINDING_NARRATIVE.md` in this directory for
the full write-up, and `param_survey.py` / `geo_fingerprints.png` for the
underlying evidence — most importantly, that seven *independently*
optimized per-index searches converge on the same base winding frequency
(`M(NM)`) to within 0.2%, with no cross-index constraint in the fitting
procedure to force that.

This document isn't the evidence. It's the recipe for getting the evidence
yourself, and an argument for why doing so costs almost nothing.

## Why the compute/training bar here is close to zero

This matters enough to be explicit about, because it's the whole reason
this is worth other people's time:

- **No neural network is trained anywhere in this pipeline.** The
  "model" for each climate index is a small number of physically-motivated
  free parameters (winding numbers/wave-numbers, their amplitudes and
  phases, a handful of independent calendar terms) found by ordinary
  multivariate OLS regression, wrapped in a random-descent search over
  integer harmonic choices. It's classical numerical methods, implemented
  in Ada, running on a CPU.
- **The datasets are tiny.** Each climate index is a single time series —
  a few hundred to roughly fifteen hundred monthly or annual values. This
  is not big-data territory.
- **The search space per index is low-dimensional.** A handful of base
  modes, a handful of integer harmonic multipliers, and the regression
  coefficients that follow deterministically from them. Runs converge in
  minutes on ordinary hardware, not distributed training jobs measured in
  GPU-days.
- **The astronomical calibration step has ground truth.** `GEM.dLOD`
  (`src/gem-dlod.adb`) fits the lunar tidal constituents against the
  *measured* Earth delta-LOD record before that basis is ever reused for
  a climate index. That's a well-posed problem with a known, checkable
  answer — you can verify this step in isolation, on its own, before
  trusting anything downstream of it. That's what makes the whole approach
  falsifiable rather than a black box.
- **The LLM's role so far has been software engineering, not statistical
  learning on the climate data.** No model was fine-tuned or trained on
  any of this data. An off-the-shelf coding-assistant LLM was used to
  implement the Ada optimizer's regularization mechanisms (ridge
  regression, coverage-based cross-validation, enclosing/held-out-window
  validation), fix a bug in the harmonic search, build the Python survey
  and visualization tooling, and help write up the results — the same way
  it would help with any other software project. That's the entire "LLM
  training requirement": zero. Any researcher with access to a capable
  coding assistant already has everything needed to reproduce and extend
  this.

Put together: this is a weekend-scale replication, not a resourcing
problem.

## The recipe

**Step 1 — Calibrate the manifold against something you can check.**
Fit the lunar tidal constituents to measured Earth delta-LOD. Use the
**monthly (Mm, ~27.55d) and fortnightly (Mf, ~13.66d) lunar terms as the
primary basis** — not the 18.6-year nodal / 8.85-year perigee periods
that dominate the standard tidal literature. Those long apparent periods
are *emergent*: beat/envelope periods produced by Mm/Mf interacting with
the annual cycle (and with the slow lunar node/perigee precession rates),
not independent forcing terms in their own right. `src/gem-lte.ads`'s
`Doodson_Args` basis is built exactly this way — fast monthly/fortnightly
periods (`Draconic`/`Tropical`/`Anomalistic`) combined via Doodson
arguments, with the long periods appearing only as slow argument
coefficients riding on top. Get this step right and checked against real
dLOD data before doing anything else — it's the one part of the pipeline
with an independently verifiable answer.

**Step 2 — Reuse that manifold unmodified as forcing input.** Don't refit
the tidal basis per climate index. The whole point is that the manifold is
shared and externally calibrated, not tuned per target.

**Step 3 — Search for a low-DOF winding-number fit.** Regress the shared
manifold onto the climate index: a small number of base wave-numbers, a
small number of integer harmonics of them, plus independent annual/
semi-annual/trend terms fit directly on the calendar. Order 5-10 free
winding parameters total is typical here — see the per-index tables in
`WINDING_NARRATIVE.md`.

**Step 4 — Guard against overfitting cheaply, without adding compute.**
This project's optimizer includes a small toolkit for this that's worth
copying rather than reinventing: `RIDGE` (L2 penalty on regression
coefficients), `COVERAGE` (fit on a leading fraction of the training
window, score the whole window), `ENCLOSING`/`ENCLOSED` (fit only on data
*outside* a held-out interval, score on the held-out interval, the held-out
gap or held-out-plus-fit-region). None of these need more compute than the
plain fit — they just change what slice of data the accept/reject
criterion sees.

**Step 5 — Check for cross-index convergence before trusting anything
else.** If the same base winding frequency shows up across independently
optimized fits to unrelated indices, that's the signature the manifold is
real and shared rather than an artifact of a flexible-enough fitting
procedure. This is the single most important sanity check in the whole
recipe, and it's cheap: run the search per index, compare the base
`M(NM)` values.

**Step 6 — Only then, interpret the differences.** Once the shared
backbone is established, the *differences* between indices (which
harmonics are energetically significant, how flat or peaked the spectrum
is, geographic clustering of similar signatures) become physically
readable — see `WINDING_NARRATIVE.md` and `geo_fingerprints.png` for what
that looks like in practice.

## Why this suits a multi-agent LLM effort specifically

This isn't a generic "AI could help with science" claim — the shape of
this particular problem happens to line up unusually well with what
multi-agent LLM coordination is actually good at:

- **The task decomposes into independent, parallelizable units with
  objective, checkable success criteria.** One agent per climate index,
  one dedicated to the dLOD calibration/verification step, one doing
  cross-index synthesis, one on physical/geographic interpretation — each
  with a quantitative score (correlation coefficient, DTW distance,
  held-out validation score) to report back, not a subjective judgment
  call. Coordinating agents is easy when "did this work" has a number
  attached.
- **The iteration loop is fast enough for many agents to explore variants
  cheaply.** A single index's search converges in minutes on ordinary
  hardware. That means an agent swarm can try different harmonic search
  ranges, different regularization settings, different held-out windows,
  or entirely different indices/basins in parallel, and compare results
  head-to-head, without anyone waiting on a training run.
- **Verification is built into the method at two levels.** The
  astronomical calibration step (Step 1) has ground truth an agent can
  check before trusting anything downstream — a natural gate any
  replication effort should hit first. The climate-index fits (Step 3-4)
  have explicit held-out validation built in — a natural objective for
  agents to optimize against or debate over, rather than in-sample fit
  quality alone.
- **Everything is open and inspectable.** This is a small Ada codebase
  plus Python tooling, not a proprietary model. Any agent with code-reading
  and execution tools can pick this up directly — there's no API to
  request access to, no weights to download.

## The actual ask

Replicate first. Build the Ada tool (`lte.gpr`), run the `GEM.dLOD`
calibration and check it against measured delta-LOD, pick one climate
index, run the search, and see whether your fitted base winding frequency
lands anywhere near the `M(NM)≈0.207` this project keeps finding across
seven independent runs. That's an afternoon, not a proposal. If it
replicates, the natural next steps are obvious and cheap to parallelize:
more indices, more basins, alternative regularization schemes, and a
harder look at tying specific winding numbers to specific named tidal
constituents beyond Mm/Mf. None of it requires more than what's already
sitting in this repository.
