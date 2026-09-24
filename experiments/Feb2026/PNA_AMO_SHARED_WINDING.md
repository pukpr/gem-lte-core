# PNA's dominant winding (M ≈ 2.065) and its link to AMO

Written 2026-09-24, following up on a real optimizer run: `pna` (Pacific/
North American pattern, monthly index, 1950.0-2026.33, 917 months),
manifold search seeded directly from AMO's own converged `lt.exe.p`,
`EXCLUDE=TRUE TRAIN_START=1970 TRAIN_END=1985 TREND=TRUE F9=1`. The
`VALIDATE=TRUE` run's own held-out score for the excluded 1970-1985 window
was weak (0.145 — see the immediately preceding session exchange, where
two real Ada bugs in the new `VALIDATE` lockbox were found and fixed), but
the *training* fit itself (0.66, over the ~61 remaining years of record)
converged cleanly on one dominant term:

```
M = 2.065414   amp = 0.536644   (largest of 6 windings, ~2.7x the next)
```

This note asks: is that real, and if so, what would it mean for an
atmospheric index whose own dynamics are not well understood.

## 1. It replicates from a cold-ish restart

Re-ran the identical config (`EXCLUDE=TRUE TRAIN_START=1970 TRAIN_END=1985
TREND=TRUE F9=1`), `VALIDATE=FALSE`, fresh from the same AMO seed, 8
threads, only 4 minutes:

```
M = 2.064714   amp = 0.386755   (tied for largest, with M=1.606, amp=0.386779)
```

Same frequency, same harmonic slot (harm=9, i.e. `9 x k2` where
`k2 ≈ 0.2294` is `ltep[1]`), independently re-found from scratch in a
fraction of the time the original run took. `k2` itself barely moved from
AMO's own seed value in either run (AMO: `0.229413` -> pna original:
`0.229488` -> pna replicate: `0.229413`, i.e. back to AMO's own value
exactly) — the search keeps landing on the same carrier frequency AMO
already had.

## 2. It is independently confirmed in PNA's own raw data — not just the fit

`winding_rank.py` tests candidate winding numbers against the *data*
column (never the model), so this isn't circular: does `pna`'s own
monthly series really correlate with `sin(2*pi*M*Forcing(t))`, for the
current manifold's own phase trajectory `Forcing(t)`, more than an
AR(1)-matched noise floor would predict?

```
$ python3 winding_rank.py --index pna --m-max 12 --dm 0.01
M=2.07   peak=4.19 bits   FWHM=0.06   continuity=1.00   PASS   <- tied for #1 of 22 scanned ridges
M=10.93  peak=4.19 bits   FWHM=0.06   continuity=1.00   PASS   <- tied #1, unexplained (see below)
M=7.40   peak=4.03 bits   FWHM=0.08   continuity=0.92   PASS
...(6 more passing ridges, all M >= 4.06)
```

`M=2.07` — within one grid step (`dm=0.01`) of the optimizer's own
2.065 — ties for the single strongest, sharpest, most stationary ridge in
the whole 0-12 sweep. That is real, load-bearing evidence: it isn't merely
"the optimizer kept a term it inherited," the raw data itself, run through
an independent statistical test, picks the same frequency out as its best
candidate.

**Important caveat, stated plainly:** this test is run against the
*current* (AMO-inherited) manifold's own phase trajectory, not some
manifold-agnostic absolute spectrum of `pna.dat`. A `winding_rank` sweep
against `pna`'s own *earlier*, independently-fit (non-AMO-seeded) manifold
five days ago (see `project_pna_manifold_decomposition` memory) found a
completely different set of 13 ridges (1.46, 2.97, 4.73, 5.40, 5.81, ...),
with nothing near 2.065-2.07. So the honest reading is: **given AMO's own
nonlinear tidal/Bessel manifold as the "clock," PNA's data locks onto
M≈2.07 better than any other candidate in range** — not "PNA demands this
frequency independent of any manifold choice." Which manifold-phase basis
is the *right* one to test against is exactly the open question this whole
project keeps running into (see `WINDING_NARRATIVE.md`'s "shared reference
frequency across all seven" section) — this result is a new data point
for AMO's manifold specifically being a good clock for PNA too, not proof
of an absolute truth about PNA's spectrum.

## 3. AMO already carries the same term, independently

AMO's own currently-converged `lt.exe.p` (fit on real Atlantic SST, with
no PNA involvement) has:

```
M = 2.064720   amp = 0.631977   (2nd-largest of AMO's own 6 windings,
                                  behind only its k1≈0.016 near-DC mode)
```

This is not a coincidence of the seeding — AMO had this term at
comparable relative strength before `pna` ever borrowed its parameters.

## 4. But AMO's own raw data does *not* independently demand it

Running the same `winding_rank.py` check against `--index amo`:

```
0 of 22 candidate ridges pass in the M=0-12 range.
Nearest to 2.065: M=2.55 (2.42 bits, fails), M=2.91 (2.43 bits, fails).
```

So the asymmetry is real and worth stating precisely: **PNA's raw data
gives this frequency much stronger, cleaner, independent support than
AMO's own raw data does**, even though the fitted *parameter* lives in
both models. AMO's optimizer evidently found the term compatible enough
to keep (nonzero, real amplitude) as part of its own joint multi-term fit,
but it isn't AMO's own standout signal the way it is PNA's.

## Summary table

| check | result |
|---|---|
| PNA optimizer (23+ min, EXCLUDE 1970-85) | M=2.0654, amp=0.537 (dominant, 2.7x next) |
| PNA optimizer (4 min replicate, cold restart) | M=2.0647, amp=0.387 (tied #1) |
| PNA raw-data ridge test (winding_rank) | M=2.07, 4.19 bits, FWHM 0.06, cont 1.00 — tied best of 22 |
| AMO fitted manifold (independent origin) | M=2.0647, amp=0.632 (AMO's own 2nd-largest) |
| AMO raw-data ridge test (winding_rank) | nothing near 2.065 passes |

## Possible mechanisms

PNA (Pacific/North American pattern) is a 500-hPa geopotential height
teleconnection over the North Pacific/North America, driven by
extratropical Rossby wave trains — not itself a directly-forced tidal
quantity the way an ocean SST index might more plausibly be. That a
manifold built around long-period lunar/solar tidal beats (the whole
premise of GEM-LTE — see `feedback_lte_tidal_terms.md`'s Mm/Mf framing)
fits it this cleanly is exactly the kind of result that needs a mechanism,
not just a fitted number. Four candidates, in decreasing order of how much
direct support they have from mainstream literature vs. this project's
own internal findings:

**(a) A common Length-of-Day / atmospheric-angular-momentum pathway.**
This is the mechanism GEM-LTE's own architecture is already built to
check — every fit reports a `:dLOD:` correlation (`CompareRef` against a
reference LOD series) precisely because long-period tidal torques modulate
Earth's rotation (LOD), and LOD is mechanistically linked to atmospheric
angular momentum (AAM) on sub-decadal timescales in the mainstream
literature (mostly discussed there as an ENSO-AAM-LOD chain). PNA's own
circulation pattern is fundamentally a redistribution of atmospheric
angular momentum via the jet stream's position and amplitude. If the same
tidal beat that paces AMO's slow ocean-mixing-driven SST evolution also
imprints on the solid-Earth/ocean/atmosphere AAM budget, both indices
could inherit the same M≈2.065 harmonic from one upstream driver without
one causing the other. (This note deliberately does *not* cite this run's
own `:dLOD:` numbers as evidence for this — they were uniformly ~0.996-
0.999 across every region tested this session regardless of index, which
looks like an artifact of the shared manifold's dominant terms dominating
that comparison rather than a per-index confirmation; the mechanism is
plausible on independent grounds, the number itself isn't being used to
argue for it here.)

**(b) An AMO -> tropical Pacific -> PNA teleconnection chain.** This has
real support in the mainstream literature independent of GEM-LTE: AMO's
multidecadal SST state is documented to modulate tropical Pacific
SST/convection and Walker circulation strength (the Atlantic-Pacific
trans-basin variability literature), and PNA is one of the most robustly
established ENSO teleconnections (warm ENSO -> enhanced tropical
convection -> Rossby wave train -> positive PNA). A shared M≈2.065
harmonic could be propagating from AMO into PNA via this established
ocean-atmosphere-teleconnection pathway (through the tropical Pacific as
an intermediary) rather than via a direct, independent tidal imprint on
the atmosphere.

**(c) A genuine but as-yet-unidentified rectification chain.** `M=2.065`
is the 9th harmonic of `k2≈0.2294` — itself not one of the base
13-27-day tidal periods in the shared `lpap` table, but a *downstream*,
already-nonlinear modulation frequency (the same status `0.2073-0.2077`
has across seven other indices in `WINDING_NARRATIVE.md`). What specific
beat/rectification chain of the raw lunar/solar periods produces a
~0.229-cycle carrier, and why its 9th harmonic specifically, is not
worked out here — flagged as open, not asserted.

**(d) Shared optimizer attractor, not shared physics.** The most
skeptical read: every run so far has been warm-started from AMO's exact
parameters, and even the "replicate" run in section 1 never left that
basin. §2's raw-data ridge test mitigates this concern substantially
(it's a non-parametric test of PNA's own data, not of the optimizer's
search dynamics) but does not fully retire it, because it was still run
against a manifold phase that AMO's own fit produced. **The single
strongest remaining test is a genuinely cold-started PNA search** — random
initial conditions, no AMO seed at all — to see whether M≈2.065 is
rediscovered independently or whether it only ever appears when the search
starts near it. Not yet done; flagged as the natural next step.

## 5. Significance test + AIC/BIC on the windowed correlation (added 2026-09-24)

`pnasite1970-1985.png`'s bottom-left panel (50-month running windowed
correlation, Model vs Data) stays high and consistently oscillating
across the ~61 non-excluded training years, not just briefly. Two
rigorous checks on that, via `pna_windowed_corr_significance.py`:

**AR(1)-surrogate significance test** (same discipline as `winding_rank`'s
own null model, applied to the windowed-correlation trace instead of a
ridge scan): generate 2000 AR(1) surrogates matched to the real Data's own
lag-1 autocorrelation (`phi=0.757`) and variance, run each through the
identical 50-month rolling correlation against the SAME real (fixed)
fitted Model, outside the excluded 1970-1985 gap.

```
                          real value   surrogate mean (std)   z      percentile
mean windowed corr        0.603        -0.002  (0.065)        9.37   100.0%
worst-case (min) window  -0.065        -0.562  (0.092)        5.41   100.0%
```

Both the average level AND the single worst 50-month window over the
entire 61 non-excluded years clear the AR(1) noise floor by a wide margin
(z > 5 even in the worst case). This is about as unambiguous as this kind
of test gets: the sustained windowed correlation is real structure, not
an artifact of Data's own autocorrelation.

**Nested AIC/BIC** (full 6-winding model vs. a trend+annual+semiannual-only
baseline with zero windings, both refit via plain OLS on the identical
EXCLUDE=1970-1985 training subset, so this is a clean apples-to-apples
comparison rather than reusing the stochastic search's own amp/phase):

```
                    k    RSS       delta AIC (full - base)   delta BIC
raw N=737 months    22 vs 8   217.8 vs 329.7   +277.6 (decisive)   +213.2 (decisive)
eff N=117 months*   22 vs 8   217.8 vs 329.7    +20.5 (decisive)    -18.1 (favors baseline)
```
\* effective sample size after correcting for the baseline residuals' own
lag-1 autocorrelation (`phi=0.727`) via the standard AR(1) formula
`n_eff = n*(1-phi)/(1+phi)` -- 737 nominal monthly points are only worth
about 117 independent ones given how strongly autocorrelated this series
is. Naive (raw-N) AIC/BIC on monthly climate data is close to meaningless
overconfidence; this correction is not optional.

**Honest reading, not the more flattering half of it:** AIC (which weighs
prediction error more than model-selection consistency) endorses the full
model decisively either way. BIC (heavier, `k*ln(n)`, complexity penalty,
designed for consistent model selection) flips once the effective sample
size is used -- it mildly *prefers the simpler 8-parameter baseline* over
the full 22-parameter, 6-winding model. The RSS reduction is real and
substantial (34% lower residual variance, full vs. baseline) and the
correlation-trace significance test above is completely unambiguous -- so
the windowed correlation itself is not in question. What's genuinely
ambiguous is whether *this specific 6-winding parameterization*, with its
14 extra parameters (2 of them continuous/searched, 4 of them discrete
harmonic-index choices not even counted as continuous parameters above,
which if anything makes the baseline look relatively *better*), is the
most parsimonious way to capture that real structure, versus a lower-
order manifold. That's a real, open question this result raises rather
than answers -- worth testing directly by sweeping the number of
windings (e.g. NM=2/NH=0 through NM=2/NH=4) and checking where
eff-N-corrected BIC actually bottoms out.

## Open questions

- Run a cold-started (non-AMO-seeded) PNA optimization for long enough to
  reach comparable convergence, and check whether M≈2.065 is
  re-discovered from scratch.
- `M=10.93` ties for PNA's #1 raw-data ridge alongside 2.07 and is
  currently unexplained — no clean small-integer relationship to 2.065,
  to the existing harm=9 slot, or to `k1≈0.3435` was found. Worth its own
  look.
- Work out what beat/rectification of the primitive lunar/solar tidal
  periods produces a ~0.229-cycle carrier specifically (mechanism (c)),
  the way `project_amo_rectification_study.md` did for AMO's own
  120yr/60yr doubling.
- Test mechanism (b) directly: does adding a tropical Pacific
  (NINO4/NINO3.4) intermediary term change PNA's own winding structure or
  explain the M≈2.065 carrier's amplitude better than AMO's manifold
  alone?
- Sweep winding count (NM/NH) and re-run the effective-N-corrected
  BIC/AIC comparison at each step to find where BIC actually bottoms
  out — §5 above shows the full 6-winding model isn't obviously the most
  parsimonious once autocorrelation-corrected complexity is priced in
  fairly.
