# Simulated referee pass: GEM-LTE cascaded/layered formulation (2026-10-06)

Protocol: `skills/skills/referee-sim/SKILL.md`, run in a fresh agent context (no authoring history).
No existing file was modified. Section 6 of the protocol (defensive edits) is given as PROPOSED
edits only, because the instruction was not to modify anything.

The "document" audited is the claim as stated for review:

> AMO, PDO, NINO4 and 89 Kaplan 20x20-degree quads are all described by one cascaded formulation
> model = level + k0 L + trend + annual + sum_j a_j sin(2 pi k_j L + phi_j), then F9, with
> L = A sin(2 pi k_f (F + B sin(2 pi k_b F + phi_b)) + phi_f), where F is the lag-integrated,
> annually impulsed lunisolar forcing whose ~42 constituents are calibrated to daily dLOD
> (gate CC >= 0.99) and "essentially locked". Shared: F, k_f, k_b, A, phi_f, B, phi_b. Per series:
> 1-2 windings (k ~ 4; ~6.2 for NINO4) plus harmonics and linear coefficients, about 7 nonlinear +
> 13 linear parameters. Evidence: low-winding level structure ("LTE Modulation").

together with its supporting files: DECLARATION*.md, SURROGATE_RESULT.md, levels_result.txt,
year_y1.json, year_y2.log, quads_frozen.json, Oct2026_layer_ir/GATE*.md, and the code in `src/`.

Objections the author already raised (not repeated as headline objections here): data-side phase
scrambling tests timing only weakly; high windings must be penalized; tidal parameters are locked
by LOD; moving AMO's slow phase is weak with few cycles; YEAR detuning is a cheap model-side null.
Several objections below show that two of these defences (LOD lock, YEAR null) are weaker than the
author believes, for reasons the author has not stated.

---

## 0. Checks run for this pass (new numbers, all reproducible from seeds/)

All computed from `seeds/{amo,pdo,nino4}/lt.exe.p`, `lte_results.csv`, `dlod_ref.dat`.

C1. **What the dLOD gate can see.** The gate (`gem-lte-primitives-solution.adb` lines 75-100) is
the CC between two *model* tidal sums on the dLOD dates: Tide_Sum with constituents re-regressed
on dlod_ref.dat, and Tide_Sum with the SST fit's constituents. It is not a fit to the dLOD data
of the SST fit's constituents. Using AMO's 42 constituents on the 1962-2019 dLOD dates:
| perturbation of AMO's constituent set | tidal-sum CC vs unperturbed |
|---|---|
| delete all 11 constituents with period > 200 d (18.6 yr, 8.85 yr, 9.3 yr, 4.4 yr, 3 yr, 6 yr, Ssa, ...) | 0.9968 |
| delete all 24 constituents with amplitude < 0.01 | 0.9983 |
| replace AMO's long-period terms by PDO's | 0.9982 |
| multiply long-period amplitudes by 3 and add 1.5 rad to their phases | 0.967 |
| shift the whole sum by exactly 1 day | 0.890 |
The gate is therefore blind to all decadal and all minor constituents. It is very sensitive to a
1-day timing error in Mf/Mm, and only to that.

C2. **Per-series values of the "locked" constituents differ.** Major terms agree across the three
fits (Mf 13.6608 d: amplitude 0.296/0.306/0.303, phase 0.638/0.648/0.638). Minor and long terms do
not: 1616.3 d phase -4.00/-0.75/-3.75; 2190.35 d amplitude 0.0018/0.0109/0.0088; 6167 d
amplitude 0.0031/0.0013/0.00003; 1656.3 d phase 0.08/-0.79/2.89; 27.3217 d amplitude
-0.00001/0.0086/0.0034. The set also contains 1095.175 d (3.0 yr) and 2190.35 d (6.0 yr), which
have no obvious Doodson argument.

C3. **The reference dLOD has no low-frequency content.** Power of dlod_ref.dat: below 0.2 cyc/yr
0.1%, 0.2-1.5 cyc/yr 0.5%, 20-40 cyc/yr 87%. It is a band-passed (or modelled) tidal series, so it
cannot constrain the 3-18.6 yr terms even in principle.

C4. **The latent is decadal.** The power of F2 (lte_results.csv column 4) lies 84-86% at periods
> 8 yr, 10-12% at 2-8 yr, ~4% below 2 yr; lag-12 autocorrelation 0.79-0.83.

C5. **"Cross" manifolds are near-copies.** Correlation of F2 between fits: AMO-PDO 0.971,
AMO-NINO4 0.971, PDO-NINO4 0.991. Yet S1 drops from own to cross manifold by ~2x (AMO 0.456 to
0.366, PDO 0.263 to 0.113, NINO4 0.212 to 0.105).

C6. **F2 bins are time-segregated.** Correlation ratio of (t - tbar)^2 by the 12 equal-count F2
bins used in S1: AMO 0.29, PDO 0.19, NINO4 0.15 (of t itself: 0.14/0.06/0.08). Bins of F2 partly
encode *where in the record* a month lies, so S1 is partly sensitive to any slow non-linear
residual trend left after linear detrending.

C7. **YEAR is a time lag, not a year length.** Tidal phase is 2 pi (Y/P) t with t in calendar
years from year 0 (Tide_Sum, gem-dLOD). Changing Y by delta days shifts every constituent by
delta x t days, i.e. about delta x 1950 days in this record (varying by only delta x 140 days
across it). The fitted YEAR offset +0.00405 d (Y = 365.2463 d, which is no astronomical year) is
equivalent to a common ~7.9-day lag. Y1's delta = 0.0005 d is a ~1-day shift; Y1's dLOD of
0.876-0.886 there matches C1's exact 1-day shift (0.890).

C8. **Sampling.** The tidal sum is evaluated at the instants of the monthly time stamps, then
multiplied by a two-slot comb (`Impulse_Delta`: DelA at month slot DelB, Asym at the slot half a
year later; DelA = -4.23, Asym = +6.24 for all three fits), so F is driven by the tidal value at two
instants per year, which are set by the data's time-stamp convention.

C9. **Y2 (equal-budget re-search) is complete:** delta = 0 wins P for all three, but every
delta != 0 run ends with dLOD between -0.80 and -0.15. The search never re-pinned the tidal
phases, and the delta = 0 run starts at its own long-searched optimum.

---

## 1. Audiences

| audience | why they get a row |
|---|---|
| A. Time-series statistician | methods borrowed: surrogate testing (IAAFT), BIC, CV, correlation ratio |
| B. Geodesist / Earth-rotation (LOD) | borrowed authority: "calibrated against dLOD, CC 0.99, locked" |
| C. Physical oceanographer | borrowed theory: Laplace tidal equations; a mechanism claim for SST |
| D. Climate-variability community (AMO/PDO) | incumbent explanations: forced response and stochastic red noise |
| E. SST-dataset community (Kaplan, HadSST, ERSST) | dataset borrowed: Kaplan quads, pre-1950 back-cast target |
| F. ENSO forecasters (practitioners, incumbent) | would act on forecasts; the 2026 El Nino miss |
| G. Tidal-climate literature (incumbent on novelty) | Keeling & Whorf 2000, Ray 2007, Yasuda 18.6-yr, Loder & Garrett |
| H. General reader | sees only the headline claim |

---

## 2. Objections (steelmanned, referee's voice), two axes per audience

### A. Time-series statistician

**A1 (correctness, HIGH): the parameter count is not the capacity, and the project's own file
contradicts the headline count.** "You quote ~7 nonlinear + ~13 linear parameters. Your own
SURROGATE_RESULT.md says the full search has ~150 free parameters per series, and that at that
freedom the model fits phase-scrambled surrogates almost as well (Arm B failed, 1 of 3). The small
count belongs to a restricted estimator that you have not yet built or tested. Beyond that, a count
of parameters that enter through phase (sin(k x), sin(2 pi (Y/P) t)) bounds nothing: sin(k x) with
one parameter has infinite VC dimension. Your Y1 result shows the fit is destroyed by a 1-day
shift, i.e. some parameters are effectively specified to ~1e-6 relative precision, which in MDL
terms costs many bits per parameter." What they know: GDF (Ye 1998) and MDL measure complexity by
sensitivity of the fit to the data, not by counting.
*Decisive test:* build the declared restricted estimator (six shared inner parameters searched
jointly across AMO/PDO/NINO4; tidal constituents frozen, see B1; outer sums by regression) and run
Arm B with it on surrogate triplets, cold and cross-seeded, equal budget. Also report generalized
degrees of freedom: perturb the data by small i.i.d. noise, re-run the full search, and estimate
sum d(yhat)/d(y). If GDF >> 20 the parsimony claim fails as stated.

**A2 (correctness, HIGH): the "cross-manifold" test is not out-of-sample, and the own-manifold
results are circular.** "The cross manifolds correlate 0.97-0.99 with each target's own manifold
(C5), and the shared inner parameters were developed on this same family of series. That the
remaining 1-3% of manifold difference halves eta2 (C5) is the textbook signature of in-sample
tuning. Your S2 'BIC picks k ~ 4 on all real indices' holds only on the own manifold, which was
co-tuned with a k ~ 4 winding: it is circular. On the cross manifold the BIC picks 2.84 / 3.62 /
3.68, not 4." What they know: a selection-free holdout must withhold the data from every stage of
fitting, the latent included.
*Decisive test:* temporal holdout of the whole pipeline. Fit everything (inner layer, tidal
amplitudes and phases, k) on 1950-2022 only. Freeze it. Compute S1 and the k-selection on
1880-1949 only. Repeat identically for each IAAFT surrogate (fit its 1950-2022 part, score its
1880-1949 part). Rank the real series.

**A3 (correctness, MEDIUM): S1 is partly a slow-trend detector.** "F2 is 85% decadal (C4), and
its bins encode position in the record (C6: eta2 of (t - tbar)^2 by F2 bins = 0.29 for AMO).
After linear detrending, AMO keeps a strong curvature (warming acceleration, 1960-90 aerosol
dip). Twelve free bin means on a decadal, time-segregated covariate can capture that. IAAFT
preserves the spectrum but scatters the timing of that curvature, so the real series wins for a
reason unrelated to tides." What they know: the effective sample size of a decadal statistic over
142 yr is ~15-20, and a nuisance shared by the data and the covariate must be removed first.
*Decisive test:* re-run S1 after regressing out a quadratic, or better the forced response (see
D1), and also within 30-yr blocks (bin edges fixed globally, eta2 pooled within blocks), so that
across-block level differences cannot contribute.

**A4 (correctness, MEDIUM): forking paths in the 2000-05 window.** "GATE.md declares the 2000-05
window contaminated (it entered accept rules). SURROGATE_RESULT.md still lists 'NINO4 test window
ranks first in Arm B' as a residual positive signal." *Fix:* drop it as evidence. Use only
1880-1949 (as in A2) or data after 2026-10.

**A-novelty:** "Spectrum-preserving surrogates are standard. What is new is the claim of a
low-dimensional latent driver. That claim must be compared with the generic alternative: a
decadal latent of any origin (e.g. a smoothed random walk or a CMIP forced response)
through the same 12-bin / one-winding readout." *Test:* replace F2 by 19 random decadal
latents with F2's spectrum (IAAFT of F2 itself) and by the CMIP6 forced response. Run S1/S2 on the
real data with each. If a generic decadal latent scores like F2, the tidal origin adds nothing.

### B. Geodesist / Earth-rotation (LOD) expert

**B1 (correctness, HIGH): "LOD-locked" does not hold for the constituents that make the
decadal latent.** "Your gate is a CC between two tidal *model* sums, dominated by Mf/Mm (87% of the
reference's power is at 20-40 cyc/yr, C3). Deleting every constituent longer than 200 d, or every
constituent with amplitude below 0.01 (24 of 42), leaves the gate at 0.997-0.998 (C1). The per-series
values of those constituents do differ, by up to 6x in amplitude and ~3 rad in phase (C2). Yet F2,
the object your evidence rests on, is 85% decadal (C4). So the parameters that control the
evidence are free per series and invisible to the gate. Also: the 18.6-yr zonal tide in LOD is not
small (IERS Conventions 2010, Table 8.1 zonal model: roughly a third of Mf's LOD amplitude);
your fitted 6793-d amplitude is ~1-2% of Mf's. Either your amplitudes are not LOD-calibrated in
the long-period band, or the reference was high-passed (C3) and cannot speak to it." What they know:
the IERS zonal-tide model (Yoder et al. 1981; Defraigne & Smits 1999) fixes the long-period LOD
amplitudes theoretically, and the observed LOD below ~1 cyc/yr is dominated by atmospheric angular
momentum (AAM) and the core, not tides.
*Decisive test:* (i) freeze all 42 constituents at the values of a single dLOD regression, or at
the IERS Table 8.1 values scaled to the same Mf amplitude, identical for all series. Remove
non-Doodson terms (1095.175 d, 2190.35 d) unless each is given a Doodson argument. Refit only the
inner layer and the regression. (ii) Attribute F2's > 8 yr variance by constituent group (major
zonal / minor / long-period / non-Doodson) by deleting each group. If (i) loses most of the skill,
or (ii) shows that the decadal latent comes from the gate-invisible groups, "locked" is withdrawn.

**B2 (correctness, MEDIUM): YEAR is a free ~8-day lag disguised as a calendar constant, and
the gate's sensitivity is only date bookkeeping.** "With the phase origin at year 0, YEAR acts as
a common time shift of delta x ~1950 days (C7). Your fitted 365.2463 d matches no tropical,
sidereal, anomalistic or Julian year; it is a fitted ~7.9-day group lag. The Y1 'detuning' is a
1-to-100-day shift of the tide against the data. That dLOD collapses under it (0.99 to 0.89 at 1
day) shows that Mf's phase is pinned to ~1 day. It says nothing about SST." What they know:
standard practice computes arguments from Delaunay/Doodson variables at J2000 (Simon et al. 1994)
with the standard year, and fits any lag explicitly, in days.
*Decisive test:* reparameterize. Use the standard tropical year, J2000 origin, and an explicit
lag parameter tau (days) applied to the forcing. Confirm the fit is identical. Then report tau
with its uncertainty and its physical meaning. Rerun Y1/Y2 as "tau shifted by N days", with the
tidal phases re-derived from dLOD at each tau.

**B3 (correctness, MEDIUM): provenance of the dLOD reference.** "dlod_ref.dat has almost no power
below 1.5 cyc/yr (C3). Is it observed IERS C04 LOD after filtering, or a tidal model (IERS
zonal)? If it is a model, '0.99 against dLOD' is agreement with a table, not with data. If it is
observed and filtered, state the filter, and show the CC against the *unfiltered* observed LOD in the
same band, with the AAM-corrected residual."
*Test:* document the source and filter. Report the gate against IERS C04 minus AAM (NCEP/ERA) in
the 5-40 cyc/yr band.

**B-novelty:** "Fitting zonal tides to LOD at CC ~0.99 is routine (the standard tidal LOD
correction). It is not, by itself, a validation of a geophysical forcing hypothesis." *Fix:*
present the gate as a consistency check of the tidal arithmetic, not as evidence.

### C. Physical oceanographer

**C1 (correctness, HIGH): day-precision timing is not physical for monthly SST.** "Shifting the
forcing by 1 day collapses your level statistic (P 0.48 to 0.03-0.07 for all three indices, Y1),
and shifts of 1 to 100 days never recover it. Your forcing is the tidal sum at two instants per
year (C8), set by the time-stamp convention of monthly data. No seasonal ocean process (spring
mixed-layer shoaling, westerly-wind-burst season, thermocline-feedback onset) samples the tidal
potential at a calendar instant with less than 1 day of year-to-year jitter. Such processes vary by
weeks. A mechanism this timing-fragile is a numerical alias, unless you can name the gate and show
its date is that stable." What they know: seasonal phase-locking of ENSO has a spread of
several weeks (e.g. onset and peak dates vary by about a month); monthly means average Mf over two
cycles.
*Decisive test:* (i) **jitter**: replace the fixed impulse instant by instant + e_y, e_y ~ U(-J, J)
days independently each year, J = 3, 7, 14, 30. Fully re-fit (restricted estimator), 10 jitter
draws each. (ii) **windowing**: replace the point value of the tide by its average over a window of
W = 3, 7, 15 days around the impulse instant. A physical gate must keep most of its skill at
J = 7-14 d. If skill vanishes at J = 3 d, the claimed mechanism is excluded.

**C2 (correctness, MEDIUM): the LOD weighting is the wrong transfer function for any ocean
pathway.** "LOD weights constituents by the degree-2 order-0 zonal potential (rotation). The
established tide-climate pathway is mixing modulated by the spring-neap and declinational
envelopes of the semidiurnal and diurnal tides (MSf, MNf, ...; Loder & Garrett 1978; Munk &
Wunsch 1998). Its weights are set by M2/S2/K1/O1 energetics and their dissipation sites, not by
LOD. A rotation pathway (delta Omega/Omega ~ 1e-8) is energetically negligible." What they know:
order-of-magnitude budgets. A 0.5 K basin-scale SST anomaly in a 50 m mixed layer needs ~1e21 to
1e22 J.
*Decisive test:* give an order-of-magnitude energy or momentum budget for the chosen pathway. Then
replace the LOD-weighted forcing by the astronomically computed spring-neap envelope of the
local M2+S2+N2+K1+O1 potential, or by TPXO dissipation-weighted envelopes, with the same
estimator. The physically right forcing should fit at least as well.

**C3 (correctness, MEDIUM): one winding for every ocean region.** "The LTE reduction relies on
equatorial dynamics (small f, waveguide). The same k ~ 4 (with the same A) is claimed for the
subpolar North Atlantic (AMO), the North Pacific (PDO) and boxes up to high latitudes, where the
deformation radius and baroclinic wave speeds differ by about 10x." What they know: Chelton et
al. 1998 deformation-radius climatology.
*Test:* on the quads, estimate k_j per quad with the restricted estimator, and regress it on the
local first-baroclinic gravity-wave speed / deformation radius and on basin width. LTE predicts
a dependence; "universal k" predicts none. Either result must be reported with the theory it
supports.

**C-novelty:** "What does the formulation predict that a damped oscillator forced by weather
noise (Hasselmann; Penland & Sardeshmukh 1995) does not? State one observable (a spectral line,
a phase relation to a named tidal argument, a seasonal timing) that the null cannot produce."

### D. Climate-variability community (AMO/PDO)

**D1 (correctness, HIGH for AMO and the global back-cast): the forced response is an unremoved
nuisance.** "Linearly detrended AMO keeps the anthropogenic-aerosol and volcanic signal (Booth et
al. 2012; Mann et al. 2021; Trenberth & Shea 2006). The global SST back-cast target 1880-1949
contains Krakatoa (1883), Santa Maria (1902), Katmai (1912), and the WWII bucket/engine-intake bias
peak (1939-45; Thompson et al. 2008). A decadal latent with time-segregated bins (C6) can align
with any of these." This is the same failure mode as the PNA trend and the brestexcl seasonal
cycle in this project's own history.
*Decisive test:* subtract the CMIP6 historical multi-model-mean forced response (or use the
Trenberth-Shea relative AMO, N Atlantic minus global mean), then rerun S1 versus surrogates and
the back-cast. For the back-cast, include volcanic AOD (Sato) as a regressor. Score against
bias-adjusted HadSST4 and ERSSTv5 as well as Kaplan.

**D-novelty:** "AMO has ~2 cycles in the record. Any positive result on AMO needs a statement of
how many independent cycles support it." (Partly raised by the author for the phase-shift test.
Here it applies to S1: report the effective degrees of freedom of S1 under the decadal
autocorrelation.)

### E. SST-dataset community

**E1 (correctness, MEDIUM): the 89 quads are not 89 independent confirmations, and pre-1950
Kaplan is a projection.** "Kaplan SST is a reduced-space optimal smoother on EOFs from the modern
era (Kaplan et al. 1998). Before ~1950 every grid box is a combination of a few dozen large-scale
modes, fitted to sparse observations. The 89 quads therefore span far fewer independent degrees of
freedom. A back-cast trained on 1950-2023 is extrapolating inside the same mode basis that the
reconstruction imposed."
*Test:* report the effective number of independent quads (e.g. the participation ratio of the
quad correlation matrix, 1950-2023 and 1880-1949). Repeat the back-cast against HadSST4 grid
boxes with coverage masks (no EOF infilling). Give a back-cast null: identical pipeline with F
shifted by whole years (as in the quads N2) and with a fake-Moon forcing (F1 below). The CC 0.52
/ 0.63 means nothing without that null distribution.

**E-novelty:** "Global SST back-casts at CC ~0.5 are routinely achieved by regression on a few
modes. Compare with a damped-persistence / EOF-regression baseline trained on the same 1950-2023
period."

### F. ENSO forecasters (practitioners, incumbent)

**F1 (correctness, HIGH): no lunar-specificity null at equal budget.** "Y2 shows delta = 0 wins,
but the delta = 0 run starts at its own months-long optimum, and no delta != 0 run re-pinned its
tides (dLOD stayed at -0.15 to -0.80, C9). That is an incumbency effect, not a clock test. The
null you need is a counterfactual Moon."
*Decisive test:* **fake-Moon forcing**. Scale all lunar frequencies by factors 0.97, 0.98, 0.99,
1.01, 1.02, 1.03 (giving non-astronomical months, fortnights and nodal periods). Keep the solar
terms. Re-derive the constituent phases by regression on the same dLOD reference (the gate will
fail by design, which is the point). Then fit real AMO/PDO/NINO4 from a COLD, cross-seeded start
with the restricted estimator and an equal budget, and do the same for the real Moon. The real Moon
must win clearly (for example rank 1 of 7, by more than the run-to-run spread from the ctrl2
replicate) on the 1880-1949 holdout of A2.

**F2 (correctness, MEDIUM): forecast benchmark.** "The 2026 El Nino amplitude was missed. On the
post-2023 window and on 1880-1949 NINO4, what is the skill relative to persistence, damped
persistence, a LIM (Penland & Sardeshmukh 1995), and NMME?"
*Test:* skill score versus these baselines on identical windows, with block-bootstrap intervals.

**F-novelty:** "What lead time is claimed? A deterministic tidal driver implies skill at
arbitrary lead. State the lead-time-skill curve and how it differs from the ~6-12 month barrier."

### G. Tidal-climate literature (incumbent on novelty)

**G1 (novelty, MEDIUM):** "Ray (2007, J. Climate, 'Decadal climate variability: is there a tidal
connection?') found tidal-energy modulation too weak to drive decadal SST. Keeling & Whorf (2000)
proposed the long tidal cycles. Yasuda et al. found 18.6-yr signals in the North Pacific. Where
does this work stand against Ray's energetic objection, and what is the advance over
Keeling-Whorf?" *Fix:* a prior-work paragraph that engages Ray's argument quantitatively (ties to
C2's budget).

### H. General reader

**H1 (correctness, MEDIUM):** "Describes 89 quads with one formulation" reads as 89 successes. The
quads test failed its positive control (no power), and Arm B failed. A general reader will not
reach those files. *Fix:* see section 4.

---

## 3. Verdict table

| # | audience | objection | axis | answered where (or NOT) | severity | fix (proposed, not applied) |
|---|---|---|---|---|---|---|
| B1 | geodesist | gate blind to long and minor constituents; they differ per series; F2 85% decadal | correctness | NOT (the claim says "essentially locked") | HIGH | frozen-at-dLOD/IERS forcing refit; F2 variance attribution |
| A1 | statistician | parameter count is not capacity; own file says ~150 params | correctness | contradicted by SURROGATE_RESULT.md | HIGH | restricted-estimator Arm B; GDF |
| A2 | statistician | cross manifold correlates 0.97-0.99 with own; own S2 k~4 circular | correctness | partly (DECLARATION_LEVELS says own is biased); the k~4 claim is still quoted without the caveat | HIGH | 1950-2022 fit, 1880-1949 S1 holdout incl. surrogates |
| C1 | oceanographer | 1-day timing fragility, point-sampled two-slot comb | correctness | NOT | HIGH | impulse jitter and windowing test |
| F1 | ENSO practitioner | no lunar-specificity null; Y2 is an incumbency artifact | correctness | NOT (Y2 read as support) | HIGH | fake-Moon, cold-start, equal budget |
| D1 | climate variability | forced response / volcanic / bucket bias nuisance | correctness | NOT | HIGH (AMO, back-cast); MEDIUM (PDO, NINO4) | CMIP6-MMM removal; volcanic regressor; HadSST4/ERSSTv5 |
| B2 | geodesist | YEAR is a ~8-day lag; Y1 = day shifts | correctness | partly (author noted absolute origin) | MEDIUM | J2000 + explicit lag reparameterization |
| A3 | statistician | S1 partly detects slow curvature (bins time-segregated) | correctness | NOT | MEDIUM | quadratic/forced removal; within-block eta2 |
| C2 | oceanographer | LOD weighting is the wrong transfer for mixing; energy budget | correctness | NOT | MEDIUM | budget; spring-neap envelope forcing comparison |
| C3 | oceanographer | universal k across regimes contradicts the equatorial LTE basis | correctness | NOT | MEDIUM | k versus deformation radius across quads |
| E1 | SST data | quads not independent; pre-1950 Kaplan is a projection; back-cast lacks a null | correctness | NOT | MEDIUM | participation ratio; HadSST4 boxes; shifted-F back-cast null |
| B3 | geodesist | dLOD reference provenance and filter | correctness | NOT | MEDIUM | document; gate against C04 minus AAM |
| A4 | statistician | 2000-05 window still cited as positive | correctness | contradicted by GATE.md | MEDIUM | drop it |
| F2 | ENSO practitioner | no skill relative to LIM / persistence / NMME | correctness | NOT | MEDIUM | skill scores |
| A-nov | statistician | generic decadal latent may do as well | novelty | NOT | MEDIUM | IAAFT-of-F2 latents; CMIP latent |
| G1 | tidal-climate | Ray 2007 energetics; Keeling-Whorf | novelty | NOT | MEDIUM | prior-work paragraph |
| C-nov, F-nov | oceanographer, forecaster | unique observable; lead-time curve | novelty | NOT | MEDIUM | state one discriminating prediction |
| B-nov | geodesist | LOD tidal fit is routine | novelty | NOT | LOW | reframe the gate as an arithmetic check |
| H1 | general | "89 quads described" implies 89 successes | correctness | NOT on the reading path | MEDIUM | scope sentence |

---

## 4. Headline-claim ("abstract") audit, cold

Read as the only text a reader sees, the claim implies four things the files do not support:
1. **"Essentially locked by LOD"**: implies the forcing has no per-series freedom. False for 24
   of 42 constituents and for every decadal term (C1, C2), and these control the decadal latent.
2. **"~7 nonlinear + ~13 linear parameters"**: implies the evidence was obtained at that freedom.
   The surrogate result that exists was obtained at ~150 parameters, and it failed (Arm B).
3. **"k ~ 4 on all real indices"**: holds only on each index's own, co-tuned manifold. On the
   cross manifold the BIC chooses 2.84-3.68.
4. **"89 quads described by one formulation"**: the only quad-level test failed its positive
   control. No quad result is yet evidence.
No conjecture labelling is present at headline level.

Proposed minimal edits (scope sentences, no new claims; NOT applied):
- After "essentially locked": "The dLOD gate constrains the fortnightly and monthly constituents.
  Long-period (>200 d) and minor constituents are not constrained by it and are fitted per series."
- After the parameter count: "This count refers to the restricted estimator, which has not yet
  been tested against surrogates. The full search (~150 parameters) fits phase-scrambled
  surrogates nearly as well as the real indices."
- Replace "picks k ~ 4.0 on all real indices" with "picks k ~ 4.0 on each index's own manifold
  (which was tuned with a k ~ 4 winding), and k = 2.8-3.7 on another index's manifold."
- Headline label: "a candidate description", with the quads stated as "fitted, not yet tested
  against a powered null".

---

## 5. Ranking for the next experiments (by severity and by how cheaply each decides)

1. B1, frozen-forcing refit plus F2 attribution: cheap, and it directly tests the "locked" pillar.
2. C1, impulse jitter/windowing: cheap. A failure at J = 3 d would exclude a physical seasonal gate.
3. F1 + A2, fake-Moon with cold start and 1880-1949 holdout of the whole pipeline: the decisive
   specificity test. It needs the restricted estimator (A1) first.
4. D1, forced-response removal for S1 and the back-cast: cheap.
5. E1, back-cast null and HadSST4 target.
6. C2, C3, B2, B3, F2, then the novelty paragraphs.

**Fatal?** Nothing is shown to be fatal yet. B1 and C1 could each be fatal: B1 to the "few
parameters, LOD-locked" parsimony argument, C1 to any physical interpretation. Both are cheap
to run. If the frozen-at-dLOD forcing keeps the level structure and the fit survives a 7-14 day
jitter, the case becomes much stronger than anything shown so far.
