# Phase down, energy up? (`updown_lags.py` → `updown_lags.png`, `updown_lags.json`, 2026-10-27)

**The question.** If the QBO's descent is a phase pattern and the physical driving is upward, what shows the upward
part? Phase cannot: any phase lag is just a filter property. What carries energy is the wave **envelope** (amplitude).

**Data.**
- **Daily:** ERA5 0N 0E 00 UTC at 10/20/30/50/70 hPa.
- **Hourly:** ERA5 u at 30 and 50 hPa, 5S–5N, all 72 longitudes, 32 half-years in 1999–2016 (`docs/era5/data`; 2003 and 2012 still downloading).

**Sign convention.** Lag of the lower level behind the upper level; positive = downward, negative = upward.

| quantity | 10→20 | 20→30 | 30→50 | 50→70 |
|---|---|---|---|---|
| QBO phase (u low-passed 200 d), daily | +119 d (r 0.89) | +76 d (0.94) | +111 d (0.86) | +67 d (0.91) |
| 4–20 d wave phase, daily | +2.2 d (0.21) | +1.5 d (0.36) | +2.2 d (0.22) | +2.4 d (0.44) |
| 4–20 d wave envelope, daily | −0.3 d (0.11) | +0.3 d (0.21) | −0.5 d (0.16) | +1.0 d (0.30) |

Hourly, 30→50 hPa:

| quantity | pooled peak | r | half-years upward |
|---|---|---|---|
| 4–20 d phase, all waves | +61 h | 0.21 | 0/32 |
| 4–20 d envelope, all waves | +6 h | 0.13 | 11/32 |
| **Kelvin phase** (eastward k = 1–3, 4–20 d) | **+75 h** | 0.31 | 3/32 |
| **Kelvin envelope** | **−6 h** | **0.38** | **24/32** (binomial p ≈ 0.004) |

(Rerun 2026-10-27 with 32 half-years. The first run used 24 and gave Kelvin envelope −8 h, r 0.40, 18/24 upward, p ≈ 0.01.)

## Reading

- **The QBO descent is real phase descent, not a lag ambiguity.** Each step is far below half a cycle, and every
  8-year block agrees (blue, panel a). The zero-wind lines go from 20 to 70 hPa in about 250 days.
- **Wave phase also descends, about 100× faster** (2–3 days per level pair, panel b). For vertically propagating
  Kelvin and gravity waves, downward phase is the signature of upward group velocity.
- **The envelope is where "upward" shows, and only once the waves are isolated.**
  - With all 4–20 d variance at one point, the envelope lag is zero within its noise (daily and hourly).
  - With the Kelvin space-time filter, the envelope correlation triples (0.13 → 0.38) and its peak moves to −6 h:
    50 hPa amplitude leads 30 hPa, in 24 of 32 half-years.
  - In the same filtered field the phase goes the other way (+75 h). Panel c shows this opposite tilt, which is the
    visual answer: the energy rises while the phase descends.
- **Caveat on the size.** For a non-dispersive hydrostatic Kelvin wave, the group and phase speeds have equal
  magnitude, so the envelope lag should be about −75 h, not −6 h. The sign agrees with theory; the size is about a
  tenth of the expectation.

  Possible reasons, untested:
  - wave packets that arrive horizontally from the convective sources reach both levels at nearly the same time;
  - damping and Doppler shifting in the QBO shear;
  - the 60-day envelope high-pass.

  The assimilation-window seams cannot produce a −6 h lag: both levels are analysed in the same window, so a seam
  would give zero lag.

## Relevance to the LTE fit

The multi-leg fit reproduces the descent through per-level winding phases on one common forcing (see
`../legs/qbo9_kit_era_ex2015_2021/RESULTS.md`). That is a phase statement and is silent on direction. The upward
part is a separate observable. For the model to say anything about it, it would need a wave-amplitude (envelope)
output to compare with the Kelvin envelope here.

## Next

- Rerun when the hourly 10/70 hPa winds and the daily 10 hPa point series are in.
- Vertical momentum flux (u′w′, or the Eliassen–Palm form) needs hourly T and ω, which have not been requested.

# Accommodating the flux: tests (a) and (c) (2026-10-27)

**Framing.** The LTE solution is the barber pole: its apparent descent is a phase pattern. Two consequences follow.
- The momentum it needs can still be computed from it: the drag G = du/dt + u/τ, and the upward flux at each level,
  F(z) = ∫_z^top ρ G dz′.
- The energy can come either from a tidal body force or from upward wave flux whose *timing* the tide gates.

## (a) Implied drag and flux (`implied_drag.py` → `implied_drag.png/.json`)

**Data.** The six KIT monthly legs of `../legs/qbo9_kit_era70full_ex2015_2021`, observed and model, with τ = ∞ (the
τ = 1 yr run is in the json and similar).

**Results.**
- **The flux needed at 70 hPa is matched:** rms 0.62 (obs) vs 0.56 mPa (model); r = 0.65 on training years and 0.54
  held out (2015–21).
  - Order of magnitude: about 1 mPa, similar to Kelvin plus mixed Rossby-gravity flux estimates.
  - This underestimates the true requirement, because tropical upwelling (w̄ ∂u/∂z) is neglected.
- **Drag by level:** model vs observed r is 0.39–0.60 on training years and 0.37–0.64 held out, best at 30–40 hPa.
  The model's rms drag is 65–90 % of the observed.
- **Drag is less tied to shear in the model** (r(G, ∂u/∂z) 0.49–0.63) than in the observations (0.66–0.81).
  - A pure translating pattern would give r = 1, so both depart from a rigid barber pole.
  - The model's drag includes a column-wide component that ignores the shear zones: panel (b) shows near-vertical
    drag stripes from the impulse comb.
- **The model misses the westerly/easterly asymmetry.**
  - Observed: westerly accelerations are sharp and narrow, easterly ones broad. The 95th-percentile westerly over
    5th-percentile easterly ratio is 1.43.
  - Model: 1.03. A sinusoidal-phase solution is symmetric by construction.
  - The asymmetry is the signature of narrow Kelvin-wave absorption in westerly shear, the part of the momentum
    budget that is physical rather than phase.

## (c) Comb lines in the incoming wave activity (`envelope_lines.py` → `envelope_lines.png/.json`)

**Method.**
- log of the 4–20 d wave activity, monthly, daily ERA5 0N 0E, 1979–2026, at 100 hPa (added 2026-10-27 evening) and
  70–20 hPa.
- Lines q, 2q, 1±q, 2±q, 1, 2, fitted jointly; SNR measured against 400 off-line frequencies.
- **Control:** remove a multiplicative season × QBO-state model (month, local u bin × shear sign, u70 bin).

| hPa | q raw (p) | q after control (p) | 1±q, 2±q raw | state-model R² |
|---|---|---|---|---|
| 100 | 2.1 (0.05) | 1.8 (0.14) | 2+q 2.2 (p 0.04), others ≤ 1.6 | 0.23 |
| 70 | 2.5 (0.01) | **2.2 (0.04)** | none (best 2−q 1.7, p 0.10) | 0.10 |
| 50 | 4.1 (0.00) | 1.0 (0.55) | none (2−q 2.0, p 0.11) | 0.28 |
| 30 | 6.9 (0.00) | 0.5 (0.86) | 1−q 2.8 (p 0.06) | 0.42 |
| 20 | 4.3 (0.00) | 0.6 (0.78) | 1−q 2.4 (p 0.02) | 0.31 |

**Reading.**
- **Wave activity carries the QBO period at every level, but at 50–20 hPa the local QBO state fully explains it.**
  This is ordinary filtering: the mean wind Doppler-shifts and absorbs the waves.
- **The gate's signature sidebands (1±q, 2±q) are absent from the wave activity,** raw or controlled, apart from
  1−q at 20 hPa (p 0.02), which is high in the column and not at the entry level. The QBO wind itself does have them
  (70 hPa 1−q SNR 4.6, p 0.03; 2 at 5.8).
- **At 70 hPa a weak q line survives the control** (p 0.05). The residual peak sits nearer 0.44 cpy than 0.4219,
  close to the 1979–2026 mean QBO frequency. It may be filtering by the wind between 100 and 70 hPa, which the
  control does not see.
- **Verdict (before 100 hPa).** No evidence that the tide gates the incoming wave flux. The test could not separate a
  gate from QBO self-filtering where the QBO is strong, and a gate correlated with the QBO state would be absorbed by
  the control.

### 100 hPa, the entry level (2026-10-27, `get_era5_point.py 100`, 1979.0–2026.8)

The 70–20 hPa numbers moved by at most 0.1 in SNR. They are not identical because 70 hPa's shear control now uses the
100 hPa neighbour, and the off-line draws changed.

- **The QBO reaches 100 hPa weakly.** The wind itself has q at SNR 3.2 (p 0.01), against 12–30 at 70–20 hPa, plus
  strong annual (9.9) and semiannual (6.8) lines.
- **There is no q line in the wave activity.** Raw SNR is 2.1 (p 0.05), and 1.8 (p 0.14) after the control. The
  low-frequency residual peaks at **0.511 cpy**. That is neither q (0.4219) nor the observed QBO (0.456 cpy at 50 hPa
  from zero crossings; the 70 hPa residual peak is at 0.442). So the marginal q value is the skirt of an unrelated
  peak.
- **The gate's sidebands are absent.** Only 2+q is nominal (SNR 2.2, p 0.04, raw and controlled). With six non-seasonal
  lines tested per level, the chance of one at p ≤ 0.04 is about 22 %. Its partner 2−q (1.3, p 0.36) and the 1±q pair
  (0.7, 1.5) are at noise. A multiplicative gate on a seasonal carrier would put power in both sidebands and in q.
- **Verdict.** At the level where the waves enter the column, ERA5 wave activity (4–20 d) shows no tidal comb lines.
  A gate F_in = F0[1 + ε g(M(t))] has to be weak. The q line is the gate's direct signature; the off-line noise
  places that line below about 2× noise in log E. The 70 hPa q residual (p 0.04, peak 0.442 cpy) is most simply QBO
  filtering between 100 and 70 hPa.
- **Caveats.**
  - One grid point (0N 0E), 00 UTC, so the 4–20 d band mixes Kelvin and mixed Rossby–gravity waves.
  - u² is a proxy for wave activity, not the vertical flux.
  - ERA5 at 100 hPa is model-dominated over the ocean.
  - So this rules out a strong gate in this proxy, not a gate in the true flux. Hourly u′w′ (needs T and ω) is the
    stronger test.

# Fast and slow windings: where the stripes and the descent come from in the LTE fit (2026-10-27)

**The LTE picture.** In the response sin(2πk F(t) + φ_L), a fast winding on a sloping plateau (slope s) runs at k·s
cycles per year, so a per-level phase shifts it by only Δφ/(2πk s), a matter of days. Slow windings advance mainly
from step to step between plateaus, so the same Δφ becomes months. One mechanism could then give both the
near-vertical stripes and the gradual descent.

## What the nine-leg fit actually does (`winding_lags.py` → `winding_lags.png/.json`)

**Reconstruction check.** The model is rebuilt exactly (max error ~1e-12 on every leg):

  model = F9[ m(t) − IR·m(t − 1 yr) ], IR = 0.79

where m holds the windings plus per-leg annual and semiannual terms. Those per-leg terms are not saved: `lt.exe.p`
holds only leg 1's, so they are recovered by least squares.

**Manifold.** Plateau slope |dF/dt| ≈ 0.29 /yr; about 2 steps a year.

**The windings are all slow:** |k| = 0.11–2.65, so k·s ≤ 0.8 cpy. Their effective (stepped, aliased) periods run from
12 to 120 months. The k = −0.11 winding has a 28.4-month period, which is the QBO carrier.
- A per-level phase on such a winding does **not** translate it in time, because its changes are pinned to the
  calendar steps.
- Cross-correlating one winding's component between levels gives a lag of 0, or a jump of whole half-years.
- The descent appears only in the *sum* over windings, the yearly IR difference and the per-leg seasonal terms.

**Whole-model lags**, same method as the data (KIT monthly, low-pass 200 d, adjacent levels):

| | adjacent-level lags, 15→20 … 50→70 (months) | cumulative 15→70 |
|---|---|---|
| observed | 1.5 / 2.5 / 2.3 / 1.7 / 2.4 | **10.3 months** |
| model | 0.9 / 1.6 / 1.6 / 1.1 / 1.9 | **7.0 months** |

The model descends, but at about 70 % of the observed rate. This sharpens the earlier "lags match" statement in
`../legs/qbo9_kit_era_ex2015_2021/RESULTS.md`, which measured lags from 15 hPa and rounded to whole months.

**The model's stripes come from the steps, not from plateau windings.**
- 51 % of the model's high-pass (< 120 d) variance falls within ±10 d of a manifold step, which covers 8 % of the
  time.
- Model stripe lags between levels are 0.1–0.4 d (r 0.71–0.95). The observed lags are 1.3–2.3 d (r 0.16–0.36).
- So the model's stripes are vertical and too coherent because a step hits every level at the same instant.

## Adding fast windings directly (`fast_windings.py` → `fast_windings.png/.json`)

**Setup.**
- Basis: sin and cos of 2πk·F(t), for k = 15–300 (k·s spans about 3–90 d periods), with F the nine-leg manifold at
  daily resolution.
- Target: daily ERA5 stripes (u minus a 120-day running mean) at 20/30/50/70 hPa.
- Fit: ridge regression per level, trained outside 2015–21, λ chosen on 2009–15, scored on 2015–21.

| basis | train r | held-out 2015–21 r | fitted stripe lags 20→30 / 30→50 / 50→70 |
|---|---|---|---|
| **manifold F** | 0.22–0.26 | −0.06 … +0.02 | +0.2 / **+2.1 / +2.6 d** |
| RAMP (straight F, no manifold) | 0.22–0.25 | −0.05 … +0.03 | +1.4 / +2.3 / +2.8 d |
| SHIFT (F delayed 1 yr) | 0.18–0.24 | −0.01 … +0.05 | +0.8 / +1.9 / +2.5 d |
| observed | | | +1.3 / +2.0 / +2.3 d |

**Reading.**
- **Kinematically, the picture works.** Fast windings with per-level phases tilt the stripes by days: the fit
  produces 2.1–2.6 d per level pair, as observed. Nothing in the LTE form forces vertical stripes; they are vertical
  in the present fit only because it has no fast windings.
- **But the tilt is not evidence for the manifold.**
  - The controls give the same tilt, because any per-level fit inherits the data's lag.
  - None of the three bases predicts the held-out stripes (|r| ≤ 0.06). The fast stripes are not phase-locked to
    this manifold at the precision available.
  - This matches the SOI fast-harmonic finding: a planted fast term was recovered only when the manifold was
    accurate about 100× beyond dLOD.
- **The two ends behave differently.**
  - The slow end, the descent, is carried by the fitted manifold and holds out of sample: held-out r 0.6–0.8 per
    level.
  - The fast end, the stripes, can be *represented* by plateau windings but not *predicted*. To the model they are
    weather.

## Phase-free pattern tests for fast windings (`plateau_scaling.py`, `plateau_slope_inversion.py`)

**Why held-out phase is the wrong test.** At each step a fast winding's phase jumps by 2πk·ΔF (k ≈ 100, ΔF ≈ 2),
so phase cannot carry across steps. What survives is the pattern: inside plateau p, a fixed set of windings
appears at f = k·|s_p|. Plateau slopes vary 0.13–0.42 /yr (90 usable plateaus) and have the same distribution in
both seasons. So the stripe spectra should line up across plateaus in the coordinate k = f/|s_p|.

**Scaling test (manifold slopes imposed).**
- Per plateau and level: a whitened periodogram over 3–90 cpy, mapped to k = 23–212.
- Consistency C = mean pairwise correlation across plateaus.
- Null: 3000 permutations of slopes, stratified by season × sign of the plateau-mean u, so that neither the seasonal
  wave cycle nor QBO Doppler shifting can create alignment.

| hPa | C (own slopes) | C (permuted) | z | p | C in plain f |
|---|---|---|---|---|---|
| 10 | 0.029 | 0.023 ± 0.005 | +1.09 | 0.15 | 0.031 |
| 20 | 0.033 | 0.031 ± 0.006 | +0.29 | 0.36 | 0.037 |
| 30 | 0.016 | 0.019 ± 0.005 | −0.78 | 0.77 | 0.027 |
| 50 | 0.031 | 0.032 ± 0.005 | −0.33 | 0.62 | 0.044 |
| 70 | 0.006 | 0.016 ± 0.005 | −1.89 | 0.98 | 0.029 |

**Power check** (planting k = 60, 110, 170 on the manifold into the real 30 hPa stripes):

| planted share of stripe variance | z |
|---|---|
| 1 % | −0.6 |
| 2 % | +0.3 |
| 5 % | +4.6 |
| 10 % | +12.4 |
| 25 % | +28 |

With 10 % planted, slope errors between the planted signal and the manifold matter: 5 % error gives z +7.3, 10 % gives
+2.6, and 20 % gives +1.2.

**Reading.**
- No level shows the slope scaling.
- At every level the spectra align better in plain frequency than in k. The stripes are organised by physical wave
  frequency (Kelvin and gravity waves) more than by plateau slope.
- So a fixed-k plateau pattern carries less than about 5 % of the stripe variance, **or** the manifold's plateau slopes
  are off by more than about 15 %.

**Inversion (each plateau picks its own stretch; frequency-domain warping).** This answers the slope-accuracy
loophole.
- Power check: it recovers planted slopes only at 50 % planted (ρ = +0.52); at 10 % it fails (ρ = +0.13). Its
  freedom fits noise: alignment on pattern-destroyed spectra is already 0.37–0.41.
- Real data: the inferred stretches do not track the manifold slopes (ρ −0.25 … +0.15, p 0.10–0.99), nor the
  plateau-mean wind (ρ −0.23 … +0.04).
- Alignment exceeds the shuffled-spectrum null only marginally, at 30 and 70 hPa.
- This is a weak test; it does not close the loophole, only fails to find anything.

**Verdict.** Kinematically, fast plateau windings can produce the observed stripe tilt (`fast_windings.py`). But
the stripe spectra show no sign of the slope scaling that windings on *this* manifold would imprint. Above about 5 %
of the variance, that rules out such a pattern unless the plateau slopes are more than about 15 % wrong.

## Multi-altitude co-inversion of the plateau profile (`joint_plateau_inversion.py`, `joint_plateau_null.sh`)

**Idea.** The manifold's plateau slopes may simply be wrong; fast windings are too sensitive to trust them. But
every altitude shares one manifold, so a real fast-winding pattern would stretch each plateau's stripe spectrum
by the same factor at every level.

**Method.**
- Estimate each plateau's stretch from two disjoint groups of levels, then compare the two estimates.
- Spectra: sine multitaper, 5 tapers.
- Refinement: either local (|s_p|·(1 + ε_p), with |ε| ≤ 0.3) or blind.

**Results.**
- **Power (local refinement, all five levels, planted k = 60/110/170 with the true slopes 20 % off the model's).**
  The recovered corrections correlate with the true corrections at +0.00 for a 10 % plant, +0.19 for 25 % and +0.23
  for 50 %. Per-plateau free stretches are poorly determined by ~170-day single-point spectra.
- **Real data.** The groups' corrections agree: ρ = +0.33 between the far groups ({10,20} vs {50,70}) and +0.50
  between the interleaved groups.
- **Null.** The identical procedure on data circularly shifted against the manifold (24 shifts of 2–40 yr, the same
  shift at every level) gives +0.26 ± 0.10 and +0.42 ± 0.10. 21 % of the shifts reach the real value (z +0.6, +0.9).
  The agreement is therefore method bias plus shared windows, not a shared plateau profile.
- The stratified permutation null, by season × length, does not remove that bias. Only the time-shift surrogate
  does.
- In local mode, a correlation between σ and the model's |s_p| holds by construction and means nothing.

**Reading.** The nonparametric route (90 free stretches estimated from daily single-point spectra) cannot use the
multi-altitude leverage. Its noise per plateau is too high, and its bias is shared across levels.

**More focused route.** Keep the manifold *parametric*, a handful of tidal and impulse parameters shared by every
leg. The multi-level data then constrains ~10 numbers instead of 90. That is an N-leg `lt.exe` run:
- legs: ERA5 daily at 10/20/30/50/70 hPa (and 100 when it arrives);
- terms: fast windings added;
- metric: phase-free, e.g. `METRIC=FB` (band-limited spectral correlation, in cycles per year) over the stripe band, or
  `DTW`/`CTW`.
Validate it first with a planted fast winding on a perturbed manifold, as here. Cleaner input would also help: the
Kelvin-filtered hourly fields (5S–5N, all longitudes) carry far less local noise than a single grid point.

## Is IR an annual feedback or a filter? (`ir_criticality.py`, 2026-10-27)

The model is F9[m(t) − IR·m(t−1 yr)]. A constant, t, t² and the annual harmonics map back into their own span under a 1-yr lag, so IR can act only through F(t−1 yr). The test fixes the manifold F and the k set, then refits all linear coefficients by least squares for IR from 0 to 0.99. It uses five ERA 5-day levels with 2015–21 held out.

- **Not formally degenerate.** The 1-yr manifold advance varies (ΔF = +0.03 ± 2.6). For the lagged copy of each winding with |k| ≥ 0.5, R² on its own sin/cos span is ≤ 0.06. Only the slowest winding (k −0.11) is mostly absorbed (R² 0.72).
- **Not critical.** On the 9-leg manifold, held-out r is flat from IR ≈ 0.4 to 0.99: 10 hPa 0.59–0.61, 20 hPa 0.66–0.67, 30 hPa 0.69–0.71, 50 hPa 0.66–0.70, 70 hPa 0.49–0.52. The only loss is towards IR = 0, where held-out r drops by 0.03–0.16 and DTW worsens. Training r rises steadily to IR → 1, with no interior optimum. The best held-out IR differs from level to level (0.40–0.99).
- **One shared IR is better than freeing the lag.** With IR = 0 plus free lagged windings, training r rises by about 0.1. Held-out r is the same or worse (70 hPa: 0.31 vs 0.52).
- **IR did not cause the CC search's held-out loss.** On the CC-searched 5-level manifold, held-out r is 0.10–0.33 for every IR. On the 9-leg manifold it is 0.51–0.71 for every IR. The loss came from the manifold and k set: the k −0.55 winding was replaced by −2.95.

Reading: IR acts as a regularizing band filter. Its gain 1 − IR·z⁻¹ʸʳ is 1.5–1.9 in the QBO band and 0.06–0.5 at 0 and 1 cpy. Any IR ≳ 0.5 does the job, so its exact value (0.79 vs 0.94) is not identified by the data. These numbers say nothing for or against a physical annual feedback. Caveat: this is a profile with the manifold held fixed. In lt.exe the search moves F together with IR.
