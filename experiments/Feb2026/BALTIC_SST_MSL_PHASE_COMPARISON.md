# Baltic SST vs. MSL: a shared tidal anchor, two different responses

## 1. Summary

Two independent measurements of the same enclosed basin — Kaplan gridded
SST (`kN_baltic`, 54-66N/10-30E, narrowed from the standard 20x20 sweep
grid to the Baltic's actual extent) and the project's own 211-gauge PSMSL
sea-level composite (`baltic`) — both resolve a genuine, sharp winding
ridge at the **same frequency, M = 1.245** (six times the shared 0.2076
backbone; the established Baltic-family harmonic). Both series find it
independently, using their own separately-optimized forcing pipelines,
over overlapping and non-overlapping eras alike (the ridge holds with
perfect continuity even before 1890, decades before the MSL record
begins). That agreement is a strong, cross-validated fact.

Despite sharing that anchor exactly, the two raw time series correlate
only weakly (R = 0.137, 1880.5-2023 overlap; best lagged value ~0.20
detrended, no lag beyond noise). This is not a contradiction. Extracting
each series' own amplitude and phase *at* M = 1.245 (regressed against
the same shared forcing, for a fair comparison) shows why:

| quantity | SST (`kN_baltic`) | MSL (`baltic`) |
|---|---|---|
| phase at M=1.245 | +28.3° | +85.2° |
| phase difference (SST − MSL) | **−56.9°** (SST lags) | — |
| lag in physical time (period ≈2.24 yr) | **≈0.35 yr (~4 months)** | — |
| amplitude, normalized by own std | **0.388** (38.8% of SST's variance) | **0.214** (21.4% of MSL's variance) |
| phase gap at a second mode (M=2.5) | −15.1° (much smaller) | — |

Two findings, not one: a **real phase lag** (SST behind MSL, and it is
mode-specific — the gap shrinks at M=2.5, so it is not a single constant
time-delay applied uniformly across all frequencies), and a **real
difference in relative importance** — this mode is proportionally almost
twice as significant to SST's own variability as it is to MSL's.

The raw-series trend comparison over the same window adds a third,
separate fact: SST shows a modest *positive* trend (+1.02 units/century)
while the MSL composite shows a *negative* one (−3.01 units/century) —
opposite signs. This is very likely a distinct, lower-frequency artifact
(glacial isostatic rebound dominating many of the 211 PSMSL gauges,
particularly in the northern Baltic/Gulf of Bothnia, where land is
rising faster than the ocean itself) rather than part of the same
winding-phase story, and should be kept separate from the M=1.245
discussion below.

## 2. Candidate physical mechanisms for the phase/weight split

None of these are confirmed — they are physically motivated hypotheses
for the next round of investigation, ranked by how directly they follow
from established geophysical fluid dynamics.

### 2.1 Barotropic vs. thermal response timescales (leading candidate)

Sea level in an enclosed, tidally-forced basin is a **barotropic**
response: mass/volume redistribution and basin seiche/Kelvin-wave
dynamics adjust on the order of days, essentially "instantaneously" on
the multi-year timescale of M=1.245's ~2.24-year period. SST is not a
direct mechanical response to the same forcing at all — it is set by the
balance of solar/atmospheric heat flux and how efficiently that heat is
**mixed** through the water column, and heat content has genuine
thermal inertia (a mixed-layer heat budget integrates forcing over
weeks to months rather than responding to it instantaneously). A
forcing that modulates mixing intensity would only show up in SST after
that inertial lag has been absorbed. The observed ~4-month SST-behind-
MSL lag is the right order of magnitude for a mixed-layer thermal
adjustment time, and the direction (SST lagging, not leading) is exactly
what this mechanism predicts.

### 2.2 Indirect coupling path: forcing → mixing → heat content

Related to 2.1, but specifically about *why* the relative amplitude
differs. MSL couples to the tidal forcing directly and (to first order)
linearly. SST couples through an extra intermediate step — the forcing
modulates *mixing/straining*, and mixing intensity then modulates how
much of the ambient seasonal heat flux stays near the surface. That is
a **multiplicative**, not additive, coupling (tidal-mixing strength ×
whatever heat flux happens to be locally available), so its imprint on
SST is not simply a scaled-down copy of its imprint on MSL — the
transfer function is a different, higher-order one, which is
independently consistent with the relative-amplitude finding (SST
leans on this mode *more*, not less, than a simple linear pass-through
would predict, since concentrated mixing events land on top of whatever
seasonal thermal gradient exists at the time).

### 2.3 Consistency with the already-established Coulomb-friction result

Last session's finding — that this whole model family's IIR/friction
term is decisively better described by constant-magnitude (Coulomb/dry)
friction than by any proportional/exponential alternative — is itself a
classic signature of **bottom friction in shallow straits** (the
Danish straits connecting the Baltic to the North Sea are exactly the
kind of narrow, shallow channel where dry-friction-dominated barotropic
exchange flow is the standard textbook picture). That mechanism acts
directly on volume flux/sea level. It is not obviously the same
mechanism that would set SST's response, which reinforces the idea that
MSL and SST are reading two different physical processes off the same
shared astronomical clock, rather than one process observed twice.

### 2.4 Seasonal stratification modulation (secondary candidate)

The Baltic's stratification is strongly seasonal — sharply stratified
in summer, closer to well-mixed in winter. If the tidal-mixing forcing's
efficiency at reaching the surface heat budget itself depends on the
season it lands in, the effective SST response would carry a seasonal
amplitude/phase modulation on top of the base winding, which could
further explain both the phase offset and its non-uniformity across
different M values (2.1's simple thermal-lag story predicts a lag that
grows with frequency; 2.4 predicts a lag that also depends on which
calendar phase each M's cycle happens to test).

## 3. Path forward

1. **Extend the phase/amplitude comparison across the full winding set.**
   Only M=1.245 and M=2.5 have been checked. Building the full
   phase-vs-M and normalized-amplitude-vs-M curves for both channels
   would show whether the lag grows smoothly with frequency (favors 2.1,
   a single thermal time constant) or is mode-specific with no clear
   trend (favors 2.4, or a richer multi-timescale response).
2. **Test the thermal-lag hypothesis against known mixed-layer
   timescales.** The Baltic's published mixed-layer heat-content
   adjustment time is an independent, literature-derived number; check
   whether it's consistent with the observed ~4-month lag before relying
   on this mechanism further.
3. **Season-stratify the coupling.** Split the SST regression by calendar
   season (or fit a seasonally-modulated amplitude term) and check
   whether the mode's amplitude/phase genuinely varies by season, which
   would be direct evidence for 2.4.
4. **Build a two-channel coupled model.** Rather than fitting SST and
   MSL as two independent single-channel LTE pipelines, test whether a
   single shared forcing, split into a fast barotropic branch and a
   slow thermally-lagged branch (one new time constant, not two
   independently-optimized full parameter sets), can reproduce *both*
   channels' amplitude and phase at once. If it can with only one extra
   physical parameter, that would be strong, parsimonious support for
   mechanism 2.1/2.2 over treating the two as unrelated fits.
5. **Generalize beyond Baltic.** Check whether other basins with both a
   Kaplan SST proxy and an independent second physical measurement (sea
   level, salinity, whatever is available) show the same qualitative
   pattern — SST lagging, and carrying proportionally more of the shared
   mode's variance than the more "primary" physical channel. If this
   generalizes, it's a property of the barotropic/thermal split, not a
   Baltic-specific quirk.

## 4. Caveats

- The phase/amplitude extraction used `kN_baltic`'s own forcing for
  *both* regressions (not each series' independently-optimized one), to
  keep the comparison fair. Re-running with `baltic`'s own forcing as
  the shared reference would be a useful cross-check that the phase gap
  isn't an artifact of which of the two forcings was chosen as common.
- Only two modes have been checked; "mode-specific, not a uniform
  delay" is based on a two-point comparison and should be treated as
  suggestive, not established.
- The opposite-sign century trends (SST warming, MSL falling) are very
  likely a separate GIA/land-uplift artifact in the PSMSL composite, not
  part of the same winding-phase mechanism — flagged here so it isn't
  conflated with the M=1.245 findings above in any follow-up work.
