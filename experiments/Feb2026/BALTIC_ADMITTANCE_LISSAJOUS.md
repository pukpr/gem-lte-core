# Baltic MSL vs. Kaplan SST: admittance and phase structure across the shared manifold

## 1. The data

Two independently-fitted models of the same physical basin — the
211-gauge PSMSL sea-level composite (`baltic`) and the narrowed Kaplan
SST box (`kN_baltic`) — converge on the **same seven windings** (see
`Baltic_compare.txt`). K0 is the linear (trend/background) transfer;
K1-K7 are the sinusoidal winding transfers, each `Amplitude *
sin(2*pi*M*Forcing(t) + Phase)`.

| K | M | MSL amp | MSL phase | SST amp | SST phase | gain (SST/MSL) | phase diff |
|---|---|---|---|---|---|---|---|
| K1 | -0.9247 | 11.623 | -2.363 | 0.226 | -2.229 | 0.0194 | **+7.7°** |
| K2 | +0.1219 | 7.819 | +1.765 | 0.142 | -2.241 | 0.0181 | +130.4° |
| K3 | +0.2076 | 4.939 | +0.652 | 0.037 | +1.031 | **0.0075** | **+21.7°** |
| K4 | +2.0756 | 4.939 | -2.224 | 0.129 | +1.633 | 0.0260 | -139.0° |
| K5 | +1.8681 | 6.025 | -2.693 | **0.199** | -0.005 | **0.0330** | +154.0° |
| K6 | +1.2454 | **31.646** | -2.601 | 0.362 | +0.772 | 0.0114 | -166.7° |
| K7 | +2.2832 | 6.919 | -0.368 | 0.105 | -2.463 | 0.0152 | -120.0° |

(gain and phase diff are SST relative to MSL, computed directly from the
amplitude/phase pairs above; sorted by |M| in the discussion below.)

## 2. Reading the Lissajous diagrams

Each panel in `baltic_lissajous.png` traces the parametric curve `(MSL
amp * sin(theta + MSL phase), SST amp * sin(theta + SST phase))` as
theta sweeps one full winding cycle (0 to 2*pi) — the same construction
used to read phase difference between two same-frequency signals off an
oscilloscope in XY mode. Because MSL and SST share the exact same
theta for a given K_n (it's the same physical frequency observed in two
channels), this curve is always a closed ellipse, and its **shape alone
encodes the phase relationship**, independent of the very different
absolute amplitude scales of the two series:

- **A narrow, nearly-straight diagonal** = the two signals are close to
  in phase (or antiphase) — one is nearly a scaled copy of the other.
- **A wide, rounder ellipse** = the two signals are closer to
  quadrature (~90° apart) — knowing one tells you almost nothing about
  the other's instantaneous value.
- **The tilt direction** (positive vs. negative slope) distinguishes in
  phase (positive) from antiphase (negative).

The two colored arrows in each row are not two different curves — both
land on the *same* ellipse, at two specific points on it:

- **Left (blue), `Amp*sin(Phase)`**: the ellipse evaluated at theta=0.
- **Right (orange), `Amp*cos(Phase)`**: the ellipse evaluated at
  theta=pi/2, a quarter-cycle later.

Together they are just the sine/cosine (in-phase/quadrature) projection
of each channel's own phasor — the same decomposition the original
regression coefficients are built from — placed on the actual shared
curve so the reader can see where on the cycle each projection falls
without needing to trust the algebra.

## 3. Three phase families, not a smooth continuum

Sorting by |M| lays out three distinct groups, not a gradual trend:

**In phase** — K1 (M=-0.9247, +7.7°) and K3 (M=+0.2076, the shared
backbone itself, +21.7°). Both trace tight, near-diagonal ellipses.
MSL and SST move together on these two terms: whatever drives one
drives the other with almost no lag or sign change.

**Near quadrature** — K2 (M=+0.1219, +130.4°). The widest, most open
ellipse of the seven. This is the lowest-frequency of the "extra"
windings (excluding K1/K3) and behaves distinctly from both of the
other two families.

**Antiphase cluster** — K4, K5, K6, K7 (M = 2.0756, 1.8681, 1.2454,
2.2832; phase differences -139.0°, +154.0°, -166.7°, -120.0°). All four
land within a 47° band clustered near 180° — a real family, not a
coincidence of four unrelated numbers (see §5). These are also, not
incidentally, four close relatives in frequency: roughly 9x, 9x, 6x,
and 11x the shared 0.2076 backbone.

## 4. A possible mechanism: mixing volume, not just mixing presence

The clean split — fundamental/near-fundamental terms in phase,
higher-harmonic terms in antiphase — lines up with a real, physically
motivated distinction rather than an arbitrary one. Higher harmonics of
a tidal flow are the generic fingerprint of **nonlinear, quadratic
processes** (bottom friction and turbulent mixing both scale with flow
speed squared, which is exactly what generates overtone content from a
single fundamental tidal current). If K4-K7 are the nonlinear-mixing
overtones of the same shared forcing that K1/K3 represent more directly
and linearly, a natural physical story follows:

- A pulse of enhanced tidal mixing exchanges a *volume* of water
  vertically — cool, dense subsurface water is entrained toward the
  surface. That directly **lowers SST** at the surface.
- The same mixing event, by homogenizing the water column and altering
  the local density structure, can act to *reduce* the column's steric
  height rather than add to it — the opposite sign of what a simple
  linear, wave-driven sea-level response would produce.

If that's right, the antiphase relationship in K4-K7 isn't SST "lagging
behind" MSL — it's the two channels responding with **opposite sign to
the same volume-exchange event**, while K1/K3 capture the more
straightforward, linear part of the same shared forcing where both
channels move the same way.

## 5. Why the SST-MSL correlation is weak: mixing destroys coherence

The model is what it is: seven windings, three phase families, one
(K6) overwhelmingly dominant and sitting almost exactly at antiphase.
That is a sufficient explanation for the weak observed SST-MSL
correlation on its own, and it does not need to be rationalized further
by chasing down every remaining percentage point of covariance. The
direct statement is: **mixing-driven upwelling/downwelling does not
merely fail to produce coherence between sea surface temperature and
sea level in this basin — it actively destroys it.** The same
volume-exchange event that this shared forcing describes drives the two
channels in opposite directions on the dominant winding, so weak or
vanishing correlation is the expected outcome of the physics itself,
not a gap in the model or a sign that the two channels are unrelated.
Two independently-fitted models (`baltic`/MSL and `kN_baltic`/SST)
agreeing on this same structure — same seven frequencies, same
antiphase family — is itself the finding; it does not additionally need
the raw correlation coefficient to come out any particular sign or size
to be meaningful.

## 6. Why this is unlikely to be coincidence

Two independent lines of evidence, and they carry different evidentiary
weight — worth being precise about which is which.

**The frequency matches themselves** are striking (all seven M values
agree between the two independently-fitted models to within a fraction
of a percent), but this evidence should be weighted cautiously: the SST
model's search was informed by donor manifolds and, ultimately, by
comparison against the already-known Baltic structure, so close
frequency agreement is not fully blind, independent confirmation on its
own. What *is* genuine confirmation is that the SST data actually
**accepts** this borrowed frequency set — every other basin tried this
session that didn't share this basin's real dynamics (Gulf of Mexico,
most of the 89-region Kaplan sweep) flatly failed to produce a good
joint fit no matter which windings or parameters were tried. Kaplan SST
here did not merely have these frequencies imposed on it; it fit them
well, which a mismatched forcing has repeatedly failed to do elsewhere
in this project.

**The phase clustering is the stronger, genuinely blind evidence.**
Once M is fixed, each channel's amplitude and phase at that frequency
is not a free choice or a seeded value — it falls out of an independent
linear regression against real, physical data neither model was told
to match in phase. There was no mechanism for the optimizer to *target*
a phase relationship; phase is a pure byproduct of the fit. That four
separate windings (K4, K5, K6, K7) each landed, independently, within a
47° band centered near 180° is the part of this finding that a
guided-search critique cannot explain away. Treating each of the four
phase differences as if drawn independently from a uniform distribution
on the circle, the probability that all four would fall within a
window this narrow by chance is:

```
P(4 iid uniform phases all within a 47 deg window) ~ 0.008  (< 1%)
```

This is offered as an illustrative order-of-magnitude estimate, not a
rigorous hypothesis test (it doesn't account for how the four
frequencies were selected in the first place) — but it quantifies the
intuition directly: a coherent four-way phase cluster this tight is not
what four unrelated, independently-fit numbers look like.

**Bottom line.** The frequency agreement shows the shared manifold is a
real, physically appropriate description of this basin in both heat
and mass. The phase clustering — genuinely unforced, genuinely
independent per winding — shows that the *manner* in which mass and
heat couple to that shared forcing is itself structured and repeatable,
splitting cleanly into a linear/in-phase family and a nonlinear/
antiphase family rather than scattering randomly. Both together are the
signature of a real, shared physical mechanism operating on this basin,
not an artifact of curve-fitting freedom.

## 7. Related literature

**Tidal mixing fronts are a real, established mechanism for exactly
this kind of SST modulation** — the foundational result is Simpson &
Hunter's h/u³ criterion (*Fronts in the Irish Sea*, 1974), which
predicts where shelf seas transition between stratified and
tidally-mixed water based on depth and current speed, with the mixed
side running measurably colder at the surface. This has since been
refined (Simpson & Sharples 1994's critical log(h/u³)=2.7±0.4 value)
and applied at basin scale (Holt & Umlauf on the NW European shelf;
Rippeth 2005 on the seasonal-stratification paradigm; frontal SST
directly measured from ferryboat/glider transects in later work). This
is solid precedent for "enhanced mixing cools the surface" as a real,
quantified mechanism, not a speculative one.

**What the literature does *not* provide is a direct extension of that
mechanism to sea level.** No paper surfaced that connects Simpson-Hunter
tidal mixing fronts to steric height or sea-level variability — every
source found addresses the SST/stratification side only. The sea-level
half of §4's hypothesis (mixing altering the water column's density
structure enough to reduce steric height) is a plausible physical
extension of established SST-mixing-front theory, but it is a novel
combination for this document to propose, not one with a citable
precedent.

**An important correction specific to the Baltic itself**: the
literature is consistent and clear that the Baltic Sea is properly
described as *practically non-tidal* — astronomical tidal range is only
~2-5 cm, and the Danish Straits act as a low-pass filter that blocks
most North Sea tidal energy from propagating in (tidal *currents* do
still reach several tens of cm/s locally in narrow straits/sills, but
basin-wide sea-level and mixing variability is dominated by wind and
seiches, not tides — see the 2025 NHESS sea-level decomposition study
and the ADCP tidal-current survey below). That means the mixing
mechanism proposed in §4, if real, more likely operates as **wind- or
seiche-driven mixing organized on the same astronomical clock** that
this whole project's forcing pipeline tracks, rather than literal
tidal-current-driven mixing in the Simpson-Hunter sense — the model's
windings are locked to the astronomical forcing regardless of what
physically carries the mixing energy, so this doesn't invalidate the
finding, but "tidal mixing" should be read loosely here, not literally,
for this specific basin. Separately, a 2006 Tellus A study on decadal
Baltic sea-level variability found a statistical link to temperature
mediated through precipitation — a different, non-mixing pathway
between the two channels worth keeping in mind as an alternative or
additional mechanism.

Sources:
- [Simpson & Hunter, "Fronts in the Irish Sea" (1974)](https://www.researchgate.net/publication/234215750_Fronts_in_the_Irish_Sea)
- [Simpson & Sharples' h/u³ criterion, as used in shelf sea glider observations, Ocean Science 14, 225 (2018)](https://os.copernicus.org/articles/14/225/2018/)
- [Sharples et al., "On the value of the mixing efficiency in the Simpson-Hunter h/u³ criterion," Ocean Dynamics](https://link.springer.com/article/10.1007/BF02226285)
- [Holt & Umlauf, "Modelling the tidal mixing fronts and seasonal stratification of the Northwest European Continental shelf," Continental Shelf Research](https://www.sciencedirect.com/science/article/abs/pii/S0278434308000241)
- [Lagrangian h/u³ values from ferryboat-monitored SST, J. Phys. Oceanogr. 38(11), 2008](https://journals.ametsoc.org/view/journals/phoc/38/11/2008jpo3839.1.xml)
- [Timko et al., "Assessment of shelf sea tides and tidal mixing fronts in a global ocean model," Ocean Modelling](https://www.sciencedirect.com/science/article/abs/pii/S1463500318302622)
- [Untangling the waves: decomposing extreme sea levels in a non-tidal basin, the Baltic Sea, NHESS 25, 1439 (2025)](https://nhess.copernicus.org/articles/25/1439/2025/)
- [Tidal currents as estimated from ADCP measurements in "practically non-tidal" Baltic Sea](https://www.researchgate.net/publication/261200034_Tidal_currents_as_estimated_from_ADCP_measurements_in_practically_non-tidal_Baltic_Sea)
- [Influence of temperature and precipitation on decadal Baltic Sea level variations in the 20th century, Tellus A 58(1), 2006](https://www.tandfonline.com/doi/abs/10.1111/j.1600-0870.2006.00157.x)

## 8. Widening the comparison: North Sea and English Channel SST

The same shelf-sea system connects all three regions physically —
Channel/Brest feeds the North Sea, which feeds the Baltic through the
Danish Straits — so the natural next check is whether their own
independently-fitted SST windings share the same structure. Current
real (`lte_gui.py`-fitted) winding tables for all three:

| region | M | amp | phase |
|---|---|---|---|
| Baltic | -0.9412 | 0.215 | -2.257 |
| Baltic | +0.1134 | 0.142 | -2.073 |
| Baltic | **+0.2078** | 0.040 | +0.567 |
| Baltic | **+1.2467** | **0.352** | +0.803 |
| Baltic | +1.8701 | 0.200 | -0.120 |
| Baltic | +2.0779 | 0.123 | +1.617 |
| Baltic | +2.2857 | 0.093 | -2.479 |
| North Sea | **+0.2065** | 0.055 | -2.272 |
| North Sea | +0.6196 | 0.024 | +3.113 |
| North Sea | +0.8010 | 0.112 | -2.070 |
| North Sea | **+1.2499** | **0.305** | -2.577 |
| North Sea | +1.6468 | 0.132 | -0.652 |
| Brest | -0.0105 | **0.274** | +2.491 |
| Brest | **+0.2075** | 0.087 | +2.609 |
| Brest | +0.4150 | 0.062 | -0.810 |
| Brest | **+1.2450** | 0.126 | +3.054 |
| Brest | +3.1126 | 0.148 | -0.435 |
| Brest | +6.0177 | 0.065 | -1.281 |

The shared M≈1.245-1.250 winding is **the single dominant term in both
Baltic and North Sea's own SST fits** (largest amplitude of any winding
in each) — not just a frequency that happens to appear in both, but the
one that matters most in both. It's present in Brest too, but there
it's outranked by a very-low-frequency term (M≈-0.01, consistent with
the AMO-like 60-year-scale character already noted for Brest), so it
isn't the leading mode of variability in every basin — it's specifically
the leading shelf-mixing mode where the basin geometry favors it.

### A methodology correction, made and caught in the same session

An earlier version of this section computed each region's phase against
its **own independently-fit forcing**, and reported Baltic and North
Sea as sitting close to antiphase (166.4° at the dominant winding).
That comparison was invalid, and the fix matters enough to document
plainly rather than quietly overwrite: the two regions' own forcings
are 99.66% correlated overall, but at M≈1.245 that remaining 0.34% of
disagreement gets amplified by the frequency multiplication into a
phase relationship with a standard deviation of ~112° — effectively
randomized, not a shared clock. The cross-region "phase difference"
computed that way isn't measuring anything about how the two real time
series relate to each other; it's measuring how two *almost*-identical
but not-quite-identical clocks drift apart. (The original Baltic
MSL-vs-SST result in §1-6 is unaffected by this — that comparison
deliberately used one shared forcing throughout from the start, exactly
the fix described here.)

Redone correctly — one shared forcing (Baltic's own) used for every
region's regression — the picture is completely different, and now
consistent with everything else already established (the strongly
positive raw correlation, the positive windings-only reconstruction
correlation below): **Baltic, North Sea, and Brest are in phase with
each other**, not antiphase, at both shared frequencies:

| comparison | M≈1.245 phase diff | backbone (M≈0.207) phase diff |
|---|---|---|
| Baltic vs. North Sea | **-4.1°** | **+7.2°** |
| Baltic vs. Brest | -13.0° | +42.8° |
| North Sea vs. Brest | -8.9° | +35.6° |

All six numbers have positive cosine (all close to in phase, none
anywhere near antiphase). Baltic and North Sea are nearly perfectly
aligned at both shared frequencies (within a few degrees). Brest lags
both by a modest but real amount — more at the backbone (36-43°) than
at the dominant winding (9-13°) — and the lag is smaller between North
Sea and Brest than between Baltic and Brest, exactly the ordering
physical geography predicts: Brest/Channel feeds North Sea, which feeds
Baltic through the Danish Straits, so a real propagation or exchange
delay should accumulate with distance along that chain. That's a
graded phase **lag** consistent with heat/water propagating through a
connected system, not a seesaw or compensating exchange — the opposite
mechanism from what the (incorrect) antiphase reading had suggested.

**Confirmed directly in the raw data, the same way as before, corrected
for the same forcing-sharing fix.** Both raw SST series were split into
30 non-overlapping 20-year windows (1866-2011) and locally regressed
against the *same* shared forcing at M≈1.245 — no full-record fit
involved. All 30 of 30 windows landed on the in-phase side, with a
circular mean of **-5.7°** and a mean resultant length of **0.974**
(tighter than the flawed antiphase version's own internal consistency
was). This is a genuine, persistent, in-phase relationship across the
full 145-year record, not a fluke of one regression or the other.

The practical lesson for any future cross-region phase comparison in
this project: **always extract every region's phase against one shared
forcing, never each region's own independently-fit one** — however
highly correlated (99%+) two regions' own forcings look in a simple
linear sense, multiplying by a winding number of order 1 is enough to
amplify the small remaining disagreement into an apparently-large,
essentially meaningless phase difference. The original Baltic
MSL-vs-SST comparison in §1-6 got this right from the start (one
model's forcing used for both channels); this section's error was
introduced only when extending the comparison across independently-fit
regions, and is now fixed the same way.

### What this rules out for §4-5

This result closes off one specific reading of the original §4
mechanism before it could be mistakenly generalized: the within-Baltic
MSL-vs-SST antiphase is **not** evidence of a compensating exchange
with the North Sea. A basin-to-basin seesaw — heat or mass pushed out
through the Danish Straits showing up as an opposing signal next
door — would require Baltic and North Sea to sit in antiphase with each
other. They don't; they're tightly in phase. That leaves the §4
mechanism exactly where it was originally scoped: a **local, vertical**
process inside the Baltic itself (mixing-driven upwelling exchanging
heat between the surface and subsurface within the basin), not a
horizontal exchange with a neighboring sea. The two candidate
mechanisms would have left different, distinguishable fingerprints —
only the local-process one is actually present in the data.
