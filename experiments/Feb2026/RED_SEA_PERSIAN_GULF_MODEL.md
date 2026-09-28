# The Red Sea / Persian Gulf Basin (kN020_E050): a GEM-LTE Case Study

Region: Kaplan SST quadrangle centered at 20°N, 50°E (10°N–30°N, 40°E–60°E),
covering the Red Sea, the Gulf of Aden approach, the Persian Gulf, and the
Gulf of Oman. Two narrow, sill-bounded straits gate exchange with the open
ocean at either end of this box: Bab-el-Mandeb (Red Sea ↔ Gulf of Aden) and
the Strait of Hormuz (Persian Gulf ↔ Gulf of Oman). This geometry — two
young, shallow, evaporation-dominated seas each choked to the open ocean
through a single narrow gap — is what makes the region interesting both to
GEM-LTE's forcing-domain winding analysis and, independently, to
first-principles basin hydrodynamics.

This note documents the region in two passes — first the production
1950–2023 fit (`kN020_E050/`), then the full 1856–2023 span
(`kN020_E050_/`, the trailing-underscore "full record" companion directory)
— and closes by relating both to Johannes Lawen's basin-dynamics modeling
work, which happens to include a validated hydrodynamic model of Doha Bay,
inside this exact quadrangle.

## Part I — The 1950+ model

### Setup

`kN020_E050/` fits monthly Kaplan SST (877 points, 1950.0–2023.0) against
GEM-LTE's tidal forcing manifold with `NM=2` (two independently-searched
winding rates plus their harmonics), `DECAY=48`, `DTW`-then-`CC` two-phase
optimization, and `EXCLUDE=true` (the optimizer never sees 2000–2005;
that span is a held-out block used only for post-hoc validation). The
current production fit reaches **r ≈ 0.74** data-vs-model correlation over
the full 877-month span.

### Winding structure

Every GEM-LTE index in this project shares one universal, spatially
uniform "backbone" term — and the Red Sea/Persian Gulf box is no
exception: its fitted `ltep` backbone is **M = 0.20750**, matching the same
constant found independently in AMO, PDO, NAO, PNA, and Baltic (all
0.207–0.208). What distinguishes this region is the *density* of
harmonics riding on top of that backbone — the `lpap` table carries
`harm = [2, 5, 6, 7, 8, 11, 15, 17, 31]`, and several of the resulting
winding rates (≈1.66, 1.87, 2.28, 2.90) show large, well-fit amplitudes
with no analog in most other quads. Earlier basin-classification work in
this project flagged these as plausible candidates for the region's *own*
basin eigenmode (the same phenomenon that gives the Gulf of
Mexico/Aegean/Black Sea cluster its shared order-6 term), pending a
second, comparably shaped silled basin as a cross-check — which, as
discussed below, is exactly the geometry the Persian Gulf provides
relative to the Red Sea.

Also structurally important is the once-a-year *impulsive* gate,
`Impulse_Delta`, keyed on `delB = 1.126`: rather than a smooth annual
sinusoid, this fires as a sharp, single-month-per-year kick, and its
amplitude and month-of-year slot are both free-fit parameters. That this
region's SST is "uniformly consistent with continuous winding ridges" (the
observation that motivated this whole investigation) coexists with a
*non-smooth* seasonal gate turns out to matter for two of the tests below.

### Temporal power spectrum: what the model actually reproduces

A standard calendar-time periodogram (`power_spectrum.py`, AR(1) red-noise
floor, 200 surrogates) tells a different — and complementary — story from
the forcing-domain winding scalogram. The winding scalogram shows dense,
persistent ridge structure across essentially the whole record; the plain
frequency-domain view shows the model's fidelity to the real data is
*sharply band-dependent*:

| band (cyc/yr) | data's share of its own variance | model/data power ratio |
|---|---|---|
| 0.05–0.30 (multi-year) | 0.225 | 0.52 |
| 0.30–0.70 | 0.146 | 0.37 |
| 0.70–1.30 (annual) | 0.271 | 0.53 |
| 1.30–2.50 | 0.158 | 0.36 |
| 2.50–4.00 | 0.072 | 0.08 |
| 4.00–6.00 (sub-5-month) | 0.041 | 0.06 |

Only 4 of the data's peaks clear the AR(1) significance floor at all
(annual, semi-annual, and two faster harmonics), and together they account
for under 10% of total periodogram power — the real signal here is
broadband, not a clean annual cycle. The model reproduces roughly half the
real power at annual/interannual timescales but only ~6–8% of it at the
fastest bands. **Lesson reinforced from this region specifically:** a rich,
continuous forcing-domain winding ridge does not imply uniform fidelity
across calendar-time frequency — the two are genuinely different axes.

### The smooth seasonal terms are a minor player

Isolating just the model's explicit seasonal regression terms (`ann1`,
`ann2`, `sem1`, `sem2` — a fixed 1×/2× per-year sin/cos basis) and
comparing their own peak-to-peak excursion against the real data's
1.0-to-(-1.0)-scale excursion confirms directly what the user flagged from
inspecting `lt.exe.p`: **these four terms together account for only
roughly 10% of the real seasonal-cycle amplitude.** Whatever gives this
basin its strong, regular seasonal *shape* is not primarily these smooth
harmonics — which is exactly why the impulsive `Impulse_Delta` gate (and
the `Annual_Impulse` test below) matter more here than in most other
regions.

### The Annual_Impulse (ImpC) probe

Motivated directly by the finding above, a new mechanism was added to
`gem-lte-primitives-solution.adb`: `Annual_Impulse`, a genuine once-a-year
*delta* (not a sinusoid) with its own independent amplitude carried by the
previously-inert `ImpC` field. It deliberately reuses `Impulse_Delta`'s own
`DelB`-derived month-of-year slot rather than adding a second free phase,
so the test is specifically "does a sharp kick *coincident with* the
existing tidal-gating impulse explain variance the smooth terms can't
reach" — not an independently-phased new cycle. `ImpC = 0.0` is a no-op by
construction, so every pre-existing fit is untouched.

Implementing this caught a real latent bug worth recording: `DelB`'s own
raw `DPos` slot is not naturally range-limited to a valid month index, and
this region's own fitted `DelB` value lands outside that range —
`Impulse_Delta`'s second branch already guards against exactly this with a
`mod` wraparound, but the naive first draft of `Annual_Impulse` omitted it,
which would have made the feature a silent, permanent no-op for this
region specifically. Fixed and rebuilt clean.

Empirically, testing this against kN020_E050: at `DelB`'s own coupled
phase, no signal (p ≈ 0.72). Scanning all 12 calendar months
independently found one marginal, inconclusive lead at August (p ≈ 0.005,
ΔR ≈ +0.0028 — borderline against Bonferroni correction for testing 12
months) and a weaker echo in December; every other month, null. **Verdict:
mostly negative, one unconfirmed lead.** Interestingly, once this feature
was live and the region was re-optimized through further production
passes, the fitted `ImpC` did settle on small nonzero values (currently
−0.019 in the 1950+ fit, −0.059 in the full-record fit) — consistent with
a real but modest effect, not with the null finding being simply
overturned. This is exactly the kind of ambiguous, small-effect-size
result that the "decompose nuisance before reporting a match" discipline
says to report honestly rather than round up.

### The 2000–2005 holdout: confirmed non-leaking

`EXCLUDE=true` blocks 2000–2005 from training. A rolling-window
correlation over the full record dips specifically and worst at the
*center* of that untrained gap (furthest from any training anchor:
r ≈ 0.40 at year 2003.25) and recovers immediately on both flanks
(r ≈ 0.75 before, r ≈ 0.84 after) — the textbook signature of a genuine,
non-leaking blocked holdout, and a match to the GUI's own
`kN020_E050site2000-2005.png` diagnostic.

## Part II — The full 1856–2023 span

### Why pre-1950 data exists here at all

Kaplan SST's ship-based (pre-satellite, pre-1950) coverage is patchy
worldwide, but the Red Sea/Persian Gulf corridor is an exception: it has
carried dense shipping traffic — and therefore dense ship-log SST
sampling — since long before 1950, plausibly reinforced by Suez Canal
traffic since 1869. That gives this region one of the longest trustworthy
pre-satellite records in the whole 89-quad grid, and made it a natural
candidate for a stationarity check that most other quads can't support.

### Companion directory and current numbers

A full-span companion directory, `kN020_E050_/` (trailing underscore,
matching this project's naming convention for "same region, full
historical record" companions, with its own `kN020_E050__loc.png` map),
was seeded from the production fit and has since been iterated on directly
through `lte_gui.py` (a series of dated checkpoints — 1950-focused passes,
full-1856 passes, and a couple of runs through 1970 — are preserved
alongside the current `lt.exe.p`). As of the current fit:

| span | r (data vs model) | n |
|---|---|---|
| Full record, 1856–2023 | 0.628 | 2005 |
| Pre-1950 only (1856–1949.9) | **0.666** | 1128 |
| Post-1950 only (1950–2023) | 0.406 | 877 |
| 2000–2005 window | 0.402 | 60 |

The pre-1950 era fits *at least as well as*, and by this snapshot somewhat
better than, the post-1950 era it was never specifically tuned against
independently — a non-trivial stationarity result: whatever the winding
model is capturing here does not appear to be an artifact of the
better-instrumented modern era. (Because this directory is being actively
refit through the GUI, these exact numbers will drift fit-to-fit; the
qualitative finding — pre-1950 fidelity is comparable to or better than
post-1950 — has held up across more than one iteration now.)

### Scalogram persistence

`winding_scalogram.py --m-max 5` run against the full 1870–2010 span shows
several bands running essentially unbroken the entire length — most
notably the universal M=0.500 comb, with M=0.207 (the backbone itself),
1.660, 2.282, and 2.905 all showing persistent structure well into the
pre-1950 era rather than only appearing once modern data density kicks in.
That visual claim was independently re-verified by direct inspection of
the saved scalogram image before being reported.

## Part III — Connecting to Johannes Lawen's basin-dynamics work

Two things pointed at Johannes Lawen's research for this write-up:
`environment.report` (his modeling portal — a single-page app whose live
content couldn't be scraped directly, so the summary below leans on his
published work and institutional profile instead) and `wavedyne.com`,
his own resolved coastal-ocean hydrodynamic model (the specific `cv_jl.pdf`
link 404'd, but his CV content is mirrored via his TU Hamburg-Harburg GHI
profile and publication record).

**What Wavedyne is.** Lawen's own model is a Voronoi-mesh-based, fully
parallelized finite-volume coastal ocean code — shallow-water Navier–Stokes
with Coriolis, viscous stress, wave-resolving dynamics, and sediment
transport via the Rouse number — built as an alternative to
triangle-mesh models like FVCOM, and fast enough to exceed real-time
performance on commodity hardware. Critically, **the model's most recent
published validation (2025, *Ocean Science*) is against tidal
observations in Doha Bay, Qatar** — inside the Persian Gulf half of this
exact kN020_E050 quadrangle — with tidal-station RMSEs under 7% across two
seasons.

**Earlier tidal-dissipation work.** Before Wavedyne, Lawen was a
co-author on Yu et al. (2017, *Continental Shelf Research*), "Tidal
propagation and dissipation in the Taiwan Strait," an FVCOM study of how
the M2, S2, K1, and O1 tidal constituents propagate through and dissipate
within a narrow, strait-bounded shelf sea. That is precisely the class of
problem the Red Sea and Persian Gulf pose: both are shallow seas gated to
the open ocean through a single narrow sill (Bab-el-Mandeb, Hormuz), the
exact geometry this project's own grid-wide harmonic classification
flagged as producing *discrete basin eigenmodes* — quantized by
connectivity geometry, not a smooth function of distance — rather than the
smoothly-decaying or globally-uniform behavior seen in open-ocean quads.
Lawen's strait-dissipation framework gives a first-principles physical
mechanism for exactly that empirical distinction: a narrow sill acts as an
impedance-mismatched gate that partially reflects incoming tidal energy,
setting up a basin-scale standing structure with its own resonant
frequencies, distinct from and generally higher-Q than open-shelf tidal
response.

**Heat exchange and evaporation.** Lawen's current doctoral work (from
2026, at TUHH's Institute of Geo-Hydroinformatics) is explicitly focused
on surface heat exchange and evaporation. This is directly relevant here:
the Red Sea and Persian Gulf are among the most strongly evaporative,
negative-freshwater-budget semi-enclosed seas on Earth, each running a
persistent two-layer exchange flow through its bounding strait (fresher,
lighter inflow at the surface; saltier, denser outflow at depth) driven
by basin-wide evaporative heat and salt loss. That overturning has its own
natural annual rhythm, set by the seasonal cycle of evaporation and surface
heating/cooling — a mechanistically *impulsive*, once-a-year-forced
process, not a smooth sinusoidal insolation cycle. That is a plausible
physical candidate for exactly the gap this document's Part I identified:
the finding that this region's real seasonal SST *shape* is only ~10%
explained by smooth annual/semi-annual harmonics, with the rest evidently
riding on the sharper `Impulse_Delta`/`Annual_Impulse` gate. The marginal
August lead found in the `ImpC` probe above is far too weak to call this
confirmed, but it is at least consistent in kind with an evaporation-driven
exchange-flow pulse rather than a smooth radiative seasonal cycle.

**How the two approaches relate.** It's worth being explicit that
GEM-LTE's kN020_E050 fit and a Wavedyne-style simulation are different in
kind, not degree. The GEM-LTE model here is a *reduced-order, closed-form
regression* against a shared tidal-forcing manifold — a small number of
winding rates and impulsive gates fit by direct linear/nonlinear
regression, with no explicit representation of basin geometry, bathymetry,
or momentum balance. Wavedyne is a *resolved, first-principles finite
volume solver* of the actual momentum and continuity equations on a real
mesh of the basin. This project has run exactly this kind of pairing
before — resolved shallow-water PDE side-studies alongside the empirical
winding approach for AMO, QBO, PDO, NAO, PNA, and Brest — and the general
pattern has been that the two lenses are complementary rather than
redundant: the resolved PDE can test *whether a proposed physical
mechanism is even dynamically plausible* (e.g., whether Coriolis alone
can produce a given resonance), while the winding regression is far
cheaper to fit and can consume 150+ years of sparse, gappy real data
directly. A genuine next step — not undertaken here — would be a
Wavedyne-style resolved run of the Red Sea/Persian Gulf system with a
realistic seasonal evaporation forcing, specifically to check whether its
exchange-flow response through Bab-el-Mandeb/Hormuz is closer to an
impulsive once-a-year gate or a smooth sinusoid, which would give an
independent, first-principles read on the `Annual_Impulse` question this
document leaves only marginally resolved.

## Summary

- The 1950+ fit reaches r≈0.74, sits on the same universal M=0.2075
  backbone as every other index in this project, and carries an unusually
  dense set of high-order harmonics plausibly tied to a basin-specific
  eigenmode.
- Calendar-time power-spectrum fidelity is strongly band-dependent
  (~50% at annual/interannual, ~6–8% at sub-5-month) despite dense
  forcing-domain winding structure throughout — a reminder that the two
  analyses answer different questions.
- Smooth seasonal harmonics explain only ~10% of the real seasonal
  amplitude; the `Annual_Impulse`/`ImpC` mechanism built to test an
  impulsive alternative found a mostly-negative, one-lead-inconclusive
  result.
- The 2000–2005 holdout is a confirmed genuine (non-leaking) blocked
  validation window.
- The full 1856–2023 record shows pre-1950 fidelity comparable to or
  better than post-1950 — real stationarity, not a modern-data artifact
  — and persistent scalogram ridges reaching back to 1870.
- Johannes Lawen's Wavedyne model has already been validated inside this
  exact quadrangle (Doha Bay), and his strait-dissipation and
  evaporation/heat-exchange research both offer plausible first-principles
  mechanisms for two of this region's more distinctive empirical findings
  (basin eigenmode harmonics; the impulsive-vs-smooth seasonal shape
  question) — worth a resolved-PDE follow-up study, not yet performed.

## Sources

- [Johannes Lawen — GHI, TU Hamburg-Harburg](https://www.tuhh.de/ghi/people/mr-johannes-lawen)
- [Wave-resolving Voronoi model of the Rouse number for sediment entrainment (Ocean Science, 2025)](https://os.copernicus.org/articles/21/877/2025/)
- [EGUsphere preprint of the above](https://egusphere.copernicus.org/preprints/2024/egusphere-2024-1213/)
- [Yu, Yu, Wang, Kuang, Wang, Ding, Ito & Lawen (2017), "Tidal propagation and dissipation in the Taiwan Strait," Continental Shelf Research 136:57-73](https://ui.adsabs.harvard.edu/abs/2017CSR...136...57Y/abstract)
- [environment.report](https://www.environment.report/) (Lawen's modeling portal; single-page app, content not directly scrapable — cited for reference, not directly quoted)
- `wavedyne.com/cv_jl.pdf` — link provided by the user 404s as of 2026-09-21; not directly accessible
