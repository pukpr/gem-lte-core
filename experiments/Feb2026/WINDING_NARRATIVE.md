# Winding-Frequency Signatures, by Index

An analytic read of the `Forcing -> Model` winding-frequency fingerprints
computed by `param_survey.py` (CV window 1880-1885, `EXCLUDE=true`, the
project default, per-index resp override where set — see the script's own
printed header for exact per-index windows). All wave-numbers below are
`|M|`; harmonic order means the integer multiplier of each index's own
`M(NM)` reference (i.e. `M_harmonic / M(NM)`, rounded).

## The manifold comes first

The headline result here isn't any one index's winding fingerprint — it's
that a *shared latent forcing manifold* makes low-DOF fits possible for
seven climate indices that have each, individually, been genuinely hard to
model precisely (ENSO/PDO/AMO in particular have decades of literature
behind exactly that difficulty). This is a SINDy-style result in the sense
Brunton and Kutz mean it: the hard part of a fluid/climate dynamics problem
is rarely the regression once you have the right coordinates — it's
discovering the coordinates (the manifold, the library of candidate terms)
in the first place. Sparse, parsimonious dynamics are a *consequence* of
having the right basis, not something you get by throwing a bigger basis at
noisy data.

The forcing manifold used here isn't discovered from the climate indices
themselves — it's the LTE tidal-forcing manifold, calibrated against the
*measured* Earth delta-LOD (length-of-day) record (`GEM.dLOD`,
`src/gem-dlod.adb`), a well-understood astronomical signal with a known
answer, before that basis is ever reused as forcing input for a climate
index. The terms doing the physical work there are the **monthly (Mm,
anomalistic, ~27.55d) and fortnightly (Mf, ~13.66d) lunar tidal
constituents** (`src/gem-lte.ads`'s `Draconic`/`Tropical`/`Anomalistic`
periods, combined via Doodson arguments `s,h,p,N`) — not the 18.6-year
nodal / 8.85-year perigee periods that dominate the standard tidal
literature. Those long apparent periods are *emergent*: they're the slow
precession rates of the lunar node (`N`) and perigee (`p`) modulating the
fast Mm/Mf terms, and in particular arise from Mm/Mf **beating against the
annual cycle** — the legacy constants in `gem-lte.ads` even construct them
literally that way (`Nodal := Annual*18.6`, `Perigee := Annual*8.85`).
The integrated manifold is the product of that interaction, not a basis
built from 18.6yr/8.85yr as independent drivers. That ordering matters:
the manifold is validated on a tractable problem first, then reused,
rather than being fit ad hoc per index. Once that's in place, each index's
model reduces to a small number of free winding numbers (harmonics of a
common base) plus independent annual/trend terms — which is exactly the
parsimony you'd expect if the manifold is real and shared rather than an
artifact of a flexible enough fitting procedure.

The direct evidence for this is the M(NM) convergence below: seven
*independently* optimized per-index searches — no cross-index constraint
in the fitting procedure at all — land on the same base winding frequency
to within 0.2%. If the manifold weren't real and shared, there would be no
reason for that convergence; each index's optimizer is free to wander
anywhere in its own search space. Everything in the rest of this document —
the per-index differences, the AMO/TNA kinship, the PDO/NINO4 beat pair,
baltic's single dominant mode — is second-order structure layered on top of
that shared backbone, and is only interpretable as physically meaningful
*because* the backbone itself is independently calibrated rather than
curve-fit.

**A note on `baltic`**: its `lt.exe.p` was deliberately calibrated with
`RIDGE`/`COVERAGE`/`ENCLOSING` regularization — a genuine tuning campaign
(visible in its snapshot history), not incidental churn. Its raw amplitudes
aren't in the same physical units as the other six indices (`baltic`'s
`.dat` is mean sea level in mm-scale units; the others are dimensionless-ish
SST/pressure anomaly indices, ~17-63x smaller in raw std dev), so don't
read magnitude comparisons across indices literally. But the specific
*within-index* finding below — one winding mode dominating the rest by
roughly 3x — was checked against baltic's actual calibration windows
(1955-2000, 1979-2000), not just the survey's default probe window, and
holds up in all of them. The one loose end: the live `lt.exe.p` differs in
some regression coefficients from the last labeled snapshot in its tuning
history, so exact figures may drift slightly on a future re-run, though the
dominant-mode structure itself has been stable across every window tested.

## A shared reference frequency across all seven

Before getting into what's different, the most striking thing in the data
is what's the *same*: every single index's `M(NM)` — the base frequency its
own harmonics are built from — lands within `0.2073` to `0.2077`:

| index | M(NM) |
|---|---|
| nino4 | 0.2073 |
| pdo | 0.2073 |
| iode | 0.2073 |
| baltic | 0.2076 |
| amo | 0.2075 |
| tna | 0.2075 |
| nao | 0.2077 |

Seven independently-fit climate indices, spanning the Atlantic, Pacific,
Indian Ocean, and Baltic Sea, converging on the same reference winding
frequency to within 0.2%. This is the same "solid manifold" result this
whole investigation kept running into from other directions — here it shows
up directly in the fitted wave-numbers themselves.

## (1) AMO: a low-frequency winding carrying the multidecadal signal

AMO's dominant mode — by a wide margin — sits at `M = -0.0134`
(`amp = 1.9155`), the lowest-magnitude winding frequency of any base mode
in the whole survey (excluding `k0`). Its second mode, the shared
`M(NM) ≈ 0.2075` reference, carries real but smaller weight (`amp = 1.0528`).
A winding frequency this close to zero corresponds to the longest apparent
period in the reconstructed `Model` — directly consistent with a
multidecadal (~60-year) signal riding on top of the shared, faster
tidal-manifold response. This reads clearly in the data.

## (2) TNA: a vestige of AMO's low-frequency mode

TNA — geographically the closest index to AMO in this set (tropical N.
Atlantic vs. the basin-wide Atlantic) — has a base mode at `M = +0.0167`,
nearly the same location as AMO's dominant `M = -0.0134`. But where that
mode *dominates* AMO's spectrum, in TNA it's roughly co-equal with the
shared `M(NM) ≈ 0.2075` mode (`amp = 0.1328` vs. `0.1359` — a near-tie).
So TNA carries the same low-frequency component AMO is built around, just
demoted from lead role to a supporting one — "vestige" is the right word
for it.

## (3) PDO vs. NINO4: nearly identical dominant windings

PDO's dominant mode: `M = 0.4489` (`amp = 0.8496`). NINO4's dominant mode:
`M = 0.4490` (`amp = 0.4122`). PDO's is indeed marginally lower, confirming
the observation — but by a razor-thin margin (`ΔM ≈ 0.0001`), close enough
that I'd read this less as "PDO is measurably slower" and more as "PDO and
NINO4 are locked to essentially the same basin-scale Pacific winding,"
which is itself a notable result given how differently these two indices
are usually described (PDO decadal/basin-wide vs. NINO4 the faster
equatorial ENSO index).

## (4) NINO4's high harmonic carries real weight; PDO's doesn't

This one needed the exact harmonic orders to read correctly, and the
straightforward version of the claim (NINO4 reaches a numerically higher
harmonic than PDO) turns out to be backwards — PDO's highest harmonic order
is actually **20** (`M = 4.146`) vs. NINO4's **15** (`M = 3.109`). But look
at the amplitude each one carries:

| index | harmonic order | M | amplitude | vs. that index's peak amplitude |
|---|---|---|---|---|
| nino4 | 15 | 3.109 | 0.2405 | ~58% |
| pdo | 20 | 4.146 | 0.0064 | ~0.75% |

PDO's higher-order harmonic is numerically real but energetically
negligible — noise-level compared to its dominant mode. NINO4's 15th
harmonic, by contrast, carries more than half the weight of its own
dominant mode — a genuinely significant, not incidental, feature of the
fit. So the observation holds once "higher harmonic" means *energetically
significant* reach rather than raw frequency: NINO4 has real structure out
near `M≈3.1` that PDO simply doesn't. Tropical instability waves are a
plausible physical candidate — they're a real, well-documented
high-frequency (~20-40 day) phenomenon riding on the equatorial Pacific
that has no PDO analogue, which is consistent with what shows up here.

## Deeper dive: why PDO reads as more decadal than NINO4, despite sharing a dominant winding

(3) and (4) leave a puzzle: PDO and NINO4 share essentially the same
dominant winding (`M≈0.449`), yet PDO is the one that reads as more
decadal in character. AMO's characterization above gives the tool to
resolve it — but the mechanism turns out to be more subtle than AMO's.

AMO's decadal character comes from something literal: one mode sits at
`M≈-0.013`, and a winding frequency that close to zero corresponds
directly to a very long period. PDO doesn't have that — its own
lowest-`M` base mode (the shared `M(NM)≈0.207` reference every index
carries) is nearly silent for PDO, `amp=0.0149`, essentially noise. So
PDO's decadal-ish character can't be coming from a direct low-`M` term
the way AMO's does. Something else is producing it.

What's actually there is a **near-degenerate beat pair**. Both PDO and
NINO4's two strongest components sit at almost the same two frequencies:

| index | dominant | 2nd | ΔM | amplitude ratio (2nd/dominant) |
|---|---|---|---|---|
| PDO | M=0.4489 (amp 0.8496) | M=0.4146 (amp 0.4425) | 0.0343 | **0.52** |
| NINO4 | M=0.4490 (amp 0.4122) | M=0.4146 (amp 0.1612) | 0.0344 | **0.39** |

That second mode in both cases is the order-2 harmonic of the shared
`M(NM)≈0.207` reference — so both indices have access to the *identical*
near-degenerate pair, with essentially the same frequency gap
(`ΔM≈0.034`). This isn't a coincidence particular to one index; it's
built into the shared manifold both are fit against.

Two sinusoids this close in frequency don't just add — they beat,
producing a slow amplitude-modulation envelope at the *difference*
frequency (`ΔM≈0.034`, much lower than either component). That's a
mechanism for synthesizing an effectively low-frequency, long-period
signal out of two mid-frequency components, without needing a literal
low-`M` mode at all — the same end result as AMO's dominant low winding,
reached by interference instead of directly.

The depth of that beat envelope scales with how close the two amplitudes
are to each other. PDO's pair sits at a 0.52 ratio — closer to matched,
so the interference is deeper and the resulting slow modulation is more
pronounced. NINO4's same pair sits at 0.39 — present, but weaker, so the
beat is shallower and less able to dominate NINO4's character the way it
does PDO's.

NINO4 also has somewhere else to put its energy that PDO doesn't: its
order-15 harmonic (`M=3.109`, `amp=0.2405`, ~58% of its own dominant —
the same harmonic from (4)) is genuinely significant, while PDO's
equivalent reach (order 20, negligible `amp=0.0064`) isn't. So on top of
having a weaker beat, NINO4 also carries real fast structure layered on
— pulling its overall character toward the shorter, ~2-7 year ENSO
timescale — while PDO's spectrum stays confined to the low-order band
where its (deeper) beat envelope dominates unopposed.

Same underlying logic as AMO, in other words — low effective frequency
drives decadal-looking behavior — just reached by a different route. AMO
gets there directly, with one mode that's genuinely near zero. PDO gets
there indirectly: two mid-frequency modes, shared with NINO4 almost
exactly, interfering closely enough in both frequency and amplitude to
beat down to an emergent slow envelope — while NINO4, holding the same
ingredients, doesn't lean on them as hard and instead carries real
high-frequency content that keeps it fast.

## (5) IODE resembles NINO4's structure, with extended harmonic reach

IODE and NINO4 share the same base-mode skeleton (both `nm=3`, anchored on
the shared `M(NM)≈0.207` reference, with a second base mode in the
`0.44-0.56` range). Where they diverge is harmonic reach: NINO4's four
harmonics land at orders `{2, 3, 4, 15}`; IODE's four land at
`{3, 8, 15, 21}` — IODE keeps NINO4's order-15 harmonic (`M≈3.11`, IODE's
own is `M=3.109`, essentially identical) but adds an order-8 harmonic NINO4
lacks and pushes further out to order 21 (`M≈4.35`). So "similar to NINO4,
plus additional winding harmonics" is a precise description of what's in
the data — same core, more high-frequency structure layered on top,
plausibly reflecting the Indian Ocean Dipole's own additional dynamics on
top of whatever ENSO-linked component it shares with NINO4.

## (6) NAO: energy spread across several comparable windings

NAO's seven modes carry amplitudes ranging `0.096` to `0.317` — no single
mode dominates the way AMO's or baltic's does. The largest-to-smallest
ratio is about 3.3x, compared to ~6.8x for AMO and ~7.2x for baltic. That's
a genuinely flatter profile — NAO's winding "fingerprint" (visible directly
in `geo_fingerprints.png`) looks like a comb of similar-height teeth rather
than one or two dominant spikes, which is a fair characterization of "a
number of similar winding numbers." (NINO4's profile is comparably flat by
this same ratio measure, worth noting if you want a more rigorously
"flattest of all seven" claim — NAO isn't uniquely flat, but it is
genuinely un-peaked.)

## (7) Baltic: one dominant winding, an order of magnitude above everything else

Baltic's order-6 harmonic (`M = 1.245`) carries `amp = 30.95` — the largest
value anywhere in this seven-index survey by more than an order of
magnitude. Some of that gap is the raw-units difference noted above (baltic
is MSL in mm-scale units, not a normalized anomaly index), so the absolute
number isn't comparable to the other six directly. What is comparable is
the *within-index* dominance: that mode outweighs baltic's own next-largest
(`10.90`) by roughly 2.8-3.1x, a ratio that holds across baltic's actual
calibration windows (1955-2000, 1979-2000), not just the survey's default
probe window — so this isn't a fitting artifact. One dominant mode this
cleanly separated from the rest, distinct from the flatter multi-peak
profiles of NAO or IODE, is a plausible genuine feature of a mean-sea-level
aggregate: averaging dozens of tide-gauge sites should cancel local,
instrumental, and short-period noise while retaining a basin-scale coherent
signal — exactly the setting where a single dominant resonance would read
out cleanly rather than getting buried in site-specific noise the way it
might in any one station's record.

## Taken together

The argument this document makes is really about the manifold, not the
fingerprints — the fingerprints are what the manifold argument predicts you
should see if it's true. The shared `M(NM)≈0.207` reference is the direct
evidence: it's not that these indices merely correlate with similar tidal
forcing, their *fitted* base winding frequency — arrived at by
independent, unconstrained per-index optimization — converges to the same
value, to within 0.2%, across seven indices spanning four ocean basins.
That's the signature of a real shared latent coordinate system, calibrated
externally against measured dLOD rather than fit to any one index, in
exactly the sense Brunton's SINDy framing predicts: once you have the right
coordinates, the dynamics on top of them come out low-DOF and sparse almost
for free.

Read against that backbone, the differences that show up between indices
stop looking like independent curve-fitting choices and start looking like
real, physically-motivated structure: AMO/TNA's shared low-frequency mode,
NINO4/IODE's shared mid-spectrum skeleton, PDO's near-but-not-quite match to
NINO4 via a shared near-degenerate beat pair, NAO's flatter multi-mode
spread, and baltic's single dominant mode standing apart from all of it —
plausibly the signature of a basin-scale resonance surviving in an
averaged, noise-suppressed aggregate where it wouldn't survive as cleanly
in any single site's record. None of that differential structure would be
readable as *physical* rather than *coincidental* without the shared,
externally-calibrated manifold underneath it to judge it against.

## (8) Cross-validation across physical quantities: Baltic SST confirms the MSL winding, back to before MSL data even exists

Section (7) established that Baltic MSL's `M=1.24538` (order-6 harmonic)
dominates its own spectrum. `kN_baltic` (Kaplan reconstructed SST for the
same basin) is fit completely independently — different physical quantity,
different `.dat` file, its own free optimizer run — and lands on
`M=1.24674` at the same harmonic order, agreement to 0.1%, consistent with
the shared-manifold argument above. That much is just another entry in the
convergence table.

What makes this pairing special is that `kN_baltic`'s own record starts in
**1856**, twenty-four years before Baltic MSL's record begins in **1880**.
That gap is a genuine, temporally disjoint natural experiment: any real
signature this winding leaves in the pre-1880 stretch of SST data cannot
possibly have been influenced by MSL data, because no MSL data exists yet
at those dates. This is a strictly stronger test than the cross-index
convergence above — it isn't just "two independently-optimized fits agree",
it's "one fit's own signature shows up, unprompted, in a stretch of a
*different* series that predates the other series' entire existence."

`winding_scalogram.py` (a windowed Gabor-style transform,
`G_X(t0,M) = < window(t-t0) · X(t) · exp(-i·2π·M·Forcing(t)) >`, scored
against a matched-length AR(1) red-noise floor) makes this directly
checkable: does `kN_baltic`'s own fitted `M=1.24674` carry power that
persists, above the noise floor, into the pre-1880 window nothing else in
the two-index comparison touches? Swept across window half-width `sigma`:

| sigma (yr) | earliest window center | bits above AR(1) floor |
|---|---|---|
| 15 | 1863.5 | +3.43 |
| 10 | 1861.0 | +2.70 |
| 7 | 1859.5 | +1.97 |
| 5 | 1858.5 | +1.23 |

At every window size tested, power at this M stays above the red-noise
floor right up to within a couple of years of the actual start of the SST
record itself (1856) — weakening gracefully as the window shrinks (less
data, more noise, exactly as expected) but never crossing zero. Restricted
to the pre-1880 stretch specifically (`sigma=15yr`, n=5 windows spanning
1871.5-1879.5), mean power there is **+3.47 bits above floor — actually
stronger than the post-1880 average (+2.50 bits)**, and essentially flat
across those years (`log2(power)` = -3.35, -3.35, -3.36, -3.37, -3.39). A
control frequency one unit away (`M=1.7467`, not a fitted winding of
either index) shows the contrast: only +1.55 bits pre-1880, and *negative*
(-0.43, below the floor) post-1880 — this kind of cross-era persistence is
specific to the real winding, not a generic property of any nearby M.

## A note on method: why this needed a stationarity test, not a residual-minimizing one

This result came out of a deliberate methodological detour worth recording.
Two general-purpose discovery tools were tried first, both scored by
aggregate residual reduction against a dictionary of `sin(M·manifold)`
candidates: a PySR symbolic-regression search (`pysr_optimizer.py`), and a
grid + LASSO sparse-selection search (`grid_lasso_optimizer.py`). Both were
validated on a *clean, noise-free synthetic* target built from Baltic's own
known windings and — under enough tuning — could recover `M=1.245` almost
exactly there. Neither could recover it reliably on the real, noisy Baltic
MSL data itself; the LASSO search instead settled on other, less coherent
combinations that reduced the fold's aggregate residual just as well.

The reason isn't a tuning failure, it's a mismatch between the scoring
criterion and the hypothesis under test. LASSO (and STLSQ, the sparse
solver behind SINDy-style discovery) asks one question of a candidate term:
*does including it reduce total squared error over the fitted stretch?* A
term that fits beautifully in one decade and contributes nothing (or
actively hurts) in another can satisfy that criterion exactly as well as a
term that is genuinely, persistently present throughout — aggregate
residual reduction has no way to distinguish the two. That is a reasonable
prior for the kind of problem SINDy is built for: recovering the
right-hand side of an **autonomous or forced-autonomous ODE**,
`dy/dt = F(y(t), u(t))`, where `y`'s own current state is central and the
underlying assumption is a single, time-invariant, pointwise functional
relationship — the same `F` should hold at every sampled `(y, dy/dt)` pair,
regardless of which era it came from. Even SINDy variants that explicitly
acknowledge an external forcing/manifold `u(t)` keep that same
autonomous-ODE scaffolding: `u(t)` is folded in as one more library term
feeding a self-referential state equation, still scored by aggregate fit.

The model underneath this whole document is not that kind of object at
all. `y(t) = trend(t) + Σ Amp_i·sin(M_i·Forcing(t) + phase_i)` is a pure,
memoryless **transfer function** of an externally-calibrated, non-autonomous
manifold — `y`'s own past values never appear anywhere on the right-hand
side. There is no ODE state to discover, autonomous or otherwise; all of
the apparent complexity comes from `Forcing(t)`'s own rich, independently-
determined structure (real astronomical periods, calibrated against dLOD
before ever touching a climate index — see "The manifold comes first"
above), not from any self-referential memory in `y`. For a hypothesis of
that shape, the right validation question isn't "does this term reduce
aggregate error" — it's "is this term's relationship to the external clock
*stationary*", i.e. does it hold up independently in disjoint eras,
including ones the fit never saw. `winding_scalogram.py`'s windowed
transform asks exactly that question, per candidate `M`, as a function of
time — which is precisely why it could find, and the pre-1880 SST result
above could confirm, a signature that both aggregate-residual-based search
methods missed on the same real data.

The practical takeaway: a researcher's toolkit built around autonomous or
SINDy-style dynamical discovery — however readily it accommodates an
external forcing term — carries an implicit residual-minimization prior
that is a poor match for a genuinely non-autonomous transfer-function
claim like this one. That mismatch, not the underlying physics, is
plausibly why a result like the one above can read as unfamiliar to
researchers coming from that tradition: the discovery machinery they're
used to was built to reward a different kind of structure than the one
being claimed here.
