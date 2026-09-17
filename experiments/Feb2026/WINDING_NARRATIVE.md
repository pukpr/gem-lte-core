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
