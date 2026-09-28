# SAM: a narrative, a history, and its connective tissue to QBO

## 1. What SAM is, and where it came from

The Southern Annular Mode (SAM, also called the Antarctic Oscillation,
AAO) is the leading mode of atmospheric variability in the extratropical
Southern Hemisphere: a north-south seesaw of atmospheric mass between the
mid-latitudes (~40°S) and the Antarctic margin (~65°S). In its positive
phase the SH polar vortex and jet stream contract poleward and
strengthen; in its negative phase they relax equatorward and weaken. It
is not a regional pattern tied to any ocean basin or continent — it is
a **zonally symmetric, wavenumber-0 annular ring** that encircles the
entire hemisphere at once.

The index used in this project (`sam/sam.dat`) is the Marshall (2003)
station-pressure index, built from 12 sub-Antarctic and mid-latitude
surface pressure stations and extending back to 1957 — the same
"space-age reliable" starting point the project has already leaned on
for polar-motion (EOP) data. SAM's modern significance in the climate
literature is hard to overstate: it is now understood as a primary
driver of Antarctic sea-ice variability, Southern Ocean carbon uptake,
and drought/rainfall patterns across Chile, southern Australia, and
South Africa. Its most famous single finding (Thompson et al. 2011 and
predecessors) is that stratospheric **ozone depletion over Antarctica
is a major forcing of its multi-decadal trend** — the ozone hole has
pushed SAM toward its positive phase since the 1970s, strengthening and
poleward-shifting the SH jet. That trend is directly visible in this
project's own per-winding amplitude analysis (`sam_ws_params.png`,
several turns back): most of SAM's windings show the same qualitative
shape — high 1957-65, a dip through the 1960s-80s, and a renewed rise
from ~1990 onward — lining up with the ozone hole's own timeline
almost exactly.

## 2. The finding: near-perfect ridge striations

Re-checked against the properly conservative AR1-floor significance
test (not the anti-conservative white-noise floor that produced a false
positive on pure red noise a few turns back — see §4), SAM's real,
jointly-optimized forcing shows **three ridges with perfect continuity
(1.00)**:

| M | peak (bits) | FWHM | continuity |
|---|---|---|---|
| 4.25 | **3.90** | 0.03 | **1.00** |
| 2.33 | 2.61 | 0.09 | **1.00** |
| 2.16 | 2.36 | 0.05 | **1.00** |

This is the cleanest, most rigorously confirmed result of any index or
region examined in this entire investigation — stronger than the
Baltic's M=1.245 (continuity 0.77), stronger than the North Sea's
independently-verified M=1.245 (once the donor-forcing bug was fixed).
No other case this session has produced *three* ridges at perfect
continuity simultaneously.

## 3. The connective tissue: why QBO is the right comparison

This project's own index taxonomy (documented in the sibling
`winding_scalogram` repo) already separates two categories:

| category | forcing anchor | example |
|---|---|---|
| "ocean" (wavenumber > 0) | dLOD-calibrated tidal comb, region-specific | AMO, PDO, NINO4, Baltic, kN040_E030, ... |
| **"k=0" (wavenumber = 0)** | draconic declination torque, zonally symmetric | **QBO**, Chandler wobble |

SAM belongs structurally in the second category, not the first — and
this is not a superficial resemblance. QBO is the equatorial
stratosphere's own zonally symmetric (wavenumber-0) oscillation, and
this project has already established, independently and rigorously,
that QBO's winding is not a free-running oscillator but a **phase-locked
response to the draconic (lunar nodal) torque**, with its winding number
set dynamically by the reciprocal of the forcing's own local plateau
slope (`QBO_FORMULATION.md`). SAM's own near-perfect ridge striations —
continuity=1.00, not just "passing" — are exactly the signature that
finding predicts: a genuinely forced, phase-locked zonally symmetric
mode, not a regionally-confined standing wave competing with local
noise the way an ocean basin does.

There is also a real, independently-published physical mechanism
connecting the two, not just a shared category label: QBO phase is
documented to modulate polar-vortex strength (the Holton-Tan effect,
established for the NH and increasingly documented for the SH), via
wave propagation and refraction in the stratosphere. If the same
draconic-tidal locking mechanism that drives QBO in the equatorial
stratosphere also reaches the extratropical/polar troposphere — whether
directly, or relayed through QBO's own well-documented modulation of the
polar vortex — that would be a coherent, physically motivated pathway
for why SAM shows the *same kind* of clean, forced ridge structure as
QBO, rather than the messier standing-wave-competing-with-noise picture
seen everywhere in the "ocean" category.

## 4. Honest caveats

- **The white-noise floor lesson applies here too, and was already
  worth re-checking before writing this.** A few turns back, pure red
  noise (zero real signal) was shown to produce a false "bright ridge"
  against `winding_scalogram.py`'s white-noise floor. SAM's result has
  now been re-verified against the proper AR1-floor test specifically
  *because* of that lesson — the perfect-continuity result above is the
  AR1-floor number, not a reading off the visual scalogram, and it
  survives the more skeptical test cleanly.
- **The k=0/QBO connection is a structural and mechanistic hypothesis,
  not a proven causal chain.** No cross-correlation or phase-comparison
  between QBO's own winding and SAM's has been run yet — the case above
  is built on shared category membership, shared "phase-locked not
  free-running" behavior, and an independently-published physical
  coupling pathway (Holton-Tan), not a direct empirical test between
  the two series.
- **The Marshall SAM index used here starts in 1957** (satellite/modern
  reanalysis era); earlier reconstructions exist in the literature back
  to ~1905 but were not sourced for this project (an earlier attempt to
  fetch a longer monthly reconstruction failed — see the
  conversation's own record of failed URL guesses). The continuity=1.00
  result is therefore validated only within the 1957-2026 window, the
  same caveat that applies to the North Sea/Baltic/Brest post-1950
  results discussed earlier.

## 5. Path forward

1. **Direct QBO-SAM comparison.** Run the same phase/amplitude
   extraction used for the Baltic-vs-MSL comparison, but between SAM's
   and QBO's own windings — do they share a phase relationship, or only
   a shared mechanism class?
2. **Check whether SAM's M=4.25 ridge (the strongest single ridge found
   this session) has a draconic-family origin**, the same way the
   Chandler wobble's 0.844 c/yr line was traced to the rectified
   draconic 2N torque and QBO's 0.422 c/yr line to the unrectified one.
   If M=4.25 converts to a real calendar period matching a known
   draconic alias, that would be a much sharper, falsifiable version of
   the hypothesis in §3 than category membership alone.
3. **Extend the ozone-hole trend observation.** The per-winding
   amplitude timeline (§1) is suggestive but was never formally tested
   against the published ozone-hole SAM-trend literature's own timing —
   worth a direct comparison rather than an eyeballed match.
