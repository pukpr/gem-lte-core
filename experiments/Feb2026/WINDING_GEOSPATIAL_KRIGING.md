# Geospatial Structure of Winding Parameters, and What It Means for Kriging the Kaplan SST Grid

The winding triplets `(M, Amplitude, Phase)` that `lt.exe.windings.json`
now records (`GEM.LTE.Primitives.Shared.Save_Windings`, added this
session) are, in this model, the *only* thing that distinguishes one
climate index from another — the forcing manifold itself is shared
(`WINDING_NARRATIVE.md`). If that's true, winding-parameter coherence
between two sites should behave like a real spatial field: nearby sites
should share more of it than remote ones, in a way well-behaved enough to
support kriging/interpolation across the Kaplan SST grid
(`build_k_sst.py`). This document tests that directly, across every one
of the 19 indices with a `lt.exe.windings.json` as of this session.

## 1. A methodology bug had to be fixed first

The very first version of this analysis compared each index's own
`lt.exe.windings.json` phase directly against another index's own phase.
That is unsound: every index here has been through its own independent
Ada optimizer run since creation (confirmed directly — even the 42-term
LPAP Doodson tables have drifted apart from any common donor, not just
the Bessel-step `impA`/`impB`/`bg`/`offs`), so any two indices' manifolds
are typically only ~99.6% correlated overall — and that residual 0.4%
disagreement decorrelates **catastrophically** once multiplied by a
winding number `M` inside `sin(2*pi*M*Forcing)`. Caught directly on
`kN_baltic` vs. `kN_northsea`: their manifolds correlate 0.9966 overall,
but `sin(2*pi*M*F)` between them collapses to +0.31 at the shared
backbone (`M~0.21`) and goes **negative** at the dominant winding
(`M~1.25`) — manufacturing a spurious "antiphase" reading with no
counterpart in the real SST data at all (which shows R=0.974, strongly
in-phase, over 30 independent 20-year windows — see
`BALTIC_ADMITTANCE_LISSAJOUS.md` §8).

**Fix**: for every pair (A, B), rebuild A's own manifold once and use it
as the ONE shared reference forcing for both regressions — refit *both*
A's and B's raw F9-filtered data via the same single-mode OLS against
that one forcing, at the M values that match within tolerance (0.05).
Verified against a hand calculation before trusting it at scale: corrected
`kN_baltic`-`kN_northsea` phase agreement is `+0.999` (backbone) and
`+0.998` (dominant winding) — both essentially perfectly in-phase, fully
consistent with the raw SST correlation, once the artifact is removed.

**Also fixed for future use**: `lt.exe.windings.json` now carries a
`"manifold"` sub-object (`delA`, `delB`, `asym`, `ma`, `mp`, `shiftT`,
`impA`, `impB`, `offs`, `bg`, `year`) specifically so a future comparison
can check two files' manifold compatibility before trusting a direct
phase comparison, instead of rediscovering this bug a third time.

All results below use the corrected, shared-forcing methodology
(`geospatial_winding_analysis.py`, `winding_variogram.py`).

## 2. Case study: NINO4 vs. PDO

A detailed, single-pair worked example before the full-scale analysis,
since it's small enough to read term-by-term and it reproduces a result
`WINDING_NARRATIVE.md` had already documented independently (the
NINO4/PDO "near-degenerate beat pair," `M~0.449`/`M~0.415`):

| M | amp(nino4) | amp(pdo) | phase diff | relationship |
|---|---|---|---|---|
| 0.2073 (shared backbone) | 0.063 | 0.080 | +26° | in-phase |
| 0.4146 (2×backbone, beat-pair partner) | 0.133 | 0.104 | -46° | quadrature |
| 0.4491 (dominant mode, both indices) | 0.277 | 0.661 | -17° | strongly in-phase |
| 0.5604 (not a clean harmonic) | 0.126 | 0.106 | +105° | quadrature/weak antiphase |
| 0.6218 (3×backbone) | 0.210 | 0.292 | -31° | in-phase |

Unique to each: nino4 carries `M=1.385, 0.829, 3.109` (the last one
already documented as the tropical-instability-wave-scale harmonic); PDO
carries `M=1.866, 2.280, 4.146` (its own order-20 harmonic, also already
documented). Distance: 4477km.

Reading: the *low-order, clean-harmonic* windings (backbone and its low
multiples, plus the shared dominant mode) are strongly in-phase between
two indices 4477km apart — behaving like a real shared background field.
Each index's own *high-order* windings are entirely unshared — genuinely
local fine structure. The one exception, `M=0.560` (shared in frequency
but only weakly coherent, and notably *not* a clean multiple of the
backbone), hints that coherence tracks harmonic-family membership, not
just raw frequency agreement.

## 3. Systematic classification: the harmonic ladder, at full scale

Every winding across all 19 indices, classified as harmonic-family
(within tolerance of an integer multiple of the shared backbone,
`M(NM)~0.2075`) or not:

| harmonic order | M | n indices carrying it (of 19) |
|---|---|---|
| 1 (backbone) | 0.207 | **19** (universal) |
| 2 | 0.415 | 10 |
| 3 | 0.622 | 11 |
| 4 | 0.830 | 6 |
| 5 | 1.037 | 6 |
| 6 | 1.245 | 7 |
| 7 | 1.452 | 3 |
| 8 | 1.660 | 4 |
| 9 | 1.867 | 4 |
| 10 | 2.075 | 2 |
| 11 | 2.282 | 6 |
| 12 | 2.490 | 2 |
| 15 | 3.112 | 7 |
| 16, 18, 20, 21, 26, 29 | 3.32-6.02 | 1-3 each |

142 total winding instances, 68% harmonic-family / 32% non-harmonic. The
order-1 backbone is literally universal (found in every single index that
has a `lt.exe.windings.json`); membership thins out steadily with order,
and several of the highest orders are singletons — genuinely
index-specific fingerprints, not shared field.

**Practical rule this gives the Kaplan SST kriging plan**: winding order
is a principled way to separate "global background field" (low order,
worth interpolating spatially) from "local fine structure" (high order,
index-specific, should not be interpolated — treat as a separate
per-cell term).

## 4. The empirical variogram: real near-field decay, noisy far field

Testing whether coherence decays smoothly with distance (the property a
real kriging variogram needs), across all 249 matched-winding instances
from 91 point-like index pairs, split by harmonic-family membership:

| distance bin (km) | harmonic: n, mean cos(dφ) | non-harmonic: n, mean cos(dφ) |
|---|---|---|
| 0-2000 | 6, **+0.990** | 0, — |
| 2000-5000 | 48, +0.417 | 13, +0.439 |
| 5000-8000 | 37, +0.282 | 8, +0.298 |
| 8000-12000 | 56, +0.397 | 15, +0.101 |
| 12000-16000 | 42, +0.593 | 14, +0.190 |
| 16000-20000 | 9, -0.164 | 1, +0.911 (n=1, noise) |

Overall correlation(distance, coherence) is weak for both groups
(harmonic: r=-0.037; non-harmonic: r=-0.119) — **not** because there's no
spatial structure, but because the relationship isn't a simple smooth
global decay. The plot shows why directly: a sharp drop from ~+1.0 to
~+0.4 in the first ~3000km (a real nugget-to-sill transition), then a
noisy, non-monotonic plateau from 3000-16000km where individual pairs
scatter across the *entire* range from -1 to +1 regardless of exact
distance. The 16000-20000km tail is too thin (n=9 and n=1) to read.

## 5. What actually explains the far-field scatter: ocean-basin membership, not raw distance

Individual case studies already hinted at this — `kN_okhotsk` is
strongly *anti-phase* with `nino4` (6891km) yet only *quadrature* with
the much more distant `tna` (15212km) — geometric distance alone doesn't
predict this, but oceanographic disconnection does (Okhotsk is a
sub-arctic NW-Pacific marginal sea; nino4 is equatorial open Pacific;
neither shares a current system with the other or with the tropical
Atlantic in any direct way).

Grouping all 14 point-like indices into ocean-basin categories (N.
Atlantic, European shelf seas, Mediterranean/Black Sea, Equatorial
Pacific, N. Pacific, NW Pacific marginal seas, Indian Ocean, Southern
Ocean) and comparing same-basin vs. different-basin pairs directly:

```
SAME-basin pairs (n=7):  mean coherence = +0.820   median = +0.979
DIFF-basin pairs (n=84): mean coherence = +0.383   median = +0.387
```

More than double. And critically, this survives controlling for
distance — restricting to the 2000-9000km range where both groups have
enough pairs to compare fairly:

```
MATCHED DISTANCE (2000-9000km):
  SAME-basin (n=4): mean = +0.692
  DIFF-basin (n=37): mean = +0.312
```

Same-basin pairs remain more than double the cross-basin mean even at
matched distance — basin membership carries real explanatory power
*beyond* raw geographic distance, not merely a proxy for it. (One
exception worth flagging honestly: `kN040_E130`-`kN_okhotsk`, both
classified NW-Pacific-marginal, scores essentially incoherent, +0.008 at
2229km — basin membership is a strong predictor, not a perfect one.)

## 6. Verdict and practical recommendation

| Claim | Status |
|---|---|
| Winding parameters are the only thing distinguishing indices in this model | Assumed, per `WINDING_NARRATIVE.md`; consistent with everything below |
| A shared-forcing correction is required before comparing any two indices' winding phases | Confirmed — naive comparison manufactures spurious antiphase readings |
| Windings decompose cleanly into a universal low-order "background" and index-specific high-order "fine structure" | Confirmed — order-1 backbone in 19/19 indices, several high orders in exactly 1 |
| Coherence decays smoothly with distance globally (a simple isotropic variogram) | **Not confirmed** — real short-range decay (<3000km) but noisy, non-monotonic beyond that |
| Ocean-basin membership explains coherence better than raw distance | Confirmed — same-basin pairs score >2x cross-basin coherence, even at matched distances |

**For the Kaplan SST kriging plan specifically**: don't fit a single
global isotropic variogram to raw great-circle distance — the far-field
scatter shown in §4 means that would mostly be fitting noise. Two
things are load-bearing instead: (1) restrict the "global background
field" being interpolated to the low harmonic orders identified in §3 —
that's the part with real spatial coherence; (2) use ocean-basin
membership (or a real oceanographic-connectivity measure — current
systems, strait connections — rather than Euclidean distance) as a
categorical covariate, not just raw distance, the way §5 demonstrates is
necessary. Given `build_k_sst.py`'s own 20°×20° grid, adjacent cells will
mostly fall in the strongly-coherent <2000-3000km regime documented in
§4, which is exactly where this approach should work well; interpolating
confidently across ocean-basin boundaries is not supported by anything
found here.
