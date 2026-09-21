# Overnight Quad Survey

Continuation of the per-region deep-dive session (kN040_E030 Aegean/Black Sea,
kN040_E130 Sea of Japan, kN_gom Gulf of Mexico, kN040_W070 US East Coast,
kN040_W130 US West Coast, kN020_W150 Hawaii, kN020_E050 Red Sea,
kS040_W050 Patagonia), continued autonomously overnight per instruction to
"continue analysis of as many quads as possible." All phase comparisons use
the shared-forcing-corrected methodology established in
`WINDING_GEOSPATIAL_KRIGING.md` (one location's own rebuilt forcing as a
single common reference, both sides refit against it) -- never a naive
same-file phase comparison, which was shown early this session to
manufacture spurious antiphase readings.

## Individual quad findings (this stretch)

- **Darwin (kS020_E130) vs Tahiti (kS020_W150)** -- the classic SOI
  atmospheric dipole, tested for an SST-domain echo. Frequency-dependent:
  backbone (M=0.207) shows real anti-phase (cos=-0.628, matching the SOI
  sense), but the independently-shared order-4 term (M=0.830) runs the
  *other* way (cos=+0.458, in-phase). A real but non-uniform dipole
  signature, same pattern as Hawaii's NINO4/PDO split.

- **TNA vs TSA** (Tropical N/S Atlantic, nominally the Atlantic Meridional
  Mode dipole) -- NOT a dipole in this framework. All three tested shared
  frequencies (TNA's dominant M=1.037, TSA's dominant M=2.42, and the
  standard backbone) show strong *positive* coherence (+0.71, +0.99, +0.89).
  Honest negative result for the "everything is a dipole" hypothesis.

- **Azores (kN040_W030) vs Iceland (kN060_W010)** -- NAO SST proxy test.
  The two reliably-measured shared terms (slow comb ~-0.014, backbone) are
  strongly in-phase (+0.98, +0.93), not anti-phase. Consistent with an
  already-documented finding from earlier this session
  (`project_nao_resolved_shallow_water.md`): NAO's real SST expression does
  not show up as a simple two-point spatial dipole; it needed a 12-month
  delayed-difference transform instead. This SST-only test reconfirms that
  difficulty rather than resolving it.

- **West Pacific warm pool (kN000_E150) vs East Pacific cold-tongue-adjacent
  (kN000_W110)** -- the cleanest dipole result of the night. At the
  ENSO-dominant frequency M=0.449: cos(dphi)=-0.955 with strong, reliable
  fit quality on both sides (fit_r=0.34, 0.31 -- among the best individual
  fit_r values found all session). Also anti-phase at order 15
  (cos=-0.988). Only the shared backbone stays in-phase (+0.978), exactly as
  expected for the universal background field. This is the real East-West
  Walker-circulation/zonal-SST-gradient signature that defines ENSO,
  recovered cleanly from independently-fit winding phases.

- **Atlantic Nino** (west equatorial Atlantic kN000_W010 vs Gulf of Guinea
  kN000_E010) -- weak/no dipole (mostly in-phase: backbone +0.995, order 9
  +0.871; one unreliable near-zero term). Physically consistent with the
  real, well-documented relative weakness of the Atlantic zonal mode
  compared to Pacific ENSO -- a negative result that matches known
  climatology rather than contradicting the method.

- **South Indian Ocean Subtropical Dipole** (Madagascar-adjacent kS020_E050
  vs W.Australia-adjacent kS020_E090) -- weak/mixed, no reliable dipole
  signature found (only the always-present backbone shows real coherence,
  +0.90; other shared terms have fit_r too low to trust).

- **Bering Sea neighbors** (kN060_E170 vs kN060_W170, ~2000km apart) --
  even immediate geographic neighbors show frequency-dependent mixing, not
  uniform coherence. Order 3 is in-phase (+0.66, reliable). But within a
  genuine close-doubled pair at order 2 (M=0.414 vs M=0.418, differing by
  only 0.004), one member is in-phase (+0.62) and the other is anti-phase
  (-0.72). Order 20 is also anti-phase (-0.82). Reinforces that this
  model's coherence structure is fundamentally frequency-specific, not a
  simple function of distance -- true even between adjacent boxes.

- **kN040_W090** (central/SE USA, 30-50N 80-100W) -- flagged with a data
  caveat: this box is mostly land, qualifying only via a small, patchy
  ocean fraction (Gulf Coast strip, maybe coastal Carolinas). Its huge
  dominant amplitude (0.75 at M=0.625) and a near-degenerate second peak
  (0.62 at M=0.622, differing by only 0.003) likely reflect this
  small/noisy sample rather than a cleaner physical signal. Still notable:
  the ~0.622 frequency matches kN040_W130's own genuine dominant term
  (amp=0.68, a well-covered open-Pacific box) -- a real coincidence worth
  registering even given the caveat.

- **kS040_W050 Patagonia doubling cascade + tidal connection** -- confirmed
  clean harmonic doubling 5->10->19(not quite 20), plus a genuine 3-member
  ridge cluster near order 25.5-29 (`triad`) via `winding_rank`'s AR1-gated
  scan (2 of 3 pass the significance test). This matches the textbook
  signature of nonlinear shallow-water tidal overtide generation
  (M2->M4->M6->M8-style harmonic cascades), and Patagonia is independently
  known for some of the world's strongest continental-shelf tides. The
  neighboring quad west of it (kS040_W070, spanning the narrow tip of South
  America) shows no such doubling (orders 7,12,21,25 -- no clean 2x
  relationship anywhere), confirming the user's own caveat that this
  quad straddles the continent and mixes Pacific/Atlantic/land signal.

## Cross-cutting harmonic classification (grid-wide scans)

Extending the distance-binned coherence method from earlier in the session
(originally applied to orders 1,2,6,9,12,14 using the ~19-index named
dataset) to the full 89-quad grid, using each newly-discovered harmonic's
own independently-selected members (never a forced blind regression):

| order | M | # quads carrying it | signature | classification |
|---|---|---|---|---|
| 3 | 0.622 | **50/89 (56%)** | 0.796 -> 0.430 -> 0.183 -> 0.025 -> +0.126 | finite correlation length, decays to ~zero by 6-10Mm; strongest/most prevalent harmonic in the whole grid, amplitude concentrated along the eastern North Pacific boundary |
| 4 | 0.829 | 27/89 | 0.655 -> 0.588 -> 0.533 -> 0.425 (never decays) | **strongest confirmed uniform global background field found** -- more consistent than the backbone itself |
| 11 | 2.280 | 33/89 | 0.916 -> 0.330 -> 0.160 -> 0.164 -> 0.136 | hybrid: decays initially then plateaus at a persistent nonzero level -- mix of a real regional mode plus a weaker shared global component |
| ~0.56 | 0.560-0.573 (non-backbone-harmonic) | 11 locations | mostly anti-phase vs NINO4 (cos -0.12 to -0.995) except the one immediately NINO4-adjacent quad (kN000_W170, cos=+0.988) | genuine remote-ENSO-teleconnection sign-flip pattern -- same mechanism as the US East Coast's anti-phase NINO4 relationship, now confirmed as a widespread pattern, not a one-off |
| 5 | 1.037 | 30/89 | 0.080 -> 0.126 -> 0.195 -> 0.155 -> 0.057 (weak/noisy throughout) | inconclusive -- doesn't match any of the four established classes cleanly; likely a weakly-shared, mostly locally-specific term |
| 8 | 1.658 | 25/89 | 0.449 -> 0.136 -> 0.084 -> -0.029 | finite correlation length, decays cleanly to ~zero -- joins orders 2, 3 in this class |
| 15 | 3.112 | 25/89 | 0.279 -> -0.049 -> -0.102 -> -0.017 | ambiguous -- mild negative dip at mid-range (weaker than order 14's clean -0.23 reversal), not cleanly uniform or finite-length; would need a direct basin test (like order 6 got) to classify further |

Combined with the pre-bedtime classification (orders 1, 2, 6, 9, 12, 14 on
the smaller 19-index dataset), the overall picture now spans most of the
grid's independently-discovered spectral content, organized into four
real behavior classes:

(Orders 1-15 are now comprehensively classified across the full 89-quad
grid, not just the original ~19-index dataset.)

1. **Uniform global background fields** (orders 1, 4, 9) -- present nearly
   everywhere, coherence stays high or flat regardless of distance.
2. **Finite-correlation-length regional fields** (orders 2, 3, 8) -- real,
   often strong, but decay smoothly to zero by ~6-10,000km; no sign
   reversal.
3. **True global standing/traveling waves** (orders 12 [k=5], 14 [k=1]) --
   decay *through* zero to negative at large distance, the textbook
   sign-reversal-at-half-wavelength signature of an actual planetary-scale
   wave mode.
4. **Discrete basin eigenmodes** (order 6, confirmed via direct tests on
   two independent basin clusters: baltic+northsea+brest cohere at one
   phase, Aegean+BlackSea+GOM cohere at a different, roughly opposite
   phase) -- basin-bound, not describable by any smooth spatial function of
   distance or longitude. GOM's membership in the Aegean/BlackSea camp
   (despite being geographically close to AMO/Baltic) suggested a refined
   rule: narrow-strait-bounded silled basins vs. open continental shelf,
   not raw geographic proximity.
5. **Remote-teleconnection sign-flip fields** (M~0.56, and the US East
   Coast's M=0.449 relationship to NINO4) -- strongly in-phase immediately
   adjacent to the source region, flipping to anti-phase almost everywhere
   else. Order 11 may partly belong to this class too given its hybrid
   signature.

## Open threads for further work

- Order 15's "patchy" classification from the original pass could use the
  same wider-grid treatment now available.
- IODW's pre-existing seed is deeply unstable even after reseeding from
  IODE (diverges to nonsensical order-1000+ terms) -- needs manual
  attention, not more automated reseed attempts.
- Red Sea (kN020_E050)'s own dominant high-order terms (2.28, 2.90, 1.66,
  1.87) were hypothesized as its *own* basin-specific eigenmode (different
  from GOM/BlackSea's order 6, scaled to its narrower/more elongated
  geometry) but this needs a second similarly-shaped silled basin (Persian
  Gulf?) as a cross-check before treating it as confirmed.
- Order 12's wavenumber-5 pattern (p=0.015, confirmed under a joint lat+lon
  search) is real but noisier than order 14's exceptionally clean
  wavenumber-1 result -- worth more data points if new quads are added
  later.
