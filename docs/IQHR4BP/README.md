# Invariant tori and stroboscopic maps: the IQHR4BP and the GEM-LTE equatorial ladder

Revised 2026-10-27 against the paper itself. The first draft is kept as `README_original.md`; §6 lists what was
corrected and why.

**Paper:** G. Grossi, N. Baresi, F. Topputo, *Quasi-Periodic Invariant Manifolds in the In-Plane Quasi-Hill
Restricted Four-Body Problem*, AAS 26-622 (2026). Local copy:
[GrossiBaresiTopputo2026…pdf](GrossiBaresiTopputo2026Quasi-PeriodicInvariantManifoldsintheIn-PlaneQuasi-HillRestrictedFour-BodyProblem.pdf).
(Politecnico di Milano / University of Surrey. The draft's "Grossi & Ross 2024" is wrong.)

## 1. What the paper actually does

The paper is a **numerical-methods paper**. It computes invariant tori for spacecraft orbits in the Earth–Moon system
and continues them through a hierarchy of models. It makes no claim about geophysics.

| model | Earth–Moon motion | external frequencies | spacecraft torus D = Di + De |
|---|---|---|---|
| CR3BP | circular, fixed | none | 1 = 1 + 0 (periodic orbit) |
| HR4BP | Hill's variational orbit (HVO) | solar tide, period **½ synodic month** | 2 = 1 + 1 |
| IQHR4BP | quasi-Hill variational orbit (QHVO), **planar** | ½-synodic + **anomalistic** (27.56 d) | 3 = 1 + 2 |

- **The Sun enters through Hill's approximation as a uniform tidal field.** In the Earth–Moon rotating frame it
  repeats every **π TU_S, half a synodic month (14.77 d)**: the MSf period, not an annual one. The half-period
  comes from symmetry: the solar tidal field is symmetric about the barycentre, and the HVO is invariant under a
  half-turn.
- **Hill's variational orbit has two center manifolds.**
  - The in-plane one has a period close to the **anomalistic month**.
  - The out-of-plane one has a period close to the **draconic month** (Aydin 2023, refs 25, 55).
  - So Hill's linearization alone produces the Moon's perigee and node rates. Their slow differences from the
    sidereal month are the 8.85-yr and 18.6-yr cycles: emergent beats, not separate forcings.
- **The IQHR4BP excites only the in-plane manifold.** The QHVO is the member of that family with e = 0.055, so the
  model has lunar eccentricity (Earth–Moon pulsation) but **no inclination and no draconic month**. The paper says
  so explicitly: the HVO "does not capture the inclination of the lunar orbit and the associated draconic month",
  and the model is named *In-plane*.
- **Continuation.**
  - *External* continuation adds a forcing frequency: m from 0 to m_SEM = 0.08085 (synodic month / year), then
    e_M from 0 to 0.055.
  - *Internal* continuation moves along a family within one model: HVO → QHVO, parameterized by the internal
    frequency.
- **Method.** GMOS-type two-point boundary-value problem for the **invariant section of a stroboscopic map**.
  - Adaptive Chebyshev–Gauss pseudospectral collocation propagates each point of the section for one
    stroboscopic period.
  - The Fourier order of the section adapts as the torus deforms.
  - Unfolding parameters square the Newton system.
  - Pseudo-arclength continuation, with up to ~200 000 unknowns.

## 2. The stroboscopic map: the central equivalence

The paper's definitions (Eqs. 4–7):
- Torus function τ(θ) with invariance φ(τ(θ), t) = τ(θ + ωt).
- Pick one angle θ_s as the strobe. Sampling every T_s = 2π/ω_s reduces the D-torus to a (D−1)-dimensional
  invariant section.
- On that section, the remaining angles advance by the **rotation vector ρ = ω T_s** per stroboscopic step.

GEM-LTE does the same thing with the **annual impulse comb** as the strobe:

| IQHR4BP / torus language | GEM-LTE language |
|---|---|
| strobe angle θ_s, period T_s | annual phase, T_s = 1 yr (the `Impulse_Delta` comb) |
| reduced angle advanced each strobe | draconic phase F, advanced 2π·13.4220 per year |
| rotation number ρ/2π = ω·T_s mod 1 | **q = 13.4220 − 13 = 0.4220 /yr → 2.370 yr, the QBO** |
| invariant section σ(ϑ) | the forcing manifold M(t) read at the impulses |
| torus function τ(θ) = Σ c_j e^{i j·θ} (Eq. 10) | LTE response sin(k M(t)): a Fourier series on the lunar–annual angles; by Jacobi–Anger its coefficients are Bessel functions of k × the tidal amplitudes |
| frequencies j·ω of the torus spectrum | the alias "comb": q, 2q, 1 ± q, 2 ± q, and the anomalistic alias 13.2552 − 13 = 0.2552 /yr (3.92 yr) |
| choice of strobe changes the reduced description, not the torus | changing the impulse timing changes which aliases appear, not the tidal input |

So in torus terms the QBO period is **not a beat between an internal mode and the forcing**. It is the rotation
number of the draconic angle on the annual stroboscopic section. That is a Fourier mode j = (−13, 1) of a 2-torus
(annual, draconic). It is not a third frequency.

**Half-wave symmetry, in both.**
- *Paper:* the QHVO's Fourier coefficients with even index along θ₁ are exactly zero, inherited from the
  half-synodic symmetry of the HR4BP. Only odd coefficients are stored.
- *GEM-LTE:* the QBO comb is odd-harmonic: twice-yearly impulses of alternating sign (`asym`), so annual
  fundamental with odd harmonics only.
- *Why it matters:* 13 is odd, so the draconic line folds to q = 0.4220 /yr. A plain, same-sign semiannual comb
  would fold it to |13.4220 − 14| = 0.578 /yr (1.73 yr), which is not the QBO.

## 3. Where the two ladders correspond, and where they don't

**Correspondences that hold:**
- **The lunar orbit itself.** The QHVO (D = 2: synodic, anomalistic), extended by the out-of-plane (draconic)
  manifold, *is* the quasi-periodic lunar motion that generates the tidal forcing. The GEM-LTE forcing manifold is
  a function on the same angles (draconic, anomalistic, tropical, annual), calibrated to dLOD. The strongest
  equivalence is therefore **QHVO ↔ forcing manifold**, not spacecraft torus ↔ QBO.
- **Center manifolds ↔ tidal species:**
  - *In-plane (anomalistic, Earth–Moon distance):* Mm, 27.55 d, and its annual alias 0.2552 /yr.
  - *Out-of-plane (draconic, ecliptic latitude):* the draconic month, whose annual alias q is the QBO.
  - *Missing in the IQHR4BP:* the paper's model has no out-of-plane part. The atmospheric analogue of the QBO
    forcing needs the **spatial** extension (a "QHR4BP" with inclination), which the paper doesn't treat.
- **External continuation ↔ adding tidal constituents.** Each added incommensurate forcing frequency raises De, and
  with it the torus dimension. That is the paper's De ≤ D ≤ De + n.
- **Solar forcing at ½ synodic ↔ MSf.** The paper's single solar frequency is the 14.77-d MSf period in the
  rotating frame. In ERA5 the ±MSf sidebands were traced to 4D-Var assimilation seams (`docs/era5/README.md`), so
  that line has not been found in the atmosphere as a physical forcing.

**Differences that matter:**
- **Hamiltonian vs dissipative.**
  - *Paper:* the IQHR4BP is conservative. Every torus keeps an internal frequency (Di = 1) and comes in a Cantor
    family.
  - *GEM-LTE:* the atmosphere and ocean are damped (the IIR stage). Free internal oscillations decay, and the
    forced response settles onto a **slaved torus with Di = 0**: a graph over the external angles alone. That is
    why one autonomous manifold can be shared across indices.
  - Counting by the paper's rule: annual alone gives D = 1, De = 1; adding the draconic month gives D = 2, De = 2,
    with the QBO as a Fourier mode of that 2-torus.
- **Response depth.** The spacecraft torus is in the paper's 6-D state space. In GEM-LTE the analogue is the
  wavenumber k in sin(k M(t)), a spatial or modal parameter, not a frequency. The "windings" are harmonics of k on
  a fixed manifold, which corresponds most closely to the Fourier order M_F of the torus function, not to new torus
  dimensions.
- **Resonance.** The paper assumes incommensurate frequencies; tori break near resonances. q = 0.4220 is not near a
  low-order rational (5/12 = 0.4167, 3/7 = 0.4286), consistent with a persistent QBO period.

## 4. Corrected ladder

```
Orbital (Grossi, Baresi, Topputo 2026)          Equatorial atmosphere (GEM-LTE reading)
--------------------------------------          ---------------------------------------
CR3BP periodic orbit   D=1 (Di=1, De=0)          free equatorial waves: damped, no torus survives (Di→0)
  │ external: Sun, period ½ synodic (MSf)          │ external: annual (odd-harmonic, alternating-sign comb)
HR4BP 2-torus          D=2 (Di=1, De=1)          annual/semiannual response (SAO)   D=1 (Di=0, De=1)
  │ external: anomalistic (e_M 0 → 0.055)          │ external: draconic month
IQHR4BP 3-torus        D=3 (Di=1, De=2)          QBO: rotation number q of draconic on annual strobe   D=2 (Di=0, De=2)
  │ (not in paper) out-of-plane / draconic         │ external: anomalistic, tropical (Mm, Mf) …
spatial model, then ephemeris                    full tidal manifold (ENSO, AMO … share it)   D=De=n
```

The ladders run parallel in the **external** frequencies, not row by row. The atmospheric column has no internal
frequency. Its key step, the draconic month, is the one the paper leaves out; the paper's anomalistic step comes
later in the atmospheric ladder.

## 5. Precise statements of the facts the draft used

- **Draconic month vs declination.**
  - The draconic month (27.2122 d) is the period of the Moon's **ecliptic latitude** (node to node).
  - Declination is the sum of an obliquity term (ε = 23.44°, tropical month 27.3216 d) and an inclination term
    (i = 5.14°, argument of latitude, draconic).
  - So the draconic month modulates declination but does not by itself "control" it.
  - The equatorially symmetric (k = 0) tidal term goes as δ², which puts it at fortnightly periods (Mf). The draconic
    fundamental, linear in latitude, sits in the antisymmetric part.
  - Hourly ERA5 shows no Mf / Mm / MSf line at 30–50 hPa (`docs/era5/README.md`, strat_tides). That is consistent
    with the forcing acting through the impulse comb rather than as a direct spectral line.
- **18.6 yr and 8.85 yr** are the node and perigee rates: the beats of draconic and anomalistic against the
  sidereal/tropical month (18.6 and 8.85 yr from the month lengths above). They are not independent forcings. The
  primary terms are the monthly and fortnightly ones (Mm, Mf).
- **Chandler wobble.** The repo reads F_CW = 2q (2·13.4220 − 26 = 0.8440 /yr), the same stroboscopic folding at
  twice the rotation number (`docs/qbo_six_year_beat.md`). This is the repo's hypothesis, not something in the
  paper.

## 6. Corrections to the original draft (`README_original.md`)

| draft said | correction |
|---|---|
| Grossi & Ross (2024), ResearchGate 412048007 | Grossi, Baresi & Topputo, AAS 26-622 (2026); the ResearchGate link was not checked |
| IQHR4BP adds "eccentricity + inclination", "anomalistic + draconic" | in-plane only: eccentricity / anomalistic month; inclination and the draconic month are explicitly excluded |
| HR4BP adds "solar forcing" as an annual-type frequency | the solar forcing period in the model is ½ synodic month (MSf, 14.77 d) |
| in-plane manifold ↔ SAO; out-of-plane ↔ draconic forcing | in-plane ↔ anomalistic (Mm); out-of-plane ↔ draconic; the SAO is the annual strobe and its forced harmonic, not a center manifold |
| ω_QBO = \|p ω_SAO − q ω_drac\| ≈ \|13 − ω_drac\| | mixes units and uses the wrong comb: q = 365.2422/27.2122 − 13 = 0.4220 cycles/yr; the 13 is the 13th *annual* harmonic, reached through an odd-harmonic comb; a semiannual comb gives 0.578 /yr |
| QBO is a 2-torus "beat mode" from SAO (internal) + draconic | QBO is the rotation number of the draconic angle on the annual stroboscopic section: a Fourier mode of a 2-torus with no internal frequency |
| table "mapping is exact"; "exactly the same mathematical structure" | the shared structure is torus + stroboscopic map + continuation; the ladders differ in Di (Hamiltonian vs damped) and in order (draconic is absent from the paper) |
| full ephemeris / "HFEM" step and "ELFO draconic sidebands" as paper content | not in the paper; the full-ephemeris step is stated as future use, and ELFOs are not discussed |
| full atmospheric ephemeris adds "18.6-yr nodal, 8.85-yr perigee" forcing | these are beat rates of the monthly terms (Mm, Mf, draconic), not separate forcings |
| draconic cycle "controls lunar declination" | the draconic cycle sets ecliptic latitude; declination's main monthly cycle is tropical |
