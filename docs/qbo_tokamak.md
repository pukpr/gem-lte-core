# The QBO as an axisymmetric toroidal device with a reversing poloidal drive (2026-10-22)

Companion to `landreman_lte_synthesis.md` section 8. Code: `qbo_tokamak.py` (figure `qbo_tokamak.png`, numbers
`qbo_tokamak.json`). Data: KIT radiosonde equatorial zonal wind at 30, 50, 70 hPa (`experiments/Oct2026_generator/legs`).

## 1. Geometry and correspondence

| Device | QBO |
|---|---|
| toroidal direction (around the major axis) | longitude around the equator; the QBO is the zonal mean (m = 0) |
| poloidal cross-section (minor radius r, poloidal angle) | latitude–height plane |
| toroidal field B_φ (B ↔ u) | zonal wind u |
| edge of the plasma | top of the QBO domain (~10 hPa) |
| inward penetration of a field change from the edge | downward propagation of a wind regime |
| external axisymmetric coils | zonal (m = 0) tidal triggers |

Axisymmetric equilibria obey the Grad–Shafranov equation; the fluid equivalent with swirl is the Bragg–Hawthorne
equation, with the angular momentum C(ψ) = r·u_φ in place of the toroidal-field function F(ψ) = R·B_φ. The state
reduces to the poloidal flux ψ and two flux functions, so no 3D (stellarator) shaping is needed for the QBO.

## 2. The device

**Drive.** A poloidal drive D(t) that reverses sign, applied at the edge. The tidal force has no zonal component
(section 8.3 of the synthesis), so the toroidal flow is not driven directly; it is generated from the poloidal drive.

**Edge relaxation (reversed-field-pinch threshold).** A tokamak's toroidal field is fixed by its coils and never
reverses; a reversed-field pinch generates its toroidal field and reverses it at the edge once the drive passes a
threshold (Taylor state B_φ ∝ J₀(μr), reversal for μa > 2.405). The device uses the simplest version, a threshold with
hysteresis on the normalized drive:

  u_edge = +1 when D > +θ, −1 when D < −θ, unchanged in between   (θ = 0: u_edge = sign D)

with a "linear" variant u_edge = D for comparison.

**Penetration (current penetration with an inward pinch).** The edge state moves inward (downward) by
advection–diffusion with damping:

  ∂u/∂t = −w ∂u/∂z + η ∂²u/∂z² − u/τ_d,   z = depth below the edge,

solved exactly by convolving u_edge with the first-passage Green's function
G(s; z) = z / √(4πηs³) · exp(−(z − ws)²/(4ηs) − s/τ_d). w is the analogue of the tokamak's inward (Ware) pinch, η of the
resistive diffusivity; a pure skin effect (w = 0) would give equal phase and amplitude decay lengths, which the QBO
does not have. Levels: 30, 50, 70 hPa at depths d, d + 3.5 km, d + 5.9 km (log-pressure heights).

## 3. The two drives

**Simplified cycle.** D = cos(2π(t − t_φ)/T) with T = 2.370 yr fixed in advance (below); only the phase t_φ is free.

**Comb impulse + modulation.** Impulses twice a year, at a trigger day of year and half a year later, with weight

  a_n = sin β_moon(t_n) · sin δ_sun(t_n) · (solar distance factor)³,

followed by a leaky integrator (time constant τ_i). β_moon is the Moon's **ecliptic latitude**, computed from Meeus'
mean elements, so its phase is astronomical, not fitted.

**Alias algebra (why 2.37 yr).** The draconic month is 27.212221 d, i.e. f_D = 13.4219 cycles/yr.
- Sampling twice a year with the **same** sign: a phase step of π·f_D = π·13.4219 per half year, i.e. 0.711 cycles,
  equivalent to −0.289 → an alias at **0.578 cycles/yr (1.73 yr)**.
- Sampling twice a year with **alternating** sign (the solar declination changes sign between the two triggers):
  multiplying by (−1)ⁿ adds π per step, π·(f_D + 1) = π·14.4219 → 0.211 cycles per half year → **0.422 cycles/yr
  (2.370 yr)**.
- **Modulation** sidebands: the 18.6-yr nodal cycle at 0.422 ± 0.054 cycles/yr (2.10, 2.72 yr) and the annual
  solar-distance asymmetry, which leaks the 1.73-yr line.

So the QBO's period follows from Doodson numbers alone, and **which** lunar argument matters is decided by it:

| lunar term in the comb | dominant month | alias with alternating semiannual comb |
|---|---|---|
| ecliptic latitude β (argument F) | draconic, 27.2122 d | **2.370 yr** |
| declination δ (dominated by ecliptic longitude λ) | tropical, 27.3216 d | 2.72 yr |

The declination's draconic part is only its inclination (±5.1°) modulation, i.e. the nodal sideband. The QBO driver
must therefore be the Moon's position relative to the **ecliptic** (the Sun's plane), not relative to the equator.

## 4. Results (fit 1953–2005, held out 2006–2026)

| Drive, edge | free parameters | train r (30/50/70 hPa) | held-out r | mean period at 30 hPa | descent lag 30→70 |
|---|---|---|---|---|---|
| simplified, threshold | 6 | 0.67 / 0.62 / 0.59 | +0.01 / +0.14 / +0.03 | 2.370 yr | 8 months |
| simplified, linear | 5 | 0.67 / 0.62 / 0.58 | −0.01 / +0.13 / +0.03 | 2.370 yr | −21 months |
| comb (β), threshold | 7 | 0.65 / 0.61 / 0.57 | −0.06 / +0.07 / −0.01 | 2.373 yr | 6 months |
| comb (β), linear | 6 | 0.67 / 0.63 / 0.59 | +0.06 / +0.17 / +0.04 | 2.367 yr | 8 months |
| comb (declination), threshold | 7 | 0.36 / 0.37 / 0.33 | −0.12 / −0.26 / −0.21 | 2.84 yr | 5 months |
| **observed** | | | | 2.36 yr (1953–2005) | 7 months |

- **Period and descent are reproduced.** The astronomically phased comb (only the trigger day free) fits 1953–2005
  as well as a sinusoid with a free phase, and gives the observed descent lag (6–8 months against 7).
- **The declination drive fails** (r ≈ 0.35), as the alias table predicts.
- **The threshold does not improve the correlation** over the linear edge; it shapes the regimes (flat-topped,
  sharp transitions) and kept the simplified model's descent the right way round (the linear simplified fit found a
  reversed-lag solution).
- **No held-out skill after 2006.** Against a 2.370-yr carrier, the observed QBO phase stays within about ±0.2 cycles
  from 1955 to 2005 (mean period 2.36 yr), then slips by about half a cycle across the 2016 and 2019–20 disruptions
  (phase +0.30 in 2005, +0.13 in 2020, +0.75 in 2025). A phase-locked device cannot follow a slip.
- **The disruptions are not explained by a weak drive.** In the threshold reading, a missed or extra reversal needs the
  drive near threshold. The comb drive's envelope in 2016–2019 is average to high (1.42–1.49 against a mean of 1.34;
  its lowest years are 2022, 1984, 2003). So the slips need another ingredient: a second trigger, a change in the
  edge threshold, or the QBO's known sensitivity to extratropical intrusions (the 2016 disruption is usually attributed
  to Rossby-wave momentum from the Northern Hemisphere).
- The fitted trigger days (235 and 53, late August and late February) trade off against the integrator lag τ_i, so
  they are not a clean measure of the physical trigger date.

## 5. What the design establishes

- An axisymmetric (tokamak-family) device suffices: the m = 0 triggers, an edge reversal and inward penetration
  reproduce the QBO's period, phase (to 2005), alternation and descent with 6–7 parameters, three of them physical
  transport coefficients.
- The **period is a prediction**, not a fit: 2.370 yr from the draconic month and an alternating semiannual comb,
  against an observed 2.36 yr. The same algebra rules out the declination (2.72 yr) and an unsigned comb (1.73 yr).
- The **reversal rule** (RFP threshold) is not yet tested by the data: it neither raises the correlation nor explains
  the post-2006 slips. Next tests: an edge threshold that depends on the drive's slow envelope (the RFP pinch
  parameter), and the 5-day ERA5 winds, which resolve the transition timing far better than monthly KIT data.

## 6. The QBO record read as an MHD shot (`qbo_tokamak_device.py`, `qbo_tokamak_device.png`, `qbo_tokamak_shot.png`)

**Rescaling.** 1 year of the record = 1 second of discharge; zonal wind = edge toroidal field, 1 m/s → 1 mT
(`qbo50.dat` is stored in units of 100 m/s and is rescaled accordingly).

| QBO | Shot |
|---|---|
| record 1953–2026 | a 73-s discharge (long pulse, but within reach of present long-pulse devices) |
| 30 / 50 / 70 hPa | magnetic probes at r/a = 0.90 / 0.80 / 0.75 |
| semiannual trigger | 2-Hz pulsed supply, alternating polarity |
| draconic month | 13.42-Hz amplitude modulation |
| QBO period 2.37 yr | field reversal every 2.37 s (0.42 Hz) |
| descent 30 → 70 hPa, 7 months | inward penetration 0.90 → 0.75 in ≈ 0.6 s |

**The reversing poloidal circuit** (panel c):
1. **Solar clock, the strong carrier comb.** It fires an H-bridge of thyristors twice per second and flips its polarity
   every firing, (−1)ⁿ: the sign of the solar declination at the two semiannual triggers.
2. **Lunar modulator.** It sets the capacitor bank's charge for each pulse, V_n ∝ sin β_moon(t_n) × sin δ_sun(t_n) ×
   (solar distance)³.
3. **L/R integrator** (τ_i): the solenoid's loop voltage drives plasma current and so the poloidal field. The pulse train
   integrates to the slow 2.37-s beat (inset).
4. **Edge relaxation:** past ±θ the edge toroidal field flips (RFP-like), creating a new reversal surface.
5. **Inward penetration** (pinch w, diffusion η, damping τ_d) carries each reversal surface past the probes.

**Why the ecliptic and not the equator.** The carrier comb is the strong impulse, and its timing and polarity are set
by the Sun. Each pulse samples the Moon at a moment fixed by the Sun's own plane, so what the comb measures is the
Moon's position relative to that plane: its ecliptic latitude, whose period is the draconic month. The declination
is measured from the equator, which is not the clock's plane; it is dominated by the tropical month and aliases to
2.72 s instead of 2.37 s. In circuit terms, the modulator must be referenced to the clock, not to an independent
frame.

**What the shot shows** (`qbo_tokamak_shot.png`, first 25 s). The measured probe map has the tilted reversal fronts
of inward penetration, the device reproduces their spacing and tilt, and over 1953–2005 its reversals stay in step
with the measured ones (r ≈ 0.6). After t ≈ 53 s (2006) the measured shot slips by about half a reversal cycle and the
device does not follow (section 4).

## 7. Shoulders: semiannual stalls from a gated (stroboscopic) penetration (`qbo_shoulders.py`, `.png`, `.json`)

**Observation.** The ragged square wave has shoulders: the wind holds at an intermediate level for months (10% of
months at 30 hPa, 29% at 50 hPa, 25% at 70 hPa), and month-to-month changes cluster by season. At 30 hPa the
semiannual component of |du/dt| is at the 100th percentile of surrogates (stronger May–July, weaker Nov–Jan); at
50 hPa the May–June peak dominates.

**Devices compared** (each refit on 1953–2005 with the objective [corr(u) + corr(3-month change)]/2):

| Device | Train r (30/50/70) | Train r, 3-month change | Seasonal \|du/dt\| profile r | Shoulder months | Held-out r |
|---|---|---|---|---|---|
| switch (threshold edge) | 0.62 / 0.62 / 0.57 | 0.51 / 0.49 / 0.42 | 0.53 / 0.04 / 0.10 | 0 / 0 / 0% | ≈ 0 to +0.2 |
| saturating core, unequal triggers | 0.66 / 0.63 / 0.60 | 0.51 / 0.50 / 0.41 | −0.09 / 0.71 / 0.16 | 2 / 0 / 0% | ≈ 0 |
| + satellite pulses (±31 d, 0.18) | 0.67 / 0.63 / 0.58 | 0.55 / 0.48 / 0.39 | 0.44 / 0.58 / 0.02 | 1 / 0 / 0% | ≈ 0 |
| switch + gated penetration | 0.64 / 0.62 / 0.56 | 0.53 / 0.53 / 0.42 | 0.78 / 0.93 / 0.67 | 0 / 6 / 1% | ≈ 0 |
| **saturating core + gated penetration** | **0.67 / 0.64 / 0.61** | 0.53 / 0.51 / 0.45 | 0.68 / 0.83 / 0.59 | 3 / 6 / 5% | ≈ 0 |
| observed | | | | 10 / 29 / 25% | |

- **Steps at the edge alone do not make shoulders.** The saturating core produces a staircase (volt-second steps on a
  saturable core, the weak trigger at about 0.6–0.7 of the strong one), but the inward diffusion washes the steps out
  before the probes. Satellite pulses add little.
- **Gating the penetration does.** When the inward transport (pinch and diffusion together) is strobed through the
  year, rate ∝ exp(m₁ cos 2π(t − p₁) + m₂ cos 4π(t − p₂)), the fronts stall between strobes and lurch forward at
  them. The fitted gate is mainly semiannual (m₂ ≈ 1.1) with an annual asymmetry (m₁ ≈ 0.3): a strong strobe around
  June and a weak one around December, matching the observed seasonal clustering (profile r 0.6–0.9 against ≈ 0.1 for
  the switch). The model now shows shoulders at the times the data do (e.g. 1989, 1991–92, 2001; figure).
- **Still short of the data.** The model's seasonal profile is too peaked, and it makes 3–6% shoulder months against
  10–29% observed: real stalls are longer and less regular than one fixed gate. Held-out skill after 2006 is still about
  zero (the phase slip of section 4).

**In circuit terms.** The strobe acts twice: it fires the H-bridge (the edge drive) and it gates the inward pinch,
strongly at one semiannual trigger and weakly at the other. The weak trigger's smaller advance is what leaves a
shoulder: a partial penetration that waits half a year for the next strobe.

## 8. The stellarator-like part: helical coils for the fast, high-winding ripple (`qbo_stellarator.py`, `.png`)

The aim here is the QBO's characteristics and overall timing, not a point-by-point match.

**What the 5-day ERA5 winds at 30 hPa show.**
- About 20% of the variance is faster than a year. The two-leg QBO fit (`legs/qbo30_point0E`) captures it with high
  windings of the common manifold, up to 26× the base (8×, 11×, 14×, 17×, 24×, 26×). Its sub-annual part correlates
  0.61 with the data's.
- The fast ripple (periods under about 4 months) is strongest near reversals (|u| < 5 m/s: 2.6 m/s) and while
  westerlies accelerate (2.2 m/s), and weakest inside settled westerlies (1.8 m/s).

**Device.** The axisymmetric device (section 7) plus a helical (m > 0) coil set, the stellarator-like part. The coils
are driven by the high windings of the fitted manifold, with amplitude gated by the slow state,
G = 1 + g₀·exp(−(u/5)²) + g₁·max(du/dt, 0), and g₀, g₁ taken from the data's state ratios (not fitted). In a
reversed-field pinch the helical modes likewise do their work at and around the reversal surface. Panel (a) of
`qbo_tokamak_device.png` draws the helical winding; panel (c) adds its driver.

**Characteristics, data vs emulator (30 hPa, 5-day, 1979–2026):**

| | log-spectrum r (0.1–30 cycles/yr) | fast ripple: westerly / easterly / transition / W-accel. / E-accel. (m/s) |
|---|---|---|
| observed | | 1.82 / 1.98 / 2.63 / 2.23 / 1.96 |
| axisymmetric only | 0.84 | 0.51 / 1.09 / 1.29 / 1.12 / 1.06 |
| + helical, no gate | 0.89 | 1.81 / 2.54 / 2.12 / 2.37 / 2.22 |
| + helical, gated | 0.90 | 1.77 / 2.47 / 2.27 / 2.43 / 2.17 |

Mean reversal period: 2.37 yr in the emulator against 2.43 yr in the 5-day record.

- **The helical coils supply the missing fast band.** The axisymmetric device has almost no power above about 2 cycles
  per year; the helical part fills 2–10 cycles/yr to the observed level and brings the ripple amplitude to the observed
  1.8–2.6 m/s.
- **The gate moves the state dependence the right way**, toward reversals and westerly onset, but only part way: the
  emulator's easterly ripple is still too strong.
- **The fastest band is still short.** Above about 10 cycles/yr (periods under ~5 weeks, where the equatorial
  Kelvin waves sit) the emulator's power falls an order of magnitude below the data; windings up to 26× do not reach it.
- **Artifact.** The gated penetration of section 7 is impulsive (a rate ratio of about 9 between strobe and off), so the
  emulator's reversals have sharper jumps than observed.
- **Not tested: the wavenumber itself.** The data are a single point at 0° E, so "m > 0" is an interpretation of the
  fast content, not a measurement. ERA5 winds at several longitudes along the equator would show whether the ripple
  travels eastward (Kelvin-like) as a helical mode should.
