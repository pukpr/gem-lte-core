# The QBO as a switched toroidal plasma (2026-10-22)

Figure: `qbo_switching_plasma.png` (script `qbo_switching_plasma.py`). The detailed working is in `qbo_tokamak.md`;
this note is the short version.

## The claim

The QBO behaves like a ring of conducting fluid driven by a **polarity-switching pulsed supply**: the flow around
the ring reverses every 2.37 cycles of the switch, and each reversal sinks from the top of the ring into its interior.
Read with 1 year = 1 second, the equatorial stratosphere is a 73-second discharge whose edge field flips every
2.37 s.

## Three ingredients

1. **A switch (the Sun).** It fires twice per cycle and flips polarity each time (+, −, +, −): the semiannual
   trigger, signed by which hemisphere the Sun is over.
2. **A modulator (the Moon).** It sets each pulse's height: how far the Moon is above or below the Sun's plane (its
   ecliptic latitude), which repeats every draconic month, 13.42 times per cycle.
3. **A ring that reverses.** The ring's flow (or field) flips sign at its edge once the integrated drive passes a
   threshold, and each new reversal moves inward.

## Why the period is 2.37, with nothing tuned

This is stroboscopic beating, the wagon-wheel effect. Flipping the sign at every pulse shifts the 13.42 modulation to
14.42; sampling it twice per cycle folds 14.42 down to 14.42 − 14 = **0.42 cycles, a 2.37-cycle beat**: the QBO's
28-month period. Two controls fail, which is what makes the result specific:
- without the polarity flip the beat is 1.73;
- with the Moon's declination (its height above the equator) instead of its latitude it is 2.72.

The switch sets the reference plane, so the modulator has to be measured from the Sun's plane.

## Systems that already do each part

| Ingredient | Where it exists |
|---|---|
| Polarity-switched (AC) drive of a toroidal plasma | **AC tokamak operation**: the plasma current is reversed from one half-cycle to the next with the plasma sustained. Demonstrated on JET (Tubbing et al., *Nucl. Fusion* 1992) and on smaller tokamaks (STOR-M, CT-6B, ISTTOK, the last running many AC cycles in a row). |
| A ring whose edge field reverses by itself past a threshold | **Reversed-field pinch** (RFX-mod in Padova, MST in Madison): the toroidal field reverses near the edge on its own once the driving current is large enough; the relaxed state is Taylor's minimum-energy state (Taylor, *PRL* 1974), with reversal when μa > 2.405. |
| Sustaining a toroidal plasma with oscillating, sign-reversing voltages | **Oscillating-field current drive (OFCD)**: AC poloidal and toroidal loop voltages, 90° apart, inject magnetic helicity and drive net current. Proposed by Bevir and Gray (1980), tested on ZT-40M and on MST (McCollam et al., *PRL* 2006). |
| Descending flow reversals driven by an oscillating boundary | **The QBO laboratory analogue**: Plumb and McEwan (*J. Atmos. Sci.* 1978) produced descending, periodically reversing mean flow in a stratified annulus forced by a standing wave; Semin, Garaud and Fauve (*PRL* 2018) revisited it. This is a fluid, not a plasma, but it is the same ring geometry and the same reversals. |

So every element is realized somewhere. What has not been built is the combination: an AC- or OFCD-driven
reversed-field pinch whose switching is amplitude-modulated off-resonance, so that the edge reversal follows the
stroboscopic beat rather than the switch itself.

## What a lab test would look like

- Drive an RFP (or AC tokamak) with a pulsed supply that flips polarity twice per cycle, at, say, 2 Hz.
- Modulate each pulse's height at 13.42 Hz.
- Expected: edge-probe reversals every 2.37 s, and later at probes further in. Remove the polarity flip and the
  reversals should come every 1.73 s instead.
- The test is cheap in principle, because only the timing of the supply changes; the plasma is standard.

## Limits

- The fit to the real QBO holds its phase from 1953 to 2005 (correlation about 0.6 at 30–70 hPa, descent lag 6–8 months
  against 7 observed) but not across the 2016 and 2019–20 disruptions (`qbo_tokamak.md`, section 4).
- In the real QBO the mechanism is usually attributed to upward-propagating waves depositing momentum (Lindzen and
  Holton 1968; Baldwin et al., *Rev. Geophys.* 2001). The switched-plasma picture replaces the waves' timing with the
  Sun–Moon triggers; it does not show that the waves are absent.
- The details (shoulders from strobed penetration, a fast helical ripple) are in `qbo_tokamak.md`, sections 7–8.

## Phase-locking test (`qbo_phase_lock.py`, `qbo_phase_lock.png`)

Reversal dates from the KIT winds at 30, 50 and 70 hPa (3-month mean, ±3 m/s hysteresis), 1953–2024: 59–64
reversals per level. Each reversal's phase in the 2.370-yr beat (easterly onsets shifted by half a cycle); Rayleigh R
measures how tightly they cluster. Null: the same intervals between reversals in random order (same period spread,
phase free to wander), 5000 draws.

| Level | R at 2.370 yr | Shuffled median (95th pct) | p | Best period in a scan (R) | Scan max, p | R at 1.73 / 2.72 yr |
|---|---|---|---|---|---|---|
| 30 hPa | 0.61 | 0.31 (0.57) | 0.027 | 2.345 yr (0.70) | 0.041 | 0.05 / 0.13 |
| 50 hPa | 0.52 | 0.19 (0.39) | 0.002 | 2.350 yr (0.57) | 0.020 | 0.09 / 0.09 |
| 70 hPa | 0.44 | 0.14 (0.28) | 0.0006 | 2.350 yr (0.48) | 0.0008 | 0.04 / 0.08 |

- **The reversals are phase-locked**, beyond what the same reversal intervals in random order give, at all three
  levels (p from 0.03 at 30 hPa to 0.0006 at 70 hPa). The QBO's phase is regulated: it does not drift like a free
  oscillator with its own period spread.
- **"2.35 yr" is not a distinct period.** Over the record the difference from 2.370 amounts to about a quarter of a
  cycle (30.4 against 30.7 cycles for 1953–2024; 31.2 against 31.4 to late 2026), comparable to the phase jitter
  itself, so where the cycle is counted from and to decides it. Restricted to 1953–2005, the best period is
  2.356 / 2.363 / 2.362 yr at 30 / 50 / 70 hPa, with R at 2.370 essentially equal to the best (0.71 vs 0.73,
  0.66 vs 0.67, 0.57 vs 0.57): no measurable departure from 2.370 over five decades. The pull toward 2.345–2.35 comes
  entirely from the reversals after 2006, the decade with the 2016 and 2019–20 disruptions, i.e. a phase step, not a
  different steady period.
- **The lock weakens but survives after 2006** at 30 hPa (R 0.70, phase shifted from 0.54 to 0.32 cycles); at 50 and
  70 hPa it falls to R 0.44 and 0.33.
- **What the controls do and do not show.** R at 1.73 and 2.72 yr is near zero, but any oscillator with a period near 2.37 yr
  would fail those; they rule out those periods, not other mechanisms.
- **What this does not yet show.** Any phase-regulated oscillator near 2.37 yr would pass this test, including a
  wave-driven QBO whose phase is regulated by the seasonal cycle. Seasonal regulation with half-year steps and a
  stroboscopic beat are close relatives: a 2.37-yr average built from half-year steps is exactly what the beat
  produces. The discriminating test is the **absolute phase**: with the trigger dates fixed by astronomy (not fitted),
  the strobe model predicts *where* in the beat the reversals fall. That prediction, checked against the observed mean
  phases above, is the next step.

## Absolute-phase test (`qbo_abs_phase.py`, `qbo_abs_phase.png`)

The drive is fixed by the ephemeris alone: impulses at the June and December solstices, heights from the Moon's ecliptic
latitude, the Sun's declination and its distance. Nothing is fitted.

**1. The astronomical beat is almost a perfect clock.** Its phase wobbles by only 0.026 cycles over 1953–2024
(sd 0.006) and its amplitude by ±4%. The observed reversal phase wanders by about 0.2–0.3 cycles. So:
- the lag-independent test (do reversals follow the ephemeris wobble?) has no power: the ephemeris-aligned phase
  gives R 0.617 / 0.523 / 0.448 at 30 / 50 / 70 hPa, the same as time-shifted ephemerides (34th–52nd percentile);
- the QBO's phase wander is not astronomical. In the gating picture it is the internal (wave) part: jitter in when
  the system reaches its threshold.

**2. Implied lag.** Taking the drive's beat maximum as the edge trigger, the observed reversals follow it by:

| | 30 hPa | 50 hPa | 70 hPa |
|---|---|---|---|
| drive sign as written | 21.8 months | 27.0 | 30.1 |
| drive sign reversed | 7.6 months | 12.8 | 15.9 |

(beat = 28.4 months). The lags increase downward, by 5 and 3 months, matching the observed descent. A descent of
about 1 km per month puts the trigger near 10 hPa for the reversed sign (7.6 months above 30 hPa) or near 2–3 hPa
for the sign as written (22 months). Both are plausible places for a QBO to start, and the drive's sign is not yet
fixed by physics. So the absolute phase is **consistent** with the strobe, but it does not discriminate: a random
phase would also land on a plausible lag with good probability.

**3. What would discriminate.** A wave-driven QBO regulated by the seasonal cycle and a lunar-gated QBO make the same
half-year-stepped pattern and the same spectral images (f₀ ± 1, 2 − f₀). They differ in the **long-run mean period**:
the gated QBO must average exactly 2.3702 yr indefinitely; a wave-driven one has no reason to, and climate models
project its period to change. Over the record the difference between 2.370 yr and the best-fitting period is about a quarter of a cycle, all of it
from a phase step after 2006; over 1953–2005 the record is indistinguishable from 2.370. Not a contradiction, and not
yet a confirmation. Two ways to
sharpen it:
- fix the drive's sign from the physics (which sign of the Moon–Sun product pushes westerly), which turns the lag into
  a prediction of the trigger altitude;
- every further decade of reversals tightens the mean-period test; a QBO that keeps drifting from 2.370 would rule the
  strobe out.

## The SAO above: the clock made visible (`qbo_sao.py`, `qbo_sao.png`)

ERA5 monthly equatorial zonal wind (zonal mean, 5°S–5°N), 1–70 hPa, 1979–2025.

**1. The SAO is calendar-locked, so it is forced, not resonant.** At 1–3 hPa the mean annual cycle carries 77 / 67 / 48%
of the wind variance, and the semiannual harmonic's phase repeats from year to year with R = 0.98 / 0.97 / 0.94. The
lock fades below 5 hPa (0.69) and is gone by 10–30 hPa (0.16–0.24).

**2. Its two cycles are unequal, and that inequality is the polarity flip.** At 1 hPa the January easterlies reach
−39 m/s and the July easterlies −20; the April westerlies +32 and the October westerlies +23. A same-sign semiannual
strobe would give a 1.73-yr beat, which the data reject; the inequality between the two half-years is an annual
component, and an annual component has opposite sign half a year apart. That is the alternating polarity the 2.37-yr
beat needs. Physically it is the hemispheric (and perihelion, January) asymmetry of solar heating, which the drive
already carries through sin δ_sun and the solar distance.

**3. The 2.37-yr lock is a lower-altitude property.** Reversal clustering at the beat (Rayleigh R) by level:
1 hPa 0.03 (none), 2–3 hPa 0.33–0.39, 5–10 hPa 0.44–0.47, 20–70 hPa 0.47–0.58. Above, the calendar dominates; below,
the beat. (ERA5's 47 years give about 40 reversals per level, so the shuffle tests here are weaker, p 0.03–0.19, than
for the 71-year KIT record.)

**4. Descent and the trigger.** Reversal timing relative to 30 hPa (months), continuous from 70 hPa upward:

| hPa | 2 | 3 | 5 | 7 | 10 | 20 | 30 | 50 | 70 |
|---|---|---|---|---|---|---|---|---|---|
| months | −16.8 | −16.4 | −12.2 | −9.7 | −6.9 | −2.8 | 0 | +5.1 | +7.5 |

Reversals first appear at 2–3 hPa, at about the same time at both levels, then descend at roughly 1 km per month.
This settles the drive's sign from the absolute-phase test:
- **reversed sign:** the trigger would fall 7.6 months above 30 hPa, between 7 and 10 hPa, after the reversals have
  already appeared above it, so it cannot be the trigger;
- **sign as written:** the trigger falls 21.8 months before 30 hPa, about 5 months before the reversals first appear at
  2–3 hPa, i.e. at or just above the top of the QBO, beneath the SAO layer.

So the data favour the written sign, with the drive's beat maximum acting at the SAO–QBO interface and the reversal
emerging about 5 months later at 2–3 hPa (the integration time of the gated medium) before descending.

**5. SAO westerly phases seed QBO westerly onsets.** At 5–7 hPa, QBO westerly onsets cluster in April–June (5 of 21
in May at 5 hPa), just after the SAO's strong April westerly peak above. Easterly onsets cluster less.

**Summary.** The vertical structure matches the switched-plasma picture: a calendar-locked, solar-driven switch at the
top (the SAO, too thin and too strongly damped to integrate anything), its unequal half-cycles supplying the
alternating polarity, and below it a medium with memory where the lunar-modulated beat takes hold, enters at
2–3 hPa and descends.
