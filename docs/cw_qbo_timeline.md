# The Chandler wobble collapse and the QBO anomalies: one timeline, two devices (2026-10-26)

Code: `cw_qbo_timeline.py` (figure `cw_qbo_timeline.png`, numbers `cw_qbo_timeline.json`). Inputs, nothing refitted:
- **CW:** `ChandlerWobble/docs/cw_suppression.json`. This holds the envelope and phase of raw IERS polar motion
  demodulated at 0.8440 /yr (Gaussian σ = 1 yr), and the required forcing per epoch.
- **QBO:** the KIT radiosonde winds, demodulated at 0.4219 /yr with the same smoothing, and the required pulse train of
  the inverted QBO device (`qbo_control.json`).

The two devices share the semiannual comb and the draconic month (`cw_suppression.md` §5):

| | drive | comb | alias | comb harmonic |
|---|---|---|---|---|
| CW | B² (sign-blind) | same polarity | 2 × 13.4219 − 26 = 0.8438 /yr | 26 |
| QBO | B (sign-sensitive) | alternating | 13.4219 − 13 = 0.4219 /yr | 13 |

## 1. When: the observed timeline

Phases are relative to each alias clock (QBO relative to 1973–90; 30 hPa shown, 50 and 70 hPa agree within ~20°).

| | Chandler wobble | QBO |
|---|---|---|
| to 2002 | 160–210 mas (peak 211 in 1993) | phase within ±50° of the 2.370-yr clock |
| **2003–07** | first step down, 159 → 104 mas | phase advances: −2° (2005) → +48° (2007) |
| 2008–11 | plateau 105–120 mas | phase plateau about **+60°** at all levels |
| **2010–15** | **collapse 112 → 16 mas** (half of the 2010 level by about 2012.3; e-fold 2.6 yr; minimum 2015) | phase eases back, +62° → +29°; nothing abrupt |
| **2016.1–2016.2** | phase swing begins (−2° in 2015 → +35° in 2016 → +58° in 2018) | **first disruption**: 40 hPa easterly inside a 30/50 hPa westerly (Feb–Mar 2016); phase **jumps to about −100°** (−30° / −96° / −114° at 30 / 50 / 70 hPa) |
| 2016–19 | anomalous phase, amplitude 20–27 mas | anomalous phase, −80° to −113° |
| **2020.0** | phase returns (+29° → −3°); **recovery begins** (21 → 27 → 36 → 47 mas by 2022) | **second disruption** (Jan 2020); phase jumps back to about +43° |
| 2021–25 | (record in the json ends 2022) | phase advances steadily, +52° → +161°: cycles of about 2.0 yr |

- **The two QBO disruptions bracket the CW's anomalous-phase interval.** The QBO jumps out of phase in Feb–Mar 2016
  and back in Jan 2020. The CW's phase swing starts in 2015–16 and ends in 2020, when its recovery begins.
- **Both records have an earlier, smaller step around 2003–07:** the first CW drop, and the QBO's +50–60° phase
  advance.
- **The CW amplitude collapse (2010–15) has no abrupt QBO counterpart.** In those years the QBO phase only drifts.
- **Smoothing caveat:** σ = 1 yr blurs dates by about ±1 yr. The QBO "amplitude" dips in 2016 and 2020 come mainly from
  averaging across abrupt phase jumps, not from a weaker QBO. The CW phase is poorly defined near its 2015 minimum.
  The disruption months themselves come from the monthly 40 hPa winds, not the smoothed demodulation.

## 2. How: each device's required forcing

Relative to 1973–90 (amplitude / phase):

| epoch | QBO device (pulses) | CW, Q = 85 | CW, Q = 32 |
|---|---|---|---|
| 1990–2005 | 1.06 / −39° (30) | 0.49 / −108° | 0.46 / −42° |
| 2005–10 | 1.26 / **+50°** (10) | 0.65 / −121° | 0.45 / −90° |
| 2010–13 | 1.41 / +29° (6) | **2.31 / +122°** | 1.04 / +111° |
| 2013–16 | 1.06 / +76° (6) | **0.98 / +128°** | 0.51 / +123° |
| 2016–19 | 1.55 / **−122°** (6) | 0.36 / −52° (quiet) | 0.18 / −34° (quiet) |
| 2019–27 | 0.75 / +90° (16) | 0.87 / +1° (resumed) | 0.56 / 0° |

- **CW device:** a counter-forcing episode. The forcing reverses in 2010–16, which free decay cannot mimic (it would
  need Q ≈ 7). It nearly vanishes in 2016–19 and resumes near its original phase after 2019.
- **QBO device:** phase-locked to the astronomical comb until about 2005. It then needs a forcing rotated about +50°
  (2005–10), stays moderately rotated through 2016, **flips to −122° in 2016–19**, and settles at about +90° after
  2019, with a noisier input (`qbo_control.md` §3). At the pulse level, the 2016 disruption is one out-of-schedule pulse
  (2016.15, +2.7σ).
- **The devices' anomalies are offset:** the CW forcing reverses in 2010–16 and goes quiet in 2016–19, which is exactly
  when the QBO forcing flips.

## 3. One cause? The comb-timing test

A shift Δt in the seasonal comb timing rotates the CW alias forcing by 360·26·Δt and the QBO's by 360·13·Δt. So the
QBO rotation should be half the CW one, known only modulo 180° because the CW rotation is known modulo 360°.

| epoch | 1990–2005 | 2005–10 | 2010–13 | 2013–16 | 2016–19 | 2019–27 |
|---|---|---|---|---|---|---|
| miss, Q = 85 | 15° | 70° | 32° | 12° | 83° | 90° |
| miss, Q = 32 | 18° | 85° | 27° | 14° | 75° | 90° |

- **2010–16 is consistent:** the CW reversal (+111° to +128°) halves to the QBO's +29° to +76° within 9–32°. A week-scale
  comb shift can account for that episode in both devices at once, as `cw_suppression.md` §5 proposed (7 days: CW 179°,
  QBO 90°).
- **The other epochs are not:** 2005–10, 2016–19 and 2019–27 miss by 65–90°.
- **Overall, chance level.** With the 180° ambiguity, a random phase lands within 30° of a prediction with probability
  ≈ 1/3, so about 2 of 6 epochs are expected to "hit"; 3 are observed. A single, persistent comb-timing change is
  therefore **not** supported as the common cause of the whole decade. At most it fits the 2010–16 episode.

## 4. Reading

- **The coincidence is real and specific.** The QBO's two disruptions (Feb–Mar 2016, Jan 2020) bound the CW's
  anomalous-phase interval and the start of its recovery. Both records also show a smaller step in 2003–07. This
  agrees with the published attribution of the CW's lost atmospheric excitation to the QBO disruptions (through the
  meridional angular-momentum budget; cited in `cw_suppression.md`).
- **The two models register the decade differently.** The CW device needs a reversal and then silence. The QBO device
  needs a drift, a flip and a reset. They are offset by one epoch: the CW reverses in 2010–16; the QBO flips in
  2016–19. A shared astronomical comb that simply re-timed itself does not produce that pattern.
- **A shared atmospheric event fits better than a shared astronomical one.** The 2016 disruption was driven by
  extratropical Rossby-wave momentum intruding into the equatorial stratosphere, and that same atmosphere is the CW's
  excitation source. In control terms (`qbo_control.md`), both devices received an out-of-schedule input from the
  atmosphere in this decade. It is not a change in the lunar or seasonal schedule.

## Caveats

- The epoch tables are descriptive: each epoch has its own complex forcing value. The QBO's short epochs (6 pulses)
  give noisy phases.
- The QBO's 1973–2005 agreement with the comb is in-sample with respect to its plant fit.
- The CW json ends in 2022. Extending the CW envelope to 2026 (`fetch_eop.py` in the ChandlerWobble repo) would show
  whether its recovery tracks the QBO's post-2020 phase advance.
