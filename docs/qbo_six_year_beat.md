# The CW's 6.4-yr / 6.0-yr beat in the QBO (2026-10-26)

Code: `qbo_six_year_beat.py` (figure `qbo_six_year_beat.png`, numbers `qbo_six_year_beat.json`). CW inputs, read and not
refitted: `ChandlerWobble/docs/cw_suppression.json` (two-line beat minima) and `cw_longrecord.json` (IERS C01
envelope, 12-yr windows). QBO data: KIT radiosondes, 1953–2026, 20–70 hPa.

## 1. Why the beat carries over

The CW alias is exactly twice the QBO alias: 2F − 26 = 2(F − 13), so F_CW = 2q, with q = 0.4219 /yr. In polar motion
(`cw_suppression.md` §6.1), two 6-year signals meet:
- **annual − CW = 1 − 2q = 0.1560 /yr (6.41 yr);**
- **F − M = 0.1668 /yr (6.00 yr),** the lunar perigee–node cycle (draconic − anomalistic).

Folded by the annual (odd-harmonic) impulse, F − M sits at 1 − (F − M) = 0.8332 /yr beside the CW at 0.8438. The two
lines beat with a ~93-yr envelope: the fitted minima are 1924–29 and 2017–22, and out of sample (1890–1972) r = 0.42–0.49.

F − M is linear in the ecliptic latitude, the QBO's kind of drive, and the QBO comb is odd-harmonic. So the QBO has the
same pair twice, each separated by 2q + (F − M) − 1 = 0.0107 /yr (93 yr):

| QBO band | QBO line | lunar F − M line | F − M beat in the self-referenced phase |
|---|---|---|---|
| 0.84 /yr | 2q = 0.8438 (2nd harmonic, = F_CW) | 1 − (F − M) = 0.8332 | **−39° per decade** |
| 0.58 /yr | 1 − q = 0.5781 (annual sideband) | q + (F − M) = 0.5887 | **+39° per decade** |

(Also q − (F − M) = M − 13 = 0.2551 /yr, the anomalistic alias.)

**What a 74-yr record can show.** It cannot resolve either pair, since that needs 93 yr. It can show the beat as a slow
drift of each band's phase **relative to the QBO line itself**: phase(2q) − 2·phase(q), and phase(1 − q) + phase(q).
These are constant if the band is only the QBO's own waveform or seasonal modulation, and cancel the QBO's period
wandering exactly. The record can also show at most one **node**: an amplitude minimum with a phase jump of up to 180°.

## 2. Results (12-yr windows, 1-yr step; fixed-line least squares, as in `cw_longrecord.py`)

**Band 0.84 (2q), where the CW's own pair lives:**

| relative-phase drift (° per decade) | 20 hPa | 30 hPa | 50 hPa | 70 hPa |
|---|---|---|---|---|
| 1966–2002 (F − M beat: −39) | −79 ± 43 | **−41 ± 16** | **−38 ± 11** | **−46 ± 12** |
| 1963–2020 | −86 | −85 | −78 | −84 |

- **1960s to about 2003: drift at the F − M beat rate.** At 30–70 hPa it matches within 1σ, with the right sign.
- **Then a node-like event.** The 2q amplitude reaches its minimum (50 hPa: 0.20 m/s in the window centred 2004), the
  relative phase steps by about −200° between the 2002 and 2008 windows, and the amplitude then recovers (2q/q:
  0.05 → 0.32 at 30 hPa by the 2014 window, before the disruptions enter). That is the sequence a two-line beat
  produces at a node: decline, phase jump, recovery. The amplitude also declines gradually toward it, from about
  0.2–0.4 of the QBO line in the late 1970s to about 0.05–0.2 by 2002–08.
- **The symmetric-waveform alternative fails.** A square wave's 2nd harmonic vanishes when the westerly and easterly
  phases last equally long, but the westerly duty cycle stays at 0.58–0.63 (50 hPa) and about 0.66 (70 hPa) through
  2004, and does not track the 2q ratio (r = −0.04 to +0.14).

**Band 0.58 (1 − q): no F − M signature.** The relative-phase drift is −90 to −250° per decade at 20, 30 and 70 hPa,
the wrong sign for the beat (+39). At 50 hPa it is +11 ± 18. This band is dominated by the drift of the QBO's seasonal
synchronization (`era5/qbo_sideband_rules.py`).

**Against the CW:** the QBO's off-clock bands strengthened relative to the QBO line from the 2000s, and sharply once the
2016 and 2020 disruptions enter the windows, while the CW envelope declined. The descriptive correlations of the band
envelopes with the CW envelope are −0.14 to −0.64, with only about 5 independent windows.

## 3. Reading

- **In the 0.84 band, the QBO behaves like the CW's two-line beat.** It drifts at the F − M rate for four decades and
  then passes a node. That band is the QBO's version of the Chandler frequency (2q = F_CW).
- **The QBO node falls around 2004–06, not 2017–22.** The CW's own envelope has its first step down in 2003–07
  (`cw_qbo_timeline.md`), and the QBO's phase advance also begins in 2005–07. The CW's fitted node (2017–22) depends on
  its resonance: the Q = 32–85 response shifts the relative phase of two lines only 0.011 /yr apart, and the node year
  moves by 5 yr across Q. The QBO has no sharp resonance at 0.84 /yr, so its node is set closer to the forcing's own
  phases. An offset between the two node years is therefore expected; its size has not been predicted.
- **The 1950s start does not allow more.** With 74 yr against a 93-yr beat, the QBO can show at most one node, and the
  one it shows (about 2004) is it. The previous node would be about 1911 and the next about 2097, both outside any
  equatorial wind record. The 1920s CW minimum and reversal have no QBO counterpart to check.

## Caveats

- The 1966–2002 slope window was chosen after seeing the 2004–08 step (post hoc). Over the whole record the slope is
  about twice the beat rate.
- The 2F + M line (17.6-yr beat with the CW) has minima at 2004/06 too, and 12-yr windows cannot resolve a 17.6-yr
  modulation. The node could belong to either beat.
- Other waveform changes than the duty cycle (e.g. transition sharpness) can change the 2nd harmonic and were not tested.
- Window uncertainties assume one independent sample per 12 yr. The 2008–2020 windows include the disruptions.
