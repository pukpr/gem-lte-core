# ERA5 tropical surface pressure: eastward wavenumber 1 and its "13.66-day" sidebands

Replication of the spectrum discussed in
[The just-so story narrative](https://geoenergymath.com/the-just-so-story-narrative/), which comments on
Sakazaki & Hamilton (2020), *An Array of Ringing Global Free Modes Discovered in Tropical Surface Pressure Data*,
JAS 77, 2519 (their Fig. 4). The blog makes four claims:
1. Symmetric double sidebands appear only about the **even** harmonics of the daily tide.
2. They come from the **13.66-day** fortnightly lunar tide (Mf).
3. The odd harmonics are weaker.
4. The pattern comes from an asymmetric, sawtooth-like daily carrier, or a half-day modulation.

**Status:** numbers are from the contiguous 1999–2016 hourly record (18 yr) unless stated; 1979–1986 and 1989–1998
reproduce them (table below). The 1987b–1988 downloads stalled in the CDS queue and were not needed.

## Files

| file | what |
|---|---|
| `fetch_sp.py` | ERA5 hourly `surface_pressure`, 20S–20N, 5° grid, half-year NetCDF files into `data/` (not in git) |
| `wn1_spectrum.py` | k=1 coefficient c(t) and the paper's method (Fig. 4 replica, 1/365 cpd bins, AR(1) background); whole-record DFT line table → `wn1_spectrum.png`, `wn1_lines.json`, `wn1_coef.npy` |
| `sideband_test.py` | joint least-squares fit of S_n (n = 1..11) with ±Mf, ±MSf and annual sidebands; sideband/carrier ratio vs n → `sideband_test.png/.json` |
| `ringing_model.py` | corrected time-domain model of the blog's "atmospheric ringing" figure, simulated and compared with the data → `ringing_model.png/.json` |
| `k0_longperiod.py` | wavenumber-0 band-mean daily pressure: long-period lunar lines, nodal satellites, draconic/tropical controls → `k0_longperiod_<Y0>_<Y1>.json` |
| `antisym_lunar.py` | antisymmetric (0–20N minus 20S–0) k=1 diurnal lunar lines and nodal satellites; k=0 antisymmetric draconic control |
| `l2_check.py` | migrating lunar L2 (westward k=2 M2): amplitude vs classical value, N2/M2, per-UTC-hour window dependence |
| `kmix_test.py` | S2 x M2 products by wavenumber (westward k = 0, 2, 4): physical quadratic mixing vs analysis-window multiplication |
| `strat_tides.py` | same tidal-line analysis for surface pressure and hourly u, v at 70/50/30/10 hPa (partial results below: sp, u30, u50; u70/u10/v downloading) |
| `fetch_uv.py` | ERA5 hourly u, v at 70/50/30/10 hPa, 20S–20N, 5°, one level per request |
| `qbo_sideband_rules.py` | QBO 2 ± q sideband selection by level (KIT radiosondes 70–10 hPa, 1953–2026): seasonal amplitude vs timing modulation, bootstrap |
| `qbo_sideband_descent.py` | why the 2 ± q mix swaps between 30 and 70 hPa: calendar-locked timing gate + descending response (delay from lag correlation) |
| `qbo_gate_model.py` | equinox sample-and-hold of the draconic factor (Meeus ephemeris): does it produce the observed column-wide QBO timing gate? mixed-gating fit + lunar-phase control |
| `qbo_spectral_fit.py` | QBO model (equinox-impulsed draconic → leaky integrator → sin(kL+d) → response) fitted ONLY to the amplitude spectrum; time domain generated for visual comparison |
| `seam_test.py` | jumps by UTC hour and lunar M2 amplitude by UTC hour (assimilation-window test) → `seam_test.png/.json` |

Method: cos(lat)-weighted 20S–20N mean p(t, lon). Its wavenumber-1 coefficient is c(t) = ⟨p e^{-i lon}⟩.
The time FFT of the complex c(t) separates eastward waves (negative numpy frequency) from westward ones; a synthetic test
in the script checks the sign. Exact line frequencies are built from the mean-longitude rates s, h, p.

## Results

**Fig. 4 reproduces.** Averaging yearly spectra gives lines at every integer cpd up to 11. It also shows the broad
free-mode bumps of the paper: the ~33-h Kelvin wave and the ~9.4-h mode.

**Claim 1 (even harmonics only): confirmed.**

**Claim 2 (13.66 days): true in the sidereal frame; the solar-frame spacing is MSf** (see *Sidereal frame* below).
Relative to the solar harmonics S_n, the record resolves Mf (13.661 d) from MSf (14.765 d) by 36 Rayleigh bins:

| eastward, joint LS | n=1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| −MSf amp / noise | 1.1 | 5.3 | 1.1 | **15.9** | 1.7 | **15.7** | 0.7 | **11.6** | 1.0 | **11.2** |
| +MSf amp / noise | 0.7 | 5.5 | 2.0 | 5.8 | 1.5 | 8.4 | 2.2 | **10.9** | 2.4 | **10.9** |
| −Mf amp / noise | 0.7 | 1.3 | 1.2 | 0.9 | 0.9 | 0.8 | 0.5 | 0.9 | 0.2 | 0.9 |
| +Mf amp / noise | 0.4 | 1.2 | 1.2 | 1.5 | 0.8 | 1.5 | 1.4 | 1.6 | 1.1 | 1.7 |

About S_n the pairs sit at ±MSf (the synodic, spring–neap fortnight): S_n − MSf = M2 + S_{n−2}, S_n + MSf = S_{n+2} − M2.
Since MSf = Mf − 2h, the same lines are exact **Mf (13.66 d) offsets of the sidereal/semiannual carriers**:
- lower: S_n − MSf = (S_n + 2h) − Mf, with M2 = K2 − Mf;
- upper: S_n + MSf = (S_n − 2h) + Mf.

Below 1 cpd, **O1 = K1 − Mf**. So "13.66 d" and "14.77 d" are the same lines referred to the sidereal or the solar
day, aliased through the 2h (semiannual) shift of the annual impulse comb. That comb is present: the eastward diurnal
band has P1, S1, K1 of comparable size plus S1 ± 2h.

**Claim 3 (odd harmonics weaker): partly.** True from n = 4 up. There the even carriers are nearly **flat**
(~0.002 hPa, n = 4..10) while the odd ones fall roughly as 1/n².

**Claim 4 (sawtooth daily carrier): not as stated.** A fortnightly modulation of any one daily waveform, sawtooth or
not, gives the same sideband/carrier ratio at every harmonic. Observed: the median ratio is 0.22 at even n and 0.05
(noise) at odd n. A uniform modulation would predict a ~0.006 hPa sideband at n = 1, about 10× what is observed.
A sawtooth also has all harmonics, decaying as 1/n; only a 12-h-periodic factor gives even-only lines.

**But there *is* a 12-h sawtooth, and it is ERA5's assimilation cycle** (`seam_test.png`):
- ERA5 runs 4D-Var in two 12-h windows (09–21 and 21–09 UTC). After removing the mean solar day, the hour-to-hour
  step in c(t) is **5.7–5.8× background at exactly 09→10 and 21→22 UTC**, and ~1.0× at every other hour.
- ERA5's forecast model has no lunar gravitational forcing (S&H 2020 §2a), so M2 enters only through the analysis.
  Fitted separately by UTC hour, the westward k=1 M2 amplitude **collapses at 10 and 22 UTC (0.017 → 0.006 hPa) and
  ramps linearly through each window**: a sawtooth reset at the seams. The eastward M2 jumps at the same hours.
- Multiplying M2 by a 12-h-periodic sawtooth puts lines at M2 + 2j cpd and 2j − M2, i.e. at S_{2k} ± MSf. That gives
  even harmonics only, a 14.77-d spacing, unequal sideband pairs, and a flat tail to n = 10 because the resets are
  sharp. The reconstruction of the line amplitudes from the per-hour fits agrees to <1%, but that agreement is a
  mathematical identity, a consistency check rather than an independent prediction. The evidence is the time-domain
  location of the discontinuities at fixed UTC hours. A physical modulation of a global-scale pressure wave would be
  smooth and could not be locked to 09|10 and 21|22 UTC.
- The flat even solar carriers at n ≥ 4 likely have the same origin. They have not been tested separately yet.

**Other lines.** The eastward diurnal is strongly seasonally modulated: P1 = S1 − h and K1 = S1 + h are comparable to
S1 itself. The paper's AR(1) background fit hits r → 1 here (cosmetic; it only draws the red curve).

## Corrected ringing model (`ringing_model.png`)

This recreates the blog's two-panel figure (forcing in the time domain, model vs data spectrum) with the forcing corrected.
Every part is estimated from the 1999–2016 eastward k=1 data:

| component | blog | corrected |
|---|---|---|
| daily forcing | impulse train | fitted S_n, n = 1..11 with ±1, ±2 cycle/yr neighbours. At 0E it looks like a sharp-rise / slow-decay daily sawtooth |
| lunar | whole waveform × (1 + m cos 2πt/13.66 d) | M2 (east and west) × the measured per-UTC-hour analysis-window function (12-h sawtooth, reset 10 & 22 UTC), plus O1 |
| ringing modes | 33 h, 9.4 h, 7.2 h | noise-driven damped oscillators from Lorentzian fits (S&H Eq. 5): **32.5 h** (half-width 0.018 cpd, e-fold 9 d), **9.3 h** (0.029 cpd, 5 d), **7.2 h** (at the 0.01 cpd minimum width; barely visible at k=1 E) |
| pedestals | — | two "wandering daily pattern" families, each one fixed waveform × a slowly varying random amplitude: even (12-h; harmonics 2..12) and odd (24-h; 3..11) |
| background | — | running median of the data spectrum, random phase |

- **The pedestal families are seam-locked.** Demodulated envelopes are coherent across harmonics: 0.73–0.91 among
  even n = 4..10, 0.56–0.81 among odd n = 3..7, ≤0.16 between the two groups. Each family's waveform is the leading
  eigenvector, reconstructed with no reference to UTC timing, and both change abruptly at 09|10 and 21|22 UTC (figure
  inset). Two alternatives were rejected:
  - a 12-h *multiplicative* modulation of the slow field, which is ~10⁴ too weak;
  - independent per-window ramps, which give broadband noise rather than pedestals.
- **Discriminating score** (rms log₁₀ error of the 40 sideband/carrier ratios, n = 1..10, ±MSf and ±Mf):
  **0.29 corrected vs 0.80 blog-style**. The blog model puts ~0.2 Mf sidebands on every harmonic, odd and even;
  the data and the corrected model do not. A whole-spectrum log correlation does *not* discriminate
  (both ≈ 0.98, dominated by the red background), so it is not used.
- **Triad:** f(32.5 h) + f(9.3 h) = 3.329 cpd vs f(7.2 h) = 3.313 cpd. That is within the 9.3-h mode's half-width,
  consistent with the blog's triad remark but not a demonstration (that needs a bicoherence test).
- **Shortfalls:**
  - Pedestal contrast at even n = 4, 8, 10 is ~1.4 in the model vs 1.6–2.2 in the data.
  - n = 3 is too high (1.28 vs 0.76).
  - The model reproduces about half the eastward seam-step excess (1.4–1.5× vs 3.0–3.1×). The rest is presumably jumps
    that don't repeat from window to window (measured seam-to-seam correlation 0.14); the model carries that energy in its
    stationary background rather than at the seams.

## Independent eras

The same analysis on non-overlapping data (run from scratch copies of the scripts):

| | 1979–1986 | 1989–1998 | 1999–2016 |
|---|---|---|---|
| step excess at 09→10 / 21→22 UTC | 7.1× / 6.6× | 6.0× / 6.6× | 5.7× / 5.8× |
| westward M2 drop at 10 UTC (hPa) | 0.029 → 0.019 | 0.022 → 0.011 | 0.017 → 0.006 |
| ±MSf at even n = 4..10 (× noise) | 3.5–5.7 | 5.3–8.2 | 5.8–15.9 |
| ±Mf, any n (× noise) | ≤2.3 | ≤2.1 | ≤1.7 |

## Sidereal frame, wavenumber 0 and the draconic month

- **Selection rule.** If a 13.66-d factor modulated the whole seasonally pulsed carrier, the largest sideband would sit
  on the largest comb tooth (P1 − Mf). Instead only the tidal-potential lines appear: O1 = K1 − Mf and M2 = K2 − Mf.
  S1 − Mf, P1 − Mf and OO1 are at noise, because the lunar hour angle is the sidereal one minus s.
  Only through the annual comb does the 13.66-d spacing become 14.77 d about the solar harmonics.
- **k = 0 (zonal mean, the QBO's wavenumber)** (`k0_longperiod.py`, 1979–2016, 38 yr): none of Mm, Mf, its nodal
  satellites 2s + n and 2s + 2n, MSf, Mtm, MSm, or the control months (draconic 27.2122 d, tropical 27.3216 d) is
  detectable. All are at 0.5–2.3× noise; noise is 0.01–0.03 hPa in daily band-mean pressure. At k > 0 the
  diurnal/semidiurnal carriers move the lunar content up to 1–2 cpd, where the background is ~10³ lower. That is why it
  is visible only there.
- **Where the draconic month lives** (`antisym_lunar.py`). In the degree-2 potential the k = 0 term goes as δ² and is
  equatorially symmetric, so its draconic content is 2(s + n), 13.606 d. The k = 1 diurnal term goes as sin 2δ and is
  **antisymmetric**, so the declination envelope δ(t), including its inclination (draconic) part, appears linearly only
  there. In the antisymmetric (0–20N minus 20S–0) westward k = 1 band:
  - O1 is 11.0× noise, ~5× its symmetric-band amplitude.
  - Its 18.6-yr nodal satellite **O1′ = O1 − n is 3.1× noise, with O1′/O1 = 0.225 against the equilibrium 0.189**
    (≈ i/ε = 5.15°/23.44°). This is marginal, but both the frequency and the ratio were predicted beforehand.
  - Non-potential controls (S1 − Mf, P1 − Mf, OO1, M1) are ≤1.2×.
  - The k = 0 antisymmetric draconic fundamental (degree 3) is 0.8×.
  - In the symmetric band, the nodal satellites K1′, O1′ and M2′ are all at noise (their expected sizes are ~1×).
- **k = 0 (QBO) forcing: ecliptic only.** In the QBO model the forcing is the longitude-independent (k = 0) part of
  the lunisolar pull, the same at all longitudes at once. It depends only on ecliptic geometry: the Moon's ecliptic
  latitude β = i sin(λ − Ω) (draconic month, 27.2122 d) and the Sun's equinox crossings (semi-annual, as for the SAO).
  The hour-angle (k ≥ 1) tides analysed above are not that forcing. A 400-yr numerical check of the alias:

  | forcing | strongest sub-annual line |
  |---|---|
  | β × semi-annual impulse, same sign each equinox | 1.732 yr |
  | β × semi-annual impulse, **alternating sign** (northward vs southward crossing) | **2.367 yr** |
  | β × annual impulse | 2.367 yr |
  | β² × same-sign semi-annual impulse | 1.183 yr |

  An alternating semi-annual comb has only the odd annual harmonics: 13.4219 − 13 = 0.4219 cpy. So the SAO → QBO link
  needs the equinox impulse to carry the *direction* of the solar crossing (antisymmetric N/S). Its product with β
  (also antisymmetric) is symmetric and k = 0. ERA5 registers the ecliptic (inclination) geometry only marginally,
  through O1′ above. The k = 0 response itself belongs in the stratospheric zonal wind, not in surface pressure, where
  k = 0 lunar lines are below noise. ERA5's own 09/21 UTC analysis updates are also simultaneous at all longitudes
  (k = 0-coherent); they are 12/24-h features, far from the 2.37-yr band.

## Lunar signals in ERA5 are observational

ERA5's forecast model has no lunar gravitational forcing, and the observing system runs on solar/UTC clocks. Lunar
lines at exact astronomical frequencies can therefore enter only through the assimilated observations: they are real
atmospheric signals, not model fabrications. Evidence:
- M2 sits exactly on frequency (0.0-bin offset).
- O1 is at 11× noise, and its 18.6-yr nodal satellite appears near the equilibrium ratio.
- **L2** (`l2_check.py`): **0.068 hPa at the equator** (0.062 hPa averaged over 20S–20N), close to the classical
  barometric ~0.07 hPa.
- **N2/M2 = 0.156** against the equilibrium 0.191. N2 needs the lunar perigee (anomalistic month), which nothing in
  ERA5 knows about.

Two qualifications:
1. **Source vs spectral placement.** The even-harmonic sideband comb (MS4, S_n − MSf, …) is real M2 multiplied by the
   12-h analysis-window function. Its lunar content is genuine, but those lines are not atmospheric tides.
2. **Amplitude is window-dependent.** L2 fitted per UTC hour ranges over 0.053–0.080 hPa (max/min 1.5, phase spread
   ~40°), jumping at the 09|10 and 21|22 seams. Long-record averages are close to the classical value; instantaneous
   lunar amplitudes carry about ±20% analysis-cycle error.

## Are the 14.77-day (MSf) sidebands nonlinear mixing? (`kmix_test.py`)

By frequency they are the classical combination tides of S2 and M2 (MSf = S2 − M2, MS4, 2SM2, …). Wavenumber separates
the mechanisms. Quadratic mixing of the migrating S2 and M2 (both westward k = 2) puts the sum MS4 at **k = 4** and the
difference MSf at **k = 0**. A UTC-locked analysis window multiplies all longitudes alike and keeps the wavenumber.

| westward, sym, 1999–2016 (hPa, ×noise) | k = 0 | k = 2 | k = 4 |
|---|---|---|---|
| S2 | 0.074 (36) | 1.231 (169) | 0.021 (31) |
| M2 | 0.006 (20) | 0.062 (53) | 0.0019 (8) |
| MS4 = S2 + M2 | 0.0008 (10) | 0.0097 (83) | 0.0024 (16) |
| M4 = 2 M2 | — | 0.0002 (2) | 0.0001 (0.9) |
| MSf = S2 − M2 | 0.004 (0.2) | — | — |

- **k = 2 (and the eastward k = 1 comb): analysis window.** The wavenumber is preserved; MS4/M2 = 0.16.
- **k = 4: a k-changing process.** MS4/M2 = 1.25 rules out window leakage of the k = 4 M2, and the degree-2 lunar
  potential has no k = 4 component. The ratios fit quadratic S2 × M2 mixing: MS4/S4 = 0.07 vs 2·M2/S2 = 0.10, and M4 is
  at noise as (M2/S2)² predicts. This is suggestive, not proven: S4 is also directly heated, and wavenumber mixing by
  uneven observing networks is argued against, not excluded.
- **The difference branch** (MSf at k = 0, the zonal-mean forcing relevant to QBO/ENSO-type behaviour) should be of
  comparable size, ~0.002 hPa. That is about 10× below the k = 0 noise at fortnightly periods in surface pressure, so it
  is untestable here rather than absent. The stratospheric winds are the test (watch for M2 aliasing to MSf through
  00/12 UTC radiosonde sampling).
- **To annual scales:** under an annual impulse MSf and Mf alias to the same fractional frequency (they differ by
  2 cpy): 24.7366 and 26.7366 cpy → 0.7366 / 0.2634 cpy, i.e. **3.80 yr** and 1.36 yr. The QBO's 2.37 yr instead needs
  the draconic fundamental with a signed (odd-harmonic) impulse, as above.

## QBO sideband selection rules by altitude (`qbo_sideband_rules.py`)

The requirement for an astronomical-tide model of the QBO: explain why 2 + q survives at some levels and 2 − q at others
(q = 0.4219 cpy, the draconic alias). Data: KIT/FU Berlin monthly radiosonde winds, 7 levels, 1953–2026, 74 yr.

Write the QBO as Re[a(month) e^{i2πqt}], with a(month) fitted per calendar month over all years, and a = ā(1 + α + iφ):
α is the seasonal **amplitude** modulation and φ the seasonal **timing** modulation. Then the upper sideband 2 + q goes as
|α₂ − iφ₂| and the lower 2 − q as |α₂ + iφ₂|. One sideband cancels when the two semiannual modulations have equal depth
and are ±90° apart, and the sign picks which survives. This is an algebraic **restatement of the requirement**, not
independent evidence: the agreement between the phase sign and the dominant sideband is automatic.

| hPa | 70 | 50 | 40 | 30 | 20 | 15 | 10 |
|---|---|---|---|---|---|---|---|
| upper/lower (2+q : 2−q) | 0.82 | 0.55 | 0.36 | **9.1** | 0.64 | 0.81 | **10.6** |
| 95% bootstrap | 0.26–3.5 | 0.13–2.6 | 0.09–1.4 | 0.73–8.9 | 0.19–2.0 | 0.24–4.0 | 0.88–14.9 |
| AM − PM phase (± sd) | −143 ± 84 | −48 ± 60 | −61 ± 36 | +98 ± 56 | −135 ± 56 | −22 ± 68 | +81 ± 54 |

- **The asymmetry is not significant at any single level** (every 95% interval includes 1). The 30-vs-70 hPa flip is
  suggestive only. The QBO's wandering period broadens all its lines.
- **The informative part is the altitude structure.** The timing modulation is nearly level-independent: the QBO phase
  is most retarded around April and advanced July–November at every level, like one seasonal clock for the column. The
  seasonal **amplitude** peak descends: ~Feb at 20 hPa, May–Jun at 30, Jun–Jul at 40, Jul at 50, Aug at 70 hPa. Their
  relative semiannual phase therefore winds monotonically with depth: unwrapped from 10 hPa, ≈ −500° by 40 hPa,
  ≈ −50°/km, with per-level sd 36–68° below the 70 hPa end. The surviving sideband flips whenever the winding crosses
  quadrature.
- **Model requirement in physical terms:** a seasonally gated timing modulation acting on all levels at once (as an
  equinox-comb sample-and-hold would), plus an amplitude modulation that propagates downward with the QBO. The sideband
  pattern at each level then follows from where the descending amplitude peak sits relative to the fixed timing clock.
- **Not yet reconciled:** ERA5 u70 (1979–2020) gives a much stronger lower sideband (power ratio 57.6 vs 8.8) than KIT
  70 hPa (amplitude ratio 1.2), though the direction agrees.

## Reconciling the 30 → 70 hPa sideband swap (`qbo_sideband_descent.py`)

The objection runs: an astronomical forcing has fixed spectral lines, so the sideband mix should change only
infinitesimally over ~5 km. But the mix is not a property of the forcing alone. It is the interference of a
**calendar-locked** timing gate (measured: the same at every level, most retarded around April) with a **descending**
response. A pure delay τ cannot change |2 ± q|. A calendar-locked term plus a delayed term rotate against each other by
720°·τ (semiannual), however. The descent from 30 to 70 hPa takes ~8 months, which rotates the semiannual relative
phase by more than a full turn, while the 2.37-yr carrier shifts by only ~0.28 cycle. The sidebands act as a sensitive
interferometer on the descent delay, so large changes over a few km are expected, not anomalous.

Test with no fitted constant. The delay τ is taken from lag correlation with 10 hPa (independent of the seasonal
analysis), and the predicted winding is −720°·τ:

| hPa | 10 | 15 | 20 | 30 | 40 | 50 | 70 |
|---|---|---|---|---|---|---|---|
| delay (mo) | 0 | 2 | 4 | 6 | 9 | 11 | 14 |
| winding predicted / measured (sd) | 0 / 0 | −120 / −103 (68) | −240 / −216 (56) | −360 / −343 (56) | −540 / −502 (36) | −660 / −489 (60) | −840 / −584 (84) |
| lower/upper observed | 0.09 | 1.23 | 1.57 | **0.11** | 2.80 | 1.81 | **1.22** |
| model A (measured AM depth, phase descended from 10 hPa) | 0.06 | 1.26 | 1.63 | **0.09** | 5.07 | 0.71 | **1.26** |
| model B (10 hPa AM, delayed only) | 0.06 | 1.55 | 1.94 | **0.06** | 16.7 | 0.65 | **1.55** |

- From 10 to 40 hPa the measured winding matches −720°·τ within 1σ at every level. A high r alone would be weak, since
  both quantities grow monotonically with depth; the slope near 1 with no free constant is the informative part.
- Model B (one common gate plus the 10 hPa AM delayed by τ) gets the dominant sideband right at 6 of 7 levels,
  including the 30 hPa (upper) → 70 hPa (lower) swap. Model A also matches the magnitudes at 10, 15, 20, 30 and 70 hPa.
- **Misfit below 40 hPa:** the winding stalls (50 and 70 hPa ~3σ short of −720°·τ) and model B misses 50 hPa.
  Plausibly the weak lower-stratosphere QBO (|ā| = 4 m/s at 70 hPa) picks up an AM that is partly calendar-locked
  (annual upwelling) rather than descended. This is untested.
- **Scope:** this removes the objection; it does not by itself discriminate astronomical forcing from internal
  seasonal synchronization. Any column-wide calendar gate plus a descending response gives the same interference. What
  the astronomical model must supply is that gate: its calendar phase (April timing extreme) and its depth (~0.1 rad
  semiannual).

## Does an equinox sample-and-hold produce the gate? (`qbo_gate_model.py`)

Forcing: the Moon's ecliptic latitude β = sin F (Meeus mean argument of latitude, real draconic phase), sampled at the
actual equinox instants (Meeus; March → September 186.4 d) and held, run raw or through a 60-d response, averaged
monthly over 1953–2026, and analysed exactly like the radiosondes. Target: the observed column-mean timing gate, most
retarded in **March** (–14° at month 3, +10° in July), semiannual depth 0.10 rad, annual 0.15 rad.

| comb (60-d response) | most retarded month | semiannual / annual timing depth (rad) | r with observed curve |
|---|---|---|---|
| observed | **3** | 0.10 / 0.15 | — |
| annual, March (northward crossing) only | **3** | 0.21 / 0.85 | **0.82** |
| unsigned semiannual | 4 | 0.20 / 2.33 | 0.72 |
| signed semiannual (+Mar, −Sep) | 9 | 0.21 / 0.51 | −0.41 |
| annual, September only | 9 | 0.21 / 0.85 | −0.44 |

- **Calendar phase: matched.** Pure sample-and-hold gives a timing sawtooth that resets in the month of the
  sampling equinox. A comb dominated by the March (S→N) crossing reproduces the observed March retardation;
  September-dominated combs, including the signed comb in this sign convention, anti-correlate.
- **Depth: pure sample-and-hold is 2–15× too strong.** The QBO is only partly gated. A two-parameter mix fitted to the
  12-month curve (comb asymmetry ρ = b/a, gated fraction g, the rest an ungated q line) picks **ρ ≈ 0 (March-only),
  g ≈ 0.15** and gives rms 4.2° against 7.3° with no gate (observational scatter ~2.7°/month). With a 90-d response,
  g ≈ 0.20 and rms 3.2°.
- **Not reproduced: the semiannual gate** (model 0.03 rad vs observed 0.10). No ρ reaches the observed
  semiannual/annual ratio (0.68; best pure comb 0.41). Yet this semiannual part is what drives the 2 ± q sideband
  selection. It needs another ingredient: a semiannual response (e.g. coupling to the SAO), a second gate, or a
  nonlinear hold.
- **Control: the gate is a property of the comb, not of the Moon.** Shifting the draconic phase by arbitrary offsets
  leaves the gate unchanged (r = 1.000). So the calendar match supports the **equinox sample-and-hold** part; the
  **lunar** part is tested by the carrier frequency (the observed 74-yr peak at 2.36 yr vs the 2.370 yr draconic alias),
  not by the gate.

## Spectrum-only fit, time domain for visual comparison (`qbo_spectral_fit.py`)

The only fitting step matches the model's **amplitude spectrum** to the observed one (KIT 30 hPa, 1953–2026,
deseasonalised, 0.1–3 cpy). The time series is then generated and drawn against the observations, with no time-domain
fitting. u and −u have the same spectrum, so both are drawn. Model: draconic factor sin F (Meeus) impulsed at the real
equinoxes (weight 1 March, ρ September) → leaky integrator (T_i) → A sin(kL + d) → one-pole response (τ).

| run (figure) | criterion | ρ | T_i (yr) | k | τ (mo) | main peak (obs 0.420 cpy) | misfit |
|---|---|---|---|---|---|---|---|
| `qbo_spectral_fit_30_log.png` | log amplitude, all bins equal | −0.50 | 2.2 | 52.8 | 4.6 | 0.468 | 0.230 dex |
| `qbo_spectral_fit_30_lin.png` | linear (energy-weighted) | −0.76 | 5.3 | 52.9 | 4.9 | 0.439 | 0.446 |
| `qbo_spectral_fit_30_lin_s1` (seed 1) | linear | −0.33 | 12.5 | 49.4 | 3.7 | 0.408 | 0.427 |
| `qbo_spectral_fit_30_lin_k5.png` | linear, k ≤ 5 | −0.78 | 0.9 | 2.5 | 4.2 | **0.422** | 0.605 |

What the visual comparison shows:
- **Log criterion: unusable.** Equal weighting of all bins fits the broadband floor and misses the 2.37-yr peak.
- **Low winding (k ≤ 5), the near-pure sample-and-hold regime:**
  - The spectrum is a rigid comb of lines at q, 1 ± q, 2 ± q, 2q, 3q. Several coincide with observed bumps (e.g. 2 + q),
    but the observed peak is broad.
  - In time the model is strictly periodic at 2.370 yr. It sits in phase with the observations over some stretches,
    then drifts through quadrature to antiphase (the −model tracks other decades). **The mean period is right, but the
    QBO's phase wanders by up to half a cycle over decades, and a fixed astronomical clock cannot follow it.**
  - Waveform: sharp rise then slow decay (the leaky hold), against observed sharp transitions both ways with flat
    plateaus.
- **High winding (k ≈ 50):** the large phase jumps at each impulse broaden the peak and give QBO-like cycles with
  wandering phase. But the alignment with the observations is no better, semiannual ripples appear that the
  observations lack, and the parameters are **not identified**: different seeds give ρ −0.33…−0.76 and T_i 5…12.5 yr
  at equal misfit, and seed 2 found the seed-7 solution with d + 180° (the sign ambiguity).
- **Consistent across all energy-weighted fits: ρ < 0.** The September impulse enters with the opposite sign to March,
  i.e. an alternating-sign comb (N→S vs S→N crossing).

Mismatches to target: (1) phase wandering. Real QBO cycles last 20–36 months, while the forcing clock is rigid. That is
the phase-slip / dynamic-winding element of the production model, which this standalone version lacks. (2) The
square-wave shape with plateaus, which needs a saturating response rather than a leaky hold.

**Saturating response added** (5th argument `sat`: y = tanh(G(sin(kL+d) + b))/tanh(G) before the low-pass;
same spectrum-only fitting):

| run (figure) | ρ | T_i (yr) | k | τ (mo) | G | b | main peak | misfit |
|---|---|---|---|---|---|---|---|---|
| `qbo_spectral_fit_30_lin_k5_sat.png` (k ≤ 5) | −0.66 | 4.9 | 3.48 | 2.9 | 20 (bound) | +0.21 | **0.422** | 0.445 (was 0.605) |
| `qbo_spectral_fit_30_lin_sat.png` (k ≤ 60) | −0.33 | 15.5 | 49.8 | 4.3 | 20 (bound) | −0.29 | 0.419 | 0.376 |
| `qbo_spectral_fit_30_lin_s1_sat` (seed 1, k ≤ 60) | −0.66 | 4.9 | 3.48 | 2.9 | 20 | −0.21 | 0.422 | 0.445 (= k ≤ 5 solution, sign-flipped) |

- **The fit wants a hard switch** (G at the bound), i.e. a sign comparator: square wave with flat plateaus.
- **Low winding + saturation (`_k5_sat`), the clearest view:**
  - Shape is now QBO-like: plateaus, sharp switches, and stronger easterlies than westerlies via b (≈ +19 / −25 m/s),
    as observed.
  - Remaining mismatches: (i) the rigid 2.370-yr clock still drifts in and out of phase over decades; (ii) every switch
    is abrupt-then-exponential, the same in both directions, whereas at 30 hPa the switch into easterlies is visibly
    slower than into westerlies; (iii) the spectrum has a strong **2 − q** line (1.578 cpy) that the 30 hPa
    observations lack (the selection-rule problem, now in the model), and a deep notch near 2 cpy the observations
    lack.
  - Reproducible: seed 1 with k free reached the same solution.
- **High winding + saturation (`_sat`):** best spectrum (broadened peak), but the time series shows phase slips at the
  wrong times. One example is a 2000–2004 easterly stall the observations do not have, plus spurious short excursions.
  The flexibility is spent on spectral breadth, not on timing.
- ρ < 0 again in all runs.

Next targets in order: timing (phase wandering / slips at the right times); direction-dependent transition speed
(hysteresis rather than a symmetric switch); the 2 − q suppression at 30 hPa.

## Tides from the surface into the QBO levels (`strat_tides.py`, partial run, 2026-10-27)

This run covers surface pressure and hourly u at 30 and 50 hPa, 1999–2016, 20S–20N. For u, 32 of 36 half-years are in; 2003 and 2012 are still downloading. Output is `strat_tides_1999_2016.json`. Each value below is amplitude/noise (SNR) in the cos(lat)-weighted equatorially symmetric component.

| line | sp | u30 | u50 |
|---|---|---|---|
| k = 0 Mf / Mm / MSf | 0.9 / 1.4 / 0.2 | 0.8 / 1.0 / 0.3 | 0.2 / 0.3 / 0.4 |
| k = 0 Mtm / MSm | 1.0 / 0.5 | 0.9 / 1.3 | 1.8 / 0.7 |
| westward k = 2 S2 / M2 | 487 / 53 | 141 / 51 | 166 / 35 |
| westward k = 2 M2/S2 amplitude | 0.050 | 0.087 | 0.074 |
| eastward k = 2 S2+MSf / S4+MSf | 11.7 / 20.1 | 4.6 / 2.5 | 3.5 / 1.6 |
| seam step excess (09\|10, 21\|22 UTC) | 4.0× / 4.3× | 5.1× / 5.0× | 6.5× / 6.2× |

- **No zonal-mean lunar long-period line at the QBO levels.** Mf, Mm and MSf at k = 0 are at noise (SNR ≤ 1.4) at the surface, 30 and 50 hPa. The draconic and tropical controls are also at noise. So hourly ERA5 does not carry the k = 0 Mf forcing as a direct spectral line. Mtm at 50 hPa (1.8, and 2.0 antisymmetric) is marginal, given the many candidate lines.
- **The migrating lunar semidiurnal tide does reach the stratosphere.** Westward k = 2 M2 has SNR 51 at 30 hPa and 35 at 50 hPa. Its ratio to S2 rises from 0.05 at the surface to 0.07–0.09 aloft. IFS has no lunar forcing, so this M2 comes from the assimilated observations (see "Lunar signals in ERA5 are observational").
- **The ±MSf sidebands weaken with height.** The eastward k = 2 S2+MSf / S4+MSf family that marks the 4D-Var seams in surface pressure (12–20) drops to 2–5 at 30/50 hPa. The seam step excess is larger aloft (5–6.5×), so these sidebands should still be treated as assimilation artifacts.
- **Diurnal lunar O1 is weak everywhere** (SNR ≤ 2 in u). K1 and S1 are strong, as expected for solar thermal tides.

Rerun when u at 70/10 hPa and v at all levels are in: `python3 strat_tides.py 1999 2016`.

## Caveats and follow-ups

- SNRs and amplitude/noise ratios are not corrected for searching several candidate frequencies. The ±MSf lines at
  even n ≥ 4 (11–16× noise) are unaffected; values of 5–8 are marginal.
- Controls not yet run:
  - a reanalysis with a different window layout, e.g. MERRA-2 (6-h windows), where the comb should move to multiples
    of 4 cpd;
  - raw hourly station barometers, where any lunar sidebands would be physical;
  - ERA5 split into its separately produced streams, to see whether seam strength changes with the observing system.

Fixed 2026-10-26: the lunar perigee rate was entered as 1/8847 d⁻¹ instead of 1/3232.6 d⁻¹ (8.85 yr). That moved Mm,
Q1, N2 and L2 by ~1–2 Rayleigh bins (N2 SNR 6.3 → 7.0, Q1 4.8 → 3.4 after the fix); Mf, MSf, M2, O1 and K1 do not
depend on it.
