# Winding Scalogram — Qualifications & 17-Index Comparative Ridge Survey

Companion to `WINDING_NARRATIVE.md` and `winding_scalogram.py`. All numbers
below were produced by `/tmp/winding_check/ridge_survey.py` (2026-09-13
session; candidate for promotion into this directory), reading each index's
own `lte_results.csv` (Data = col 3, the non-circular panel),
`dlod_compare.csv`, and `compute_winding()` from `param_survey.py`.

## Qualifications brought forward (now verified, not assumed)

1. **The manifold is external-calibrated before it is ever climate-fit.**
   Pearson r between the two columns of each index's `dlod_compare.csv`
   (measured dLOD vs model forcing): **+0.9946 to +0.9989** for all 14
   winding-number indices (nino4 .9951, pdo .9960, amo .9963, tna .9964,
   nao .9947, baltic .9964, iode .9951, emi .9951, nino34 .9989,
   darwin1880 .9980, npi .9951, pna .9946, noi .9948, tpi .9958).
   Column 4 of every `lte_results.csv` is therefore not tuned per-index to
   that index's climate data — the CC>0.99 tie to LOD survives the
   simultaneous multi-index fine-tune. (Caveat: qbo30/qbo50's
   `dlod_compare.csv` pairs correlate at r≈−0.02 — that file compares a
   different quantity under their wavenumber-0 manifold; worth re-checking
   what its columns are before citing it.)
2. **The manifold is fixed and shared.** nino4 vs pdo column-4 correlation
   on the overlapping record: **0.9995**. It is also genuinely a
   *non-monotone* clock (162 dF sign changes; linear-in-time R²=0.023), so
   the transform is not secretly a calendar-wavelet relabeling — calendar
   oscillators land at M≈0.03, never at the winding ridges (verified with a
   matched-apparent-period control).
3. **Circularity is now properly scoped.** The Model panel remains
   self-confirming by construction (it is regressed onto this basis). The
   qualification does not change that — but because the manifold is
   externally calibrated (point 1) and held fixed, **Data-panel ridges are
   genuine evidence of phase-locked (non-autonomous) structure**, and the
   Data-vs-Model ridge comparison measures how much of that locked power
   the fit actually captured. That is the right way to present the tool.
4. **Categorical control for wavenumber-0.** qbo30/qbo50 were fit against a
   wavenumber-0 manifold (their forcing spans only ±2.4, 168–180 sign
   changes, |slope| 0.001–0.003/yr — a different, near-zero-drift clock).
   Their scalograms validate the method's discrimination from the opposite
   side (see survey below: qbo30 still shows a persistent, sharp-ish ridge
   — 7.4 bits peak, 100% window continuity).
5. **Negative control.** `pdo_iaaft_detuned_surrogate` — an arbitrary
   IAAFT-drawn time series through the identical pipeline — shows **no
   horizontal ridge anywhere**: its best-M structure has FWHM 0.93 in M
   (vs 0.04–0.08 for real windings) and its argmax wanders across windows.
   Visible in `pdo_iaaft_detuned_surrogate/winding_scalogram.png`.

## Survey table (Data panel; white-noise floor; σ=15yr windows)

pkSNR = max log2(power/noise-floor) on the ridge; FWHM = width of the mean
|G| profile at half-max around Mdom; cont = fraction of time windows whose
local ridge stays on Mdom above 2 bits. AR1 = matched pipeline run on AR(1)
ρ=0.97 red noise, at *its own* most-favorable M (an upper bound for the
null).

    index                           n  dLOD_CC   M(NM)    Mdom  pkSNR   FWHM  cont |  AR1fwhm AR1cont
    nino4                        1731  +0.9951  0.2073  0.4491    6.2  0.050  1.00 |    0.090    0.96
    pdo                          1706  +0.9960  0.2073  0.4488    6.4  0.070  1.00 |    0.090    0.92
    amo                          1712  +0.9963  0.2075  0.0134    8.1  0.040  1.00 |    0.150    0.73
    tna                          1741  +0.9964  0.2075  0.2075    6.4  0.690  0.58 |    0.140    0.54
    nao                          1722  +0.9947  0.2077  0.8307    5.8  0.040  0.85 |    0.140    1.00
    baltic                       1742  +0.9964  0.2076  1.2454    6.7  0.060  1.00 |    0.050    0.74
    iode                         1705  +0.9951  0.2073  0.4455    4.9  0.470  0.69 |    0.150    0.88
    emi                          1667  +0.9951  0.2073  0.4488    5.8  0.300  0.76 |    0.090    1.00
    nino34                       1732  +0.9989  0.2079  0.8166    7.2  0.040  1.00 |    0.140    0.85
    darwin1880                   1731  +0.9980  0.4377  8.0838    4.4  0.080  1.00 |    0.240    0.88
    npi                          1530  +0.9951  0.2073  4.1457    4.0  0.060  0.91 |    0.070    0.78
    pna                           917  +0.9946  0.2076  1.4532    4.6  0.050  1.00 |    0.120    1.00
    noi                           942  +0.9948  0.2081  3.7453    4.6  0.040  1.00 |    0.060    1.00
    tpi                          1758  +0.9958  0.2072  0.4473    6.5  0.060  1.00 |    0.100    0.74
    qbo30                         856  -0.0210  0.1105  0.1105    7.4  0.250  1.00 |    0.400    0.58
    qbo50                         856  -0.0205  0.1115  0.1115    5.3  0.260  0.92 |    0.420    0.58
    pdo_iaaft_detuned_surrogate  1706  +0.9964  0.2076  0.9247    5.3  0.930  0.77 |    0.140    1.00

## Automated tooling delivered alongside this report (2026-09-13)

- **`winding_rank.py`** — goal (2): ranks winding numbers objectively from
  any index subdir (`--index nino4`) or a bare .dat remapped onto a chosen
  manifold (`--dat X.dat --manifold-from nino4`). Null model upgraded to
  AR(1) ρ=0.97 persistence surrogates (fixes the white-noise
  under-calibration above); a ridge must be strong (≥2 bits), sharp
  (FWHM ≤0.12) and steady (continuity ≥0.70). Validated: nino4 → [3.11,
  0.45]; amo --m-max 9 → [3.11, 0.02, 1.04];
  `pdo_iaaft_detuned_surrogate` → **NONE** (clean rejection).
- **`psl_extract.py`** — goal (3): programmatic NOAA PSL timeseries form
  (Kaplan SST V2 `--ntype 4`, NCEP Reanalysis `--ntype 1`, etc.); builds a
  ready analysis subdir (.dat + year×month .txt). Verified to reproduce the
  saved `kap10-10-20-30` extract exactly (corr 1.000000, max diff 0.0000).
  The month convention is **start-of-month** (year + m/12), and PSL reports
  the *cell-centered* grid actually used — e.g. requesting 20..30N returns
  the 22.5..32.5 box (nearest cells) — so compare grids, not just inputs.
  End-to-end demo: fresh `atlbox5n15s` (tropical E. Atlantic box, Kaplan)
  → ranked against the nino4 manifold → single passing ridge at M=0.01
  (2.4 bits, FWHM 0.04, continuity 0.85) — AMO-like low-frequency locking;
  record is short (1948–2023 → only 13 windows), so treat as a pilot.

## Reading

- **11 of 14 winding indices show a textbook ridge**: FWHM 0.04–0.08 and
  cont ≥ 0.85 — a sharp, horizontally persistent ridge at (or adjacent to)
  their fitted Mdom. The AR(1) null at its own best M never combines sharp
  width with high SNR *away from M≈0*: red-noise ridges pile up near
  M≈0.02–0.05, whereas these ridges sit at 0.449, 0.82, 0.83, 1.25, 1.45,
  3.75, 4.15. AMO at 0.0134 is the deliberate exception — its claim rests
  on FWHM 0.040 + cont 1.00 + 8.1 bits (the sharpest, strongest ridge in
  the survey), which the AR1 null at low M does not reproduce (fwhm 0.15,
  cont 0.73).
- **Weak/borderline cases**: tna (FWHM 0.69, cont 0.58), iode (0.47/0.69),
  emi (0.30/0.76). These three need the red-noise-surrogate upgrade before
  their ridges can be asserted — at those widths they are barely better
  than the null's wandering low-M ridge.
- **Surrogate behaves exactly as advertised**: no horizontal ridge
  (FWHM 0.93 ≈ 93 grid cells); its 5.3-bit "peak" is broadband red-like
  power with a wandering argmax. Clean rejection.
- **QBO wavenumber-0 pair**: on their own near-zero-drift manifold,
  both still show continuous ridges at their fitted base M (qbo30: 7.4 bits,
  cont 1.00; AR1 control on that manifold is *also* fuzzy — fwhm 0.40,
  cont 0.58), i.e. the method is not merely detecting the LOD manifold's
  spectral texture; it discriminates signal from null on a categorically
  different clock too.
- **darwin1880 needs a scalogram re-run**: its Mdom = 8.0838 = order-19
  harmonic (NM=1, harm=[4]; base M(NM)=0.4377) sits ABOVE the plot's
  default `--m-max 5.0`, so its persistent ridge (cont 1.00, FWHM 0.080)
  is off the top of the current PNG. Regenerate:
  `python3 winding_scalogram.py darwin1880 --m-max 9`. Also flag: its
  M(NM)=0.4377 departs from the ~0.207 backbone — worth understanding
  before citing it alongside the other thirteen.

## Revised viability verdict

With the external dLOD calibration, the near-identity of column 4 across
indices, the wavenumber-0 QBO contrast, and the surrogate rejection, the
tool's central claim holds up as executed: **localized winding-power against
a fixed, externally calibrated, non-monotone manifold separates
non-autonomous from autonomous content, and recovers ground-truth M values
on synthetics**. Remaining caveats that survive the qualifications:
(i) swap the white-noise floor for AR(1)/persistence surrogates per index
    before publishing tna/iode/emi ridge claims (the 11 sharp ones are safe
    as-is; AMO's low-M claim is defensible on width+continuity but the
    surrogate upgrade makes it airtight);
(ii) significance should be quoted as ridge width + temporal continuity,
    not peak height alone;
(iii) present Model-panel panels as fit-adequacy diagnostics, not evidence.

Prior-art framing unchanged: position against order tracking (Fyfe et al.)
and lock-in demodulation; the novelty is the application — a time-localized
winding-number scalogram against a non-monotone tidal phase manifold — not
the transform's skeleton.
