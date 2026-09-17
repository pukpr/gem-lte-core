# Methods blurb — winding-number scalogram (draft for public release)

Draft status: text only; references verified 2026-09-13 against Fyfe et al.
(order tracking), Antoni (cyclostationarity reviews), and the lock-in /
tidal-harmonic lineage. Aim: ~350 words, suitable for a preprint Methods
section or software-paper summary. A condensed 1-paragraph version follows
the full text.

---

## Winding-number scalograms: phase-locked spectral analysis on a
## non-monotone manifold

Consider a scalar climate record x(t) and a known latent phase coordinate
F(t) — the tidal forcing manifold of the Laplace's Tidal Equation (LTE)
framework, calibrated independently against the measured length-of-day
(dLOD) record (|r| > 0.995 retained after simultaneous multi-index
fine-tuning). A dynamical component "winds" around this manifold at
integer or fractional wave-number M if it contributes
A sin(2*pi*M*F(t) + phi) to the record. Because F(t) is non-monotone — it
reverses direction hundreds of times per century of monthly sampling and
is uncorrelated with any linear clock (linear-fit R^2 = 0.02) — calendar
time and winding "time" are incommensurate clocks, and no reparameterized
wavelet or Fourier transform on t can substitute for the winding axis.

We therefore compute the winding scalogram

    G(t0, M) = < w(t - t0) x(t) exp(-i 2*pi*M*F(t)) >,

a Gaussian-windowed correlation of the record against the manifold's own
native basis, localized in time t0. This is the fixed-manifold,
time-localized analogue of two established techniques. In rotating-machinery
diagnostics, computed order tracking (Fyfe & Munill 1977; Fyfe, Fritzen &
Rand 1997) resolves vibration spectra against shaft angle rather than time,
removing speed-fluctuation smearing; classical implementations resample the
signal at constant angular increments, which presupposes a monotone phase.
When the reference phase is non-monotone and multiply periodic — as for a
tidal manifold built from Doodson arguments (s, h, p, N) — resampling is
undefined and direct projection onto exp(-i 2*pi*M*F) is required, which is
what G(t0, M) performs. In experimental physics, lock-in (phase-sensitive)
demodulation measures amplitude and phase at a known drive phase; winding
scalograms generalize this in two directions: to a continuum of
simultaneous "orders" M, and to time-localized windows, yielding ridge
structures — persistent, sharp ridges in (t0, M) — as the signature of
genuine phase locking rather than incidental correlation. The closest
climate-domain ancestor is classical tidal harmonic analysis (Doodson 1928),
which regresses against known astronomical arguments but globally in time;
the scalogram's windowed form makes non-stationarity of the locked components
visible, and the surrogately-normalized power (log2 |G|^2 divided by an
AR(1)-matched noise floor through the identical basis) calibrates
significance against persistence, not white noise.

Applied to spatially distributed records, the procedure inverts the
conventional dimension-reduction order. EOF/PCA analysis diagonalizes the
temporal covariance and reports whatever directions maximize variance —
an assumption-free-in-form but assumption-laden-in-fact basis, since
empirical modes of non-stationary or phase-modulated oscillations split and
mix across several unrelated EOFs. A systematic sweep over sampling boxes
(lat, lon) — e.g. gridded Kaplan SST V2 or NCEP Reanalysis extracts — and
over M instead yields a continuum: for each geographic point, a ranked set
of winding numbers with amplitude, phase, and persistence statistics,
defined relative to a single externally calibrated clock. Structure is read
off directly — dominant winding versus multi-mode spread, inter-index
frequency locks, near-degenerate beat pairs — as coordinates with physical
meaning fixed a priori, rather than variance-ordered mixtures re-derived
separately for every dataset.

## One-paragraph condensed version

Winding-number scalograms measure how strongly a climate record contains
components phase-locked to a known, externally calibrated manifold F(t):
G(t0,M) = <w(t-t0) x(t) exp(-i 2 pi M F(t))>, a windowed correlation against
the manifold's native basis. Mathematically the technique is a fixed-manifold
generalization of computed order tracking (Fyfe et al. 1977, 1997) and lock-in
demodulation to a non-monotone, multiply periodic phase reference — the
turbine-shaft-angle resampling of classical order tracking is undefined for
such a clock, so direct projection is used — and of Doodson tidal harmonic
analysis, made time-local and surrogate-calibrated (AR(1) noise floor through
the identical basis). Sharp, persistent (t0,M) ridges then separate
non-autonomous, forcing-locked dynamics from autonomous (internally generated)
variability — a discrimination that variance-maximizing decompositions such as
EOF/PCA cannot make, since their basis is re-derived per dataset and has no
reference clock to be locked or un-locked to.

## References to cite (verified)

- Doodson, A.T. (1921). "The harmonic development of the tide-generating
  potential." Proc. R. Soc. Lond. A 99, 570-577. [the 38 constituents and
  five fundamental arguments behind every Doodson number]
- Doodson, A.T. (1928). "The reduction of observations of monthly mean sea
  level." (published for the International Association of Geodesy
  sea-level commission; see also Doodson & Warburg 1941 for the handbook
  version.) [the tide-to-MSL "reduction" direction the baltic result
  inverts]
- Doodson, A.T. (1928). "The tidal topic." Proc. R. Soc. Lond. A 121,
  571-602. [phase-argument regression framework]
- PSMSL (Permanent Service for Mean Sea Level), founded 1928; baltic
  composite = PSMSL ID 10006, "Baltic 211 sites composite" (monthly MSL,
  1880-2025). [data lineage of the baltic index]
- Fyfe, K.R. & Munill, P. (1977). "Computed Order Tracking Obsoletes Older
  Methods." Sound and Vibration (reprinted SAE). [origin of computed order
  tracking: resample at constant shaft-angle increments]
- Fyfe, K.R., Fritzen, C.-P., Rand, A.G. (1997). "Analysis of computed order
  tracking." Mechanical Systems and Signal Processing 11(2), 187-205.
- Antoni, J. (2009). "Cyclostationarity by examples." Mechanical Systems and
  Signal Processing 23(4), 987-1036. [frame: order analysis = spectral
  analysis against known phase, within cyclostationary theory]
- Meade, M.L. (1983). "Lock-in amplifiers: principles and practice."
  IEEE Trans. Education 26(4), 173-175. (or Butterworth & Harmon, Rev. Sci.
  Instrum. 49, 1978) [phase-sensitive demodulation at a known reference phase]
- Doodson, A.T. (1928). "The tidal topic." Proc. R. Soc. Lond. A 121,
  571-602. [the direct ancestor: regression on known astronomical phase
  arguments via Doodson numbers]
- (project) LTE tidal manifold & dLOD calibration: gem-lte-primitives.adb,
  src/gem-dlod.adb; WINDING_NARRATIVE.md.

## Lineage and ancestry: the Doodson arguments, and the inverted operation

The winding-number scalogram sits on a pedigree older than any of the
signal-processing ancestors above. In 1921 Doodson expanded the
tide-generating potential into 38 primary constituents whose phases are
linear combinations of five fundamental astronomical arguments — Moon's
mean longitude (s), Sun's (h), lunar perigee (p), lunar node (N), and
solar perigee (p1) — codified as the Doodson numbers that still name every
tidal constituent today (Mm, Mf, M2, S2, N2, K1…). In 1928 Doodson
published the *reduction method*: tables allowing the tides to be removed
from monthly mean sea level, developed for the international geodetic
sea-level program that founded PSMSL the same year — a program built on
precisely the kind of multi-gauge monthly MSL composites analyzed here.
For a century, the discipline's use of these arguments has run in one
direction: tidal terms are *nuisance, subtracted from MSL* so that
"climate" (or geodesy, or land uplift) can be studied in what remains.

The baltic index turns the operation around, and that inversion is the
historically intuitive headline. `baltic` is a PSMSL composite of 211
coastal gauge sites (ID code 10006); averaging the basin suppresses
site-specific, instrumental, and short-period noise, and what survives —
per the scalogram survey above — is one dominant, sharp, whole-record ridge
at M = 1.2454 = **6.00 x 0.2076**, the exact order-6 harmonic of the same
base winding frequency independently converged on by six other basins
(6.7 bits over the AR(1) null, FWHM 0.06, continuity 1.00). That is:
non-tidal sea level of the Baltic, after everything classical reduction
would strip away, is *built from* the Doodson manifold — windings at
integer harmonics of a frequency calibrated against measured dLOD. The
predictive claim follows in the generative direction: an MSL time series
synthesized from the calibrated forcing alone (LTE Model vs baltic data
CC = 0.70 on the survey CV window — the ridge finding is solid; the
full-fidelity generation claim should be stated at this strength, not
beyond it).

"Conceivable in 1925" is therefore nearly true, and the near-miss is
instructive: Doodson had the arguments, the Baltic gauge networks, the
composite means. What the era lacked was (i) the measured dLOD record to
calibrate the manifold against a *known answer* rather than assume it,
(ii) the non-autonomous framing — classical harmonic analysis presumes
stationary constants over the analysis window, and 18.6-yr/8.85-yr
precession terms enter as constituents in their own right, whereas here
they emerge as beating *modulations* of the fast Mm/Mf windings against
the annual cycle, and (iii) the regression machinery to search a
continuum of winding numbers rather than a fixed constituent library.
Item (iii) was computationally hard in 1925; items (i) and (ii) were
conceptional, which is why the inversion took a century.

## Notes for Paul — positioning the EOF/PCA contrast honestly

The claim above is defensible as written; two cautions to keep it so:

1. Don't say EOFs are "blind" or "wrong" — say they answer a different
   question. EOF/PCA asks "which directions maximize variance?" Winding
   scalograms ask "which components are phase-locked to a calibrated
   external clock, where and when?" A non-stationary locked oscillation
   genuinely does smear across several EOFs (well documented — e.g. ENSO
   modes distributed across EOF pairs), so the two answers legitimately
   differ; the winding answer has the advantage that its coordinates are
   fixed a priori and comparable across datasets and sites without
   re-derivation. That last property is what makes the geospatial winding
   continuum (your M(lat,lon) field) possible in a way an EOF map is not:
   EOF labels only mean something within the dataset that produced them.

2. Anticipate the obvious reviewer question: "the manifold must be right."
   Answer is already in hand — dLOD calibration (independent validation
   target with a known astronomical answer), the fixed cross-index
   manifold (r = 0.9995 between indices), the IAAFT surrogate rejection
   (no ridge for arbitrary series), and the wavenumber-0 QBO contrast.
   State those as the method's falsification battery; the surrogate
   rejection is the strongest single card because it shows the instrument
   can answer "no."
