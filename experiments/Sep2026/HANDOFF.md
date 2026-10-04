# GEM-LTE grid-cell sweep — handoff notes (2026-09-27, updated 2026-09-29)

## NEW: MANIFOLD quadrature-flip regularizer added at the Ada source (2026-09-29)

Direct answer to the subharmonic quadrature metastability finding
documented below ("Subharmonic quadrature metastability -- a real,
unresolved structural degeneracy"): the user's own proposed fix --
"regularize the metric via a correlation coefficient that compares to a
mean or median manifold ... If CC~1.0, no penalty, but if CC~0, large
penalty" -- is now implemented in `GEM.LTE.Primitives.Solution`'s local
`Metric` function (`src/gem-lte-primitives-solution.adb`), the single
choke point every accept/reject decision in the search already routes
through (CorrCoeff, CorrCoeffP, CorrCoeffTest, Exclude_Metric,
Validate_Metric all call it).

**New env var `MANIFOLD`** (default unset = feature off, zero behavior
change): points at a reference manifold `.dat` file (same two-column
`date value` format as any other `.dat` file this project reads, on the
SAME date grid as the cell's own data -- meant to be a mean/median
manifold across the 89-quad collective, written by a sweep script, never
by `lt.exe` itself). When set and the file exists, every `Metric(X,Y,Z)`
call now also computes `Manifold_CC := CC(Z, Manifold_Ref(Z'First ..
Z'Last))` (Z is always some slice of Forcing/the manifold already, so no
separate date-alignment step is needed) and multiplies the raw score by
`Long_Float'Max(0.0, Manifold_CC)` -- CC~1.0 (this run's own manifold
slice agrees with the collective reference) leaves the score untouched;
CC~0.0 (the classic quadrature-swap signature -- sin/cos orthogonality,
not a sign flip) or negative (anti-correlated) crushes it toward/through
zero, steering the search away from the wrong quadrature without a
second, separate search objective.

**Safety note (important if reusing this pattern elsewhere)**:
`Make_Data`'s own failure path for a missing/bad file is
`GNAT.OS_Lib.OS_Exit(0)` -- silent process termination, NOT a catchable
exception -- so `MANIFOLD` is deliberately checked with
`Ada.Directories.Exists` BEFORE ever calling `Make_Data`, rather than
wrapping the call in a handler. Confirmed live: a `MANIFOLD` pointing at
a nonexistent file produces byte-identical output to not setting it at
all, no crash.

**Verified via a new standalone self-test**
(`src/manifold_regularizer_test.adb`, built via `gprbuild -P lte.gpr
manifold_regularizer_test`, same pattern as the existing
`iir_invariant_test.adb`) rather than the full threaded search harness
(background/non-PTY runs of the real binary hit an unrelated,
pre-existing environment flake -- `STORAGE_ERROR` from GNAT's interrupt
manager, reproduced identically with AND without `MANIFOLD` set,
confirming it's not caused by this change, just not a usable test
signal here). The self-test confirms all 5 properties the real penalty
logic depends on: self-vs-self CC=+1.0 (no penalty), self-vs-negated
CC=-1.0 (multiplier clamps to 0.0), both clamp values exactly, and --
the one genuinely subtle point -- that a misaligned/too-short reference
correctly raises (and is caught as) `Constraint_Error` rather than
propagating uncaught. That last check caught a real Ada semantics trap
during development: an exception raised while ELABORATING a block's own
declarative part is NOT caught by that same block's own handler in Ada,
only by an ENCLOSING one -- the production code already had the correct
double-nested `begin declare ... begin ... end; exception ... end;`
shape, but the test's first draft didn't and crashed uncaught until
fixed to match.

**Not yet done**: no sweep script writes a `MANIFOLD` file yet (per the
user's own framing, "the contents of the MANIFOLD file will be set by
the sweep script and stored at this level" -- i.e. one shared file at
the Sep2026 experiment root, likely the mean or median of all 89
cells' own current manifolds). No cell's `.resp` sets `MANIFOLD` yet
either, so this is purely a capability landing, not yet exercised
against a real mean-manifold file or a real ambiguous cell.

**First live test on kN040_W050 (2026-09-29 evening, resumed after a power
outage)**: `medianManifold.dat` (built by `build_manifold_reference.py`)
flags kN040_W050 (CC=-0.04) and kS020_W010 (-0.30) as the only
off-quadrature cells. Ran all 4 cascade attempts on a copy
(`_manifold_test_kN040_W050/`, TIMEOUT=400, each from the SAME fresh
seed = kN040_W050's accepted lt.exe.p, per-attempt outputs in
`attempts/`). Backbone stays 0.2294 in all four. Results:

| attempt | manifold CC vs median | raw fit CC all / excl 2000-05 | dLOD | gate |
|---|---|---|---|---|
| original (no MANIFOLD) | -0.041 | 0.775 / 0.785 | 0.9991 | -- |
| VALIDATE=T DTW | +0.663 | 0.684 / 0.691 | 0.947 | fail (dLOD) |
| VALIDATE=F DTW | +0.230 | 0.737 / 0.762 | 0.9989 | pass |
| VALIDATE=T CC  | +0.372 | 0.767 / 0.791 | 0.992 | fail (dLOD) |
| VALIDATE=F CC  | +0.290 | 0.748 / 0.776 | 0.987 | fail (dLOD) |

The pre-outage run (attempt 3 chained on attempts 1-2, so ~3x400s of
search) reached manifold CC +0.80 with raw fit 0.788/0.804, BETTER than
the original on both counts -- so the regularizer does rotate the
quadrature, but 400s from a wrong-quadrature seed isn't enough. Fresh
400s runs are part-way there and currently cost some fit/dLOD.

**Gotcha found**: the penalty multiplies EVERY Metric call, including the
short 2000-2005 test/validate window reports. Over that 5yr slice the
manifold-vs-median CC is often <=0 (e.g. -0.10, -0.74), so the PRINTED
test/val1 scores come out exactly 0.0 / -0.0 -- they're penalized scores,
not comparable to unregularized ledger scores (and they can trip
`stiff()`/acceptance logic spuriously). Probably the penalty should apply
only to the search objective (full/training Z), or use the full-record
manifold CC, and the reported scores should stay raw. Not changed yet.

**trend/accel bookkeeping fixed in the Ada (2026-09-29, uncommitted).**
Two defects made lt.exe.windings.json (and stdout :trend:/:accel:)
misdescribe the model's secular term; the MLR itself was always fine:
1. Accel reference mismatch. Regression_Factors fits accel*(t-ref)**2
   with ref = first date of its fit array (1950 under EXCLUDE, where
   TRAIN_START..TRAIN_END=2000..2005 is the excluded TEST gap), but LTE
   evaluated around Forcing(First).Date = 2000 (adadb00's fix assumed
   TRAIN_START = training start). Fix: Regression_Factors now returns
   `Accel_Ref` (out param; gem-dlod.adb caller updated) and LTE uses it.
   Verified: identical trend/accel, column 2 changes by exactly
   accel*((t-1950)**2-(t-2000)**2) (IR delay applied), CC 1.0000.
2. SECULAR (default TRUE) took abs() of trend/accel AFTER the model was
   built -- only the saved/printed values changed sign. Removed.
   kS020_W010's MLR accel is -5.6e-6 (saved as +5.6e-6 before).
windings.json now also has `secular_context`: accel_ref, fit_span,
train_start_end + role, record span, rate at start/end, rise,
lowest_winding. All 89 quads were re-evaluated in place with the fixed
lt.exe (TEST_ONLY, own lte_run.sh settings, lte_run.sh/resp restored;
pre-fix copies in the session scratchpad only). Every JSON now has
secular_context, accel_ref=1950 everywhere; 42/89 MLR accels and 9/89
trends are negative (all had been saved positive). No dLOD fell below
0.994, none deadlocked. lt.exe.p changed only in ann1/ann2/sem1/sem2
(write-only copies of the MLR's annual cycle, refreshed). 27 quads'
lte_results.csv changed by more than the tilt; an old-binary re-run of
each showed data/manifold columns identical to the new ones, i.e. those
CSVs had been stale relative to lt.exe.p, and the model difference is
the predicted tilt to within 6%. Prorated result (global_sst_from_quads.py):
trend +0.084 C/decade, accel -1.16e-6 C/yr^2, rate +0.084 -> +0.082
C/decade, secular rise +0.605 C; stage 2 model CC 0.950 raw / 0.855
detrended. Printed CC of
existing parameters drops very slightly under the corrected evaluation
(kS020_W010 train 0.6895->0.6884, kN040_W050 0.7840->0.7838).
Curvature budget (stage 3, corrected JSON): model total 1.19e-4
C/yr^2 = data 1.19e-4; accel term -1%, lowest winding 31%, other
windings 29%, k0*F 25% -- the curvature is carried by the sinusoids and
the manifold term, not by accel.

**ACCEL opt-out + overnight Feb2026 sweep (2026-09-30 ~01:00).**
Python re-fit of all 89 Sep2026 quads (own manifold + windings, with vs
without the (t-1950)**2 column): the other regressors explain a median
97% of that column; training CC gain +0.0002; decade-blocked CV WORSE in
63/89; 2000-05 holdout 46/89 better (coin toss); prorated global CC
0.9471 without vs 0.9468 with. So `accel` adds nothing and destabilizes
the trend/accel split. New env/resp var ACCEL (in gem.adb Options enum;
default TRUE = byte-identical output, verified) drops the column in
Regression_Factors (accel := 0, trend = plain linear rate); JSON
secular_context records accel_enabled.
Overnight: experiments/Feb2026_sweep/ = copy of the 89 Feb2026 quads +
frozen lt.exe (00:49 build) + dlod3.dat; run_overnight_sweep.py runs
each quad (6 parallel): resp standardized to EXCLUDE TRUE /
TRAIN_START 2000 / TRAIN_END 2005 / ACCEL FALSE, TEST_ONLY baseline
(DTW, VALIDATE=TRUE), sweep.py 4-attempt cascade (TIMEOUT 300),
TEST_ONLY re-score, keep only if cascade gates pass and min(train,
validate,test) >= baseline. Log sweep.log, per-quad results.jsonl,
SWEEP_DONE + global_sst_from_quads.txt/.png at the end. Resumable: just
re-launch `nohup setsid python3 run_overnight_sweep.py >> sweep.log 2>&1 &`
(quads already in results.jsonl are skipped; an interrupted quad may
need its folder re-copied from Feb2026 first). Feb2026 is untouched.

**AERO volcanic-aerosol regression term (2026-09-30, uncommitted).**
Per the CSALT post (geoenergymath.com/2014/03/12/csalt-volcanic-aerosols):
AERO=<file> (two-column "date value") adds one MLR column
Aero_Value(t) = series(t - AERO_LAG) and LTE adds aero*Aero_Value(t).
Unset = off; verified byte-identical model output vs the previous build.
Files in experiments/ext_data/: aero_saod.dat (GISS CMIP7 2.2.1 global
cos-lat mean stratospheric AOD, 1850-2023) and aero_sparse.dat (VEI>=5
sparse profile, the post's intensities /100, fast-rise/1-yr-decay pulse
peaking ~0.3 yr after each eruption); CC 0.875 with each other 1950-2023.
Both are aerosol LOADING; AERO_LAG (years) adds the temperature response
delay (post: cooling minimum 8-16 months after eruption -> ~0.5-1.0).
JSON: top-level "aero", secular_context aero_enabled/aero_file/aero_lag;
stdout ":aero:". Coefficient < 0 = cooling. First check, kN000_W130 with
its current params and 2000-2005 gap, lag 0.5: aero +1.25 (SAOD) /
+1.18 (sparse), i.e. WARMING -- the El Nino coincidence again; not yet
tried with the 1980-1985 gap or inside a search.
Also fixed in the same edit: after a search, trend/accel/aero/M/MAP
reported and saved were the LAST evaluated candidate's (possibly
rejected), not the kept one's; now stored at accept and restored after
the loop (Keep_Trend etc.). All JSON files produced this week came from
single TEST_ONLY evaluations, so they were not affected.

**Sep2026 -> Feb2026_sweep same-metric merge (2026-09-30 afternoon).**
All 89 Sep2026 parameter sets (.p + .resp) scored exactly like the
overnight sweep (TEST_ONLY, DTW, VALIDATE=TRUE, EXCLUDE 2000-2005,
ACCEL=FALSE; 3 controls reproduced the overnight scores exactly). Sep2026
beat the sweep on 25 quads by min(train,validate,test); the 5 hand-tuned
equatorial quads (AERO/NONLIN) were not touched. Installing all 24 lowered
the global detrended CC (0.8710 -> 0.8698), so only the 8 that raise the
global fit on their own AND keep manifold CC vs median >= 0.5 were kept:
kN000_E150 kN020_W090 kN060_E170 kN060_W050 kS020_E010 kS040_E050
kS040_E070 kS040_W130 (global detrended 0.8710 -> 0.8754 vs prorated
data; 0.856 -> 0.860 vs Kaplan global). kN040_W050's Sep set helps
globally but is the known off-quadrature one (-0.04), excluded. Backups
and scoring in Feb2026_sweep/_sepcmp/ (final.json lists kept/reverted).

**Collective R^2 on Feb2026_sweep + quadrature fixes (2026-09-30).**
mlr_shared_windings.py / collective_r2_log.py / sweep_manifold.py now take
--root DIR (default unchanged: Sep2026 reproduces 0.9464 exactly);
sweep_manifold.py also takes --env KEY=VALUE and --donor QUAD.
Collective R^2: Sep2026 0.9464, Feb2026 0.9425, Feb2026_sweep 0.9373 --
two inherited off-quadrature quads (kN060_W170 -0.245, kS020_W090 -0.269
manifold CC vs median). kN060_W170 fixed via neighbour seed kN060_E170 +
MANIFOLD (CC +0.996, fit 0.643 -> 0.644, global SST det CC 0.8754 ->
0.8751): collective 0.9382, per-quad min 0.502 -> 0.853. kS020_W090 NOT
fixed: seeds kS020_W110 (90 s and 300 s) and kS020_W070 (300 s) all reach
quadrature (+0.96) but fit only 0.59-0.64 vs its current 0.721; still
off-quadrature. Backups in Feb2026_sweep/_qfix_backup/.

**1880-1949 back-cast of the Feb2026_sweep quads (2026-09-30/10-01).**
Each quad's Kaplan box series rebuilt 1880-2023 (matches the existing
.dat exactly from 1950); TEST_ONLY with EXCLUDE=FALSE, fit 1950-2023
only, model evaluated from 1880 (Feb2026_sweep/_backcast/, plot
backcast_1880.png). Anchoring: non-extended quads (IDATE 1880, data from
1950) were seeded at the first sample, so the back-cast sets IDATE 1949.99
(IIR seeds Jan 1950, exact backward pass to 1880); IDATE >= 1950 kept;
INIT_DATE/STRICT_IDATE quads unchanged. NOT INIT_DATE=1950: since adadb00
("Option A") init always means the value at IDATE in strict mode; the
comment block near solution.adb ~707-755 still describes the old meaning
(stale). Result: all 89 manifolds reproduce 1950-2023 (81 exact, 8 to
1e-5). Composite vs Kaplan global 1880-1949: raw CC +0.445 is mostly the
backward-extrapolated linear trend term (+0.076 C/decade gives -0.75 C by
1880 vs data ~-0.2); without it raw CC +0.489 and trends match (+0.021 vs
+0.026 C/decade). Band CC 1880-1949 (fit period in brackets): >20 yr +0.77
(0.99), 7-20 yr +0.17 (0.91), 2-7 yr +0.08 (0.92), 1-2 yr +0.01 (0.73).
Multidecadal shape transfers; interannual/ENSO timing does not.
ALSO FOUND: lt.exe TEST_ONLY finishes in ~5 s but never exits (shutdown
hang); every "2-minute" TEST_ONLY this week was the watcher's own timeout
killing a finished process. Workaround in the back-cast runner: wait for
the final dLOD lines + settled lte_results.csv, then kill the process
group. Not fixed in Ada.

**Alternate long-period quadrature sweep launched (2026-10-01 01:00).**
experiments/Oct2026_altquad/run_altquad_sweep.py. Seeds: the full-record
Sep2026 index fits amo, nino4, nino34, nino3, nino12 -- one shared manifold
(pairwise CC +1.000), CC -0.28 vs Sep2026 medianManifold, backbone 0.1736
(Feb2026's amo/nino4/nino34 are in the CURRENT quadrature, +0.95).
All 89 quads freshly seeded (index .p + .resp, STRICT_IDATE TRUE so the
1880-anchored init carries over -- verified: seeded quad reproduces the
index manifold, CC 0.999998). MANIFOLD lock: altMedianManifold_seed.dat
(z-median of the 5 index manifolds) for the 20 North Atlantic quads
first, then altMedianManifold.dat rebuilt from accepted NA quads + indices
for the rest (outward by distance). Per quad: 5 seeds x 300 s CC search,
unpenalized TEST_ONLY DTW re-score; accept best min score with dLOD >=
0.994, backbone in range, manifold CC >= 0.5 vs reference (else best kept
+ flagged). Results: results.jsonl, sweep.log, SWEEP_DONE,
global_sst_from_quads.txt. Note: mlr_shared_windings' k-set assumes the
0.207 backbone, so the collective R^2 does not apply to this set as is.

**FIXED 2026-10-01: saved lt.exe.p did not reproduce the displayed/kept
optimum (harmonics).** Root cause: harmonic multipliers live in the local
Harms/Harms_Keep arrays; D.C is only read at start-up, so DKeep := D at
accept carried the ORIGINAL file's harmonics, and the 09-29 save-time
rebuild of M from DKeep.C (commit 5d66099) wrote those back as "harm".
Any search that accepted a harmonic change saved a file whose reload gave
a different regression (same manifold, same B params, different harm) and
a lower score than displayed. Reproduced on kN000_W130 (reload model off
by 0.2 C, 0.7938 displayed vs 0.7918 reloaded); unrelated to AERO/NONLIN.
Fix: at accept, D.C(1..NH) := Harms before DKeep := D. Also the
non-VALIDATE final report now uses the kept candidate (D := DKeep; Model
:= KeepModel), as the VALIDATE branch already did -- the final score line
was computed from the last-tried candidate. Verified: 8 replicates
(3 with changed harmonics) reload to 1e-11; VALIDATE=FALSE displayed ==
re-run exactly. With VALIDATE=TRUE the middle (validate) number can still
differ slightly: it is the live running maximum across threads by design.
Consequence for earlier sweeps (Feb2026_sweep overnight, Oct2026_altquad):
their files are self-consistent (every kept result was re-scored with
TEST_ONLY from the saved file) but harmonic improvements found during
those searches were dropped at save. DEBUG_KEEP=TRUE (env) prints
full-precision fingerprints at each accept and at save.
Oct2026_altquad/lt.exe replaced with the fixed build (old one kept as
lt.exe.sweep_build_0930).

**CC refinement sweep of Oct2026_altquad (re)launched 2026-10-01 11:23.**
Sep2026/sweep_cc_refine.py --root ../Oct2026_altquad --ref
altMedianManifold.dat --ALL --timeout 250 --threads 4 --workers 4 (~1.8 h).
Recipe = the user's post-save-fix run on kS040_W010 (fit CC 0.55 -> 0.72,
manifold still locked): own current params, METRIC=CC, VALIDATE=TRUE, F9=1.
Accept: train not lower, validate not lower by >0.01, 2000-05 test not
lower by >0.10 (5-yr CC is noisy), dLOD >= 0.994, backbone in range,
manifold CC >= 0.5 vs the alt reference; rejection reasons logged.
Scratch copies in _refine/<quad>, backups _refine/backup/, results
_refine/results.jsonl, refine.log, _refine/REFINE_DONE. Resumable.
First launch (1000 s, 2 threads x 7) was stopped: it exposed two more
lt.exe problems, both fixed 11:20 build:
 - File_Lines/Make_Data opened reference files without "shared=yes"; two
   tasks reading the dLOD file at once (TEST_ONLY: both threads save) hit
   USE_ERROR "reopening shared file" -> 0 records.
 - Regression_Factors' unconditional Accel_Ref := Forcing(First).Date
   (added 09-30) then raised CONSTRAINT_ERROR on that empty array, killing
   the thread before its dLOD line; now guarded.
 The altquad overnight sweep was slowed by the same thing (runner waited
 for two dLOD lines). Runners now wait for one dLOD line + quiet output,
 retry reading lte_results.csv until complete, and kill lt.exe by cwd.

**CC refinement sweep of Oct2026_altquad DONE (2026-10-01 13:28).**
49/89 quads improved and installed, 40 unchanged (38 rejected only for the
2000-05 test CC dropping >0.10, 1 validate, 1 dLOD 0.9931). Accepted: fit
CC +0.06 median, train +0.07, validate +0.12, test -0.01. Global vs
Kaplan: monthly 0.927/0.703 -> 0.940/0.768 (raw/detrended), 12-mo
0.954/0.776 -> 0.963/0.838. The per-quad test gate matters: taking the 38
test-only rejections too would raise the global detrended CC (0.782 ->
0.802 vs prorated data) but drop the GLOBAL held-out 2000-05 CC from
+0.31 to +0.12. Held-out 2000-05 global window: Oct2026_altquad CC +0.31
(detrended -0.18, RMSE 0.105) vs Feb2026_sweep +0.83 (+0.52, RMSE 0.053)
-- but Feb2026_sweep's 2000-05 scores have been used for accept decisions
over several rounds, so that window is no longer pristine there.
ANOTHER lt.exe FIX (13:25 build): Save_Mutex serializes the final save
block. In TEST_ONLY both threads save; with shared reads one thread read
dlod_ref.dat while the other rewrote it -> END_ERROR, process aborted
mid-save, 0-byte lte_results.csv (26/89 re-scores). sweep_cc_refine.py
--rescore-failed redid just those re-scores (all 26 reported cleanly).

**ALPHA piecewise manifold-minimum gain (2026-10-02, uncommitted; build
01:37).** ALPHA (env/resp, default 0 = code skipped, output identical):
winding sin/cos terms are multiplied by ALPHA_GAIN (default 2.0) while
F - M_min <= ALPHA*|M_min| (manifold within an ALPHA fraction of its
minimum negative excursion), unchanged otherwise; same function
(Alpha_Gain) in Regression_Factors and LTE; M_min = min of the candidate's
whole manifold; JSON secular_context has alpha, m_min. (A first version,
exp(-ALPHA*|F-M_min|) on all windings, only removed response away from the
minimum and was replaced.) Purpose: isolate conditions for super El Ninos.
Findings (alt-quadrature kN000_W090/W110/W130 mean, and Sep2026/nino12):
manifold min -40.9 at 1978.8, within 10% of it 1976-2003. Events switch on
as ALPHA passes their manifold position: 1976 (6% above min) and 1987-88
(5%) at 0.1 -- landing on the data (nino12 1976: 0.92 -> 1.53, data 1.53;
quads 1987-88: 0.83 -> 1.29, data 1.18); 1997-98 (15%) at 0.2; 1982-83
(20%) at 0.3; 1972 (34%) at 0.4. With gain 2 the overall fit is ~flat
(nino12 CC 0.545 -> 0.550, skew 0.13 -> 0.39 vs data 1.21) but 1982-83 /
1997-98 only reach ~1.2-1.5 vs ~4 observed. Raising ALPHA_GAIN (3,4,6)
lifts them (to ~2.0-2.3) but overshoots 1976 (2.1 vs 1.53), zeroes events
outside the window (2015-16, 1957) and lowers CC (0.47 at gain 6). So
"manifold near its minimum" is necessary-looking for 1982-83/1997-98 but
not sufficient (1976 shares it) and not necessary for 2015-16, 2023,
1957, 1888 (manifold 75-105% above min).
NOTE: Sep2026/nino12's saved lt.exe.p re-evaluates to fit CC 0.545, not
the 0.636 of its lte_results.csv (pre-fix harmonic-save bug).

**ALPHA gain now POSITIVE-ONLY (2026-10-02, build 02:35).** Per the user:
inside the ALPHA window the gain applies only when the summed winding (LTE
modulation) term is positive; negative excursions are unchanged. The sign
depends on the fitted amplitudes, so Regression_Factors does two passes
(pass 1 without gain -> sign of the winding sum per row -> scale the
winding columns of qualifying rows -> pass 2 is the fit); LTE applies the
same rule to its own winding sum (WS). ALPHA=0 verified identical
(kN000_W090, diff 0.0). Sep2026/nino12 re-run in place with ALPHA 0.1 (the user had edited the resp to ALPHA 0.1, ALPHA_GAIN 2.0, VALIDATE off)
(300 s, user's lte_run.sh settings): plot CC 0.649 (train 0.656,
validation 0.41) vs 0.636 without ALPHA; model max 2.13 (was 1.68), model
min unchanged (-2.01); skew 0.49 (0.38; data 1.21). Peaks: 1976 1.97
(data 1.53), 1982-83 1.64 (3.91), 1997-98 1.01 (4.05). Pre-ALPHA folder
kept as Sep2026/nino12_pre_alpha/.

**Overnight ALPHA update of Oct2026_altquad launched 2026-10-02 03:07.**
User: positive-only ALPHA "works excellent" on nino12 at ALPHA 0.25,
ALPHA_GAIN 3.0. Oct2026_altquad/run_alpha_overnight.sh: two back-to-back
passes (A, B) of Sep2026/sweep_cc_refine.py --ALL, 300 s, 4 threads x 4
workers, with --force-set ALPHA=0.25 --force-set ALPHA_GAIN=3.0 on 15
El Nino pathway quads only: equatorial Pacific kN000_E170/W170/W150/W130/
W110/W090; N. America west coast kN020_W090 (box also covers part of the
Gulf of Mexico), kN020_W110, kN040_W130, kN060_W130, kN060_W150; S.
America west coast kS020_W090, kS020_W070, kS040_W070, kS060_W070. Forced
quads always end with the setting (searched fit if it beats the baseline
scored WITH ALPHA, else that baseline: action kept_old_with_settings).
The other 74 quads: plain CC refinement. Outputs: refine_alpha_{A,B}.log,
_refine_alpha_{A,B}/ (results.jsonl, backups), global_sst_alpha_{A,B}.txt,
global_sst_from_quads_alpha_{A,B}.png, ALPHA_OVERNIGHT_DONE. Binary:
Oct2026_altquad/lt.exe = 02:35 build (previous kept as
lt.exe.build_1001_1325). Earlier pass-2 records in _refine_pass2/.
Feb2026_sweep not touched.

**Which quadrature? evidence + path (2026-10-02 morning).**
ALPHA overnight on Oct2026_altquad finished 07:45 (passes A,B: 61 and 65
improved); pathway El Nino peaks inside the manifold-min window up (eq
Pacific 6-quad mean: 1982-83 0.72->0.95 of 1.47, 1987-88 0.68->0.98 of
1.02, 1997-98 0.64->1.26 of 1.71), outside it down (1957, 1972, 2015).
Global alt: monthly 0.949/0.811, held-out 2000-05 CC +0.16 (was +0.31,
+0.25: falling with each refinement pass).
Comparison original (Feb2026_sweep) vs alternate (Oct2026_altquad):
per-quad median detrended fit 0.620 vs 0.668 (alt better 55/89), per-quad
held-out 0.29 vs 0.35; windings 8 vs 6; BUT global detrended CC 0.875 vs
0.825 and global held-out 2000-05 +0.83 vs +0.16. Confounded: alt has had
4 CC refinement passes with the fixed binary, original none. Same-quad
residual CC between sets 0.66, global residual CC 0.75 (misfit mostly
common). Basins: alt better in the Atlantic (+0.08..+0.15), original
better in eq/N Pacific. Weak in both: kS040_E130, kS040_W070, kS040_W150,
kN040_E150, kN040_E130, kN040_E010, kN020_W010, kN000_W010, kS020_E030.
Pathway skew (data/orig/alt): eq Pacific +0.27/+0.05/+0.52, N America
coast +0.17/-0.14/+0.26, S America coast +0.03/+0.11/+0.43 -> ALPHA
removed from the 4 S. America coast quads (re-fit, ALPHA 0.0, all 4
improved); ALPHA stays on 11 quads. Manifolds: CC(orig, d/dt alt)=+0.50,
CC(alt, d/dt orig)=-0.04; the 4 largest El Ninos sit at low |slope| of
the original manifold (<=28th pct), moderate ones at high slope (n=10).
1880-1949 back-cast, alt (Oct2026_altquad/_backcast, all 89 manifolds
reproduce exactly) vs original: raw +0.514 vs +0.445; no-trend raw +0.339
vs +0.489; lin-detrended +0.42 vs +0.37; 12-mo +0.50 both; bands 2-7 yr
+0.19 vs +0.08, 7-20 yr -0.52 vs +0.17, >20 yr +0.85 vs +0.77. Plots
Oct2026_altquad/backcast_1880_detrended_both.png, backcast_1880_both.png.
IN PROGRESS (launched 09:17): Feb2026_sweep/run_refine_overnight.sh --
two CC refinement passes of the original set (equal effort); outputs
refine_{A,B}.log, _refine_{A,B}/, global_sst_refine_{A,B}.txt,
REFINE_AB_DONE. Then compare on global + held-out + back-cast.

**Equal-effort quadrature comparison (2026-10-02 15:10) -- DONE.**
Feb2026_sweep refinement finished 13:52 (pass A 63, pass B 64 improved;
backups in _refine_A/, _refine_B/). Both back-casts refreshed.
Global (Kaplan, 1950-2023) original vs alternate: monthly raw/det
0.960/0.880 vs 0.949/0.813; 12-mo 0.977/0.938 vs 0.971/0.879. Held-out
2000-05 global CC +0.84 (RMSE 0.047) vs +0.18 (0.107); excluding
kN000_W150/W170 (1980-85 gap in original) +0.74 vs +0.13. Per-quad
median det fit 0.676 vs 0.668 (alt better 35/89; held-out 47/89).
Basins: alt better S/eq Atlantic; original better N/eq Pacific.
Pathway skew data/orig/alt: eq Pac +0.27/+0.05/+0.52; N Am coast
+0.17/-0.11/+0.26; S Am coast (ALPHA off) +0.03/+0.08/+0.17.
1880-1949 back-cast lin-detr monthly: original +0.32 (was +0.37 BEFORE
refinement), alternate +0.42; bands 2-7 yr -0.01 (was +0.08) vs +0.19;
7-20 yr -0.01 (was +0.17) vs -0.52; >20 yr +0.79 vs +0.84. So CC
refinement cost the original its small pre-1950 interannual skill.
Weak in both: kS040_E130, kN020_W010, kS040_W150, kN000_W010, kN000_E070.
OPEN (user to decide): step 5 = low-|dM/dt| gain for the original
manifold (analogue of ALPHA); which quadrature to carry forward.

**Slope-keyed ALPHA + forward prediction (2026-10-02 evening).**
Ada: ALPHA_SLOPE (months; window tested on signed manifold slope) and
ALPHA_SIGN (-1 falling, default) added, build 15:38 (previous kept as
Feb2026/enso_opt.build_1002_0235). Scratch refits of kN000_W090/W110/W130
(Feb2026_sweep/_slope/root_*): falling side is right; best W=13, ALPHA
0.7, gain 2 (train CC 0.817->0.842, 0.806->0.825, 0.804->0.807; skew
matches data; 1957/82/87/97 peaks near observed; 1972 and 2015 shrink;
2000-05 held-out drops). Not installed in real quads. ALPHA ON HOLD per
user. Forward prediction: experiments/Oct2026_forecast (forecast.py,
verify.py, forecast_vs_observed.png, sstoi.indices). Observed Aug 2026:
Nino3.4 +2.52, Nino3 +3.13, Nino1+2 +4.08. No model reaches it (Jun-Aug
2026 model -0.4..+0.7); alt quads W110/W130/W150 + alt_nino3 forecast CC
+0.5..+0.8, original quads ~0; orig_nino4 best single (CC +0.62).

**LAYER two-stage manifold (2026-10-02 evening, Feb2026/_layer/).**
Env/resp LAYER (0 off; enum added): lowest winding M(1) leaves the MLR and
builds a new manifold F2 from the tidal F after its first-order impA/impB
Bessel correction: LAYER=1 F2 = offs*sin(2pi k1 Fb) + bg*cos(2pi k1 Fb),
LAYER=2 F2 = Fb + that; higher windings regressed on F2; offs/bg (formerly
only 2nd-order Bessel factors, dropped in this mode) hold the layer
amplitude (user's suggestion; v1 used impA/impB and lost the Bessel
correction). Build 18:16. AMO results: baseline/control (300 s search)
0.780/0.624/0.572 -> 0.781/0.619/0.583 (train/val/test). v2 grid (168
pts, fixed params): best 0.705 (A 10, k1 x0.75) vs 0.780. 300 s searches
from the 3 best: 0.759/0.566/0.516, 0.755/0.420/0.397, 0.786/0.536/0.331
-- one beats the control in training, none out of sample. The search
barely moves offs/bg/k1 (one-of-~150 small relative steps). LAYER=2 not
yet tested in v2. Scripts layer_test.py, layer_seed_scan.py,
layer_grid.py (v1), layer_grid2.py (v2); results grid2_amo.json.

**LAYER sweep of all quads launched 2026-10-03 01:05 (experiments/Oct2026_layer).**
run_layer_sweep.py: seeds = user's Feb2026 LAYER=1 fits amo, pdo, nino4 (.p AND
.resp; ltep = [layer k1 ~0.021, Bessel k ~0.232, regressed ~3.8-4.1 (+6.1 nino4)]).
Per quad: 200 s CC search per seed, TEST_ONLY re-score, usable if dLOD >= 0.99,
choose best mean(train, validate); then stage 2 adds one harmonic of the last
regressed winding (kept if train +0.02, validate -0.03 max, dLOD >= 0.99).
4 workers x 4 threads, ~6 h. Frozen binary Oct2026_layer/lt.exe. Outputs
results.jsonl, sweep.log, SWEEP_DONE, installed quads in Oct2026_layer/<quad>.
Also: AMO/PDO common-manifold round-robin in Feb2026/_rr (joint 0.490 -> 0.690;
PDO searches push dLOD < 0.994; round 4 the balanced state). Smoke test
kN040_W050: 0.782/0.755/0.712 with nino4 seed + harmonic.

**Warm-side IR paired pass launched 2026-10-03 12:13 (experiments/Oct2026_layer_ir).**
New Ada options since the LAYER commit (uncommitted, all off by default): LAYER_MIX
(impC as tidal-manifold mixing weight w under LAYER=1; opt-in because existing fits
carry small impC impulses), LAYER_HP (high-pass the mixed Fb), IR_SIDE (+1 warm-side
only / -1 cold-side only 12-month delay differential, "warm" = above own IR_BASE=121
month running mean, data and model sides). Findings: w mixing never helps the eastern
Pacific quads (raw or high-passed); symmetric IR weak; warm-only IR seeded -0.2 beat
IR=0 searches on 2000-05 test in 10/10 targets (amo, pdo, 8 Pacific quads).
run_ir_pass.py: per quad, from Oct2026_layer, a 300 s control search (IR 0) and a
300 s warm-only search (IR_SIDE 1, IR -0.2), installed in ctrl/ and warm/; Pacific
quads first; ~4.3 h; results.jsonl, pass.log, PASS_DONE. Then compare both sets
globally (held-out + 1880 back-cast) at equal effort.

**ANNUAL_DITHER set aside (2026-10-03).** Dummy annual signal added to the manifold
(ANNUAL_DITHER a, _PHASE, _STAGE 1 after layer / 0 before) with climatology removal of model
output and regression columns. Grid (Oct2026_layer/_dither/grid.log): all a > 0 worse;
the modulated annual sidebands survive climatology removal. Code stays, off by default.
Focus: warm-only IR paired pass (Oct2026_layer_ir).

## NEW GOVERNING PRINCIPLE: judge every future change by the COLLECTIVE measure, not per-cell scores (2026-09-28 ~21:30)

Explicit reframing from the user: the objective from here on is to move
the whole 89-quad collective toward a better shared parameterization --
"the FIRST time 90 spatially separated SST time-series have been
modelled with a handful of standing-wave parameters." A few individually
weak cells are fine. A uniform shared bias in any parameter (dLOD,
year-length-in-days, any other `.p` scalar) is fine too, as long as it's
uniform across the collective -- what matters is the collective measure,
not each cell hitting its own independently-best CC/DTW/dLOD.

**The collective measure**: `mlr_shared_windings.py`'s overall geo-smooth
R^2 (89 quads' real manifolds explained by the confirmed 20-winding
shared harmonic/subharmonic set + a degree-1 geographic smoothness
constraint on amplitude/phase -- 338 total parameters for ~78,000 data
points). New wrapper `collective_r2_log.py "<label>"` runs it and appends
one line to `collective_r2_log.jsonl` -- **run this before and after any
future batch, LOCKW change, or manual-fix round, and compare consecutive
entries** as the actual go/no-go signal, instead of counting how many
individual cells pass some per-cell threshold.

**Baseline logged 2026-09-28 ~21:30** (after 89/89 quads solved +
manifold-outlier fixes + LOCKW batch #1): overall geo-smooth R^2 =
**0.9341**, independent-per-quad mean = 0.9346 (min 0.6125, max 0.9658).
Notably this is DOWN from 0.9461 (measured earlier today, before LOCKW
batch #1) -- a real, concrete illustration of the tension the user is
describing: LOCKW batch #1 improved several individual cells' own CC/DTW
scores but slightly reduced cross-quad geo-smooth coherence. Worth
digging into which of the 5 LOCKW-batch fixes caused the drop before
running further batches blind -- possible that one or more landed on a
locally-good-scoring but collectively-incoherent harmonic (a smaller-
scale echo of the earlier "metastable"/`kS040_W030` pattern).

## LOCKW batch #1: 5/8 weak cells improved (2026-09-28 ~16:00-19:30)

First real-world validation of `LOCKW` (see below), per the user's
"more sweeps, improve collective fit/CV AND reinforce canonical
backbone" instruction. `sweep_lockw_batch.py` targeted the 8 weakest-
scoring currently-accepted cells (lowest min(train/validate/test) or
min(pair)), donor-seeded each from its nearest already-solved neighbor,
ran with `LOCKW=TRUE` at a longer 400s internal timeout (LOCKW needs
more iterations to re-adapt around a locked winding than the usual
120s), and kept the result only if its minimum score beat the original.

**5/8 improved and committed, all with confirmed canonical backbones**:
- `kS020_W090`: -1.807 -> 0.120 (this was the session's one documented
  pathological case -- validate=0.0 exactly, test~-1.8 -- finally
  resolved)
- `kN020_W090`: 0.014 -> 0.694
- `kS040_W010`: 0.043 -> 0.384 (backbone landed on 0.414 = 2x canonical,
  a valid harmonic)
- `kS040_W070`: 0.072 -> 0.416
- `kN060_W130`: 0.135 -> 0.244

**3/8 failed and rolled back cleanly** (no data lost, kept original):
`kN040_W010` (-0.200, unchanged), `kN020_W030` (0.017, unchanged),
`kS040_E170` (0.164, unchanged) -- each after a very long cascade
(~3200-3600s across all 4 DTW/CC x VALIDATE variants + retry) that never
found a canonical-locked fit beating the original. These 3 may need a
different donor, a manual fit, or may simply be cells whose real local
data doesn't fit any canonical harmonic well (worth investigating
individually if pursued further, same as the earlier hard-cell pattern
this session).

Script (`sweep_lockw_batch.py`) is reusable for further batches --
change `TARGET_CELLS` to whatever's next (e.g. re-run the remaining
weak-score cells from `quad_status_map.png`'s amber/salmon set, or the
3 still-flagged `kS040_W130/150/170` once the user's pending Feb2026
update lands).

## NEW: canonical-backbone winding lock (`LOCKW`) added to lt.exe (2026-09-28 ~15:30)

Real Ada engineering, not a Python workaround. Goal (per the user): support
"more sweeps" whose search is biased toward a canonical fundamental
backbone with harmonics/subharmonics, for parsimony and easy
parameterization -- following up on the MLR investigation that found the
shared 20-winding harmonic set explains ~94% of manifold variance per
quad (see below).

**What it does**: constrains the `ltep` (winding/backbone) slice of the
random-descent search's parameter vector to snap to the nearest value in
a caller-supplied canonical set (backbone x harmonics/subharmonics),
every time Markov perturbs one of those slots -- without touching the
generic `Markov` procedure (which stays fully unconstrained for every
other parameter) and without changing default behavior at all (opt-in,
off by default).

**Checked first**: the old `~/refactor/gem-lte-core` branch's
`Phase_Lock`/`Enhanced_Nonlinear_Structure` (flagged in
`REFACTOR_RECOVERY_PLAN.md` as unwired/never call-sited) looked
superficially related but isn't -- it's a fixed seasonal amplitude
modulation and a quadratic overtide term, neither of which touches which
winding *frequencies* the search explores. Not reusable for this; built
fresh instead.

**New env vars** (all opt-in, default off/matching the pre-existing
canonical set from the MLR work):
- `LOCKW` (default False) -- enable the lock.
- `CANON_BACKBONE` (default 0.207) -- the fundamental.
- `CANON_HARMONICS` (default `"1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 18 19 23 27"`,
  space-separated like the existing `NH` convention) -- integer
  multiples of the backbone to allow.
- `CANON_SUBHARMONIC` (default 11) -- adds `backbone/N` as one more
  candidate (this default matches the confirmed ~0.0188 slow "120/60yr
  rectification" subharmonic).

**Implementation** (`src/gem-random_descent.ads/.adb`,
`src/gem-lte-primitives-solution.adb`): new `Nearest_In_Set` utility
function in `GEM.Random_Descent` (unit-tested standalone, confirmed
correct incl. sign handling). `ltep` occupies exactly the last `NM`
entries of the search's `Set` array (`Size_Shared - NM + 1 ..
Size_Shared`, matches `Param_B_Overlay.First_LT_Index`). Since `Markov`
perturbs one parameter per call without telling its caller which index
changed, the snap is applied by the CALLER (the model-specific search
loop) right after each `Walker.Markov (Set, Keep, Spread, Set0)` call: it
scans the `ltep` index range for whichever entry differs from `Keep`
(the pre-perturbation backup) and snaps just that one. Compiles clean
(`gprbuild -p -P lte.gpr`, only 2 pre-existing unrelated warnings).

**Validated behavior** (see full test writeup in this session's
conversation, not reproduced here): works correctly as a **drift-
prevention guardrail on an already-reasonable (e.g. donor-seeded from a
near-canonical neighbor) starting point** -- touched slots snap to exact
canonical values, fit quality holds. Does **NOT** reliably pull an
arbitrary/uncorrelated cold-start seed onto canonical values within
practical run times (confirmed at 30s/240s/1000s, all landing at the
same off-canonical local optimum) -- physically expected, since winding
frequency error integrates into large phase drift over a ~73-year
record, and single-parameter hill-climbing can't jointly re-adapt the
correlated scalars (impA/impB/delA/delB/asym/etc.) fast enough to escape
a well-adapted off-canonical optimum. **Recommended usage pattern for
the next round of sweeps**: donor-seed from an already-canonical
neighbor (there are many now -- 0.207-backbone gold standards are the
majority of the 89 solved quads) as usual, then run with `LOCKW=TRUE` to
prevent the kind of alias-drift that caused this whole session's
backbone-invalidation problem in the first place, rather than expecting
it to fix a cell that's never been close to canonical before.

**Also fixed in passing**: `.gitignore` now excludes `enso_opt` (the
actual file `lt.exe` symlinks to -- confirmed by rebuild; was previously
untracked-but-not-ignored, a gap since `lt.exe` itself has always been
gitignored), `iir_invariant_test`, `node_modules/`, and `obj/`.

## FLAGGED FOR RE-EVALUATION: kS040_W130, kS040_W150, kS040_W170 (2026-09-28 ~13:15)

After the manifold-outlier fixes, re-running `manifold_comparison_all_quads.py`
(RAW corr 0.937, DETRENDED 0.591 -- both way up) surfaced a new, different
kind of anomaly: at year 2020.0 specifically, these 3 adjacent South
Pacific quads (all 40S) sit ~4 std below the ensemble mean (-0.577 vs
mean 0.829), while every other quad is clustered tightly positive. Not a
data-edge artifact (all 3 have full records to 2023). Traced to a direct
donor chain: `kS040_W130` -> `kS040_W150` -> `kS040_W170` (each seeded
from the previous, all from the ORIGINAL automated sweep pipeline
pre-dating this session) -- so the dip is either a genuine shared
regional signal (physically plausible, adjacent quads) or an inherited
artifact from the chain's root (`kS040_W130`) that never got
independently re-checked against each cell's own data near the record's
end. Ledger entries for all 3 now carry a `flag: needs-re-evaluation`.
**User has better results pending in Feb2026** -- re-verify these 3 (same
process as every other manual fix this session: clear stale checkpoint
first, TEST_ONLY verify, sync, check manifold coherence AND the 2020.0
value specifically) once available.

## INCIDENT: manual TEST_ONLY verification corrupted a fresh update (2026-09-28 ~12:45)

The stale-per-CLIMATE_INDEX-checkpoint bug (documented earlier, fixed in
`sweep.py`'s own `clear_stale_checkpoint`) bit again -- this time in my
own ad-hoc MANUAL verification workflow, which every gold-standard fix
this whole session used (`cp lte_run.sh ...; sed TEST_ONLY=true; expect
watch_run.exp; mv back`) and which never cleared `lt.exe.<cell>.dat.p`
first. The user gave a fresh `kS020_W030`-seeded update for
`kS040_W030`; I verified it via TEST_ONLY without deleting the stale
checkpoint left over from an earlier attempt, and the 2-thread run's
`Report_Final` wrote back the STALE state over the fresh one -- silently
destroying the update in both `Sep2026/` and (via my own Feb-sync)
`Feb2026/` before I'd recorded its true scores. **Recovered**, not lost:
the user's `lte_gui.py` had already auto-saved timestamped backups
(`lt.exe.p.09-28-2026-2-2000-2005`), which I restored from -- very close
to but NOT byte-identical to the true lost final state (scores 0.532/
0.536 recovered vs 0.707/0.661 seen briefly before the corruption).

**Lesson: any manual TEST_ONLY (or real) verification of a cell must
delete `lt.exe.<cell>.dat.p` first, every time, not just when
`sweep.py`'s own automated code paths run.** This was already true in
principle from the earlier finding, but hadn't been consistently applied
to the manual workflow used for essentially every "gold standard" fix
recorded this session -- those may be fine (most cells' checkpoint and
`lt.exe.p` were written together by the same run, so no staleness), but
this is now a confirmed live risk for ANY cell where a NEW update lands
while an OLDER per-index checkpoint is still sitting in the same
directory. If in doubt about a past verification, re-check for a stale
checkpoint before trusting the recorded scores.

## MANIFOLD COHERENCE CHECK + OUTLIER REFIT (2026-09-28 ~11:00-12:00)

Post-completion sanity check, per the user: compare `lte_results.csv`
column 4 (the "hidden latent forcing" / autonomous manifold) across all
89 quads, adapting `Feb2026/plot_hidden_latent_forcing.py` (originally
hardcoded to 7 flagships) to enumerate all of them. New script:
`Sep2026/plot_hidden_latent_forcing_all_quads.py`, output
`Sep2026/manifold_comparison_all_quads.png` (spaghetti overlay + raw
correlation heatmap + detrended correlation heatmap, geographically
ordered).

**Key finding, corrected once:** raw mean pairwise correlation looked
very high (0.764) but a strong shared secular trend across all 89
manifolds was inflating it -- per this project's own documented lesson
(memory: decompose a nuisance trend before crediting a correlation),
detrending each series individually before correlating dropped the mean
to 0.525. Still clearly real coherence, just not as dramatic as the raw
number suggested. **6 clear outlier cells** (low coherence in both raw
and detrended views): `kS040_W030`, `kS020_E150`, `kN040_E150`,
`kN000_W150`, `kS020_E010`, `kN020_E150` -- notably NOT weak fits by
their own train/validate/test scores (all solid, 0.5-0.87), so this
wasn't a data-quality tell on its own.

**User's call**: these 6 show nothing physically distinct from their
neighbors, so refit them (donor-seed from a clean proximal quad) rather
than treat the divergence as real local physics. `Sep2026/
refit_manifold_outliers.py` did this, judging success on BOTH
criteria -- not just CC/DTW score (already fine before!) but whether the
detrended correlation to a leave-these-6-out ensemble reference actually
improves by a real margin.

**Results: 3/6 fixed, 3/6 couldn't be** (kept at their original good
fit, not forced):
- `kS020_E150`: coherence 0.044->0.487, fit 0.634/0.745/0.659 -- committed
- `kS020_E010`: coherence 0.173->0.499, fit 0.745/0.849/0.768 -- committed
- `kN020_E150`: see bug note below -- final coherence 0.147->0.649, fit
  0.416/0.420/0.415 (moderate but balanced) -- committed
- `kS040_W030`: donor-seeded refit succeeded on score (0.698/0.848/0.732)
  but coherence barely moved (0.009->0.054) -- rolled back, kept original
- `kN040_E150`, `kN000_W150`: donor-seeded cascade failed outright (even
  `kN000_W150` seeded from `nino34` at just 556km) -- rolled back, kept
  original. These 2 plus `kS040_W030` remain genuine outliers.

Net effect: RAW mean corr 0.764->0.836, DETRENDED mean corr 0.525->0.546.

**Bug found and fixed in `sweep.py`**: `solve_cell_triangulated_with_retry`
had no self-exclusion -- when refitting `kN020_E150` (which already had
an accepted ledger entry), it silently picked `kN020_E150` itself as one
of its own 3 triangulation donors, producing a near-degenerate fit
(train=0.16/validate=0.13/test=0.08) that I almost left committed because
its COHERENCE looked great (reinforcing its own already-bad state looks
"coherent" with a coherence check that doesn't verify fit quality). Caught
this by eyeballing the fit scores, not automatically -- **the refit
script's automated commit criterion checked coherence improvement only,
with no fit-quality floor, which is a real gap if this pattern gets
reused.** Fixed `solve_cell_triangulated_with_retry` to always exclude
the cell itself (`exclude = set(exclude) | {cell}`) regardless of caller.
Re-ran cleanly afterward for the final committed `kN020_E150` result.

## SWEEP COMPLETE (2026-09-28 ~10:40): 89/89 QUADS SOLVED

Every single k_sst grid quad now has an accepted fit. Summary of how the
final 36 (previously untouched by `sweep.py` at all) got there, per the
user's instruction to try each cell's own prior Feb2026 work before
falling back to donor-seeding:

- **34/36** solved directly from their own pre-existing Feb2026
  parameter set (predating this sweep session) — `verify_untouched_36.py`
  batch-tested both DTW and CC (TEST_ONLY, ~20s/check) against each
  cell's untouched prior fit and kept whichever cleared the dLOD floor
  with the better minimum score; `finalize_verified_36.py` then wrote
  the winning config, synced to `Sep2026/`, and recorded the ledger in
  one batch pass. **DTW won 25/34, CC won 9/34** — reinforcing this
  session's earlier DTW-favoring finding, though CC still wins often
  enough it's worth checking both, not skipping straight to DTW.
- **2/36** (`kN060_W150`, `kN080_E010`) had decent correlation scores
  from their own prior fit but failed the dLOD floor (0.745 and 0.978
  vs the 0.994 floor) — fell back to `sweep.py`'s normal donor-seeding
  cascade. `kN060_W150` succeeded that way (dLOD=0.99543).
  `kN080_E010` did NOT (donor-seeded from `kN060_E010`, stuck at
  dLOD=0.936 across all 4 cascade attempts + retry) — this one
  genuinely needed the user's own manual fix. **User's insight**: this
  cell (80N/10E, Svalbard/Arctic) is physically unique, with a strong
  annual component that comes and goes (intermittent/amplitude-modulated
  seasonal cycle, plausible for Arctic sea-ice-covariant SST) —
  accounting for that got dLOD to 0.99609 cleanly. This was the very
  last cell solved, completing the full 89/89.

**Verify current live count** (should read 89):
`python3 -c "import sys; sys.path.insert(0,'.'); import sweep;
print(len([e for e in sweep.latest_entry_per_cell(sweep.load_sweep_ledger()).values() if e.get('accepted')]))"`
from `Sep2026/`.

**What's NOT necessarily done**: "accepted" here means dLOD floor +
backbone_ok + a real score reported — it does NOT mean every fit is
equally strong. Several (especially among the last 36) carry the
familiar weak-test or weak-pair caveats (see `quad_status_map.png`,
regenerate via `python3 make_status_map.py`, for the amber/salmon
cells). If further quality passes are wanted, that map plus the ledger
notes are the place to find which ones to revisit.

## ALL 6 REMAINING CELLS SOLVED (2026-09-28 ~09:00): 53/89 quads accepted

After the backlog run finished with 6 still bad (see below), the user
gave a key instruction that cleared all 6: **for a hard cell, also try
its OWN prior Feb2026 parameter set (predates this sweep session) as an
alternative to a nearby-quad donor seed, and keep whichever gives the
better fit/CV.** Every single one of the 6 remaining cells
(`kS020_E170`, `kN020_E090`, `kN020_E070`*, `kN020_E110`, `kS040_E150`,
`kN000_E150`) turned out to have real, pre-existing Feb2026 work
(backup snapshots dated Sep 20-26, well before this sweep started) that
was **immediately better than every donor-seeded automated attempt** —
in 4/6 cases just a `TEST_ONLY` check of the untouched existing state
was enough, no new search needed. Where both DTW and CC were tried
(4 cells), **DTW won every time**, several by a wide margin (e.g.
`kS040_E150` test=0.853 DTW vs 0.302 CC).

(*`kN020_E070` was fixed slightly differently, via `winding_scalogram.py`
identifying a real M=1.15 ridge — see below — not simply reusing a
pre-existing fit untouched.)

**Practical lesson for the 36 still-untouched quads and any future
backlog**: check each cell's own Feb2026 directory for a pre-existing
`lt.exe.p` FIRST (and reconstruct a missing `lte_run.sh` from the
standard template if needed — 4 of these 6 had none) before reaching for
donor-seeding. `sweep.py`'s automated flood-fill never tries this at
all — it only ever seeds from OTHER cells, never from a cell's own
history. Also worth checking: whether `Sep2026/<cell>/` is currently
holding a corrupted leftover from an old automated attempt that needs
overwriting from the good Feb2026 state (this was true for all 4 of the
"just verify" cases here: backbones found were 200.0, 18.68, 150.0,
10.0, 9.996 -- all now overwritten).

`python3 -c "import sys; sys.path.insert(0,'.'); import sweep;
print(len([e for e in sweep.latest_entry_per_cell(sweep.load_sweep_ledger()).values() if e.get('accepted')]))"`
from `Sep2026/` reconfirms 53 solved live.

## RUN COMPLETE (2026-09-28 04:16): 47 solved, 6 still bad

The `resume_refix.py` backlog run (started after the checkpoint-bug fix,
see below) finished cleanly. Final: **47 total cells accepted, 6 genuinely
still bad**: `kS020_E170`, `kN020_E090`, `kN020_E070`, `kN020_E110`,
`kS040_E150`, `kN000_E150` — all clustered in the eastern Indian
Ocean/Maritime Continent region (N/S 020, E 070-150), suggesting that
specific area may be a real trouble spot worth a deliberate manual look,
same as the shelf/channel gold standards. Up from 25 "genuinely
trustworthy" cells at the start of this handoff.

**UPDATE 2026-09-28 ~08:15: `kN020_E070` manually fixed, 5 remain in the
cluster.** The user started from this session's own automated flood-fill
seed, then used `winding_scalogram.py` to visually identify a strong,
continuous ridge candidate at **M=1.15** and set it as `ltep[0]` —
correctly placed *before* the backbone slot `ltep[NM-1]`, so the shared
0.2073 tidal backbone is preserved there untouched (`ltep=[1.150, 0.231,
0.207]`, NM=3). DTW, VALIDATE not set (pair). Result: pair=(0.725,
0.693), dLOD=0.99608, backbone=0.2073 — clean, no caveat. Recorded as
`manual-gold-standard`. Remaining cluster: `kS020_E170`, `kN020_E090`,
`kN020_E110`, `kS040_E150`, `kN000_E150`.

`python3 -c "import sys; sys.path.insert(0,'.'); import
sweep; print(len([e for e in sweep.latest_entry_per_cell(sweep.load_sweep_ledger()).values() if e.get('accepted')]))"`
from `Sep2026/` reconfirms the live solved-count (48 as of the above fix).

**Status map**: `Sep2026/make_status_map.py` generates
`Sep2026/quad_status_map.png`, a full-world Mercator plot of all 89 quads
color-coded by fit status (good/weak/negative/unsolved/untouched) with
train/validate/test (or pair) labels on each box — re-run it after any
further manual fixes to keep it current (it's a point-in-time snapshot,
not auto-regenerated). Mercator can't reach the poles: keep
`ax.set_extent`'s northern bound at 80, NOT 85 or 90 — both broke the
whole projection (a Mercator polar singularity), confirmed directly this
session.

Context for a fresh session picking this up. The prior session's context is
very long; this file plus `sweep.py`'s own inline comments should be
sufficient to continue without re-reading the whole history.

## What this is

An automated sweep that fits the ~89 Kaplan SST grid cells (`kN020_E050`-style
directories in `experiments/Feb2026/`) using the real GEM-LTE stochastic
optimizer (`lt.exe`), seeding each cell's `lt.exe.p` from a nearby
already-solved donor (a "flagship" named index, or another solved grid cell),
flood-filling outward. Fully isolated in `experiments/Sep2026/` — never
writes into `experiments/Feb2026/` except to *read* donor `.p`/`.dat`/
`lte_results.csv` files.

Core files:
- `sweep.py` — all the machinery (geo helpers, seeding, the attempt cascade,
  the backbone sanity gate, the ledger).
- `watch_run.exp` — TCL Expect script that spawns `lte_run.sh` under a real
  PTY. **Required**, not optional — Ada's Text_IO fully block-buffers when
  stdout isn't a TTY, so a plain `subprocess.PIPE` sees nothing until the
  process exits (or never, if killed first). Emits tagged lines
  (`MATCHED_TRIPLET:`, `MATCHED_PAIR:`, `MATCHED_DEADLOCK:`, `MATCHED_DLOD:`)
  that `sweep.py` parses.
- `sweep_ledger.json` — this sweep's own append-only log. **Read it with
  `sweep.latest_entry_per_cell(sweep.load_sweep_ledger())`**, not raw — a
  cell can have multiple entries over time (retries, invalidations), and
  only the *last* one is authoritative.

## THE *REAL* ONE THING TO FIX FIRST (found 2026-09-27 ~20:15, supersedes the N=1-vs-N=2 story below)

**`lt.exe` writes a per-CLIMATE_INDEX checkpoint file, `lt.exe.<cell>.dat.p`
(e.g. `lt.exe.kN040_W170.dat.p`), alongside the generic `lt.exe.p`, and
PREFERS LOADING THAT FILE OVER `lt.exe.p` on its next run if both exist**
(confirmed directly by the user, who knows the program's actual behavior).
`sweep.py`'s `setup_cell_dir`/`setup_cell_dir_triangulated` only ever wrote
a fresh donor seed to `lt.exe.p` — they never touched or deleted the stale
per-index checkpoint. That means: **the very first-ever automated attempt
on a cell was the only one that could possibly have used the intended
donor seed.** Every subsequent attempt on that same cell (a retry, a
different donor, N=1 vs N=2, single-donor vs triangulated) would have
silently reloaded whatever `lt.exe.<cell>.dat.p` was left behind by the
FIRST attempt and just continued from there — completely ignoring the
new seed. Confirmed empirically: every one of 5 spot-checked backlog
cells (`kN060_W050`, `kN060_W170`, `kN040_W150`, `kN020_E170`,
`kN020_W050`) already had a stale per-index checkpoint sitting in its
`Sep2026/` directory.

**This almost certainly explains the ~0% success rate of the entire
24-cell backbone-refix backlog across BOTH N=1 and N=2** — the "N=1 vs N=2
processor-count" story below was investigating a difference that may
never have mattered, because neither run was actually re-seeding most of
these cells at all after their first attempt. Treat every "STILL BAD"
verdict in `refix.log` for a cell's *second-or-later* attempt (i.e.
everything except each cell's very first-ever attempt, back in the
original 46-cell sweep) as **not a real test of that donor/N combination**
— the backlog needs to be re-run from here with the fix below in place
before any conclusion about "genuinely hard cells" can be trusted.

**Fix applied** (`sweep.py`): new `clear_stale_checkpoint(target, cell)`
helper, called from `setup_cell_dir`, `setup_cell_dir_triangulated`, AND
`attempt()` (the last one matters just as much — within a single 4-attempt
cascade, DTW+TRUE/DTW+FALSE/CC+TRUE/CC+FALSE, each attempt after the first
would have silently continued from the PREVIOUS attempt's own ending
checkpoint rather than being an independent fresh trial from the same
donor seed; now every attempt in the cascade genuinely starts from the
same `lt.exe.p`). No cleanup script was needed for the existing stale
files scattered across the backlog — the fix deletes each cell's stale
checkpoint automatically the next time that cell is attempted.

**The manually-fixed gold-standard cells this session are NOT affected**
by this bug — the user was updating `lt.exe.p` and its matching per-index
checkpoint together in the same `lte_gui.py` run each time (confirmed by
inspecting both files' contents: identical core parameters, just written
moments apart), so there was never a stale/mismatched pair for any of
those. This bug is specific to `sweep.py`'s own cross-directory donor
seeding.

## THE ONE THING TO FIX FIRST (ORIGINAL, now superseded by the above)

`sweep.py`'s `BASE_OVERRIDES["NUMBER_OF_PROCESSORS"]` was just reverted from
`"1"` back to `"2"` (line ~65-82). It had been dropped to 1 out of caution
about a memory-kill issue that turned out to be caused by something else
entirely (`ulimit -s unlimited` — see below, already fixed separately). A
24-cell backbone-refix batch run at N=1 had a **near-0% success rate** (9
cells processed, 0 refixed) versus every earlier N=2 batch succeeding on the
first attempt most of the time. Single-threaded random-descent search likely
just can't escape a bad local optimum within `TIMEOUT` as reliably as 2+
parallel threads exploring independently. **Re-run everything below at N=2**
before concluding a cell is genuinely hard to fit.

**UPDATE (later same day, 2026-09-27 ~16:00-17:30):** N=2 re-run is in
progress (see "State of the 24-cell backbone-invalidation refix" below) —
so far 3/24 processed, all 3 still failing (`kN040_W170`, `kN020_W030`,
and `kN060_W050` in progress). N=2 is *not* turning out to be the clean
fix N=1-vs-N=2 evidence originally suggested; these may just be genuinely
harder cells (see `kN040_W170` note just below). Don't over-index on the
"N=1 vs N=2" framing from earlier in the day — treat per-cell difficulty
as the more likely explanation until more data comes in.

**BUG FOUND (2026-09-27 ~15:55): `invalidate_and_refix.py` cannot be
safely re-run as-is.** Its `find_bad_cells()` only catches cells whose
*latest ledger entry* is still `accepted: true` with a bad backbone. Once
a cell has already been invalidated (latest entry `accepted: false`,
reason `"invalidated-bad-backbone"` or a subsequent failed refix attempt
`"unsolved-needs-review"`), a bare re-run of that script finds **zero**
bad cells and does nothing — contrary to this file's own earlier claim
that "it's safe to just re-run it directly." Confirmed directly: running
it fresh logged `invalidating 0 bad-backbone cells: []` and exited
immediately claiming "25 total solved, 0 still bad" even though the
24-cell backlog was untouched. **Fix:** use the new `resume_refix.py`
instead (same flood-fill loop, but resumes from `accepted: false` cells
whose reason is `"invalidated-bad-backbone"` or `"unsolved-needs-review"`,
recomputed fresh from the ledger each run — safe to stop/restart). This is
what's currently running in the background (see below).

## Key fixes made this session (all already in `sweep.py`)

1. **Cross-date seeding via Bessel-inversion** (`derive_seed_init`,
   `setup_cell_dir`). A donor's raw `init` value means "value at the
   donor's own record start date" — copying it verbatim into a target
   whose record starts on a different date is wrong, sometimes
   catastrophically (confirmed: sent dLOD from +0.99 to **-0.48** on the
   exact same cell/shape-params, just because of this). Fix: read the
   donor's own `lte_results.csv` at the target's start date (column 4,
   post-Bessel Forcing) and invert the Bessel nonlinearity
   (`invert_bessel`) to recover the correct pre-Bessel seed value. Only
   skipped when donor and target start dates already match (~0.5yr
   tolerance) — then it's a safe verbatim copy (the established
   `kS040_W050`→`kS040_W050_` pattern).

2. **DTW-first, not CC-first** (`_run_cascade`'s attempt order: DTW+VALIDATE
   TRUE → DTW+FALSE → CC+TRUE → CC+FALSE). CC repeatedly overfits train/
   validate while collapsing on the true held-out test window (seen
   directly: CC test scores of 0.26, 0.155, -1.8 despite fine train/
   validate); DTW stays consistent across all three splits. This reversed
   the *original* design's CC-first assumption, based on direct evidence.

3. **Bounded stack ulimit, not "unlimited"** (`write_lte_run_sh`, `ulimit -s
   65536` not `unlimited`). Every attempt using `unlimited` (the project's
   own long-standing convention in every hand-written `lte_run.sh`) got
   killed by the harness's own memory-safety monitor within seconds, with
   zero output — even though actual RSS usage was a modest ~50-60MB. A
   side-by-side test with a bounded 64MB stack ran fine in the identical
   system state. Do not revert this back to "unlimited".

4. **Backbone sanity gate** (`backbone_ok`, `BACKBONE_MIN=0.1`,
   `BACKBONE_MAX=2.5`) — **the most important and most recently-added
   fix, not yet fully exercised**. dLOD and CC/DTW scores alone do NOT
   validate that the fitted winding frequency (`ltep[nm-1]`) is physically
   sensible. A real optimizer run, given the resp's wide `SPREAD_MAX=0.45`
   perturbation range, can drift a perfectly good seed to a spurious
   high-frequency alias (values like 18.68, 100, 150, 200, 600, 678
   recurred *identically* across unrelated cells) that still clears the
   dLOD/CC gates by coincidence. Discovered by comparing the user's
   manually-tuned `kS020_E050` (backbone 0.207, correct) against its
   *automated* neighbors (backbone 18.68, wrong) — **24 of the 46 cells
   accepted before this fix (52%!) had a nonsensical backbone** and were
   invalidated. `backbone_ok(cell_dir)` is now ANDed into every acceptance
   check in `_run_cascade`. Every fit accepted going forward is
   backbone-sane by construction; anything accepted *before* this fix was
   added should be treated with suspicion unless it's been re-verified.

5. **`latest_entry_per_cell()` / ledger correctness** (`solved_cells`,
   `current_solved_nodes` both use this now). Needed so an explicit
   invalidation entry (`accepted: false, reason: "invalidated-bad-backbone"`)
   correctly overrides a stale earlier `accepted: true` entry for the same
   cell, instead of the donor pool still trusting it.

## Manually-verified "gold standard" cells (NOT sweep.py outputs — live in Feb2026/)

The user hand-tuned two cells directly in their `Feb2026/<cell>/` working
directories (visible via the many `lt.exe.p.MM-DD-2026-N-...` snapshot
backups there). Both independently verified this session via `TEST_ONLY`
through `watch_run.exp` (not just trusting the live directory's own claim),
and both have sane backbones:

- **`kS020_E050`** (Madagascar / Mozambique Channel box): CC, VALIDATE=FALSE,
  F9=0, UNCOMPENSATED=TRUE. train=0.844, validate=0.708, test=0.655,
  dLOD=0.99513, backbone=0.2072. Genuinely excellent, no caveats.
- **`kN040_W070`** (Georges Bank / Gulf of Maine mouth): CC, VALIDATE=FALSE,
  F9=0, UNCOMPENSATED=TRUE. train=0.824, validate=0.803, **test=0.205**
  (the specific 2000-2005 window is notably weaker than the rest of the
  record — not disqualifying, but worth knowing), dLOD=0.99589,
  backbone=0.2066.
- **`kN040_W170`** (North Pacific, ~40N/170W): added 2026-09-27 ~17:20.
  This was one of the 24 backbone-invalidated cells — the automated search
  (both N=1 and N=2 attempts) drove it into a stiff, high-frequency-alias
  state and never accepted it (the `backbone_ok` gate correctly rejected
  it both times; ledger integrity was never compromised, only the scratch
  `Sep2026/kN040_W170/lt.exe.p` was left in a bad state as a side effect
  of repeated failed attempts). The user manually re-fit it directly (via
  `lte_gui.py`) and reset it as a gold standard. CC, VALIDATE=TRUE, F9=1.
  train=0.784, validate=0.802, test=0.533, dLOD=0.99443 (note: close to
  the 0.994 floor, not much margin), backbone=0.2073 — matching the
  shared ~0.207 tidal backbone almost exactly, same as `kS020_E050`'s
  0.2072. Independently verified this session via `TEST_ONLY` through
  `watch_run.exp` (train/validate/test/dLOD above are from that
  verification run, not just trusted from the live directory). Ledger
  entry recorded with `attempt_tag: "manual-gold-standard"`, superseding
  the earlier automated `"unsolved-needs-review"` entry. Unlike the other
  two, `Sep2026/kN040_W170/` was *also* already synced to match (backbone
  identical to 15 decimal places) — no separate donor-dir special-casing
  needed for this one, but **verify this before trusting any future
  manually-fixed cell as a donor**: `donor_dir_for()` / the inline
  duplicate of that logic in `invalidate_and_refix.py`/`resume_refix.py`
  only resolves to `FEB /` for cells in `FLAGSHIPS` — anything else
  resolves to `HERE /` (Sep2026) even if the real fit lives in Feb2026,
  so a manually-fixed grid cell's `Sep2026/` copy must be kept in sync or
  donor lookups will silently pull a stale/wrong fit.
- **`kN000_W170`** and **`kN020_W170`** (added 2026-09-27 ~17:50): two more
  cells stacked directly below `kN040_W170` at the same longitude (W170),
  forming a N000/N020/N040 @ W170 gold-standard column. Neither had ever
  been touched by `sweep.py` before (no prior ledger entry either way).
  Both independently verified this session via `TEST_ONLY` through
  `watch_run.exp`, both already synced Feb2026↔Sep2026 (identical
  backbones), both recorded with `attempt_tag: "manual-gold-standard"`.
  - `kN000_W170`: CC, VALIDATE=TRUE, F9=1. train=0.712, validate=0.796,
    test=0.579, dLOD=0.99475, backbone=0.2073. **User notes this cell's
    behavior is very close to NINO4** — worth following up on (same
    longitude band as the equatorial Pacific, N000 = equator).
  - `kN020_W170`: CC, VALIDATE=TRUE, F9=1. train=0.835, validate=0.807,
    **test=0.280** (same weak-test-window caveat pattern as `kN040_W070`'s
    0.205 — not disqualifying), dLOD=0.99424, backbone=0.2073.

  This column gives the automated `resume_refix.py` backlog fresh nearby
  donors for `kN060_W170` (directly above, in the backlog), `kN020_W150`,
  and `kN020_W130` (same latitude band) — `current_solved_nodes()` is
  recomputed fresh every loop iteration so these take effect automatically
  without restarting the background job.
- **`kN060_E170`** (added 2026-09-27 ~20:00): DTW, VALIDATE not set (pair,
  not triplet). pair=(0.701, 0.648), dLOD=0.99542, backbone=0.2069.
  Never previously touched by sweep.py at all — **no `Sep2026/` directory
  existed for it**, so I created one and copied `lt.exe.p`/`lt.exe.resp`/
  `lte_run.sh`/`kN060_E170.dat`/`lte_results.csv`/`dlod_ref.dat` over from
  Feb2026, otherwise a future donor lookup for this cell would have hit a
  `FileNotFoundError` in `setup_cell_dir`/`derive_seed_init` (both read
  `donor_dir/{cell}.dat`, `lt.exe.p`, and `lte_results.csv`). **Not** in
  the 22-cell automated backlog, so no restart of the background job was
  needed for this one. It went on to seed a genuine automated success:
  `kN060_W170` REFIXED on its very first post-checkpoint-bugfix attempt
  (DTW, VALIDATE=TRUE, train=0.690, validate=0.805, test=0.412,
  dLOD=0.99545, backbone=0.2069) — the first hard evidence the checkpoint
  fix actually works.
- **`kN020_W030`** (manual-user-update, ~20:10): pair=(0.017, 0.791),
  dLOD=0.99623, backbone=0.2073. **Caveat: one pair value (0.017) is very
  weak** compared to every other manual fix this session — flagged in the
  ledger note, not treated as a clean gold standard.
- **`kN060_W050`** (manual-user-update, ~21:05): pair=(0.235, 0.824),
  dLOD=0.99610, backbone=0.2080. Solid result, no caveats. Automated
  search (both pre- and post-checkpoint-fix) had failed on this cell
  every time — genuinely needed the manual fit.
- **`kS020_E070`** and **`kS040_E070`** (~22:55-23:00): a real manual
  *donor chain*, per the user: `kS020_E050` (Madagascar gold standard) →
  `kS020_E070` → `kS040_E070`, each seeded from the previous. Both CC,
  VALIDATE=TRUE, F9=1.
  - `kS020_E070`: train=0.849, validate=0.736, **test=0.204** (weak,
    same pattern as `kN040_W070`'s 0.205), dLOD=0.99489 (passes floor
    with very little margin, ~0.0009). Never previously touched by
    sweep.py — created `Sep2026/kS020_E070/` and synced files over
    (same reason as `kN060_E170`).
  - `kS040_E070`: train=0.740, validate=0.747, **test=0.264** (also
    weak), dLOD=0.99494 (also thin margin above floor). Supersedes the
    automated backlog's `unsolved-needs-review` entry (which had tried
    triangulating from `kS040_E090`/`kS020_E050`/`iodw`, not from the
    directly-north `kS020_E070`. Already excluded from the running job
    (it had already failed and moved past this cell before the manual
    fix landed), so no restart was needed for either of these two.
- **`kS040_E050`** (~23:08): fourth link in the same donor chain:
  `kS020_E050` → `kS020_E070` → `kS040_E070` → `kS040_E050`, each seeded
  from the previous ("like dominoes", per the user) — a real manual
  flood-fill, propagating the shared backbone across the whole
  Mozambique-Channel/Madagascar-shelf region one neighbor at a time. DTW,
  VALIDATE=TRUE. train=0.701, validate=0.797, **test=0.792** — strong
  across all three splits, no weak-test caveat unlike the two links
  before it. dLOD=0.99477, backbone=0.2071. **This one WAS still in the
  running job's live backlog** (unlike the previous two), so it needed
  the full kill-and-restart treatment — done immediately, backlog now
  down to 13 cells.

I've recorded ledger entries for both (`attempt_tag: "manual-gold-standard"`)
so `latest_entry_per_cell()` reflects them correctly. **They still physically
live in `Feb2026/`, not `Sep2026/`** — if you want to use either as a donor
for further seeding, pass `sweep.FEB / "kS020_E050"` (etc.) as `donor_dir`
directly, not `sweep.HERE / "kS020_E050"`.

### The emerging physical pattern (user's own hypothesis, worth continuing to test)

Three "gold standard" locations found so far — Patagonian shelf (earlier
project history, not this session), the Mozambique Channel/Madagascar shelf
break, and Georges Bank/Gulf of Maine — are all places where shelf or channel
geometry is known to resonantly amplify tidal currents. The hypothesis: these
locations show an unusually clean, strong expression of the shared ~0.207
tidal-forcing backbone (the project's central finding across nino4, pdo,
etc.), and are *also* where the automated search is most prone to getting
distracted by a spurious alias instead of finding the true signal — i.e. shelf/
channel sites are simultaneously the best places to see the real effect and
the hardest for blind search to find it. Worth deliberately testing more
shelf/channel locations if any are identified.

**Fourth data point, 2026-09-28: `kS040_E170`** (near New Zealand). The user
independently drew the same connection unprompted: this cell shows "very
strong continuous windings" and is physically comparable to the Madagascar
site — a strong tidal-current *differential* on both sides of New Zealand,
same shelf/channel-resonance mechanism as the Mozambique Channel. It also
fits the second half of the hypothesis: automated search failed on it every
single time (both pre- and post-checkpoint-bugfix), and it needed a manual
fit (CC, VALIDATE=TRUE: train=0.804, validate=0.781, test=0.164 — weak test
again, dLOD=0.99511, backbone=0.2072). Worth actively looking for more
strait/channel/shelf-break locations around other island arcs and
continental margins as deliberate candidates, rather than waiting for them
to turn up in the generic grid sweep.

## Cross-testing the Madagascar donor onto its immediate neighbors

- **`kN000_E050`** (Somalia coast, directly north of Madagascar): seeded from
  `Feb2026/kS020_E050` → **SOLVED**, first attempt (DTW, VALIDATE=TRUE).
  train=0.592, validate=0.749, test=0.598, dLOD=0.99513, backbone=-0.210
  (sane). Clean success — much better than its earlier `iodw`-seeded
  automated attempt, which failed.
- **`kS020_E030`** (Mozambique Channel proper, directly west of Madagascar):
  seeded from `Feb2026/kS020_E050` → **FAILED**, all 4 attempts (both
  retries), at N=1. Notably, dLOD was actually fine every time (~0.9951,
  above the 0.994 floor) — so this was rejected on `backbone_ok`, `stiff()`,
  or a triplet/pair parsing miss, not the dLOD gate. **Re-run this one first
  at N=2** — given it already looked close (good dLOD every time), it may
  just need the extra search breadth N=2 provides. Check
  `sweep.HERE / "kS020_E030" / "lt.exe.p"`'s backbone after re-running to see
  which gate actually rejected it.

## State of the 24-cell backbone-invalidation refix (UPDATE: now running at N=2 via `resume_refix.py`)

All 24 originally-bad-backbone cells were invalidated in the ledger. The
original refix driver (`invalidate_and_refix.py`) got through 9 of them at
N=1 (all 9 failed — see above). **`invalidate_and_refix.py` itself cannot
be bare-re-run to resume this** (see the bug writeup near the top of this
file) — use `resume_refix.py` instead, which correctly re-derives the
still-bad list from the ledger (`accepted: false` with reason
`"invalidated-bad-backbone"` or `"unsolved-needs-review"`) every time it's
run, so it's safe to stop/restart at any point.

**`kN040_W170` and `kN000_W030` have been removed from this backlog** —
`kN040_W170` because the user manually re-fit it (`manual-gold-standard`
ledger entry); `kN000_W030` because the user manually updated it too
(`manual-user-update` ledger entry — passes dLOD/backbone gates cleanly,
though it reports a VALIDATE=FALSE train/test *pair* rather than a
triplet, and the user separately flagged its raw 1950-1965 window as
suspiciously quiet — investigated and traced to a real, non-corrupted
low-variance stretch in the source Kaplan SST data itself, not a fit
bug). `kN000_W030` was re-verified a second time after the user made a
further small tweak and explicitly synced it to `Sep2026/kN000_W030`
themselves (confirmed identical backbones both places): final numbers
CC pair (0.259, 0.817), dLOD=0.99605, backbone=0.207218. As of the
checkpoint-bug fix (below), this list has been changing fast — **don't
trust the snapshot below as current**; run
`python3 -c "import sys; sys.path.insert(0,'.'); import resume_refix;
print(resume_refix.find_still_bad_cells())"` from `Sep2026/` for the
live truth. Snapshot at ~19:25 (22 cells, before the checkpoint-bug fix):

`kN020_W030`, `kN060_W050`, `kN060_W170`, `kN040_W150`,
`kS040_E170`, `kN040_W010`, `kS040_E070`, `kN020_E170`, `kN020_W050`,
`kN020_W150`, `kN020_W130`, `kS020_E170`, `kN020_E090`, `kN020_E070`,
`kN020_E110`, `kS040_E150`, `kS020_E150`, `kN000_E150`, `kS040_E050`,
`kS020_E030`, `kS020_E010`, `kN040_W090`

Since then: `kN020_W030` and `kN060_W050` also got manual fixes (see
gold-standard section above); `kN060_W170` was genuinely REFIXED by the
automated pipeline post-bugfix (the first confirmed automated success
under corrected conditions). Roughly **~18 cells** likely remain as of
this writing but check live.

**LESSON LEARNED (2026-09-28 ~00:40):** when the user says "cell X
updated," verify/record/sync it IMMEDIATELY, before doing anything else
(including answering a side question they ask in the same breath) — I
got pulled into investigating a mechanistic question about `kN040_W010`
right after the user flagged it as updated, and by the time I came back
to actually record it, the background job had already re-attempted
`kN040_W010` from scratch with its own donor seed (ignoring the user's
fix entirely, since nothing in the ledger yet marked it as done) and
failed again, overwriting `Sep2026/kN040_W010/` with its own bad
leftover state. Recovered fine (re-synced from Feb2026, no data lost),
but the wasted ~20min automated attempt and confusion could have been
avoided. `kN040_W010`'s final recorded scores (manual-user-update):
CC, VALIDATE=TRUE, train=0.833, validate=0.684, **test=-0.1995
(NEGATIVE — worse than the usual "weak" caveat pattern)**, dLOD=0.99591,
backbone=0.2073.

**IMPORTANT — restart discipline:** whenever a cell in this backlog gets
a manual fix (new `accepted: true` ledger entry) *while the background
driver is already running*, you must kill and restart `resume_refix.py`
before it reaches that cell — its to-do list (`bad_cells`/`remaining`) is
computed once at startup from the ledger, so a mid-run ledger update
won't remove an already-queued cell, and the driver will happily
overwrite the manual fix in `Sep2026/<cell>/` if it gets there first
(confirmed this needed doing twice this session: once for the `pna`
removal, once for `kN000_W030`). No ledger corruption risk from
killing it — `solve_cell_with_retry`/`solve_cell_triangulated_with_retry`
only call `append_ledger` after a cell's full cascade finishes, so a
`kill -9` mid-cell just loses that one cell's in-progress attempt, not
any ledger data.

**Live status as of 2026-09-27 ~19:25** (this session): `resume_refix.py`
is running in the background (restarted twice, most recently ~19:25 —
PID varies, check `pgrep -f resume_refix.py`), processing this backlog
nearest-first via the same flood-fill donor logic, at N=2 and now without
`pna` in the donor pool, ~17-35 min/cell (faster if the single-donor
attempt succeeds, slower if it falls through to the triangulated k=3
fallback). Progress so far (before the two restarts above; unaffected by
them): `kN040_W170` (now excluded, see above), `kN020_W030`,
`kN060_W050`, `kN060_W170`, `kN000_W030` (now excluded, see above) all
still bad. **So far N=2 is not showing the clean recovery the N=1-vs-N=2
comparison predicted** — treat that framing with suspicion; these may
just be genuinely hard cells needing manual fits like the ones above.
Log: `Sep2026/refix.log` (append-only, safe to `tail -f`). Completion
marker: `Sep2026/RESUME_REFIX_DONE` (written only when the whole backlog
is processed — check `find_still_bad_cells()`'s output in the log's
final `=== resume-refix complete ===` line for the true up-to-date
remaining list, since more cells may get fixed by later flood-fill
donors even after early failures). A `Monitor` watching `refix.log` for
`REFIXED|STILL BAD|EXCEPTION|resume-refix complete` lines (see
`watch_refix.sh`) is the recommended way to track it rather than
polling.

## Currently solved (23 cells, all backbone-sane, all from `sweep.py`'s own pipeline)

`kN000_E050, kN000_E090, kN000_E110, kN000_W150, kN020_E050, kN020_E150,
kN040_E150, kN040_W030, kN040_W050, kN060_E010, kN060_W010, kN060_W030,
kS020_E110, kS020_W090, kS020_W110, kS020_W130, kS040_E090, kS040_W050,
kS040_W130, kS040_W150, kS040_W170, kS060_W050, kS060_W070`

Plus the 2 manually-verified gold standards (`kS020_E050`, `kN040_W070`,
recorded separately, see above) = **25 genuinely trustworthy solved cells**
as of this handoff.

One remaining known-suspicious case: `kS020_W090` was accepted twice (once
via single-donor, once via a self-excluded triangulated re-check) with
near-identical pathological scores both times (`validate=0.0` exactly,
`test≈-1.8`) despite a sane backbone and passing dLOD. This wasn't caused by
the backbone-alias bug (backbone is fine) — worth investigating the cell's
own raw data directly (possible degenerate/near-constant segment landing in
the VALIDATE window) before trusting its CC/validate/test numbers.

## Flagships / donor pool (`sweep.py`'s `FLAGSHIP_COORDS`)

Original 6: `nino4, amo, pdo, pna, baltic, sam` (pna/baltic/sam coordinates
are best-guesses, no formal reference). Added this session: `nino34, nao,
tna, iode, iodw` — all already mature, independently-fitted indices (12-92
backup snapshots each in their own `Feb2026/` dirs), used purely as
additional donor nodes, never re-fit themselves. `iode`/`iodw` gave the
Indian Ocean basin its first donor coverage at all (previously 5,000-8,800km
from any flagship; now 1,200-4,000km) and enabled 8 solved cells there.

**`pna` REMOVED 2026-09-27 ~19:20** — the user flagged it as a bad donor
for Pacific cells (it produced the seed for `kN060_W170`'s and
`kN040_W150`'s failed refix attempts, and it was always just a
best-guess coordinate, `(50.0, -165.0)`, with no real anchor). Real
verified Pacific donors exist now instead: `kN000_W170`, `kN020_W170`,
`kN040_W170` (see gold-standard section above). Removed from both
`FLAGSHIPS` and `FLAGSHIP_COORDS` in `sweep.py`. **Important:** the
already-running `resume_refix.py` background process had the old
`sweep` module cached in memory (Python doesn't hot-reload), so editing
the file alone didn't change its behavior — it had to be killed and
restarted to pick up the fix (no ledger data was lost; the interrupted
attempt hadn't logged anything yet). If you edit `sweep.py` again while
a sweep/refix driver is running in the background, remember to restart
it too.

## `lte_forward.py` bitrot fixed + real Ada Forcing-save bug found (2026-09-29)

Per the user's instruction to fix `experiments/Feb2026/lte_forward.py`'s
bitrot, verifying the manifold (`lte_results.csv` column 4, "Forcing")
BEFORE the Model (column 2), three real bugs were found and fixed, plus one
genuine pre-existing bug in the live Ada `enso_opt`/`lt.exe` binary itself.

**Python bugs fixed in `lte_forward.py`:**
1. `harms` was parsed from the resp's `NH` key (the *initial* harmonic-count
   template, e.g. `"1 1 1 1"`) instead of `lt.exe.p`'s actual fitted `harm`
   array (e.g. `[6.0, 2.0, 9.0, 4.0]`). This duplicated the backbone winding
   several times in the regression design matrix, making it singular — the
   literal `--verify` crash. Fixed to read `params["harm"]`.
2. `STRICT_IDATE`/explicit `INIT_DATE` (which implies it, per
   `gem-lte-primitives-solution.adb` lines ~669-719) were never implemented.
   This activates a different IIR-seeding path (extend the template
   backward, seed at `IDATE - 0.1/sampling` instead of plain `IDATE`) —
   `amo`'s resp sets `INIT_DATE 1880.0` explicitly, silently taking this
   path. Added `extend_backward()` + the branch. Also fixed the non-strict
   path to seed at `dates[0]` (the record's own first date), not `IDATE` —
   the two are conflated in Ada only when they happen to coincide.
3. The IIR backward-reconstruction pass implemented the exact *deprecated*
   heuristic (`y[i-1] = -x[i-1] + Mem*y[i] + copysign(lag_c, y[i])`) that
   the Ada source's own comment says was replaced after
   `iir_invariant_test.adb` caught a ~2.5-unit round-trip error. Replaced
   with the current sign-candidate exact-inverse algorithm.

After these three fixes, a controlled single-threaded (`NUMBER_OF_PROCESSORS
=1`), single-iteration (`Counter=0`, no Markov perturbation yet) live Ada
run — instrumented with temporary debug dumps, since removed — showed
`lte_forward.py`'s forcing (both pre- and post-Bessel) matches the live Ada
`Calc_Forcing`/`Bessel` output to **~1e-13**, i.e. exact to float64 noise.

**Real Ada bug found and fixed** (`gem-lte-primitives-solution.adb`): the
still-remaining ~1.8e-4 relative gap against the *saved* `lte_results.csv`
was not a Python bug at all. Under `VALIDATE=TRUE` (the setting used by
`amo` and most Sep2026 cells), the cross-thread "lockbox" reporting path
resyncs `D` and `Model` to the winning thread's `DKeep`/`KeepModel` right
before saving (`D := DKeep; Model := KeepModel;`, ~line 1968) — but
`Forcing` was never included in that resync. Since `Forcing` is recomputed
unconditionally every iteration regardless of accept/reject, by the time
`Save` fires it can reflect a different (possibly rejected, possibly a
different thread's own trial) candidate than the one whose `Model`/`DKeep`/
`lt.exe.p` is actually being saved — so `lte_results.csv` column 4 (the
manifold) could silently disagree with its own column 2 and its own
`lt.exe.p`. **As the user pointed out**, this is not just a cosmetic
CSV-column mismatch: it means the cross-thread lockbox could report/save
Model/JSON from the winning thread while the accompanying Forcing came from
a non-optimal thread, i.e. the historically-saved parameter sets for
`VALIDATE=TRUE` cells may not be as internally consistent as intended.

Fixed by caching `KeepForcing := Forcing;` at the same accept-time
assignment as `KeepModel`/`DKeep` (`Keep := Set; KeepModel := Model;` ~line
1758), then unconditionally restoring `Forcing := KeepForcing;` whenever
`Ever_Accepted` is true, right before the `Save_Now` decision — covering
both the `VALIDATE=TRUE` lockbox path and the legacy `Best_Client=ID` path
uniformly. Verified on a fresh, real (multi-threaded, `VALIDATE=TRUE`,
15s) `amo` run after rebuilding: the freshly-saved `lte_results.csv`
column 4 now matches `lte_forward.py`'s independent computation of the
same freshly-saved `lt.exe.p` to **1.7e-10 relative** (float64 noise).

**Second Ada bug found and fixed, same root cause:** the same "`M` is
recomputed unconditionally every iteration, never resynced to `DKeep`"
pattern also affects the harmonic-multiplier report block right before
`Save_Windings`/`Shared.Save(DKeep)` (~line 1980): `DKeep.C(I-NM) :=
Integer(M(I)/M(NM))` derives the harmonic multipliers FROM the live
(possibly stale) `M`, and `DKeep.C` is exactly what gets serialized as
`lt.exe.p`'s `"harm"` field via the immediately-following
`Shared.Save(DKeep)`. Checked with the user whether this mattered as much
as the Forcing bug: amplitude/phase (`MAP`) is always freshly
regression-derived and `lt.exe.windings.json` is confirmed genuinely
write-only (grepped the whole source — no `Read_Windings`/`Load_Windings`
function exists anywhere), so that part is low-stakes. But `lt.exe.p`
itself is NOT write-only — it's the checkpoint reloaded by future runs and
by `lte_forward.py`'s own `harm`-array fix (bug 1 above) — so a stale `M`
corrupting `DKeep.C` here would propagate a real inconsistency into the
checkpoint itself. Fixed by rebuilding `M` from `DKeep.B.LT`/`DKeep.C`
directly at the top of the `Save_Now` block (before the report loop that
writes `DKeep.C` back from `M`), making that write-back a harmless
round-trip. Verified: `harm` in a fresh run's saved `lt.exe.p` stays
correctly matched to a fresh `lte_forward.py` replay (manifold match still
1.7e-10 after this second fix, unchanged from before it — confirming no
regression).

**Not yet done:** re-verifying/fixing the Model (column 2) path in
`lte_forward.py` itself — it still needs at minimum the missing
`Annual_Impulse` (`ImpC`) term (added to Model, not Forcing, right after
the LTE response and before the IR delay-differential, in
`gem-lte-primitives-solution.adb` ~line 1567) which `lte_forward.py`
doesn't implement at all. Per the user's explicit ordering instruction,
this was deferred until the manifold match above was confirmed.

## Post-fix audit: how many cells were actually hit (2026-09-29)

Per the user's own recollection ("a VALIDATE=TRUE was stopped at a TEST CC
of 0.6 but when rerun only started at 0.3... that never happened with
VALIDATE=FALSE") — exactly the symptom the Forcing/harm resync bug above
would produce — wrote `audit_validate_bug.py`: for every cell whose resp
sets `VALIDATE=TRUE` (67 of 91; falls back to the Feb2026 companion dir
when a Sep2026 cell has no own `lt.exe.p`/`lte_results.csv` yet), replays
its CURRENTLY SAVED `lt.exe.p` through `lte_forward.py`'s `forward()`
directly (no `lt.exe` launch) and compares the resulting Forcing against
that cell's own saved `lte_results.csv` column 4.

First pass (`audit_validate_bug_results.json`, before a 4th `lte_forward.py`
bug below was found): **30 OK**, **37 SUSPECT** — of which 8 looked severe
(corr 0.82 down to -0.25): `kN000_E150, kN040_W050, kS020_E150,
kS040_E170, kS040_E070, kS040_E050, kN020_E150, kS040_E150`.

**A 4th `lte_forward.py` bug was found mid-re-sweep**, inflating that
first-pass SUSPECT count: `Calc_Forcing`'s non-strict branch actually
seeds `IIR` at `Start => Initial_Conditions_Date` (IDATE) directly, not at
`dates[0]` as the earlier fix wrongly assumed -- a misreading that
confused `Impulse_Amplify`'s own, separate, numerically-irrelevant (since
Offset=Ramp=0 always here) `Start` argument for the IIR seed. This was
silently correct only for cells where IDATE happens to equal `dates[0]`
(true for `amo`, false for e.g. `kN000_E150`, whose `IDATE=1880` but data
starts `1950`) -- so cells with IDATE genuinely predating their own
record were misdiagnosed as SUSPECT by this Python bug, not real Ada
corruption. Fixed by seeding at plain `idate` in the non-strict path too
(the natural `Start_Index` clamping when IDATE precedes the record, not a
separate `dates[0]` seed, is what makes non-strict often *behave* as if
seeded at `dates[0]`). Caught by first running `resweep_corrupted.py`
against the original 8 (2 succeeded and were committed --
`kN000_E150`/`kN040_W050` -- before the residual pattern on
`kN000_E150`'s post-fix audit, ~1e-3 instead of ~1e-13, prompted digging
further), then re-verified against `amo` (still matched, 1.7e-10, once
compared to a FRESH Ada regen rather than the repo's still-pre-fix
`amo/lte_results.csv` snapshot) and `kN000_E150` (now 6e-11) directly.

Re-running the full audit with this fix: **SUSPECT dropped from 37 to
17** -- confirms most of the original "mild" 29-cell tier was this
Python bug's own false positives, not genuinely corrupted Ada saves. Of
the original "severe 8", 6 remain genuinely SUSPECT after the fix:
`kS020_E150, kS040_E170, kS040_E070, kS040_E050, kN020_E150,
kS040_E150` -- these are the real, confirmed casualties of the Ada
Forcing/harm resync bug. `kN040_W050` (IDATE already equalled `dates[0]`,
so immune to the Python bug) was genuinely Ada-corrupted and is now fixed
by its re-sweep; `kN000_E150`'s status is more ambiguous (may have been a
pure Python-bug false positive, or a mix) but its re-sweep produced a
valid, self-consistent, good-scoring result regardless.

**Session lesson:** `pkill -f resweep_corrupted.py` (to abort a batch
mid-run) only kills the Python driver, not the `expect`/`lt.exe`
descendants `subprocess.run` spawned -- found one orphaned `lt.exe -j`
still running under `kS020_E150/` afterward (harmless here, hadn't
written anything yet, killed directly by PID). Kill the process group
(or find and kill descendants explicitly) next time, not just the
driver's own name.

`resweep_corrupted.py`'s `TARGET_CELLS` now holds the 6 confirmed-real
cells (`kN000_E150`/`kN040_W050` removed, already committed). Re-sweep
of the 11 remaining milder SUSPECT cells is the natural follow-up after
these 6 finish.

## Third, more fundamental Ada bug: Forcing missing from the cross-thread lockbox entirely (2026-09-29)

Running the 6-cell re-sweep surfaced a THIRD Ada bug, deeper than the
first two. `kS040_E170` came out of its re-sweep perfectly clean
(post-fix audit rel=2.6e-10, corr=1.0) but `kS020_E150` did not
(rel=0.81, corr=0.77) -- **even though both were re-swept with the
already-`KeepForcing`-fixed binary**. Confirmed this was not
`lte_forward.py` imprecision by re-instrumenting the real Ada binary
(temporary `DUMP_FORCING` hooks, same technique as before) and running a
completely fresh, single-threaded (`NUMBER_OF_PROCESSORS=1`), zero-
perturbation (`Counter=0`) replay of `kS020_E150`'s exact saved
`lt.exe.p`: the pre-Bessel forcing matched `lte_forward.py`'s own
computation to 7e-10 (so both Ada's live replay and the Python port
agree, and are self-consistent) -- but that SAME correct Ada output,
run through the SAME live Ada `Bessel` step with the SAME `lt.exe.p`
parameters, still did not match the file that same overall process had
saved as `lte_results.csv`. Pure Ada vs. its own saved output, zero
Python involved -- proof the bug was still on the Ada side.

Root cause: `Dipole_Model`'s local variables (`D`, `Set`, `M`,
`KeepForcing`, etc.) ARE genuinely task-local (each `Thread` task calls
`Dipole_Model` with its own independent call frame -- not a data race).
But `Monitor.Report_Final`/`Monitor.Winner` (the protected-object
lockbox VALIDATE=TRUE uses to pick the actual best candidate ACROSS all
threads, not just within one) never included `Forcing` in the
handshake at all -- only `D`, `Model`, `M`, `MAP`, `Trend`, `Accel`. So
when `Winner` hands the reporting thread a DIFFERENT thread's winning
`D`/`Model`/`M`/`MAP` (the common case whenever the reporting thread --
whichever happens to call last -- isn't also the actual best thread),
the `Forcing := KeepForcing` fix from earlier in the day still only
resynced within-thread (this thread's OWN last accept), never picking up
the ACTUAL WINNING thread's own `Forcing`. This explains the earlier
"kN040_W050 mild (9.8e-4), kS020_E150 catastrophic (0.81)" split
directly: it depends on whether the reporting thread happened to also be
the true winner (no swap needed, first fix already sufficient) or not
(swap needed, first fix insufficient).

Fixed properly this time: added `Forcing` as an `in`/`out` parameter to
`Monitor.Report_Final`/`Monitor.Winner`, added a
`Best_Validate_Forcing : Final_Model_P` field to the protected object's
private lockbox state (reallocated in `Report_Final` alongside
`Best_Validate_Model`, nulled in `Reset`, returned in `Winner`), and
updated both call sites (the normal per-thread report, and the
exception-handler safety-net report) to pass `Forcing => KeepForcing`.

**Verified 2026-09-29** (after a mid-session Bash-classifier outage
blocked the rebuild for a while -- resolved by switching to manual
approval mode): rebuilt cleanly, re-swept `kS020_E150` a third time
(`sweep._run_cascade`, 408s, `validate_true_dtw`, train/validate/test =
0.646/0.757/0.678) and re-audited -- **rel=9.3e-10, corr=1.0** (was
rel=0.81, corr=0.77 after the first two fixes alone). Appended a new
ledger entry noting it supersedes the earlier, still-inconsistent commit
for this cell. Removed the temporary `DUMP_FORCING` debug hooks
afterward and rebuilt clean again (both hooks were gated behind an env
var and inert for every real sweep run in the meantime, so no real runs
were affected by their presence).

**All 4 remaining cells re-swept and verified 2026-09-29**
(`kS040_E070, kS040_E050, kN020_E150, kS040_E150`), all clean
(rel 5.5e-11 to 1.2e-10, corr=1.0 on every one). All 8 of the original
"severe 8" cells are now fully fixed and self-consistent.

**Final full 67-cell audit** (`audit_v3.log`) with all three Ada fixes
in place: **56 OK** (up from 30), **11 SUSPECT** (down from 37, then
17) -- and critically, none of the remaining 11 are severe: worst is
`kN040_W030` at corr=0.994, the rest are corr>=0.997, six of them
corr>=0.9999. This confirms the hypothesis noted above -- the earlier
17-cell count WAS partly inflated by the cross-thread gap (9 of those 17
are now clean) -- and shows the three fixes together account for the
overwhelming majority of real corruption. Remaining 11:
`kN040_W030, kN060_E010, kN020_W130, kS060_W050, kS020_E010,
kN040_W090, kN020_W090, kN060_W130, kN020_E170, kN000_W170, kN040_W170`.
None were in the original "severe 8"; re-sweeping them is a natural,
lower-urgency follow-up (same `resweep_corrupted.py` pattern, new
`TARGET_CELLS` list) whenever picked back up -- no further Ada-side
investigation should be needed, this looks like ordinary remaining
cross-thread-lockbox cases the fix already handles correctly, just not
yet re-run.

**All 11 re-swept 2026-09-29 (user going to bed, ran unattended):** 10 of
11 committed clean (rel ~1e-10, corr=1.0 on every one):
`kN040_W030, kN020_W130, kS060_W050, kS020_E010, kN040_W090, kN020_W090,
kN060_W130, kN020_E170, kN000_W170, kN040_W170`. One,
`kN060_E010`, failed its full 4-attempt cascade (1628s, deadlock/no
gate pass) and was rolled back to its pre-batch state unchanged --
still SUSPECT (corr=0.998, mild) but not regressed.

**Final full 67-cell audit (`audit_v4.log`): 66 OK, 1 SUSPECT
(`kN060_E010` only).** Started this investigation at 30 OK/37 SUSPECT.
All three Ada bugs (KeepForcing, DKeep.C/harm, and the
Report_Final/Winner cross-thread Forcing gap) plus all four
`lte_forward.py` bugs are confirmed fixed project-wide. `kN060_E010` is
the only remaining open item -- try a longer TIMEOUT or a donor-seeded
reset (`solve_cell_with_retry`, discarding its own possibly-poisoned
starting state entirely) rather than another in-place resume next time,
since in-place resume has now failed for it once already.

## ImpC pilot + cheap Feb2026-vs-Sep2026 best-of pass (2026-09-29)

Per user direction ("what I want to see is progress... spend only a
minute on each"), pivoted away from expensive re-optimization sweeps
toward two cheap, mostly-search-free improvements:

1. **`lte_forward.py` gained the missing Annual_Impulse (ImpC) term**
   (added directly to Model, not Forcing, right after the LTE response).
   Confirmed via `gem-random_descent.adb` that ImpC has been stuck at
   0.0 in every prior fit project-wide: `Markov`'s `Set(I)=0.0` branch
   never perturbs a parameter sitting at exactly zero (recurses to pick
   a different one instead, unless `FLIP<0.0`, which no resp sets) -- a
   structural dead zone, not evidence the term is useless. Piloted
   seeding `impC=0.01` (`seed_impc.py`, `pilot_impc.py`) on 3 diverse
   cells: all three improved, with `impC` growing meaningfully away
   from the seed on the biggest gainer (`kN060_W050`: 0.235->0.666,
   impC 0.01->0.061) -- real signal. A full 86-cell seeding sweep was
   started then explicitly stopped per user redirect toward the faster
   approach below; re-seeding is still available as a lower-priority
   follow-up (`sweep_impc_full.py`, TARGET_CELLS already excludes the
   3 piloted cells).
   - Also found/fixed: the first seeding attempt wrote compact single-
     line JSON via `json.dumps`, which Ada's JSON reader rejected
     ("2:9: comma expected"), corrupting two cells' checkpoints
     (recovered via `git checkout`). Fixed by editing the file as text
     (single targeted regex substitution), preserving Ada's own
     multi-line pretty-printed format exactly.
   - Also found/fixed a real, unrelated gap while running the sweep: 29
     Sep2026 cells had a checkpoint but no local `.dat` file at all
     (silently missing, inherited from Feb2026 companion dirs) --
     copied from Feb2026.

2. **`pick_best_of_feb_sep.py`**: a genuinely cheap, no-search pass --
   for every one of the 89 cells, compares Feb2026 vs Sep2026
   `lte_results.csv` directly (Pearson CC between Model and Data
   columns, zero Ada invocations), and where Feb2026 is meaningfully
   better, gates the import on manifold/backbone range (instant) and a
   dLOD reading via Ada's native `TEST_ONLY=TRUE` mode (evaluate the
   loaded parameters exactly once, no search, ~4s). Result: 43
   imported, 44 kept Sep2026, 2 initially rejected on the dLOD gate.
   Status map: 77 good/11 warning/1 serious -> 85 good/4 warning/0
   serious/0 critical.

   **Important bug found in `TEST_ONLY` itself, not just this script**:
   `TEST_ONLY=TRUE` is NOT a pure passive "evaluate as-loaded, don't
   mutate anything" mode, contrary to what its use here assumed. The
   harmonic-collision-avoidance code (`Walker.Force_Harmonic`,
   gem-lte-primitives-solution.adb ~line 1440) runs BEFORE the
   `if Test_Only then exit;` check (~line 1820), in the same (only)
   iteration -- so if the loaded harmonics happen to collide, Force_
   Harmonic silently redraws a NEW RANDOM harmonic value before saving,
   even under TEST_ONLY. Found by re-investigating a user report that
   `kN060_E010` "was never optimized earlier in Feb2026" (the user was
   concurrently, manually re-fitting it in Feb2026 -- 6 snapshot saves
   08:42-08:54 -- while this pass ran, initially causing a race); a
   direct re-check showed the SAVED `lt.exe.p`'s windings (`ltep`,
   `harm`) had changed between two successive `TEST_ONLY` calls on the
   identical loaded state, non-reproducibly. Cross-checked the other 41
   originally-"imported" cells' current `lte_results.csv` CC against
   the CC recorded at import time: **9 of 43 had drifted**, some
   severely (`kN020_E050` 0.726->0.304, `kN060_W010` 0.612->0.156) --
   real corruption from this same mechanism, not noise.

   `dLOD` itself is unaffected (computed from `Forcing`, via
   `Calc_Forcing`, BEFORE the harmonic-collision code runs each
   iteration) -- so the dLOD readings already gated on remain valid.
   Fixed by re-copying `lt.exe.p`/`lte_results.csv`/`lt.exe.windings.json`
   fresh from Feb2026 for all 43 imported cells (no further Ada
   invocation), restoring exact consistency; verified all 43 now match
   their recorded import-time CC exactly. `kN060_E010` and
   `kS040_W170` (the 2 originally dLOD-rejected cells) were re-added
   using dLOD sourced from their own real (non-TEST_ONLY) historical
   fit logs (`cc_refit.log`) instead of a fresh `TEST_ONLY` re-check.

   **Do not reuse `TEST_ONLY` for a "peek without touching anything"
   check without also gating on whether the loaded harmonics might
   collide** -- it is a check-and-possibly-mutate operation, not a pure
   read. A real fix would move the `Test_Only` exit check before the
   collision-avoidance block, or add a `Test_Only`-aware guard around
   `Force_Harmonic`'s call inside it; not done this session.

## Subharmonic quadrature metastability -- a real, unresolved structural degeneracy (2026-09-29)

User finding: the LOWEST-value winding (~0.01-0.02, the `BACKBONE/11`
subharmonic family responsible for the 60-120yr AMO-type frequency-
doubling modulation) is genuinely METASTABLE. Its amplitude and phase
are decided purely by `Regression_Factors`' OLS fit against
`[sin(k*F), cos(k*F)]` basis columns, per-quad, with no physical
constraint steering the outcome. Because this winding's period is
comparable to or longer than the ~140yr record, `sin(k*F)` and
`cos(k*F)` are nearly collinear over the observed range -- the
regression can land on either a "cosine-dominant" or "sine-dominant"
solution (a ~90 degree quadrature swap, not a simple 180 degree sign
flip) with no way to algebraically choose between them after the fact.
Per the user: "there is no way to automatically invert the winding
other than determining which fits the data better" (the "Mach-Zehnder
encryption paradigm" analogy) -- resolving it per-cell requires a real
A/B optimization (lock each quadrature, refit, compare), not a cheap
linear-algebra trick. **Relevant to long-term-trend interpretation**:
since this winding's period rivals the record length, an arbitrarily-
selected phase can partially alias into the fitted trend/level terms,
muddying whether an individual quad's apparent SST trend is physical or
partly a fitting artifact of this degeneracy.

**Direct evidence gathered** (`experiments/Feb2026/amo` vs
`experiments/Sep2026/amo`, two independent fits of the identical AMO
data -- Feb2026 is described by the user as the "cosine" branch,
Sep2026 as the newer "sine" branch):
- Fit quality nearly tied: CC(model,data) = 0.776 (Feb2026) vs 0.782
  (Sep2026).
- Sep2026 slightly more parsimonious: `impC=0.0` (Annual_Impulse
  inactive) vs Feb2026's `impC=-0.028` (active) -- same fit quality
  with one fewer working mechanism.
- The two full manifolds show a real but PARTIAL mirror-image
  relationship: `corr(Feb2026, -Sep2026)` on the raw post-Bessel
  forcing is ~0.26, stable across every smoothing scale from 5 to 40
  years (not noise, but far from a clean identity). Backbone values
  also differ genuinely (0.229 vs 0.174), and harmonic sets are
  entirely different (`[6,2,9,4]` vs `[2,5,12,4]`). Conclusion: these
  are two genuinely DIFFERENT local optima, not the same solution in
  two equivalent conventions -- though both are plausibly anchored to
  the same real AMO regime-shift timing (see `project_amo_resolved_
  shallow_water.md` in memory).

**A population-wide coherence test was attempted and is INCONCLUSIVE --
do not re-run this exact approach expecting a different answer**
(`check_subharmonic_coherence.py`): tested whether the subharmonic's
fitted (amp*cos(phase), amp*sin(phase)) vector is spatially LESS smooth
across the 89 quads than faster "control" windings (k=0.207, 0.414,
0.828), and whether a sign-canonicalization improves that. Result:
noisy in both directions for every winding tested including controls
(subharmonic raw R^2=0.122 actually HIGHER than the 0.207 backbone
control's 0.059; canonicalization sometimes helped, sometimes hurt,
by comparable amounts for subharmonic and controls alike). Root cause
of the inconclusive result: the canonicalization tested was a 180-
degree sign flip, but the real ambiguity (per the amo evidence, ~81
degree phase difference between the two fits) is closer to a genuine
90-degree cos/sin quadrature swap -- which is NOT fixable by a sign
flip at all (cos and sin are orthogonal, not sign-related), so this
script was testing the wrong transform. No automatic fix exists per
the user's own diagnosis above; don't try to invent one algebraically.

**Decision (2026-09-29): not pursuing a full resolution now** -- doing
the real per-cell A/B optimization needed to resolve this properly
across even a handful of affected cells cuts against the fast-
turnaround priority established this session. Documented here as a
known, real, structural degeneracy for whenever it's picked back up;
`check_subharmonic_coherence.py` remains in the repo as a (currently
inconclusive) starting point, not a working diagnostic.

## Real state-recovery bug found and fixed at the source: stale per-CLIMATE_INDEX secondary checkpoint (2026-09-29)

User report: for `VALIDATE=FALSE`, re-running a cell didn't recover the
previously-saved optimal state -- e.g. `pdo`'s "last training CC=0.6398"
became "starts at CC=0.62731" on rerun, a real regression, not a display
artifact ("usually works reliably" -- this was new). Initial hypothesis
(that `Monitor.Check`'s cross-thread `Best_Client` selection uses each
iteration's live, possibly-rejected trial score rather than each
thread's own kept/accepted best) was plausible and independently
verified against the numbers (`lte_forward.py`, no Ada involved,
reproduced 0.6287 for the then-current saved state -- matching the
user's "0.62731" almost exactly) -- but the user pointed at the REAL,
simpler root cause instead: `lt.exe.<CLIMATE_INDEX>.p`, a secondary
per-index checkpoint file, still existed from this project's much
longer history (a "known use at one time" per the user) and had gone
stale relative to `lt.exe.p`.

Root cause, confirmed by reading `GEM.LTE.Primitives.Shared.Read`/
`Write_JSON` directly: `Read(D)` loads the PRIMARY file first (full
state, including `lpap`), then unconditionally loads the SECONDARY file
`lt.exe.<CLIMATE_INDEX>.p` INTO THE SAME `D` -- and `Read_JSON`'s
`Include_LPAP` flag (`False` for the secondary read) only ever gates
whether `lpap` gets parsed; every scalar, `ltep`, AND `harm` field gets
applied regardless of that flag. The secondary file never even contains
an `lpap` key. So the secondary read did nothing but silently OVERWRITE
the just-loaded primary's scalars/`ltep`/`harm` with the secondary
file's own (frequently stale) values whenever the two diverged -- and
`Write_JSON` writes both files together on every save, so they only
stay in sync as long as NOTHING external ever touches `lt.exe.p` alone
(e.g. this session's own Feb2026-import copying, `pick_best_of_feb_sep.py`,
never touched the secondary file for any of the 43 imported cells --
a real, live risk that existed until this fix).

**Fixed by removing the secondary file entirely** (both the read in
`Read` and the write in `Write_JSON`, `src/gem-lte-primitives-shared.adb`)
-- per the user's own diagnosis, `lt.exe.p` already contains everything
needed. Rebuilt clean (no warnings after also removing the now-dead
`FN2_JSON`/`FN2`/`CI` declarations). Verified directly on `pdo`: deleted
its stale `lt.exe.pdo.dat.p`, reran, confirmed the file is no longer
recreated and the recovered/saved state no longer reverts.

**Project-wide cleanup**: deleted all matching `lt.exe.*.dat.p` files
(the precise pattern Ada builds via `Base & "." & CI & ".p"` where CI
includes `.dat`, NOT the broader `lt.exe.*.p` glob, which also matches
unrelated files like the user's own `lt.exe.pysr.*.p` symbolic-
regression artifacts and numerous `lt.exe.p.<label>` dated snapshots --
careful with the exact pattern if repeating this) -- 89 in Sep2026, 139
in Feb2026, all now removed (gitignored already, so no tracked diff).
These files are also now fully inert even if any reappear (e.g. from an
old un-rebuilt binary, or copied in from elsewhere) -- the current
`enso_opt` no longer reads or writes them at all.

## Full ImpC=0.01 sweep, short (90s) runs: 35/53 cells improved (2026-09-29)

Follow-up to the stale-checkpoint fix above, per the user's explicit
direction: "continue on but try to proceed with short training runs.
Since the stale .p bug has been fixed, progress will be maintained
from run-to-run so that short durations will not cause problems
regressing to non-optimal sets."

`sweep_impc_full.py` re-run against all 53 cells whose `impC` was still
exactly `0.0` (i.e. never yet touched by the earlier 3-cell pilot, the
Feb2026-vs-Sep2026 import pass, or any organic escape from `0.0`), with
`TIMEOUT` dropped from 400s to 90s. Rationale confirmed correct in
practice: no cell regressed to a worse-than-starting state mid-sweep,
and per-cell wall time dropped roughly 4x with no loss of hit rate
(35/53 improved here vs the pilot's 3/3 and the general expectation
from the seed's own structural-bug rationale -- `Markov` never perturbs
a parameter sitting at exactly `0.0`, so `ImpC` was permanently
inaccessible to search until manually seeded).

**Result: 35/53 committed, 18 reverted to their pre-sweep state** (either
no improvement or one cascade failure, `kS040_E030`). Several large
single-attempt gains: `kN020_W170` +0.379, `kS040_W030` +0.307,
`kS020_E050` +0.201, `kN000_W050` +0.273, `kS040_E110` +0.245,
`kN000_W130` +0.147, `kN060_E170` +0.147, `kS020_W070` +0.142,
`kS020_W130` +0.075, plus many smaller positive deltas. Full per-cell
detail in `sweep_impc_full.log` and `SWEEP_IMPC_DONE`.

**One implementation wrinkle worth remembering**: `seed_impc()` mutates
`lt.exe.p` on disk (writing `impC=0.01`) BEFORE `backup()` snapshots the
cell, so a "rolled back" cell's `restore()` puts back the post-seed,
pre-cascade state, not pristine git HEAD -- i.e. the on-disk seed
survives a no-improvement rollback unless something else discards it.
Handled manually this round by identifying the 18 non-committed cells
from the log and `git checkout --`-ing them before committing, so the
repo stays consistent with the ledger (no orphaned uncommitted `impC`
seeds). Worth fixing in the script itself (move `backup()` before
`seed_impc()`) before the next run, rather than repeating the manual
cleanup step.

**Status map: 85 good/4 warning -> 87 good/2 warning/0 serious/0
critical/0 unattempted** -- the best collective state reached this
session. Remaining 2 warnings: `kS020_W090` (weak, 0.106-0.120, this
sweep tried it and made it worse, reverted) and `kS020_E070` (untouched
-- `impC` already nonzero from an earlier pass, so outside this sweep's
target set; a manual look is the natural next incremental step).

Per the user's own reasoning validated: short bursts + the fixed
checkpoint = compounding progress with fast turnaround, no need to
return to long single-cell searches for this kind of exploratory sweep.

## Suggested order of operations for a fresh session

1. Read this file, then skim `sweep.py` itself (well-commented, ~700 lines).
2. Re-run `kS020_E030` at N=2 (quick win candidate, dLOD was already fine).
3. Re-run `invalidate_and_refix.py` for the remaining ~23 invalidated cells
   at N=2 (should be much faster/more successful than the N=1 attempt).
4. Consider testing the Madagascar donor against its other neighbors
   (`kS040_E050`, `kS020_E070` if it exists, etc.) now that `kS020_E030` and
   `kN000_E050` have been tried.
5. Once the invalidated-cell backlog is cleared, consider widening beyond
   the original "poor-holdout" 43-cell candidate pool (`load_poor_cells()`)
   to the remaining ~40 grid cells never touched at all this sweep.
