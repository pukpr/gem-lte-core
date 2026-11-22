# Leg groups: one search, two fits

`LEGS 2` divides the worker tasks of `lt.exe` (`enso_opt`) into two leg groups that share one
parameter file. Each leg can have its own index file, sampling, training interval and any other
response-file setting. The search accepts or rejects every candidate on a composite of the two
legs' scores.

Added 2026-10-12. Without `LEGS`, the program behaves exactly as before (output bit-identical).

**More than two legs.** `LEGS N` (2 to 9) runs groups of N tasks: one leader on leg 1 and a follower
on each of legs 2 to N, with settings from `lt.exe.resp.1` to `lt.exe.resp.N` and outputs
`lte_results.1.csv` to `lte_results.N.csv`. The composite is the weighted mean over all N legs, and
`NUMBER_OF_PROCESSORS` is rounded down to a multiple of N. Everything below written for two legs
applies unchanged.

## Files

| File | Purpose |
|---|---|
| `lt.exe.resp` | Common settings, plus `LEGS 2` |
| `lt.exe.resp.1` | Settings for leg 1 only |
| `lt.exe.resp.2` | Settings for leg 2 only |
| `lt.exe.p` | The one shared parameter file (read at start, saved at the end) |
| `lte_results.1.csv`, `lte_results.2.csv` | Results per leg, same columns as `lte_results.csv` |
| `lt.exe.windings.1.json`, `lt.exe.windings.2.json` | Regression amplitudes per leg |

## Settings

- **`LEGS 2`** in `lt.exe.resp` or the environment turns the mode on.
- **Leg files** use the same `KEY value` format as `lt.exe.resp`. Any key the response file accepts
  may appear: `CLIMATE_INDEX`, `SAMPLING`, `TRAIN_START`, `TRAIN_END`, `EXCLUDE`, `F9`, `METRIC`,
  `DTW`, and so on.
- **`LEG_WEIGHT`** (default 1.0), set in a leg file, is that leg's weight in the composite.
- **Precedence:** leg file, then `lt.exe.resp`, then environment variables, then the default. An
  `export` in `lte_run.sh` therefore does not override a value given in a leg file.
- **Threads:** `NUMBER_OF_PROCESSORS` is rounded down to an even number; each pair of tasks is one
  independent search.

## Choosing which search to save: `SELECT_TOP`

Works with or without legs. Default 0: the search with the best training score saves, as before.

- **`SELECT_TOP N`**: when every search has stopped, each reports its kept candidate's training and
  held-out scores. Of the N best training scores lying within `SELECT_TIE` of the best, the one with
  the highest held-out score is saved.
- **`SELECT_TIE`** (default 0.01): how close to the best a training score must be to count as tied.
- The log lists the tied candidates with both scores and the one selected.
- With legs the scores are the pair's composites; with `EXCLUDE TRUE` the held-out score is the
  excluded interval.
- Searches are not independent under `THRESHOLD_ACTION CATCHUP`, and choosing the best of N on a
  short held-out interval raises that score by selection (about 0.85 of its standard deviation for
  N = 3). Treat the saved held-out score as a tie-breaker, not as a test.

## Choosing by a validation interval: `VALIDATE` with legs

`VALIDATE TRUE` never steers the search. When every search has stopped, it picks which one to save
by its score on a validation interval. Without legs that was already so. With legs (2026-11-19):

- **Scoring:** each leg scores its kept candidate on the interval. The group's score is the
  `LEG_WEIGHT`-weighted mean over its legs, and the group with the best one is saved, instead of
  the best training composite.
- **Log:** one line per group (`VALIDATE group N: leg 1 validate ..., composite validate ...,
  training composite ...`), then the `SELECT_TOP` lines with the validate score in the first slot.
  The `V:` in the status line is leg 1's score only.
- **`VALIDATE_START`** (date; default `TRAIN_END`, as before): where the interval begins, so a test
  gap can sit between training and validation. Train before 1960, test 1960-1970, validate from 1970:

      EXCLUDE FALSE  TRAIN_START 1880  TRAIN_END 1960  VALIDATE TRUE  VALIDATE_START 1970

  `leg_scores2.py --held 1960 1970` reports the test decade.
- **Groups:** the choice is between groups, so it needs several: threads / legs, e.g. 18 threads
  for 3 groups of 6 legs.
- **`MANIFOLD_PRE TRUE`** works with it: column 4 is the shared pre-Bessel manifold of the saved
  group.
- **Checks** (2026-11-19, the 6-leg `.leg6` fit):
  - With `VALIDATE` unset, and with `VALIDATE TRUE` in `TEST_ONLY`, the new binary's scores and
    files are identical to the previous one's.
  - 90 s search, 18 threads, `VALIDATE_START 1970`: group 3 was saved (validate 0.624, training
    0.762) over group 2 (best training, 0.766). The six saved `lte_results` reproduce 0.62434 as
    the mean of the legs' 1970-on CCs, and their column 4 agrees across legs to 0.9999.

## Legs that are the same index: `LEG_COMMON`

For two different indices each leg regresses its own amplitudes and phases, which is right. For two
samplings or two intervals of the SAME index it leaves open whether one response fits both. These
settings close that gap.

- **`LEG_COMMON H`** (common file; default 0 = off). A follower does not regress the harmonics above
  `H` itself. It applies the amplitudes and phases the leader regressed on its own leg, and regresses
  only the remaining terms (level, slope, annual, slow winding, base, harmonics up to `H`) on what is
  left of its data.
- **`LEG_SCALE`** (leg file; default 1.0) multiplies the leader's amplitudes on that leg, for series
  in different units.
- **`LEG_SCORE FAST`** (leg file; default `FULL`). The leg is scored on the leader's fast part alone:
  the metric between that part of the model and the data less the leg's own remaining model, over the
  leg's interval. The composite then rewards a manifold on which the fast response calibrated on leg 1
  is present on the other legs.
- At save each leg prints its own score beside `COMPOSITE`.
- With two identical legs the follower reproduces the leader's model to 1e-14 (a least-squares
  identity), which is the built-in check.
- Use `ALPHA 0` and `NONLIN 1.0` with it.

## Each leg's own windings: `LEG_WINDINGS`

By default every leg uses the same winding wavenumbers (`ltep`, the `harm` multipliers) and `IR`
from `lt.exe.p`, and only the regressed amplitudes and phases differ. With **`LEG_WINDINGS TRUE`**
(common file; default off) the manifold stays shared, but the parameters after it are each leg's own:

- **Files.** Each leg reads `lt.exe.p.<leg>` (`lt.exe.p.1`, `lt.exe.p.2`, ...) after `lt.exe.p`, and
  takes from it whichever of `ltep`, `harm` and `IR` it holds. A missing file means the leg starts
  from the shared values. At save every leg writes its own `lt.exe.p.<leg>`. The leader still writes
  `lt.exe.p` (with leg 1's values), so that file stays complete on its own.
- **Counts.** `NM` and `NH` come from the main `lt.exe.resp`, so every leg has as many `ltep` and
  `harm` values. Do not set `NM` or `NH` in a leg file.
- **Search.** The leader walks the shared parameters and its own leg's values as before. Each follower
  puts its own `ltep`, `IR` and multipliers back over the leader's candidate and steps them itself:
  one `Markov` step on one of its `ltep` / `IR` values at the rate the leader's walk gives any one
  parameter, and the usual `Random_Harmonic` step on its multipliers. Accept and reject stay joint
  (the composite), and a rejected step reverts each leg's own values.
- **Zeros never move**, as in the walk: a leg whose `IR` should be searched needs a nonzero start.
  `LEG_IR` still fixes a leg's `IR`; `LOCKW` snaps followers' `ltep` too.
- **The manifold stays shared.** Its Bessel correction used to take its wavenumber from `ltep`
  (`ltep(NM)`, or `ltep(NM-1)` with `LOCKF`), which made each leg's own `ltep` move its manifold.
  It is now the shared `lt.exe.p` parameter **`besk`** (2026-11-17), searched with the other
  manifold parameters; a `.p` without it is seeded from that `ltep`, so older fits reproduce.
  `LAYER` still reads `ltep(1)` and `ltep(2)` and so does not keep the manifold shared.
- **`LEG_P`** (common file or environment, with `LEG_WINDINGS`; any case) sets how much more
  than `ltep`, `harm` and `IR` each leg owns. The extra keys go into `lt.exe.p.<leg>`, and each
  follower steps them with its `ltep` and `IR`:
  - **`NONE`** (default): nothing more.
  - **`BESSEL`**: the Bessel correction, `impA`, `impB`, `besk`, `offs`, `bg` and `bg2` (`bg2` is
    stepped only with `NM 1`, the one case that uses it). Only the tidal forcing and its impulse
    integration stay shared.
  - **`PARAMS`**: every `.p` parameter except `lpap` and `year`. That is the Bessel set plus `impC`,
    `delA`, `delB`, `delC`, `asym`, `ann1`, `ann2`, `sem1`, `sem2`, `ma`, `mp`, `mb`, `shfT`, `init`,
    `rho`, `ztau`, `zdist`, `ztilt`, `ssaA`, `ssaP`, `cw1` and `cw2`. Only the tidal constituents
    and the year length stay shared. `ann*`/`sem*` (with `TREND`) and `cw*` are regression outputs,
    so they are each leg's own in any mode; `PARAMS` only carries them in the leg file.

  With `BESSEL` or `PARAMS`, each leg's manifold column is its own. A leg file without some of
  these keys starts those from the shared values. `LEG_BESSEL TRUE`, the earlier switch, still
  means `LEG_P BESSEL`.
- **`MANIFOLD_PRE TRUE`** (common file or environment; default off; also without legs) puts the
  manifold as it is *before* the Bessel step (or its `LAYER` / `FRICTION` variant) into column 4
  of `lte_results.csv` / `lte_results.<leg>.csv`. That is the tidal forcing and its impulse
  integration only, after a stage-0 `ANNUAL_DITHER`. It is cached with each accepted candidate, so
  it matches the saved fit. Columns 1–3 and all scores are unchanged, and `leg_scores2.py`'s
  manifold line and bottom panel then compare the legs' manifolds before Bessel. It is not
  applied with `VALIDATE TRUE` without legs, where the save can take another thread's candidate;
  the run then prints a note and writes the usual post-Bessel column.
  - Check (2026-11-18, the 4-leg `amo_pdo_nao_nino4_feb_overnight.wow` fit, `LEG_P BESSEL`):
    applying each leg's Bessel step to the new column reproduces the post-Bessel column to 4e-13.
    The pre-Bessel manifolds of the four legs agree to 1.00000, against 0.963 / 0.941 / 0.898
    after Bessel, so the whole spread between them comes from each leg's Bessel step.
- **Checks for `LEG_P`** (2026-11-17, `TEST_ONLY`, 2-leg AMO + NAO):
  - Unset gives the same scores as the binary before `LEG_P`.
  - `BESSEL`, `LEG_BESSEL TRUE` and `PARAMS` with leg files holding only the Bessel keys all give
    the same scores as the earlier `LEG_BESSEL TRUE`.
  - With NAO's full fit as `lt.exe.p` and every non-lpap, non-year key of each single-leg fit in
    `lt.exe.p.<leg>`, the NAO leg reproduces its single-leg 0.8505 in both modes. The AMO leg goes
    from 0.483 (`BESSEL`) to 0.585 (`PARAMS`).
  - A `PARAMS` search writes all 31 keys to each leg file, and never `lpap` or `year`.
- **Use `STRICT_IDATE TRUE`** for a strictly common manifold: without it each leg's forcing is
  integrated from its own first data date (Baltic July 1880, NAO January 1880 gives a manifold
  CC of 0.9994 between the legs, in every run before this feature too).
- Off with `LEG_COMMON`, which applies the leader's amplitudes to windings it assumes are common.
- `leg_scores2.py` lists each leg's wavenumbers when they differ.
- Checks (2026-11-15): off, and on without per-leg files, is bit-identical to before; a hand-set
  `lt.exe.p.2` changes leg 2 only; reloading a search's saved files reproduces every leg's score.

Other per-leg settings added for particular indices: **`LEG_IR`** (the leg's own one-year lag
coefficient; `1.0` with `UNCOMPENSATED TRUE` makes the leg the one-year difference of its response, as
for NAO against AMO) and **`ALPHA_GAIN`** (positive-side gain, now readable per leg).

## Examples

One series, two intervals:

```
lt.exe.resp.1                  lt.exe.resp.2
TRAIN_START  1992              TRAIN_START  1992
TRAIN_END    2005              TRAIN_END    2005
EXCLUDE      FALSE             EXCLUDE      TRUE
```

Two samplings of the same index:

```
lt.exe.resp.1                        lt.exe.resp.2
CLIMATE_INDEX  soi_lag_5.dat         CLIMATE_INDEX  soi_lag_m.dat
SAMPLING       73.0                  SAMPLING       12.0
TRAIN_START    1992                  TRAIN_START    1992
TRAIN_END      2005                  TRAIN_END      2005
EXCLUDE        FALSE                 EXCLUDE        TRUE
```

For two samplings to share a manifold, keep `IMPULSE_EXACT TRUE`, `INIT_DATE` and `STRICT_IDATE`
in the common file.

## Running

```
cd <run directory>          # holds lt.exe.p, lt.exe.resp, lt.exe.resp.1, lt.exe.resp.2, data files
export TIMEOUT=1800 NUMBER_OF_PROCESSORS=16 TEST_ONLY=FALSE
ulimit -s 65536
../lt.exe
```

- Keys while running: `q` or `s` stop and save; `x` exits without saving.
- `TEST_ONLY=TRUE` scores the current `lt.exe.p` on both legs, writes the leg files and exits.
- The `Status:` lines show the best composite and the composite of the two legs' out-of-interval
  scores. At save, each leg prints its own two scores and a `COMPOSITE` line.

## How it works

- **Per-task settings.** `GEM.Getenv` checks the calling task's leg (a task attribute set when the
  worker starts) and looks in that leg's file first. Settings read inside a worker are therefore
  the leg's own without any change to the code that reads them.
- **Pairs.** Tasks are started in pairs: odd-numbered tasks are leaders on leg 1, even-numbered
  ones followers on leg 2.
- **Exchange.** Each pair shares a protected `Exchange`. The leader posts a candidate, both tasks
  evaluate it at the same time on their own legs, and the exchange returns
  `(w1*s1 + w2*s2) / (w1 + w2)`. Both tasks accept or reject on that number, so they keep the same
  candidate.
- **Regression amplitudes are per leg.** Only the parameters in `lt.exe.p` are shared.
- **Saving.** The choice of pair is made once, centrally, when every search has stopped (see
  `SELECT_TOP`; 0 acts as 1 with legs). That pair saves: the leader writes `lt.exe.p` and leg 1's
  files, the follower leg 2's files.
- **`SAMPLING` lag.** The one-year lag for `IR` is now read where it is used, so each leg gets its
  own.

## Tests

`experiments/Oct2026_generator/legs/legs_test.py`, from the `soi_daily/rr/t_daily` fit:

| Test | Result |
|---|---|
| No legs, new binary against the previous one (5-day and monthly) | Identical |
| Two legs with identical settings | Both leg files identical to the single run |
| One series, leg 1 trains 1992–2005, leg 2 trains outside it | Each leg file identical to its single run; composite equals the mean of the two single scores |
| 5-day leg plus monthly leg | Each leg file identical to its single run; composite equals the mean |
| 90 s paired search, 4 pairs | Exits by itself, composite rises, no errors |
| Three legs, one series, three intervals | Each leg file identical to its single run; composite equals the mean of the three |
| 60 s three-leg search, 3 groups | Exits by itself, three leg files written, no errors |
| Re-run of the saved `lt.exe.p` | Both leg files reproduced to 8e-7 |

## Limits

- `ALTERNATE` and `SPLIT_TRAINING` are not supported with legs. `VALIDATE` is (see above).
- `TIDES ZONAL` is untested with legs; its forcing cache is shared between tasks.
- With `TREND TRUE` the annual terms are regressed per leg and the saved `lt.exe.p` holds the
  leader's.
- Plot scripts that read `lte_results.csv` need pointing at the `.1` and `.2` files.

## Source

`src/gem.ads`, `src/gem.adb` (leg lookup, leg files, `LEGS` and `LEG_WEIGHT` keys);
`src/gem-lte-primitives-solution.ads/.adb` (pairs, `Exchange`, leader and follower paths);
`src/gem-lte-primitives.adb` (`SAMPLING` lag); `src/gem-lte-primitives-shared.adb` (windings file
name). Backups: `src/*.10-12-2026-1-pre-legs`; previous binary `experiments/Feb2026/enso_opt.pre-legs`.
