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

- `VALIDATE`, `ALTERNATE` and `SPLIT_TRAINING` are not supported with legs.
- `TIDES ZONAL` is untested with legs; its forcing cache is shared between tasks.
- With `TREND TRUE` the annual terms are regressed per leg and the saved `lt.exe.p` holds the
  leader's.
- Plot scripts that read `lte_results.csv` need pointing at the `.1` and `.2` files.

## Source

`src/gem.ads`, `src/gem.adb` (leg lookup, leg files, `LEGS` and `LEG_WEIGHT` keys);
`src/gem-lte-primitives-solution.ads/.adb` (pairs, `Exchange`, leader and follower paths);
`src/gem-lte-primitives.adb` (`SAMPLING` lag); `src/gem-lte-primitives-shared.adb` (windings file
name). Backups: `src/*.10-12-2026-1-pre-legs`; previous binary `experiments/Feb2026/enso_opt.pre-legs`.
