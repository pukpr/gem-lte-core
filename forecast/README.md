# GEM-LTE forecast ledger

A weekly, self-scoring record of what every committed GEM-LTE fit says
about data it was never fitted to. It runs as a GitHub Action; the git
history provides the timestamps.

```
python3 forecast/ledger.py all                 # observe + issue + score + report
python3 forecast/ledger.py backfill            # rolling-origin hindcasts (LOD, Chandler)
python3 forecast/ledger.py list                # registered targets
python3 forecast/ledger.py all -t soi qbo30    # subset (target or group names)
```

Results: [`REPORT.md`](REPORT.md), per-issue scores in `scores/`, frozen
projections in `ledger/<target>/<issue_id>/{projection.csv,manifest.json}`.

## What is scored

Frozen LTE fits (`engine: lte`). The committed `lt.exe.p`, `lt.exe.resp`,
`<idx>.dat` and `lte_run.sh` are run by the Ada solver itself (`lt.exe`, in
`TEST_ONLY` mode, with `FORECAST_YEARS` set to the horizon). `lt.exe` extends
the date grid past the last sample on the data's own step and computes the
manifold and model there. The extension takes no part in the regression, the
data smoother or the scores. A second run without the extension checks that
the fit on the real samples is unchanged (`extension_max_dev_over_std` in the
manifest). The scores in the manifest are the ones `lt.exe` prints for the
same fit. Nothing is refit, and the tidal forcing is deterministic, so every
sample past the data end is out of sample. The issue id is the first 12 hex
digits of a SHA-256 over those files, the solver source (`src/`, `lte.gpr`), the engine code and
any `env:` override in `targets.yml`.

The solver is found as `$LTE_EXE`, else `experiments/Feb2026/lt.exe`, else
`experiments/Feb2026/enso_opt` (the `gprbuild -P lte.gpr enso_opt` output).
`env:` in a target block adds `lt.exe` environment, e.g. `TIDE_LEGACY: "TRUE"`
for fits made before the 2026-10-06 tide fix that do not set it in `lte_run.sh`.

Issues before 2026-10-15 21:00 UTC used `lte_forward.py`, a Python
re-implementation that does not reproduce the current fits (NINO4: +0.26
training against +0.71 from `lt.exe`). They stay in the ledger as history;
the report shows the most recent issue per target. A new fit therefore creates
a new issue, and old issues keep being scored.

Each projection is scored over three segments:

| segment | meaning |
|---|---|
| `holdout` | the TRAIN_START..TRAIN_END interval excluded from the committed fit (EXCLUDE=true) |
| `post_data` | after the last sample in the committed `.dat`; never seen by the fit |
| `post_issue` | after the ledger timestamp; strictly ex-ante |

Baselines are issued from the same origin using only data available then.
They are persistence, AR(1)-damped persistence on anomalies (30-year window),
and calendar-month climatology. `ss_x = 1 - MSE(model)/MSE(x)`. Observations
get the same `F9` (1,2,1)/4 passes as the fit. The last `F9` months are
held back until the centred smoother is final.

Refit engines (`tidal_lod`, `alias_cw`) are refit once a month from live
IERS EOP C04, and `backfill` adds rolling-origin hindcasts every 3 months
over the last 10 years. Each one uses only data available at its origin.

* `lod`: the Doodson constituents in a donor `lt.exe.p` (default `amo`, the
  LOD-calibrated `lpap` set), Rayleigh-de-duplicated for the training span,
  plus annual/semiannual terms and a quadratic. Ridge least squares on 12 years of daily LOD.
* `chandler_x`, `chandler_y`: prograde complex fit at the forced-alias
  period 365.242 / (2 frac(365.242/27.2122)) = 432.7 d, plus annual and pole drift, on 24 years.

## Deferred targets

Disabled in `targets.yml` (`enabled: false`), 2026-10-15: `amo` (NOAA's AMO
series ends in 2023), the Kaplan quads (Kaplan Extended SST v2 ends 2023-01)
and `brest` (set aside for now).

## Live observation sources (`sources.py`)

| target | feed |
|---|---|
| soi, soi_daily | Queensland LongPaddock monthly / daily SOI (Tahiti, Darwin, SOI columns) |
| qbo30, qbo50, qbo70 | FU Berlin / KIT Canton-Gan-Singapore monthly winds (scored on the monthly leg of the two-leg runs) |
| amo (deferred) | NOAA PSL `amon.us.data`, falling back to `amon.us.long.data` |
| pdo | NOAA NCEI ERSSTv5 PDO |
| nino4, nino34 | NOAA CPC `sstoi.indices` |
| kaplan quads, deferred (all `k[NS]ddd_[EW]ddd` dirs) | ERSSTv5 20-degree box means (Kaplan v2 ends 2023-01; `kaplan_box` is kept for history) |
| baltic | PSMSL RLR composite of every gauge in 53.5-66N, 9.5-30.5E |
| brest (set aside), other MSL | PSMSL RLR by station id |
| lod, chandler | IERS EOP 20 C04 daily |

Live feeds are mapped onto the committed `.dat` units by a linear fit over
their overlap before the data end. MSL sites are de-seasoned first, because
their `.dat` has annual and semiannual removed upstream. `REPORT.md` shows
`obs cal r` for each target. A LOW value (< 0.9) means the live feed is not
yet the same series as the `.dat`. Commit the builder, or point the target
at the exact feed, before reading those scores.

## To do: running it as a GitHub Action

Exercised locally so far (2026-10-15): `python3 forecast/ledger.py all` runs end to end in
about 4 minutes, and a second run issues nothing new. Steps before the scheduled
workflow should run this version:

1. **Commit the fits the ledger reads.** The working-tree `lt.exe.p` / `lt.exe.resp` of
   `soi`, `pdo`, `nino4`, `nino34` and `baltic` in `experiments/Feb2026/` were refit for the
   current tides and are uncommitted; the committed versions predate the 2026-10-06 tide fix
   and would score far lower. The QBO leg runs (`legs/qbo30_point0E`, `qbo50_point0E`,
   `qbo70_point0E`) are already committed.
2. **Check from a clean checkout** (`git worktree add`) that `ledger.py issue` gives the same
   issue ids as the local tree. Ids hash the run files, `src/`, `lte.gpr` and `engines.py`,
   not the binary, so a CI build of the same commit should match.
3. **Workflow** (`.github/workflows/forecast-ledger.yml`, rewritten 2026-10-15): installs GNAT
   13 on ubuntu-24.04, builds `enso_opt` (cached on `src/`, `lte.gpr`), sets `LTE_EXE`, runs the
   ledger and commits the results. It is manual-only (`workflow_dispatch`); after step 1, enable
   the schedule and push triggers listed in its header comment. Pushing it needs a token with the
   `workflow` scope.
4. **First manual run on GitHub**: check that the build succeeds, that `ulimit -s 65536` is
   accepted, and that the issue ids match the local ones (no NEW issues).
5. **Live feeds that do not match their `.dat`**: `pdo` (obs cal r 0.65) and `baltic` (0.41).
   Point them at the exact series used for the fit, or mark their post-data scores as
   provisional.
6. **Starting point**: the issues made with `lte_forward.py` were removed on 2026-10-15 (git
   history keeps them), with the scores and observations of the deferred targets. The ledger
   starts from the Ada-engine issues.
7. **Old workflow**: its schedule is removed by commit `2203fa4` (branch
   `disable-forecast-schedule`, on top of `origin/master` 1be4348); push it with
   `git push origin disable-forecast-schedule:master`, or disable the workflow in the GitHub
   Actions page.
8. Optional: a newer ERA5 70 hPa point series (the current one ends 2020.45) would let the
   `qbo70` 5-day leg run to the present like the 30 and 50 hPa legs.

## Extending

* New index: add one block to `targets.yml`, plus a fetcher in `sources.py`
  if the feed is new (it returns `(dates, values)`).
* New MSL site: `<dir>: {engine: lte, dir: <dir>, deseasonalize: true,
  source: {kind: psmsl_rlr, stations: [<PSMSL id>, ...]}}`.
* Multi-leg runs: `legs_dir: <repo path>` and `leg: N` instead of `dir:`. All legs run
  together as fitted; leg N is scored, and each leg's projection is saved as `legN.csv`.
* Daily fits: the `soi_daily` slot activates
  once `experiments/Feb2026/<dir>/lt.exe.p` with `SAMPLING` set for daily
  data is committed. The projection grid follows `SAMPLING`.
* New Kaplan quad: commit the fit dir; the glob picks it up.

## Engine fidelity

`lte_forward.py` is the Python mirror of the Ada `lt` binary. On PDO, Ada in
`TEST_ONLY` mode gives r = 0.217 and the Python forward gives r = 0.204 for the
same committed parameters. Each manifest records `r_train`, `r_holdout`,
`r_all` and `extension_max_dev_over_std`. The last checks that extending the
date grid leaves the in-sample model unchanged; it is about 0 because
zero-valued padding is dropped from the regression and the data smoother is
edge-held.
