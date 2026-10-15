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
`<idx>.dat` and `lte_run.sh` are run through `experiments/Feb2026/lte_forward.py`
on a date grid that extends 5 years past the end of the data. Nothing is
refit, and the tidal forcing is deterministic, so every month past the
data end is out of sample. The issue id is the first 12 hex digits of a
SHA-256 over those files plus the engine code. A new fit therefore creates
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

## Live observation sources (`sources.py`)

| target | feed |
|---|---|
| soi, soi_daily | Queensland LongPaddock monthly / daily SOI (Tahiti, Darwin, SOI columns) |
| qbo30, qbo50 | FU Berlin / KIT Canton-Gan-Singapore monthly winds |
| amo | NOAA PSL `amon.us.data`, falling back to `amon.us.long.data` |
| pdo | NOAA NCEI ERSSTv5 PDO |
| nino4, nino34 | NOAA CPC `sstoi.indices` |
| kaplan quads (all `k[NS]ddd_[EW]ddd` dirs) | ERSSTv5 20-degree box means (Kaplan v2 ends 2023-01; `kaplan_box` is kept for history) |
| baltic | PSMSL RLR composite of every gauge in 53.5-66N, 9.5-30.5E |
| brest, other MSL | PSMSL RLR by station id |
| lod, chandler | IERS EOP 20 C04 daily |

Live feeds are mapped onto the committed `.dat` units by a linear fit over
their overlap before the data end. MSL sites are de-seasoned first, because
their `.dat` has annual and semiannual removed upstream. `REPORT.md` shows
`obs cal r` for each target. A LOW value (< 0.9) means the live feed is not
yet the same series as the `.dat`. Commit the builder, or point the target
at the exact feed, before reading those scores.

## Extending

* New index: add one block to `targets.yml`, plus a fetcher in `sources.py`
  if the feed is new (it returns `(dates, values)`).
* New MSL site: `<dir>: {engine: lte, dir: <dir>, deseasonalize: true,
  source: {kind: psmsl_rlr, stations: [<PSMSL id>, ...]}}`.
* Daily fits: the `soi_daily` slot (and `qbo70_daily`, disabled) activate
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
