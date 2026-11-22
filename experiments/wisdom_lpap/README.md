# Crowd-median ("wisdom of the crowds") tidal set

`../lpap_wisdom.py` reads the `lpap` of every saved `lt.exe.p` / `lt.exe.p.*` on the machine and
writes one median triple (period, amplitude, phase) per tidal line. Its docstring lists every
option. The steps:

- Copies of the same lpap count once.
- QBO-related sets are dropped.
- Lines with |period| < 2 d are dropped.
- Negative amplitudes are made positive, with their phase moved by π.
- Phases are folded into 0..2π and averaged as a circular mean.

| File | Contents |
|---|---|
| `wisdom_raw.lt.exe.p` | Median amplitudes. dLOD CC 0.9907 vs 0.985 for the single AMO/NAO fits. The benchmark. |
| `wisdom_wow.lt.exe.p` | Median of the fits saved as exceptional: any `lt.exe.p*` whose file or folder name holds "wow"/"WOW" (261 distinct sets, QBO excluded). dLOD CC 0.9913. |
| `wisdom_norm.lt.exe.p` | Each set scaled to unit variance before the median. dLOD CC 0.965. |
| `*_stats.csv` | Per line: count, period spread, amplitude median/trimmed/mean, phase, R (phase agreement, 0..1). |
| `dlod_calibrated.lt.exe.p` | The dLOD calibration itself: the solver's fit of the dLOD record at the 52 tidal periods (dLOD CC 1.000). Made 2026-11-20 by a `TEST_ONLY`, `DLOD_REF FALSE` run on a copy of `../Feb2026/amo`, which copies that fit into lpap; depends slightly on the fit's year correction. The reference for `plot_lpap_multi_bar.py --dlod`. |
| `wisdom.lt.exe.p`, `wisdom_lpap_stats.csv` | Earlier sorted 60-line version. Its row order mis-assigns lines in the model, so don't use it. |

All three sets are written with `--template ../Feb2026/amo/lt.exe.p`, which follows the model's
52-period table. Rebuild them with:

    python3 -I ../lpap_wisdom.py --template ../Feb2026/amo/lt.exe.p --no-normalize \
        -o wisdom_raw.lt.exe.p --csv wisdom_raw_stats.csv

and the wow set with:

    python3 -I ../lpap_wisdom.py / --include wow --skip '^/tmp/|leg5\.(bessel|ridge_wis)(/|$)' \
        --template ../Feb2026/amo/lt.exe.p --no-normalize \
        -o wisdom_wow.lt.exe.p --csv wisdom_wow_stats.csv

`--skip` leaves out scratch copies under /tmp and the two 2026-11-19 reruns of the 5-leg fit (one
was itself seeded from `wisdom_raw`). The wow set is almost the raw set: amplitude CC 0.9997 over the
42 lines, median phase difference 9°. Seed a run from it with `WISDOM_REF TRUE` and
`WISDOM_P=<path>/wisdom_wow.lt.exe.p`.

## Seeding a run: `WISDOM_REF TRUE`

`WISDOM_REF TRUE` (environment or `lt.exe.resp`; default off) is the counterpart of
`DLOD_REF FALSE`. Instead of starting lpap from the dLOD fit, a run starts it from this set:

- **File:** `WISDOM_P` when set; otherwise the first `experiments/wisdom_lpap/wisdom_raw.lt.exe.p`
  found in the run folder or one of its parents. With neither, the run stops with a message
  (exit 1).
- **Nonzero rows:** each is copied onto the table period within 1e-4 of it (relative). The match
  is tighter than the 1% used when reading `lt.exe.p`, so a row can't take a neighbor's slot.
- **Zero rows:** the 10 diurnal lines the set leaves out keep the loaded `lt.exe.p` values. A
  zero would never be searched again.
- **Combined with `DLOD_REF FALSE`:** the wisdom seed is applied last. The dLOD reference used by
  the `:dLOD:` check is unchanged either way.
- The log prints `Seeding lpap from wisdom <file>: 42 lines seeded, 10 rows left as loaded`.

Checks (2026-11-18, single-leg AMO, `TEST_ONLY`):

- With the option unset, scores are identical to the binary before it was added.
- `WISDOM_REF TRUE` gives the same scores as an `lt.exe.p` with the 42 rows pasted in by hand.
  dLOD goes from 0.9848 to 0.9907.
- The default search finds the file from a folder inside the repo.

## Seeding a run by name: `DLOD_REF <name>`

`DLOD_REF` (environment or `lt.exe.resp`, case-insensitive; 2026-11-20) takes:

| Value | lpap at the start of the run |
|---|---|
| `TRUE` or `REF` (default; also `T`, `YES`, `1`) | as loaded from `lt.exe.p`; the dLOD calibration is kept as the reference for the `:dLOD:` check |
| `FALSE` or `CAL` (also `F`, `NO`, `0`) | the dLOD calibration |
| any other `<name>` | seeded from `<name>.lt.exe.p` in the first `experiments/wisdom_lpap` at or above the run folder, the name matched case-insensitively (e.g. `wisdom_raw`, `wisdom_wow`, `dlod_calibrated`, `index`); the dLOD calibration is kept as the reference |

- **Named seed:** only the file's lpap rows are used; its other parameters are ignored. Each nonzero row goes to the
  table period within 1e-4 of it, and zero rows keep the loaded values (as `WISDOM_REF`). The log prints
  `DLOD_REF <name>: seeding lpap from <file>: N lines seeded, M rows left as loaded`.
- **Unknown name:** the run stops (exit 1) and lists the names available.
- **Checks** (`TEST_ONLY`, `Feb2026/amo`): `true`/`REF` and `false`/`Cal` give the previous binary's `TRUE`/`FALSE` scores
  exactly; `Wisdom_Raw` gives dLOD 0.9907 as `WISDOM_REF`; `index` and `INDEX` agree.
