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
| `wisdom_norm.lt.exe.p` | Each set scaled to unit variance before the median. dLOD CC 0.965. |
| `*_stats.csv` | Per line: count, period spread, amplitude median/trimmed/mean, phase, R (phase agreement, 0..1). |
| `wisdom.lt.exe.p`, `wisdom_lpap_stats.csv` | Earlier sorted 60-line version. Its row order mis-assigns lines in the model, so don't use it. |

Both sets are written with `--template ../Feb2026/amo/lt.exe.p`, which follows the model's
52-period table. Rebuild them with:

    python3 -I ../lpap_wisdom.py --template ../Feb2026/amo/lt.exe.p --no-normalize \
        -o wisdom_raw.lt.exe.p --csv wisdom_raw_stats.csv

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
