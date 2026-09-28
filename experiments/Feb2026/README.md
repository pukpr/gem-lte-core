# February 2026 Experiment Workspace

This directory is the main GEM-LTE core working area for running fits, checking
saved parameter sets, and exploring the cross-index research tooling.

## Workspace layout

Most runnable index/site directories in this folder follow the same convention:

- `<index>/<index>.dat` — target series
- `<index>/lt.exe.p` — saved parameter state
- `<index>/lt.exe.resp` — run-time option file

The Python tools in this directory assume the Ada runtime executable is present
here as `lt` or `lt.exe`.

## Primary tools

### GUI launcher

`lte_gui.py` is the interactive front end for browsing index directories,
running the solver, and previewing plots.

Typical usage:

```bash
python3 lte_gui.py
```

### Command-line runner

`lte_run.py` is the command-line equivalent of the GUI's **Run lt** action. It
loads defaults from the selected directory's `lt.exe.resp`, then launches the
staged `lt` / `lt.exe` binary with matching environment variables.

Typical usage:

```bash
python3 lte_run.py amo
python3 lte_run.py amo --cv 1880 1885
python3 lte_run.py nino4 --timeout 3600 --metric DTW
python3 lte_run.py baltic --no-exclude --cv 1940 2020 --ridge 0.1
```

Outputs are written back into the selected index directory, including
`lte_results.csv` and updated parameter state.

### Deterministic forward model

`lte_forward.py` keeps the forward computation but removes the Ada
random-descent search. Use it when you want to inspect or verify what a saved
parameter set is doing.

Typical usage:

```bash
python3 lte_forward.py amo --verify
python3 lte_forward.py amo --cv 1880 1885
```

### Local sweep helper

`lte_sweep.py` runs isolated parameter sweeps around a saved fit by cloning an
index directory into temporary trial folders, invoking `lte_run.py`, and
preserving the trial outputs.

Typical usage:

```bash
python3 lte_sweep.py amo
python3 lte_sweep.py amo --target 0.9 --timeout 300
python3 lte_sweep.py amo --drop-low-modulation --year-offsets=-.001,0,.001
```

## Analysis and plotting helpers

Some of the most directly useful inspection tools are:

- `param_survey.py` — compare persisted annual/semi-annual parameters and IR
  across a default or explicit set of indices
- `plot_lpap_multi_bar.py` — compare absolute LPAP tidal amplitudes across many
  index directories
- `plot_lte_parameters_multi_bar.py` — compare LTE parameter magnitudes across
  saved runs
- `ws.py`, `winding_rank.py`, `winding_scalogram.py`, and related scripts —
  winding-number and forcing diagnostics

Typical usage:

```bash
python3 param_survey.py
python3 param_survey.py nino4 pdo amo
python3 plot_lpap_multi_bar.py amo nao pdo nino4 --mode grouped --no-show
```

## Research notes in this directory

This folder also contains several Markdown write-ups that explain the modelling
ideas, experimental findings, and replication workflow, especially:

- `RECIPE_FOR_RESEARCHERS.md`
- `OPTIMIZATION_GUIDE.md`
- `WINDING_NARRATIVE.md`
- `WINDING_SCALOGRAM_FEASIBILITY.md`

Use those for methodological context, and use the tools above for the actual
runtime workflow.
