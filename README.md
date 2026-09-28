# GEM-LTE Core

Core Ada solver and experiment workspace for the GEM-LTE modelling pipeline.
This repository keeps the LTE engine and the `experiments/Feb2026` research
tools together, but without the broader packaging/setup scripts from
[`pukpr/GEM-LTE`](https://github.com/pukpr/GEM-LTE).

## Build the Ada solver

This repository builds with GNAT/GPRbuild using `lte.gpr`.

```bash
gprbuild -P lte.gpr enso_opt
```

The project file writes executables into `experiments/Feb2026/`. The Python
helpers there expect the runtime binary to be named `lt` or `lt.exe`, so after
building, stage the produced `enso_opt` executable under that name in
`experiments/Feb2026/`.

Unlike the full `GEM-LTE` repository, this core repo does **not** bundle the
Windows setup/staging scripts, so building and staging the runtime executable is
done directly.

## Essential operation

The main Ada driver is `src/enso_opt.adb`. Runtime fits are controlled by the
per-index files stored under `experiments/Feb2026/<index>/`:

- `<index>.dat` — target time series
- `lt.exe.p` — saved parameter state
- `lt.exe.resp` — runtime options / defaults

Typical workflow:

1. Build `enso_opt` with `gprbuild -P lte.gpr enso_opt`.
2. Stage the binary as `experiments/Feb2026/lt` or
   `experiments/Feb2026/lt.exe`.
3. Run the solver from the experiment workspace with either the GUI or the CLI
   helper wrappers.

## Main usage paths

### GUI workflow

From `experiments/Feb2026/`:

```bash
python3 lte_gui.py
```

The GUI discovers subdirectories that contain `lt.exe.p`, lets you choose an
index/site, and launches the staged `lt` / `lt.exe` binary with the matching
environment variables.

### Command-line run helper

Use `lte_run.py` when you want the same operation as the GUI's **Run lt**
button, but from the shell:

```bash
cd experiments/Feb2026
python3 lte_run.py amo --cv 1880 1885
```

Example with a different metric and timeout:

```bash
python3 lte_run.py nino4 --timeout 3600 --metric DTW
```

### Deterministic forward reproduction

`lte_forward.py` runs the forward model only, without the Ada random-descent
optimizer:

```bash
cd experiments/Feb2026
python3 lte_forward.py amo --verify
```

This is the easiest way to inspect the saved forcing/manifold/LTE pipeline for a
given parameter set.

### Local parameter sweeps

`lte_sweep.py` runs isolated sweeps around a saved fit:

```bash
cd experiments/Feb2026
python3 lte_sweep.py amo --target 0.9 --timeout 300
```

## Repository layout

- `src/` — Ada packages and drivers
- `lte.gpr` — GNAT project file
- `experiments/Feb2026/` — primary research workspace, index datasets, GUI, and
  analysis helpers
- `experiments/rr/` — round-robin / pared-down helper workspace

## See also

- `experiments/README.md`
- `experiments/Feb2026/README.md`
- `experiments/Feb2026/RECIPE_FOR_RESEARCHERS.md`
