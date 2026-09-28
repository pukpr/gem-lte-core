# Experiments

This directory contains the working experiment spaces that sit alongside the
core Ada solver.

## February 2026 workspace

`Feb2026/` is the main interactive and research-oriented workspace. It contains:

- the staged runtime executable location expected by the Python wrappers
- per-index subdirectories with `<index>.dat`, `lt.exe.p`, and `lt.exe.resp`
- the Tk GUI launcher `lte_gui.py`
- command-line helpers such as `lte_run.py`, `lte_forward.py`, and
  `lte_sweep.py`
- supporting analysis and plotting tools

Typical usage from the repository root:

```bash
cd experiments/Feb2026
python3 lte_gui.py
```

Command-line equivalent of the GUI run button:

```bash
cd experiments/Feb2026
python3 lte_run.py amo --cv 1880 1885
```

For a fuller tool summary, see `Feb2026/README.md`.

## Round-robin workspace

`rr/` is a smaller helper workspace for round-robin and Pareto-style parameter
experiments. It includes its own copies of the lightweight Python helpers plus a
small set of index directories (`amo`, `nao`, `nino4`, `pdo`).

Typical usage from `experiments/rr`:

```bash
python3 round_robin_forcing.py --rounds 3 --indices amo nao pdo nino4
```

Or to inspect saved fits without running the optimizer:

```bash
python3 lte_forward.py amo
```
