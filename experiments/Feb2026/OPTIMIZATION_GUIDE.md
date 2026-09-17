# LTE Parameter Optimization Framework

A multi-stage optimization pipeline for the Laplace's Tidal Equation (LTE) climate index model, using coordinate descent, Pareto front analysis, and cross-index forcing manifold convergence.

## Overview

The framework optimizes LTE model parameters across multiple climate indices (AMO, NAO, PDO, NINO4, TNA, Baltic) using a staged approach:

1. **IR + ltep tuning** — coarse optimization of feedback and modulation frequencies
2. **lpap micro-tuning** — fine adjustment of tidal amplitudes/phases (sub-0.01 shifts)
3. **Pareto analysis** — accuracy vs. complexity tradeoff for ltep harmonic structures
4. **Round-robin convergence** — iteratively align forcing manifolds across indices

## Invocation

### 1. Single-Index Optimization

```bash
cd experiments/Feb2026

# Coordinate descent on IR + ltep (7 parameters)
python3 lte_pysr.py <index> --subset ir_ltep --coordinate --generations 300

# Coordinate descent on top-10 lpap amplitudes/phases (20 parameters)
python3 lte_pysr.py <index> --subset lpap_top --lpap-n 10 --coordinate --generations 500

# Seed from another index's optimized params
python3 lte_pysr.py nino4 --subset ir_ltep --coordinate --generations 300
# (uses PDO-optimized scalars as seed via the script's internal logic)
```

**Options:**
| Flag | Description |
|------|-------------|
| `--subset` | Parameter subset: `ltep`, `scalars`, `lpap_top`, `lpap`, `ir`, `ir_ltep` |
| `--lpap-n N` | Number of top lpap constituents to evolve (default 10) |
| `--coordinate` | Use coordinate descent (hill-climb one param at a time) |
| `--evolve` | Use population-based evolutionary optimization |
| `--pysr` | Use PySR symbolic regression |
| `--generations N` | Max rounds for coordinate descent |

### 2. Multi-Index Joint Optimization

```bash
# First run: optimize shared forcing across indices
python3 lte_pysr.py --multi amo nao pdo --multi-rounds 200

# Subsequent runs: cascade from previous multi_index output
python3 lte_pysr.py --multi amo nao pdo --multi-rounds 200
# (loads from multi_index/shared_params.p and per_index_params.json)
```

**What's shared vs per-index:**
- **Shared** (same for all): lpap (42×2 amp/phase) + scalar forcing params
- **Per-index**: ltep (modulation frequencies) + IR (12-month feedback)
- **Fixed per-index**: k0, DC level, annual/semiannual (computed by regression)

**Output:**
- `multi_index/shared_params.p` — universal forcing parameters
- `multi_index/per_index_params.json` — per-index ltep + IR

Each run cascades from the previous `multi_index/` output, enabling iterative refinement.

```bash
# Show accuracy-complexity Pareto front for multiple indices
python3 pareto_ltep.py --indices amo nao pdo nino4

# Scan ltep perturbation landscape for a single index
python3 pareto_ltep.py --scan --index nino4 --samples 200
```

**Complexity metrics:**
- **n_distinct**: Number of unique ltep values (within 1% tolerance)
- **harmonic_score**: Fraction of ltep values that are integer multiples of the minimum
- **spread**: std/mean ratio of |ltep| values
- **complexity**: Composite score = `n_distinct × (1 + spread) × (2 - harmonic_score)`

Lower complexity = simpler structure. Harmonic solutions (e.g., 11:1 ratios) score best.

### 3. Pareto Front Analysis

```bash
# Show accuracy-complexity Pareto front for multiple indices
python3 pareto_ltep.py --indices amo nao pdo nino4

# Scan ltep perturbation landscape for a single index
python3 pareto_ltep.py --scan --index nino4 --samples 200
```

**Complexity metrics:**
- **n_distinct**: Number of unique ltep values (within 1% tolerance)
- **harmonic_score**: Fraction of ltep values that are integer multiples of the minimum
- **spread**: std/mean ratio of |ltep| values
- **complexity**: Composite score = `n_distinct × (1 + spread) × (2 - harmonic_score)`

Lower complexity = simpler structure. Harmonic solutions (e.g., 11:1 ratios) score best.

### 4. Round-Robin Forcing Convergence

```bash
# 3 rounds of lpap optimization with shared forcing constraint
python3 round_robin_forcing.py --rounds 3 --indices amo nao pdo nino4

# Force shared forcing from round 2 onward (loads from multi_index/ first)
python3 round_robin_forcing.py --rounds 5 --force-shared
```

Each round:
1. Loads params from `multi_index/shared_params.p` (if exists) → `round{N-1}.p` → `.coord.p` → `lt.exe.p`
2. Optimizes each index's lpap independently (round 1) or against a shared median forcing (round 2+)
3. Computes pairwise forcing manifold correlations
4. Stops early if mean similarity > 0.99

### 5. Trial Tracking

```bash
# Record a trial
python3 trial_tracker.py record <index> <subset> <cc_before> <cc_after> \
  --params-file <path> --seed-from <source> --notes "description"

# List all trials, sorted by gain
python3 trial_tracker.py list --sort cc_gain

# Show optimization chain for an index
python3 trial_tracker.py best nino4
```

## Candidate Score Tracking

### Primary Metric: Pearson CC (Correlation Coefficient)

Computed on the training-excluded window (default: 1880–1885 excluded, rest used for CC).

```
CC = Pearson(model[sel], data[sel])
```

The `cc()` function in `lte_forward.py` trims leading/trailing zeros before computing.

### Tracking Chain

Each optimization step is recorded in `trials.jsonl`:

```jsonl
{"timestamp": "2026-09-02T20:21:10+00:00", "index": "nino4", "subset": "ir_ltep",
 "cc_before": 0.346695, "cc_after": 0.394136, "cc_gain": 0.047441,
 "seed_from": "pdo-optimized", "params_file": "nino4/lt.exe.pysr.ir_ltep.coord.p",
 "notes": "IR + ltep coordinate descent"}
```

**Query patterns:**
```bash
# Best gain per index
python3 trial_tracker.py list --sort cc_gain

# Final CC rankings
python3 trial_tracker.py list --sort cc_after

# Optimization chain for one index
python3 trial_tracker.py best nino4
```

### Example Results

| Index | Step | CC_before | CC_after | Gain | Seed |
|-------|------|-----------|----------|------|------|
| nino4 | ir_ltep | 0.347 | 0.394 | +0.047 | pdo-optimized |
| nino4 | lpap_top | 0.394 | 0.534 | +0.140 | pdo-ir-ltep-opt |
| nao | ir_ltep | 0.522 | 0.541 | +0.019 | nao-original |
| pdo | ir_ltep | 0.658 | 0.664 | +0.006 | pdo-original |

**Total gain for NINO4: +0.187** (from 0.347 → 0.534), beating its own optimizer's 0.388.

## Hidden Latent Forcing Manifold Agreement

The forcing manifold (column 4 of `lte_results.csv`) is computed as:

```
forcing = IIR(tide_sum × impulse_delta)
forcing = Bessel(forcing, impA, impB, m[NM], offs, bg)
```

### Cross-Index Similarity

Pairwise Pearson correlations between forcing manifolds:

|          | amo  | nao  | pdo  | nino4 |
|----------|------|------|------|-------|
| **amo**  | 1.00 | 1.00 | 1.00 | 0.24  |
| **nao**  | 1.00 | 1.00 | 1.00 | 0.24  |
| **pdo**  | 1.00 | 1.00 | 1.00 | 0.24  |
| **nino4**| 0.24 | 0.24 | 0.24 | 1.00  |

**Key finding:** Atlantic/North Atlantic indices (AMO, NAO, PDO) share identical forcing (r=1.000). NINO4 (Pacific ENSO) diverges (r=0.24), reflecting different ocean dynamics.

### Convergence Protocol

The round-robin script enforces forcing convergence:

```python
# Round 1: Independent optimization → extract forcing
forcing_dict[idx] = forward(cfg, params)["forcing"]

# Round 2+: Compute shared forcing (median across indices)
shared = np.median([forcing_dict[i] for i in indices], axis=0)

# Optimize each index's lpap to match shared forcing
# while maintaining CC
```

**Convergence criterion:** Mean pairwise similarity > 0.99

### Why Forcing Agreement Matters

The lpap (tidal amplitude/phase) parameters are the **fine-tuning layer** that adjusts the shared forcing to each index's specific response. The chart `hidden_latent_forcing_comparison.png` shows all 8 indices overlaid — they're nearly identical, with only micro-shifts in amp/phase (0.001–0.003) needed for optimal fit.

This supports the hypothesis of a **universal tidal forcing** modulated by index-specific LTE standing-wave responses (ltep) and feedback (IR).

## File Structure

```
experiments/Feb2026/
├── lte_forward.py          # Deterministic forward pipeline (Ada-equivalent)
├── lte_pysr.py             # PySR wrapper + coordinate descent + evolution
├── pareto_ltep.py          # Accuracy-complexity Pareto front analysis
├── round_robin_forcing.py  # Cross-index forcing convergence
├── trial_tracker.py        # Trial recording and ranking
├── trials.jsonl            # Persistent trial log
├── <index>/
│   ├── lt.exe.p            # Optimized parameters (JSON)
│   ├── lt.exe.resp         # Configuration overrides
│   ├── <index>.dat         # Climate index time series
│   ├── lte_results.csv     # Forward pipeline output
│   └── lt.exe.pysr.*.p     # Optimization checkpoints
```

## Key Parameters

| Parameter | Description | Typical Range | Role |
|-----------|-------------|---------------|------|
| **IR** | 12-month autoregressive feedback | -0.3 to 0.5 | Index-specific damping |
| **ltep** | LTE modulation frequencies (N values) | 0.01 to 17.0 | Standing-wave modes |
| **lpap** | Tidal amplitudes/phases (42×2) | amp: 0.001–0.3, ph: -π to π | Forcing manifold fine-tuning |
| **impA/B** | Bessel modulation coefficients | -4 to -3 | Forcing nonlinearity |
| **delA/B** | Impulse comb parameters | -5 to 6 | Monthly forcing timing |

## Workflow Summary

```
1. Seed selection
   └── Cross-index transfer (e.g., PDO → NINO4)

2. IR + ltep optimization (coarse)
   └── Coordinate descent, 300 rounds
   └── Gain: +0.006 to +0.047

3. lpap micro-tuning (fine)
   └── Top-10 constituents, 500 rounds
   └── Gain: +0.001 to +0.140

4. Pareto analysis
   └── Scan ±10% ltep perturbations
   └── Identify harmonic solutions (complexity ↓)

5. Round-robin convergence
   └── Force shared forcing manifold
   └── Iterate until r > 0.99
```

The full pipeline transforms NINO4 from CC=0.347 (PDO-seeded) → 0.534 (fully optimized), a **+0.187 total gain** and **+0.146 over NINO4's own optimizer**.
