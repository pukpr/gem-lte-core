#!/usr/bin/env python3
"""PySR wrapper for LTE parameter evolution.

Uses PySR to discover symbolic expressions that mutate key parameters
of the LTE forward model. The fitness function runs lte_forward.forward()
and returns the CC (correlation coefficient) to maximize.

Strategy: evolve a small subset of parameters first (ltep = 6 values),
then expand to scalar B params and lpap amplitudes/phases.

Usage:
    ./lte_pysr.py amo                  # run PySR search on ltep
    ./lte_pysr.py amo --subset scalars # evolve scalar B params instead
    ./lte_pysr.py amo --evaluate "expression_string"  # evaluate one expression
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Callable

import numpy as np

# Add parent dir so we can import lte_forward
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from lte_forward import Config, read_resp, forward


# ---------------------------------------------------------------------------
# Parameter subsets — which params to evolve
# ---------------------------------------------------------------------------

# ltep: LTE modulation frequencies (6 values for AMO)
# These control the standing-wave response modes.
LTEP_KEYS = ["ltep"]

# Scalar B params: the most impactful forcing-manifold knobs
SCALAR_KEYS = ["impA", "impB", "delA", "delB", "asym", "ma", "mp", "shfT", "init"]

# lpap amplitudes and phases (42 constituents × 2 = 84 values, but we
# pick the top-N by amplitude to keep the search tractable)
LPAP_TOP_N = 10  # evolve top-10 constituents' amp+phase

# IR + ltep: the NAO-relevant combo (IR feedback + modulation freqs)
IR_KEY = ["IR"]


def get_param_vector(params: dict, subset: str, lpap_n: int = 10) -> np.ndarray:
    """Extract the flat parameter vector to evolve."""
    if subset == "ltep":
        return np.array(params["ltep"], dtype=float)
    elif subset == "scalars":
        return np.array([float(params[k]) for k in SCALAR_KEYS])
    elif subset == "lpap_top":
        lpap = np.array(params["lpap"], dtype=float)
        amps = np.abs(lpap[:, 1])
        top_idx = np.argsort(amps)[-lpap_n:]  # top-N by |amplitude|
        vals = []
        for idx in sorted(top_idx):
            vals.extend([lpap[idx, 1], lpap[idx, 2]])  # amp, phase
        return np.array(vals, dtype=float)
    elif subset == "lpap":
        # All 42 constituents: amp+phase for each
        lpap = np.array(params["lpap"], dtype=float)
        vals = []
        for idx in range(len(lpap)):
            vals.extend([lpap[idx, 1], lpap[idx, 2]])
        return np.array(vals, dtype=float)
    elif subset == "ir":
        return np.array([float(params["IR"])], dtype=float)
    elif subset == "ir_ltep":
        return np.concatenate([[float(params["IR"])], params["ltep"]])
    else:
        raise ValueError(f"Unknown subset: {subset}")


def apply_deltas(params: dict, deltas: np.ndarray, subset: str, lpap_n: int = 10) -> dict:
    """Apply evolved deltas to a copy of params."""
    mutated = json.loads(json.dumps(params))  # deep copy

    if subset == "ltep":
        old = np.array(mutated["ltep"], dtype=float)
        mutated["ltep"] = (old + deltas).tolist()
    elif subset == "scalars":
        for k, d in zip(SCALAR_KEYS, deltas):
            mutated[k] = float(mutated[k] + d)
    elif subset == "lpap_top":
        lpap = np.array(mutated["lpap"], dtype=float)
        amps = np.abs(lpap[:, 1])
        top_idx = np.argsort(amps)[-lpap_n:]
        sorted_idx = sorted(top_idx)
        for i, idx in enumerate(sorted_idx):
            lpap[idx, 1] += deltas[2 * i]       # amp delta
            lpap[idx, 2] += deltas[2 * i + 1]   # phase delta
        mutated["lpap"] = lpap.tolist()
    elif subset == "lpap":
        # All 42 constituents
        lpap = np.array(mutated["lpap"], dtype=float)
        for i in range(len(lpap)):
            lpap[i, 1] += deltas[2 * i]
            lpap[i, 2] += deltas[2 * i + 1]
        mutated["lpap"] = lpap.tolist()
    elif subset == "ir":
        mutated["IR"] = float(mutated["IR"] + deltas[0])
    elif subset == "ir_ltep":
        mutated["IR"] = float(mutated["IR"] + deltas[0])
        old = np.array(mutated["ltep"], dtype=float)
        mutated["ltep"] = (old + deltas[1:]).tolist()
    else:
        raise ValueError(f"Unknown subset: {subset}")

    return mutated


# ---------------------------------------------------------------------------
# Feature generation — inputs to the symbolic expression
# ---------------------------------------------------------------------------

def make_features(params: dict, subset: str, lpap_n: int = 10) -> np.ndarray:
    """Build feature matrix (n_params × n_features) for the expression.

    Each row corresponds to one parameter being evolved. Features include
    the parameter's base value, its index, and global context.
    """
    vec = get_param_vector(params, subset)
    n = len(vec)

    # Features per parameter:
    #   - base value (normalized)
    #   - index / n
    #   - index² / n²  (captures position-dependent patterns)
    #   - global: mean of all params, std of all params
    mean_val = float(np.mean(vec))
    std_val = float(np.std(vec)) + 1e-30

    X = np.zeros((n, 4))
    for i in range(n):
        X[i, 0] = vec[i] / (std_val)           # normalized value
        X[i, 1] = i / n                          # normalized index
        X[i, 2] = (i / n) ** 2                   # index²
        X[i, 3] = mean_val / std_val             # global context (same for all)
    return X


# ---------------------------------------------------------------------------
# Fitness function — runs forward pipeline, returns CC
# ---------------------------------------------------------------------------

def make_fitness(run_dir: Path, cfg: Config, base_params: dict, subset: str):
    """Returns a callable: expression_callable(X) -> np.ndarray of CC values.

    PySR expects a function that takes feature matrix X and returns
    predictions y. We repurpose this: for each row of X (one param),
    the expression produces a delta, we apply all deltas, run forward(),
    and return CC broadcast to each row.
    """

    def fitness_func(X: np.ndarray) -> np.ndarray:
        """X has shape (n_params, n_features).
        The PySR expression maps each row to a scalar delta.
        We collect all deltas, mutate params, run forward, return CC."""
        n = X.shape[0]
        # PySR evaluates the expression row-by-row through its internal
        # machinery. Here we just return a placeholder — the real evaluation
        # happens in evaluate_population below.
        return np.zeros(n)

    return fitness_func


def evaluate_params(cfg: Config, mutated_params: dict) -> float:
    """Run forward pipeline, return CC (the metric to maximize)."""
    try:
        result = forward(cfg, mutated_params)
        return float(result["cc"])
    except Exception as e:
        import sys
        print(f"evaluate_params ERROR: {e}", file=sys.stderr)
        return -1.0  # penalty for invalid parameter sets


# ---------------------------------------------------------------------------
# PySR search wrapper — direct fitness maximization
# ---------------------------------------------------------------------------

def run_pysr_search(run_dir: Path, cfg: Config, base_params: dict,
                    subset: str, niterations: int = 50,
                    population_size: int = 100, ncyclesperiteration: int = 500,
                    lpap_n: int = 10):
    """Run PySR to evolve parameter-delta expressions.

    Strategy: generate a dataset of (features, delta) pairs where each
    delta was randomly sampled. Then train PySR to predict the delta that
    maximizes CC. We do this by training on the top-K performers only,
    teaching PySR "what good deltas look like."

    Alternative (better): use PySR's ability to accept a custom loss that
    directly maximizes CC. Here we use the regression approach but on a
    carefully constructed dataset.
    """
    from pysr import PySRRegressor

    vec = get_param_vector(base_params, subset, lpap_n)
    n_params = len(vec)
    base_cc = evaluate_params(cfg, base_params)

    # Generate exploration data
    n_samples = max(200, n_params * 20)
    rng = np.random.RandomState(42)

    print(f"Generating {n_samples} random parameter samples...")
    cc_scores = []
    delta_sets = []

    for i in range(n_samples):
        if i % 50 == 0:
            print(f"  sample {i}/{n_samples} ...")
        # Small perturbations (0.5% of base or 0.0005)
        scale = np.maximum(np.abs(vec) * 0.005, 0.0005)
        deltas = rng.normal(0, 1, size=n_params) * scale

        mutated = apply_deltas(base_params, deltas, subset, lpap_n)
        cc = evaluate_params(cfg, mutated)
        cc_scores.append(cc)
        delta_sets.append(deltas)

    cc_arr = np.array(cc_scores)
    delta_arr = np.array(delta_sets)
    features_base = make_features(base_params, subset, lpap_n)

    print(f"CC range: {cc_arr.min():.6f} .. {cc_arr.max():.6f} "
          f"(base={base_cc:.6f})")
    print(f"Top 10 CC: {np.sort(cc_arr)[-10:]}")

    # Build training data: for each param, features → delta
    # We only use samples that improved over baseline (or are close to it)
    threshold = base_cc * 0.99  # within 1% of baseline
    good_mask = cc_arr >= threshold

    if good_mask.sum() == 0:
        # No samples beat baseline — use top 5% anyway
        threshold = np.percentile(cc_arr, 95)
        good_mask = cc_arr >= threshold
        print(f"No samples beat baseline, using top 5% (CC >= {threshold:.6f})")
    else:
        print(f"{good_mask.sum()} samples at or above baseline*0.99")

    X_rows = []
    y_rows = []
    for i in range(n_samples):
        if good_mask[i]:
            for j in range(n_params):
                # Features: param index, base value, normalized position
                row = [
                    j / n_params,                      # normalized index
                    features_base[j, 0],               # normalized value
                    vec[j],                            # raw base value
                    float(cc_arr[i]),                  # achieved CC
                ]
                X_rows.append(row)
                y_rows.append(delta_arr[i, j])

    X_train = np.array(X_rows)
    y_train = np.array(y_rows)

    feature_names = ["param_idx", "norm_val", "raw_val", "cc"]

    print(f"Training on {len(X_train)} samples from good performers")

    # Run PySR
    print("Starting PySR search...")
    model = PySRRegressor(
        niterations=niterations,
        population_size=population_size,
        ncycles_per_iteration=ncyclesperiteration,
        binary_operators=["+", "-", "*"],
        unary_operators=["sin", "cos", "abs", "sqrt"],
        maxdepth=5,
        maxsize=15,
        populations=8,
        parsimony=0.005,
        progress=True,
        verbosity=0,
    )

    model.fit(X_train, y_train)

    print("\n=== Top PySR Expressions ===")
    print(model)

    # Evaluate best expression
    best_expr = model.get_best()["equation"]
    print(f"\nBest expression: {best_expr}")
    cc = evaluate_single_expression(run_dir, cfg, base_params, subset, best_expr)

    return model


# ---------------------------------------------------------------------------
# Single expression evaluation
# ---------------------------------------------------------------------------

def evaluate_single_expression(run_dir: Path, cfg: Config, base_params: dict,
                               subset: str, expr: str, lpap_n: int = 10):
    """Evaluate a single symbolic expression for parameter mutation."""
    import sympy
    from sympy import sin, cos, exp, log, Abs, sqrt

    vec = get_param_vector(base_params, subset, lpap_n)
    n = len(vec)
    features = make_features(base_params, subset, lpap_n)

    # Parse the expression (PySR format uses sympy-compatible strings)
    # Expected variables: x0, x1, x2, x3 (the 4 feature columns)
    try:
        # Create a lambda from the expression string
        x0, x1, x2, x3 = sympy.symbols("x0 x1 x2 x3")
        # Try to evaluate it
        expr_obj = eval(expr, {"sin": sin, "cos": cos, "exp": exp,
                               "log": log, "Abs": Abs, "sqrt": sqrt,
                               "x0": x0, "x1": x1, "x2": x2, "x3": x3})
        f = sympy.lambdify([x0, x1, x2, x3], expr_obj, modules="numpy")

        # Compute deltas
        deltas = np.array([f(features[i, 0], features[i, 1],
                             features[i, 2], features[i, 3])
                          for i in range(n)])

        mutated = apply_deltas(base_params, deltas, subset, lpap_n)
        cc = evaluate_params(cfg, mutated)

        print(f"Expression: {expr}")
        print(f"Deltas: {np.round(deltas, 6)}")
        print(f"Base CC: {evaluate_params(cfg, base_params):.6f}")
        print(f"Mutated CC: {cc:.6f}")
        return cc
    except Exception as e:
        print(f"Error evaluating expression: {e}")
        return -1.0


# ---------------------------------------------------------------------------
# Iterative evolution loop (alternative to PySR, simpler)
# ---------------------------------------------------------------------------

def iterative_evolution(run_dir: Path, cfg: Config, base_params: dict,
                        subset: str, generations: int = 100,
                        population_size: int = 50, mutation_rate: float = 0.1,
                        lpap_n: int = 10):
    """Simple evolutionary optimization without PySR."""
    rng = np.random.RandomState(42)
    vec = get_param_vector(base_params, subset, lpap_n)
    n_params = len(vec)

    base_cc = evaluate_params(cfg, base_params)
    print(f"Base CC: {base_cc:.6f} ({n_params} parameters)")

    # Initialize population: very small perturbations around base
    population = []
    for _ in range(population_size):
        scale = np.maximum(np.abs(vec) * 0.01, 0.001)
        deltas = rng.normal(0, scale, size=n_params)
        mutated = apply_deltas(base_params, deltas, subset, lpap_n)
        cc = evaluate_params(cfg, mutated)
        population.append((cc, deltas, mutated))

    population.sort(key=lambda x: x[0], reverse=True)

    # Track best ever (not just last gen's best)
    best_ever_cc, best_ever_deltas, best_ever_params = population[0]
    current_mutation_rate = mutation_rate

    print(f"{'Gen':>4} {'Best CC':>10} {'Mean CC':>10} {'Delta':>10} {'MutRate':>8}")

    for gen in range(generations):
        # Select top 20%
        n_elite = max(population_size // 5, 2)
        elites = population[:n_elite]

        # Generate new population from elites + mutation
        new_population = list(elites)

        while len(new_population) < population_size:
            # Pick random elite
            elite_cc, elite_deltas, elite_params = elites[rng.randint(len(elites))]

            # Mutate with adaptive rate
            scale = np.maximum(np.abs(vec) * current_mutation_rate, 0.0001)
            new_deltas = elite_deltas + rng.normal(0, scale, size=n_params)
            mutated = apply_deltas(base_params, new_deltas, subset, lpap_n)
            cc = evaluate_params(cfg, mutated)
            new_population.append((cc, new_deltas, mutated))

        population = sorted(new_population, key=lambda x: x[0], reverse=True)

        best_cc = population[0][0]
        mean_cc = np.mean([p[0] for p in population])
        improvement = best_cc - base_cc

        # Track best ever
        if best_cc > best_ever_cc:
            best_ever_cc, best_ever_deltas, best_ever_params = population[0]
            current_mutation_rate = min(mutation_rate, current_mutation_rate * 1.2)
        else:
            current_mutation_rate *= 0.97  # shrink if no improvement

        current_mutation_rate = max(current_mutation_rate, 1e-6)

        if gen % 10 == 0 or gen == generations - 1:
            print(f"{gen:>4} {best_cc:>10.6f} {mean_cc:>10.6f} "
                  f"{improvement:>+10.6f} {current_mutation_rate:>8.2e}")

    # Report best ever
    print(f"\nBest ever CC: {best_ever_cc:.6f} (improvement: {best_ever_cc - base_cc:+.6f})")
    print(f"Best deltas: {np.round(best_ever_deltas, 6)}")

    # Save best
    out_path = run_dir / f"lt.exe.pysr.{subset}.p"
    with open(out_path, "w") as f:
        json.dump(best_ever_params, f, indent=2)
    print(f"Saved best params to {out_path}")

    # Verify
    verify_cc = evaluate_params(cfg, best_ever_params)
    print(f"Verified CC: {verify_cc:.6f}")

    return best_ever_params, best_ever_cc


# ---------------------------------------------------------------------------
# Coordinate-wise hill climbing — one param at a time
# ---------------------------------------------------------------------------

def coordinate_descent(run_dir: Path, cfg: Config, base_params: dict,
                       subset: str, max_rounds: int = 20,
                       step_sizes: list = None, lpap_n: int = 10):
    """Hill-climb one parameter at a time."""
    vec = get_param_vector(base_params, subset, lpap_n)
    n_params = len(vec)
    base_cc = evaluate_params(cfg, base_params)

    print(f"Base CC: {base_cc:.6f} ({n_params} parameters)")
    print(f"{'Round':>6} {'Param':>6} {'Step':>12} {'CC':>10} {'Delta':>10}")

    current_params = dict(base_params)
    current_cc = base_cc

    if step_sizes is None:
        # Default: 1% of each parameter's magnitude, min 0.001
        step_sizes = np.maximum(np.abs(vec) * 0.01, 0.001).tolist()

    improved = True
    round_num = 0
    while improved and round_num < max_rounds:
        improved = False
        for j in range(n_params):
            if round_num >= max_rounds:
                break
            round_num += 1

            step = step_sizes[j]
            # Try positive direction
            pos_deltas = np.zeros(n_params)
            pos_deltas[j] = step
            pos_mutated = apply_deltas(current_params, pos_deltas, subset, lpap_n)
            pos_cc = evaluate_params(cfg, pos_mutated)

            # Try negative direction
            neg_deltas = np.zeros(n_params)
            neg_deltas[j] = -step
            neg_mutated = apply_deltas(current_params, neg_deltas, subset, lpap_n)
            neg_cc = evaluate_params(cfg, neg_mutated)

            if pos_cc > current_cc and pos_cc >= neg_cc:
                current_params = pos_mutated
                current_cc = pos_cc
                improved = True
                print(f"{round_num:>6} {j:>6} {step:>12.6f} {pos_cc:>10.6f} "
                      f"{step:>10.6f}")
            elif neg_cc > current_cc:
                current_params = neg_mutated
                current_cc = neg_cc
                improved = True
                print(f"{round_num:>6} {j:>6} {step:>12.6f} {neg_cc:>10.6f} "
                      f"{-step:>10.6f}")

            # If neither direction improves, try smaller step
            if not improved:
                step_sizes[j] *= 0.5

    print(f"\nFinal CC: {current_cc:.6f} (improvement: {current_cc - base_cc:+.6f})")

    # Save
    out_path = run_dir / f"lt.exe.pysr.{subset}.coord.p"
    with open(out_path, "w") as f:
        json.dump(current_params, f, indent=2)
    print(f"Saved to {out_path}")

    return current_params, current_cc


# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description="PySR parameter evolution for LTE model.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("index", nargs="?", default="amo",
                    help="index subdirectory (default amo)")
    ap.add_argument("--subset", default="ltep",
                    choices=["ltep", "scalars", "lpap_top", "lpap", "ir", "ir_ltep"],
                    help="parameter subset to evolve")
    ap.add_argument("--lpap-n", type=int, default=10,
                    help="number of top lpap constituents to evolve (default 10)")
    ap.add_argument("--evaluate", type=str, default=None,
                    help="evaluate a single expression string")
    ap.add_argument("--evolve", action="store_true",
                    help="run simple evolutionary optimization")
    ap.add_argument("--coordinate", action="store_true",
                    help="run coordinate-wise hill climbing (one param at a time)")
    ap.add_argument("--pysr", action="store_true",
                    help="run PySR symbolic regression")
    ap.add_argument("--generations", type=int, default=100,
                    help="generations for evolution (default 100)")
    ap.add_argument("--population", type=int, default=50,
                    help="population size (default 50)")
    ap.add_argument("--niterations", type=int, default=50,
                    help="PySR iterations (default 50)")
    ap.add_argument("--cv", nargs=2, metavar=("START", "END"),
                    default=("1880", "1885"))
    args = ap.parse_args()

    run_dir = (ROOT / args.index).resolve()
    if not run_dir.is_dir():
        print(f"error: {run_dir} not found", file=sys.stderr)
        return 2
    os.chdir(run_dir)

    # Load params
    params = json.loads((run_dir / "lt.exe.p").read_text())

    # Build config
    resp = read_resp(run_dir / "lt.exe.resp")
    resp.setdefault("CLIMATE_INDEX", f"{args.index}.dat")
    overrides = {"TRAIN_START": args.cv[0], "TRAIN_END": args.cv[1]}

    # Sync NH from .p file's harm field (resp often has stale NH values)
    if "harm" in params:
        overrides["NH"] = " ".join(str(int(float(h))) for h in params["harm"])
    cfg = Config(resp, overrides)

    # Report baseline
    try:
        base_cc = evaluate_params(cfg, params)
    except Exception as e:
        import traceback
        print(f"ERROR computing baseline: {e}", file=sys.stderr)
        traceback.print_exc()
        return 2
    vec = get_param_vector(params, args.subset, args.lpap_n)
    print(f"Index: {args.index}, Subset: {args.subset} (lpap_n={args.lpap_n})")
    print(f"Parameters to evolve: {len(vec)}")
    print(f"Baseline CC: {base_cc:.6f}")
    print(f"CWD: {os.getcwd()}")
    print(f"Data file: {(run_dir / f'{args.index}.dat').exists()}")

    if args.evaluate:
        cc = evaluate_single_expression(run_dir, cfg, params, args.subset,
                                        args.evaluate, args.lpap_n)
        return 0 if cc > base_cc else 1

    if args.evolve:
        iterative_evolution(run_dir, cfg, params, args.subset,
                           generations=args.generations,
                           population_size=args.population,
                           lpap_n=args.lpap_n)
        return 0

    if args.coordinate:
        coordinate_descent(run_dir, cfg, params, args.subset,
                          max_rounds=args.generations,
                          lpap_n=args.lpap_n)
        return 0

    if args.pysr:
        run_pysr_search(run_dir, cfg, params, args.subset,
                       niterations=args.niterations,
                       lpap_n=args.lpap_n)
        return 0

    # Default: run evolution
    print("\nNo mode specified, running evolutionary optimization.")
    print("Use --pysr for symbolic regression or --evolve for simple evolution.")
    iterative_evolution(run_dir, cfg, params, args.subset,
                       generations=args.generations,
                       population_size=args.population,
                       lpap_n=args.lpap_n)
    return 0


if __name__ == "__main__":
    sys.exit(main())
