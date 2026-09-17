#!/usr/bin/env python3
"""Pareto front analysis for ltep complexity vs CC accuracy.

Measures ltep "complexity" via:
1. Harmonic content: how many ltep values are integer multiples of a base frequency
2. Value spread: std/mean ratio (tighter = stiffer)
3. Distinct clusters: number of unique values (rounded to tolerance)

Usage:
    ./pareto_ltep.py --indices amo nao pdo nino4
    ./pareto_ltep.py --scan --index nino4  # scan ltep perturbations
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from lte_forward import Config, read_resp, forward
from lte_pysr import evaluate_params


def ltep_complexity(ltep: list[float], tol: float = 0.01) -> dict:
    """Measure ltep complexity.
    
    Returns:
        n_distinct: number of unique values (within tol)
        harmonic_score: fraction of values that are integer multiples of the min
        spread: std/mean ratio (higher = more spread out)
        min_val: smallest absolute value (the "base frequency")
    """
    vals = np.array([abs(v) for v in ltep if abs(v) > 1e-10])
    if len(vals) == 0:
        return {"n_distinct": 0, "harmonic_score": 0, "spread": 0, "min_val": 0}
    
    # Distinct clusters
    sorted_vals = np.sort(vals)
    n_distinct = 1
    for i in range(1, len(sorted_vals)):
        if sorted_vals[i] - sorted_vals[i-1] > tol * sorted_vals[i-1]:
            n_distinct += 1
    
    # Harmonic score: fraction that are near-integer multiples of the minimum
    base = vals.min()
    multiples = vals / base
    harmonic_count = sum(1 for m in multiples if abs(m - round(m)) < tol)
    harmonic_score = harmonic_count / len(multiples)
    
    # Spread
    spread = vals.std() / vals.mean() if vals.mean() > 0 else 0
    
    return {
        "n_distinct": n_distinct,
        "harmonic_score": round(harmonic_score, 3),
        "spread": round(spread, 3),
        "min_val": round(base, 6),
        "n_terms": len(vals),
    }


def evaluate_index(idx: str, params: dict, label: str = "") -> dict:
    """Evaluate an index's ltep complexity and CC."""
    ltep = params.get("ltep", [])
    harm = params.get("harm", [])
    
    # Compute CC
    run_dir = ROOT / idx
    os.chdir(run_dir)
    resp = read_resp(run_dir / "lt.exe.resp")
    resp["CLIMATE_INDEX"] = f"{idx}.dat"
    
    nh = ""
    if "harm" in params:
        nh = " ".join(str(int(float(h))) for h in params["harm"])
    
    cfg = Config(resp, {"TRAIN_START": "1880", "TRAIN_END": "1885", "NH": nh})
    cc = evaluate_params(cfg, params)
    
    # Complexity
    comp = ltep_complexity(ltep)
    
    # Combined complexity metric: lower is simpler
    # Penalize high spread, reward harmonics, penalize many distinct values
    complexity = comp["n_distinct"] * (1 + comp["spread"]) * (2 - comp["harmonic_score"])
    
    return {
        "index": idx,
        "label": label,
        "cc": round(cc, 6),
        "ltep": ltep,
        "harm": harm,
        "complexity": round(complexity, 3),
        **comp,
    }


def print_pareto_front(results: list[dict]):
    """Print Pareto front: high CC, low complexity."""
    # Sort by CC descending
    results.sort(key=lambda r: r["cc"], reverse=True)
    
    print(f"{'Index':<10} {'Label':<20} {'CC':>10} {'Complexity':>12} "
          f"{'Distinct':>10} {'Harmonic':>10} {'Spread':>8} {'Min':>10}")
    print("-" * 95)
    
    # Track Pareto frontier
    best_complexity = float("inf")
    pareto = []
    for r in results:
        if r["complexity"] < best_complexity:
            best_complexity = r["complexity"]
            pareto.append(r)
            marker = "  <-- Pareto"
        else:
            marker = ""
        
        print(f"{r['index']:<10} {r['label']:<20} {r['cc']:>10.6f} "
              f"{r['complexity']:>12.3f} {r['n_distinct']:>10} "
              f"{r['harmonic_score']:>10.3f} {r['spread']:>8.3f} "
              f"{r['min_val']:>10.6f}{marker}")
    
    print(f"\nPareto front: {len(pareto)} solutions")
    return pareto


def scan_ltep_perturbations(idx: str, base_params: dict, n_samples: int = 100):
    """Scan random ltep perturbations to map the accuracy-complexity landscape."""
    run_dir = ROOT / idx
    os.chdir(run_dir)
    resp = read_resp(run_dir / "lt.exe.resp")
    resp["CLIMATE_INDEX"] = f"{idx}.dat"
    nh = ""
    if "harm" in base_params:
        nh = " ".join(str(int(float(h))) for h in base_params["harm"])
    cfg = Config(resp, {"TRAIN_START": "1880", "TRAIN_END": "1885", "NH": nh})
    
    base_ltep = np.array(base_params["ltep"], dtype=float)
    base_cc = evaluate_params(cfg, base_params)
    base_comp = ltep_complexity(base_params["ltep"])
    
    print(f"Base: CC={base_cc:.6f}, complexity={base_comp}")
    print(f"ltep: {base_ltep}")
    print(f"\nScanning {n_samples} perturbations...")
    
    results = []
    rng = np.random.RandomState(42)
    
    for i in range(n_samples):
        # Perturb ltep by ±10%
        scale = np.abs(base_ltep) * 0.1 + 0.01
        perturbed = base_ltep + rng.normal(0, scale, size=len(base_ltep))
        
        params = dict(base_params)
        params["ltep"] = perturbed.tolist()
        
        cc = evaluate_params(cfg, params)
        comp = ltep_complexity(perturbed.tolist())
        complexity = comp["n_distinct"] * (1 + comp["spread"]) * (2 - comp["harmonic_score"])
        
        results.append({
            "index": idx,
            "label": f"perturb_{i}",
            "cc": round(cc, 6),
            "complexity": round(complexity, 3),
            **comp,
            "ltep": perturbed.tolist(),
        })
    
    # Add base
    results.append({
        "index": idx,
        "label": "base",
        "cc": round(base_cc, 6),
        "complexity": round(ltep_complexity(base_params["ltep"])["n_distinct"] * 
                          (1 + ltep_complexity(base_params["ltep"])["spread"]) * 
                          (2 - ltep_complexity(base_params["ltep"])["harmonic_score"]), 3),
        **ltep_complexity(base_params["ltep"]),
        "ltep": base_params["ltep"],
    })
    
    return results


def main():
    ap = argparse.ArgumentParser(description="Pareto front analysis for ltep")
    ap.add_argument("--indices", nargs="+", default=["amo", "nao", "pdo", "nino4"])
    ap.add_argument("--scan", action="store_true", help="Scan ltep perturbations")
    ap.add_argument("--index", default=None, help="Index to scan (with --scan)")
    ap.add_argument("--samples", type=int, default=100, help="Perturbation samples")
    
    args = ap.parse_args()
    
    if args.scan:
        idx = args.index or args.indices[0]
        run_dir = ROOT / idx
        params = json.loads((run_dir / "lt.exe.p").read_text())
        results = scan_ltep_perturbations(idx, params, args.samples)
        print_pareto_front(results)
    else:
        # Evaluate each index's current best params
        results = []
        for idx in args.indices:
            run_dir = ROOT / idx
            # Try optimized params first, fall back to original
            for suffix in ["lpap_top.coord.p", "ir_ltep.coord.p", "lt.exe.p"]:
                pfile = run_dir / suffix
                if pfile.exists():
                    params = json.loads(pfile.read_text())
                    label = suffix.replace(".coord.p", "").replace("lt.exe.p", "original")
                    break
            
            r = evaluate_index(idx, params, label)
            results.append(r)
        
        print_pareto_front(results)


if __name__ == "__main__":
    main()
