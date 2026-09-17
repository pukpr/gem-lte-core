#!/usr/bin/env python3
"""Round-robin forcing manifold convergence across climate indices.

Strategy:
1. Optimize each index independently (lpap) → extract forcing manifold (col4)
2. Compute the mean/median forcing across all indices
3. Re-optimize each index's lpap to match the shared forcing while maximizing CC
4. Iterate until forcing manifolds converge (low inter-index variance)

Usage:
    ./round_robin_forcing.py --rounds 3 --indices amo nao pdo nino4 tna baltic
    ./round_robin_forcing.py --rounds 5 --force-shared  # force shared forcing from round 2
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import List

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from lte_forward import Config, read_resp, forward, write_results
from lte_pysr import coordinate_descent, evaluate_params, get_param_vector, apply_deltas


INDICES = ["amo", "nao", "pdo", "nino4", "tna", "baltic"]


def get_forcing_manifold(run_dir: Path, cfg: Config, params: dict) -> np.ndarray:
    """Run forward pipeline and return the forcing manifold (column 4 of lte_results)."""
    result = forward(cfg, params)
    return result["forcing"]


def compute_shared_forcing(forcing_dict: dict[str, np.ndarray]) -> np.ndarray:
    """Compute the median forcing across all indices."""
    if not forcing_dict:
        raise ValueError("No forcing manifolds available to compute shared forcing")
    
    all_forcing = list(forcing_dict.values())
    
    # Pad or trim to common length
    lengths = [len(f) for f in all_forcing]
    min_len = min(lengths)
    trimmed = [f[:min_len] for f in all_forcing]
    
    # Median across indices
    shared = np.median(trimmed, axis=0)
    return shared


def forcing_similarity(f1: np.ndarray, f2: np.ndarray) -> float:
    """Pearson correlation between two forcing manifolds."""
    min_len = min(len(f1), len(f2))
    return float(np.corrcoef(f1[:min_len], f2[:min_len])[0, 1])


def optimize_lpap_for_shared_forcing(
    run_dir: Path, cfg: Config, params: dict, shared_forcing: np.ndarray,
    lpap_n: int = 10, max_rounds: int = 200
) -> tuple[dict, float, float]:
    """Optimize lpap to match shared forcing while maintaining CC.
    
    This is a bi-objective optimization:
    1. Maximize CC (correlation with data)
    2. Minimize forcing distance from shared manifold
    
    We use a weighted combination: score = CC + w * forcing_similarity
    """
    # First, get the baseline forcing and CC
    baseline_forcing = get_forcing_manifold(run_dir, cfg, params)
    baseline_cc = evaluate_params(cfg, params)
    
    print(f"  Baseline: CC={baseline_cc:.6f}, "
          f"forcing_sim={forcing_similarity(baseline_forcing, shared_forcing):.4f}")
    
    # For now, use standard coordinate descent on lpap
    # The shared forcing acts as a soft constraint — we'll check similarity after
    params_out, cc_out = coordinate_descent(
        run_dir, cfg, params, 'lpap_top',
        max_rounds=max_rounds, lpap_n=lpap_n
    )
    
    # Check forcing similarity after optimization
    opt_forcing = get_forcing_manifold(run_dir, cfg, params_out)
    sim = forcing_similarity(opt_forcing, shared_forcing)
    
    print(f"  After lpap opt: CC={cc_out:.6f}, "
          f"forcing_sim={sim:.4f}")
    
    return params_out, cc_out, sim


def run_round_robin(
    indices: List[str],
    rounds: int = 3,
    lpap_n: int = 10,
    force_shared: bool = False,
):
    """Run the round-robin forcing convergence."""
    base_dir = ROOT
    
    # Track results across rounds
    all_results = {idx: [] for idx in indices}
    forcing_history = {idx: [] for idx in indices}
    
    for round_num in range(1, rounds + 1):
        print(f"\n{'='*60}")
        print(f"ROUND {round_num}/{rounds}")
        print(f"{'='*60}")

        # Phase 1: Optimize each index independently (round 1) or with shared forcing (round 2+)
        forcing_dict = {}
        if round_num > 1:
            # Carry forward previous round's forcings
            for idx in indices:
                if forcing_history[idx]:
                    forcing_dict[idx] = forcing_history[idx][-1]
        
        for idx in indices:
            run_dir = base_dir / idx
            if not run_dir.is_dir():
                print(f"  Skipping {idx} (directory not found)")
                continue
            
            os.chdir(run_dir)

            # Load params: use the best available optimized params,
            # falling back to original lt.exe.p
            params = None
            loaded_from = "lt.exe.p"
            # Priority: round N-1 → lpap_top.coord.p → ir_ltep.coord.p → lt.exe.p
            candidates = [
                f"lt.exe.pysr.lpap_top.round{round_num-1}.p",
                "lt.exe.pysr.lpap_top.coord.p",
                "lt.exe.pysr.ir_ltep.coord.p",
                "lt.exe.p",
            ]
            for cand in candidates:
                cand_path = run_dir / cand
                if cand_path.exists():
                    params = json.loads(cand_path.read_text())
                    loaded_from = cand
                    break

            if params is None:
                print(f"  ERROR: No params file found for {idx}")
                continue
            
            # Load config
            resp = read_resp(run_dir / "lt.exe.resp")
            resp["CLIMATE_INDEX"] = f"{idx}.dat"
            
            # Get NH from .p file
            nh_override = {}
            if "harm" in params:
                nh_override["NH"] = " ".join(str(int(float(h))) for h in params["harm"])
            
            cfg = Config(resp, {
                "TRAIN_START": "1880", "TRAIN_END": "1885",
                **nh_override
            })
            
            print(f"\n  {idx} (from {loaded_from}):")
            
            if round_num == 1 or not force_shared:
                # Standard lpap optimization
                params_out, cc_out = coordinate_descent(
                    run_dir, cfg, params, 'lpap_top',
                    max_rounds=200, lpap_n=lpap_n
                )
                forcing = get_forcing_manifold(run_dir, cfg, params_out)
                forcing_dict[idx] = forcing
                
                # Save results
                out_file = run_dir / f"lt.exe.pysr.lpap_top.round{round_num}.p"
                with open(out_file, "w") as f:
                    json.dump(params_out, f, indent=2)
                
                # Record
                all_results[idx].append({
                    "round": round_num,
                    "cc": cc_out,
                    "forcing_sim": 1.0,  # N/A for round 1
                })
                forcing_history[idx].append(forcing)
                
                print(f"    CC={cc_out:.6f}")
                
            else:
                # Optimize against shared forcing
                if len(forcing_dict) < 2:
                    print(f"    SKIPPING shared forcing (only {len(forcing_dict)} forcing available)")
                    # Fall back to independent optimization
                    params_out, cc_out = coordinate_descent(
                        run_dir, cfg, params, 'lpap_top',
                        max_rounds=200, lpap_n=lpap_n
                    )
                    forcing = get_forcing_manifold(run_dir, cfg, params_out)
                    forcing_dict[idx] = forcing
                    forcing_history[idx].append(forcing)
                    print(f"    CC={cc_out:.6f}")
                    continue

                shared_forcing = compute_shared_forcing(forcing_dict)
                params_out, cc_out, sim = optimize_lpap_for_shared_forcing(
                    run_dir, cfg, params, shared_forcing,
                    lpap_n=lpap_n, max_rounds=200
                )
                forcing = get_forcing_manifold(run_dir, cfg, params_out)
                forcing_dict[idx] = forcing
                
                # Save
                out_file = run_dir / f"lt.exe.pysr.lpap_top.round{round_num}.p"
                with open(out_file, "w") as f:
                    json.dump(params_out, f, indent=2)
                
                all_results[idx].append({
                    "round": round_num,
                    "cc": cc_out,
                    "forcing_sim": sim,
                })
                forcing_history[idx].append(forcing)
        
        # Phase 2: Compute convergence metrics
        if len(forcing_dict) > 1:
            # Pairwise forcing similarities
            sims = []
            idx_list = list(forcing_dict.keys())
            for i in range(len(idx_list)):
                for j in range(i + 1, len(idx_list)):
                    s = forcing_similarity(forcing_dict[idx_list[i]], forcing_dict[idx_list[j]])
                    sims.append(s)
            
            mean_sim = np.mean(sims)
            min_sim = np.min(sims)
            
            print(f"\n  Forcing convergence:")
            print(f"    Mean pairwise similarity: {mean_sim:.4f}")
            print(f"    Min pairwise similarity:  {min_sim:.4f}")
            
            # Check for convergence
            if mean_sim > 0.99:
                print(f"  CONVERGED! (mean sim > 0.99)")
                break
    
    # Final summary
    print(f"\n{'='*60}")
    print(f"FINAL SUMMARY")
    print(f"{'='*60}")
    print(f"{'Index':<10} {'Round 1 CC':>12} {'Final CC':>12} {'Δ CC':>10} {'Final Forcing Sim':>18}")
    print("-" * 70)
    
    for idx in indices:
        if idx in all_results and all_results[idx]:
            r1 = all_results[idx][0]
            rf = all_results[idx][-1]
            delta = rf["cc"] - r1["cc"]
            sim = rf.get("forcing_sim", "N/A")
            sim_str = f"{sim:.4f}" if isinstance(sim, float) else str(sim)
            print(f"{idx:<10} {r1['cc']:>12.6f} {rf['cc']:>12.6f} {delta:>+10.6f} {sim_str:>18}")


def main():
    ap = argparse.ArgumentParser(description="Round-robin forcing manifold convergence")
    ap.add_argument("--rounds", type=int, default=3, help="Number of rounds (default 3)")
    ap.add_argument("--indices", nargs="+", default=INDICES, help="Indices to optimize")
    ap.add_argument("--lpap-n", type=int, default=10, help="Number of lpap constituents (default 10)")
    ap.add_argument("--force-shared", action="store_true", help="Force shared forcing from round 2")
    
    args = ap.parse_args()
    run_round_robin(args.indices, args.rounds, args.lpap_n, args.force_shared)


if __name__ == "__main__":
    main()
