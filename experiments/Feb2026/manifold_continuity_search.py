#!/usr/bin/env python3
"""manifold_continuity_search.py -- search for a manifold (tidal factors +
mA/mP/etc integration factors) that makes the REAL DATA's own winding
scalogram show the most continuous, unbroken ridges at a chosen handful of
winding rates M.

Motivation: comparing kN020_E050's pre-1950 and 1950+ fits (see
plot_composite_windings.py) found 9 of 11 fitted M values matching to 4
decimal places between two independently-optimized eras -- "tantalizingly
close". But checking those M values against the full 1856-2023 record's
winding scalogram (winding_transform + winding_rank.continuity, the same
machinery winding_rank.py already uses to rank ridges) shows NONE of them
clear the usual 0.70 continuity bar on the CURRENT manifold -- the best is
only ~0.52. The manifold (how Forcing(t) is built from the tidal factors
via Tide_Sum, then Impulse_Delta + IIR's mA/mP/etc integration, then
Bessel) determines what "winding around Forcing" even means, so a
different manifold could turn some of these into genuinely continuous
ridges -- or reveal that they can't be, on this data.

Design choices driven directly by the two things this was explicitly
asked to account for:

  "it will work on tidal factors and the mA, mP, etc integration factors"
    -- the search perturbs exactly two parameter groups, nothing else:
    lpap's own (amplitude, phase) per tidal constituent ("tidal factors"),
    and the scalar delA/delB/asym/shfT/ma/mp/offs/bg/impA/impB/init
    manifold parameters ("integration factors", named for mA/mP's role in
    IIR's own feedback/lag integration, generalized to their scalar
    siblings). The winding rates M themselves (lpap's own backbone/
    harmonic multipliers) are held FIXED -- moving them would just be
    re-doing the original frequency search, not asking "does a different
    manifold make THESE rates more continuous."

  "computing a full scalogram is expensive, unless it focuses on just a
  few windings" -- the objective never calls winding_rank.rank_series'
  full 0..m_max scan. Each evaluation only runs winding_transform (plus
  its own small AR1-surrogate floor) on a narrow +/-0.15 window around
  each of the (few) target M values -- the same transform winding_rank.py
  itself uses, just restricted in M-extent, which is what actually makes
  each iteration cheap enough to search at all. The floor still has to be
  recomputed every iteration (it depends on the candidate Forcing, not
  just M), so this remains real work per step -- expect slow, incremental
  improvement over many iterations, not a fast convergence.

Objective: mean of winding_rank.continuity(...) across the target M list,
computed for the REAL DATA against the CANDIDATE forcing (never the
Model) -- matching winding_rank.py's own anti-circularity rule: this asks
whether the manifold reveals structure already in the data, not whether a
regression can be made to fit it.

Search: greedy random-descent (propose one small perturbation to one
randomly chosen parameter, keep it only if the objective improves),
mirroring the project's own Ada Walker.Markov accept-if-better rule and
its cosine-annealed step-size schedule, not full simulated annealing.

Usage:
    manifold_continuity_search.py --index kN020_E050_ \\
        --targets 0.2075 1.245 2.2825 3.1125 6.4325 \\
        --seconds 120 --out kN020_E050_/manifold_search_best.json
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ws                                            # noqa: E402
from wavelet_scalogram import standardize            # noqa: E402
from winding_scalogram import winding_transform      # noqa: E402
from winding_rank import continuity as ridge_continuity  # noqa: E402

SCALAR_NAMES = ("delA", "delB", "asym", "shfT", "ma", "mp", "offs", "bg",
                "impA", "impB", "init")


def local_ar1_floor(t: np.ndarray, forcing: np.ndarray, m_grid: np.ndarray,
                     t0_grid: np.ndarray, sigma: float, n_reps: int,
                     rho: float, seed: int) -> np.ndarray:
    """Same AR(1) surrogate convention as winding_rank.ar1_floor, but
    callable with an arbitrary (small) m_grid -- this is what keeps a
    single evaluation cheap: the floor is real work (it depends on the
    CANDIDATE forcing, so it can't be precomputed once), but only over a
    handful of M values' narrow neighborhoods, not the full spectrum."""
    rng = np.random.default_rng(seed)
    n = len(t)
    floor = np.zeros(len(m_grid))
    s = np.sqrt(max(1e-12, 1.0 - rho * rho))
    for _ in range(n_reps):
        e = rng.standard_normal(n)
        z = np.zeros(n)
        for i in range(1, n):
            z[i] = rho * z[i - 1] + s * e[i]
        Gn, _, _ = winding_transform(t, standardize(z), forcing, m_grid, t0_grid, sigma)
        floor += np.mean(np.abs(Gn) ** 2, axis=1)
    return floor / n_reps


def continuity_at(t: np.ndarray, data_std: np.ndarray, forcing: np.ndarray,
                   M: float, sigma: float, t0_grid: np.ndarray, halfband: float,
                   n_reps: int, rho: float, seed: int) -> float:
    m_grid = np.linspace(M - halfband, M + halfband, 31)
    G, _, _ = winding_transform(t, data_std, forcing, m_grid, t0_grid, sigma)
    floor = local_ar1_floor(t, forcing, m_grid, t0_grid, sigma, n_reps, rho, seed)
    lp = np.log2(np.maximum(np.abs(G) ** 2 / np.maximum(floor[:, None], 1e-12), 1e-9))
    return ridge_continuity(lp, m_grid, M)


def objective(t: np.ndarray, data_std: np.ndarray, forcing: np.ndarray,
              targets: list[float], sigma: float, t0_grid: np.ndarray,
              halfband: float, n_reps: int, rho: float, seed: int
              ) -> tuple[float, list[float]]:
    scores = [continuity_at(t, data_std, forcing, M, sigma, t0_grid, halfband,
                             n_reps, rho, seed + i)
              for i, M in enumerate(targets)]
    return float(np.mean(scores)), scores


def rebuild_forcing(dates: np.ndarray, lpap: np.ndarray, periods: np.ndarray,
                     yl: float, manifold: dict, k1: float, idate: float) -> np.ndarray:
    """Same four-stage pipeline as ws.build_forcing_at (Tide_Sum ->
    Impulse_Delta -> IIR -> Bessel), parameterized so every "tidal factor"
    (lpap's own amplitude/phase columns) and "integration factor" (the
    manifold dict) can be perturbed independently. k1 (the Bessel carrier
    / backbone winding rate) is held fixed -- see module docstring."""
    tf = ws.tide_sum(dates, lpap[:, 1:3], periods, yl, 0.0, manifold["shfT"])
    comb = ws.impulse_delta(dates, manifold["delA"], manifold["delB"], manifold["asym"], 12)
    R = ws.iir(tf * comb, lag_a=1.0 - manifold["ma"], lag_c=manifold["mp"],
               init=manifold["init"], start_date=idate, dates=dates)
    return ws.bessel(R, manifold["impA"], manifold["impB"], k1, manifold["offs"], manifold["bg"])


def search(prep: dict, targets: list[float], seconds: float, sigma: float,
           t0_step: float, halfband: float, n_reps: int, rho: float,
           seed: int, report_every: int) -> dict:
    dates = prep["dates"]
    data_std = standardize(prep["data"])
    periods = ws.doodson_periods(prep["year_startup"], prep["year_cand"])
    t0_grid = np.arange(dates[0] + sigma / 2, dates[-1] - sigma / 2 + 1e-9, t0_step)
    k1 = float(np.array(prep["params"]["ltep"], dtype=float)[prep["nm"] - 1])

    manifold0 = {k: float(prep["params"].get(k, 0.0)) for k in SCALAR_NAMES}
    lpap0 = np.array(prep["params"]["lpap"], dtype=float)

    def forcing_of(manifold: dict, lpap: np.ndarray) -> np.ndarray:
        return rebuild_forcing(dates, lpap, periods, prep["yl"], manifold, k1, prep["idate"])

    rng = np.random.default_rng(seed)

    best_manifold, best_lpap = dict(manifold0), lpap0.copy()
    f0 = forcing_of(best_manifold, best_lpap)
    best_score, best_per_m = objective(dates, data_std, f0, targets, sigma, t0_grid,
                                        halfband, n_reps, rho, seed)
    print(f"[{prep['index']}] targets: {targets}")
    print(f"start: mean continuity={best_score:.3f}  "
          f"per-M={[round(c, 2) for c in best_per_m]}")

    n_lpap = len(lpap0)
    t_start = time.time()
    it = 0
    accepted = 0
    while time.time() - t_start < seconds:
        it += 1
        # Cosine-annealed step size, same shape as the Ada search's own
        # Spread_Min + Spread_Max*(1-cos(cycle)) schedule -- oscillates
        # rather than monotonically shrinking, to keep some chance of
        # escaping a local plateau in this discrete, saturating objective.
        spread = 0.03 + 0.20 * (1.0 - np.cos(it / 40.0))
        cand_manifold, cand_lpap = dict(best_manifold), best_lpap.copy()
        if rng.random() < 0.35:
            name = SCALAR_NAMES[rng.integers(len(SCALAR_NAMES))]
            scale = abs(cand_manifold[name]) if cand_manifold[name] != 0.0 else 1.0
            cand_manifold[name] += rng.normal(0.0, spread * scale * 0.3 + 1e-4)
            what = name
        else:
            row = int(rng.integers(n_lpap))
            cand_lpap[row, 1] += rng.normal(0.0, spread * 0.05)
            cand_lpap[row, 2] += rng.normal(0.0, spread * 0.5)
            what = f"lpap[{row}]"

        f = forcing_of(cand_manifold, cand_lpap)
        score, per_m = objective(dates, data_std, f, targets, sigma, t0_grid,
                                  halfband, n_reps, rho, seed + it)
        if score > best_score:
            best_score, best_per_m = score, per_m
            best_manifold, best_lpap = cand_manifold, cand_lpap
            accepted += 1
            print(f"[{it:5d} t={time.time() - t_start:6.1f}s] IMPROVED "
                  f"mean_cont={best_score:.3f}  via {what}  "
                  f"per-M={[round(c, 2) for c in best_per_m]}")
        elif it % report_every == 0:
            print(f"[{it:5d} t={time.time() - t_start:6.1f}s] best={best_score:.3f}  "
                  f"accepted={accepted}/{it}")

    return dict(index=prep["index"], targets=targets, seconds=seconds, iterations=it,
                accepted=accepted, start_score=None, best_score=best_score,
                best_per_m=best_per_m, manifold=best_manifold,
                lpap=best_lpap.tolist(), k1=k1)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--index", required=True, help="index directory (needs lt.exe.p + .dat)")
    ap.add_argument("--root", default=None)
    ap.add_argument("--targets", type=float, nargs="+", required=True,
                     help="the (few) winding rates M to optimize continuity for")
    ap.add_argument("--seconds", type=float, default=120.0,
                     help="search time budget (default 120s)")
    ap.add_argument("--sigma", type=float, default=15.0, help="window width, years")
    ap.add_argument("--t0-step", type=float, default=5.0)
    ap.add_argument("--halfband", type=float, default=0.15,
                     help="per-target local m_grid half-width (keeps each "
                          "evaluation cheap -- this is the 'few windings' "
                          "restriction, not a full scalogram)")
    ap.add_argument("--n-reps", type=int, default=8,
                     help="AR1 surrogates per floor estimate (lower = faster "
                          "but noisier objective; winding_rank.py's own "
                          "default is 24 for final reporting)")
    ap.add_argument("--rho", type=float, default=0.97)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--report-every", type=int, default=20)
    ap.add_argument("--out", type=Path, default=None,
                     help="save the best manifold+lpap found to this JSON path")
    args = ap.parse_args()

    prep = ws.load_index(ws.find_index_root(args.index, args.root), args.index)
    result = search(prep, args.targets, args.seconds, args.sigma, args.t0_step,
                     args.halfband, args.n_reps, args.rho, args.seed, args.report_every)

    print(f"\nfinal: mean continuity {result['best_score']:.3f} after "
          f"{result['iterations']} iterations ({result['accepted']} accepted) "
          f"in {args.seconds:.0f}s")
    print(f"per-M: {dict(zip(result['targets'], [round(c, 2) for c in result['best_per_m']]))}")

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(result, indent=2))
        print(f"saved best manifold to {args.out}")


if __name__ == "__main__":
    main()
