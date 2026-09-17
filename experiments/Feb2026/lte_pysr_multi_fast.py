#!/usr/bin/env python3
"""lte_pysr_multi_fast.py — performance-focused clone of lte_pysr.py's
--multi coordinate descent (coordinate_descent_multi / evaluate_multi_index).

`--multi` runs correctly but slowly at high --multi-rounds (e.g. 10000):
every trial in Phase 1/2 only changes ONE index's per-index params, yet the
original evaluate_multi_index recomputes forward() — including re-reading
and re-parsing that index's .dat file from disk — for ALL indices on every
call, including the redundant "recompute just to print stats" call right
after each accepted move. This clone adds a value-signature cache
(IndexCache) keyed on the actual shared_params/per-index values that feed
forward() for a given index: if neither has changed since the index was
last evaluated, the cached (cc, forcing) is reused instead of recomputing
(and re-reading its .dat file) — numerically identical results, no
optimizer-behavior change, just skipped duplicate work. It also bakes each
index's CLIMATE_INDEX path as absolute at setup, removing the os.chdir()
in/out of the index directory on every single evaluation. The baked-in
per-index Config still skips a rebuild on every call for the parts that
never change (CLIMATE_INDEX, TRAIN_START/END); only its NH override is
rebuilt per trial (via _cfg_with_harm), since coordinate_descent_multi's
harmonic-search phase does mutate "harm" per index to maximize fit.

Only the --multi path is touched; --subset/--evolve/--coordinate/--pysr are
unchanged from lte_pysr.py (kept only so this file still runs standalone;
prefer lte_pysr.py for those modes).

Usage (same CLI as lte_pysr.py):
    ./lte_pysr_multi_fast.py --multi amo tna pdo --multi-rounds 10000
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

# Non-lpap manifold params that --multi normally shares across indices
# (SCALAR_KEYS plus these). Used by --lpap-only-shared to keep them
# per-index instead, so each index's forcing manifold is free to diverge.
MANIFOLD_KEYS = SCALAR_KEYS + ["offs", "bg", "ann1", "ann2", "year"]

# lpap amplitudes and phases (42 constituents × 2 = 84 values, but we
# pick the top-N by amplitude to keep the search tractable)
LPAP_TOP_N = 10  # evolve top-10 constituents' amp+phase

# IR + ltep: the NAO-relevant combo (IR feedback + modulation freqs)
IR_KEY = ["IR"]

DEFAULT_MAXH = 100  # ceiling for auto-generated harmonic multipliers


def sync_harm_with_resp(idx_p: dict, resp: dict, idx_dir: Path, idx: str) -> dict:
    """resp's NH is authoritative for the HARMONIC COUNT — its own tokens
    are placeholders (e.g. NH="1 1" just means "2 slots"; using literal 1s
    as multipliers would duplicate an existing mode and produce a singular
    regression design matrix). lt.exe.p's "harm" array holds the actual
    per-slot multiplier VALUES that get used in Config's NH override.

    If resp wants more slots than lt.exe.p currently has, extend "harm"
    with new values — each the smallest positive integer <= MAXH not
    already present in "harm" — and persist the addition back to
    lt.exe.p on disk. If resp wants fewer, truncate "harm" (also
    persisted). No-op, no write, if the counts already match."""
    desired_n = len(str(resp.get("NH", "")).split())
    existing = [float(h) for h in idx_p.get("harm", [])]

    if desired_n == len(existing):
        return idx_p

    if desired_n < len(existing):
        new_harm = existing[:desired_n]
    else:
        maxh = int(float(resp.get("MAXH", DEFAULT_MAXH)))
        used = {int(round(h)) for h in existing}
        new_harm = list(existing)
        candidate = 2  # 1 == the base mode itself; skip to avoid a duplicate
        while len(new_harm) < desired_n:
            if candidate > maxh:
                raise ValueError(
                    f"{idx}: resp NH wants {desired_n} harmonics but ran out "
                    f"of unique values <= MAXH={maxh} (have {new_harm})")
            if candidate not in used:
                new_harm.append(float(candidate))
                used.add(candidate)
            candidate += 1

    idx_p = dict(idx_p)
    idx_p["harm"] = new_harm
    (idx_dir / "lt.exe.p").write_text(json.dumps(idx_p, indent=2))
    print(f"  {idx}: resp NH wants {desired_n} harmonic(s) (lt.exe.p had "
          f"{len(existing)}) -> harm={new_harm} (saved to lt.exe.p)")
    return idx_p


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
# Multi-index joint optimization — shared forcing, per-index modulation
# ---------------------------------------------------------------------------

def _shared_signature(shared_params: dict) -> tuple:
    """Hashable snapshot of every shared_params field forward() reads."""
    scal = tuple(float(shared_params.get(k, 0.0)) for k in SCALAR_KEYS)
    lpap = tuple(np.asarray(shared_params["lpap"], dtype=float).ravel().tolist())
    extra = tuple(float(shared_params.get(k, 0.0))
                  for k in ("offs", "bg", "ann1", "ann2", "year"))
    return scal + lpap + extra


def _pip_signature(pip: dict) -> tuple:
    """Hashable snapshot of the per-index fields forward() reads (ltep/IR/harm,
    plus any MANIFOLD_KEYS held per-index under --lpap-only-shared). harm and
    the manifold keys never change during optimization, but it's cheap to
    include them and it keeps the cache correct if that ever stops being true."""
    ltep = tuple(float(v) for v in pip["ltep"])
    ir = float(pip.get("IR", 0.0))
    harm = tuple(float(h) for h in pip.get("harm", []))
    extra = tuple(float(pip[k]) for k in MANIFOLD_KEYS if k in pip)
    return ltep + (ir,) + harm + extra


def _cfg_with_harm(cfg: Config, harm: list) -> Config:
    """A Config whose NH override reflects `harm`, without touching the
    baked-in base Config. Cheap (two dict copies) — used so the harmonic
    search in coordinate_descent_multi can retry NH every trial."""
    overrides = dict(cfg.overrides)
    overrides["NH"] = " ".join(str(int(round(float(h)))) for h in harm)
    return Config(cfg.resp, overrides)


class IndexCache:
    """Per-index memo of (signature -> (cc, forcing)).

    Phase 1/2 trials only change ONE index's per-index params, so the other
    indices' forward() result — including the .dat file read — is exactly
    what's already cached, every time. Phase 3 changes shared_params, which
    every index depends on, so it invalidates all of them (a genuine
    recompute, not a caching gap)."""

    def __init__(self):
        self._store: dict[str, tuple[tuple, float, np.ndarray]] = {}
        self.hits = 0
        self.misses = 0

    def get(self, idx: str, sig: tuple, compute):
        cached = self._store.get(idx)
        if cached is not None and cached[0] == sig:
            self.hits += 1
            return cached[1], cached[2]
        self.misses += 1
        cc, forcing = compute()
        self._store[idx] = (sig, cc, forcing)
        return cc, forcing


def evaluate_multi_index(
    indices: list[str],
    shared_params: dict,
    per_index_params: dict[str, dict],
    configs: dict[str, Config],
    cache: IndexCache,
    objective: str = "mean",
) -> tuple[float, dict[str, float], dict[str, np.ndarray]]:
    """Run forward pipeline for multiple indices with shared forcing.

    configs[idx] must already carry an absolute CLIMATE_INDEX and the final
    NH override baked in (done once in main()) — no os.chdir(), no per-call
    Config() rebuild here, unlike the original evaluate_multi_index."""
    ccs = {}
    forcings = {}
    shared_sig = _shared_signature(shared_params)

    for idx in indices:
        pip = per_index_params[idx]
        sig = (shared_sig, _pip_signature(pip))

        def compute(idx=idx, pip=pip):
            try:
                params = dict(shared_params)
                params["ltep"] = pip["ltep"]
                params["IR"] = pip["IR"]
                cfg = configs[idx]
                if "harm" in pip:
                    params["harm"] = pip["harm"]
                    # configs[idx]'s NH override is baked in at setup time;
                    # the harmonic search below mutates pip["harm"] per
                    # trial, so build a per-trial Config to match it.
                    cfg = _cfg_with_harm(cfg, pip["harm"])
                # Under --lpap-only-shared, shared_params holds only "lpap";
                # each index's own MANIFOLD_KEYS values (saved per-index in
                # pip) fill in the rest instead of a synced shared value.
                for k in MANIFOLD_KEYS:
                    if k in pip:
                        params[k] = pip[k]
                result = forward(cfg, params)
                return float(result["cc"]), result["forcing"]
            except Exception as e:
                print(f"  evaluate_multi_index ERROR for {idx}: {e}", file=sys.stderr)
                return -1.0, np.zeros(1)

        ccs[idx], forcings[idx] = cache.get(idx, sig, compute)

    mean_cc = float(np.mean(list(ccs.values())))
    min_cc = float(np.min(list(ccs.values())))
    score = mean_cc if objective == "mean" else min_cc
    return score, ccs, forcings


def forcing_agreement(forcings: dict[str, np.ndarray]) -> float:
    """Compute mean pairwise correlation between forcing manifolds."""
    idxs = list(forcings.keys())
    if len(idxs) < 2:
        return 1.0
    sims = []
    for i in range(len(idxs)):
        for j in range(i + 1, len(idxs)):
            f1, f2 = forcings[idxs[i]], forcings[idxs[j]]
            min_len = min(len(f1), len(f2))
            if min_len < 10:
                sims.append(0.0)
            else:
                r = np.corrcoef(f1[:min_len], f2[:min_len])[0, 1]
                sims.append(r)
    return float(np.mean(sims))


def coordinate_descent_multi(
    indices: list[str],
    shared_params: dict,
    per_index_params: dict[str, dict],
    configs: dict[str, Config],
    max_rounds: int = 200,
    cache: IndexCache | None = None,
    lpap_only_shared: bool = False,
) -> tuple[dict, dict[str, dict], float]:
    """Coordinate descent on shared params + per-index ltep/IR/manifold.

    Two categories only, matching --lpap-only-shared's intent: shared (all
    indices forced identical) vs. per-index (each index optimizes its own
    copy independently, which is what lets them diverge — not frozen).

    Alternates between:
    1. Optimize shared lpap/scalars (affects all indices) — just lpap when
       lpap_only_shared is set, since MANIFOLD_KEYS then live per-index
       (see per_index_vecs below) instead.
    2. Optimize per-index ltep/IR[/MANIFOLD_KEYS under lpap_only_shared]
       (affects only that index)
    """
    # Flatten the parameter vector for coordinate descent
    # Shared: lpap amp/phase (84 values) [+ 9 scalars = 93, unless
    # lpap_only_shared moves scalars into the per-index vector below]
    # Per-index: ltep (varies per index) + IR (1) [+ MANIFOLD_KEYS (14)
    # under lpap_only_shared]

    # Build initial vector
    scalar_dim = 0 if lpap_only_shared else len(SCALAR_KEYS)
    if lpap_only_shared:
        shared_vec = np.array(shared_params["lpap"])[:, 1:3].flatten()
    else:
        shared_vec = np.concatenate([
            np.array([shared_params.get(k, 0.0) for k in SCALAR_KEYS]),
            np.array(shared_params["lpap"])[:, 1:3].flatten(),  # amp, phase for 42 constituents
        ])

    # ltep length varies by index (1 for tsa, up to 6 elsewhere), so the
    # per-index vector layout is [IR, *ltep, *MANIFOLD_KEYS] with MANIFOLD_KEYS
    # only appended under lpap_only_shared; n_ltep is tracked per index so
    # slot j can always be mapped back to the right field/label.
    per_index_vecs = {}
    per_index_n_ltep = {}
    for idx in indices:
        pip = per_index_params[idx]
        n_ltep = len(pip["ltep"])
        per_index_n_ltep[idx] = n_ltep
        parts = [[pip.get("IR", 0.0)], np.array(pip["ltep"], dtype=float)]
        if lpap_only_shared:
            parts.append(np.array([float(pip.get(k, 0.0)) for k in MANIFOLD_KEYS]))
        per_index_vecs[idx] = np.concatenate(parts)

    if cache is None:
        cache = IndexCache()

    # Initial evaluation
    score, ccs, forcings = evaluate_multi_index(
        indices, shared_params, per_index_params, configs, cache
    )
    print(f"  Initial score: {score:.6f}")
    for idx in indices:
        print(f"    {idx}: CC={ccs[idx]:.6f}")
    print(f"  Forcing agreement: {forcing_agreement(forcings):.4f}")
    
    # Coordinate descent: cycle through shared + per-index params
    round_num = 0
    improved = True
    best_score = score
    
    # Shared step sizes (small for lpap, larger for scalars)
    lpap_steps = np.full(84, 0.001)  # tiny steps for lpap
    if lpap_only_shared:
        shared_steps = lpap_steps
    else:
        scalar_steps = np.maximum(np.abs(shared_vec[:scalar_dim]) * 0.01, 0.001)
        shared_steps = np.concatenate([scalar_steps, lpap_steps])
    
    per_index_steps = {}
    for idx in indices:
        pip = per_index_params[idx]
        ir_step = max(abs(pip.get("IR", 0.0)) * 0.1, 0.05)
        # ltep steps: 50% of value, min 0.30 (wide-ranging exploration)
        ltep_vals = np.array(pip["ltep"], dtype=float)
        ltep_steps = np.maximum(np.abs(ltep_vals) * 0.5, 0.30)
        # Extra aggression on ltep[0]: at least 50% or 0.30
        if len(ltep_steps) > 0:
            ltep_steps[0] = max(ltep_steps[0], 0.30, abs(ltep_vals[0]) * 0.5)
        parts = [[ir_step], ltep_steps]
        if lpap_only_shared:
            # Same step formula as the old shared scalar_steps (1% of
            # magnitude, min 0.001) — now hill-climbed per index instead of
            # jointly, which is what lets indices' manifolds diverge.
            manifold_vals = np.array([float(pip.get(k, 0.0)) for k in MANIFOLD_KEYS])
            manifold_steps = np.maximum(np.abs(manifold_vals) * 0.01, 0.001)
            parts.append(manifold_steps)
        per_index_steps[idx] = np.concatenate(parts)

    def _vec_to_pip(base_pip: dict, trial_vec: np.ndarray, n_ltep: int) -> dict:
        """Rebuild a per-index params dict from a per_index_vecs-layout
        vector: [IR, *ltep, *MANIFOLD_KEYS (only under lpap_only_shared)]."""
        trial_pip = dict(base_pip)
        trial_pip["IR"] = float(trial_vec[0])
        trial_pip["ltep"] = trial_vec[1:1 + n_ltep].tolist()
        if lpap_only_shared:
            manifold_vec = trial_vec[1 + n_ltep:]
            for mi, k in enumerate(MANIFOLD_KEYS):
                trial_pip[k] = float(manifold_vec[mi])
        if "harm" in base_pip:
            trial_pip["harm"] = base_pip["harm"]
        return trial_pip

    print(f"  {'Round':>6} {'Target':<20} {'Param':>6} {'Step':>12} {'Min CC':>10} {'Forcing Sim':>12}")

    while improved and round_num < max_rounds:
        improved = False

        # Track which indices have converged on ltep[0] THIS PASS only.
        # Reset every outer pass: Phase 2 (IR) and Phase 3 (shared lpap/
        # scalars, which shifts every index's forcing manifold) can change
        # the fitness landscape for ltep[0] again, so a lock from an earlier
        # pass must not persist — that was silently freezing indices (e.g.
        # amo never moved its ltep[0] at all; pdo took exactly one step and
        # then never got reconsidered).
        ltep0_converged = {idx: False for idx in indices}

        # Phase 0: Harmonic (NH) search. harm is per-index only (each index's
        # own resp/MAXH), so a trial here only recomputes that one index —
        # IndexCache still serves every other index from cache — making a
        # full sweep of the allowed range cheap enough to run every pass.
        # Values must stay integers in [2, MAXH] (1 is the base mode itself)
        # and mutually unique within one index's harm list.
        for idx in indices:
            if round_num >= max_rounds:
                break
            pip = per_index_params[idx]
            harm = list(pip.get("harm", []))
            if not harm:
                continue
            maxh = int(float(configs[idx].get("MAXH", DEFAULT_MAXH)))

            for slot in range(len(harm)):
                if round_num >= max_rounds:
                    break
                round_num += 1

                current = int(round(harm[slot]))
                used = {int(round(v)) for i, v in enumerate(harm) if i != slot}
                best_val = current
                best_trial_cc = best_score

                for cand in range(2, maxh + 1):
                    if cand == current or cand in used:
                        continue
                    trial_harm = list(harm)
                    trial_harm[slot] = float(cand)
                    trial_pip = dict(pip)
                    trial_pip["harm"] = trial_harm

                    trial_cc, _, _ = evaluate_multi_index(
                        indices, shared_params,
                        {k: (trial_pip if k == idx else per_index_params[k]) for k in indices},
                        configs, cache,
                    )

                    if trial_cc > best_trial_cc:
                        best_trial_cc = trial_cc
                        best_val = cand

                if best_val != current:
                    harm[slot] = float(best_val)
                    pip = dict(pip)
                    pip["harm"] = harm
                    per_index_params[idx] = pip
                    best_score = best_trial_cc
                    improved = True
                    _, ccs, forcings = evaluate_multi_index(
                        indices, shared_params, per_index_params, configs, cache
                    )
                    label = f"{idx} harm[{slot}]"
                    move = f"{current}->{best_val}"
                    print(f"  {round_num:>6} {label:<20} {slot:>6} {move:>12} "
                          f"{best_score:>10.6f} {forcing_agreement(forcings):>12.4f}")

        # Phase 1: Aggressively tune each index's ltep[0] (the big lever)
        # Try multiple step sizes to escape local minima, keep trying until converged
        for idx in indices:
            if round_num >= max_rounds:
                break
            if ltep0_converged.get(idx, False):
                continue
            vec = per_index_vecs[idx]
            steps = per_index_steps[idx]
            pip = per_index_params[idx]

            j = 1  # ltep[0] position
            if j < len(vec):
                round_num += 1
                base_step = steps[j]
                # Try 5 step sizes: 200%, 100%, 50%, 25%, 10% of the base step
                found = False
                for mult in [2.0, 1.0, 0.5, 0.25, 0.1]:
                    if round_num >= max_rounds:
                        break
                    step = base_step * mult
                    for sign in [1, -1]:
                        trial_vec = vec.copy()
                        trial_vec[j] += sign * step

                        trial_pip = _vec_to_pip(pip, trial_vec, per_index_n_ltep[idx])

                        trial_cc, _, _ = evaluate_multi_index(
                            indices, shared_params,
                            {k: (trial_pip if k == idx else per_index_params[k]) for k in indices},
                            configs, cache,
                        )

                        if trial_cc > best_score:
                            per_index_vecs[idx] = trial_vec
                            per_index_params[idx] = trial_pip
                            best_score = trial_cc
                            improved = True
                            found = True
                            _, ccs, forcings = evaluate_multi_index(
                                indices, shared_params, per_index_params, configs, cache
                            )
                            print(f"  {round_num:>6} {idx} ltep[0] x{mult:.2f}     {j:>6} {sign*step:>12.6f} "
                                  f"{best_score:>10.6f} {forcing_agreement(forcings):>12.4f}")
                            break  # Found improvement, move to next index
                    if found:
                        break
                if not found:
                    ltep0_converged[idx] = True  # No step size worked, mark converged

        # Phase 2: Tune remaining per-index params (IR, ltep[1:], and under
        # lpap_only_shared each index's own MANIFOLD_KEYS — this is what
        # actually lets manifolds diverge; they are per-index state that
        # gets optimized independently, not frozen).
        for idx in indices:
            if round_num >= max_rounds:
                break
            vec = per_index_vecs[idx]
            steps = per_index_steps[idx]
            pip = per_index_params[idx]
            n_ltep = per_index_n_ltep[idx]

            for j in range(len(vec)):
                if j == 1:  # Skip ltep[0], already tried in Phase 1
                    continue
                # Skip ltep[1] — shared Bessel frequency, frozen. Only
                # applies when ltep actually has a second element (n_ltep
                # varies 1-6 across indices, e.g. tsa has n_ltep=1); guard
                # it so a short ltep doesn't accidentally freeze what's
                # really the first MANIFOLD_KEYS slot at position 2.
                if j == 2 and n_ltep >= 2:
                    continue
                if round_num >= max_rounds:
                    break
                round_num += 1

                step = steps[j]
                for sign in [1, -1]:
                    trial_vec = vec.copy()
                    trial_vec[j] += sign * step

                    trial_pip = _vec_to_pip(pip, trial_vec, n_ltep)

                    trial_cc, _, _ = evaluate_multi_index(
                        indices, shared_params,
                        {k: (trial_pip if k == idx else per_index_params[k]) for k in indices},
                        configs, cache,
                    )

                    if trial_cc > best_score:
                        per_index_vecs[idx] = trial_vec
                        per_index_params[idx] = trial_pip
                        best_score = trial_cc
                        improved = True
                        _, ccs, forcings = evaluate_multi_index(
                            indices, shared_params, per_index_params, configs, cache
                        )
                        if j == 0:
                            label = "IR"
                        elif j <= n_ltep:
                            label = f"ltep[{j-1}]"
                        else:
                            label = MANIFOLD_KEYS[j - n_ltep - 1]
                        print(f"  {round_num:>6} {idx} {label:<19} {j:>6} {sign*step:>12.6f} "
                              f"{best_score:>10.6f} {forcing_agreement(forcings):>12.4f}")
                        break

        # Phase 3: Fine-tune shared params
        for j in range(len(shared_vec)):
            if round_num >= max_rounds:
                break
            round_num += 1

            step = shared_steps[j]
            for sign in [1, -1]:
                trial_vec = shared_vec.copy()
                trial_vec[j] += sign * step

                trial_shared = dict(shared_params)
                if not lpap_only_shared:
                    trial_shared.update(zip(SCALAR_KEYS, trial_vec[:scalar_dim]))
                trial_lpap = np.array(shared_params["lpap"])
                trial_lpap[:, 1:3] = trial_vec[scalar_dim:].reshape(42, 2)
                trial_shared["lpap"] = trial_lpap.tolist()

                trial_cc, _, _ = evaluate_multi_index(
                    indices, trial_shared, per_index_params, configs, cache
                )

                if trial_cc > best_score:
                    shared_vec = trial_vec
                    shared_params = trial_shared
                    best_score = trial_cc
                    improved = True
                    _, ccs, forcings = evaluate_multi_index(
                        indices, shared_params, per_index_params, configs, cache
                    )
                    print(f"  {round_num:>6} {'shared':<20} {j:>6} {sign*step:>12.6f} "
                          f"{best_score:>10.6f} {forcing_agreement(forcings):>12.4f}")
                    break

    print(f"\n  Final score: {best_score:.6f}")
    _, final_ccs, final_forcings = evaluate_multi_index(
        indices, shared_params, per_index_params, configs, cache
    )
    for idx in indices:
        print(f"    {idx}: CC={final_ccs[idx]:.6f}")
    print(f"  Final forcing agreement: {forcing_agreement(final_forcings):.4f}")
    print(f"  Cache: {cache.hits} hits / {cache.hits + cache.misses} lookups "
          f"({100.0 * cache.hits / max(cache.hits + cache.misses, 1):.1f}% hit rate)")

    return shared_params, per_index_params, best_score

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

    print(f"{'Gen':>4} {'Best CC':>10} {'Min CC':>10} {'Delta':>10} {'MutRate':>8}")

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
    ap.add_argument("--multi", nargs="+", default=None,
                    help="multi-index joint optimization: list of indices")
    ap.add_argument("--multi-rounds", type=int, default=200,
                    help="max rounds for multi-index optimization (default 200)")
    ap.add_argument("--lpap-only-shared", action="store_true",
                    help="--multi: only share lpap (amp/phase of tidal "
                         "constituents) across indices; scalars/offs/bg/"
                         "ann1/ann2/year become per-index instead of shared "
                         "and are optimized independently per index, so "
                         "forcing manifolds are free to diverge")
    ap.add_argument("--reset", action="store_true",
                    help="clear multi_index/ and start from lt.exe.p")
    ap.add_argument("--from-original", action="store_true",
                    help="start from original lt.exe.p, skipping .coord.p files")
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
    ap.add_argument("--metric", default="cc", choices=["cc", "dtw"],
                    help="fit metric to maximize: Pearson's cc (default) or "
                         "dtw_distance (a normalized, banded DTW score — "
                         "less sensitive to time-series alignment than cc)")
    ap.add_argument("--dtw-window", type=int, default=3,
                    help="Sakoe-Chiba band half-width for --metric dtw "
                         "(default 3)")
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
    overrides = {"TRAIN_START": args.cv[0], "TRAIN_END": args.cv[1],
                 "METRIC": args.metric.upper(), "DTW_WINDOW": args.dtw_window}

    # resp's NH is authoritative for harmonic COUNT; reconcile lt.exe.p's
    # "harm" values to match (adding unique values up to MAXH if resp wants
    # more than lt.exe.p currently has) before deriving the NH override from
    # those (real) values. See sync_harm_with_resp.
    params = sync_harm_with_resp(params, resp, run_dir, args.index)
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

    # Multi-index joint optimization
    if args.multi:
        indices = args.multi
        print(f"\nMulti-index optimization: {indices}")

        multi_dir = ROOT / "multi_index"

        # Reset: clear previous multi_index output
        if args.reset and multi_dir.exists():
            import shutil
            shutil.rmtree(multi_dir)
            print(f"Cleared {multi_dir}")

        multi_dir.mkdir(exist_ok=True)

        # Load shared params: try previous multi_index output first
        shared_params = None
        loaded_from = None

        # Build candidate list based on --from-original flag
        if args.from_original or args.reset:
            # Only load from original lt.exe.p
            candidates = [
                (ROOT / indices[0] / "lt.exe.p", f"{indices[0]}/lt.exe.p"),
            ]
        else:
            # Full cascade
            candidates = [
                (multi_dir / "shared_params.p", "multi_index/shared_params.p"),
                (ROOT / indices[0] / "lt.exe.pysr.lpap_top.coord.p", f"{indices[0]}/lt.exe.pysr.lpap_top.coord.p"),
                (ROOT / indices[0] / "lt.exe.pysr.ir_ltep.coord.p", f"{indices[0]}/lt.exe.pysr.ir_ltep.coord.p"),
                (ROOT / indices[0] / "lt.exe.p", f"{indices[0]}/lt.exe.p"),
            ]
        for cand_path, cand_name in candidates:
            if cand_path.exists():
                shared_params = json.loads(cand_path.read_text())
                loaded_from = cand_name
                break

        if shared_params is None:
            print(f"ERROR: No params file found for any index")
            return 2
        print(f"Loaded shared params from: {loaded_from}")

        # Full dict as loaded, before filtering down to what stays shared —
        # kept so --lpap-only-shared can seed each index's own MANIFOLD_KEYS
        # from the value that WAS shared (i.e. the previous run's result),
        # not from each index's possibly-stale lt.exe.p on disk. That's what
        # guarantees the initial score under --lpap-only-shared is at least
        # the previous run's final score: every index starts with the exact
        # params that produced it, so the first move can only do as well or
        # better.
        loaded_shared_params = dict(shared_params)

        # Keep only shared params: lpap [+ scalars, unless --lpap-only-shared
        # keeps those per-index so the manifolds don't converge exactly]
        if args.lpap_only_shared:
            shared_params = {k: v for k, v in shared_params.items() if k == "lpap"}
            print("--lpap-only-shared: only lpap is shared across indices; "
                  f"{MANIFOLD_KEYS} stay per-index (seeded from the "
                  "previously-shared values)")
        else:
            shared_keys = set(SCALAR_KEYS) | {"lpap", "offs", "bg", "ann1", "ann2", "year"}
            shared_params = {k: v for k, v in shared_params.items() if k in shared_keys}

        # Per-index params: ltep + IR + harm from each index's .p file
        # Also try previous multi_index output
        per_index_params = {}
        configs = {}

        # Load previous per-index params if available (and not --from-original/--reset)
        prev_per_index = None
        if not args.from_original and not args.reset:
            if (multi_dir / "per_index_params.json").exists():
                prev_per_index = json.loads((multi_dir / "per_index_params.json").read_text())
                print(f"Loaded previous per-index params from: multi_index/per_index_params.json")

        for idx in indices:
            idx_dir = ROOT / idx
            idx_p = json.loads((idx_dir / "lt.exe.p").read_text())
            resp = read_resp(idx_dir / "lt.exe.resp")
            # Absolute path: forward() no longer needs cwd == idx_dir, so
            # evaluate_multi_index can skip the os.chdir() dance entirely.
            resp["CLIMATE_INDEX"] = str(idx_dir / f"{idx}.dat")

            # resp's NH is authoritative for harmonic COUNT; reconcile
            # lt.exe.p's "harm" values (add unique ones up to MAXH, or
            # truncate) to match before deriving the NH override from them.
            idx_p = sync_harm_with_resp(idx_p, resp, idx_dir, idx)

            # NH from harm
            nh_override = {}
            if "harm" in idx_p:
                nh_override["NH"] = " ".join(str(int(float(h))) for h in idx_p["harm"])

            cfg = Config(resp, {
                "TRAIN_START": "1880", "TRAIN_END": "1885",
                "METRIC": args.metric.upper(), "DTW_WINDOW": args.dtw_window,
                **nh_override,
            })
            configs[idx] = cfg

            # Use previous multi_index output if available, else original .p
            if prev_per_index and idx in prev_per_index:
                per_index_params[idx] = dict(prev_per_index[idx])
            else:
                per_index_params[idx] = {
                    "ltep": idx_p["ltep"],
                    "IR": idx_p.get("IR", 0.0),
                }

            # "harm" is now actively optimized per index (coordinate_descent_
            # multi's harmonic-search phase), so a continued run keeps
            # whatever it last found instead of resetting to idx_p's
            # original values every time. Only fall back to idx_p's "harm"
            # on a fresh start (no prior value) or when its length no longer
            # matches — resp's NH count changed and sync_harm_with_resp
            # above already reconciled idx_p's "harm" to the new count.
            prev_harm = per_index_params[idx].get("harm")
            if "harm" in idx_p:
                if prev_harm is None or len(prev_harm) != len(idx_p["harm"]):
                    per_index_params[idx]["harm"] = idx_p["harm"]
            else:
                per_index_params[idx].pop("harm", None)

            # --lpap-only-shared: each index keeps its own scalars/offs/bg/
            # ann1/ann2/year instead of a value synced from the shared
            # cascade going forward. Seed them from loaded_shared_params
            # (the value that WAS shared, i.e. the previous run's result) —
            # not idx_p, which can be stale — so every index starts from
            # exactly the params that produced the previous run's score.
            # Only fill in keys not already present (e.g. from a continued
            # prior lpap-only-shared run) so a resumed run doesn't silently
            # reset its own per-index progress.
            if args.lpap_only_shared:
                for k in MANIFOLD_KEYS:
                    if k in per_index_params[idx]:
                        continue
                    if k in loaded_shared_params:
                        per_index_params[idx][k] = loaded_shared_params[k]
                    elif k in idx_p:
                        per_index_params[idx][k] = idx_p[k]
        
        # Single cache shared across the initial baseline and the whole
        # coordinate descent run, so every reused (shared, per-index)
        # combination anywhere is only ever computed once.
        cache = IndexCache()

        # Initial baseline
        score, ccs, forcings = evaluate_multi_index(
            indices, shared_params, per_index_params, configs, cache
        )
        print(f"Initial score: {score:.6f}")
        for idx in indices:
            print(f"  {idx}: CC={ccs[idx]:.6f}")
        print(f"Initial forcing agreement: {forcing_agreement(forcings):.4f}")
        
        # Run optimization
        shared_out, per_index_out, final_cc = coordinate_descent_multi(
            indices, shared_params, per_index_params, configs,
            max_rounds=args.multi_rounds, cache=cache,
            lpap_only_shared=args.lpap_only_shared,
        )
        
        # Save results
        out_dir = ROOT / "multi_index"
        out_dir.mkdir(exist_ok=True)
        with open(out_dir / "shared_params.p", "w") as f:
            json.dump(shared_out, f, indent=2)
        with open(out_dir / "per_index_params.json", "w") as f:
            json.dump(per_index_out, f, indent=2)
        print(f"\nSaved shared params to {out_dir / 'shared_params.p'}")
        print(f"Saved per-index params to {out_dir / 'per_index_params.json'}")
        
        return 0

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
