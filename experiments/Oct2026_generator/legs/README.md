# Deterministic (Non-Stochastic) Multi-Leg Fitting with `lpap_refine.py`

`lpap_refine.py` provides a **100% deterministic, zero-Monte-Carlo** gradient + discrete-screening optimizer for single-leg and multi-leg (`LEGS N`, `LEG_WINDINGS TRUE`, `LEG_P PARAMS`) runs. Given a run folder's `lt.exe.p` (and `lt.exe.p.<leg>`), `lt.exe.resp*`, and `lte_run.sh`, it replaces the long stochastic random-walk tail of `lt.exe` (`Walker.Markov` / `FLIP`) with trust-region Gauss-Newton and greedy harmonic screening, verified at every round by the real `lt.exe` solver (`TEST_ONLY=TRUE`).

## How the Non-Stochastic Solver Works

1. **Exact Python surrogate verification**: Rebuilds each leg's post-Bessel response (`impA`, `impB`, `besk`, `offs`, `bg`, `ltep`, `harm`, `IR`, `F9`, `TREND`, `ACCEL`) and aborts if any leg differs from `lt.exe` by more than `1e-6`.
2. **Analytic pre-Bessel Jacobian (~100 ms vs ~15–35 min of FD probes)**: Computes the exact linear tangent of `Extend_Backward → Exact_Template → Tide_Sum_Diff → Impulse_Amplify → Shape_Impulses → IIR` with respect to all `2 × 52 = 104` tidal `(A cos φ, A sin φ)` coefficients for every leg (supporting `STRICT_IDATE`, `LEG_P PARAMS`, `delC`, `mB`, `shfT`, `JERK`, `LOD_TERM`, `IMPULSE_MOD12`, `IMPULSE_EXACT`, `ALIAS`, and non-monthly `SAMPLING`). A single-column `TEST_ONLY` spot-check (`h = 1e-4`) verifies agreement to `< 1e-4` (typically `~1e-9`) concurrently with `dlod_ref`, falling back to finite differences (`--fd-jac`) only if an uncovered mode is active.
3. **Coarse-to-fine frequency continuation (`--homotopy`)**: Solves cycle-skipping on `sin(2π k F)` from cold starts (e.g., `DLOD_REF=FALSE`) by optimizing `lpap` progressively across winding cutoffs `|k| <= 0.35 → 0.90 → 2.50` before turning on the full harmonic spectrum.
4. **Block-alternating Gauss-Newton (`--per-leg`)**: Alternates:
   - **Block A**: Shared `lpap` trust-region `least_squares` with a daily-dLOD L2 penalty (`--lam`).
   - **Block B**: Per-leg continuous trust-region `least_squares` over each leg's post-manifold parameters (`impA`, `impB`, `besk`, `offs`, `bg`, `ltep`, `IR`), with `|ltep| >= MIN_LT` bounds and mild L2 drift regularization (`--leg-lam`). Because Block B is post-manifold, the pre-Bessel Jacobian stays valid across all rounds.
5. **Deterministic harmonic screening & expansion (`--harm-screen`, `--harm-expand N`)**: Replaces stochastic `FLIP` mutations with greedy 1-swap coordinate ascent over `2..MAXH` (collision-checked against `ltep`, gated by `--harm-min-gain` and `--harm-max-swaps`). `--harm-expand N` greedily appends `N` new harmonic slots and updates `NH` in `lt.exe.resp` on `--apply`.

## Usage Instructions

Run from inside a leg run folder (e.g., `experiments/Oct2026_generator/legs/atlantic_family`):

### A. Refining an existing fit (warm basin)
```bash
# 1. Refine shared lpap + per-leg continuous params + discrete harmonics (~2-3 min)
python3 ../lpap_refine.py . --full --lam 1.0 --rounds 4 --apply

# 2. Score and regenerate lte_results.<leg>.csv, windings JSON, and legs.png (~20 s)
TEST_ONLY=TRUE ./lte_run.sh
```
- Backups are saved automatically as `lt.exe.p.pre_refine`, `lt.exe.p.<leg>.pre_refine` (and `lt.exe.resp.pre_refine` if `--harm-expand` changes `NH`).
- To refine shared `lpap` only (keeping per-leg files fixed), omit `--full`.

### B. Cold-starting `lpap` from measured dLOD (`DLOD_REF=FALSE`, zero stochastic search)
```bash
# 1. Seed lt.exe.p from measured dLOD (dlod3.dat) without running Monte Carlo search
DLOD_REF=FALSE TEST_ONLY=TRUE ./lte_run.sh

# 2. Coarse-to-fine frequency continuation + alternating refinement (~4-5 min)
#    Use --lam 0.5 --rounds 2 (or --accept validate with VALIDATE_START) to guard held-out generalization
python3 ../lpap_refine.py . --homotopy --full --lam 0.5 --rounds 2 --apply

# 3. Regenerate final CSVs and legs.png
TEST_ONLY=TRUE ./lte_run.sh
```

### Key CLI Flags
| Flag | Default | Purpose |
|---|---|---|
| `--lam` | `1.0` | Daily-dLOD L2 regularization weight (higher = closer to measured dLOD, stronger generalization) |
| `--rounds` / `--iters` | `8` / `48` | Outer solver-verified trust-region rounds and total surrogate evaluations |
| `--homotopy` | off | 3-stage coarse-to-fine frequency continuation (`\|k\| <= 0.35 → 0.90 → 2.50`) before full rounds |
| `--per-leg` | off | Enable Block B per-leg continuous GN (`--leg-step 0.05`, `--leg-lam 1e-3`) |
| `--harm-screen` | off | Greedy 1-swap harmonic screening over `2..MAXH` (`--harm-max-swaps 2`, `--harm-min-gain 1e-4`) |
| `--harm-expand N` | `0` | Append `N` greedy harmonic slots to `NH` (updates `NH` in `lt.exe.resp` on `--apply`) |
| `--full` | off | Shorthand for `--per-leg --harm-screen` |
| `--accept validate` | `train` | Accept a round only if validation CC (`t >= VALIDATE_START`) does not drop by more than `--vtol` |
| `--fd-jac` | off | Force legacy finite-difference solver probes instead of the analytic pre-Bessel Jacobian |
| `--apply` | off | Write accepted parameters into `lt.exe.p` (and `lt.exe.p.<leg>`) with `*.pre_refine` backups |

## Benchmarks (2026-11-24, held out `1960–1970`)

### 1. Warm `atlantic_family` (`NAO, Baltic, AMO, TNA`, `STRICT_IDATE TRUE`, `IDATE 1880.0`, `NH=7`, `MAXH=30`)
| Run | Wall time | NAO (tr/**ho**) | Baltic (tr/**ho**) | AMO (tr/**ho**) | TNA (tr/**ho**) | Mean (tr / **ho**) | dLOD |
|---|---|---|---|---|---|---|---|
| Baseline (`STRICT_IDATE` + `baltic_test`) | — | `0.827` / `+0.443` | `0.618` / `+0.337` | `0.769` / `+0.473` | `0.819` / `+0.396` | `0.7588` / `+0.412` | `0.9840` |
| `lpap_refine.py --lam 1.0 --rounds 4` (analytic J) | **182 s** | `0.827` / `+0.447` | `0.622` / `+0.358` | `0.772` / **`+0.493`** | `0.821` / **`+0.445`** | `0.7607` / **`+0.436`** | `0.9819` |
| `lpap_refine.py --full --lam 1.0 --rounds 4` | **158 s** | `0.828` / **`+0.461`** | `0.622` / **`+0.366`** | `0.773` / `+0.474` | `0.822` / `+0.437` | **`0.7616`** / `+0.435` | `0.9818` |

*(Analytic Jacobian cut the 4-leg Jacobian build from `918 s` to `23 s` with `3.9e-9` spot-check rel diff; total runtime dropped `5.9×` from `1075 s` to `182 s` with identical results.)*

### 2. Pure Cold-Start `lpap` from Measured dLOD (`DLOD_REF=FALSE`, zero stochastic search)
Starting from `DLOD_REF=FALSE` (`dLOD = 1.0000`) with `--homotopy --full --lam 0.3 --rounds 4` (non-`lpap` scalars initialized from the run's `.p` files):

- **4-leg `atlantic_family`** (`266 s`, spot-check `9.47e-10`):
  - `baseline (DLOD_REF=FALSE)`: mean `0.4742` train / `+0.1155` held-out (`Baltic 0.323/+0.046`)
  - `after homotopy` (`60 s`): mean `0.7277` train / **`+0.1838`** held-out (`NAO 0.757/+0.473`, `Baltic 0.625/+0.260`, `AMO 0.782/+0.310`)
  - `after round 1`: mean `0.7506` train / `+0.1719` held-out (**`NAO 0.774/+0.520`**, `AMO 0.794/+0.334`)
  - `after round 4`: mean **`0.7661`** train / `+0.1064` held-out (`NAO 0.775`, **`Baltic 0.693`**, **`AMO 0.809`**, `TNA 0.787`, `dLOD 0.9709`)
- **6-leg `6indices`** (`NAO, NINO4, AMO, PDO, IODE, TNA`, `334 s`, spot-check `2.83e-09`):
  - `baseline (DLOD_REF=FALSE)`: mean `0.4934` train / `+0.0775` held-out
  - `after homotopy` (`64 s`): mean `0.7114` train / `+0.1086` held-out (5 of 6 legs positive held-out)
  - `after round 1`: mean `0.7364` train / **`+0.1303`** held-out (`NINO4 0.741/+0.349`, `AMO 0.793/+0.275`, `IODE 0.714/+0.449`)
  - `after round 4`: mean **`0.7528`** train / `+0.0686` held-out (`NAO 0.650`, `NINO4 0.759`, `AMO 0.795`, `PDO 0.767`, `IODE 0.741`, `TNA 0.804`, `dLOD 0.9781`)

**Generalization note**: On cold starts, held-out CC peaks right after `homotopy` / `round 1`, while additional rounds at low `--lam 0.3` push training CC higher (e.g., Baltic to `0.693`) at the expense of held-out overfitting. For best generalization, use `--lam 0.5–1.0` and stop after `1–2` rounds (or use `--accept validate` with `VALIDATE_START`).

### 3. Hybrid Basin-Hopping (`basin_hop.py`) on 7-Leg `7indices` (`NAO, Baltic, AMO, TNA, NINO4, PDO, IODE`)
For longer / overnight runs, `basin_hop.py` alternates **Phase A** (`LOCKT=TRUE` non-LPAP exploration, with optional `--perturb` kicks on stalls) and **Phase B** (`lpap_refine.py --full` deterministic polish), checkpointing best parameters after every accepted cycle:
```bash
python3 -u ../basin_hop.py . --cycles 12 --nonlpap-sec 180 --threads 14 --lam 1.0 --rounds 3 --iters 24 --perturb 0.03 --min-dlod 0.98
```

| Stage (`7indices`, `LEGS 7`, `STRICT_IDATE TRUE`, `IDATE 1880.0`) | Mean train | dLOD | NAO | Baltic | AMO | TNA | NINO4 | PDO | IODE |
|---|---|---|---|---|---|---|---|---|---|
| Seed (`atlantic_family` + `6indices` per-leg) | `0.7040` | `0.9818` | `0.828` | `0.622` | `0.773` | `0.822` | `0.560` | `0.596` | `0.725` |
| Initial `lpap_refine.py --full --lam 1.0` | `0.7401` | `0.9828` | `0.820` | `0.575` | `0.769` | `0.815` | `0.688` | `0.752` | `0.760` |
| `basin_hop.py` (3 + 12 cycles, `LOCKT=TRUE` + `--full`) | `0.7518` | `0.9824` | `0.816` | `0.591` | `0.767` | `0.820` | `0.714` | `0.779` | `0.775` |
| `--lam 0.5` polish + 20 min joint `LOCKT=FALSE` + `--full --lam 1.0` | **`0.7531`** | **`0.9820`** | `0.816` | `0.592` | **`0.773`** | **`0.820`** | `0.713` | **`0.781`** | **`0.776`** |

*(Pre-Bessel manifold agreement remains `1.00000` across all 7 legs; 2000–05 CCs include `NINO4 +0.704`, `TNA +0.609`, `PDO +0.579`, `Baltic +0.541`.)*

