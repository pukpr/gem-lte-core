#!/usr/bin/env bash
# Standalone single-leg & multi-leg GEM-LTE launcher (Meeus-free TIDES=POLY default).
# Runs in the caller's working directory (PWD), so BOTH `./lte_run.sh` (from
# experiments/Dec2026) and `../lte_run.sh --legs ...` (from any legs subdir) work.
#
# Usage:
#   ./lte_run.sh                                          # single-leg pdo.dat (POLY Z10, cold + homotopy)
#   ../lte_run.sh --legs nao.dat baltic_3pt_na.dat \
#       --start 1970 --end 1980 --refine --lam 1.0 --rounds 1 --iters 4 --threads 2
#   ./lte_run.sh --tier Z10_S5_N5 --plot
#   ./lte_run.sh --test --plot
#
# Options:
#   --legs f1.dat [f2.dat ...]  Multi-leg climate indices (auto-links from data/, writes lt.exe.resp.1..N)
#   --dat FILE                  Single climate index file (default: pdo.dat when --legs not given)
#   --tides MODE                TIDES mode: POLY (default) or TABLE
#   --tier TIER                 POLY_TIER: Z9, Z10, Z15, Z10_D6, Z10_N5, Z10_D6_N5, Z10_S5_N5, Z10_S5_D6_N5
#   --n5-ridge VAL              L2 ridge on multi-year N5 terms (default: 0.05)
#   --start YEAR                TRAIN_START (default: 1960)
#   --end YEAR                  TRAIN_END   (default: 1970)
#   --exclude true|false        EXCLUDE hold-out window (default: true)
#   --cold                      Cold start from dlod3.dat calibration (DLOD_REF=FALSE)
#   --warm                      Warm start from existing lt.exe.p (DLOD_REF=TRUE)
#   --refine                    Run native Ada REFINE=TRUE, then score with TEST_ONLY=TRUE
#   --homotopy                  Enable REFINE_HOMOTOPY=TRUE (implies --refine)
#   --search                    Run stochastic enso_opt search instead of REFINE
#   --lam VAL                   REFINE_LAM dLOD prior weight (default: 1.0)
#   --hlam VAL                  REFINE_HOMOTOPY_LAM (default: same as --lam)
#   --rounds N                  REFINE_ROUNDS (default: 4)
#   --iters N                   REFINE_ITERS  (default: 60)
#   --test                      Score current lt.exe.p only (TEST_ONLY=TRUE)
#   --plot                      Also run plot_lpap_year.py (--year 2000 --years 3 --measured)
#   --threads N                 NUMBER_OF_PROCESSORS (default: 4)
#   --timeout SEC               TIMEOUT for --search mode (default: 300)
#   -h, --help                  Show this help
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ulimit -s unlimited 2>/dev/null || true

LT_EXE="$SCRIPT_DIR/lt.exe"
DATA_DIR="$SCRIPT_DIR/data"
DLOD_PATH="$DATA_DIR/dlod3.dat"
if [[ ! -e "$DLOD_PATH" && -e "$SCRIPT_DIR/dlod3.dat" ]]; then
  DLOD_PATH="$SCRIPT_DIR/dlod3.dat"
fi
SCORES_PY="$SCRIPT_DIR/leg_scores2.py"
PLOT_PY="$SCRIPT_DIR/plot_lpap_year.py"

if [[ ! -x "$LT_EXE" ]]; then
  echo "ERROR: executable not found: $LT_EXE" >&2
  exit 1
fi

DAT="${CLIMATE_INDEX:-pdo.dat}"
TIDES_MODE="${TIDES:-POLY}"
TIER="${POLY_TIER:-Z10}"
N5_RIDGE="${POLY_N5_RIDGE:-0.05}"
T_START="${TRAIN_START:-1960}"
T_END="${TRAIN_END:-1970}"
EXC="${EXCLUDE:-TRUE}"
DLOD_MODE=""
MODE="homotopy"   # homotopy | refine | search | test
R_LAM="${REFINE_LAM:-1.0}"
R_HLAM="${REFINE_HOMOTOPY_LAM:--1.0}"
R_ROUNDS="${REFINE_ROUNDS:-4}"
R_ITERS="${REFINE_ITERS:-60}"
NPROC="${NUMBER_OF_PROCESSORS:-4}"
TOUT="${TIMEOUT:-300}"
DO_PLOT="false"
LEGS_FILES=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --legs)
      shift
      while [[ $# -gt 0 && "$1" != -* ]]; do
        LEGS_FILES+=("$1")
        shift
      done
      ;;
    --dat)                 DAT="$2"; shift 2 ;;
    --tides)               TIDES_MODE="${2^^}"; shift 2 ;;
    --tier)                TIER="${2^^}"; shift 2 ;;
    --n5-ridge)            N5_RIDGE="$2"; shift 2 ;;
    --start|--train-start) T_START="$2"; shift 2 ;;
    --end|--train-end)     T_END="$2"; shift 2 ;;
    --exclude)             EXC="${2^^}"; shift 2 ;;
    --cold)                DLOD_MODE="FALSE"; shift ;;
    --warm)                DLOD_MODE="TRUE"; shift ;;
    --refine)              MODE="refine"; shift ;;
    --homotopy)            MODE="homotopy"; shift ;;
    --search)              MODE="search"; shift ;;
    --test)                MODE="test"; shift ;;
    --lam)                 R_LAM="$2"; shift 2 ;;
    --hlam|--homotopy-lam) R_HLAM="$2"; shift 2 ;;
    --rounds)              R_ROUNDS="$2"; shift 2 ;;
    --iters)               R_ITERS="$2"; shift 2 ;;
    --threads)             NPROC="$2"; shift 2 ;;
    --timeout)             TOUT="$2"; shift 2 ;;
    --plot)                DO_PLOT="true"; shift ;;
    -h|--help)
      sed -n '2,35p' "$0" | sed 's/^# \{0,1\}//'
      exit 0
      ;;
    *)
      echo "Unknown option: $1 (use --help)" >&2
      exit 1
      ;;
  esac
done

if [[ -z "$DLOD_MODE" ]]; then
  if [[ "$MODE" == "test" ]]; then DLOD_MODE="TRUE"; else DLOD_MODE="FALSE"; fi
fi

link_data_file() {
  local fname="$1"
  if [[ -e "$fname" ]]; then return 0; fi
  local src=""
  for cand in "$DATA_DIR/$fname" "$SCRIPT_DIR/$fname"; do
    if [[ -e "$cand" ]]; then src="$cand"; break; fi
  done
  if [[ -z "$src" ]]; then
    echo "ERROR: could not locate $fname in PWD, $DATA_DIR, or $SCRIPT_DIR" >&2
    exit 1
  fi
  local rel
  rel="$(python3 -c 'import os,sys; print(os.path.relpath(sys.argv[1], os.getcwd()))' "$src")"
  ln -sf "$rel" "$fname"
}

seed_leg_p() {
  local dat_f="$1"
  local dst="$2"
  if [[ -f "$dst" ]]; then return 0; fi
  local base
  base="$(basename "$dat_f" .dat)"
  local atl="$SCRIPT_DIR/../Oct2026_generator/legs/baltic_3pt_overnight/repro_0693/R5_atl4_cold_H10_l03"
  local ip3="$SCRIPT_DIR/../Oct2026_generator/legs/indopacific_family"
  local cand=""
  case "${base,,}" in
    nao*)     cand="$atl/lt.exe.p.1.pre_refine" ;;
    baltic*)  cand="$atl/lt.exe.p.2.pre_refine" ;;
    amo*)     cand="$atl/lt.exe.p.3.pre_refine" ;;
    tna*)     cand="$atl/lt.exe.p.4.pre_refine" ;;
    nino4*)   cand="$ip3/lt.exe.p.1.pre_refine" ;;
    pdo*)     cand="$ip3/lt.exe.p.2.pre_refine" ;;
    iode*)    cand="$ip3/lt.exe.p.3.pre_refine" ;;
  esac
  if [[ -n "$cand" && -f "$cand" ]]; then
    cp "$cand" "$dst"
  elif [[ -f "$SCRIPT_DIR/lt.exe.p" ]]; then
    cp "$SCRIPT_DIR/lt.exe.p" "$dst"
  elif [[ -f "$SCRIPT_DIR/pdo/lt.exe.p" ]]; then
    cp "$SCRIPT_DIR/pdo/lt.exe.p" "$dst"
  elif [[ -f "$atl/lt.exe.p.1.pre_refine" ]]; then
    cp "$atl/lt.exe.p.1.pre_refine" "$dst"
  fi
}

leg_label() {
  local base
  base="$(basename "$1" .dat)"
  case "${base,,}" in
    nao*)     echo "NAO" ;;
    baltic*)  echo "Baltic" ;;
    amo*)     echo "AMO" ;;
    tna*)     echo "TNA" ;;
    nino4*)   echo "NINO4" ;;
    nino34*)  echo "NINO34" ;;
    pdo*)     echo "PDO" ;;
    iode*)    echo "IODE" ;;
    soi*)     echo "SOI" ;;
    *)        echo "${base^^}" ;;
  esac
}

# Bootstrap lt.exe.resp and lt.exe.p in PWD if absent
if [[ ! -f lt.exe.resp ]]; then
  if [[ -f "$SCRIPT_DIR/lt.exe.resp" ]]; then
    cp "$SCRIPT_DIR/lt.exe.resp" lt.exe.resp
  elif [[ -f "$SCRIPT_DIR/pdo/lt.exe.resp" ]]; then
    cp "$SCRIPT_DIR/pdo/lt.exe.resp" lt.exe.resp
  else
    cat > lt.exe.resp <<'RESP'
-FSTART                        0.1
FLIP                           0.0015
SPREAD_MIN                     0.000001
SPREAD_MAX                     0.45
MAXH                           30
NH                             "1 1 1 1 1 1 1"
NM                             3
YEAR                           0.00405305
DECAY                          48.0
DTW                            3
IDATE                          1880.0
INIT_DATE                      1880.0
STRICT_IDATE                   TRUE
THRESHOLD_ACTION               CATCHUP
PART                           FALSE
FILTER                         0.0
JERK                           0.0
RESP
  fi
fi
if [[ ! -f lt.exe.p ]]; then
  for p_cand in \
    "$SCRIPT_DIR/lt.exe.p" \
    "$SCRIPT_DIR/pdo/lt.exe.p" \
    "$SCRIPT_DIR/../Oct2026_generator/legs/baltic_3pt_overnight/repro_0693/R5_atl4_cold_H10_l03/lt.exe.p.pre_refine" \
    "$SCRIPT_DIR/../Oct2026_generator/legs/indopacific_family/lt.exe.p.pre_refine"
  do
    if [[ -f "$p_cand" ]]; then
      cp "$p_cand" lt.exe.p
      break
    fi
  done
  if [[ ! -f lt.exe.p ]]; then
    echo "ERROR: could not find a seed lt.exe.p to bootstrap $PWD" >&2
    exit 1
  fi
fi

N_LEGS=0
if [[ ${#LEGS_FILES[@]} -gt 0 ]]; then
  N_LEGS=${#LEGS_FILES[@]}
  # Remove stale lt.exe.resp.<k> beyond N_LEGS
  for old in lt.exe.resp.[0-9]*; do
    [[ -e "$old" ]] || continue
    idx="${old##*.}"
    if [[ "$idx" =~ ^[0-9]+$ ]] && (( idx > N_LEGS )); then
      rm -f "$old"
    fi
  done
  for i in "${!LEGS_FILES[@]}"; do
    k=$((i + 1))
    dat_f="${LEGS_FILES[$i]}"
    link_data_file "$dat_f"
    seed_leg_p "$dat_f" "lt.exe.p.${k}"
    printf "CLIMATE_INDEX                  %s\n" "$dat_f" > "lt.exe.resp.${k}"
  done
elif [[ -f lt.exe.resp.1 ]]; then
  for k in 1 2 3 4 5 6 7 8 9; do
    if [[ -f "lt.exe.resp.${k}" ]]; then
      N_LEGS=$k
      dat_f="$(awk '$1=="CLIMATE_INDEX"{print $2; exit}' "lt.exe.resp.${k}")"
      [[ -n "$dat_f" ]] && link_data_file "$dat_f"
    else
      break
    fi
  done
else
  link_data_file "$DAT"
fi

# Sync lt.exe.resp cleanly for single-leg vs multi-leg + TIDES/POLY_TIER
python3 - "$N_LEGS" "$DAT" "$TIDES_MODE" "$TIER" <<'PY'
import re, sys
from pathlib import Path
n_legs, dat, tides, tier = int(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4]
p = Path("lt.exe.resp")
drop = r"^\s*(CLIMATE_INDEX|TIDES|POLY_TIER|LEGS|LEG_WINDINGS|LEG_P)\b"
lines = [l for l in p.read_text().splitlines() if not re.match(drop, l, re.I)]
if n_legs > 0:
    lines.append(f"{'LEGS':<30s} {n_legs}")
    lines.append(f"{'LEG_WINDINGS':<30s} TRUE")
    lines.append(f"{'LEG_P':<30s} PARAMS")
else:
    lines.append(f"{'CLIMATE_INDEX':<30s} {dat}")
lines.append(f"{'TIDES':<30s} {tides}")
if tides == "POLY":
    lines.append(f"{'POLY_TIER':<30s} {tier}")
p.write_text("\n".join(lines) + "\n")
PY

COMMON_ENV=(
  "DLOD_DAT=$DLOD_PATH"
  "TIDES=$TIDES_MODE"
  "POLY_TIER=$TIER"
  "POLY_N5_RIDGE=$N5_RIDGE"
  "TRAIN_START=$T_START"
  "TRAIN_END=$T_END"
  "EXCLUDE=$EXC"
  "NUMBER_OF_PROCESSORS=$NPROC"
)
if (( N_LEGS == 0 )); then
  COMMON_ENV+=("CLIMATE_INDEX=$DAT")
fi

# Collect labels for reporting
LABELS=()
if (( N_LEGS > 0 )); then
  for ((k=1; k<=N_LEGS; k++)); do
    dat_f="$(awk '$1=="CLIMATE_INDEX"{print $2; exit}' "lt.exe.resp.${k}")"
    LABELS+=("$(leg_label "$dat_f")")
  done
  echo "=== Multi-Leg GEM-LTE (${N_LEGS} legs: ${LABELS[*]}) | TIDES=$TIDES_MODE ($TIER) | window=$T_START-$T_END (EXCLUDE=$EXC) | mode=$MODE ==="
else
  echo "=== Single-Leg GEM-LTE: $DAT | TIDES=$TIDES_MODE ($TIER) | window=$T_START-$T_END (EXCLUDE=$EXC) | mode=$MODE ==="
fi

if [[ "$MODE" == "homotopy" || "$MODE" == "refine" ]]; then
  HOMO=$([[ "$MODE" == "homotopy" ]] && echo "TRUE" || echo "FALSE")
  env "${COMMON_ENV[@]}" \
    DLOD_REF="$DLOD_MODE" \
    REFINE=TRUE \
    REFINE_HOMOTOPY="$HOMO" \
    REFINE_LAM="$R_LAM" \
    REFINE_HOMOTOPY_LAM="$R_HLAM" \
    REFINE_ROUNDS="$R_ROUNDS" \
    REFINE_ITERS="$R_ITERS" \
    REFINE_APPLY=TRUE \
    "$LT_EXE" | tee refine.log
  echo "--- Scoring refined lt.exe.p (TEST_ONLY=TRUE) ---"
  env "${COMMON_ENV[@]}" DLOD_REF=TRUE REFINE=FALSE TEST_ONLY=TRUE "$LT_EXE" > score.log
elif [[ "$MODE" == "search" ]]; then
  env "${COMMON_ENV[@]}" DLOD_REF="$DLOD_MODE" REFINE=FALSE TIMEOUT="$TOUT" "$LT_EXE" | tee search.log
  env "${COMMON_ENV[@]}" DLOD_REF=TRUE REFINE=FALSE TEST_ONLY=TRUE "$LT_EXE" > score.log
else
  env "${COMMON_ENV[@]}" DLOD_REF=TRUE REFINE=FALSE TEST_ONLY=TRUE "$LT_EXE" | tee score.log
fi

# Concise Train / Held-Out / Full CC summary for both single-leg and multi-leg runs
python3 - "$T_START" "$T_END" "$EXC" "$N_LEGS" "${LABELS[@]}" <<'PY'
import sys, numpy as np
from pathlib import Path
t0, t1 = float(sys.argv[1]), float(sys.argv[2])
exc = sys.argv[3].upper() == "TRUE"
n_legs = int(sys.argv[4])
labels = sys.argv[5:]

def summarize(csv_path, name):
    arr = np.loadtxt(csv_path, delimiter=",")
    t, y, m = arr[:, 0], arr[:, 1], arr[:, 2]
    v = y != 0.0
    win = v & (t >= t0) & (t <= t1)
    tr = (v & (~win)) if exc else win
    ho = win if exc else (v & (~win))
    cc_tr = float(np.corrcoef(y[tr], m[tr])[0, 1]) if tr.sum() > 2 else float("nan")
    cc_ho = float(np.corrcoef(y[ho], m[ho])[0, 1]) if ho.sum() > 2 else float("nan")
    cc_all = float(np.corrcoef(y[v], m[v])[0, 1])
    print(f"  {name:<10s} | Train CC = {cc_tr:+.4f} (N={int(tr.sum())}) | Held-Out CC = {cc_ho:+.4f} (N={int(ho.sum())}) | Full CC = {cc_all:+.4f}")
    return cc_tr, cc_ho

print("\nSUMMARY:")
if n_legs > 0:
    trs, hos = [], []
    for k in range(1, n_legs + 1):
        tr, ho = summarize(Path(f"lte_results.{k}.csv"), labels[k - 1] if k - 1 < len(labels) else f"Leg{k}")
        trs.append(tr); hos.append(ho)
    print(f"  {'MEAN':<10s} | Train CC = {np.mean(trs):+.4f}          | Held-Out CC = {np.mean(hos):+.4f}")
else:
    summarize(Path("lte_results.csv"), "SINGLE")
PY

if [[ "$DO_PLOT" == "true" && -e "$PLOT_PY" ]]; then
  python3 "$PLOT_PY" . --year 2000 --years 3 --measured --no-show --out lpap_year.png
fi
