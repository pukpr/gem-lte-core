#!/usr/bin/env bash
# Universal multi-leg LTE launcher with native Ada REFINE support and
# automatic lt.exe.resp.<k> / LEGS bookkeeping.
#
# Usage:
#   ./lte_run.sh [timeout_sec] [threads]
#   ./lte_run.sh --legs nao.dat baltic_3pt_na.dat amo.dat tna.dat --start 1960 --end 1970 --refine --homotopy --lam 1.0
#   ./lte_run.sh --start 1970 --end 1980 --test
#
# Options:
#   --legs f1.dat [f2.dat ...]   Set leg climate indices; auto-writes lt.exe.resp.1..N and updates LEGS in lt.exe.resp
#   --start YEAR                 Set TRAIN_START (default: ${TRAIN_START:-1960})
#   --end YEAR                   Set TRAIN_END   (default: ${TRAIN_END:-1970})
#   --exclude true|false         Set EXCLUDE     (default: ${EXCLUDE:-true})
#   --refine                     Run native Ada REFINE=TRUE, then score with TEST_ONLY=TRUE
#   --homotopy                   Enable REFINE_HOMOTOPY=TRUE (implies --refine)
#   --lam VAL                    Set REFINE_LAM (default: ${REFINE_LAM:-1.0})
#   --hlam VAL                   Set REFINE_HOMOTOPY_LAM (default: same as --lam)
#   --rounds N                   Set REFINE_ROUNDS (default: ${REFINE_ROUNDS:-4})
#   --iters N                    Set REFINE_ITERS  (default: ${REFINE_ITERS:-8})
#   --cold                       Cold start from measured dLOD (DLOD_REF=FALSE)
#   --test                       Score only (TEST_ONLY=TRUE)
#   --timeout SEC                Solver timeout (default: 18000)
#   --threads N                  Worker threads (default: 16)
#   -h, --help                   Show this help
set -euo pipefail
cd "$(dirname "$0")"

# Defaults (preserve env overrides if set)
T_START="${TRAIN_START:-1960}"
T_END="${TRAIN_END:-1970}"
EXC="${EXCLUDE:-true}"
DLOD_MODE="${DLOD_REF:-TRUE}"
RUN_REFINE="${REFINE:-FALSE}"
RUN_HOMOTOPY="${REFINE_HOMOTOPY:-FALSE}"
R_LAM="${REFINE_LAM:-1.0}"
R_HLAM="${REFINE_HOMOTOPY_LAM:--1.0}"
R_ROUNDS="${REFINE_ROUNDS:-4}"
R_ITERS="${REFINE_ITERS:-8}"
TOUT="${TIMEOUT:-18000}"
NPROC="${NUMBER_OF_PROCESSORS:-16}"
T_ONLY="${TEST_ONLY:-FALSE}"
LEGS_FILES=()
POS_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --legs|--dat)
      shift
      while [[ $# -gt 0 && "$1" != -* ]]; do
        LEGS_FILES+=("$1")
        shift
      done
      ;;
    --start|--train-start) T_START="$2"; shift 2 ;;
    --end|--train-end)     T_END="$2"; shift 2 ;;
    --exclude)             EXC="$2"; shift 2 ;;
    --refine)              RUN_REFINE=TRUE; shift ;;
    --homotopy)            RUN_REFINE=TRUE; RUN_HOMOTOPY=TRUE; shift ;;
    --lam)                 R_LAM="$2"; shift 2 ;;
    --hlam|--homotopy-lam) R_HLAM="$2"; shift 2 ;;
    --rounds)              R_ROUNDS="$2"; shift 2 ;;
    --iters)               R_ITERS="$2"; shift 2 ;;
    --cold)                DLOD_MODE=FALSE; shift ;;
    --test)                T_ONLY=TRUE; shift ;;
    --timeout)             TOUT="$2"; shift 2 ;;
    --threads)             NPROC="$2"; shift 2 ;;
    -h|--help)
      sed -n '2,24p' "$0" | sed 's/^# \{0,1\}//'
      exit 0
      ;;
    *)
      POS_ARGS+=("$1"); shift ;;
  esac
done

# Legacy positional compatibility: ./lte_run.sh [timeout] [threads]
if [[ ${#POS_ARGS[@]} -ge 1 ]]; then TOUT="${POS_ARGS[0]}"; fi
if [[ ${#POS_ARGS[@]} -ge 2 ]]; then NPROC="${POS_ARGS[1]}"; fi

# Helper: derive nice display label from .dat filename
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

# If --legs was provided, auto-generate lt.exe.resp.1..N and update LEGS in lt.exe.resp
if [[ ${#LEGS_FILES[@]} -gt 0 ]]; then
  N_NEW=${#LEGS_FILES[@]}
  # Remove stale lt.exe.resp.<k> beyond N_NEW
  for old in lt.exe.resp.[0-9]*; do
    [[ -e "$old" ]] || continue
    idx="${old##*.}"
    if [[ "$idx" =~ ^[0-9]+$ ]] && (( idx > N_NEW )); then
      rm -f "$old"
    fi
  done
  for i in "${!LEGS_FILES[@]}"; do
    k=$((i + 1))
    printf "CLIMATE_INDEX                  %s\n" "${LEGS_FILES[$i]}" > "lt.exe.resp.${k}"
  done
  if [[ -f lt.exe.resp ]]; then
    if grep -qE '^LEGS[[:space:]]+' lt.exe.resp; then
      sed -i -E "s/^LEGS[[:space:]]+[0-9]+/LEGS                           ${N_NEW}/" lt.exe.resp
    else
      printf "LEGS                           %d\n" "${N_NEW}" >> lt.exe.resp
    fi
  fi
fi

# Collect leg labels from lt.exe.resp.1..N (auto-discovered, zero manual list needed)
LABELS=()
for k in 1 2 3 4 5 6 7 8 9; do
  if [[ -f "lt.exe.resp.${k}" ]]; then
    dat_f="$(awk '$1=="CLIMATE_INDEX"{print $2; exit}' "lt.exe.resp.${k}")"
    if [[ -n "$dat_f" ]]; then
      LABELS+=("$(leg_label "$dat_f")")
    fi
  else
    break
  fi
done

export DLOD_REF="$DLOD_MODE" EXCLUDE="$EXC" TRAIN_START="$T_START" TRAIN_END="$T_END"
export F9="${F9:-1}" LOCKA="${LOCKA:-FALSE}" LOCKT="${LOCKT:-FALSE}" METRIC="${METRIC:-CC}" TREND="${TREND:-true}" ACCEL="${ACCEL:-TRUE}"
export MANIFOLD_PRE="${MANIFOLD_PRE:-TRUE}"
export TIMEOUT="$TOUT" NUMBER_OF_PROCESSORS="$NPROC" TEST_ONLY="$T_ONLY"

ulimit -s 65536 2>/dev/null || true
[ -x ../lt.exe ] && [ "$(wc -c < ../lt.exe)" -gt 1024 ] || ln -sf ../../Feb2026/enso_opt ../lt.exe 2>/dev/null || cp -f ../../../Feb2026/enso_opt ../lt.exe
[ -s ../dlod3.dat ] && [ "$(wc -c < ../dlod3.dat)" -gt 1024 ] || ln -sf ../../Feb2026/dlod3.dat ../dlod3.dat 2>/dev/null || cp -f ../../../Feb2026/dlod3.dat ../dlod3.dat

if [[ "${RUN_REFINE^^}" == "TRUE" || "${RUN_REFINE}" == "1" ]]; then
  export REFINE=TRUE REFINE_HOMOTOPY="$RUN_HOMOTOPY" REFINE_LAM="$R_LAM" REFINE_HOMOTOPY_LAM="$R_HLAM"
  export REFINE_ROUNDS="$R_ROUNDS" REFINE_ITERS="$R_ITERS" REFINE_APPLY="${REFINE_APPLY:-TRUE}"
  ../lt.exe
  # Follow up with TEST_ONLY=TRUE to write lte_results.*.csv and score held-out window
  export REFINE=FALSE TEST_ONLY=TRUE DLOD_REF=TRUE
  ../lt.exe
else
  export REFINE=FALSE
  ../lt.exe
fi

if [[ ${#LABELS[@]} -gt 0 && -f ../leg_scores2.py ]]; then
  python3 ../leg_scores2.py . "${LABELS[@]}" --held "$T_START" "$T_END" --post
fi
