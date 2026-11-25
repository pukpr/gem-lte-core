#!/usr/bin/env bash
# 4 legs on one manifold: NAO (leg 1), Baltic (leg 2), AMO (leg 3), TNA (leg 4).
# 1960-1970 is held out on every leg. See LEGS.md.
# usage: ./lte_run.sh [seconds, default 18000] [threads, default 16]       TEST_ONLY=TRUE ./lte_run.sh  scores only
cd "$(dirname "$0")"
export DLOD_REF=TRUE EXCLUDE=true TRAIN_START=1960 TRAIN_END=1970 F9=1 LOCKA=FALSE LOCKT=FALSE METRIC=CC TREND=true ACCEL=TRUE
# export DLOD_MIX=MEAN DLOD_SHIFT=0.5
export MANIFOLD_PRE=TRUE
# export EXCLUDE=FALSE TRAIN_START=1880 TRAIN_END=1960 VALIDATE=TRUE VALIDATE_START=1970
export TIMEOUT=${1:-18000} NUMBER_OF_PROCESSORS=${2:-16} TEST_ONLY=${TEST_ONLY:-FALSE}
ulimit -s 65536 2>/dev/null || true
[ -x ../lt.exe ] && [ "$(wc -c < ../lt.exe)" -gt 1024 ] || ln -sf ../../Feb2026/enso_opt ../lt.exe 2>/dev/null || cp -f ../../../Feb2026/enso_opt ../lt.exe
[ -s ../dlod3.dat ] && [ "$(wc -c < ../dlod3.dat)" -gt 1024 ] || ln -sf ../../Feb2026/dlod3.dat ../dlod3.dat 2>/dev/null || cp -f ../../../Feb2026/dlod3.dat ../dlod3.dat
../lt.exe
python3 ../leg_scores2.py . NAO Baltic AMO TNA --held 1960 1970 --post
