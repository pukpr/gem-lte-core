#!/usr/bin/env bash
# Two legs on one manifold: monthly AMO (leg 1) and PDO (leg 2); see LEGS.md. 2000-2005 is held out on both.
# usage: ./lte_run.sh [seconds, default 1800] [threads, default 16]       TEST_ONLY=TRUE ./lte_run.sh  scores only
cd "$(dirname "$0")"
export DLOD_REF=TRUE EXCLUDE=true TRAIN_START=2000 TRAIN_END=2005 F9=1 LOCKA=FALSE LOCKT=FALSE METRIC=CC TREND=true
export TIMEOUT=${1:-1800} NUMBER_OF_PROCESSORS=${2:-16} TEST_ONLY=${TEST_ONLY:-FALSE}
ulimit -s 65536
../lt.exe
python3 ../leg_scores2.py . AMO PDO
