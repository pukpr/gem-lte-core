#!/usr/bin/env bash
# Two legs: CRU monthly SOI (1866 on) and the 5-day SOI from Long Paddock daily data (1991 on, inverted), 2000-2005 held out.

# Seed: the Darwin/Tahiti result (legs/darwin_tahiti), deseasonalized, integration from 1860. usage: ./lte_run.sh [seconds] [threads]; TEST_ONLY=TRUE to score only
cd "$(dirname "$0")"
export DLOD_REF=TRUE F9=1 IDATE=1860.0 LOCKA=FALSE LOCKT=FALSE METRIC=CC TREND=true EXCLUDE=true TRAIN_START=2000 TRAIN_END=2005
export TIMEOUT=${1:-1800} NUMBER_OF_PROCESSORS=${2:-16} TEST_ONLY=${TEST_ONLY:-FALSE}
ulimit -s 65536
../lt.exe
python3 ../leg_scores2.py . "CRU SOI monthly" "SOI 5-day" --from 1950
python3 high_harm.py
