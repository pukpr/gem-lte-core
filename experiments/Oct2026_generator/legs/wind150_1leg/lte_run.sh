#!/usr/bin/env bash
# 1 legs on one manifold: WIND150 (leg 1); seed wind_soi_src/seed.p.
# 2011-2016 is held out on every leg. Made by make_leg_run.py; see LEGS.md.
# usage: ./lte_run.sh [seconds, default 1800] [threads, default 16]       TEST_ONLY=TRUE ./lte_run.sh  scores only
cd "$(dirname "$0")"
export DLOD_REF=TRUE EXCLUDE=true TRAIN_START=2011 TRAIN_END=2016 F9=1 LOCKA=FALSE LOCKT=FALSE METRIC=CC TREND=true 
export TIMEOUT=${1:-1800} NUMBER_OF_PROCESSORS=${2:-16} TEST_ONLY=${TEST_ONLY:-FALSE}
ulimit -s 65536
../lt.exe
python3 ../leg_scores2.py . WIND150 --held 2011 2016
