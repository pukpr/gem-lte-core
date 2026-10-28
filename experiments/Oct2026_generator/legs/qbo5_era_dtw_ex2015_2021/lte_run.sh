#!/usr/bin/env bash
# 5 legs on one manifold: ERA10 (leg 1), ERA20 (leg 2), ERA30 (leg 3), ERA50 (leg 4), ERA70 (leg 5); seed qbo9_kit_era70full_ex2015_2021/result.p.
# 2015-2021 is held out on every leg. Made by make_leg_run.py; see LEGS.md.
# usage: ./lte_run.sh [seconds, default 1800] [threads, default 16]       TEST_ONLY=TRUE ./lte_run.sh  scores only
cd "$(dirname "$0")"
export DLOD_REF=TRUE EXCLUDE=true TRAIN_START=2015 TRAIN_END=2021 F9=1 LOCKA=FALSE LOCKT=FALSE METRIC=DTW TREND=true BZONE=TRUE IDATE=1920.9 TIDE_LEGACY=TRUE
export TIMEOUT=${1:-1800} NUMBER_OF_PROCESSORS=${2:-16} TEST_ONLY=${TEST_ONLY:-FALSE}
ulimit -s 65536
../lt.exe
python3 ../leg_scores2.py . ERA10 ERA20 ERA30 ERA50 ERA70 --held 2015 2021
