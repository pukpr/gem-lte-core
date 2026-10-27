#!/usr/bin/env bash
# 9 legs on one manifold: KIT70 (leg 1), KIT50 (leg 2), KIT40 (leg 3), KIT30 (leg 4), KIT20 (leg 5), KIT15 (leg 6), ERA30 (leg 7), ERA50 (leg 8), ERA70 (leg 9); seed qbo9_kit_era_ex2015_2021/result.p.
# 2015-2021 is held out on every leg. Made by make_leg_run.py; see LEGS.md.
# usage: ./lte_run.sh [seconds, default 1800] [threads, default 16]       TEST_ONLY=TRUE ./lte_run.sh  scores only
cd "$(dirname "$0")"
export DLOD_REF=TRUE EXCLUDE=true TRAIN_START=2015 TRAIN_END=2021 F9=1 LOCKA=FALSE LOCKT=FALSE METRIC=CC TREND=true BZONE=TRUE IDATE=1920.9 TIDE_LEGACY=TRUE
export TIMEOUT=${1:-1800} NUMBER_OF_PROCESSORS=${2:-16} TEST_ONLY=${TEST_ONLY:-FALSE}
ulimit -s 65536
../lt.exe
python3 ../leg_scores2.py . KIT70 KIT50 KIT40 KIT30 KIT20 KIT15 ERA30 ERA50 ERA70 --held 2015 2021
