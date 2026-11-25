#!/usr/bin/env bash
# 4 legs on one manifold: AMO (leg 1), PDO (leg 2), IODE (leg 3), NINO4 (leg 4); seed amo_pdo_feb_delCm/lt.exe.p.
# 1960-1970 is held out on every leg. Made by make_leg_run.py; see LEGS.md.
# usage: ./lte_run.sh [seconds, default 1800] [threads, default 16]       TEST_ONLY=TRUE ./lte_run.sh  scores only
cd "$(dirname "$0")"
export DLOD_REF=TRUE EXCLUDE=true TRAIN_START=1960 TRAIN_END=1970 F9=1 LOCKA=FALSE LOCKT=FALSE METRIC=CC TREND=true BZONE=TRUE ZONE=TRUE IDATE=1920.9 ACCEL=TRUE
# export DLOD_MIX=MEAN DLOD_SHIFT=0.5
export MANIFOLD_PRE=TRUE
# export EXCLUDE=FALSE TRAIN_START=1880 TRAIN_END=1960 VALIDATE=TRUE VALIDATE_START=1970
export TIMEOUT=${1:-18000} NUMBER_OF_PROCESSORS=${2:-16} TEST_ONLY=${TEST_ONLY:-FALSE}
ulimit -s 65536
../lt.exe
python3 ../leg_scores2.py . NAO NINO4 AMO PDO IODE TNA --held 1960 1970 
