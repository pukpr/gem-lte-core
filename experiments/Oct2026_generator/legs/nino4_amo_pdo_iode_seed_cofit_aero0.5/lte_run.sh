#!/usr/bin/env bash
# 4 legs on one manifold: NINO4 (leg 1), AMO (leg 2), PDO (leg 3), IODE (leg 4); seed nino4_amo_pdo_seed_cofit/s2_fit.p.
# 2000-2005 is held out on every leg. Made by make_leg_run.py; see LEGS.md.
# usage: ./lte_run.sh [seconds, default 1800] [threads, default 16]       TEST_ONLY=TRUE ./lte_run.sh  scores only
cd "$(dirname "$0")"
export AERO=/home/paul/eval/gem-lte-core/experiments/ext_data/aero_saod.dat AERO_LAG=0.5 DLOD_REF=TRUE EXCLUDE=true TRAIN_START=2000 TRAIN_END=2005 F9=1 LOCKA=FALSE LOCKT=FALSE METRIC=CC TREND=true 
export TIMEOUT=${1:-1800} NUMBER_OF_PROCESSORS=${2:-16} TEST_ONLY=${TEST_ONLY:-FALSE}
ulimit -s 65536
../lt.exe
python3 ../leg_scores2.py . NINO4 AMO PDO IODE --held 2000 2005
