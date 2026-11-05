#!/usr/bin/env bash
# 1 legs on one manifold: NINO4 (leg 1); seed nino4_slow_then_harm/benchmark_enso_2026-11-01/lt.exe.p.
# 2000-2005 is held out on every leg. Made by make_leg_run.py; see LEGS.md.
# usage: ./lte_run.sh [seconds, default 1800] [threads, default 16]       TEST_ONLY=TRUE ./lte_run.sh  scores only
cd "$(dirname "$0")"
export DLOD_REF=TRUE EXCLUDE=true TRAIN_START=2000 TRAIN_END=2005 F9=1 LOCKA=FALSE LOCKT=FALSE METRIC=CP TREND=true FB_LO=0.125 FB_HI=0.6667
export TIMEOUT=${1:-1800} NUMBER_OF_PROCESSORS=${2:-16} TEST_ONLY=${TEST_ONLY:-FALSE}
ulimit -s 65536
../lt.exe
ln -sf lte_results.csv lte_results.1.csv; ln -sf lt.exe.windings.json lt.exe.windings.1.json
python3 ../leg_scores2.py . NINO4 --held 2000 2005
