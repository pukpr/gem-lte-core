#!/usr/bin/env bash
# 5-day ERA leg + monthly KIT leg, QBO 70 hPa; 2000-2005 held out on both (see LEGS.md). Leg settings are in lt.exe.resp.1 and lt.exe.resp.2.
# usage: ./lte_run.sh [seconds, default 1800] [threads, default 16]       TEST_ONLY=TRUE ./lte_run.sh  scores only
cd "$(dirname "$0")"
export DLOD_REF=TRUE F9=1 IDATE=1880.0 LOCKA=FALSE LOCKT=FALSE METRIC=DTW TREND=false ZONE=FALSE
export TIMEOUT=${1:-1800} NUMBER_OF_PROCESSORS=${2:-16} TEST_ONLY=${TEST_ONLY:-FALSE}
ulimit -s 65536
../lt.exe
python3 ../leg_scores2.py . "ERA 5-day" "KIT monthly" --from 1950
