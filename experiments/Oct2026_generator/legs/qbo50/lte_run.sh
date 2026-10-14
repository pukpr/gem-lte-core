#!/usr/bin/env bash
# QBO 50 hPa, two legs on the Feb2026/qbo50 configuration (see LEGS.md):
#   leg 1  KIT monthly (qbo50.dat)            leg 2  ERA5 5-day means (era_qbo_u50_5day.dat, from ../../QBO/getERA5_daily_50.py)
# 2000-2005 held out on both. TIDE_LEGACY is needed: the stored parameters predate the 2026-10-06 tidal-sum change.
# usage: ./lte_run.sh [seconds, default 1800] [threads, default 16]       TEST_ONLY=TRUE ./lte_run.sh  scores only
cd "$(dirname "$0")"
[ -f era_qbo_u50_5day.dat ] || { echo "era_qbo_u50_5day.dat is missing: run ../../QBO/getERA5_daily_50.py and copy it here"; exit 1; }
export BZONE=TRUE DLOD_REF=TRUE EXCLUDE=true TRAIN_START=2000 TRAIN_END=2005 F9=1 IDATE=1920.9 LOCKA=FALSE LOCKT=FALSE METRIC=CC TREND=true TIDE_LEGACY=TRUE
export TIMEOUT=${1:-1800} NUMBER_OF_PROCESSORS=${2:-16} TEST_ONLY=${TEST_ONLY:-FALSE}
ulimit -s 65536
../lt.exe
python3 ../leg_scores2.py . "KIT monthly 50 hPa" "ERA 5-day 50 hPa" --from 1950
python3 ../harmonic_spectrum.py . "KIT monthly 50 hPa" "ERA 5-day 50 hPa"
