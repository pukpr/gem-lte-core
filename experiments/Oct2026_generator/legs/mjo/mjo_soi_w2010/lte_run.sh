#!/usr/bin/env bash
# Two legs: MJO pentad index at 140E (CPC) and 5-day inverted SOI (Long Paddock), same 73-per-year grid, 2010-2015 held out; SOI seed. Measured lead/lag between them at 140E is zero (strongest at 0 days), so neither is shifted.
cd "$(dirname "$0")"
export DLOD_REF=TRUE F9=1 IDATE=1880.0 LOCKA=FALSE LOCKT=FALSE METRIC=CC TREND=false EXCLUDE=true TRAIN_START=2010 TRAIN_END=2015 CLIMATE_INDEX=mjo_140E_5.dat
export TIMEOUT=${1:-1800} NUMBER_OF_PROCESSORS=${2:-16} TEST_ONLY=${TEST_ONLY:-FALSE}
ulimit -s 65536
../../lt.exe
python3 ../../leg_scores2.py . "MJO 140E" "SOI 5-day" --from 1975 --held 2010 2015
