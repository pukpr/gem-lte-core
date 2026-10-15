#!/usr/bin/env bash
# MJO pentad index at 140E (CPC), 5-day grid; seed: soi
cd "$(dirname "$0")"
export DLOD_REF=TRUE F9=1 IDATE=1880.0 LOCKA=FALSE LOCKT=FALSE METRIC=CC TREND=false EXCLUDE=true TRAIN_START=2010 TRAIN_END=2015 CLIMATE_INDEX=mjo_140E_5.dat
export TIMEOUT=${1:-1800} NUMBER_OF_PROCESSORS=${2:-16} TEST_ONLY=${TEST_ONLY:-FALSE}
ulimit -s 65536
../../lt.exe
