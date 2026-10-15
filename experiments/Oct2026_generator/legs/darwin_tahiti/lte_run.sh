#!/usr/bin/env bash
# Two legs: Darwin and Tahiti monthly sea-level pressure (CRU soi_dar.dat, soi_tah.dat), each with its own monthly climatology (outside 2000-2005) removed, full records (Darwin from 1866, Tahiti from 1861, after its 1858-1860 gaps); the manifold is integrated from 1860 with init still defined at 1880 --
# 2000-2005 held
# out. Seed: the SOI two-leg result (legs/soi_5_m). usage: ./lte_run.sh [seconds] [threads]; TEST_ONLY=TRUE to score only
cd "$(dirname "$0")"
export DLOD_REF=TRUE F9=1 IDATE=1860.0 LOCKA=FALSE LOCKT=FALSE METRIC=CC TREND=true EXCLUDE=true TRAIN_START=2000 TRAIN_END=2005
export TIMEOUT=${1:-1800} NUMBER_OF_PROCESSORS=${2:-16} TEST_ONLY=${TEST_ONLY:-FALSE}
ulimit -s 65536
../lt.exe
python3 ../leg_scores2.py . "Darwin" "Tahiti" --from 1950
python3 deseason.py
