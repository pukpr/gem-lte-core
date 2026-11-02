#!/usr/bin/env bash
# NINO4 alone, tides locked (LOCKT). Stage 1: fast windings removed (no "harm", NH "") so delA/delB/asym/delC/impC/mb/IR...
# can move on the slow winding alone. Stage 2: harmonics re-added (harm from harm_saved.json, NH "1 1 1") and refined.
cd "$(dirname "$0")"; ulimit -s 65536
E="CLIMATE_INDEX=nino4.dat DLOD_REF=TRUE EXCLUDE=true TRAIN_START=2000 TRAIN_END=2005 F9=1 LOCKT=TRUE LOCKA=FALSE METRIC=CC TREND=true"
S1=${1:-900}; S2=${2:-900}
cp stage1_seed.p lt.exe.p
env $E TEST_ONLY=TRUE TIMEOUT=60 NUMBER_OF_PROCESSORS=2 ../lt.exe > s1_seed_test.out 2>&1; cp stage1_seed.p lt.exe.p
env $E TEST_ONLY=FALSE TIMEOUT=$S1 NUMBER_OF_PROCESSORS=16 ../lt.exe > s1_run.out 2>&1; cp lt.exe.p stage1_fit.p
python3 -c "
import json;p=json.load(open('stage1_fit.p'));p['harm']=json.load(open('harm_saved.json'))['harm'];json.dump(p,open('stage2_seed.p','w'),indent=1)"
sed -i 's/^NH .*/NH                             "1 1 1"/' lt.exe.resp
cp stage2_seed.p lt.exe.p
env $E TEST_ONLY=TRUE TIMEOUT=60 NUMBER_OF_PROCESSORS=2 ../lt.exe > s2_seed_test.out 2>&1; cp stage2_seed.p lt.exe.p
env $E TEST_ONLY=FALSE TIMEOUT=$S2 NUMBER_OF_PROCESSORS=16 ../lt.exe > s2_run.out 2>&1; cp lt.exe.p stage2_fit.p
for f in s1_seed_test s1_run s2_seed_test s2_run; do printf "%-14s %s\n" $f "$(grep '^CC' $f.out | tail -1)"; done
