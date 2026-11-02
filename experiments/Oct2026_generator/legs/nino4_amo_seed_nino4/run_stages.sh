#!/usr/bin/env bash
# Stage 1: slow winding only (no harm, NH ""), 45 min. Stage 2: harm [15,2] + NH "1 1", 30 min.
cd "$(dirname "$0")"
cp lt.exe.resp resp_harm
python3 - <<'PY'
import json
p = json.load(open('seed.p')); p.pop('harm', None)
json.dump(p, open('lt.exe.p', 'w'), indent=2)
PY
sed 's/^NH .*/NH                             ""/' resp_harm > lt.exe.resp
./lte_run.sh 2700 16 > s1_run.out 2>&1
cp lt.exe.p s1_fit.p; cp legs.png s1_legs.png 2>/dev/null
python3 - <<'PY'
import json
p = json.load(open('s1_fit.p')); p['harm'] = [15.0, 2.0]
json.dump(p, open('lt.exe.p', 'w'), indent=2)
PY
cp resp_harm lt.exe.resp
./lte_run.sh 1800 16 > s2_run.out 2>&1
cp lt.exe.p s2_fit.p
echo STAGES_DONE
