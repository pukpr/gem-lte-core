#!/usr/bin/env bash
# Stage 1: slow winding only (no harm, NH ""), 45 min. Stage 2: seed harm + NH "1 1", 30 min. 15 threads (3 legs).
cd "$(dirname "$0")"
cp lt.exe.resp resp_harm
python3 - <<'PY'
import json
p = json.load(open('seed.p')); p.pop('harm', None)
json.dump(p, open('lt.exe.p', 'w'), indent=2)
PY
sed 's/^NH .*/NH                             ""/' resp_harm > lt.exe.resp
./lte_run.sh 2700 15 > s1_run.out 2>&1
cp lt.exe.p s1_fit.p; mkdir -p s1_analysis; cp legs.png lte_results.?.csv lt.exe.windings.?.json s1_analysis/ 2>/dev/null
python3 - <<'PY'
import json
p = json.load(open('s1_fit.p')); p['harm'] = json.load(open('seed.p'))['harm']
json.dump(p, open('lt.exe.p', 'w'), indent=2)
PY
cp resp_harm lt.exe.resp
./lte_run.sh 1800 15 > s2_run.out 2>&1
cp lt.exe.p s2_fit.p
echo STAGES_DONE
