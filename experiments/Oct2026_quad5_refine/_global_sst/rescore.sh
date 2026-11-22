#!/usr/bin/env bash
# usage: rescore.sh QUAD F9 OUTDIR -- TEST_ONLY of the quad's centre leg alone (no LEGS) from its refined lt.exe.p
q=$1; f9=$2; out=$3; src=$(dirname "$0")/../$q
mkdir -p "$out"; cp "$src/lt.exe.p" "$src/$q.dat" "$out/"
grep -vE '^(LEGS|-?LEG_)' "$src/lt.exe.resp" > "$out/lt.exe.resp"; echo "CLIMATE_INDEX                  $q.dat" >> "$out/lt.exe.resp"
ln -sfn /home/paul/eval/gem-lte-core/experiments/Feb2026/enso_opt "$out/lt.exe"; [ -e "$out/../dlod3.dat" ] || ln -s /home/paul/eval/gem-lte-core/experiments/Feb2026/dlod3.dat "$out/../dlod3.dat"
cd "$out"; env DLOD_REF=TRUE EXCLUDE=true TRAIN_START=2000 TRAIN_END=2005 F9=$f9 LOCKA=FALSE LOCKT=FALSE METRIC=CC TREND=true \
  TEST_ONLY=TRUE TIMEOUT=60 NUMBER_OF_PROCESSORS=2 timeout 300 bash -c 'ulimit -s 65536; ./lt.exe' > run.out 2>&1; echo "$q exit $?"
