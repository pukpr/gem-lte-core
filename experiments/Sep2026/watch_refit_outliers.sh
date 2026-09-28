#!/bin/bash
cd /home/paul/eval/gem-lte-core/experiments/Sep2026
last_line=$(wc -l < refit_outliers.log)
while true; do
  total=$(wc -l < refit_outliers.log)
  if [ "$total" -gt "$last_line" ]; then
    tail -n +"$((last_line+1))" refit_outliers.log | grep -E "^===|coherence (BEFORE|AFTER)|cascade (FAILED|succeeded)|COMMITTED|rolling back|SUMMARY"
    last_line=$total
  fi
  if ! pgrep -f refit_manifold_outliers.py > /dev/null; then
    echo "PROCESS EXITED (check refit_outliers.log tail for SUMMARY / crash)"
    break
  fi
  sleep 15
done
