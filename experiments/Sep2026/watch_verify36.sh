#!/bin/bash
cd /home/paul/eval/gem-lte-core/experiments/Sep2026
last_line=$(wc -l < verify_untouched_36.log)
while true; do
  total=$(wc -l < verify_untouched_36.log)
  if [ "$total" -gt "$last_line" ]; then
    tail -n +"$((last_line+1))" verify_untouched_36.log | grep -E "winner=|needs_donor_seed|=== done"
    last_line=$total
  fi
  if [ -f VERIFY_36_DONE ]; then
    echo "DONE: $(cat VERIFY_36_DONE)"
    break
  fi
  if ! pgrep -f verify_untouched_36.py > /dev/null; then
    echo "PROCESS EXITED without writing VERIFY_36_DONE -- check verify_untouched_36.log for a crash"
    break
  fi
  sleep 20
done
