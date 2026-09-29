#!/bin/bash
cd /home/paul/eval/gem-lte-core/experiments/Sep2026
last_line=$(wc -l < lockw_batch.log)
while true; do
  total=$(wc -l < lockw_batch.log)
  if [ "$total" -gt "$last_line" ]; then
    tail -n +"$((last_line+1))" lockw_batch.log | grep -E "^===|current min|donor:|cascade (FAILED|succeeded)|new min|COMMITTED|rolling back|SUMMARY|EXCEPTION"
    last_line=$total
  fi
  if [ -f LOCKW_BATCH_DONE ]; then
    echo "DONE"
    break
  fi
  if ! pgrep -f sweep_lockw_batch.py > /dev/null; then
    echo "PROCESS EXITED without LOCKW_BATCH_DONE -- check lockw_batch.log"
    break
  fi
  sleep 20
done
