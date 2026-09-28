#!/bin/bash
cd /home/paul/eval/gem-lte-core/experiments/Sep2026
last_line=$(wc -l < refix.log)
while true; do
  total=$(wc -l < refix.log)
  if [ "$total" -gt "$last_line" ]; then
    tail -n +"$((last_line+1))" refix.log | grep -E "REFIXED|STILL BAD|EXCEPTION|resume-refix complete"
    last_line=$total
  fi
  if [ -f RESUME_REFIX_DONE ]; then
    echo "DONE: $(cat RESUME_REFIX_DONE)"
    break
  fi
  if ! pgrep -f resume_refix.py > /dev/null; then
    echo "PROCESS EXITED without writing RESUME_REFIX_DONE -- check refix.log for a crash"
    break
  fi
  sleep 30
done
