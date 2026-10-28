#!/bin/bash
# Surrogate null for joint_plateau_inversion.py: same procedure on the daily data circularly shifted by 2-40 years
# (+ a random number of days), identical shift at every level. Results: null/joint_plateau_inversion_shift*.json
cd "$(dirname "$0")"; mkdir -p null
for i in $(seq 1 ${1:-24}); do
  s=$(( (RANDOM % 39 + 2) * 365 + RANDOM % 365 ))
  SHIFT_DAYS=$s python3 joint_plateau_inversion.py > /dev/null && mv joint_plateau_inversion_shift$s.json null/
done
