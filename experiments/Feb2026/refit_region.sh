#!/usr/bin/env bash
# refit_region.sh <region> -- redo a single k_sst grid region's DTW+CC fit
# with CLIMATE_INDEX correctly set (fixing the bug where it silently
# defaulted to nino4.dat because the region's .dat filename was only
# ever passed as an ignored positional argument, never as CLIMATE_INDEX).
set -e
region="$1"
cd "/home/paul/eval/gem-lte-core/experiments/Feb2026/$region"
ulimit -s unlimited
export CLIMATE_INDEX="$region.dat"
METRIC=DTW timeout -s INT 45 ../lt.exe > dtw_refit.log 2>&1 || true
METRIC=CC timeout -s INT 30 ../lt.exe > cc_refit.log 2>&1 || true
echo "done: $region"
