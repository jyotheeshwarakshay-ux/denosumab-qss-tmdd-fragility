#!/bin/bash
# =============================================================================
# OFV-AT-TRUTH, PILOT DATASET, REPRODUCIBILITY DRIVER: 3 repeats, each a
# SEPARATE `Rscript` process launch, mirroring run_ofv_at_truth_repro.sh.
# Cross-dataset check: does the reproducibility dataset's OFV-at-truth
# result (truth falls in the high-OFV bucket) replicate on the independent
# pilot N=30 draw?
#
# Isolated -- writes ONLY to ofv_at_truth_pilot.csv.
#
# LAUNCH (when authorized): caffeinate -i -s ./run_ofv_at_truth_pilot.sh > <logfile> 2>&1
# =============================================================================
set -uo pipefail

REPO="$HOME/pharmacometrics/denosumab-tmdd-qss"
SCRIPT="$REPO/ofv_at_truth_pilot.R"

cd "$REPO" || { echo "Cannot cd to $REPO"; exit 1; }

REPS="1 2 3"
total=3
i=0

for rep in $REPS; do
  i=$((i + 1))
  echo "[$i/$total] rep=$rep"

  Rscript "$SCRIPT" "$rep"
  rc=$?

  if [ $rc -ne 0 ]; then
    echo "  !! Rscript exited with code $rc for rep=$rep"
  fi

  echo "  (pausing 10s to let memory settle)"
  sleep 10
done

echo ""
echo "ALL $total REPS PROCESSED."
echo "Results: $REPO/ofv_at_truth_pilot.csv"
