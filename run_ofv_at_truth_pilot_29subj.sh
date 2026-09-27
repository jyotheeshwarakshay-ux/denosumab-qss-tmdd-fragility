#!/bin/bash
# =============================================================================
# OFV-AT-TRUTH, PILOT DATASET, 29 SUBJECTS (ID=9 excluded), REPRODUCIBILITY
# DRIVER: 3 repeats, each a SEPARATE `Rscript` process launch, mirroring
# run_ofv_at_truth_pilot.sh. Closes the reproducibility gap for the paper's
# reported 29-subject OFV = 4006.04 (originally computed ad-hoc, never
# scripted).
#
# Isolated -- writes ONLY to ofv_at_truth_pilot_29subj.csv.
#
# LAUNCH (when authorized): caffeinate -i -s ./run_ofv_at_truth_pilot_29subj.sh > <logfile> 2>&1
# =============================================================================
set -uo pipefail

REPO="$HOME/pharmacometrics/denosumab-tmdd-qss"
SCRIPT="$REPO/ofv_at_truth_pilot_29subj.R"

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
echo "Results: $REPO/ofv_at_truth_pilot_29subj.csv"
