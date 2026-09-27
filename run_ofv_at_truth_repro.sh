#!/bin/bash
# =============================================================================
# OFV-AT-TRUTH REPRODUCIBILITY DRIVER: 4 repeats, each a SEPARATE `Rscript`
# process launch (no shared in-memory state between repeats), matching
# run_seed_determinism_forensic.sh's design. Evaluates the FOCEi EBE-only
# objective at the fully-true parameter vector on the reproducibility
# dataset, 4 times, to confirm the number (4746.271, from an initial single
# run) is actually deterministic before it's used to interpret Experiment C.
#
# Isolated -- writes ONLY to ofv_at_truth_repro.csv.
#
# LAUNCH (when authorized): caffeinate -i -s ./run_ofv_at_truth_repro.sh > <logfile> 2>&1
# =============================================================================
set -uo pipefail

REPO="$HOME/pharmacometrics/denosumab-tmdd-qss"
SCRIPT="$REPO/ofv_at_truth_repro.R"

cd "$REPO" || { echo "Cannot cd to $REPO"; exit 1; }

REPS="1 2 3 4"
total=4
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
echo "Results: $REPO/ofv_at_truth_repro.csv"
