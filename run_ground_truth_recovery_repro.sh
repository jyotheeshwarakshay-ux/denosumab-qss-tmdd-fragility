#!/bin/bash
# =============================================================================
# GROUND-TRUTH RECOVERY DRIVER (Experiment C): 5 fits, Kss FIXED at its true
# simulation value (1.56), same 5 standard seeds used throughout this repo
# (42, 123, 999, 7, 2024), on the reproducibility N=30 corrected dataset,
# nBurn=150/nEm=200 (matching the main profile analysis).
#
# Isolated -- writes ONLY to ground_truth_recovery_repro.csv and
# ground_truth_recovery_repro_logs/. Does not touch any existing reported
# results file.
#
# LAUNCH (when authorized): caffeinate -i -s ./run_ground_truth_recovery_repro.sh > <logfile> 2>&1
# =============================================================================
set -uo pipefail

REPO="$HOME/pharmacometrics/denosumab-tmdd-qss"
SCRIPT="$REPO/ground_truth_recovery_repro.R"
LOGDIR="$REPO/ground_truth_recovery_repro_logs"
mkdir -p "$LOGDIR"

SEEDS="42 123 999 7 2024"

cd "$REPO" || { echo "Cannot cd to $REPO"; exit 1; }

total=5
i=0

for seed in $SEEDS; do
  i=$((i + 1))
  LOG="$LOGDIR/seed${seed}.log"
  echo "[$i/$total] seed=$seed -> $LOG"

  Rscript "$SCRIPT" "$seed" > "$LOG" 2>&1
  rc=$?

  tail -n 2 "$LOG"

  if [ $rc -ne 0 ]; then
    echo "  !! Rscript exited with code $rc for seed=$seed -- see $LOG"
  fi

  if grep -q "^SKIP" "$LOG"; then
    :
  else
    echo "  (pausing 15s to let memory settle)"
    sleep 15
  fi
done

echo ""
echo "ALL $total SEEDS PROCESSED."
echo "Results: $REPO/ground_truth_recovery_repro.csv"
echo "Per-fit logs: $LOGDIR/"
