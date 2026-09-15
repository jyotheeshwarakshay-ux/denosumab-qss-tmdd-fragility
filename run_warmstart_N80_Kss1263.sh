#!/bin/bash
# =============================================================================
# N=80 WARM-START DIAGNOSTIC DRIVER: 3 fits, Kss=1.263 fixed, N=80 corrected
# data, nBurn=500/nEm=1000, each warm-started from a different known N=30
# good-basin vector (repro seeds 42, 999, 2024).
#
# Isolated test -- writes ONLY to warmstart_N80_Kss1263.csv and
# warmstart_N80_Kss1263_logs/. Reads gate_test_N30_Kss1263_repro.csv but never
# writes to it. Does not touch any other existing results file.
#
# LAUNCH: caffeinate -i -s ./run_warmstart_N80_Kss1263.sh > <logfile> 2>&1
# =============================================================================
set -uo pipefail

REPO="$HOME/pharmacometrics/denosumab-tmdd-qss"
SCRIPT="$REPO/warmstart_N80_Kss1263.R"
LOGDIR="$REPO/warmstart_N80_Kss1263_logs"
mkdir -p "$LOGDIR"

SOURCE_SEEDS="42 999 2024"
PAUSE_SECONDS=15

cd "$REPO" || { echo "Cannot cd to $REPO"; exit 1; }

total=3
i=0

for seed in $SOURCE_SEEDS; do
  i=$((i + 1))
  LOG="$LOGDIR/source_seed${seed}.log"
  echo "[$i/$total] source_seed=$seed -> $LOG"

  Rscript "$SCRIPT" "$seed" > "$LOG" 2>&1
  rc=$?

  tail -n 2 "$LOG"

  if [ $rc -ne 0 ]; then
    echo "  !! Rscript exited with code $rc for source_seed=$seed -- see $LOG"
  fi

  if grep -q "^SKIP" "$LOG"; then
    :
  else
    echo "  (pausing ${PAUSE_SECONDS}s to let memory settle)"
    sleep "$PAUSE_SECONDS"
  fi
done

echo ""
echo "ALL $total WARM-START FITS PROCESSED."
echo "Results: $REPO/warmstart_N80_Kss1263.csv"
echo "Per-fit logs: $LOGDIR/"
