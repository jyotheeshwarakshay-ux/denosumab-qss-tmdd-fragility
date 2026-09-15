#!/bin/bash
# =============================================================================
# REPRODUCIBILITY PROFILE DRIVER (remaining 6 grid points): completes the
# 7-point x 5-seed = 35-fit profile on the independent repro dataset (seed
# 2025), same SAEM settings as the pilot (nBurn=150/nEm=200). Grid points
# 1,2,3,5,6,7 x seeds 42,123,999,7,2024 = 30 fits. Grid point 4 (Kss=1.263) is
# already done (the 5-fit gate test) and is skipped by the worker's own
# resume-check against gate_test_N30_Kss1263_repro.csv.
#
# APPENDS to gate_test_N30_Kss1263_repro.csv (already has 5 rows, backfilled
# with grid_point_index=4) so the final file has all 35 rows.
#
# Does NOT touch multistart_results_corrected_pilot.csv,
# multistart_results_corrected_full.csv, convergence_test_N80_Kss1263.csv, or
# multistart_results_corrected_repro.csv.
#
# LAUNCH: caffeinate -i -s ./run_gate_test_N30_repro_remaining.sh > <logfile> 2>&1
# =============================================================================
set -uo pipefail

REPO="$HOME/pharmacometrics/denosumab-tmdd-qss"
SCRIPT="$REPO/gate_test_N30_repro.R"
LOGDIR="$REPO/gate_test_N30_repro_logs"
mkdir -p "$LOGDIR"

GRID_POINTS="1 2 3 5 6 7"
SEEDS="42 123 999 7 2024"
PAUSE_SECONDS=15

cd "$REPO" || { echo "Cannot cd to $REPO"; exit 1; }

total=30
i=0

for gp in $GRID_POINTS; do
  for seed in $SEEDS; do
    i=$((i + 1))
    LOG="$LOGDIR/gp${gp}_seed${seed}.log"
    echo "[$i/$total] grid_point=$gp seed=$seed -> $LOG"

    Rscript "$SCRIPT" "$gp" "$seed" > "$LOG" 2>&1
    rc=$?

    tail -n 2 "$LOG"

    if [ $rc -ne 0 ]; then
      echo "  !! Rscript exited with code $rc for grid_point=$gp seed=$seed -- see $LOG"
    fi

    if grep -q "^SKIP" "$LOG"; then
      :
    else
      echo "  (pausing ${PAUSE_SECONDS}s to let memory settle)"
      sleep "$PAUSE_SECONDS"
    fi
  done
done

echo ""
echo "ALL $total REMAINING (grid_point, seed) COMBINATIONS PROCESSED."
echo "Results (appended): $REPO/gate_test_N30_Kss1263_repro.csv"
echo "Per-fit logs: $LOGDIR/"
