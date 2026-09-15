#!/bin/bash
# =============================================================================
# DRIVER: 7-point Kss profile x 3 seeds = 21 fits, PILOT on the CORRECTED model.
# One fresh Rscript process at a time, 15s pause between fits.
#
# Per FLAG_STATUS.md pilot-gate design: 30 subjects, all 7 grid points kept,
# 3 seeds/point (shape check, not final numbers).
#
# RESUMABLE: multistart_profile_Kss_corrected_pilot.R checks
# multistart_results_corrected_pilot.csv before doing any work and skips a
# (grid_point, seed) combo already recorded.
# =============================================================================
set -uo pipefail

REPO="$HOME/pharmacometrics/denosumab-tmdd-qss"
SCRIPT="$REPO/multistart_profile_Kss_corrected_pilot.R"
LOGDIR="$REPO/multistart_logs_corrected_pilot"
mkdir -p "$LOGDIR"

SEEDS="42 123 999"
N_POINTS=7
PAUSE_SECONDS=15

cd "$REPO" || { echo "Cannot cd to $REPO"; exit 1; }

total=$((N_POINTS * 3))
i=0

for gp in $(seq 1 "$N_POINTS"); do
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
echo "ALL $total (grid_point, seed) COMBINATIONS PROCESSED."
echo "Results: $REPO/multistart_results_corrected_pilot.csv"
echo "Per-fit logs: $LOGDIR/"
