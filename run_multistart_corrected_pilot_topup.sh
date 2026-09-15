#!/bin/bash
# =============================================================================
# TOP-UP DRIVER: adds seeds 7 and 2024 (the two the ORIGINAL study used, not
# yet run on the corrected model) at all 7 Kss grid points = 14 additional fits.
# Combined with the 21 already done (seeds 42,123,999), this brings the
# corrected pilot to the IDENTICAL 5-seed set as the original: 42,123,999,7,2024.
#
# Does NOT touch seeds 42/123/999 -- SEEDS below lists only 7 and 2024.
# Same worker script (multistart_profile_Kss_corrected_pilot.R), same output
# file (multistart_results_corrected_pilot.csv) -- the worker's existing
# resume-check appends new rows and skips anything already recorded, so this
# is additive/crash-safe by the same mechanism as the original driver.
# =============================================================================
set -uo pipefail

REPO="$HOME/pharmacometrics/denosumab-tmdd-qss"
SCRIPT="$REPO/multistart_profile_Kss_corrected_pilot.R"
LOGDIR="$REPO/multistart_logs_corrected_pilot"
mkdir -p "$LOGDIR"

SEEDS="7 2024"
N_POINTS=7
PAUSE_SECONDS=15

cd "$REPO" || { echo "Cannot cd to $REPO"; exit 1; }

total=$((N_POINTS * 2))
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
echo "ALL $total TOP-UP (grid_point, seed) COMBINATIONS PROCESSED."
echo "Results (appended): $REPO/multistart_results_corrected_pilot.csv"
echo "Per-fit logs: $LOGDIR/"
