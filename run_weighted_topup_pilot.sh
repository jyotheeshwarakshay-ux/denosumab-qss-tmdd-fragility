#!/bin/bash
# =============================================================================
# WEIGHTED SEED TOP-UP (PILOT) DRIVER: 12 fits, mechanically identical
# recipe to the reproducibility dataset's top-up (Experiment B design,
# option (a) -- verified against gate_test_N30_Kss1263_repro.csv's actual
# seed/grid-point pattern, not reconstructed from an inferred criterion).
#
# grid point 1 (Kss=0.282): seeds 1, 2025
# grid point 2 (Kss=0.465): seeds 1, 314, 777, 2025
# grid point 3 (Kss=0.766): seeds 1, 314, 777, 2025
# grid point 5 (Kss=2.082): seeds 1, 2025
#
# Isolated -- writes ONLY to weighted_topup_pilot.csv and
# weighted_topup_pilot_logs/. Does not touch multistart_results_corrected_pilot.csv
# or any other existing results file.
#
# LAUNCH (when authorized): caffeinate -i -s ./run_weighted_topup_pilot.sh > <logfile> 2>&1
# =============================================================================
set -uo pipefail

REPO="$HOME/pharmacometrics/denosumab-tmdd-qss"
SCRIPT="$REPO/weighted_topup_pilot.R"
LOGDIR="$REPO/weighted_topup_pilot_logs"
mkdir -p "$LOGDIR"

cd "$REPO" || { echo "Cannot cd to $REPO"; exit 1; }

# (grid_point_index seed) pairs -- the exact 12-fit recipe.
PAIRS=(
  "1 1" "1 2025"
  "2 1" "2 314" "2 777" "2 2025"
  "3 1" "3 314" "3 777" "3 2025"
  "5 1" "5 2025"
)

total=${#PAIRS[@]}
i=0

for pair in "${PAIRS[@]}"; do
  read -r gp seed <<< "$pair"
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
    echo "  (pausing 15s to let memory settle)"
    sleep 15
  fi
done

echo ""
echo "ALL $total (grid_point, seed) COMBINATIONS PROCESSED."
echo "Results: $REPO/weighted_topup_pilot.csv"
echo "Per-fit logs: $LOGDIR/"
