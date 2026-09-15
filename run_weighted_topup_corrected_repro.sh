#!/bin/bash
# =============================================================================
# WEIGHTED SEED TOP-UP DRIVER: 12 fits, fresh seeds ONLY (1, 2025, 314, 777 --
# none of the existing 42/123/999/7/2024), at 4 under-resolved grid points,
# weighted by how under-resolved each was in the repro profile's hit rates:
#
#   grid_point=2 (Kss=0.4646317, hit rate 1/5): seeds 1, 2025, 314, 777
#   grid_point=3 (Kss=0.7660482, hit rate 1/5): seeds 1, 2025, 314, 777
#   grid_point=1 (Kss=0.2818134, hit rate 2/5): seeds 1, 2025
#   grid_point=5 (Kss=2.0823350, hit rate 3/5): seeds 1, 2025
#
# Explicit pair list, not a nested loop over one seed set, since the seed
# count differs per grid point.
#
# APPENDS to gate_test_N30_Kss1263_repro.csv -- the file that already holds the
# full 35-row repro profile -- so the 12 new fits combine into a topped-up
# profile. Logs to weighted_topup_corrected_repro_logs/. Does NOT touch
# multistart_results_corrected_pilot.csv, multistart_results_corrected_full.csv,
# convergence_test_N80_Kss1263.csv, or multistart_results_corrected_repro.csv
# (an unrelated, unused file from an earlier, different worker).
#
# LAUNCH: caffeinate -i -s ./run_weighted_topup_corrected_repro.sh > <logfile> 2>&1
# =============================================================================
set -uo pipefail

REPO="$HOME/pharmacometrics/denosumab-tmdd-qss"
SCRIPT="$REPO/weighted_topup_corrected_repro.R"
LOGDIR="$REPO/weighted_topup_corrected_repro_logs"
mkdir -p "$LOGDIR"

PAUSE_SECONDS=15

# Explicit (grid_point, seed) pairs -- weighted allocation per grid point.
PAIRS=(
  "2 1" "2 2025" "2 314" "2 777"
  "3 1" "3 2025" "3 314" "3 777"
  "1 1" "1 2025"
  "5 1" "5 2025"
)

cd "$REPO" || { echo "Cannot cd to $REPO"; exit 1; }

total=${#PAIRS[@]}
i=0

for pair in "${PAIRS[@]}"; do
  set -- $pair
  gp=$1
  seed=$2
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

echo ""
echo "ALL $total WEIGHTED TOP-UP (grid_point, seed) COMBINATIONS PROCESSED."
echo "Results (appended): $REPO/gate_test_N30_Kss1263_repro.csv"
echo "Per-fit logs: $LOGDIR/"
