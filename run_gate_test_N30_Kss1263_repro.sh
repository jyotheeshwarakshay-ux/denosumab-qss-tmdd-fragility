#!/bin/bash
# =============================================================================
# REPRODUCIBILITY GATE DRIVER: 5 fits ONLY, at Kss=1.263 (fixed), on the fresh
# independent repro dataset (seed 2025), same SAEM settings as the pilot
# (nBurn=150/nEm=200). Seeds 42, 123, 999, 7, 2024.
#
# Does NOT run the other 6 Kss grid points. Does NOT touch
# multistart_results_corrected_pilot.csv, multistart_results_corrected_full.csv,
# convergence_test_N80_Kss1263.csv, or multistart_results_corrected_repro.csv --
# writes only to gate_test_N30_Kss1263_repro.csv and
# gate_test_N30_Kss1263_repro_logs/.
#
# LAUNCH: caffeinate -i -s ./run_gate_test_N30_Kss1263_repro.sh > <logfile> 2>&1
# =============================================================================
set -uo pipefail

REPO="$HOME/pharmacometrics/denosumab-tmdd-qss"
SCRIPT="$REPO/gate_test_N30_Kss1263_repro.R"
LOGDIR="$REPO/gate_test_N30_Kss1263_repro_logs"
mkdir -p "$LOGDIR"

SEEDS="42 123 999 7 2024"
PAUSE_SECONDS=15

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
    echo "  (pausing ${PAUSE_SECONDS}s to let memory settle)"
    sleep "$PAUSE_SECONDS"
  fi
done

echo ""
echo "ALL $total SEEDS PROCESSED."
echo "Results: $REPO/gate_test_N30_Kss1263_repro.csv"
echo "Per-fit logs: $LOGDIR/"
