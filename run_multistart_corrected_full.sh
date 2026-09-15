#!/bin/bash
# =============================================================================
# FULL REVALIDATION DRIVER: 7-point Kss profile x 5 seeds = 35 fits, on the
# CORRECTED model, fresh 80-subject dataset.
#
# LAUNCH SAFEGUARDS:
#   - Wrapped in `caffeinate -i -s` by the caller (see bottom of this file) so
#     the Mac cannot sleep for the duration of the run -- fixes the confirmed
#     stall-2 cause (a software-sleep event followed by a 14h Power Nap cycle
#     killed the pilot mid-run). This does NOT protect against the
#     app/harness-level cause behind stall 1, which had no system-level
#     explanation -- that's why the crash-safe resume below still matters.
#   - One fresh Rscript process per fit, fully exits before the next starts.
#   - Each result is appended to multistart_results_corrected_full.csv and each
#     per-fit log saved to multistart_logs_corrected_full/ the MOMENT that fit
#     finishes (see multistart_profile_Kss_corrected_full.R's write.table with
#     append=TRUE) -- a stall mid-run costs at most the one in-flight fit.
#   - RESUMABLE: the worker script checks multistart_results_corrected_full.csv
#     for an existing (grid_point, seed) row before doing any work and skips it
#     if already recorded -- rerunning this driver after a stall picks up
#     exactly where it left off, same mechanism proven across all 3 pilot
#     stalls.
#   - FAILED fits are recorded with status=FAILED (via the worker's own
#     tryCatch), never silently dropped or retried.
# =============================================================================
set -uo pipefail

REPO="$HOME/pharmacometrics/denosumab-tmdd-qss"
SCRIPT="$REPO/multistart_profile_Kss_corrected_full.R"
LOGDIR="$REPO/multistart_logs_corrected_full"
mkdir -p "$LOGDIR"

SEEDS="42 123 999 7 2024"
N_POINTS=7
PAUSE_SECONDS=15

cd "$REPO" || { echo "Cannot cd to $REPO"; exit 1; }

total=$((N_POINTS * 5))
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
echo "Results: $REPO/multistart_results_corrected_full.csv"
echo "Per-fit logs: $LOGDIR/"

# =============================================================================
# NOTE ON caffeinate: this script is meant to be LAUNCHED as:
#
#   caffeinate -i -s ./run_multistart_corrected_full.sh > pilot_full.log 2>&1
#
# rather than invoking caffeinate internally -- caffeinate needs to wrap the
# actual foreground process the shell is waiting on for its assertion to hold
# for the whole run, which means wrapping the invocation at launch time, not
# re-execing itself from inside the script. Not launched here -- build only.
# =============================================================================
