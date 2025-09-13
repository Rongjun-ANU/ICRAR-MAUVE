#!/usr/bin/env bash
# Foreground_star.sh – run Foreground_star.py sequentially for a list of MAUVE galaxies.
#
# Usage examples:
#   ./Foreground_star.sh                   # default galaxy list below
#   ./Foreground_star.sh NGC4064 NGC4192   # custom subset
#
# Per-galaxy runtime is appended to each log, and a grand-total runtime
# is printed at the end.

set -euo pipefail

# ──────────────────────────────────────────────────────────────
# 1.  Configurable variables
# ──────────────────────────────────────────────────────────────
ROOT="/arc/projects/mauve"   # MAUVE root path
SCRIPT="Foreground_star.py"              # Python driver
LOGDIR="Foreground_star_logs"            # Log directory
mkdir -p "$LOGDIR"

GALAXIES=(
  NGC4064  
  NGC4298
  NGC4694  
)

[[ $# -gt 0 ]] && GALAXIES=("$@")   # override list from CLI

# ──────────────────────────────────────────────────────────────
# 2.  Main loop
# ──────────────────────────────────────────────────────────────
all_start=$(date +%s)

for GAL in "${GALAXIES[@]}"; do
  printf "\n====================  %s  ====================\n" "$GAL"
  LOGFILE="$LOGDIR/${GAL}.log"

  start=$(date +%s)
  python "$SCRIPT" -g "$GAL" --root "$ROOT" 2>&1 | tee "$LOGFILE"
  status=${PIPESTATUS[0]}          # exit code of python, not tee
  end=$(date +%s)
  dur=$((end - start))
  mins=$((dur / 60)); secs=$((dur % 60))

  if [[ $status -eq 0 ]]; then
    msg="✅  $GAL finished in ${mins}m${secs}s"
  else
    msg="🛑  $GAL failed (exit $status) after ${mins}m${secs}s – see $LOGFILE"
  fi
  echo "$msg" | tee -a "$LOGFILE"   # append to log + print to screen
done

all_end=$(date +%s)
tot=$((all_end - all_start))
printf "\n🏁  Foreground_star.sh completed in %dh%02dm%02ds\n" \
       $((tot/3600)) $(((tot/60)%60)) $((tot%60))
