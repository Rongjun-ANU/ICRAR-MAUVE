#!/usr/bin/env bash
# Mass.sh – run Mass.py sequentially for a list of MAUVE galaxies
# --------------------------------------------------------------
# Usage examples:
#   ./Mass.sh                 # default galaxy list below
#   ./Mass.sh NGC4064 NGC4192 # custom list
# --------------------------------------------------------------

set -euo pipefail

# ──────────────────────────────────────────────────────────────
# 1.  User-configurable variables
# ──────────────────────────────────────────────────────────────
ROOT="/arc/projects/mauve"      # MAUVE root directory
SCRIPT="Mass.py"                # Python script to call
LOGDIR="mass_logs"              # Per-galaxy logs live here
mkdir -p "$LOGDIR"

GALAXIES=(
  IC3392  
  NGC4064  
  NGC4192  
  NGC4293  
  NGC4298
  NGC4330 
  NGC4383  
  NGC4396  
  NGC4419  
  NGC4457
  NGC4501 
  NGC4522  
  NGC4694  
  NGC4698
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
  status=${PIPESTATUS[0]}               # exit code of python, not tee
  end=$(date +%s)
  dur=$((end - start))
  mins=$((dur / 60)); secs=$((dur % 60))

  if [[ $status -eq 0 ]]; then
    msg="✅  $GAL finished in ${mins}m${secs}s"
  else
    msg="🛑  $GAL failed (exit $status) after ${mins}m${secs}s – see $LOGFILE"
  fi
  echo "$msg" | tee -a "$LOGFILE"      # append to log + echo to screen
done

all_end=$(date +%s)
tot=$((all_end - all_start))
printf "\n🏁  Mass.sh completed in %dh%02dm%02ds\n" \
       $((tot/3600)) $(((tot/60)%60)) $((tot%60))
