#!/usr/bin/env bash
# SFR.sh – run SFR.py in parallel for a list of MAUVE galaxies.
#
# Usage examples:
#   ./SFR.sh                   # default galaxy list below
#   ./SFR.sh NGC4064 NGC4192   # custom subset
#
# Per-galaxy runtime is appended to each log, and a grand-total runtime
# is printed at the end.

set -euo pipefail

# Fix locale settings for parallel
export LC_ALL=C
export LANG=C

# ──────────────────────────────────────────────────────────────
# 1.  Configurable variables
# ──────────────────────────────────────────────────────────────
ROOT="/arc/projects/mauve"   # MAUVE root path, for CANFAR
# ROOT="$PWD"                  # MAUVE root path, for local testing
SCRIPT="SFR+Z.py"              # Python driver
LOGDIR="sfr_logs"            # Log directory
mkdir -p "$LOGDIR"

# Detect number of cores
if command -v nproc >/dev/null 2>&1; then
  CORES=$(nproc)                # Linux
elif command -v sysctl >/dev/null 2>&1; then
  CORES=$(sysctl -n hw.ncpu)    # macOS
else
  CORES=4                       # fallback
fi

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
# 2.  Process function for each galaxy
# ──────────────────────────────────────────────────────────────
process_galaxy() {
  local GAL="$1"
  local LOGFILE="$LOGDIR/${GAL}.log"
  
  printf "\n====================  %s  ====================\n" "$GAL"
  
  start=$(date +%s)
  python "$SCRIPT" -g "$GAL" --root "$ROOT" >"$LOGFILE" 2>&1
  status=$?
  end=$(date +%s)
  dur=$((end - start))
  mins=$((dur / 60)); secs=$((dur % 60))

  if [[ $status -eq 0 ]]; then
    msg="✅  $GAL finished in ${mins}m${secs}s"
  else
    msg="🛑  $GAL failed (exit $status) after ${mins}m${secs}s – see $LOGFILE"
  fi
  echo "$msg" | tee -a "$LOGFILE"
}

# Export function and variables for parallel execution
export -f process_galaxy
export ROOT SCRIPT LOGDIR

# ──────────────────────────────────────────────────────────────
# 3.  Parallel execution
# ──────────────────────────────────────────────────────────────
all_start=$(date +%s)

printf "Running %d galaxies in parallel using %d cores...\n" "${#GALAXIES[@]}" "$CORES"

# Use GNU parallel if available, otherwise use xargs
if command -v parallel >/dev/null 2>&1; then
  printf '%s\n' "${GALAXIES[@]}" | parallel -j "$CORES" process_galaxy
else
  printf '%s\n' "${GALAXIES[@]}" | xargs -n 1 -P "$CORES" -I {} bash -c 'process_galaxy "$@"' _ {}
fi

all_end=$(date +%s)
tot=$((all_end - all_start))
printf "\n🏁  SFR.sh completed in %dh%02dm%02ds using %d cores\n" \
     $((tot/3600)) $(((tot/60)%60)) $((tot%60)) "$CORES"
