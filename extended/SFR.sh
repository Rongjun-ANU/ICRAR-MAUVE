#!/usr/bin/env bash
# SFR.sh – run SFR+Z.py in parallel for a list of MAUVE galaxies.
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
ROOT_CANFAR_BASE="/arc/projects/mauve"   # MAUVE base root; Python checks products/v0.6 and cubes/v3.0 under this
ROOT_LOCAL="$PWD"                        # Local fallback root
SCRIPT="SFR+Z.py"                        # Python driver
LOGDIR="sfr_logs"                        # Log directory
mkdir -p "$LOGDIR"

if [[ -n "${PYTHON_BIN:-}" ]]; then
  PYTHON_BIN="$(command -v "$PYTHON_BIN" 2>/dev/null || printf '%s' "$PYTHON_BIN")"
elif [[ -n "${CONDA_PREFIX:-}" && -x "${CONDA_PREFIX}/bin/python" ]]; then
  PYTHON_BIN="${CONDA_PREFIX}/bin/python"
elif command -v python >/dev/null 2>&1; then
  PYTHON_BIN="$(command -v python)"
elif command -v python3 >/dev/null 2>&1; then
  PYTHON_BIN="$(command -v python3)"
else
  echo "ERROR: could not find a usable Python executable." >&2
  exit 1
fi

if ! "$PYTHON_BIN" -c 'import numpy, astropy' >/dev/null 2>&1; then
  echo "ERROR: $PYTHON_BIN is missing one or more required Python packages for SFR+Z.py." >&2
  echo "       Activate the science environment or set PYTHON_BIN to the correct interpreter." >&2
  exit 1
fi

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
  set +e
  {
    echo "Python executable: $PYTHON_BIN"
    "$PYTHON_BIN" --version
    echo "PYTHONUNBUFFERED: 1"
    echo "CANFAR base root : $ROOT_CANFAR_BASE"
    echo "Local fallback   : $ROOT_LOCAL"
    echo
  } >"$LOGFILE" 2>&1
  PYTHONUNBUFFERED=1 "$PYTHON_BIN" -u "$SCRIPT" \
    -g "$GAL" \
    --root "$ROOT_CANFAR_BASE" \
    --fallback-root "$ROOT_LOCAL" \
    >>"$LOGFILE" 2>&1
  status=$?
  set -e
  end=$(date +%s)
  dur=$((end - start))
  mins=$((dur / 60)); secs=$((dur % 60))

  if [[ $status -eq 0 ]]; then
    msg="✅  $GAL finished in ${mins}m${secs}s"
  else
    msg="🛑  $GAL failed (exit $status) after ${mins}m${secs}s – see $LOGFILE"
  fi
  echo "$msg" | tee -a "$LOGFILE"
  return "$status"
}

# Export function and variables for parallel execution
export -f process_galaxy
export ROOT_CANFAR_BASE ROOT_LOCAL PYTHON_BIN SCRIPT LOGDIR

# ──────────────────────────────────────────────────────────────
# 3.  Parallel execution
# ──────────────────────────────────────────────────────────────
all_start=$(date +%s)

printf "Running %d galaxies in parallel using %d cores...\n" "${#GALAXIES[@]}" "$CORES"
printf "Using Python executable: %s\n" "$PYTHON_BIN"

# Use GNU parallel if available, otherwise use xargs
set +e
if command -v parallel >/dev/null 2>&1; then
  printf '%s\n' "${GALAXIES[@]}" | parallel -j "$CORES" process_galaxy
else
  printf '%s\n' "${GALAXIES[@]}" | xargs -P "$CORES" -I {} bash -c 'process_galaxy "$@"' _ {}
fi
run_status=$?
set -e

all_end=$(date +%s)
tot=$((all_end - all_start))
if [[ $run_status -eq 0 ]]; then
  printf "\n🏁  SFR.sh completed in %dh%02dm%02ds using %d cores\n" \
       $((tot/3600)) $(((tot/60)%60)) $((tot%60)) "$CORES"
else
  printf "\n🛑  SFR.sh completed with one or more failures in %dh%02dm%02ds using %d cores\n" \
       $((tot/3600)) $(((tot/60)%60)) $((tot%60)) "$CORES" >&2
fi

exit "$run_status"
