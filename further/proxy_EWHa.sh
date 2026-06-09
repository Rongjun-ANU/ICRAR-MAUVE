#!/usr/bin/env bash
# proxy_EWHa.sh - run proxy_EWHa.py in parallel for a list of MAUVE galaxies.
#
# Usage examples:
#   ./proxy_EWHa.sh                   # default galaxy list below
#   ./proxy_EWHa.sh NGC4064 NGC4192   # custom subset
#   ./proxy_EWHa.sh NGC4254           # one PHANGS-native MAUVE galaxy
#
# The script checks per-galaxy v3tk products under:
#   ${ROOT_PRODUCT_BASE}/v3tk_v7.6.8/${GAL}
# and v3tk datacubes under:
#   ${CUBE_ROOT}
# Prior post-processing outputs are expected in the current working directory.

set -euo pipefail

export LC_ALL=C
export LANG=C

# ──────────────────────────────────────────────────────────────
# 1.  Configurable variables
# ──────────────────────────────────────────────────────────────
ROOT_PRODUCT_BASE="$PWD"
ROOT_LOCAL="$PWD"
CUBE_ROOT="/arc/projects/mauve/cubes/v3tk"
SCRIPT="proxy_EWHa.py"
LOGDIR="proxy_ewha_logs"
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

if command -v nproc >/dev/null 2>&1; then
  CORES=$(nproc)
elif command -v sysctl >/dev/null 2>&1; then
  CORES=$(sysctl -n hw.ncpu)
else
  CORES=4
fi

GALAXIES=(
  "IC3392"
  "NGC4064"
  "NGC4189"
  "NGC4192"
  "NGC4254"
  "NGC4293"
  "NGC4294"
  "NGC4298"
  "NGC4302"
  "NGC4321"
  "NGC4330"
  "NGC4351"
  "NGC4383"
  "NGC4388"
  "NGC4394"
  "NGC4396"
  "NGC4402"
  "NGC4405"
  "NGC4419"
  "NGC4450"
  "NGC4457"
  "NGC4501"
  "NGC4522"
  "NGC4535"
  "NGC4567_8"
  "NGC4580"
  "NGC4606"
  "NGC4607"
  "NGC4694"
  "NGC4698"
)

[[ $# -gt 0 ]] && GALAXIES=("$@")

# ──────────────────────────────────────────────────────────────
# 2.  Process function for each galaxy
# ──────────────────────────────────────────────────────────────
process_galaxy() {
  local GAL="$1"
  local LOGFILE="$LOGDIR/${GAL}.log"

  local REDSHIFT_FILE="${ROOT_LOCAL}/new_redshifts"

  printf "\n====================  %s  ====================\n" "$GAL"

  start=$(date +%s)
  set +e
  {
    echo "Python executable: $PYTHON_BIN"
    "$PYTHON_BIN" --version
    echo "Product root        : $ROOT_PRODUCT_BASE"
    echo "Cube root           : $CUBE_ROOT"
    echo "Local fallback root : $ROOT_LOCAL"
    echo "Input resolution    : delegated to $SCRIPT"
    echo "Redshift file       : $REDSHIFT_FILE"
    echo
  } >"$LOGFILE" 2>&1

  "$PYTHON_BIN" "$SCRIPT" \
    -g "$GAL" \
    --root "$ROOT_PRODUCT_BASE" \
    --fallback-root "$ROOT_LOCAL" \
    --cube-root "$CUBE_ROOT" \
    --redshift-file "$REDSHIFT_FILE" \
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

export -f process_galaxy
export ROOT_PRODUCT_BASE ROOT_LOCAL CUBE_ROOT PYTHON_BIN SCRIPT LOGDIR

# ──────────────────────────────────────────────────────────────
# 3.  Parallel execution
# ──────────────────────────────────────────────────────────────
all_start=$(date +%s)
run_status=0

printf "Running %d galaxies in parallel using %d cores...\n" "${#GALAXIES[@]}" "$CORES"
printf "Using Python executable: %s\n" "$PYTHON_BIN"

set +e
if command -v parallel >/dev/null 2>&1; then
  printf '%s\n' "${GALAXIES[@]}" | parallel -j "$CORES" process_galaxy
  run_status=$?
else
  printf '%s\n' "${GALAXIES[@]}" | xargs -n 1 -P "$CORES" -I {} bash -c 'process_galaxy "$@"' _ {}
  run_status=$?
fi
set -e

all_end=$(date +%s)
tot=$((all_end - all_start))
if [[ $run_status -eq 0 ]]; then
  printf "\n🏁  proxy_EWHa.sh completed in %dh%02dm%02ds using %d cores\n" \
     $((tot/3600)) $(((tot/60)%60)) $((tot%60)) "$CORES"
else
  printf "\n🛑  proxy_EWHa.sh completed with one or more failures in %dh%02dm%02ds using %d cores\n" \
     $((tot/3600)) $(((tot/60)%60)) $((tot%60)) "$CORES" >&2
fi

exit "$run_status"
