#!/usr/bin/env bash
# proxy_EWHa.sh - run proxy_EWHa.py in parallel for a list of MAUVE galaxies.
#
# Usage examples:
#   ./proxy_EWHa.sh                   # default galaxy list below
#   ./proxy_EWHa.sh NGC4064 NGC4192   # custom subset
#
# The script checks CANFAR per-galaxy inputs in both:
#   ${ROOT_CANFAR_BASE}/products/v0.6/${GAL}
#   ${ROOT_CANFAR_BASE}/cubes/v3.0/${GAL}
# and falls back to the current working directory for any file missing from
# both CANFAR locations.

set -euo pipefail

export LC_ALL=C
export LANG=C

# ──────────────────────────────────────────────────────────────
# 1.  Configurable variables
# ──────────────────────────────────────────────────────────────
ROOT_CANFAR_BASE="/arc/projects/mauve"
ROOT_LOCAL="$PWD"
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

[[ $# -gt 0 ]] && GALAXIES=("$@")

resolve_first_existing() {
  local candidate
  for candidate in "$@"; do
    [[ -n "$candidate" && -e "$candidate" ]] || continue
    printf '%s\n' "$candidate"
    return 0
  done
  return 1
}

# ──────────────────────────────────────────────────────────────
# 2.  Process function for each galaxy
# ──────────────────────────────────────────────────────────────
process_galaxy() {
  local GAL="$1"
  local LOGFILE="$LOGDIR/${GAL}.log"
  local ROOT_PRODUCTS="${ROOT_CANFAR_BASE}/products/v0.6/${GAL}"
  local ROOT_CUBES="${ROOT_CANFAR_BASE}/cubes/v3.0/${GAL}"

  local BIN_FILE=""
  local GAS_FILE=""
  local CONT_FILE=""
  local REDSHIFT_FILE="${ROOT_LOCAL}/new_redshifts"

  BIN_FILE="$(resolve_first_existing \
    "${ROOT_PRODUCTS}/${GAL}_SPATIAL_BINNING_maps_extended.fits" \
    "${ROOT_CUBES}/${GAL}_SPATIAL_BINNING_maps_extended.fits" \
    "${ROOT_LOCAL}/${GAL}_SPATIAL_BINNING_maps_extended.fits")" || BIN_FILE=""

  GAS_FILE="$(resolve_first_existing \
    "${ROOT_PRODUCTS}/${GAL}_gas_BIN_maps_extended.fits" \
    "${ROOT_CUBES}/${GAL}_gas_BIN_maps_extended.fits" \
    "${ROOT_PRODUCTS}/${GAL}_gas_BIN_maps.fits" \
    "${ROOT_CUBES}/${GAL}_gas_BIN_maps.fits" \
    "${ROOT_LOCAL}/${GAL}_gas_BIN_maps_extended.fits" \
    "${ROOT_LOCAL}/${GAL}_gas_BIN_maps.fits")" || GAS_FILE=""

  CONT_FILE="$(resolve_first_existing \
    "${ROOT_PRODUCTS}/${GAL}_CONTcube.fits" \
    "${ROOT_CUBES}/${GAL}_CONTcube.fits" \
    "${ROOT_LOCAL}/${GAL}_CONTcube.fits")" || CONT_FILE=""

  printf "\n====================  %s  ====================\n" "$GAL"

  start=$(date +%s)
  set +e
  {
    echo "Python executable: $PYTHON_BIN"
    "$PYTHON_BIN" --version
    echo "CANFAR products root: $ROOT_PRODUCTS"
    echo "CANFAR cubes root   : $ROOT_CUBES"
    echo "Local fallback root : $ROOT_LOCAL"
    echo "Resolved bin file   : ${BIN_FILE:-<missing>}"
    echo "Resolved gas file   : ${GAS_FILE:-<missing>}"
    echo "Resolved cont file  : ${CONT_FILE:-<missing>}"
    echo "Redshift file       : $REDSHIFT_FILE"
    echo
  } >"$LOGFILE" 2>&1

  if [[ -z "$BIN_FILE" || -z "$GAS_FILE" || -z "$CONT_FILE" ]]; then
    {
      echo "ERROR: missing required input file(s)."
      [[ -n "$BIN_FILE" ]] || echo "  - missing bin file"
      [[ -n "$GAS_FILE" ]] || echo "  - missing gas file"
      [[ -n "$CONT_FILE" ]] || echo "  - missing continuum cube"
    } >>"$LOGFILE" 2>&1
    status=1
  else
    "$PYTHON_BIN" "$SCRIPT" \
      -g "$GAL" \
      --bin-file "$BIN_FILE" \
      --gas-file "$GAS_FILE" \
      --cont-file "$CONT_FILE" \
      --redshift-file "$REDSHIFT_FILE" \
      >>"$LOGFILE" 2>&1
    status=$?
  fi
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
}

export -f process_galaxy
export -f resolve_first_existing
export ROOT_CANFAR_BASE ROOT_LOCAL PYTHON_BIN SCRIPT LOGDIR

# ──────────────────────────────────────────────────────────────
# 3.  Parallel execution
# ──────────────────────────────────────────────────────────────
all_start=$(date +%s)

printf "Running %d galaxies in parallel using %d cores...\n" "${#GALAXIES[@]}" "$CORES"
printf "Using Python executable: %s\n" "$PYTHON_BIN"

if command -v parallel >/dev/null 2>&1; then
  printf '%s\n' "${GALAXIES[@]}" | parallel -j "$CORES" process_galaxy
else
  printf '%s\n' "${GALAXIES[@]}" | xargs -n 1 -P "$CORES" -I {} bash -c 'process_galaxy "$@"' _ {}
fi

all_end=$(date +%s)
tot=$((all_end - all_start))
printf "\n🏁  proxy_EWHa.sh completed in %dh%02dm%02ds using %d cores\n" \
     $((tot/3600)) $(((tot/60)%60)) $((tot%60)) "$CORES"
