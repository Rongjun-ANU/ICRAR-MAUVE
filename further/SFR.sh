#!/usr/bin/env bash
# SFR.sh – run SFR+Z.py in parallel for a list of MAUVE galaxies.
#
# Usage examples:
#   ./SFR.sh                   # default galaxy list below
#   ./SFR.sh NGC4064 NGC4192   # custom subset
#
# Per-galaxy runtime is appended to each log, and a grand-total runtime
# is printed at the end.
#
# Changes (2026-06-29):
#   - Product runs are resolved from /arc/projects/mauve/products first, then
#     ./projects/mauve/products, then the current working directory.
#   - Supports normal/7000 run selectors and writes logs under the selected
#     v3tk_v7.6.8 or v3tk_v7.6.8_7000 product folder.
#   - Auto-discovery only queues IC/NGC galaxy folders that contain gas-map
#     inputs for the selected run.
#
# Changes (2026-07-26):
#   - Validate every queued galaxy ID, including explicit command-line targets,
#     so product-side *_logs directories cannot be treated as galaxies.

set -euo pipefail

# Fix locale settings for parallel
export LC_ALL=C
export LANG=C

# ──────────────────────────────────────────────────────────────
# 1.  Configurable variables
# ──────────────────────────────────────────────────────────────
ROOT_LOCAL="$PWD"                        # Local fallback root
SCRIPT="SFR+Z.py"                        # Python driver
LOG_PREFIX="sfr"                         # Per-run logs live under ${PRODUCT_SUBDIR}/${LOG_PREFIX}_logs

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

discover_product_parent() {
  local candidates=(
    "${PRODUCT_PARENT:-}"
    "/arc/projects/mauve/products"
    "$PWD/projects/mauve/products"
    "$PWD"
  )
  local candidate
  for candidate in "${candidates[@]}"; do
    [[ -n "$candidate" && -d "$candidate" ]] && { printf '%s\n' "$candidate"; return; }
  done
  printf '%s\n' "$PWD"
}

run_subdir() {
  case "$1" in
    normal) printf '%s\n' "v3tk_v7.6.8" ;;
    7000) printf '%s\n' "v3tk_v7.6.8_7000" ;;
  esac
}

is_run_selector() {
  case "$1" in
    normal|7000) return 0 ;;
    *) return 1 ;;
  esac
}

is_galaxy_dir_name() {
  [[ "$1" =~ ^(IC[0-9]+|NGC[0-9]+(_[0-9]+)?)$ ]]
}

has_sfr_inputs() {
  local gal_dir="$1"
  local gal="$2"
  [[ -f "$gal_dir/${gal}_gas_BIN_maps.fits" || -f "$gal_dir/${gal}_gas_BIN_maps.fits.gz" || \
     -f "$gal_dir/${gal}_gas_bin_maps.fits" || -f "$gal_dir/${gal}_gas_bin_maps.fits.gz" || \
     -f "$gal_dir/${gal}_BIN_maps.fits" || -f "$gal_dir/${gal}_BIN_maps.fits.gz" || \
     -f "$gal_dir/${gal}_bin_maps.fits" || -f "$gal_dir/${gal}_bin_maps.fits.gz" ]]
}

run_root_for_subdir() {
  local subdir="$1"
  if [[ -d "$PRODUCT_PARENT/$subdir" ]]; then
    printf '%s\n' "$PRODUCT_PARENT"
  elif [[ -d "$PWD/$subdir" ]]; then
    printf '%s\n' "$PWD"
  else
    printf '%s\n' "$PRODUCT_PARENT"
  fi
}

discover_galaxies_for_run() {
  local run_root="$1"
  local subdir="$2"
  local run_dir="$run_root/$subdir"
  local dir gal
  [[ -d "$run_dir" ]] || return 0
  for dir in "$run_dir"/*; do
    [[ -d "$dir" ]] || continue
    gal="$(basename "$dir")"
    is_galaxy_dir_name "$gal" && has_sfr_inputs "$dir" "$gal" && printf '%s\n' "$gal"
  done | sort
}

PRODUCT_PARENT="$(discover_product_parent)"
RUN_LABELS=(normal 7000)
GALAXY_ARGS=()
if [[ $# -gt 0 ]]; then
  if is_run_selector "$1"; then
    RUN_LABELS=("$1")
    shift
    GALAXY_ARGS=("$@")
  else
    GALAXY_ARGS=("$@")
  fi
fi

TASKS=()
for RUN_LABEL in "${RUN_LABELS[@]}"; do
  PRODUCT_SUBDIR="$(run_subdir "$RUN_LABEL")"
  RUN_ROOT="$(run_root_for_subdir "$PRODUCT_SUBDIR")"
  if [[ ${#GALAXY_ARGS[@]} -gt 0 ]]; then
    GALAXIES=("${GALAXY_ARGS[@]}")
  else
    mapfile -t GALAXIES < <(discover_galaxies_for_run "$RUN_ROOT" "$PRODUCT_SUBDIR")
  fi
  for GAL in "${GALAXIES[@]}"; do
    if ! is_galaxy_dir_name "$GAL"; then
      printf 'ERROR: invalid galaxy ID %q; expected IC<digits>, NGC<digits>, or NGC<digits>_<digits>.\n' "$GAL" >&2
      exit 2
    fi
    TASKS+=("${RUN_LABEL}|${RUN_ROOT}|${PRODUCT_SUBDIR}|${GAL}")
  done
done

# ──────────────────────────────────────────────────────────────
# 2.  Process function for each galaxy
# ──────────────────────────────────────────────────────────────
process_task() {
  local TASK="$1"
  local RUN_LABEL RUN_ROOT PRODUCT_SUBDIR GAL
  IFS='|' read -r RUN_LABEL RUN_ROOT PRODUCT_SUBDIR GAL <<<"$TASK"
  local LOGDIR="$RUN_ROOT/$PRODUCT_SUBDIR/${LOG_PREFIX}_logs"
  mkdir -p "$LOGDIR"
  local LOGFILE="$LOGDIR/${GAL}.log"
  
  printf "\n====================  %s / %s  ====================\n" "$RUN_LABEL" "$GAL"
  
  start=$(date +%s)
  set +e
  {
    echo "Python executable: $PYTHON_BIN"
    "$PYTHON_BIN" --version
    echo "PYTHONUNBUFFERED: 1"
    echo "Product root     : $RUN_ROOT"
    echo "Product subdir   : $PRODUCT_SUBDIR"
    echo "Local fallback   : $ROOT_LOCAL"
    echo
  } >"$LOGFILE" 2>&1
  PYTHONUNBUFFERED=1 "$PYTHON_BIN" -u "$SCRIPT" \
    -g "$GAL" \
    --root "$RUN_ROOT" \
    --product-subdir "$PRODUCT_SUBDIR" \
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
export -f process_task
export ROOT_LOCAL PYTHON_BIN SCRIPT LOG_PREFIX

# ──────────────────────────────────────────────────────────────
# 3.  Parallel execution
# ──────────────────────────────────────────────────────────────
all_start=$(date +%s)
run_status=0

printf "Running %d run/galaxy tasks in parallel using %d cores...\n" "${#TASKS[@]}" "$CORES"
printf "Using Python executable: %s\n" "$PYTHON_BIN"

# Use GNU parallel if available, otherwise use xargs
set +e
if command -v parallel >/dev/null 2>&1; then
  printf '%s\n' "${TASKS[@]}" | parallel -j "$CORES" process_task
  run_status=$?
else
  printf '%s\n' "${TASKS[@]}" | xargs -P "$CORES" -I {} bash -c 'process_task "$@"' _ {}
  run_status=$?
fi
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
