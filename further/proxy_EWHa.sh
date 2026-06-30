#!/usr/bin/env bash
# proxy_EWHa.sh - run proxy_EWHa.py in parallel for a list of MAUVE galaxies.
#
# Usage examples:
#   ./proxy_EWHa.sh                   # default galaxy list below
#   ./proxy_EWHa.sh NGC4064 NGC4192   # custom subset
#   ./proxy_EWHa.sh NGC4254           # one PHANGS-native MAUVE galaxy
#
# Changes (2026-06-13):
#   - PHANGS_CUBE_ROOT can point NGC4254/NGC4321/NGC4535 at a separate
#     PHANGS-native datacube storage directory.
#   - PHANGS-native VOS cubes are staged temporarily with vcp when no local
#     copy exists, then removed after the proxy_EWHa.py run. The default
#     staging root is /scratch.
#
# Changes (2026-06-29):
#   - Product runs are resolved from /arc/projects/mauve/products first, then
#     ./projects/mauve/products, then the current working directory.
#   - Supports normal/7000 run selectors and writes logs under the selected
#     v3tk_v7.6.8 or v3tk_v7.6.8_7000 product folder.
#   - Auto-discovery only queues IC/NGC galaxy folders that contain a continuum
#     cube and gas-map input for the selected run.
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
ROOT_LOCAL="$PWD"
CUBE_ROOT="/arc/projects/mauve/cubes/v3tk"
PHANGS_CUBE_ROOT="${PHANGS_CUBE_ROOT:-}"
PHANGS_NATIVE_VOS_DIR="vos:phangs/RELEASES/PHANGS-MUSE/DR1.0/DATACUBES"
PHANGS_STAGING_ROOT="${PHANGS_STAGING_ROOT:-/scratch}"
SCRIPT="proxy_EWHa.py"
LOG_PREFIX="proxy_ewha"

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
  case "$1" in
    IC[0-9]*|NGC[0-9]*) ;;
    *) return 1 ;;
  esac
  case "$1" in
    *[!ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_]*|*_logs|*_log) return 1 ;;
    *) return 0 ;;
  esac
}

has_proxy_ewha_inputs() {
  local gal_dir="$1"
  local gal="$2"
  local cont_ok=1 gas_ok=1
  [[ -f "$gal_dir/${gal}_cont_cube.fits" || -f "$gal_dir/${gal}_cont_cube.fits.gz" || \
     -f "$gal_dir/${gal}_CONTcube.fits" || -f "$gal_dir/${gal}_CONTcube.fits.gz" ]] && cont_ok=0
  [[ -f "$gal_dir/${gal}_gas_bin_maps_further.fits" || -f "$gal_dir/${gal}_gas_bin_maps_further.fits.gz" || \
     -f "$gal_dir/${gal}_gas_BIN_maps_further.fits" || -f "$gal_dir/${gal}_gas_BIN_maps_further.fits.gz" || \
     -f "$gal_dir/${gal}_gas_bin_maps.fits" || -f "$gal_dir/${gal}_gas_bin_maps.fits.gz" || \
     -f "$gal_dir/${gal}_gas_BIN_maps.fits" || -f "$gal_dir/${gal}_gas_BIN_maps.fits.gz" || \
     -f "$gal_dir/${gal}_bin_maps.fits" || -f "$gal_dir/${gal}_bin_maps.fits.gz" || \
     -f "$gal_dir/${gal}_BIN_maps.fits" || -f "$gal_dir/${gal}_BIN_maps.fits.gz" ]] && gas_ok=0
  (( cont_ok == 0 && gas_ok == 0 ))
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
    is_galaxy_dir_name "$gal" && has_proxy_ewha_inputs "$dir" "$gal" && printf '%s\n' "$gal"
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
    TASKS+=("${RUN_LABEL}|${RUN_ROOT}|${PRODUCT_SUBDIR}|${GAL}")
  done
done

is_phangs_native_galid() {
  case "$1" in
    NGC4254|NGC4321|NGC4535) return 0 ;;
    *) return 1 ;;
  esac
}

phangs_native_filename() {
  printf '%s_PHANGS_DATACUBE_native.fits\n' "$1"
}

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
  local PHANGS_ARGS=()
  local STAGED_PHANGS_CUBE=""
  local status=0

  local REDSHIFT_FILE="${ROOT_LOCAL}/new_redshifts"

  printf "\n====================  %s / %s  ====================\n" "$RUN_LABEL" "$GAL"

  start=$(date +%s)
  set +e
  {
    echo "Python executable: $PYTHON_BIN"
    "$PYTHON_BIN" --version
    echo "Product root        : $RUN_ROOT"
    echo "Product subdir      : $PRODUCT_SUBDIR"
    echo "Cube root           : $CUBE_ROOT"
    echo "PHANGS cube root    : ${PHANGS_CUBE_ROOT:-<none>}"
    echo "PHANGS VOS dir      : $PHANGS_NATIVE_VOS_DIR"
    echo "PHANGS stage root   : $PHANGS_STAGING_ROOT"
    echo "Local fallback root : $ROOT_LOCAL"
    echo "Input resolution    : delegated to $SCRIPT"
    echo "Redshift file       : $REDSHIFT_FILE"
    echo
  } >"$LOGFILE" 2>&1

  if is_phangs_native_galid "$GAL"; then
    PHANGS_FILE="$(phangs_native_filename "$GAL")"
    if [[ -n "$PHANGS_CUBE_ROOT" && -f "$PHANGS_CUBE_ROOT/$PHANGS_FILE" ]]; then
      PHANGS_ARGS+=(--phangs-cube-root "$PHANGS_CUBE_ROOT")
    elif [[ -f "$CUBE_ROOT/$PHANGS_FILE" ]]; then
      PHANGS_ARGS+=(--phangs-cube-root "$CUBE_ROOT")
    else
      if ! command -v vcp >/dev/null 2>&1; then
        echo "ERROR: no local PHANGS cube found for $GAL and vcp is not available." >>"$LOGFILE"
        status=1
      else
        STAGE_DIR="${PHANGS_STAGING_ROOT%/}"
        STAGED_PHANGS_CUBE="$STAGE_DIR/$PHANGS_FILE"
        mkdir -p "$STAGE_DIR"
        echo "Staging PHANGS cube: $PHANGS_NATIVE_VOS_DIR/$PHANGS_FILE -> $STAGED_PHANGS_CUBE" >>"$LOGFILE"
        vcp "$PHANGS_NATIVE_VOS_DIR/$PHANGS_FILE" "$STAGED_PHANGS_CUBE" >>"$LOGFILE" 2>&1
        status=$?
        if [[ $status -eq 0 ]]; then
          PHANGS_ARGS+=(--phangs-cube-root "$STAGE_DIR")
        fi
      fi
    fi
  fi

  if [[ $status -eq 0 ]]; then
    "$PYTHON_BIN" "$SCRIPT" \
      -g "$GAL" \
      --root "$RUN_ROOT" \
      --product-subdir "$PRODUCT_SUBDIR" \
      --fallback-root "$ROOT_LOCAL" \
      --cube-root "$CUBE_ROOT" \
      "${PHANGS_ARGS[@]}" \
      --redshift-file "$REDSHIFT_FILE" \
      >>"$LOGFILE" 2>&1
    status=$?
  fi

  if [[ -n "$STAGED_PHANGS_CUBE" ]]; then
    rm -f "$STAGED_PHANGS_CUBE"
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
  return "$status"
}

export -f process_task
export -f is_phangs_native_galid
export -f phangs_native_filename
export ROOT_LOCAL CUBE_ROOT PHANGS_CUBE_ROOT PHANGS_NATIVE_VOS_DIR PHANGS_STAGING_ROOT PYTHON_BIN SCRIPT LOG_PREFIX

# ──────────────────────────────────────────────────────────────
# 3.  Parallel execution
# ──────────────────────────────────────────────────────────────
all_start=$(date +%s)
run_status=0

printf "Running %d run/galaxy tasks in parallel using %d cores...\n" "${#TASKS[@]}" "$CORES"
printf "Using Python executable: %s\n" "$PYTHON_BIN"

set +e
if command -v parallel >/dev/null 2>&1; then
  printf '%s\n' "${TASKS[@]}" | parallel -j "$CORES" process_task
  run_status=$?
else
  printf '%s\n' "${TASKS[@]}" | xargs -n 1 -P "$CORES" -I {} bash -c 'process_task "$@"' _ {}
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
