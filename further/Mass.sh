#!/usr/bin/env bash
# Mass.sh – run Mass.py sequentially for a list of MAUVE galaxies
# --------------------------------------------------------------
# Usage examples:
#   ./Mass.sh                 # default galaxy list below
#   ./Mass.sh NGC4064 NGC4192 # custom list
#   ./Mass.sh NGC4254         # one PHANGS-native MAUVE galaxy
# --------------------------------------------------------------
# Changes (2026-06-13):
#   - PHANGS_CUBE_ROOT can point NGC4254/NGC4321/NGC4535 at a separate
#     PHANGS-native datacube storage directory.
#   - PHANGS-native VOS cubes are staged temporarily with vcp when no local
#     copy exists, then removed after the Mass.py run. The default staging
#     root is /scratch.
#
# Changes (2026-06-29):
#   - Product runs are resolved from /arc/projects/mauve/products first, then
#     ./projects/mauve/products, then the current working directory.
#   - Supports normal/7000 run selectors and writes logs under the selected
#     v3tk_v7.6.8 or v3tk_v7.6.8_7000 product folder.
#   - Auto-discovery only queues IC/NGC galaxy folders that contain stellar
#     population inputs for the selected run.
#
# Changes (2026-07-26):
#   - Validate every queued galaxy ID, including explicit command-line targets,
#     so product-side *_logs directories cannot be treated as galaxies.

set -euo pipefail

# ──────────────────────────────────────────────────────────────
# 1.  User-configurable variables
# ──────────────────────────────────────────────────────────────
ROOT_LOCAL="$PWD"                        # Local fallback root
CUBE_ROOT="/arc/projects/mauve/cubes/v3tk"
PHANGS_CUBE_ROOT="${PHANGS_CUBE_ROOT:-}"
PHANGS_NATIVE_VOS_DIR="vos:phangs/RELEASES/PHANGS-MUSE/DR1.0/DATACUBES"
PHANGS_STAGING_ROOT="${PHANGS_STAGING_ROOT:-/scratch}"
SCRIPT="Mass.py"                # Python script to call
LOG_PREFIX="mass"               # Per-run logs live under ${PRODUCT_SUBDIR}/${LOG_PREFIX}_logs

EXTRA_ARGS=()
if [[ "${MASS_DISABLE_STAT:-0}" == "1" ]]; then
  EXTRA_ARGS+=(--disable-stat-propagation)
fi

if [[ -z "${MASS_NCPUS:-}" ]]; then
  if command -v nproc >/dev/null 2>&1; then
    MASS_NCPUS="$(nproc)"
  elif command -v getconf >/dev/null 2>&1; then
    MASS_NCPUS="$(getconf _NPROCESSORS_ONLN)"
  else
    MASS_NCPUS="0"
  fi
fi
if [[ "${MASS_NCPUS}" =~ ^[0-9]+$ ]] && (( MASS_NCPUS > 0 )); then
  EXTRA_ARGS+=(--ncpus "$MASS_NCPUS")
fi

if [[ -n "${MASS_ROW_BLOCK_SIZE:-}" ]]; then
  EXTRA_ARGS+=(--row-block-size "$MASS_ROW_BLOCK_SIZE")
fi

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

if ! "$PYTHON_BIN" -c 'import numpy, astropy, scipy, matplotlib, speclite, ppxf' >/dev/null 2>&1; then
  echo "ERROR: $PYTHON_BIN is missing one or more required Python packages for Mass.py." >&2
  echo "       Activate the science environment or set PYTHON_BIN to the correct interpreter." >&2
  exit 1
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

has_mass_inputs() {
  local gal_dir="$1"
  local gal="$2"
  local bin_ok=1 sfh_ok=1 weights_ok=1
  [[ -f "$gal_dir/${gal}_SPATIAL_BINNING_maps.fits" || -f "$gal_dir/${gal}_SPATIAL_BINNING_maps.fits.gz" || \
     -f "$gal_dir/${gal}_spatial_binning_maps.fits" || -f "$gal_dir/${gal}_spatial_binning_maps.fits.gz" ]] && bin_ok=0
  [[ -f "$gal_dir/${gal}_SFH_maps.fits" || -f "$gal_dir/${gal}_SFH_maps.fits.gz" || \
     -f "$gal_dir/${gal}_sfh_maps.fits" || -f "$gal_dir/${gal}_sfh_maps.fits.gz" ]] && sfh_ok=0
  [[ -f "$gal_dir/${gal}_sfh-weights.fits" || -f "$gal_dir/${gal}_sfh-weights.fits.gz" || \
     -f "$gal_dir/${gal}_sfh_weights.fits" || -f "$gal_dir/${gal}_sfh_weights.fits.gz" ]] && weights_ok=0
  (( bin_ok == 0 && sfh_ok == 0 && weights_ok == 0 ))
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
    is_galaxy_dir_name "$gal" && has_mass_inputs "$dir" "$gal" && printf '%s\n' "$gal"
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
# 2.  Main loop
# ──────────────────────────────────────────────────────────────
all_start=$(date +%s)
run_status=0

for TASK in "${TASKS[@]}"; do
  IFS='|' read -r RUN_LABEL RUN_ROOT PRODUCT_SUBDIR GAL <<<"$TASK"
  LOGDIR="$RUN_ROOT/$PRODUCT_SUBDIR/${LOG_PREFIX}_logs"
  mkdir -p "$LOGDIR"
  printf "\n====================  %s / %s  ====================\n" "$RUN_LABEL" "$GAL"
  LOGFILE="$LOGDIR/${GAL}.log"
  PHANGS_ARGS=()
  STAGED_PHANGS_CUBE=""
  status=0

  start=$(date +%s)
  set +e
  {
    echo "Python executable: $PYTHON_BIN"
    "$PYTHON_BIN" --version
    echo "PYTHONUNBUFFERED: 1"
    echo "Product root     : $RUN_ROOT"
    echo "Product subdir   : $PRODUCT_SUBDIR"
    echo "Cube root        : $CUBE_ROOT"
    echo "PHANGS cube root : ${PHANGS_CUBE_ROOT:-<none>}"
    echo "PHANGS VOS dir   : $PHANGS_NATIVE_VOS_DIR"
    echo "PHANGS stage root: $PHANGS_STAGING_ROOT"
    echo "Local fallback   : $ROOT_LOCAL"
    echo "MASS_DISABLE_STAT: ${MASS_DISABLE_STAT:-0}"
    echo "MASS_NCPUS       : ${MASS_NCPUS}"
    echo "MASS_ROW_BLOCK_SIZE: ${MASS_ROW_BLOCK_SIZE:-default}"
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
        echo "ERROR: no local PHANGS cube found for $GAL and vcp is not available." | tee -a "$LOGFILE" >&2
        status=1
      else
        STAGE_DIR="${PHANGS_STAGING_ROOT%/}"
        STAGED_PHANGS_CUBE="$STAGE_DIR/$PHANGS_FILE"
        mkdir -p "$STAGE_DIR"
        echo "Staging PHANGS cube: $PHANGS_NATIVE_VOS_DIR/$PHANGS_FILE -> $STAGED_PHANGS_CUBE" | tee -a "$LOGFILE"
        vcp "$PHANGS_NATIVE_VOS_DIR/$PHANGS_FILE" "$STAGED_PHANGS_CUBE" 2>&1 | tee -a "$LOGFILE"
        stage_status=${PIPESTATUS[0]}
        if [[ $stage_status -eq 0 ]]; then
          PHANGS_ARGS+=(--phangs-cube-root "$STAGE_DIR")
        else
          status=$stage_status
        fi
      fi
    fi
  fi

  if [[ $status -eq 0 ]]; then
    PYTHONUNBUFFERED=1 "$PYTHON_BIN" -u "$SCRIPT" \
      -g "$GAL" \
      --root "$RUN_ROOT" \
      --product-subdir "$PRODUCT_SUBDIR" \
      --fallback-root "$ROOT_LOCAL" \
      --cube-root "$CUBE_ROOT" \
      "${PHANGS_ARGS[@]}" \
      "${EXTRA_ARGS[@]}" \
      2>&1 | tee -a "$LOGFILE"
    status=${PIPESTATUS[0]}               # exit code of python, not tee
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
    run_status=1
    msg="🛑  $GAL failed (exit $status) after ${mins}m${secs}s – see $LOGFILE"
  fi
  echo "$msg" | tee -a "$LOGFILE"      # append to log + echo to screen
done

all_end=$(date +%s)
tot=$((all_end - all_start))
printf "Using Python executable: %s\n" "$PYTHON_BIN"
if [[ $run_status -eq 0 ]]; then
  printf "\n🏁  Mass.sh completed in %dh%02dm%02ds\n" \
       $((tot/3600)) $(((tot/60)%60)) $((tot%60))
else
  printf "\n🛑  Mass.sh completed with one or more failures in %dh%02dm%02ds\n" \
       $((tot/3600)) $(((tot/60)%60)) $((tot%60)) >&2
fi

exit "$run_status"
