#!/usr/bin/env bash
# Extract LOGSFR_SURFACE_DENSITY* HDUs from further gas-map FITS products.
#
# Usage:
#   ./SFR_cutout.sh
#   ./SFR_cutout.sh NGC4254 NGC4321
#   ./SFR_cutout.sh normal NGC4254
#   ./SFR_cutout.sh 7000 NGC4254

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="$SCRIPT_DIR/SFR_cutout.py"

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

if ! "$PYTHON_BIN" -c 'import astropy' >/dev/null 2>&1; then
  echo "ERROR: $PYTHON_BIN cannot import astropy." >&2
  echo "       Activate the science environment or set PYTHON_BIN." >&2
  exit 1
fi

discover_product_parent() {
  local candidate
  for candidate in \
    "${PRODUCT_PARENT:-}" \
    "/arc/projects/mauve/products" \
    "$PWD/projects/mauve/products" \
    "$PWD"
  do
    if [[ -n "$candidate" && -d "$candidate" ]]; then
      printf '%s\n' "$candidate"
      return
    fi
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
  [[ "$1" == "normal" || "$1" == "7000" ]]
}

is_galaxy_id() {
  [[ "$1" =~ ^(IC[0-9]+|NGC[0-9]+(_[0-9]+)?)$ ]]
}

has_cutout_input() {
  local galaxy_dir="$1"
  local galaxy="$2"
  [[ -f "$galaxy_dir/${galaxy}_gas_bin_maps_further.fits" || \
     -f "$galaxy_dir/${galaxy}_gas_bin_maps_further.fits.gz" ]]
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

discover_galaxies() {
  local run_dir="$1"
  local galaxy_dir galaxy
  [[ -d "$run_dir" ]] || return 0
  for galaxy_dir in "$run_dir"/*; do
    [[ -d "$galaxy_dir" ]] || continue
    galaxy="$(basename "$galaxy_dir")"
    if is_galaxy_id "$galaxy" && has_cutout_input "$galaxy_dir" "$galaxy"; then
      printf '%s\n' "$galaxy"
    fi
  done | sort
}

PRODUCT_PARENT="$(discover_product_parent)"
RUN_LABELS=(normal 7000)
GALAXY_ARGS=()

if [[ $# -gt 0 ]]; then
  if is_run_selector "$1"; then
    RUN_LABELS=("$1")
    shift
  fi
  GALAXY_ARGS=("$@")
fi

if [[ ${#GALAXY_ARGS[@]} -eq 0 ]]; then
  GALAXY_ARGS=(NGC4254 NGC4535)
fi

task_count=0
for run_label in "${RUN_LABELS[@]}"; do
  product_subdir="$(run_subdir "$run_label")"
  run_root="$(run_root_for_subdir "$product_subdir")"

  if [[ ${#GALAXY_ARGS[@]} -gt 0 ]]; then
    galaxies=("${GALAXY_ARGS[@]}")
  else
    galaxies=()
    while IFS= read -r galaxy; do
      [[ -n "$galaxy" ]] && galaxies+=("$galaxy")
    done < <(discover_galaxies "$run_root/$product_subdir")
  fi

  for galaxy in "${galaxies[@]}"; do
    if ! is_galaxy_id "$galaxy"; then
      printf 'ERROR: invalid galaxy ID %q.\n' "$galaxy" >&2
      exit 2
    fi
    printf '\n[%s] Extracting %s\n' "$run_label" "$galaxy"
    "$PYTHON_BIN" "$SCRIPT" \
      --galaxy "$galaxy" \
      --root "$run_root" \
      --product-subdir "$product_subdir"
    task_count=$((task_count + 1))
  done
done

if [[ $task_count -eq 0 ]]; then
  echo "No matching further gas-map FITS inputs were found." >&2
  exit 1
fi

printf '\nSFR cutout complete for %d run/galaxy task(s).\n' "$task_count"
