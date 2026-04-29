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
ROOT_CANFAR_BASE="/arc/projects/mauve"   # MAUVE base root; Python checks products/v0.6 and cubes/v3.0 under this
ROOT_LOCAL="$PWD"                        # Local fallback root
SCRIPT="Mass.py"                # Python script to call
LOGDIR="mass_logs"              # Per-galaxy logs live here
mkdir -p "$LOGDIR"

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
  set +e
  {
    echo "Python executable: $PYTHON_BIN"
    "$PYTHON_BIN" --version
    echo "PYTHONUNBUFFERED: 1"
    echo "CANFAR base root : $ROOT_CANFAR_BASE"
    echo "Local fallback   : $ROOT_LOCAL"
    echo "MASS_DISABLE_STAT: ${MASS_DISABLE_STAT:-0}"
    echo "MASS_NCPUS       : ${MASS_NCPUS}"
    echo "MASS_ROW_BLOCK_SIZE: ${MASS_ROW_BLOCK_SIZE:-default}"
    echo
  } >"$LOGFILE" 2>&1
  PYTHONUNBUFFERED=1 "$PYTHON_BIN" -u "$SCRIPT" \
    -g "$GAL" \
    --root "$ROOT_CANFAR_BASE" \
    --fallback-root "$ROOT_LOCAL" \
    "${EXTRA_ARGS[@]}" \
    2>&1 | tee -a "$LOGFILE"
  status=${PIPESTATUS[0]}               # exit code of python, not tee
  set -e
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
printf "Using Python executable: %s\n" "$PYTHON_BIN"
printf "\n🏁  Mass.sh completed in %dh%02dm%02ds\n" \
       $((tot/3600)) $(((tot/60)%60)) $((tot%60))
