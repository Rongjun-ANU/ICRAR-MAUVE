#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# Run the full image pipeline:
#  1) v3tk_observed_VRI_image.py  -> XXX_observed_VRI.png/.pdf
#  2) v3tk_get_legacy.py        -> XXX_legacy_reprojected.jpg
#  3) v3tk_combined_VRI_image.py  -> XXX_combined_VRI.png

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

timestamp="$(date +%Y%m%d_%H%M%S)"
logfile="v3tk_VRI_image_${timestamp}.log"

ncpu() {
	# Most portable first: Python.
	if command -v python >/dev/null 2>&1; then
		python -c 'import os; print(os.cpu_count() or 1)' 2>/dev/null && return 0
	fi
	# POSIX-y.
	if command -v getconf >/dev/null 2>&1; then
		getconf _NPROCESSORS_ONLN 2>/dev/null && return 0
	fi
	# GNU coreutils.
	if command -v nproc >/dev/null 2>&1; then
		nproc 2>/dev/null && return 0
	fi
	# macOS/BSD (may be missing or non-functional on some systems).
	if command -v sysctl >/dev/null 2>&1; then
		sysctl -n hw.logicalcpu 2>/dev/null && return 0
		sysctl -n hw.ncpu 2>/dev/null && return 0
	fi
	echo 1
}

NCPU="$(ncpu)"
NCPU="${NCPU:-1}"

# Legacy Survey rate-limits; throttle concurrent downloads.
# You can override: MAX_DL=1 ./v3tk_VRI_image.sh
MAX_DL="${MAX_DL:-2}"
if [[ "$MAX_DL" -lt 1 ]]; then MAX_DL=1; fi

start_epoch="$(date +%s)"

{
	echo "[$(date)] Starting pipeline in: $ROOT_DIR"
	echo "[$(date)] NCPU=$NCPU"
	echo "[$(date)] MAX_DL=$MAX_DL"
	echo "[$(date)] Log: $logfile"
	echo

	# Prefer conda activation if available; fall back to conda run.
	if command -v conda >/dev/null 2>&1; then
		# shellcheck disable=SC1090
		source "$(conda info --base)/etc/profile.d/conda.sh" || true
		conda activate ICRAR 2>/dev/null || true
	fi

	# Always prefer executing in the ICRAR conda env unless it's already active.
	PY_RUN=(python)
	if [[ "${CONDA_DEFAULT_ENV:-}" != "ICRAR" ]]; then
		if command -v conda >/dev/null 2>&1; then
			PY_RUN=(conda run -n ICRAR python)
		fi
	fi

	echo "[$(date)] Step 1/3: observed R images"
	"${PY_RUN[@]}" v3tk_observed_VRI_image.py --input-dir . --overwrite --workers "$NCPU"
	echo

	echo "[$(date)] Step 2/3: legacy download + reproject"
	"${PY_RUN[@]}" v3tk_get_legacy.py --overwrite --workers "$NCPU" --max-concurrent-downloads "$MAX_DL"
	echo

	echo "[$(date)] Step 3/3: combine observed on legacy"
	"${PY_RUN[@]}" v3tk_combined_VRI_image.py --overwrite --workers "$NCPU"
	echo

	end_epoch="$(date +%s)"
	elapsed="$((end_epoch - start_epoch))"
	echo "[$(date)] DONE in ${elapsed}s"
} 2>&1 | tee "$logfile"

