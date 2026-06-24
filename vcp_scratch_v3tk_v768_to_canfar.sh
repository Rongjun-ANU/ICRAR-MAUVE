#!/usr/bin/env bash
set -euo pipefail

# Upload selected nGIST v7.6.8 products directly from Setonix scratch to CANFAR.
#
# Usage:
#   ./vcp_scratch_v3tk_v768_to_canfar.sh [--dry-run] [GALID ...]
#
# With no GALID arguments, upload every galaxy in ALL_GALAXIES. Set JOBS to
# control concurrency; it is capped at five to limit Setonix overlay usage.

start_secs=$SECONDS

SOURCE_BASE=${SOURCE_BASE:-/scratch/pawsey1308/mauve/products/v3tk_v7.6.8}
DEST_BASE=${DEST_BASE:-arc:home/RongjunHuang/ICRAR/further/v3tk_v7.6.8}
CADC_USER=${CADC_USER:-RongjunHuang}
BASE_OVERLAY=${BASE_OVERLAY:-/software/projects/pawsey1308/containers/cadc_overlay.img}
OVERLAY_DIR=${OVERLAY_DIR:-/software/projects/pawsey1308/containers}
MAX_JOBS=5
JOBS=${JOBS:-5}
RUN_ID=$$

ALL_GALAXIES=(
	IC3392
	NGC4064
	NGC4189
	NGC4192
	NGC4216
	NGC4222
	NGC4254
	NGC4293
	NGC4294
	NGC4298
	NGC4302
	NGC4321
	NGC4330
	NGC4351
	NGC4380
	NGC4383
	NGC4388
	NGC4394
	NGC4396
	NGC4405
	NGC4402
	NGC4419
	NGC4424
	NGC4450
	NGC4457
	NGC4501
	NGC4522
	NGC4535
	NGC4548
	NGC4567_8
	NGC4569
	NGC4579
	NGC4580
	NGC4606
	NGC4607
	NGC4654
	NGC4689
	NGC4694
	NGC4698
)

PRODUCT_SUFFIXES=(
	CONFIG
	_sfh_maps.fits
	_gas_bin_maps.fits
	_sfh_weights.fits
	_gas_spaxel_maps.fits
	_spatial_binning_maps.fits
	_kin_maps.fits
	LOGFILE
	_mask.fits
)

usage() {
	cat <<'EOF'
Usage: ./vcp_scratch_v3tk_v768_to_canfar.sh [--dry-run] [GALID ...]

Upload the selected nGIST v7.6.8 products directly from Setonix scratch to
CANFAR. With no GALID arguments, upload the complete built-in galaxy list.

Options:
  -n, --dry-run  Show selected galaxies and source files without contacting CADC
  -h, --help     Show this help

Environment:
  JOBS           Concurrent workers (default 5, maximum 5)
EOF
}

format_runtime() {
	local total=$1
	printf '%02d:%02d:%02d' \
		$((total / 3600)) \
		$(((total % 3600) / 60)) \
		$((total % 60))
}

is_known_galaxy() {
	local candidate=$1
	local galaxy
	for galaxy in "${ALL_GALAXIES[@]}"; do
		if [ "$candidate" = "$galaxy" ]; then
			return 0
		fi
	done
	return 1
}

source_path() {
	local galaxy=$1
	local suffix=$2
	case "$suffix" in
		CONFIG|LOGFILE)
			printf '%s/%s/%s\n' "$SOURCE_BASE" "$galaxy" "$suffix"
			;;
		*)
			printf '%s/%s/%s%s\n' "$SOURCE_BASE" "$galaxy" "$galaxy" "$suffix"
			;;
	esac
}

wait_for_overlay() {
	local overlay=$1
	local tries=0
	while ! flock -n "$overlay" -c true 2>/dev/null; do
		tries=$((tries + 1))
		if [ "$tries" -ge 30 ]; then
			echo "ERROR: overlay still busy after 30 seconds: $overlay" >&2
			return 1
		fi
		sleep 1
	done
}

DRY_RUN=0
REQUESTED_GALAXIES=()
while [ "$#" -gt 0 ]; do
	case "$1" in
		-n|--dry-run)
			DRY_RUN=1
			;;
		-h|--help)
			usage
			exit 0
			;;
		--)
			shift
			while [ "$#" -gt 0 ]; do
				REQUESTED_GALAXIES+=("$1")
				shift
			done
			break
			;;
		-*)
			echo "ERROR: unknown option: $1" >&2
			usage >&2
			exit 2
			;;
		*)
			REQUESTED_GALAXIES+=("$1")
			;;
	esac
	shift
done

if ! [[ "$JOBS" =~ ^[0-9]+$ ]] || [ "$JOBS" -lt 1 ]; then
	echo "ERROR: JOBS must be a positive integer, got: $JOBS" >&2
	exit 2
fi
if [ "$JOBS" -gt "$MAX_JOBS" ]; then
	echo "Capping JOBS from $JOBS to $MAX_JOBS for Setonix overlay safety."
	JOBS=$MAX_JOBS
fi

if [ "${#REQUESTED_GALAXIES[@]}" -gt 0 ]; then
	GALAXIES=("${REQUESTED_GALAXIES[@]}")
else
	GALAXIES=("${ALL_GALAXIES[@]}")
fi

for galaxy in "${GALAXIES[@]}"; do
	if ! is_known_galaxy "$galaxy"; then
		echo "ERROR: unknown galaxy ID: $galaxy" >&2
		echo "Known galaxy IDs: ${ALL_GALAXIES[*]}" >&2
		exit 2
	fi
done

EFFECTIVE_JOBS=$JOBS
if [ "${#GALAXIES[@]}" -lt "$EFFECTIVE_JOBS" ]; then
	EFFECTIVE_JOBS=${#GALAXIES[@]}
fi

if [ "$DRY_RUN" -eq 1 ]; then
	echo "Requested galaxies: ${GALAXIES[*]}"
	echo "Galaxy count: ${#GALAXIES[@]}"
	echo "Effective workers: $EFFECTIVE_JOBS"
	echo "Destination: $DEST_BASE"
	for galaxy in "${GALAXIES[@]}"; do
		echo "$galaxy:"
		for suffix in "${PRODUCT_SUFFIXES[@]}"; do
			path=$(source_path "$galaxy" "$suffix")
			if [ -f "$path" ]; then
				echo "  [present] $path"
			else
				echo "  [missing] $path"
			fi
		done
	done
	exit 0
fi

command -v cadc-get-cert >/dev/null 2>&1 || {
	echo "ERROR: cadc-get-cert is not available in PATH." >&2
	exit 1
}
command -v vmkdir >/dev/null 2>&1 || {
	echo "ERROR: vmkdir is not available in PATH." >&2
	exit 1
}
command -v vcp >/dev/null 2>&1 || {
	echo "ERROR: vcp is not available in PATH." >&2
	exit 1
}
command -v flock >/dev/null 2>&1 || {
	echo "ERROR: flock is not available in PATH." >&2
	exit 1
}
if [ ! -f "$BASE_OVERLAY" ]; then
	echo "ERROR: CADC overlay not found: $BASE_OVERLAY" >&2
	exit 1
fi

worker_dir=$(mktemp -d)
worker_overlays=()

cleanup() {
	local overlay
	for overlay in "${worker_overlays[@]:-}"; do
		rm -f "$overlay"
	done
	rm -rf "$worker_dir"
}
trap cleanup EXIT

stop_workers() {
	echo "Interrupted. Stopping active transfers..." >&2
	pkill -TERM -P $$ 2>/dev/null || true
	wait 2>/dev/null || true
	exit 130
}
trap stop_workers INT TERM

cadc-get-cert -u "$CADC_USER"
wait_for_overlay "$BASE_OVERLAY"

for ((i = 0; i < EFFECTIVE_JOBS; i++)); do
	worker_overlay="$OVERLAY_DIR/cadc_overlay_v3tk_upload_${RUN_ID}_${i}.img"
	cp --reflink=auto "$BASE_OVERLAY" "$worker_overlay"
	wait_for_overlay "$worker_overlay"
	worker_overlays+=("$worker_overlay")
	: > "$worker_dir/part_${i}.txt"
done

for ((i = 0; i < ${#GALAXIES[@]}; i++)); do
	printf '%s\n' "${GALAXIES[$i]}" >> "$worker_dir/part_$((i % EFFECTIVE_JOBS)).txt"
done

upload_galaxy() {
	local galaxy=$1
	local overlay=$2
	local destination="${DEST_BASE}/${galaxy}"
	local suffix path
	local sources=()

	for suffix in "${PRODUCT_SUFFIXES[@]}"; do
		path=$(source_path "$galaxy" "$suffix")
		if [ -f "$path" ]; then
			sources+=("$path")
		else
			echo "WARNING [$galaxy]: missing source, skipping: $path" >&2
		fi
	done

	if [ "${#sources[@]}" -eq 0 ]; then
		echo "ERROR [$galaxy]: no requested products found." >&2
		return 1
	fi

	echo "Uploading $galaxy (${#sources[@]} files)..."
	CADC_OVERLAY="$overlay" vmkdir "$destination" 2>/dev/null || true
	CADC_OVERLAY="$overlay" vcp -v "${sources[@]}" "${destination}/"
	echo "Finished $galaxy"
}

pids=()
for ((i = 0; i < EFFECTIVE_JOBS; i++)); do
	part_file="$worker_dir/part_${i}.txt"
	worker_overlay="${worker_overlays[$i]}"
	(
		worker_status=0
		while IFS= read -r galaxy; do
			if ! upload_galaxy "$galaxy" "$worker_overlay"; then
				worker_status=1
			fi
		done < "$part_file"
		exit "$worker_status"
	) &
	pids+=("$!")
done

status=0
for pid in "${pids[@]}"; do
	if ! wait "$pid"; then
		status=1
	fi
done

runtime_secs=$((SECONDS - start_secs))
echo "Total runtime: $(format_runtime "$runtime_secs") (${runtime_secs}s)"

if [ "$status" -ne 0 ]; then
	echo "One or more galaxy uploads failed." >&2
	exit "$status"
fi

echo "Upload to CANFAR finished."
