#!/usr/bin/env bash
set -euo pipefail

# Upload selected nGIST v7.6.8 products directly from Setonix scratch to CANFAR.
#
# Usage:
#   ./vcp_scratch_v3tk_v768_to_canfar.sh [--dry-run] [normal|7000] [GALID ...]
#
# With no run selector, upload both normal and 7000 runs. With no GALID
# arguments, upload every available galaxy directory in the selected run(s).
# Set JOBS to control concurrency; it is capped at five because each worker
# needs a copied CADC overlay image on Setonix.

start_secs=$SECONDS

SOURCE_NORMAL=${SOURCE_NORMAL:-/scratch/pawsey1308/mauve/products/v3tk_v7.6.8}
SOURCE_7000=${SOURCE_7000:-/scratch/pawsey1308/mauve/products/v3tk_v7.6.8_7000}
DEST_NORMAL=${DEST_NORMAL:-arc:projects/mauve/products/v3tk_v7.6.8}
DEST_7000=${DEST_7000:-arc:projects/mauve/products/v3tk_v7.6.8_7000}
CADC_USER=${CADC_USER:-RongjunHuang}
BASE_OVERLAY=${BASE_OVERLAY:-/software/projects/pawsey1308/containers/cadc_overlay.img}
OVERLAY_DIR=${OVERLAY_DIR:-/software/projects/pawsey1308/containers}
MAX_JOBS=5
JOBS=${JOBS:-5}
RUN_ID=$$

find_command() {
	local command_name=$1
	local override_var
	override_var=$(printf '%s_CMD' "$command_name" | tr '[:lower:]' '[:upper:]')
	if [ -n "${!override_var:-}" ]; then
		if [ -x "${!override_var}" ]; then
			printf '%s\n' "${!override_var}"
			return 0
		fi
		echo "ERROR: $override_var is set but not executable: ${!override_var}" >&2
		return 1
	fi
	local candidate
	if command -v "$command_name" >/dev/null 2>&1; then
		command -v "$command_name"
		return 0
	fi
	for candidate in \
		"/home/$USER/bin/$command_name" \
		"/work/$USER/bin/$command_name" \
		"/cadcenv/bin/$command_name"; do
		if [ -x "$candidate" ]; then
			printf '%s\n' "$candidate"
			return 0
		fi
	done
	return 1
}

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
	_cont_cube.fits
	LOGFILE
	_mask.fits
)

usage() {
	cat <<'EOF'
Usage: ./vcp_scratch_v3tk_v768_to_canfar.sh [--dry-run] [normal|7000] [GALID ...]

Upload selected nGIST products directly from Setonix scratch to CANFAR.

Examples:
  ./vcp_scratch_v3tk_v768_to_canfar.sh
      Upload all available galaxies from normal and 7000 runs.
  ./vcp_scratch_v3tk_v768_to_canfar.sh normal
      Upload all available galaxies from /scratch/.../v3tk_v7.6.8 only.
  ./vcp_scratch_v3tk_v768_to_canfar.sh 7000
      Upload all available galaxies from /scratch/.../v3tk_v7.6.8_7000 only.
  ./vcp_scratch_v3tk_v768_to_canfar.sh NGC4254
      Upload NGC4254 from both normal and 7000 runs.
  ./vcp_scratch_v3tk_v768_to_canfar.sh 7000 NGC4254
      Upload NGC4254 from the 7000 run only.

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
	local source_base=$1
	local galaxy=$2
	local suffix=$3
	case "$suffix" in
		CONFIG|LOGFILE)
			printf '%s/%s/%s\n' "$source_base" "$galaxy" "$suffix"
			;;
		*)
			printf '%s/%s/%s%s\n' "$source_base" "$galaxy" "$galaxy" "$suffix"
			;;
	esac
}

run_source_base() {
	case "$1" in
		normal) printf '%s\n' "$SOURCE_NORMAL" ;;
		7000) printf '%s\n' "$SOURCE_7000" ;;
	esac
}

run_dest_base() {
	case "$1" in
		normal) printf '%s\n' "$DEST_NORMAL" ;;
		7000) printf '%s\n' "$DEST_7000" ;;
	esac
}

available_galaxies_for_run() {
	local source_base=$1
	local galaxy
	for galaxy in "${ALL_GALAXIES[@]}"; do
		if [ -d "$source_base/$galaxy" ]; then
			printf '%s\n' "$galaxy"
		fi
	done
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
REQUESTED_RUNS=()
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
		normal|7000)
			REQUESTED_RUNS+=("$1")
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
	echo "Capping JOBS from $JOBS to $MAX_JOBS for Setonix overlay quota safety."
	JOBS=$MAX_JOBS
fi

if [ "${#REQUESTED_RUNS[@]}" -gt 0 ]; then
	RUNS=("${REQUESTED_RUNS[@]}")
else
	RUNS=(normal 7000)
fi

WORK_ITEMS=()
for run in "${RUNS[@]}"; do
	source_base=$(run_source_base "$run")
	if [ "${#REQUESTED_GALAXIES[@]}" -gt 0 ]; then
		GALAXIES=("${REQUESTED_GALAXIES[@]}")
	else
		mapfile -t GALAXIES < <(available_galaxies_for_run "$source_base")
	fi

	for galaxy in "${GALAXIES[@]}"; do
		if ! is_known_galaxy "$galaxy"; then
			echo "ERROR: unknown galaxy ID: $galaxy" >&2
			echo "Known galaxy IDs: ${ALL_GALAXIES[*]}" >&2
			exit 2
		fi
		WORK_ITEMS+=("${run}|${galaxy}")
	done
done

if [ "${#WORK_ITEMS[@]}" -eq 0 ]; then
	echo "No available galaxies selected."
	exit 0
fi

EFFECTIVE_JOBS=$JOBS
if [ "${#WORK_ITEMS[@]}" -lt "$EFFECTIVE_JOBS" ]; then
	EFFECTIVE_JOBS=${#WORK_ITEMS[@]}
fi

if [ "$DRY_RUN" -eq 1 ]; then
	echo "Requested runs: ${RUNS[*]}"
	if [ "${#REQUESTED_GALAXIES[@]}" -gt 0 ]; then
		echo "Requested galaxies: ${REQUESTED_GALAXIES[*]}"
	else
		echo "Requested galaxies: all available"
	fi
	echo "Work item count: ${#WORK_ITEMS[@]}"
	echo "Effective workers: $EFFECTIVE_JOBS"
	for work_item in "${WORK_ITEMS[@]}"; do
		run=${work_item%%|*}
		galaxy=${work_item#*|}
		source_base=$(run_source_base "$run")
		dest_base=$(run_dest_base "$run")
		echo "$run $galaxy:"
		echo "  Destination: ${dest_base}/${galaxy}"
		for suffix in "${PRODUCT_SUFFIXES[@]}"; do
			path=$(source_path "$source_base" "$galaxy" "$suffix")
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
VCP_CMD=$(find_command vcp) || {
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
	rm -f "$worker_overlay"
	if ! cp --reflink=auto "$BASE_OVERLAY" "$worker_overlay"; then
		rm -f "$worker_overlay"
		echo "ERROR: failed to copy CADC overlay for worker $i: $worker_overlay" >&2
		exit 1
	fi
	wait_for_overlay "$worker_overlay"
	worker_overlays+=("$worker_overlay")
	: > "$worker_dir/part_${i}.txt"
done

for ((i = 0; i < ${#WORK_ITEMS[@]}; i++)); do
	printf '%s\n' "${WORK_ITEMS[$i]}" >> "$worker_dir/part_$((i % EFFECTIVE_JOBS)).txt"
done

upload_galaxy() {
	local run=$1
	local galaxy=$2
	local overlay=$3
	local source_base
	local dest_base
	source_base=$(run_source_base "$run")
	dest_base=$(run_dest_base "$run")
	local destination="${dest_base}/"
	local suffix path
	local sources=()
	local stage_root=""
	local stage_galaxy_dir=""

	for suffix in "${PRODUCT_SUFFIXES[@]}"; do
		path=$(source_path "$source_base" "$galaxy" "$suffix")
		if [ -f "$path" ]; then
			sources+=("$path")
		else
			echo "WARNING [$run $galaxy]: missing source, skipping: $path" >&2
		fi
	done

	if [ "${#sources[@]}" -eq 0 ]; then
		echo "ERROR [$run $galaxy]: no requested products found." >&2
		return 1
	fi

	echo "Uploading $run $galaxy (${#sources[@]} files)..."
	stage_root=$(mktemp -d "${source_base}/.vcp_upload_stage_${RUN_ID}_${run}_${galaxy}.XXXXXX")
	stage_galaxy_dir="${stage_root}/${galaxy}"
	mkdir -p "$stage_galaxy_dir"

	for path in "${sources[@]}"; do
		ln "$path" "$stage_galaxy_dir/" 2>/dev/null || cp -p "$path" "$stage_galaxy_dir/"
	done

	if ! CADC_OVERLAY="$overlay" "$VCP_CMD" -v "$stage_galaxy_dir" "$destination"; then
		rm -rf "$stage_root"
		echo "ERROR [$run $galaxy]: vcp failed." >&2
		return 1
	fi
	rm -rf "$stage_root"
	echo "Finished $run $galaxy"
}

pids=()
for ((i = 0; i < EFFECTIVE_JOBS; i++)); do
	part_file="$worker_dir/part_${i}.txt"
	worker_overlay="${worker_overlays[$i]}"
	(
		worker_status=0
		while IFS= read -r work_item; do
			run=${work_item%%|*}
			galaxy=${work_item#*|}
			if ! upload_galaxy "$run" "$galaxy" "$worker_overlay"; then
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
