#!/usr/bin/env bash
set -euo pipefail

# Upload selected nGIST v7.6.8 products directly from Setonix scratch to CANFAR.
#
# Improvements in this version:
#   - temporary worker overlays are stored on /scratch, not /software
#   - files are uploaded one at a time using a one-file staging directory
#   - each file is retried independently with exponential backoff
#   - one failed file does not prevent the remaining files in that galaxy from running
#   - existing identical remote files are left for vcp to detect and skip
#   - no vls/vmkdir dependency; vcp handles the remote galaxy directory as before
#
# Usage:
#   ./vcp_scratch_v3tk_v768_to_canfar.sh [--dry-run] [normal|7000] [GALID ...]
#
# Examples:
#   ./vcp_scratch_v3tk_v768_to_canfar.sh 7000 NGC4607 NGC4698
#   JOBS=1 ./vcp_scratch_v3tk_v768_to_canfar.sh 7000
#   FILE_RETRIES=8 RETRY_BASE_SLEEP=20 ./vcp_scratch_v3tk_v768_to_canfar.sh normal

start_secs=$SECONDS

SOURCE_NORMAL=${SOURCE_NORMAL:-/scratch/pawsey1308/mauve/products/v3tk_v7.6.8}
SOURCE_7000=${SOURCE_7000:-/scratch/pawsey1308/mauve/products/v3tk_v7.6.8_7000}
DEST_NORMAL=${DEST_NORMAL:-arc:projects/mauve/products/v3tk_v7.6.8}
DEST_7000=${DEST_7000:-arc:projects/mauve/products/v3tk_v7.6.8_7000}
CADC_USER=${CADC_USER:-RongjunHuang}

BASE_OVERLAY=${BASE_OVERLAY:-/software/projects/pawsey1308/containers/cadc_overlay.img}
DEFAULT_OVERLAY_DIR="${MYSCRATCH:-/scratch/pawsey1308/$USER}/cadc_upload_overlays"
OVERLAY_DIR=${OVERLAY_DIR:-$DEFAULT_OVERLAY_DIR}

# Single-stream is the safest default for the Pawsey -> CANFAR link.
# Increase JOBS manually if desired.
MAX_JOBS=3
JOBS=${JOBS:-1}

# Retry policy for each individual file.
FILE_RETRIES=${FILE_RETRIES:-5}
RETRY_BASE_SLEEP=${RETRY_BASE_SLEEP:-30}
RETRY_MAX_SLEEP=${RETRY_MAX_SLEEP:-240}


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
    cat <<'USAGE'
Usage: ./vcp_scratch_v3tk_v768_to_canfar.sh [--dry-run] [normal|7000] [GALID ...]

Upload selected nGIST products directly from Setonix scratch to CANFAR.

Examples:
  ./vcp_scratch_v3tk_v768_to_canfar.sh
      Upload all available galaxies from normal and 7000 runs.
  ./vcp_scratch_v3tk_v768_to_canfar.sh normal
      Upload all available galaxies from the normal run.
  ./vcp_scratch_v3tk_v768_to_canfar.sh 7000
      Upload all available galaxies from the 7000 run.
  ./vcp_scratch_v3tk_v768_to_canfar.sh 7000 NGC4607 NGC4698
      Upload only NGC4607 and NGC4698 from the 7000 run.

Options:
  -n, --dry-run  Show selected galaxies and source files without contacting CADC
  -h, --help     Show this help

Environment:
  JOBS              Concurrent galaxy workers (default 1, maximum 3)
  FILE_RETRIES      Attempts per individual file (default 5)
  RETRY_BASE_SLEEP  Initial retry delay in seconds (default 30)
  RETRY_MAX_SLEEP   Maximum retry delay in seconds (default 240)
  BASE_OVERLAY      Base CADC overlay image
  OVERLAY_DIR       Temporary worker overlay directory (default: $MYSCRATCH/cadc_upload_overlays)
  VCP_CMD           Optional explicit path to vcp
USAGE
}

format_runtime() {
    local total=$1
    printf '%02d:%02d:%02d' \
        $((total / 3600)) \
        $(((total % 3600) / 60)) \
        $((total % 60))
}

retry_sleep_seconds() {
    local failed_attempt=$1
    local delay=$RETRY_BASE_SLEEP
    local i

    # attempt 1 -> base, attempt 2 -> 2*base, etc.
    for ((i = 1; i < failed_attempt; i++)); do
        delay=$((delay * 2))
        if [ "$delay" -ge "$RETRY_MAX_SLEEP" ]; then
            delay=$RETRY_MAX_SLEEP
            break
        fi
    done

    printf '%s\n' "$delay"
}

is_positive_integer() {
    [[ "$1" =~ ^[0-9]+$ ]] && [ "$1" -ge 1 ]
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

if ! is_positive_integer "$JOBS"; then
    echo "ERROR: JOBS must be a positive integer, got: $JOBS" >&2
    exit 2
fi
if [ "$JOBS" -gt "$MAX_JOBS" ]; then
    echo "Capping JOBS from $JOBS to $MAX_JOBS for transfer reliability."
    JOBS=$MAX_JOBS
fi

if ! is_positive_integer "$FILE_RETRIES"; then
    echo "ERROR: FILE_RETRIES must be a positive integer, got: $FILE_RETRIES" >&2
    exit 2
fi
if ! is_positive_integer "$RETRY_BASE_SLEEP"; then
    echo "ERROR: RETRY_BASE_SLEEP must be a positive integer, got: $RETRY_BASE_SLEEP" >&2
    exit 2
fi
if ! is_positive_integer "$RETRY_MAX_SLEEP"; then
    echo "ERROR: RETRY_MAX_SLEEP must be a positive integer, got: $RETRY_MAX_SLEEP" >&2
    exit 2
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
    echo "File retries: $FILE_RETRIES"

    for work_item in "${WORK_ITEMS[@]}"; do
        run=${work_item%%|*}
        galaxy=${work_item#*|}
        source_base=$(run_source_base "$run")
        dest_base=$(run_dest_base "$run")
        echo "$run $galaxy:"
        echo "  Destination: ${dest_base}/${galaxy}/"
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

if ! mkdir -p "$OVERLAY_DIR"; then
    echo "ERROR: could not create worker overlay directory: $OVERLAY_DIR" >&2
    exit 1
fi
if [ ! -w "$OVERLAY_DIR" ]; then
    echo "ERROR: worker overlay directory is not writable: $OVERLAY_DIR" >&2
    exit 1
fi

echo "Base CADC overlay: $BASE_OVERLAY"
echo "Worker overlay directory: $OVERLAY_DIR"
echo "Worker count: $EFFECTIVE_JOBS"
echo "Per-file retries: $FILE_RETRIES"
echo "Retry backoff: ${RETRY_BASE_SLEEP}s -> max ${RETRY_MAX_SLEEP}s"

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

upload_one_file() {
    local run=$1
    local galaxy=$2
    local path=$3
    local dest_base=$4
    local overlay=$5
    local stage_galaxy_dir=$6
    local filename
    local staged_path
    local attempt delay

    filename=$(basename "$path")
    staged_path="${stage_galaxy_dir}/${filename}"

    # Keep exactly one file in the staged galaxy directory.  Because the
    # staging directory lives under source_base, ln is normally a hard link
    # and consumes essentially no additional data space.  Fall back to cp if
    # a hard link is not possible.
    rm -f "$stage_galaxy_dir"/*
    if ! ln "$path" "$staged_path" 2>/dev/null; then
        if ! cp -p "$path" "$staged_path"; then
            echo "ERROR [$run $galaxy]: could not stage $filename" >&2
            return 1
        fi
    fi

    for ((attempt = 1; attempt <= FILE_RETRIES; attempt++)); do
        echo "[$run $galaxy] $filename -- attempt $attempt/$FILE_RETRIES"

        # This deliberately uses the same directory -> remote-parent layout
        # as the original working script.  The staged galaxy directory
        # contains only this one file, so each vcp invocation is independent.
        if CADC_OVERLAY="$overlay" "$VCP_CMD" -v "$stage_galaxy_dir" "${dest_base}/"; then
            echo "[$run $galaxy] OK: $filename"
            rm -f "$staged_path"
            return 0
        fi

        if [ "$attempt" -lt "$FILE_RETRIES" ]; then
            delay=$(retry_sleep_seconds "$attempt")
            echo "WARNING [$run $galaxy]: $filename failed on attempt $attempt/$FILE_RETRIES; retrying in ${delay}s..." >&2
            sleep "$delay"
        fi
    done

    echo "ERROR [$run $galaxy]: giving up on $filename after $FILE_RETRIES attempts." >&2
    rm -f "$staged_path"
    return 1
}

upload_galaxy() {
    local run=$1
    local galaxy=$2
    local overlay=$3
    local source_base dest_base
    local suffix path
    local stage_root stage_galaxy_dir
    local sources=()
    local galaxy_status=0
    local succeeded=0
    local failed=0

    source_base=$(run_source_base "$run")
    dest_base=$(run_dest_base "$run")

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

    echo "Uploading $run $galaxy (${#sources[@]} files) one file at a time..."

    stage_root=$(mktemp -d "${source_base}/.vcp_upload_stage_${RUN_ID}_${run}_${galaxy}.XXXXXX")
    stage_galaxy_dir="${stage_root}/${galaxy}"
    mkdir -p "$stage_galaxy_dir"

    for path in "${sources[@]}"; do
        if upload_one_file "$run" "$galaxy" "$path" "$dest_base" "$overlay" "$stage_galaxy_dir"; then
            succeeded=$((succeeded + 1))
        else
            failed=$((failed + 1))
            galaxy_status=1
        fi
    done

    rm -rf "$stage_root"

    if [ "$galaxy_status" -eq 0 ]; then
        echo "Finished $run $galaxy: ${succeeded}/${#sources[@]} files successful."
    else
        echo "Finished $run $galaxy with failures: $succeeded succeeded, $failed failed." >&2
    fi

    return "$galaxy_status"
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
    echo "One or more files/galaxies failed after all retries." >&2
    exit "$status"
fi

echo "Upload to CANFAR finished."
