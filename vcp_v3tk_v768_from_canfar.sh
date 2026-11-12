#!/bin/zsh

# Download v3tk_v7.6.8 products from CANFAR
# cadc-get-cert -u RongjunHuang
#
# Usage:
#   ./vcp_v3tk_v768_from_canfar.sh normal            # all galaxies, normal run
#   ./vcp_v3tk_v768_from_canfar.sh 7000              # all galaxies, 7000 run
#   ./vcp_v3tk_v768_from_canfar.sh normal NGC4321    # one galaxy, normal run
#   ./vcp_v3tk_v768_from_canfar.sh 7000 NGC4321      # one galaxy, 7000 run
#   ./vcp_v3tk_v768_from_canfar.sh NGC4321           # one galaxy, both runs
# With no arguments, all listed galaxies from both runs are downloaded.
# TRANSFER_JOBS sets the shared file-worker limit (default 40; 0 uses all slots).
# Files are verified against VOSpace MD5 or a fresh SHA-256 computed on CANFAR.
# The fallback prefers vcp's Python. CHECKSUM_PYTHON overrides interpreter selection.

LOCAL_BASE="/Users/Igniz/Desktop/ICRAR/further"
REMOTE_PRODUCT_ROOT="arc:projects/mauve/products"
ALL_RUNS=("v3tk_v7.6.8" "v3tk_v7.6.8_7000")

ALL_GALAXIES=(
    "IC3392" "NGC4064" "NGC4189" "NGC4192" "NGC4216" "NGC4222"
    "NGC4254" "NGC4293" "NGC4294" "NGC4298" "NGC4302" "NGC4321"
    "NGC4330" "NGC4351" "NGC4380" "NGC4383" "NGC4388" "NGC4394"
    "NGC4396" "NGC4405" "NGC4402" "NGC4419" "NGC4424" "NGC4450"
    "NGC4457" "NGC4501" "NGC4522" "NGC4535" "NGC4548" "NGC4567_8"
    "NGC4569" "NGC4579" "NGC4580" "NGC4606" "NGC4607" "NGC4654"
    "NGC4689" "NGC4694" "NGC4698"
)

typeset -aU RUNS GALAXIES
add_run() {
    case "$1" in
        7000|v3tk_v768_7000|v3tk_v7.6.8_7000)
            RUNS+=("v3tk_v7.6.8_7000")
            ;;
        normal|768|v768|v3tk_v768|v3tk_v7.6.8)
            RUNS+=("v3tk_v7.6.8")
            ;;
        *)
            RUNS+=("$1")
            ;;
    esac
}

while [ "$#" -gt 0 ]; do
    case "$1" in
        --run)
            if [ "$#" -lt 2 ]; then
                echo "ERROR: --run requires a run name" >&2
                exit 1
            fi
            shift
            add_run "$1"
            ;;
        --run=*)
            add_run "${1#--run=}"
            ;;
        7000|normal|768|v768|v3tk_v768|v3tk_v7.6.8|v3tk_v768_7000|v3tk_v7.6.8_7000)
            add_run "$1"
            ;;
        *)
            GALAXIES+=("$1")
            ;;
    esac
    shift
done

if [ "${#RUNS[@]}" -eq 0 ]; then
    RUNS=("${ALL_RUNS[@]}")
fi

if [ "${#GALAXIES[@]}" -eq 0 ]; then
    GALAXIES=("${ALL_GALAXIES[@]}")
fi

TRANSFER_JOBS=${TRANSFER_JOBS:-40}
if [[ "$TRANSFER_JOBS" != <-> ]]; then
    echo "ERROR: TRANSFER_JOBS must be a non-negative integer" >&2
    exit 1
fi

queue_downloads() {
    local run galaxy product log_dir
    for run in "${RUNS[@]}"; do
        for galaxy in "${GALAXIES[@]}"; do
            for product in CONFIG LOGFILE \
                "${galaxy}_mask.fits" "${galaxy}_sfh_maps.fits" \
                "${galaxy}_sfh_weights.fits" "${galaxy}_kin_maps.fits" \
                "${galaxy}_gas_bin_maps_further.fits" \
                "${galaxy}_proxy_EW_maps_further.fits" \
                "${galaxy}_spatial_binning_maps_further.fits"; do
                printf '%s\0' "$REMOTE_PRODUCT_ROOT/$run/$galaxy/$product" "$LOCAL_BASE/$run/$galaxy/$product" "-"
            done
            for log_dir in mass_logs sfr_logs proxy_ewha_logs; do
                printf '%s\0' "$REMOTE_PRODUCT_ROOT/$run/$log_dir/$galaxy.log" "$LOCAL_BASE/$run/$log_dir/$galaxy.log" "-"
            done
        done
    done
}

run_queue() {
    xargs -0 -n 3 -P "$TRANSFER_JOBS" bash -c '
    if [ "$#" -ne 3 ]; then
        echo "ERROR: incomplete download queue entry" >&2
        exit 1
    fi
    remote_file=$1
    local_file=$2
    source_hash=$3
    algorithm=SHA-256

    local_hash() {
        local value
        if [ "$algorithm" = SHA-256 ]; then
            if command -v sha256sum >/dev/null 2>&1; then
                value=$(sha256sum < "$1") || return 1
            else
                value=$(shasum -a 256 < "$1") || return 1
            fi
            printf "%s\n" "${value%% *}"
        elif command -v md5 >/dev/null 2>&1; then
            md5 -q -- "$1"
        else
            value=$(md5sum < "$1") || return 1
            printf "%s\n" "${value%% *}"
        fi
    }

    if [ "$source_hash" = - ]; then
        algorithm=MD5
        source_hash=$(vtag "$remote_file" MD5 2>&1)
        tag_status=$?
        if [ "$tag_status" -eq 2 ]; then
            echo "Skipping $remote_file: source file is missing"
            exit 0
        fi
        if [ "$tag_status" -ne 0 ]; then
            printf "ERROR: could not read CANFAR MD5 for %s: %s\n" "$remote_file" "$source_hash" >&2
            exit 1
        fi
        source_hash=$(printf "%s" "$source_hash" | tr -d "\047\042[:space:]" | tr "A-F" "a-f")
        if [[ ! "$source_hash" =~ ^[[:xdigit:]]{32}$ ]]; then
            record=$(mktemp "$checksum_dir/pending.XXXXXXXX") || exit 1
            printf "%s\0" "$remote_file" "$local_file" > "$record" || exit 1
            exit 0
        fi
    elif [[ ! "$source_hash" =~ ^[0-9a-f]{64}$ ]]; then
        echo "ERROR: invalid SHA-256 in the download queue for $remote_file" >&2
        exit 1
    fi
    if [ -f "$local_file" ] && target_hash=$(local_hash "$local_file") && [ "$target_hash" = "$source_hash" ]; then
        echo "Verified existing file, skipping: $local_file ($algorithm)"
        exit 0
    fi

    target_dir=${local_file%/*}
    mkdir -p -- "$target_dir" || exit 1
    stage_dir=$(mktemp -d "$target_dir/.vcp_download.XXXXXXXX") || exit 1
    trap "rm -rf -- \"\$stage_dir\"" EXIT
    trap "exit 1" HUP INT TERM
    staged_file=$stage_dir/${local_file##*/}
    if ! vcp -v "$remote_file" "$staged_file"; then
        echo "ERROR: download failed for $remote_file" >&2
        exit 1
    fi
    if ! target_hash=$(local_hash "$staged_file"); then
        echo "ERROR: could not hash the staged download for $remote_file" >&2
        exit 1
    fi
    if [ "$source_hash" != "$target_hash" ]; then
        echo "ERROR: $algorithm checksum mismatch for $remote_file (CANFAR $source_hash, download $target_hash)" >&2
        exit 1
    fi
    if [ -d "$local_file" ] || ! mv -f -- "$staged_file" "$local_file"; then
        echo "ERROR: could not install verified file at $local_file" >&2
        exit 1
    fi
    echo "Downloaded and verified: $local_file ($algorithm $target_hash)"
' _
}

resolve_checksum_python() {
    local vcp_path vcp_header vcp_python candidate
    local -a candidates
    if [ -n "${CHECKSUM_PYTHON:-}" ]; then
        candidates=("$CHECKSUM_PYTHON")
    else
        if vcp_path=$(whence -p vcp) && IFS= read -r vcp_header < "$vcp_path"; then
            vcp_python=${vcp_header#\#!}
            if [[ "$vcp_header" == '#!/'* && "${vcp_python:t}" == python* && -x "$vcp_python" ]]; then
                candidates+=("$vcp_python")
            fi
        fi
        candidates+=(python3)
    fi
    for candidate in "${candidates[@]}"; do
        if "$candidate" -c 'import requests, vos' >/dev/null 2>&1; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done
    echo "ERROR: no checksum Python with requests and vos. Set CHECKSUM_PYTHON to the Python in your CADC environment." >&2
    return 1
}

set -o pipefail
checksum_dir=$(mktemp -d "${TMPDIR:-/tmp}/vcp_canfar_checksums.XXXXXXXX") || exit 1
export checksum_dir
trap 'rm -rf -- "$checksum_dir"' EXIT
trap 'exit 1' HUP INT TERM
echo "Downloading from CANFAR with TRANSFER_JOBS=$TRANSFER_JOBS..."
download_status=0
queue_downloads | run_queue || download_status=1

pending=("$checksum_dir"/pending.*(N))
if (( ${#pending[@]} )); then
    if checksum_python=$(resolve_checksum_python); then
        echo "Requesting CANFAR SHA-256 for ${#pending[@]} files without VOSpace MD5..."
        "$checksum_python" "${0:A:h}/canfar_checksum_manifest.py" download-queue \
            "$checksum_dir" "$REMOTE_PRODUCT_ROOT" > "$checksum_dir/verified.queue" || download_status=1
        if [ -s "$checksum_dir/verified.queue" ]; then
            run_queue < "$checksum_dir/verified.queue" || download_status=1
        fi
    else
        download_status=1
    fi
fi

if [ "$download_status" -ne 0 ]; then
    echo "Download finished with errors." >&2
    exit 1
fi

echo "Download finished. All available selected files verified."
