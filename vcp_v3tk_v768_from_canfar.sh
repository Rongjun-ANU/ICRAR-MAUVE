#!/bin/zsh

# Download v3tk_v7.6.8 products from CANFAR
# cadc-get-cert -u RongjunHuang

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

RUNS=()
GALAXIES=()
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

# Set the maximum number of concurrent workers
BATCH_SIZE=40

# Export variables so the xargs subprocess can see them
export LOCAL_BASE REMOTE_PRODUCT_ROOT

# Use xargs to maintain a true "rolling queue" of parallel workers
for RUN in "${RUNS[@]}"; do
echo "Downloading ${RUN} from CANFAR..."
REMOTE_GAL_BASE="${REMOTE_PRODUCT_ROOT}/${RUN}"
export RUN REMOTE_GAL_BASE

printf "%s\n" "${GALAXIES[@]}" | xargs -P "$BATCH_SIZE" -I {} bash -c '
    GALID="$1"
    SOURCE_GAL_DIR="${REMOTE_GAL_BASE}/${GALID}"
    TARGET_RUN_DIR="${LOCAL_BASE}/${RUN}"
    TARGET_GAL_DIR="${TARGET_RUN_DIR}/${GALID}"

    if ! vls "${SOURCE_GAL_DIR}" >/dev/null 2>&1; then
        echo "Skipping ${RUN}/${GALID}: source folder not found or inaccessible: ${SOURCE_GAL_DIR}"
        exit 0
    fi

    remote_md5() {
        local value
        value=$(vtag "$1" MD5 2>/dev/null) || return 1
        value=$(printf "%s" "${value}" | tr -d "\047\042[:space:]")
        [[ "${value}" =~ ^[[:xdigit:]]{32}$ ]] || return 1
        printf "%s\n" "${value}"
    }

    local_md5() {
        local value
        if command -v md5 >/dev/null 2>&1; then
            md5 -q -- "$1"
        elif command -v md5sum >/dev/null 2>&1; then
            value=$(md5sum -- "$1") || return 1
            printf "%s\n" "${value%% *}"
        else
            return 1
        fi
    }

    files_identical() {
        local remote_file="$1"
        local local_file="$2"
        local source_md5 target_md5

        [ -f "${local_file}" ] || return 1
        source_md5=$(remote_md5 "${remote_file}") || return 1
        target_md5=$(local_md5 "${local_file}") || return 1
        [ "${source_md5}" = "${target_md5}" ]
    }

    download_if_needed() {
        local remote_file="$1"
        local local_file="$2"

        if files_identical "${remote_file}" "${local_file}"; then
            echo "Skipping ${remote_file}: source and target are identical"
        elif ! vcp -v "${remote_file}" "${local_file}"; then
            echo "WARNING: failed to download ${remote_file}"
        fi
    }

    echo "Downloading ${RUN}/${GALID}..."

    mkdir -p "${TARGET_GAL_DIR}"
    mkdir -p "${TARGET_RUN_DIR}/mass_logs" "${TARGET_RUN_DIR}/sfr_logs" "${TARGET_RUN_DIR}/proxy_ewha_logs"

    PRODUCTS=(
        "CONFIG"
        "LOGFILE"
        "${GALID}_mask.fits"
        "${GALID}_sfh_maps.fits"
        "${GALID}_sfh_weights.fits"
        "${GALID}_kin_maps.fits"
        "${GALID}_gas_bin_maps_further.fits"
        "${GALID}_proxy_EW_maps_further.fits"
        "${GALID}_spatial_binning_maps_further.fits"
    )

    for PRODUCT in "${PRODUCTS[@]}"; do
        download_if_needed \
            "${SOURCE_GAL_DIR}/${PRODUCT}" \
            "${TARGET_GAL_DIR}/${PRODUCT}"
    done

    for LOG_DIR in mass_logs sfr_logs proxy_ewha_logs; do
        download_if_needed \
            "${REMOTE_GAL_BASE}/${LOG_DIR}/${GALID}.log" \
            "${TARGET_RUN_DIR}/${LOG_DIR}/${GALID}.log"
    done

    echo "Finished ${RUN}/${GALID}"
' _ {}
done

echo "Download from CANFAR finished!"
