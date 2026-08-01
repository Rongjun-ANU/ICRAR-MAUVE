# Upload local v3tk_v7.6.8 products to CANFAR
# cadc-get-cert -u RongjunHuang
#
# Usage:
#   ./vcp_v3tk_v768_to_canfar.sh normal              # all galaxies, normal run
#   ./vcp_v3tk_v768_to_canfar.sh 7000                # all galaxies, 7000 run
#   ./vcp_v3tk_v768_to_canfar.sh normal NGC4321      # one galaxy, normal run
#   ./vcp_v3tk_v768_to_canfar.sh 7000 NGC4321        # one galaxy, 7000 run
#   ./vcp_v3tk_v768_to_canfar.sh NGC4321             # one galaxy, both runs
# With no arguments, all listed galaxies from both runs are uploaded.

LOCAL_ROOT="/Users/Igniz/Desktop/ICRAR/further"
REMOTE_ROOT="arc:projects/mauve/products"
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
export LOCAL_ROOT REMOTE_ROOT

# Use xargs to maintain a true "rolling queue" of parallel workers
for RUN in "${RUNS[@]}"; do
echo "Uploading ${RUN} to CANFAR..."
LOCAL_BASE="${LOCAL_ROOT}/${RUN}"
REMOTE_BASE="${REMOTE_ROOT}/${RUN}"
export LOCAL_BASE REMOTE_BASE

printf "%s\n" "${GALAXIES[@]}" | xargs -P "$BATCH_SIZE" -I {} bash -c '
    GALID="$1"
    SOURCE_GAL_DIR="${LOCAL_BASE}/${GALID}"
    TARGET_GAL_DIR="${REMOTE_BASE}/${GALID}"

    if [ ! -d "${SOURCE_GAL_DIR}" ]; then
        echo "Skipping ${RUN}/${GALID}: source folder not found: ${SOURCE_GAL_DIR}"
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
        local local_file="$1"
        local remote_file="$2"
        local source_md5 target_md5

        source_md5=$(local_md5 "${local_file}") || return 1
        target_md5=$(remote_md5 "${remote_file}") || return 1
        [ "${source_md5}" = "${target_md5}" ]
    }

    echo "Uploading ${REMOTE_BASE}/${GALID}..."
    
    # Create the remote directory first.
    vmkdir "${TARGET_GAL_DIR}" 2>/dev/null

    PRODUCTS=(
        "CONFIG"
        "${GALID}_sfh_maps.fits"
        "${GALID}_gas_bin_maps.fits"
        "${GALID}_sfh_weights.fits"
        "${GALID}_gas_spaxel_maps.fits"
        "${GALID}_spatial_binning_maps.fits"
        "${GALID}_kin_maps.fits"
        "LOGFILE"
        "${GALID}_mask.fits"
    )

    for PRODUCT in "${PRODUCTS[@]}"; do
        LOCAL_FILE="${SOURCE_GAL_DIR}/${PRODUCT}"
        REMOTE_FILE="${TARGET_GAL_DIR}/${PRODUCT}"

        if [ ! -f "${LOCAL_FILE}" ]; then
            echo "WARNING: source file not found; skipping ${LOCAL_FILE}"
        elif files_identical "${LOCAL_FILE}" "${REMOTE_FILE}"; then
            echo "Skipping ${LOCAL_FILE}: source and target are identical"
        elif ! vcp -v "${LOCAL_FILE}" "${REMOTE_FILE}"; then
            echo "WARNING: failed to upload ${LOCAL_FILE}"
        fi
    done
        
    echo "Finished ${REMOTE_BASE}/${GALID}"
' _ {}
done

echo "Upload to CANFAR finished!"
