#!/bin/zsh

# Download v3tk_v7.6.8 products from CANFAR
# cadc-get-cert -u RongjunHuang

LOCAL_BASE="/Users/Igniz/Desktop/ICRAR/further"
REMOTE_PRODUCT_ROOT="arc:projects/mauve/products"
REMOTE_LOG_BASE="arc:home/RongjunHuang/ICRAR/further"
ALL_RUNS=("v3tk_v7.6.8" "v3tk_v7.6.8_7000")

ALL_GALAXIES=(
    "IC3392" "NGC4064" "NGC4189" "NGC4192" "NGC4293" "NGC4294" 
    "NGC4298" "NGC4302" "NGC4330" "NGC4351" "NGC4383" "NGC4388" 
    "NGC4394" "NGC4396" "NGC4402" "NGC4405" "NGC4419" "NGC4457" 
    "NGC4501" "NGC4522" "NGC4567_8" "NGC4580" "NGC4606" "NGC4607" 
    "NGC4694" "NGC4698"
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
export LOCAL_BASE REMOTE_PRODUCT_ROOT REMOTE_LOG_BASE

# Use xargs to maintain a true "rolling queue" of parallel workers
for RUN in "${RUNS[@]}"; do
echo "Downloading ${RUN} from CANFAR..."
REMOTE_GAL_BASE="${REMOTE_PRODUCT_ROOT}/${RUN}"
export RUN REMOTE_GAL_BASE

printf "%s\n" "${GALAXIES[@]}" | xargs -P "$BATCH_SIZE" -I {} bash -c '
    GALID="$1"
    echo "Downloading ${RUN}/${GALID}..."

    mkdir -p "${LOCAL_BASE}/${RUN}/${GALID}"
    mkdir -p "${LOCAL_BASE}/mass_logs" "${LOCAL_BASE}/sfr_logs" "${LOCAL_BASE}/proxy_ewha_logs"

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
        vcp -v \
            "${REMOTE_GAL_BASE}/${GALID}/${PRODUCT}" \
            "${LOCAL_BASE}/${RUN}/${GALID}/" \
            || echo "WARNING: failed to download ${RUN}/${GALID}/${PRODUCT}"
    done

    vcp -v "${REMOTE_LOG_BASE}/mass_logs/${GALID}.log" "${LOCAL_BASE}/mass_logs/" 2>/dev/null
    vcp -v "${REMOTE_LOG_BASE}/sfr_logs/${GALID}.log" "${LOCAL_BASE}/sfr_logs/" 2>/dev/null
    vcp -v "${REMOTE_LOG_BASE}/proxy_ewha_logs/${GALID}.log" "${LOCAL_BASE}/proxy_ewha_logs/" 2>/dev/null

    echo "Finished ${RUN}/${GALID}"
' _ {}

vcp -v "${REMOTE_GAL_BASE}/*_logs" "${LOCAL_BASE}/${RUN}/" \
    || echo "WARNING: failed to download ${RUN} *_logs folders"
done

echo "Download from CANFAR finished!"
