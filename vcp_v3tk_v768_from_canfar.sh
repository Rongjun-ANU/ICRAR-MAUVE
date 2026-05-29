#!/bin/zsh

# Download v3tk_v7.6.8 products from CANFAR
# cadc-get-cert -u RongjunHuang

echo "Downloading v3tk_v7.6.8 from CANFAR..."
LOCAL_BASE="/Users/Igniz/Desktop/ICRAR/further"
REMOTE_GAL_BASE="arc:home/RongjunHuang/ICRAR/further/v3tk_v7.6.8"
REMOTE_LOG_BASE="arc:home/RongjunHuang/ICRAR/further"

ALL_GALAXIES=(
    "IC3392" "NGC4064" "NGC4189" "NGC4192" "NGC4293" "NGC4294" 
    "NGC4298" "NGC4302" "NGC4330" "NGC4351" "NGC4383" "NGC4388" 
    "NGC4394" "NGC4396" "NGC4402" "NGC4405" "NGC4419" "NGC4457" 
    "NGC4501" "NGC4522" "NGC4567_8" "NGC4580" "NGC4606" "NGC4607" 
    "NGC4694" "NGC4698"
)

if [ "$#" -gt 0 ]; then
    GALAXIES=("$@")
else
    GALAXIES=("${ALL_GALAXIES[@]}")
fi

# Set the maximum number of concurrent workers
BATCH_SIZE=40

# Export variables so the xargs subprocess can see them
export LOCAL_BASE REMOTE_GAL_BASE REMOTE_LOG_BASE

# Use xargs to maintain a true "rolling queue" of parallel workers
printf "%s\n" "${GALAXIES[@]}" | xargs -P "$BATCH_SIZE" -I {} bash -c '
    GALID="$1"
    echo "Downloading ${GALID}..."

    mkdir -p "${LOCAL_BASE}/v3tk_v7.6.8/${GALID}"
    mkdir -p "${LOCAL_BASE}/mass_logs" "${LOCAL_BASE}/sfr_logs" "${LOCAL_BASE}/proxy_ewha_logs"

    vcp -v \
        "${REMOTE_GAL_BASE}/${GALID}/${GALID}_gas_bin_maps_further.fits" \
        "${REMOTE_GAL_BASE}/${GALID}/${GALID}_proxy_EW_maps_further.fits" \
        "${REMOTE_GAL_BASE}/${GALID}/${GALID}_spatial_binning_maps_further.fits" \
        "${LOCAL_BASE}/v3tk_v7.6.8/${GALID}/"

    vcp -v "${REMOTE_LOG_BASE}/mass_logs/${GALID}.log" "${LOCAL_BASE}/mass_logs/" 2>/dev/null
    vcp -v "${REMOTE_LOG_BASE}/sfr_logs/${GALID}.log" "${LOCAL_BASE}/sfr_logs/" 2>/dev/null
    vcp -v "${REMOTE_LOG_BASE}/proxy_ewha_logs/${GALID}.log" "${LOCAL_BASE}/proxy_ewha_logs/" 2>/dev/null

    echo "Finished ${GALID}"
' _ {}

echo "Download from CANFAR finished!"