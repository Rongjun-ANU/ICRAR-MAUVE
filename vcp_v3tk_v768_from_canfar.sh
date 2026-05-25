#!/bin/zsh

# Download v3tk_v7.6.8 products from CANFAR
cadc-get-cert -u RongjunHuang

echo "Downloading v3tk_v7.6.8 from CANFAR..."
LOCAL_BASE="/Users/Igniz/Desktop/ICRAR/further"
REMOTE_GAL_BASE="arc:home/RongjunHuang/ICRAR/further/v3tk_v7.6.8"
REMOTE_LOG_BASE="arc:home/RongjunHuang/ICRAR/further"

GALAXIES=(
    "IC3392" "NGC4064" "NGC4189" "NGC4192" "NGC4293" "NGC4294" 
    "NGC4298" "NGC4302" "NGC4330" "NGC4351" "NGC4383" "NGC4388" 
    "NGC4394" "NGC4396" "NGC4402" "NGC4405" "NGC4419" "NGC4457" 
    "NGC4501" "NGC4522" "NGC4567_8" "NGC4580" "NGC4606" "NGC4607" 
    "NGC4694" "NGC4698"
)

# Set a batch size to limit concurrent connections
BATCH_SIZE=5
count=0

for GALID in "${GALAXIES[@]}"; do
    (
        echo "Downloading ${GALID}..."

        mkdir -p "${LOCAL_BASE}/v3tk_v7.6.8/${GALID}"
        mkdir -p "${LOCAL_BASE}/mass_logs" "${LOCAL_BASE}/sfr_logs" "${LOCAL_BASE}/proxy_ewha_logs"

        FILES=(
            "${GALID}_gas_bin_maps_further.fits"
            "${GALID}_proxy_EW_maps_further.fits"
            "${GALID}_spatial_binning_maps_further.fits"
        )
        
        for FILE in "${FILES[@]}"; do
            vcp -v "${REMOTE_GAL_BASE}/${GALID}/${FILE}" "${LOCAL_BASE}/v3tk_v7.6.8/${GALID}/"
        done

        vcp -v "${REMOTE_LOG_BASE}/mass_logs/${GALID}.log" "${LOCAL_BASE}/mass_logs/" 2>/dev/null
        vcp -v "${REMOTE_LOG_BASE}/sfr_logs/${GALID}.log" "${LOCAL_BASE}/sfr_logs/" 2>/dev/null
        vcp -v "${REMOTE_LOG_BASE}/proxy_ewha_logs/${GALID}.log" "${LOCAL_BASE}/proxy_ewha_logs/" 2>/dev/null

        echo "Finished ${GALID}"
    ) &

    # Wait for the current batch to finish before starting the next
    count=$((count + 1))
    if (( count % BATCH_SIZE == 0 )); then
        wait
    fi
done

# Wait for any remaining background jobs to finish
wait
echo "Download from CANFAR finished!"