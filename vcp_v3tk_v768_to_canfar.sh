# Upload local v3tk_v7.6.8 products to CANFAR
cadc-get-cert -u RongjunHuang

echo "Uploading v3tk_v7.6.8 to CANFAR..."
LOCAL_BASE="/Users/Igniz/Desktop/ICRAR/further/v3tk_v7.6.8"
REMOTE_BASE="arc:home/RongjunHuang/ICRAR/further/v3tk_v7.6.8"

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
        echo "Uploading ${GALID}..."
        
        # Create the remote directory first (using vmkdir or similar if vcp doesn't auto-create it)
        vmkdir "arc:home/RongjunHuang/ICRAR/further/v3tk_v7.6.8/${GALID}" 2>/dev/null

        FILES=(
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
        
        for FILE in "${FILES[@]}"; do
            # Upload files sequentially within this galaxy batch
            vcp -v "${LOCAL_BASE}/${GALID}/${FILE}" "${REMOTE_BASE}/${GALID}/"
        done
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
echo "Upload to CANFAR finished!"