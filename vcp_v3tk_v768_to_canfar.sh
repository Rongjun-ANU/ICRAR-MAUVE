# Upload local v3tk_v7.6.8 products to CANFAR
# cadc-get-cert -u RongjunHuang

echo "Uploading v3tk_v7.6.8 to CANFAR..."
LOCAL_BASE="/Users/Igniz/Desktop/ICRAR/further/v3tk_v7.6.8"
REMOTE_BASE="arc:home/RongjunHuang/ICRAR/further/v3tk_v7.6.8"

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
export LOCAL_BASE REMOTE_BASE

# Use xargs to maintain a true "rolling queue" of parallel workers
printf "%s\n" "${GALAXIES[@]}" | xargs -P "$BATCH_SIZE" -I {} bash -c '
    GALID="$1"
    echo "Uploading ${GALID}..."
    
    # Create the remote directory first 
    vmkdir "arc:home/RongjunHuang/ICRAR/further/v3tk_v7.6.8/${GALID}" 2>/dev/null

    # Upload files sequentially within this galaxy batch using a single command
    vcp -v \
        "${LOCAL_BASE}/${GALID}/CONFIG" \
        "${LOCAL_BASE}/${GALID}/${GALID}_sfh_maps.fits" \
        "${LOCAL_BASE}/${GALID}/${GALID}_gas_bin_maps.fits" \
        "${LOCAL_BASE}/${GALID}/${GALID}_sfh_weights.fits" \
        "${LOCAL_BASE}/${GALID}/${GALID}_gas_spaxel_maps.fits" \
        "${LOCAL_BASE}/${GALID}/${GALID}_spatial_binning_maps.fits" \
        "${LOCAL_BASE}/${GALID}/${GALID}_kin_maps.fits" \
        "${LOCAL_BASE}/${GALID}/LOGFILE" \
        "${LOCAL_BASE}/${GALID}/${GALID}_mask.fits" \
        "${REMOTE_BASE}/${GALID}/"
        
    echo "Finished ${GALID}"
' _ {}

echo "Upload to CANFAR finished!"