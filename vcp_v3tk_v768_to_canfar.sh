# Upload local v3tk_v7.6.8 products to CANFAR
# cadc-get-cert -u RongjunHuang

LOCAL_ROOT="/Users/Igniz/Desktop/ICRAR/further"
REMOTE_ROOT="arc:projects/mauve/products"
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
export LOCAL_ROOT REMOTE_ROOT

# Use xargs to maintain a true "rolling queue" of parallel workers
for RUN in "${RUNS[@]}"; do
echo "Uploading ${RUN} to CANFAR..."
LOCAL_BASE="${LOCAL_ROOT}/${RUN}"
REMOTE_BASE="${REMOTE_ROOT}/${RUN}"
export LOCAL_BASE REMOTE_BASE

printf "%s\n" "${GALAXIES[@]}" | xargs -P "$BATCH_SIZE" -I {} bash -c '
    GALID="$1"
    echo "Uploading ${REMOTE_BASE}/${GALID}..."
    
    # Create the remote directory first 
    vmkdir "${REMOTE_BASE}/${GALID}" 2>/dev/null

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
        
    echo "Finished ${REMOTE_BASE}/${GALID}"
' _ {}
done

echo "Upload to CANFAR finished!"
