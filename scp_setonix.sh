#!/bin/bash

# Define the local base directory
LOCAL_DIR="/Users/Igniz/Desktop/ICRAR/further/v3tk_v7.6.8"

# Define the remote Setonix host (adjust 'setonix' to your actual SSH config host if needed)
REMOTE_HOST="rhuang@setonix.pawsey.org.au"
REMOTE_BASE_DIR="/scratch/pawsey1308/mauve/products/v3tk_v7.6.8"

# List of all 26 GALIDs in alphabetical order
GALAXIES=(
    "IC3392"
    "NGC4064"
    "NGC4189"
    "NGC4192"
    "NGC4293"
    "NGC4294"
    "NGC4298"
    "NGC4302"
    "NGC4330"
    "NGC4351"
    "NGC4383"
    "NGC4388"
    "NGC4394"
    "NGC4396"
    "NGC4402"
    "NGC4405"
    "NGC4419"
    "NGC4457"
    "NGC4501"
    "NGC4522"
    "NGC4567_8"
    "NGC4580"
    "NGC4606"
    "NGC4607"
    "NGC4694"
    "NGC4698"
)

# Loop through each galaxy
for GALID in "${GALAXIES[@]}"; do
    echo "Downloading files for ${GALID}..."
    
    # Create the local directory for the galaxy if it doesn't exist
    mkdir -p "${LOCAL_DIR}/${GALID}"

    # Define the list of files to copy for this galaxy
    FILES=(
        "CONFIG"
        "${GALID}_kin_maps.fits"
        "${GALID}_sfh_weights.fits"
        "${GALID}_gas_bin_maps.fits"
        "${GALID}_mask.fits"
        "${GALID}_spatial_binning_maps.fits"
        "${GALID}_gas_spaxel_maps.fits"
        "${GALID}_sfh_maps.fits"
        "LOGFILE"
    )

    # Because Pawsey's banner breaks scp/rsync easily, another alternative is forcing a pseudo-terminal (-t) with ssh and cat, 
    # but the easiest workaround for annoying banners using rsync is to use ssh with the -e flag to suppress banners, or just pass --ignore-missing-args
    # Wait, the best fix for "protocol version mismatch -- is your shell clean?" is to use sftp, which is inherently immune to bashrc banners.
    
    # macOS recently updated `scp` to use the strict SFTP protocol by default, which crashes instantly 
    # if the server prints any text (like the Setonix 80-line warning banner).
    # The `-O` flag forces the legacy SCP protocol, which is much more tolerant of server text banners.
    for FILE in "${FILES[@]}"; do
        scp -O -q "${REMOTE_HOST}:${REMOTE_BASE_DIR}/${GALID}/${FILE}" "${LOCAL_DIR}/${GALID}/"
    done
done

echo "Done!"
