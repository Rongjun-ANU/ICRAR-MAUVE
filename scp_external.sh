#!/bin/bash

# Define the local base directory on the external disk
LOCAL_DIR="/Volumes/Untitled/v3tk_v7.6.8"

# Define the remote Setonix host
REMOTE_HOST="rhuang@setonix.pawsey.org.au"
REMOTE_BASE_DIR="/scratch/pawsey1308/mauve/products/v3tk_v7.6.8"

# List of all 26 GALIDs in alphabetical order
ALL_GALAXIES=(
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

if [ "$#" -gt 0 ]; then
    GALAXIES=("$@")
else
    GALAXIES=("${ALL_GALAXIES[@]}")
fi

# Verify that the external disk is mounted
if [ ! -d "/Volumes/Untitled" ]; then
    echo "Error: External disk '/Volumes/Untitled' not found."
    exit 1
fi

# Ensure base directory exists
mkdir -p "${LOCAL_DIR}"

# Set the maximum number of concurrent workers. 
# You might want this a bit lower than 40 for copying entire directories to avoid overwhelming the network/disk.
BATCH_SIZE=10

export LOCAL_DIR REMOTE_HOST REMOTE_BASE_DIR

# Use xargs to maintain a true "rolling queue" of parallel workers
printf "%s\n" "${GALAXIES[@]}" | xargs -P "$BATCH_SIZE" -I {} bash -c '
    GALID="$1"
    echo "Downloading folder for ${GALID}..."
    
    mkdir -p "${LOCAL_DIR}/${GALID}"

    # Grabs the entire folder over a single SSH connection
    rsync -avz --progress \
        "${REMOTE_HOST}:${REMOTE_BASE_DIR}/${GALID}/" \
        "${LOCAL_DIR}/${GALID}/"
    
    echo "Finished ${GALID}"
' _ {}

echo "Done!"
