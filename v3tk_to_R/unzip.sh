#!/usr/bin/env bash
set -euo pipefail

SECONDS=0
trap 'printf "Runtime: %02d:%02d:%02d\n" $((SECONDS/3600)) $(((SECONDS%3600)/60)) $((SECONDS%60))' EXIT

shopt -s nullglob

# 1) Handle *.gz.zip -> unzip directly into pwd (producing *.gz)
gz_zips=( ./*.gz.zip )
if (( ${#gz_zips[@]} > 0 )); then
    for zip in "${gz_zips[@]}"; do
        echo "Found: $(basename "$zip")"
        # -j: extract into pwd (ignore any paths inside the zip)
        unzip -o -q -j "$zip" -d .
    done
fi

# 2) Existing behavior for other *.zip -> unzip into a folder named after the zip
zips=( ./*.zip )
if (( ${#zips[@]} == 0 )); then
    echo "No *.zip files found in: $(pwd)" >&2
    exit 0
fi

for zip in "${zips[@]}"; do
    [[ "$zip" == *.gz.zip ]] && continue
    echo "Found: $(basename "$zip")"
    dest="$(basename "${zip%.zip}")"
    mkdir -p "$dest"
    unzip -o -q "$zip" -d "$dest"
done
