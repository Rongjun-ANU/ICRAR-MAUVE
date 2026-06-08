#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
qc_py="${script_dir}/QC_ngist_v3tk_v768.py"
product_root="${1:-$PWD}"

if ! find "$product_root" -mindepth 2 -maxdepth 2 -name '*_sfh_maps.fits' -print -quit | grep -q .; then
    if [ -d "${product_root}/v3tk_v7.6.8" ]; then
        product_root="${product_root}/v3tk_v7.6.8"
    fi
fi

jobs="${JOBS:-${QC_JOBS:-}}"
if [ -z "$jobs" ]; then
    jobs="$(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)"
fi

find "$product_root" -mindepth 2 -maxdepth 2 -name '*_sfh_maps.fits' -print |
    sed -E 's#.*/([^/]+)/[^/]+_sfh_maps\.fits#\1#' |
    sort -u |
    xargs -P "$jobs" -I {} /opt/miniconda3/envs/ICRAR/bin/python "$qc_py" {}
