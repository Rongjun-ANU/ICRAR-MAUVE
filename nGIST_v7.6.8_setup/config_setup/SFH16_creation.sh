#!/usr/bin/env bash
# SFH-only jobs: work partition, 230G, 24 hours; preserve baseline CPU count.
set -euo pipefail
cd "$(dirname "$0")"
. ./SFH16_galaxies.sh
select_sfh16_galids "$@"
sfh16_python ./SFH16_workflow.py create "${SELECTED_GALIDS[@]}"
