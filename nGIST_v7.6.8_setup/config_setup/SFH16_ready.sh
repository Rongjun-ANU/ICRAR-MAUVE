#!/usr/bin/env bash
# Run on Setonix before sbatch; existing folders are skipped unchanged.
# --check strictly verifies a previously prepared run instead of skipping it.
set -euo pipefail
script_dir="$(cd "$(dirname "$0")" && pwd)"
. "${script_dir}/SFH16_galaxies.sh"
mode=ready
if [[ "${1:-}" == --check ]]; then
  mode=check
  shift
fi
select_sfh16_galids "$@"
sfh16_python "${script_dir}/SFH16_workflow.py" "$mode" "${SELECTED_GALIDS[@]}"
