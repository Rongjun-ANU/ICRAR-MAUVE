#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

. ./27_galaxies.sh
select_galids "$@"
GALIDS=("${SELECTED_GALIDS[@]}")

for galid in "${GALIDS[@]}"; do
  slurm_script="${galid}_v3tk_v7.6.8_setonix.slurm"
  if [[ ! -f "$slurm_script" ]]; then
    echo "ERROR: missing slurm script: $slurm_script" >&2
    exit 1
  fi

  echo "Submitting ${slurm_script}"
  sbatch "$slurm_script"
done

echo "Submitted ${#GALIDS[@]} Setonix jobs."
