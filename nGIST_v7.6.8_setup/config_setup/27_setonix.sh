#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

GALIDS=(
  IC3392
  NGC4064
  NGC4189
  NGC4192
  NGC4293
  NGC4294
  NGC4298
  NGC4302
  NGC4330
  NGC4351
  NGC4383
  NGC4388
  NGC4394
  NGC4396
  NGC4402
  NGC4405
  NGC4419
  NGC4457
  NGC4501
  NGC4522
  NGC4567_8
  NGC4580
  NGC4606
  NGC4607
  NGC4694
  NGC4698
)

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
