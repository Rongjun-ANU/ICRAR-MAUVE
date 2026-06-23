#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

. ./27_galaxies.sh
select_galids "$@"
GALIDS=("${SELECTED_GALIDS[@]}")
submitted_count=0
skipped_count=0

for galid in "${GALIDS[@]}"; do
  slurm_file="${galid}_v3tk_v7.6.8_setonix.slurm"
  slurm_file_7000="${galid}_v3tk_v7.6.8_7000_setonix.slurm"
  missing=()

  if [[ ! -f "$slurm_file" ]]; then
    missing+=("$slurm_file")
  fi

  if [[ ! -f "$slurm_file_7000" ]]; then
    missing+=("$slurm_file_7000")
  fi

  if [[ "${#missing[@]}" -gt 0 ]]; then
    echo "WARNING: skipping ${galid}: missing generated slurm file(s): ${missing[*]}" >&2
    skipped_count=$((skipped_count + 1))
    continue
  fi

  echo "Submitting ${slurm_file}"
  sbatch "$slurm_file"
  submitted_count=$((submitted_count + 1))

  echo "Submitting ${slurm_file_7000}"
  sbatch "$slurm_file_7000"
  submitted_count=$((submitted_count + 1))
done

if [[ "$submitted_count" -eq 0 ]]; then
  echo "WARNING: no generated slurm file pairs available to submit; nothing to do." >&2
  exit 0
fi

echo "Submitted ${submitted_count} Setonix jobs; skipped ${skipped_count} unavailable run ID(s)."
