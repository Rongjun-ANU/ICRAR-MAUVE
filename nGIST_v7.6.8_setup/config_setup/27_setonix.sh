#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

. ./27_galaxies.sh

run_mode="both"
if [[ $# -gt 0 ]]; then
  last_arg="${!#}"
  case "$last_arg" in
    normal|7000)
      run_mode="$last_arg"
      set -- "${@:1:$(($# - 1))}"
      ;;
  esac
fi

if [[ $# -gt 0 && "$run_mode" == "both" ]] && ! looks_like_galid "${!#}"; then
  echo "ERROR: unknown run mode: ${!#}" >&2
  echo "Usage: $0 [GALID ...] [normal|7000]" >&2
  exit 2
fi

select_galids "$@"
GALIDS=("${SELECTED_GALIDS[@]}")
submitted_count=0
skipped_count=0

for galid in "${GALIDS[@]}"; do
  slurm_file="${galid}_v3tk_v7.6.8_setonix.slurm"
  slurm_file_7000="${galid}_v3tk_v7.6.8_7000_setonix.slurm"
  missing=()

  if [[ "$run_mode" != "7000" && ! -f "$slurm_file" ]]; then
    missing+=("$slurm_file")
  fi

  if [[ "$run_mode" != "normal" && ! -f "$slurm_file_7000" ]]; then
    missing+=("$slurm_file_7000")
  fi

  if [[ "${#missing[@]}" -gt 0 ]]; then
    echo "WARNING: skipping ${galid}: missing generated slurm file(s): ${missing[*]}" >&2
    skipped_count=$((skipped_count + 1))
    continue
  fi

  if [[ "$run_mode" != "7000" ]]; then
    echo "Submitting ${slurm_file}"
    sbatch "$slurm_file"
    submitted_count=$((submitted_count + 1))
  fi

  if [[ "$run_mode" != "normal" ]]; then
    echo "Submitting ${slurm_file_7000}"
    sbatch "$slurm_file_7000"
    submitted_count=$((submitted_count + 1))
  fi
done

if [[ "$submitted_count" -eq 0 ]]; then
  echo "WARNING: no generated slurm files available to submit; nothing to do." >&2
  exit 0
fi

echo "Submitted ${submitted_count} Setonix jobs; skipped ${skipped_count} unavailable run ID(s)."
