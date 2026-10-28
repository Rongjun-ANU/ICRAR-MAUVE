#!/usr/bin/env bash
# Deliberately separate from the ordinary 27_* run selection.
SFH16_GALIDS=(NGC4189 NGC4380 NGC4450 NGC4567_8 NGC4689)
select_sfh16_galids() {
  SELECTED_GALIDS=()
  if [[ $# -eq 0 ]]; then
    SELECTED_GALIDS=("${SFH16_GALIDS[@]}")
    return
  fi
  local requested known found
  for requested in "$@"; do
    found=false
    for known in "${SFH16_GALIDS[@]}"; do
      [[ "$requested" != "$known" ]] || found=true
    done
    if [[ "$found" != true ]]; then
      echo "ERROR: select from: ${SFH16_GALIDS[*]} (got $requested)" >&2
      return 2
    fi
    SELECTED_GALIDS+=("$requested")
  done
}

sfh16_python() {
  if [[ -n "${PYTHON:-}" ]]; then
    "$PYTHON" "$@"
  elif command -v ngistenv1308 >/dev/null 2>&1; then
    ngistenv1308 python "$@"
  elif [[ -x /opt/miniconda3/envs/ICRAR/bin/python ]]; then
    /opt/miniconda3/envs/ICRAR/bin/python "$@"
  else
    python3 "$@"
  fi
}

# Sourcing this file supplies helpers only. Executing it submits SFH16 jobs.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  set -euo pipefail
  cd "$(dirname "$0")"
  select_sfh16_galids "$@"
  # Check the complete selection before submitting any jobs.
  for galid in "${SELECTED_GALIDS[@]}"; do
    slurm_file="${galid}_v3tk_v7.6.8_7000_SFH16_setonix.slurm"
    [[ -s "$slurm_file" ]] || {
      echo "ERROR: missing Slurm script: $slurm_file" >&2
      exit 1
    }
  done
  for galid in "${SELECTED_GALIDS[@]}"; do
    slurm_file="${galid}_v3tk_v7.6.8_7000_SFH16_setonix.slurm"
    echo "Submitting ${slurm_file}"
    sbatch "$slurm_file"
  done
  echo "Submitted ${#SELECTED_GALIDS[@]} SFH16 jobs."
fi
