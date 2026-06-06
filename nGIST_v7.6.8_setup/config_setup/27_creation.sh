#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

. ./27_galaxies.sh
select_galids "$@"
GALIDS=("${SELECTED_GALIDS[@]}")

SLURM_TEMPLATE="v3tk_v7.6.8_setonix.slurm"
HIGHMEM_GALIDS=(
  NGC4192
  NGC4254
  NGC4321
  NGC4501
  NGC4535
  NGC4569
)

LONG_GALIDS=(
  NGC4293
  NGC4298
  NGC4302
  NGC4330
  NGC4383
  NGC4396
  NGC4419
  NGC4457
  NGC4698
)

uses_highmem_queue() {
  local galid="$1"
  local highmem_galid

  for highmem_galid in "${HIGHMEM_GALIDS[@]}"; do
    if [[ "$galid" == "$highmem_galid" ]]; then
      return 0
    fi
  done

  return 1
}

uses_long_queue() {
  local galid="$1"
  local long_galid

  for long_galid in "${LONG_GALIDS[@]}"; do
    if [[ "$galid" == "$long_galid" ]]; then
      return 0
    fi
  done

  return 1
}

if [[ ! -x ./make_gist_config_try.py ]]; then
  echo "ERROR: ./make_gist_config_try.py is missing or is not executable" >&2
  exit 1
fi

if [[ ! -f "$SLURM_TEMPLATE" ]]; then
  echo "ERROR: missing slurm template: $SLURM_TEMPLATE" >&2
  exit 1
fi

for galid in "${GALIDS[@]}"; do
  echo "Creating YAML for ${galid}"
  ./make_gist_config_try.py "$galid"

  slurm_out="${galid}_v3tk_v7.6.8_setonix.slurm"
  echo "Creating slurm script ${slurm_out}"
  sed "s/GALID/${galid}/g" "$SLURM_TEMPLATE" > "$slurm_out"

  if uses_highmem_queue "$galid"; then
    echo "Updating ${slurm_out} for Setonix highmem queue"
    tmp_slurm="${slurm_out}.tmp"
    sed \
      -e 's/^#SBATCH --partition=work$/#SBATCH --partition=highmem/' \
      -e 's/^#SBATCH --mem=230G$/#SBATCH --mem=980G/' \
      -e 's/^#SBATCH --time=24:00:00$/#SBATCH --time=96:00:00/' \
      "$slurm_out" > "$tmp_slurm"
    mv "$tmp_slurm" "$slurm_out"
  elif uses_long_queue "$galid"; then
    echo "Updating ${slurm_out} for Setonix long queue"
    tmp_slurm="${slurm_out}.tmp"
    sed \
      -e 's/^#SBATCH --partition=work$/#SBATCH --partition=long/' \
      -e 's/^#SBATCH --time=24:00:00$/#SBATCH --time=96:00:00/' \
      "$slurm_out" > "$tmp_slurm"
    mv "$tmp_slurm" "$slurm_out"
  fi

  chmod +x "$slurm_out"
done

echo "Created ${#GALIDS[@]} YAML files and ${#GALIDS[@]} slurm scripts."
