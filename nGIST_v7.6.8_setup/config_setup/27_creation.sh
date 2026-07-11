#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

. ./27_galaxies.sh
select_galids "$@"
GALIDS=("${SELECTED_GALIDS[@]}")

SLURM_TEMPLATES=(
  "v3tk_v7.6.8_setonix.slurm"
  "v3tk_v7.6.8_7000_setonix.slurm"
)
LONG_GALIDS=(
  NGC4293
  NGC4298
  NGC4302
  NGC4383
  NGC4419
  NGC4457
)
HIGHMEM_GALIDS=(
  NGC4192
  NGC4254
  NGC4298
  NGC4321
  NGC4330
  NGC4380
  NGC4396
  NGC4501
  NGC4535
  NGC4567_8
  NGC4569
  NGC4698
)
CUBE_CENTERS_CSV="cube_centers_v3tk.csv"
CUBE_SIZES_CSV="cube_sizes_v3tk.csv"
DEFAULT_ICRAR_PYTHON="/opt/miniconda3/envs/ICRAR/bin/python"

if [[ -n "${PYTHON:-}" ]]; then
  PYTHON_BIN="$PYTHON"
elif [[ -x "$DEFAULT_ICRAR_PYTHON" ]]; then
  PYTHON_BIN="$DEFAULT_ICRAR_PYTHON"
elif command -v python >/dev/null 2>&1; then
  PYTHON_BIN="$(command -v python)"
elif command -v python3 >/dev/null 2>&1; then
  PYTHON_BIN="$(command -v python3)"
else
  echo "ERROR: no Python interpreter found; set PYTHON=/path/to/python" >&2
  exit 1
fi

cube_galid_for() {
  local galid="$1"
  echo "$galid"
}

galid_is_in_list() {
  local candidate="$1"
  shift
  local galid

  for galid in "$@"; do
    if [[ "$candidate" == "$galid" ]]; then
      return 0
    fi
  done

  return 1
}

galid_has_csv_row() {
  local candidate="$1"
  local csv_path="$2"
  local csv_id

  [[ -f "$csv_path" ]] || return 1

  while IFS=, read -r csv_id _; do
    if [[ "$csv_id" == "$candidate" ]]; then
      return 0
    fi
  done < "$csv_path"

  return 1
}

filter_galids_with_cube_csv_rows() {
  local galid
  local cube_galid
  local missing
  local filtered_galids=()

  for galid in "${GALIDS[@]}"; do
    cube_galid="$(cube_galid_for "$galid")"
    missing=()

    if ! galid_has_csv_row "$cube_galid" "$CUBE_CENTERS_CSV"; then
      missing+=("$CUBE_CENTERS_CSV")
    fi

    if ! galid_has_csv_row "$cube_galid" "$CUBE_SIZES_CSV"; then
      missing+=("$CUBE_SIZES_CSV")
    fi

    if [[ "${#missing[@]}" -gt 0 ]]; then
      echo "WARNING: skipping ${galid}: missing row in ${missing[*]}" >&2
      continue
    fi

    filtered_galids+=("$galid")
  done

  GALIDS=("${filtered_galids[@]}")
}

if [[ ! -x ./make_gist_config_try.py ]]; then
  echo "ERROR: ./make_gist_config_try.py is missing or is not executable" >&2
  exit 1
fi

if [[ ! -x "$PYTHON_BIN" ]]; then
  echo "ERROR: Python interpreter is not executable: $PYTHON_BIN" >&2
  exit 1
fi

for slurm_template in "${SLURM_TEMPLATES[@]}"; do
  if [[ ! -f "$slurm_template" ]]; then
    echo "ERROR: missing slurm template: $slurm_template" >&2
    exit 1
  fi
done

filter_galids_with_cube_csv_rows

for galid in "${GALIDS[@]}"; do
  echo "Creating YAML configs for ${galid}"
  "$PYTHON_BIN" ./make_gist_config_try.py "$galid"

  rm -f \
    "${galid}_MAUVE_MasterConfig_v7.6.8_setonix.slurm" \
    "${galid}_MAUVE_MasterConfig_v7.6.8_7000_setonix.slurm"

  for slurm_template in "${SLURM_TEMPLATES[@]}"; do
    slurm_out="${galid}_${slurm_template}"
    echo "Creating slurm script ${slurm_out}"
    sed "s/GALID/${galid}/g" "$slurm_template" > "$slurm_out"

    if galid_is_in_list "$galid" "${HIGHMEM_GALIDS[@]}"; then
      echo "Updating ${slurm_out} for Setonix highmem queue"
      tmp_slurm="${slurm_out}.tmp"
      sed \
        -e 's/^#SBATCH --partition=work$/#SBATCH --partition=highmem/' \
        -e 's/^#SBATCH --mem=230G$/#SBATCH --mem=980G/' \
        -e 's/^#SBATCH --time=24:00:00$/#SBATCH --time=96:00:00/' \
        "$slurm_out" > "$tmp_slurm"
      mv "$tmp_slurm" "$slurm_out"
    elif galid_is_in_list "$galid" "${LONG_GALIDS[@]}"; then
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
done

echo "Created $((${#GALIDS[@]} * 2)) YAML files and $((${#GALIDS[@]} * 2)) slurm scripts."
