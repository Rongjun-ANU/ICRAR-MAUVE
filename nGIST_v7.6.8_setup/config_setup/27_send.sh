#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 USERNAME [HOST]" >&2
  echo "Example: $0 rhuang" >&2
  echo "Example with explicit host: $0 rhuang setonix.pawsey.org.au" >&2
}

if [[ $# -lt 1 || $# -gt 2 ]]; then
  usage
  exit 2
fi

cd "$(dirname "$0")"

REMOTE_USER="$1"
REMOTE_HOST="${2:-setonix.pawsey.org.au}"

if [[ "$REMOTE_USER" == *@* ]]; then
  REMOTE_LOGIN="$REMOTE_USER"
else
  REMOTE_LOGIN="${REMOTE_USER}@${REMOTE_HOST}"
fi

REMOTE_RUN_DIR="/software/projects/pawsey1308/ngist_supplementary_public/ngistTutorial"
REMOTE_CONFIG_DIR="${REMOTE_RUN_DIR}/configFiles"

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

yaml_files=()
slurm_files=()

for galid in "${GALIDS[@]}"; do
  yaml_file="${galid}_MAUVE_MasterConfig_v7.6.8_setonix.yaml"
  slurm_file="${galid}_v3tk_v7.6.8_setonix.slurm"

  if [[ ! -f "$yaml_file" ]]; then
    echo "ERROR: missing YAML file: $yaml_file" >&2
    exit 1
  fi

  if [[ ! -f "$slurm_file" ]]; then
    echo "ERROR: missing slurm file: $slurm_file" >&2
    exit 1
  fi

  yaml_files+=("$yaml_file")
  slurm_files+=("$slurm_file")
done

control_files=(
  27_setonix.sh
  27_status.sh
)

for control_file in "${control_files[@]}"; do
  if [[ ! -f "$control_file" ]]; then
    echo "ERROR: missing control script: $control_file" >&2
    exit 1
  fi
done

echo "Remote login: ${REMOTE_LOGIN}"
echo "Creating remote directories if needed"
ssh "$REMOTE_LOGIN" "mkdir -p '$REMOTE_CONFIG_DIR' '$REMOTE_RUN_DIR'"

echo "Sending ${#yaml_files[@]} YAML files to ${REMOTE_CONFIG_DIR}"
COPYFILE_DISABLE=1 tar -cf - "${yaml_files[@]}" | ssh "$REMOTE_LOGIN" \
  "tar -xf - -C '$REMOTE_CONFIG_DIR' && find '$REMOTE_CONFIG_DIR' -maxdepth 1 -type f -name '._*' -delete"

echo "Sending ${#slurm_files[@]} slurm scripts plus control scripts to ${REMOTE_RUN_DIR}"
COPYFILE_DISABLE=1 tar -cf - "${slurm_files[@]}" "${control_files[@]}" | ssh "$REMOTE_LOGIN" \
  "tar -xf - -C '$REMOTE_RUN_DIR' && chmod +x '$REMOTE_RUN_DIR'/27_setonix.sh '$REMOTE_RUN_DIR'/27_status.sh '$REMOTE_RUN_DIR'/*_v3tk_v7.6.8_setonix.slurm && find '$REMOTE_RUN_DIR' -maxdepth 1 -type f -name '._*' -delete"

echo "Done."
