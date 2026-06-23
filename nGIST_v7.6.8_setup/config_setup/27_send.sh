#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 USERNAME [HOST] [GALID ...]" >&2
  echo "       $0 GALID [GALID ...]" >&2
  echo "Example: $0 rhuang" >&2
  echo "Example with default user for selected galaxies: $0 NGC4569" >&2
  echo "Example with explicit host: $0 rhuang setonix.pawsey.org.au" >&2
  echo "Example for selected galaxies: $0 rhuang NGC4383 NGC4419" >&2
  echo "Example for selected galaxies and explicit host: $0 rhuang setonix.pawsey.org.au NGC4383 NGC4419" >&2
}

if [[ $# -lt 1 ]]; then
  usage
  exit 2
fi

cd "$(dirname "$0")"

. ./27_galaxies.sh

DEFAULT_REMOTE_USER="rhuang"
REMOTE_HOST="setonix.pawsey.org.au"

if is_known_galid "$1" || looks_like_galid "$1"; then
  REMOTE_USER="$DEFAULT_REMOTE_USER"
else
  REMOTE_USER="$1"
  shift

  if [[ $# -gt 0 ]]; then
    if is_known_galid "$1" || looks_like_galid "$1"; then
      :
    else
      REMOTE_HOST="$1"
      shift
    fi
  fi
fi

select_galids "$@"
GALIDS=("${SELECTED_GALIDS[@]}")

if [[ "$REMOTE_USER" == *@* ]]; then
  REMOTE_LOGIN="$REMOTE_USER"
else
  REMOTE_LOGIN="${REMOTE_USER}@${REMOTE_HOST}"
fi

REMOTE_RUN_DIR="/software/projects/pawsey1308/ngist_supplementary_public/ngistTutorial"
REMOTE_CONFIG_DIR="${REMOTE_RUN_DIR}/configFiles"

yaml_files=()
slurm_files=()
send_galids=()
skipped_galids=()

for galid in "${GALIDS[@]}"; do
  yaml_file="${galid}_MAUVE_MasterConfig_v7.6.8_setonix.yaml"
  yaml_file_7000="${galid}_MAUVE_MasterConfig_v7.6.8_7000_setonix.yaml"
  slurm_file="${galid}_v3tk_v7.6.8_setonix.slurm"
  slurm_file_7000="${galid}_v3tk_v7.6.8_7000_setonix.slurm"
  missing=()

  if [[ ! -f "$yaml_file" ]]; then
    missing+=("$yaml_file")
  fi

  if [[ ! -f "$yaml_file_7000" ]]; then
    missing+=("$yaml_file_7000")
  fi

  if [[ ! -f "$slurm_file" ]]; then
    missing+=("$slurm_file")
  fi

  if [[ ! -f "$slurm_file_7000" ]]; then
    missing+=("$slurm_file_7000")
  fi

  if [[ "${#missing[@]}" -gt 0 ]]; then
    echo "WARNING: skipping ${galid}: missing generated file(s): ${missing[*]}" >&2
    skipped_galids+=("$galid")
    continue
  fi

  yaml_files+=("$yaml_file" "$yaml_file_7000")
  slurm_files+=("$slurm_file" "$slurm_file_7000")
  send_galids+=("$galid")
done

if [[ "${#send_galids[@]}" -eq 0 ]]; then
  echo "WARNING: no generated YAML/slurm file pairs available to send; nothing to do." >&2
  exit 0
fi

control_files=(
  27_galaxies.sh
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
echo "Selected ${#send_galids[@]} generated run ID(s) to send; skipped ${#skipped_galids[@]} unavailable run ID(s)."
echo "Creating remote directories if needed"
ssh "$REMOTE_LOGIN" "mkdir -p '$REMOTE_CONFIG_DIR' '$REMOTE_RUN_DIR'"

echo "Sending ${#yaml_files[@]} YAML files to ${REMOTE_CONFIG_DIR}"
COPYFILE_DISABLE=1 tar -cf - "${yaml_files[@]}" | ssh "$REMOTE_LOGIN" \
  "tar -xf - -C '$REMOTE_CONFIG_DIR' && find '$REMOTE_CONFIG_DIR' -maxdepth 1 -type f -name '._*' -delete"

echo "Sending ${#slurm_files[@]} slurm scripts plus control scripts to ${REMOTE_RUN_DIR}"
COPYFILE_DISABLE=1 tar -cf - "${slurm_files[@]}" "${control_files[@]}" | ssh "$REMOTE_LOGIN" \
  "tar -xf - -C '$REMOTE_RUN_DIR' && chmod +x '$REMOTE_RUN_DIR'/27_setonix.sh '$REMOTE_RUN_DIR'/27_status.sh '$REMOTE_RUN_DIR'/*_v3tk_v7.6.8*_setonix.slurm && find '$REMOTE_RUN_DIR' -maxdepth 1 -type f -name '._*' -delete"

echo "Done."
