#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 [USERNAME] [HOST]" >&2
  echo "Example: $0" >&2
  echo "Example with explicit user: $0 rhuang" >&2
  echo "Example with explicit user and host: $0 rhuang setonix.pawsey.org.au" >&2
}

if [[ $# -gt 2 ]]; then
  usage
  exit 2
fi

cd "$(dirname "$0")"

DEFAULT_REMOTE_USER="rhuang"
DEFAULT_REMOTE_HOST="setonix.pawsey.org.au"
REMOTE_USER="${1:-$DEFAULT_REMOTE_USER}"
REMOTE_HOST="${2:-$DEFAULT_REMOTE_HOST}"
REMOTE_DIR="/scratch/pawsey1308/mauve/cubes/v3tk"

if [[ "$REMOTE_USER" == *@* ]]; then
  REMOTE_LOGIN="$REMOTE_USER"
else
  REMOTE_LOGIN="${REMOTE_USER}@${REMOTE_HOST}"
fi

shopt -s nullglob
mask_files=()
for mask_file in *_mask.fits; do
  if [[ -f "$mask_file" ]]; then
    mask_files+=("$mask_file")
  fi
done
shopt -u nullglob

if [[ "${#mask_files[@]}" -eq 0 ]]; then
  echo "ERROR: no *_mask.fits files found in $(pwd)" >&2
  exit 1
fi

echo "Remote login: ${REMOTE_LOGIN}"
echo "Uploading ${#mask_files[@]} mask file(s) to ${REMOTE_DIR}"
echo "Existing remote files with the same names will be overwritten."

ssh "$REMOTE_LOGIN" "mkdir -p '$REMOTE_DIR'"
COPYFILE_DISABLE=1 tar -cf - -- "${mask_files[@]}" | ssh "$REMOTE_LOGIN" \
  "tar -xf - -C '$REMOTE_DIR' && find '$REMOTE_DIR' -maxdepth 1 -type f -name '._*' -delete"

echo "Done."
