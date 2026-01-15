#!/usr/bin/env bash
# cp_to_project.sh
set -Eeuo pipefail

SRC="/arc/home/RongjunHuang/ICRAR/v3tk_to_VRI"
DST="/arc/projects/mauve/v3tk_to_VRI"
DRYRUN=0   # use -n to preview without copying

# Parse flags
while getopts "n" opt; do
  case "$opt" in
    n) DRYRUN=1 ;;
  esac
done

# Sanity checks
if [[ ! -d "$SRC" ]]; then
  echo "ERROR: Source directory not found: $SRC" >&2
  exit 1
fi

mkdir -p "$DST"

# Prefer rsync for safe, resumable copy
if command -v rsync >/dev/null 2>&1; then
  RSYNC_OPTS=(-a --human-readable --info=stats2,progress2)
  [[ $DRYRUN -eq 1 ]] && RSYNC_OPTS+=(--dry-run)

  # Include only the targets we want
  rsync "${RSYNC_OPTS[@]}" \
    --include='*' \
    "$SRC"/ "$DST"/
else
  echo "rsync not found; falling back to cp (no progress/dry-run support)."
  # Copy FITS files
  shopt -s nullglob
  mkdir -p "$DST"
  for f in "$SRC"/*; do
    cp -av "$f" "$DST"/
  done
fi

echo "✅ Copy complete to: $DST"
