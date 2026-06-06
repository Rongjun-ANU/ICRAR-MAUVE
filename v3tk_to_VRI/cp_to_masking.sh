#!/usr/bin/env bash
set -euo pipefail

dest="../v3tk_masking_VRI"

usage() {
  printf 'Usage: %s [GALID]\n' "$0" >&2
}

copy_file() {
  local src=$1
  if [[ "${DRY_RUN:-0}" == "1" ]]; then
    printf 'cp -p %q %q\n' "./$src" "$dest/"
  else
    cp -p "./$src" "$dest/"
  fi
}

copy_galid() {
  local galid=$1
  local png="${galid}_combined_VRI.png"
  local -a fits_matches

  if [[ ! -f "$png" ]]; then
    printf 'Missing required file: %s\n' "$png" >&2
    return 1
  fi

  shopt -s nullglob
  fits_matches=("${galid}"*_VRI.fits)
  shopt -u nullglob

  if (( ${#fits_matches[@]} == 0 )); then
    printf 'Missing required file matching: %s*_VRI.fits\n' "$galid" >&2
    return 1
  fi
  if (( ${#fits_matches[@]} > 1 )); then
    printf 'Multiple FITS matches for %s:\n' "$galid" >&2
    printf '  %s\n' "${fits_matches[@]}" >&2
    return 1
  fi

  copy_file "${fits_matches[0]}"
  copy_file "$png"
}

if (( $# > 1 )); then
  usage
  exit 2
fi

if [[ ! -d "$dest" ]]; then
  printf 'Missing destination directory: %s\n' "$dest" >&2
  exit 1
fi

if (( $# == 1 )); then
  copy_galid "$1"
  exit
fi

shopt -s nullglob
combined_pngs=(*_combined_VRI.png)
shopt -u nullglob

if (( ${#combined_pngs[@]} == 0 )); then
  printf 'No *_combined_VRI.png files found.\n' >&2
  exit 1
fi

for png in "${combined_pngs[@]}"; do
  galid=${png%_combined_VRI.png}
  [[ "$galid" == "ALL" || "$galid" == "All" ]] && continue
  copy_galid "$galid"
done
