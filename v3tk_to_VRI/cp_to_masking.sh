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
  local copied=0

  shopt -s nullglob
  fits_matches=("${galid}"*_VRI.fits)
  shopt -u nullglob

  if (( ${#fits_matches[@]} > 1 )); then
    printf 'Multiple FITS matches for %s:\n' "$galid" >&2
    printf '  %s\n' "${fits_matches[@]}" >&2
    return 1
  fi

  if (( ${#fits_matches[@]} == 1 )); then
    copy_file "${fits_matches[0]}"
    copied=1
  else
    printf 'Warning: missing file matching: %s*_VRI.fits\n' "$galid" >&2
  fi

  if [[ -f "$png" ]]; then
    copy_file "$png"
    copied=1
  else
    printf 'Warning: missing file: %s\n' "$png" >&2
  fi

  if (( copied == 0 )); then
    return 1
  fi
}

list_galids() {
  local item galid

  shopt -s nullglob
  for item in *_combined_VRI.png; do
    galid=${item%_combined_VRI.png}
    [[ "$galid" == "ALL" || "$galid" == "All" ]] && continue
    printf '%s\n' "$galid"
  done

  for item in *_VRI.fits; do
    case "$item" in
      *_PHANGS_DATACUBE*_VRI.fits)
        galid=${item%%_PHANGS_DATACUBE*}
        ;;
      *_DATACUBE*_VRI.fits)
        galid=${item%%_DATACUBE*}
        ;;
      *)
        continue
        ;;
    esac
    printf '%s\n' "$galid"
  done
  shopt -u nullglob
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

galids=()
while IFS= read -r galid; do
  galids+=("$galid")
done < <(list_galids | sort -u)

if (( ${#galids[@]} == 0 )); then
  printf 'No *_combined_VRI.png or *_VRI.fits files found.\n' >&2
  exit 1
fi

for galid in "${galids[@]}"; do
  copy_galid "$galid"
done
