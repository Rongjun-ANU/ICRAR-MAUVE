#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# Locate inkscape (PATH install or macOS app bundle)
if command -v inkscape >/dev/null 2>&1; then
  INKSCAPE="inkscape"
elif [[ -x "/Applications/Inkscape.app/Contents/MacOS/inkscape" ]]; then
  INKSCAPE="/Applications/Inkscape.app/Contents/MacOS/inkscape"
else
  echo "ERROR: Inkscape not found."
  exit 1
fi

pdfs=( *.pdf )
if (( ${#pdfs[@]} == 0 )); then
  echo "No PDF files found in: $(pwd)"
  exit 0
fi

echo "Using Inkscape: $INKSCAPE"
echo "Converting ${#pdfs[@]} PDF(s) with --pdf-poppler (text -> vector paths)"
echo

for pdf in "${pdfs[@]}"; do
  base="${pdf%.pdf}"
  out="${base}.svg"
  echo "Converting: $pdf -> $out"

  "$INKSCAPE" --pdf-poppler "$pdf" \
    --export-type=svg \
    --export-plain-svg \
    --export-filename="$out" \
    >/dev/null 2>&1

  echo "  OK"
done

echo
echo "Done."
