#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./check_v3tk_cube_centers.sh
#   ./check_v3tk_cube_centers.sh NGC4254
#   ./check_v3tk_cube_centers.sh NGC4254 NGC4321 NGC4535
#   ./check_v3tk_cube_centers.sh IC3392,NGC4254,NGC4567_8
#   ./check_v3tk_cube_centers.sh [target_dir] [glob_patterns] [output_csv]
#
# Defaults:
#   target_dir    = /scratch/pawsey1308/mauve/cubes/v3tk
#   glob_patterns = *_DATACUBE*.fits,*_DATACUBE*.fits.gz
#   output_csv    = ./cube_centers_v3tk.csv

DEFAULT_TARGET_DIR="/scratch/pawsey1308/mauve/cubes/v3tk"
DEFAULT_GLOB_PATTERNS="*_DATACUBE*.fits,*_DATACUBE*.fits.gz"
DEFAULT_OUTPUT_CSV="cube_centers_v3tk.csv"

is_galaxy_token() {
    # Supports:
    #   IC3392
    #   NGC4254
    #   NGC4567_8
    [[ "$1" =~ ^([Nn][Gg][Cc]|[Ii][Cc])[0-9]+(_[0-9]+)?$ ]]
}

looks_like_galaxy_list() {
    local token="$1"
    IFS=',' read -ra parts <<< "$token"
    for part in "${parts[@]}"; do
        part="$(echo "$part" | xargs)"
        if ! is_galaxy_token "$part"; then
            return 1
        fi
    done
    return 0
}

TARGET_DIR="$DEFAULT_TARGET_DIR"
GLOB_PATTERNS="$DEFAULT_GLOB_PATTERNS"
OUTPUT_CSV="$DEFAULT_OUTPUT_CSV"
GALAXIES=()

if [[ $# -eq 0 ]]; then
    :
elif looks_like_galaxy_list "$1"; then
    for arg in "$@"; do
        IFS=',' read -ra parts <<< "$arg"
        for part in "${parts[@]}"; do
            part="$(echo "$part" | xargs | tr '[:lower:]' '[:upper:]')"
            [[ -n "$part" ]] && GALAXIES+=("$part")
        done
    done
else
    TARGET_DIR="${1:-$DEFAULT_TARGET_DIR}"
    GLOB_PATTERNS="${2:-$DEFAULT_GLOB_PATTERNS}"
    OUTPUT_CSV="${3:-$DEFAULT_OUTPUT_CSV}"
fi

if [[ ! -d "$TARGET_DIR" ]]; then
    echo "Error: target directory does not exist: $TARGET_DIR" >&2
    exit 1
fi

TMP_PY="$(mktemp)"
trap 'rm -f "$TMP_PY"' EXIT

cat > "$TMP_PY" <<'PY'
import csv
import gzip
import math
import re
import sys
from pathlib import Path


BLOCK_SIZE = 2880

FIELDNAMES = [
    "ID",
    "file",
    "shape",
    "nz",
    "ny",
    "nx",
    "center_y",
    "center_x",
    "status",
]


def round_half_up(value):
    return int(math.floor(value + 0.5))


def cube_id_from_filename(file_name):
    """
    Handles both:
      NGC4064_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk.fits
      NGC4254_PHANGS_DATACUBE_native.fits
      NGC4567_8_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk.fits
    """
    cube_id = file_name.split("_DATACUBE", 1)[0]

    if cube_id.endswith("_PHANGS"):
        cube_id = cube_id[:-len("_PHANGS")]

    return cube_id.upper()


def is_complete_fits(file_path):
    name = file_path.name
    return name.endswith(".fits") or name.endswith(".fits.gz")


def file_preference_key(file_path):
    """
    Lower is better.
    Prefer uncompressed .fits over .fits.gz.
    Ignore .part files elsewhere.
    """
    name = file_path.name

    if name.endswith(".fits"):
        suffix_rank = 0
    elif name.endswith(".fits.gz"):
        suffix_rank = 1
    else:
        suffix_rank = 9

    # Prefer ordinary MAUVE cube over PHANGS only if both exist for same ID.
    # In your current directory this should not matter.
    phangs_rank = 1 if "_PHANGS_DATACUBE" in name else 0

    return (suffix_rank, phangs_rank, name)


def deduplicate_by_cube_id(files):
    best = {}

    for file_path in files:
        if not is_complete_fits(file_path):
            continue

        cube_id = cube_id_from_filename(file_path.name)

        if cube_id not in best:
            best[cube_id] = file_path
        else:
            if file_preference_key(file_path) < file_preference_key(best[cube_id]):
                best[cube_id] = file_path

    return [best[key] for key in sorted(best.keys(), key=galaxy_sort_key_from_id)]


def galaxy_sort_key_from_id(galaxy_id):
    galaxy_id = str(galaxy_id).upper()
    match = re.match(r"^(IC|NGC)(\d+)(?:_(\d+))?$", galaxy_id)

    if not match:
        return (10**12, galaxy_id)

    prefix = match.group(1)
    number_1 = int(match.group(2))
    number_2 = int(match.group(3)) if match.group(3) else -1

    # Numeric order first, then prefix.
    # Example:
    #   IC3392 before NGC4064 because 3392 < 4064.
    #   NGC4567_8 sorted at 4567 with suffix 8.
    return (number_1, number_2, prefix)


def expand_glob_patterns(target_dir, pattern_arg):
    files = []
    seen = set()

    for part in pattern_arg.split(","):
        pattern = part.strip()
        if not pattern:
            continue

        for file_path in sorted(target_dir.glob(pattern)):
            if not is_complete_fits(file_path):
                continue

            key = str(file_path)
            if key not in seen:
                seen.add(key)
                files.append(file_path)

    return deduplicate_by_cube_id(files)


def expand_galaxy_files(target_dir, galaxies):
    files = []
    seen = set()

    for gal in galaxies:
        gal = gal.upper()

        patterns = [
            f"{gal}_DATACUBE*.fits",
            f"{gal}_DATACUBE*.fits.gz",
            f"{gal}_PHANGS_DATACUBE*.fits",
            f"{gal}_PHANGS_DATACUBE*.fits.gz",
        ]

        matches = []

        for pattern in patterns:
            for file_path in sorted(target_dir.glob(pattern)):
                if is_complete_fits(file_path):
                    matches.append(file_path)

        if not matches:
            print(
                f"Warning: no complete .fits/.fits.gz cube found for {gal} in {target_dir}",
                file=sys.stderr,
                flush=True,
            )

        for file_path in matches:
            key = str(file_path)
            if key not in seen:
                seen.add(key)
                files.append(file_path)

    return deduplicate_by_cube_id(files)


def parse_card_value(raw):
    if "/" in raw:
        raw = raw.split("/", 1)[0]
    raw = raw.strip()

    if raw.startswith("'"):
        return raw.strip("' ")

    if raw in ("T", "F"):
        return raw == "T"

    try:
        return int(raw)
    except ValueError:
        try:
            return float(raw)
        except ValueError:
            return raw


def read_header(f):
    header = {}

    while True:
        block = f.read(BLOCK_SIZE)
        if not block:
            raise EOFError("unexpected EOF before END card")

        text = block.decode("ascii", errors="replace")

        for i in range(0, BLOCK_SIZE, 80):
            card = text[i:i + 80]
            key = card[:8].strip()

            if key == "END":
                return header

            if not key:
                continue

            if len(card) > 10 and card[8] == "=":
                header[key] = parse_card_value(card[10:])


def data_size_bytes(header):
    naxis = int(header.get("NAXIS", 0))
    bitpix = abs(int(header.get("BITPIX", 8)))

    if naxis == 0:
        raw_size = 0
    else:
        n_elem = 1
        for i in range(1, naxis + 1):
            n_elem *= int(header.get(f"NAXIS{i}", 0))

        pcount = int(header.get("PCOUNT", 0))
        gcount = int(header.get("GCOUNT", 1))

        raw_size = (n_elem * gcount + pcount) * (bitpix // 8)

    if raw_size == 0:
        return 0

    return ((raw_size + BLOCK_SIZE - 1) // BLOCK_SIZE) * BLOCK_SIZE


def first_image_shape(file_path):
    opener = gzip.open if file_path.name.endswith(".gz") else open

    with opener(file_path, "rb") as f:
        hdu_index = 0

        while True:
            header = read_header(f)
            naxis = int(header.get("NAXIS", 0))

            if naxis >= 2:
                shape = tuple(int(header[f"NAXIS{i}"]) for i in range(naxis, 0, -1))
                return hdu_index, shape

            skip = data_size_bytes(header)
            if skip:
                f.seek(skip, 1)

            hdu_index += 1


def make_row(file_path):
    file_name = file_path.name
    cube_id = cube_id_from_filename(file_name)

    try:
        _, shape = first_image_shape(file_path)

        if len(shape) < 2:
            raise ValueError(f"expected at least 2D data, got shape={shape}")

        ny = int(shape[-2])
        nx = int(shape[-1])
        nz = int(shape[-3]) if len(shape) >= 3 else 1

        center_x = round_half_up(nx / 2.0)
        center_y = round_half_up(ny / 2.0)

        return {
            "ID": cube_id,
            "file": str(file_path),
            "shape": str(shape),
            "nz": nz,
            "ny": ny,
            "nx": nx,
            "center_y": center_y,
            "center_x": center_x,
            "status": "ok",
        }

    except Exception as exc:
        return {
            "ID": cube_id,
            "file": str(file_path),
            "shape": "",
            "nz": "",
            "ny": "",
            "nx": "",
            "center_y": "",
            "center_x": "",
            "status": f"error: {exc}",
        }


def read_existing_csv(output_csv):
    if not output_csv.exists():
        return []

    rows = []

    with output_csv.open("r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            clean = {name: row.get(name, "") for name in FIELDNAMES}
            rows.append(clean)

    return rows


def row_sort_key(row):
    return galaxy_sort_key_from_id(row.get("ID", ""))


def main():
    target_dir = Path(sys.argv[1])
    glob_patterns = sys.argv[2]
    output_csv = Path(sys.argv[3])
    galaxies = [g.upper() for g in sys.argv[4:]]

    specific_mode = bool(galaxies)

    if specific_mode:
        files = expand_galaxy_files(target_dir, galaxies)
    else:
        files = expand_glob_patterns(target_dir, glob_patterns)

    if not files:
        print("No matching complete cube files found.", file=sys.stderr, flush=True)
        return 1

    print(f"Found {len(files)} matching cube file(s).", flush=True)

    new_rows = []

    for i, file_path in enumerate(files, start=1):
        print(f"[{i}/{len(files)}] Checking {file_path.name}", flush=True)
        new_rows.append(make_row(file_path))

    if specific_mode:
        existing_rows = read_existing_csv(output_csv)
        requested_ids = set(row["ID"].upper() for row in new_rows)

        kept_rows = [
            row for row in existing_rows
            if str(row.get("ID", "")).upper() not in requested_ids
        ]

        rows = kept_rows + new_rows
    else:
        rows = new_rows

    rows = sorted(rows, key=row_sort_key)

    output_csv.parent.mkdir(parents=True, exist_ok=True)

    with output_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    ok_count = sum(1 for r in new_rows if r["status"] == "ok")
    err_count = len(new_rows) - ok_count

    print("", flush=True)

    if specific_mode:
        print(f"Updated {output_csv} with requested galaxies: {', '.join(galaxies)}", flush=True)
    else:
        print(f"Wrote full table to {output_csv}", flush=True)

    print(f"New checked rows: {len(new_rows)}", flush=True)
    print(f"OK: {ok_count}, Errors: {err_count}", flush=True)

    print("", flush=True)
    print("New/updated rows:", flush=True)

    for row in new_rows:
        print(
            f"- {Path(row['file']).name}: "
            f"ID={row['ID']} "
            f"shape={row['shape']} "
            f"ny={row['ny']} nx={row['nx']} nz={row['nz']} "
            f"center=(y={row['center_y']}, x={row['center_x']}) "
            f"status={row['status']}",
            flush=True,
        )

    return 0 if ok_count > 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
PY

python3 "$TMP_PY" "$TARGET_DIR" "$GLOB_PATTERNS" "$OUTPUT_CSV" "${GALAXIES[@]}"
