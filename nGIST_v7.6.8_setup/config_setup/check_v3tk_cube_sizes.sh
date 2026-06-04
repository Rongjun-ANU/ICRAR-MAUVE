#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./check_v3tk_cube_sizes.sh
#   ./check_v3tk_cube_sizes.sh NGC4254
#   ./check_v3tk_cube_sizes.sh NGC4254 NGC4321 NGC4535
#   ./check_v3tk_cube_sizes.sh IC3392,NGC4254,NGC4567_8
#   ./check_v3tk_cube_sizes.sh [target_dir] [glob_patterns] [output_csv]
#
# Defaults:
#   target_dir    = /scratch/pawsey1308/mauve/cubes/v3tk
#   glob_patterns = *_DATACUBE*.fits,*_DATACUBE*.fits.gz
#   output_csv    = ./cube_sizes_v3tk.csv
#
# Notes:
#   - No nGIST overlay is used.
#   - No astropy/numpy is required.
#   - Reads FITS headers directly.
#   - available_spaxels is counted from the central wavelength layer.
#   - .part files are ignored.
#   - Both .fits and .fits.gz are kept, because this is a size table.

DEFAULT_TARGET_DIR="/scratch/pawsey1308/mauve/cubes/v3tk"
DEFAULT_GLOB_PATTERNS="*_DATACUBE*.fits,*_DATACUBE*.fits.gz"
DEFAULT_OUTPUT_CSV="cube_sizes_v3tk.csv"

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
import struct
import sys
from pathlib import Path


BLOCK_SIZE = 2880

FIELDNAMES = [
    "ID",
    "file",
    "shape",
    "size_B",
    "size_GB",
    "total_spaxels",
    "available_spaxels",
    "status",
]


def product(values):
    out = 1
    for v in values:
        out *= int(v)
    return out


def cube_id_from_filename(file_name):
    """
    Handles:
      IC3392_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk.fits
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


def galaxy_sort_key_from_id(galaxy_id):
    galaxy_id = str(galaxy_id).upper()
    match = re.match(r"^(IC|NGC)(\d+)(?:_(\d+))?$", galaxy_id)

    if not match:
        return (10**12, 10**12, galaxy_id)

    prefix = match.group(1)
    number_1 = int(match.group(2))
    number_2 = int(match.group(3)) if match.group(3) else -1

    return (number_1, number_2, prefix)


def row_sort_key(row):
    file_name = Path(str(row.get("file", ""))).name

    if file_name.endswith(".fits"):
        suffix_rank = 0
    elif file_name.endswith(".fits.gz"):
        suffix_rank = 1
    else:
        suffix_rank = 9

    return galaxy_sort_key_from_id(row.get("ID", "")) + (suffix_rank, file_name)


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

    return sorted(files, key=lambda p: row_sort_key({"ID": cube_id_from_filename(p.name), "file": str(p)}))


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

    return sorted(files, key=lambda p: row_sort_key({"ID": cube_id_from_filename(p.name), "file": str(p)}))


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


def bytes_per_value(bitpix):
    return abs(int(bitpix)) // 8


def data_size_bytes(header):
    naxis = int(header.get("NAXIS", 0))
    bitpix = int(header.get("BITPIX", 8))

    if naxis == 0:
        raw_size = 0
    else:
        n_elem = 1
        for i in range(1, naxis + 1):
            n_elem *= int(header.get(f"NAXIS{i}", 0))

        pcount = int(header.get("PCOUNT", 0))
        gcount = int(header.get("GCOUNT", 1))

        raw_size = (n_elem * gcount + pcount) * bytes_per_value(bitpix)

    if raw_size == 0:
        return 0

    return ((raw_size + BLOCK_SIZE - 1) // BLOCK_SIZE) * BLOCK_SIZE


def first_image_info(file_path):
    """
    Return:
      hdu_index, shape, header, data_start
    where shape is in numpy/astropy order:
      (nz, ny, nx) for a 3D cube.
    """
    opener = gzip.open if file_path.name.endswith(".gz") else open

    with opener(file_path, "rb") as f:
        hdu_index = 0

        while True:
            header = read_header(f)
            data_start = f.tell()
            naxis = int(header.get("NAXIS", 0))

            if naxis >= 2:
                shape = tuple(int(header[f"NAXIS{i}"]) for i in range(naxis, 0, -1))
                return hdu_index, shape, header, data_start

            skip = data_size_bytes(header)
            if skip:
                f.seek(skip, 1)

            hdu_index += 1


def count_values_from_stream(f, offset, count, bitpix, blank_value):
    value_size = bytes_per_value(bitpix)
    total_bytes = count * value_size

    f.seek(offset)

    # Integer FITS images have no NaN. If BLANK is absent, all values are available.
    if bitpix > 0 and blank_value is None:
        return count

    if bitpix == -32:
        fmt = ">f"
    elif bitpix == -64:
        fmt = ">d"
    elif bitpix == 8:
        fmt = ">B"
    elif bitpix == 16:
        fmt = ">h"
    elif bitpix == 32:
        fmt = ">i"
    elif bitpix == 64:
        fmt = ">q"
    else:
        raise ValueError(f"unsupported BITPIX={bitpix}")

    available = 0
    remaining = total_bytes
    chunk_values = 262144
    chunk_bytes = chunk_values * value_size

    while remaining > 0:
        n_read = min(chunk_bytes, remaining)
        data = f.read(n_read)

        if not data:
            raise EOFError("unexpected EOF while reading central layer")

        # Ensure buffer length is a multiple of one value.
        usable = (len(data) // value_size) * value_size
        data = data[:usable]

        for item in struct.iter_unpack(fmt, data):
            value = item[0]

            if bitpix < 0:
                if math.isfinite(value):
                    available += 1
            else:
                if blank_value is None or value != blank_value:
                    available += 1

        remaining -= usable

    return available


def available_spaxels_central_layer(file_path, shape, header, data_start):
    """
    Mimics the old intent:
      available_spaxels = np.sum(np.isfinite(data[..., center_z, :, :]))

    For normal MUSE cubes with shape=(nz, ny, nx), this reads only the central
    wavelength layer, not the whole cube.
    """
    bitpix = int(header.get("BITPIX", 0))
    blank_value = header.get("BLANK", None)

    ny = int(shape[-2])
    nx = int(shape[-1])
    plane_values = nx * ny
    value_size = bytes_per_value(bitpix)

    opener = gzip.open if file_path.name.endswith(".gz") else open

    with opener(file_path, "rb") as f:
        if len(shape) == 2:
            return count_values_from_stream(
                f=f,
                offset=data_start,
                count=plane_values,
                bitpix=bitpix,
                blank_value=blank_value,
            )

        nz = int(shape[-3])
        center_z = nz // 2

        if len(shape) == 3:
            leading_count = 1
        else:
            leading_count = product(shape[:-3])

        available = 0

        for leading_index in range(leading_count):
            # FITS storage order:
            #   x fastest, then y, then z, then higher axes.
            plane_index = leading_index * nz + center_z
            offset = data_start + plane_index * plane_values * value_size

            available += count_values_from_stream(
                f=f,
                offset=offset,
                count=plane_values,
                bitpix=bitpix,
                blank_value=blank_value,
            )

        return available


def make_row(file_path):
    file_name = file_path.name
    cube_id = cube_id_from_filename(file_name)
    size_B = file_path.stat().st_size
    size_GB = size_B / float(1024 ** 3)

    try:
        _, shape, header, data_start = first_image_info(file_path)

        if len(shape) < 2:
            raise ValueError(f"expected at least 2D data, got shape={shape}")

        ny = int(shape[-2])
        nx = int(shape[-1])
        total_spaxels = nx * ny

        available_spaxels = available_spaxels_central_layer(
            file_path=file_path,
            shape=shape,
            header=header,
            data_start=data_start,
        )

        return {
            "ID": cube_id,
            "file": str(file_path),
            "shape": str(shape),
            "size_B": str(size_B),
            "size_GB": f"{size_GB:.3f}",
            "total_spaxels": str(total_spaxels),
            "available_spaxels": str(available_spaxels),
            "status": "ok",
        }

    except Exception as exc:
        return {
            "ID": cube_id,
            "file": str(file_path),
            "shape": "",
            "size_B": str(size_B),
            "size_GB": f"{size_GB:.3f}",
            "total_spaxels": "",
            "available_spaxels": "",
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
            f"size={row['size_GB']}GB "
            f"shape={row['shape']} "
            f"total_spaxels={row['total_spaxels']} "
            f"available_spaxels={row['available_spaxels']} "
            f"status={row['status']}",
            flush=True,
        )

    return 0 if ok_count > 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
PY

python3 "$TMP_PY" "$TARGET_DIR" "$GLOB_PATTERNS" "$OUTPUT_CSV" "${GALAXIES[@]}"
