from __future__ import annotations

from pathlib import Path
import re
import struct


_FITS_BLOCK_SIZE = 2880
_ROW_FORMAT = ">7sffdffh"


def _read_header(raw: bytes, start: int) -> tuple[list[str], int]:
    cards: list[str] = []
    pos = start
    while True:
        block = raw[pos : pos + _FITS_BLOCK_SIZE]
        if not block:
            raise ValueError("Unexpected end of FITS file while reading header")
        for i in range(0, _FITS_BLOCK_SIZE, 80):
            card = block[i : i + 80].decode("ascii", errors="ignore")
            cards.append(card)
            if card.startswith("END"):
                return cards, pos + _FITS_BLOCK_SIZE
        pos += _FITS_BLOCK_SIZE


def _parse_header_values(cards: list[str]) -> dict[str, str]:
    values: dict[str, str] = {}
    for card in cards:
        if "=" in card[:10]:
            values[card[:8].strip()] = card[10:30].strip()
    return values


def _padded_data_size(header_values: dict[str, str]) -> int:
    naxis = int(header_values.get("NAXIS", "0").split("/")[0])
    if naxis == 0:
        return 0

    bitpix = int(header_values["BITPIX"].split("/")[0])
    size = abs(bitpix) // 8
    for axis in range(1, naxis + 1):
        size *= int(header_values[f"NAXIS{axis}"].split("/")[0])

    return ((size + _FITS_BLOCK_SIZE - 1) // _FITS_BLOCK_SIZE) * _FITS_BLOCK_SIZE


def normalize_galaxy_name(name: str) -> str:
    clean = name.replace("\x00", "")
    return re.sub(r"\s+", "", clean).upper()


def load_mauve_master_properties(table_path: str | Path) -> tuple[dict[str, float], dict[str, float]]:
    """
    Load integrated log stellar masses and log SFRs from mauve_master_wiki.fits.

    The FITS file is a binary table with columns:
      Galaxy, logMstar, logSFR, logMHI, logMH2, DMS, vel
    """
    table_path = Path(table_path)
    if not table_path.exists():
        raise FileNotFoundError(f"Missing MAUVE master sample file: {table_path}")

    raw = table_path.read_bytes()

    primary_cards, primary_end = _read_header(raw, 0)
    primary_values = _parse_header_values(primary_cards)
    ext_start = primary_end + _padded_data_size(primary_values)

    ext_cards, ext_end = _read_header(raw, ext_start)
    ext_values = _parse_header_values(ext_cards)

    row_length = int(ext_values["NAXIS1"].split("/")[0])
    row_count = int(ext_values["NAXIS2"].split("/")[0])
    expected_length = struct.calcsize(_ROW_FORMAT)
    if row_length != expected_length:
        raise ValueError(
            f"Unexpected row length in {table_path}: got {row_length}, expected {expected_length}"
        )

    data_start = ext_end
    logmstar: dict[str, float] = {}
    logsfr: dict[str, float] = {}

    for idx in range(row_count):
        row = raw[data_start + idx * row_length : data_start + (idx + 1) * row_length]
        if len(row) != row_length:
            raise ValueError(f"Truncated FITS row {idx} in {table_path}")

        name_bytes, mass, sfr, *_ = struct.unpack(_ROW_FORMAT, row)
        name = normalize_galaxy_name(name_bytes.decode("ascii", errors="ignore"))
        if not name:
            continue

        if mass == mass:
            logmstar[name] = float(mass)
        if sfr == sfr:
            logsfr[name] = float(sfr)

    return logmstar, logsfr
