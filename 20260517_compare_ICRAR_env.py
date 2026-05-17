#!/usr/bin/env python3
"""Compare a Setonix ICRAR_Rongjun conda list against the local ICRAR reference.

Usage:
    python 20260517_compare_ICRAR_env.py \
      20260517_ICRAR_local_conda_list_reference.json \
      ICRAR_Rongjun_setonix_conda_list.json

The comparison is intentionally strict for Python-facing packages and Pip
packages, but tolerant for low-level platform packages because macOS arm64 and
Setonix linux-64 cannot have identical binary/runtime packages.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path


MAC_ONLY = {
    "appnope",
    "libcxx",
    "pyobjc-core",
    "pyobjc-framework-cocoa",
    "python_abi",
}

PINNED_BINARY_PYTHON_PACKAGES = {
    "astropy",
    "astropy-base",
    "astropy-healpix",
    "astroquery",
    "bokeh",
    "bqplot",
    "dask",
    "dask-core",
    "distributed",
    "h5py",
    "mamba",
    "matplotlib",
    "matplotlib-base",
    "numpy",
    "opencv",
    "pandas",
    "pillow",
    "pip",
    "py-opencv",
    "pyarrow",
    "pyarrow-core",
    "pyerfa",
    "python",
    "pyvo",
    "radio-beam",
    "reproject",
    "scikit-learn",
    "scipy",
    "shapely",
    "speclite",
    "spectral-cube",
    "synphot",
    "zarr",
}


def load_rows(path: str) -> list[dict]:
    return json.loads(Path(path).read_text())


def by_name(rows: list[dict]) -> dict[str, dict]:
    return {row["name"].lower(): row for row in rows}


def is_python_facing(row: dict) -> bool:
    name = row["name"].lower()
    if name in MAC_ONLY:
        return False
    if row.get("channel") == "pypi":
        return True
    build = row.get("build_string", "")
    return (
        row.get("platform") == "noarch"
        or build.startswith("py")
        or name in PINNED_BINARY_PYTHON_PACKAGES
    )


def main() -> int:
    if len(sys.argv) != 3:
        print(__doc__.strip(), file=sys.stderr)
        return 2

    local_rows = load_rows(sys.argv[1])
    remote_rows = load_rows(sys.argv[2])
    remote = by_name(remote_rows)

    checked = [row for row in local_rows if is_python_facing(row)]
    missing: list[str] = []
    version_mismatch: list[str] = []
    source_mismatch: list[str] = []

    for local in checked:
        name = local["name"].lower()
        got = remote.get(name)
        if got is None:
            missing.append(f"{local['name']}=={local['version']}")
            continue
        if got.get("version") != local.get("version"):
            version_mismatch.append(
                f"{local['name']}: local {local.get('version')} != Setonix {got.get('version')}"
            )
        if local.get("channel") == "pypi" and got.get("channel") != "pypi":
            source_mismatch.append(
                f"{local['name']}: local source pypi != Setonix source {got.get('channel')}"
            )

    print(f"Checked Python-facing/Pip packages: {len(checked)}")
    print(f"Local total packages: {len(local_rows)}")
    print(f"Setonix total packages: {len(remote_rows)}")

    if missing:
        print("\nMissing packages:")
        print("\n".join(f"- {item}" for item in missing))
    if version_mismatch:
        print("\nVersion mismatches:")
        print("\n".join(f"- {item}" for item in version_mismatch))
    if source_mismatch:
        print("\nPip source mismatches:")
        print("\n".join(f"- {item}" for item in source_mismatch))

    if missing or version_mismatch or source_mismatch:
        return 1

    print("\nOK: Python-facing Conda packages and Pip packages match the local reference.")
    print("Low-level system libraries were not required to match because macOS and linux-64 differ.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
