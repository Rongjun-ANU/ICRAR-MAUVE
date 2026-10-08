#!/usr/bin/env python3
"""Hash actual CANFAR product files and write a small SHA-256 manifest."""

import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import uuid

RUNS = {"normal": "v3tk_v7.6.8", "7000": "v3tk_v7.6.8_7000"}
SUFFIXES = ("CONFIG", "_sfh_maps.fits", "_gas_bin_maps.fits", "_sfh_weights.fits",
            "_gas_spaxel_maps.fits", "_spatial_binning_maps.fits", "_kin_maps.fits",
            "_cont_cube.fits", "LOGFILE", "_mask.fits")
FORMAT = "mauve-sha256-v1"


def hash_file(path):
    if not path.exists():
        return {"status": "missing", "size": None, "sha256": None}
    if not path.is_file():
        raise ValueError("Product is not a regular file")
    before = path.stat()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns, before.st_ino) != (after.st_size, after.st_mtime_ns, after.st_ino):
        raise RuntimeError("File changed while hashing; keep product writers idle")
    return {"status": "ok", "size": after.st_size, "sha256": digest.hexdigest()}


def generate(args):
    root = Path(args.root).resolve()
    run_dir = root / RUNS[args.run]
    if not run_dir.is_dir():
        raise ValueError("Run directory does not exist: " + str(run_dir))
    galaxies = args.galaxies or sorted(p.name for p in run_dir.iterdir()
                                      if p.is_dir() and re.fullmatch(r"(?:IC[0-9]+|NGC[0-9]+(?:_[0-9]+)?)", p.name))
    suffixes = ("_cont_cube.fits",) if args.cont_only else SUFFIXES
    entries = {}
    errors = 0
    for galaxy in galaxies:
        if not re.fullmatch(r"(?:IC[0-9]+|NGC[0-9]+(?:_[0-9]+)?)", galaxy):
            raise ValueError("Invalid galaxy ID: " + galaxy)
        for suffix in suffixes:
            filename = suffix if suffix in ("CONFIG", "LOGFILE") else galaxy + suffix
            key = galaxy + "/" + filename
            path = (run_dir / key).resolve()
            path.relative_to(root)
            try:
                entries[key] = hash_file(path)
            except Exception as error:
                entries[key] = {"status": "error", "error": str(error)}
                errors += 1
            print(key, entries[key]["status"], flush=True)
    manifest = {"format": FORMAT, "manifest_id": uuid.uuid4().hex,
                "run": args.run, "run_subdir": RUNS[args.run],
                "generated_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                "entries": entries}
    output = Path(args.output) if args.output else run_dir / ".checksum_manifest.json"
    temporary = output.with_name(output.name + ".tmp_" + manifest["manifest_id"])
    temporary.write_text(json.dumps(manifest, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    os.replace(str(temporary), str(output))
    print("Saved CANFAR checksum manifest:", output, flush=True)
    return 1 if errors else 0


def lookup(args):
    with Path(args.manifest).open("rb") as stream:
        data = stream.read(4 * 1024 * 1024 + 1)
    if len(data) > 4 * 1024 * 1024:
        raise ValueError("Checksum manifest exceeds 4 MiB")
    manifest = json.loads(data.decode("utf-8"))
    if manifest.get("format") != FORMAT or manifest.get("run") != args.run:
        raise ValueError("Manifest format/run mismatch")
    manifest_id = manifest.get("manifest_id")
    if not isinstance(manifest_id, str) or not re.fullmatch(r"[0-9a-f]{32}", manifest_id):
        raise ValueError("Invalid manifest ID")
    key = args.galaxy + "/" + args.filename
    if key not in manifest.get("entries", {}):
        raise ValueError("File not covered by manifest: " + key + "; regenerate with matching selections")
    entry = manifest["entries"][key]
    if entry.get("status") == "missing":
        print("MISSING -", manifest_id)
    elif entry.get("status") == "ok":
        digest, size = entry.get("sha256"), entry.get("size")
        if not isinstance(digest, str) or not re.fullmatch(r"[0-9a-f]{64}", digest):
            raise ValueError("Invalid SHA-256 in manifest")
        if not isinstance(size, int) or isinstance(size, bool) or size < 0:
            raise ValueError("Invalid byte size in manifest")
        print(size, digest, manifest_id)
    else:
        raise ValueError("CANFAR could not hash " + key + ": " + entry.get("error", "invalid status"))
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    generator = sub.add_parser("generate", help="Run inside CANFAR; reads file bytes, not VOS checksum metadata")
    generator.add_argument("--root", default="/arc/projects/mauve/products")
    generator.add_argument("--cont-only", action="store_true")
    generator.add_argument("--output")
    generator.add_argument("run", choices=RUNS)
    generator.add_argument("galaxies", nargs="*")
    reader = sub.add_parser("lookup", help="Read a manifest entry on Setonix")
    reader.add_argument("--manifest", required=True)
    reader.add_argument("--run", choices=RUNS, required=True)
    reader.add_argument("--galaxy", required=True)
    reader.add_argument("--filename", required=True)
    args = parser.parse_args()
    try:
        return generate(args) if args.command == "generate" else lookup(args)
    except Exception as error:
        print("ERROR:", error, file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
