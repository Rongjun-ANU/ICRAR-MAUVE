#!/usr/bin/env python3
"""Hash actual CANFAR product files and write a small SHA-256 manifest."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import sys
import time
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


def hash_download_request(args):
    request = Path(args.request)
    root = request.parent.parent.resolve()
    paths = json.loads(request.read_text(encoding="utf-8"))

    def hash_product(relative):
        try:
            path = (root / relative).resolve()
            path.relative_to(root)
            return hash_file(path)
        except Exception as error:
            return {"status": "error", "error": str(error)}

    with ThreadPoolExecutor(max_workers=4) as pool:
        entries = dict(zip(paths, pool.map(hash_product, paths)))
    response = request.with_suffix(".response.json")
    temporary = response.with_suffix(".tmp")
    temporary.write_text(json.dumps(entries), encoding="utf-8")
    os.replace(temporary, response)
    return 0


def fetch_download_checksums(directory, remote_root, paths):
    import requests
    from cadcutils.exceptions import NotFoundException
    from vos import Client

    cert = str(Path.home() / ".ssl" / "cadcproxy.pem")
    api = os.environ.get("CANFAR_API", "https://ws-uv.canfar.net/skaha/v1").rstrip("/")
    if not api.startswith("https://"):
        raise ValueError("CANFAR_API must use HTTPS")
    image = os.environ.get("CANFAR_CHECKSUM_IMAGE", "images.canfar.net/skaha/astroml:latest")
    timeout = int(os.environ.get("CHECKSUM_WAIT_SECONDS", "7200"))
    if timeout <= 0:
        raise ValueError("CHECKSUM_WAIT_SECONDS must be positive")
    request_id = uuid.uuid4().hex
    name = ".download_checksum_" + request_id
    request = directory / (name + ".json")
    response = request.with_suffix(".response.json")
    request.write_text(json.dumps(paths), encoding="utf-8")
    remote_base = remote_root + "/" + paths[0].split("/")[0] + "/"
    worker_node = remote_base + name + ".py"
    request_node = remote_base + request.name
    response_node = remote_base + response.name
    worker_path = "/arc/" + worker_node.removeprefix("arc:")
    request_path = "/arc/" + request_node.removeprefix("arc:")

    client = Client()
    with requests.Session() as session:
        session.cert = cert
        job_url = None
        try:
            client.copy(str(Path(__file__).resolve()), worker_node)
            client.copy(str(request), request_node)
            started = session.post(api + "/session", data={
                "name": "dl-checksum-" + request_id[:12],
                "image": image, "type": "headless", "cores": 4, "ram": 2,
                "cmd": "/usr/bin/env",
                "args": shlex.join(["python3", worker_path, "hash-download-request", request_path]),
            }, timeout=(15, 60))
            started.raise_for_status()
            job_id = started.text.strip()
            if not re.fullmatch(r"[A-Za-z0-9_-]+", job_id):
                raise ValueError("Unexpected CANFAR checksum job ID")
            job_url = api + "/session/" + job_id
            print("Started CANFAR checksum job " + job_id, file=sys.stderr)
            deadline = time.monotonic() + timeout
            finished_at = None
            while time.monotonic() < deadline:
                try:
                    client.copy(response_node, str(response))
                except NotFoundException:
                    pass
                else:
                    return json.loads(response.read_text(encoding="utf-8"))
                reply = session.get(job_url, timeout=(10, 30))
                reply.raise_for_status()
                state = str(reply.json().get("status", "Unknown")).lower()
                if state in {"failed", "terminated", "error"}:
                    raise RuntimeError("CANFAR checksum job ended with status " + state)
                if state in {"completed", "succeeded"}:
                    finished_at = finished_at or time.monotonic()
                    if time.monotonic() - finished_at >= 30:
                        raise RuntimeError("CANFAR checksum job finished without a response")
                time.sleep(5)
            raise TimeoutError("CANFAR checksum job timed out after {} seconds".format(timeout))
        finally:
            if job_url:
                try:
                    deleted = session.delete(job_url, timeout=(10, 30))
                    if deleted.status_code != 404:
                        deleted.raise_for_status()
                except requests.RequestException as error:
                    print("WARNING: could not delete CANFAR checksum job " + job_id + ": " + str(error),
                          file=sys.stderr)
            for node in (worker_node, request_node, response_node):
                try:
                    client.delete(node)
                except NotFoundException:
                    pass
                except Exception:
                    print("WARNING: could not remove " + node, file=sys.stderr)


def download_queue(args):
    directory = Path(args.directory)
    prefix = os.fsencode(args.remote_root + "/")
    pending = {}
    for record in sorted(directory.glob("pending.*")):
        remote, local, end = record.read_bytes().split(b"\0")
        if end or not remote.startswith(prefix):
            raise ValueError("Invalid pending download record")
        relative = os.fsdecode(remote[len(prefix):])
        if not re.fullmatch(r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+", relative) or \
                any(part in {".", ".."} for part in relative.split("/")):
            raise ValueError("Invalid CANFAR product path " + relative)
        pending[relative] = (remote, local)
    if not pending:
        return 0
    checksums = fetch_download_checksums(directory, args.remote_root, list(pending))
    failed = False
    for relative, (remote, local) in pending.items():
        entry = checksums.get(relative, {})
        digest = entry.get("sha256")
        if entry.get("status") == "ok" and isinstance(digest, str) and re.fullmatch(r"[0-9a-f]{64}", digest):
            sys.stdout.buffer.write(b"\0".join((remote, local, digest.encode("ascii"), b"")))
        else:
            print("ERROR: CANFAR SHA-256 unavailable for " + relative + ": " +
                  entry.get("error", entry.get("status", "missing response entry")), file=sys.stderr)
            failed = True
    return int(failed)


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
    processor = sub.add_parser("hash-download-request")
    processor.add_argument("request")
    downloader = sub.add_parser("download-queue")
    downloader.add_argument("directory")
    downloader.add_argument("remote_root")
    args = parser.parse_args()
    try:
        actions = {"generate": generate, "lookup": lookup,
                   "hash-download-request": hash_download_request, "download-queue": download_queue}
        return actions[args.command](args)
    except KeyboardInterrupt:
        return 130
    except Exception as error:
        print("ERROR:", error, file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
