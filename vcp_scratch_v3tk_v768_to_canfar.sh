#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage: ./vcp_scratch_v3tk_v768_to_canfar.sh [OPTIONS] [normal|7000] [GALID ...]

Upload selected nGIST products directly from Setonix scratch to CANFAR.

Examples:
  ./vcp_scratch_v3tk_v768_to_canfar.sh
      Upload all available galaxies from normal and 7000 runs.
  ./vcp_scratch_v3tk_v768_to_canfar.sh normal
      Upload all available galaxies from the normal run.
  ./vcp_scratch_v3tk_v768_to_canfar.sh 7000
      Upload all available galaxies from the 7000 run.
  ./vcp_scratch_v3tk_v768_to_canfar.sh 7000 NGC4607 NGC4698
      Upload only NGC4607 and NGC4698 from the 7000 run.
  ./vcp_scratch_v3tk_v768_to_canfar.sh --force-overwrite --no-checksum 7000 NGC4216 NGC4689
      Back up existing files and upload fresh copies; destination bytes UNVERIFIED.

Options:
  -n, --dry-run  Show selected galaxies and source files without contacting CADC
  --cont-only    Compare/upload only the continuum cubes
  --canfar-manifest  First option: hash files INSIDE CANFAR using this script.
                    Optional --root overrides /arc/projects/mauve/products.
  --auto-checksum   Launch CANFAR checksum jobs automatically (the default).
  --manual-checksum Use an existing manually generated manifest instead.
  --force-overwrite --no-checksum  Required together: back up each existing file,
                    then upload without CANFAR checksum jobs or manifests.
                    Local source hash/change checks remain enabled.
  -h, --help     Show this help

Environment:
  JOBS              Concurrent galaxy workers (default 5, maximum 5)
  CADC_USER         CADC username for certificate renewal (default RongjunHuang)
  TRANSFER_JOBS     Shared parallel file slots across all galaxies (default 10, maximum 10)
  CADC_READ_ONLY    auto (default), 1 to share a read-only overlay, or 0 to copy overlays
  FILE_RETRIES      Attempts per individual file (default 5)
  RETRY_BASE_SLEEP  Initial retry delay in seconds (default 30)
  RETRY_MAX_SLEEP   Maximum retry delay in seconds (default 240)
  BASE_OVERLAY      Base CADC overlay image
  OVERLAY_DIR       Temporary worker overlay directory (default: $MYSCRATCH/cadc_upload_overlays)
  VCP_CMD           Optional explicit path to vcp
  VMV_CMD           Optional path to vmv in the same CADC environment.
                    If absent, reuse a recognized vcp container wrapper.
  CHECKSUM_PYTHON    Python 3.6+ for reading checksum manifests (default python3)
  CHECKSUM_MANIFEST_NAME  Small manifest filename in each remote run directory
  CHECKSUM_RECEIPT_DIR    Persistent pending-upload receipts under OVERLAY_DIR
  CHECKSUM_WAIT_SECONDS   Maximum response wait (default 1800 seconds)
  CHECKSUM_POLL_SECONDS   Response polling interval (default 10 seconds)
  CANFAR_API             Session API (default https://ws-uv.canfar.net/skaha/v1)
  CANFAR_CHECKSUM_IMAGE  Batch image with Bash/Python 3 (default skaha/astroml:latest)
  CADC_PYTHON_CMD        Optional container-Python wrapper; normally derived
                        automatically from the existing vcp wrapper

Default checksum mode, ON SETONIX; no SSH into CANFAR or manually running session is needed:
  ./vcp_scratch_v3tk_v768_to_canfar.sh 7000 NGC4698
All selected products are checked, not only continuum cubes. Small copies of
this script and request/response JSON files are retained in each run directory.
One-core, 1 GB CANFAR jobs compute actual byte hashes and are removed afterward.
Missing/expired certificates prompt for login; valid certificates are reused.
Recognized membership/TLS/reset failures trigger one coordinated certificate
refresh per run, even for a date-valid certificate. Login may prompt again.
Safe reads retry once; uncertain submissions/moves are not replayed.
Job/API failures stop verification and never fall back to trusting metadata.

Optional manual mode: first run this INSIDE CANFAR with the same selections:
  bash vcp_scratch_v3tk_v768_to_canfar.sh --canfar-manifest --cont-only 7000 NGC4698
Then run this uploader with --manual-checksum on Setonix. Only the small .checksum_manifest.json
returns from CANFAR; FITS files are NEVER downloaded for verification.
Matching actual hashes are skipped. Mismatches are preserved as .corrupt_*
objects using vmv before fresh upload. Missing files are uploaded directly.
After any upload, regenerate the manifest INSIDE CANFAR and rerun this uploader
to verify the replacement. Reusing the pre-upload manifest cannot trigger
another repair; persistent receipts enforce this. Keep both product trees idle.
Exit codes: 0=all verified, or all transfers completed in explicit no-checksum mode;
            1=error, 2=awaiting fresh manifest in manual-checksum mode.
Forced mode preserves .overwrite_backup_* files (not proven corrupt), including
partial destinations before retry. Existing pending receipts are refreshed, not
cleared as verified. Verify forced uploads later using default automatic mode;
manual snapshots are blocked for files marked by forced receipts.
Restore/check destinations if a replacement upload fails.
Do not remove pending receipts or reuse old manifest snapshots during repair.
USAGE
}

# Upload selected nGIST v7.6.8 products directly from Setonix scratch to CANFAR.
#
# Improvements in this version:
#   - refresh date-valid certificates once on recognized CANFAR service/auth errors
#   - explicit forced overwrite without CANFAR checksums preserves destination backups
#   - report checksum job status and image-pull failures; capture events before cleanup
#   - final summaries identify each failed run, galaxy, file and processing stage
#   - up to five galaxies run concurrently, sharing ten parallel file slots by default
#   - each parallel file uses isolated staging; writable overlays stay exclusive
#   - the known CADC container wrapper shares a read-only overlay, avoiding large copies
#   - temporary worker overlays are stored on /scratch, not /software
#   - each transfer sends one file using its own staging directory
#   - each file is retried independently with exponential backoff
#   - one failed file does not prevent the remaining files in that galaxy from running
#   - compares Setonix bytes with SHA-256 computed inside CANFAR
#   - missing continuum cubes fail the galaxy instead of silently succeeding
#   - --cont-only selects only continuum cubes for targeted recovery
#   - confirmed manifest mismatches quarantine only the exact affected object
#   - never downloads FITS products to Setonix for verification
#   - a fresh CANFAR manifest verifies replacements after upload
#   - --canfar-manifest uses embedded Python; only this script is needed
#   - missing vmv wrappers reuse the known vcp container launch configuration
#   - automatic API-launched checksum jobs need no manually started watcher
#   - valid CADC certificates are reused until a recognized remote failure
#   - no vls/vmkdir dependency; vcp handles the remote galaxy directory as before
#
# Usage:
#   ./vcp_scratch_v3tk_v768_to_canfar.sh [OPTIONS] [normal|7000] [GALID ...]
#
# Examples:
#   ./vcp_scratch_v3tk_v768_to_canfar.sh 7000 NGC4607 NGC4698
#   JOBS=1 ./vcp_scratch_v3tk_v768_to_canfar.sh 7000
#   TRANSFER_JOBS=8 ./vcp_scratch_v3tk_v768_to_canfar.sh 7000 NGC4698
#   FILE_RETRIES=8 RETRY_BASE_SLEEP=20 ./vcp_scratch_v3tk_v768_to_canfar.sh normal

# The same script runs on CANFAR and Setonix; no separate helper is needed.
checksum_manifest() {
    "${CHECKSUM_PYTHON:-python3}" - "$@" <<'CANFAR_CHECKSUM_PYTHON'
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
        filenames = getattr(args, "filenames", None) or [suffix if suffix in ("CONFIG", "LOGFILE") else galaxy + suffix for suffix in suffixes]
        allowed = {suffix if suffix in ("CONFIG", "LOGFILE") else galaxy + suffix for suffix in SUFFIXES}
        for filename in filenames:
            if filename not in allowed:
                raise ValueError("Unexpected product filename: " + filename)
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
    if getattr(args, "request_id", None):
        manifest["request_id"] = args.request_id
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


def create_request(args):
    request_id = uuid.uuid4().hex
    request = {"request_id": request_id, "run": args.run,
               "galaxy": args.galaxy, "filenames": args.filenames}
    Path(args.output).write_text(json.dumps(request) + "\n")
    print(request_id)
    return 0


def response_matches(args):
    with Path(args.manifest).open("rb") as stream:
        data = stream.read(4 * 1024 * 1024 + 1)
    if len(data) > 4 * 1024 * 1024:
        raise ValueError("Response exceeds 4 MiB")
    manifest = json.loads(data.decode("utf-8"))
    return 0 if manifest.get("request_id") == args.request_id else 1


def process_request(args):
    path = Path(args.path).resolve()
    request_id = path.name[len(".checksum_request_"):-len(".json")]
    if path.name != ".checksum_request_" + request_id + ".json" or not re.fullmatch(r"[0-9a-f]{32}", request_id):
        raise ValueError("Invalid checksum request filename")
    with path.open("rb") as stream:
        data = stream.read(65537)
    if len(data) > 65536:
        raise ValueError("Request exceeds 64 KiB")
    request = json.loads(data.decode("utf-8"))
    run = request.get("run")
    if run not in RUNS or path.parent.name != RUNS[run] or request.get("request_id") != request_id:
        raise ValueError("Request ID/run/path mismatch")
    filenames = request.get("filenames")
    if not isinstance(filenames, list) or not filenames or len(filenames) > len(SUFFIXES) or not all(isinstance(f, str) for f in filenames):
        raise ValueError("Invalid filename selection")
    output = path.parent / (".checksum_response_" + request_id + ".json")
    job = argparse.Namespace(root=str(path.parent.parent), run=run,
        galaxies=[request.get("galaxy")], cont_only=False, filenames=filenames,
        output=str(output), request_id=request_id)
    return generate(job)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    # Python 3.6 has no required= keyword for add_subparsers.
    sub = parser.add_subparsers(dest="command")
    sub.required = True
    generator = sub.add_parser("generate", help="Run inside CANFAR; reads file bytes, not VOS checksum metadata")
    generator.add_argument("--root", default="/arc/projects/mauve/products")
    generator.add_argument("--cont-only", action="store_true")
    generator.add_argument("--output")
    generator.add_argument("--filename", dest="filenames", action="append")
    generator.add_argument("run", choices=RUNS)
    generator.add_argument("galaxies", nargs="*")
    reader = sub.add_parser("lookup", help="Read a manifest entry on Setonix")
    reader.add_argument("--manifest", required=True)
    reader.add_argument("--run", choices=RUNS, required=True)
    reader.add_argument("--galaxy", required=True)
    reader.add_argument("--filename", required=True)
    processor = sub.add_parser("process-request")
    processor.add_argument("path")
    request = sub.add_parser("request")
    request.add_argument("--output", required=True)
    request.add_argument("--run", choices=RUNS, required=True)
    request.add_argument("--galaxy", required=True)
    request.add_argument("--filename", dest="filenames", action="append", required=True)
    response = sub.add_parser("response-matches")
    response.add_argument("--manifest", required=True)
    response.add_argument("--request-id", required=True)
    args = parser.parse_args()
    try:
        actions = {"generate": generate, "lookup": lookup,
                   "request": create_request, "response-matches": response_matches,
                   "process-request": process_request}
        return actions[args.command](args)
    except KeyboardInterrupt:
        return 130
    except Exception as error:
        print("ERROR:", error, file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
CANFAR_CHECKSUM_PYTHON
}

# Dispatch before Setonix paths, overlays or authentication are accessed.
if [ "${1:-}" = --canfar-request ]; then
    shift
    checksum_manifest process-request "$@"
    exit $?
fi
if [ "${1:-}" = --canfar-manifest ]; then
    shift
    checksum_manifest generate "$@"
    exit $?
fi

record_issue() {
    local kind=$1 run=$2 galaxy=$3 filename=$4 stage=$5 detail=$6 entry
    [ -n "${ISSUE_DIR:-}" ] || return 0
    detail=${detail//$'\n'/ }
    detail=${detail//$'\t'/ }
    mkdir -p "$ISSUE_DIR" || return 0
    entry=$(mktemp "$ISSUE_DIR/issue.XXXXXX") || return 0
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$kind" "$run" "$galaxy" "$filename" "$stage" "$detail" > "$entry"
}

print_error_summary() {
    [ -n "${ISSUE_DIR:-}" ] && [ -d "$ISSUE_DIR" ] || return 0
    "${CHECKSUM_PYTHON:-python3}" - "$ISSUE_DIR" >&2 <<'UPLOADER_ERROR_SUMMARY'
from pathlib import Path
import sys
entries = set()
for path in Path(sys.argv[1]).glob("issue.*"):
    fields = path.read_text().rstrip("\n").split("\t", 5)
    if len(fields) == 6:
        entries.add(tuple(fields))
for kind, heading in [("ERROR", "FINAL ERROR SUMMARY"), ("PENDING", "AWAITING VERIFICATION SUMMARY")]:
    rows = sorted(row for row in entries if row[0] == kind)
    if not rows:
        continue
    print("\n" + heading + " (" + str(len(rows)) + " issue(s))")
    for _, run, galaxy, filename, stage, detail in rows:
        print("  run=" + run + " | galaxy=" + galaxy + " | file=" + filename + " | stage=" + stage)
        print("    " + detail)
UPLOADER_ERROR_SUMMARY
}

finish_reporting() {
    local result=$? summary=""
    trap - EXIT
    if [ -n "${ISSUE_DIR:-}" ] && compgen -G "$ISSUE_DIR/issue.*" >/dev/null; then
        summary=$(print_error_summary 2>&1)
    elif [ "$result" -ne 0 ]; then
        summary=$(printf '\nFINAL ERROR SUMMARY\n  run=%s | galaxy=%s | file=(not identified) | stage=%s\n    Run exited with status %s; see the diagnostic above.\n' \
            "${RUNS[*]:-${REQUESTED_RUNS[*]:-not resolved}}" "${REQUESTED_GALAXIES[*]:-all selected}" "${CURRENT_STAGE:-startup}" "$result")
    fi
    if [ -n "${worker_dir:-}" ] && declare -F cleanup >/dev/null; then cleanup; fi
    [ -z "$summary" ] || printf '%s\n' "$summary" >&2
    exit "$result"
}

trap finish_reporting EXIT
CURRENT_STAGE=argument-validation
start_secs=$SECONDS

SOURCE_NORMAL=${SOURCE_NORMAL:-/scratch/pawsey1308/mauve/products/v3tk_v7.6.8}
SOURCE_7000=${SOURCE_7000:-/scratch/pawsey1308/mauve/products/v3tk_v7.6.8_7000}
DEST_NORMAL=${DEST_NORMAL:-arc:projects/mauve/products/v3tk_v7.6.8}
DEST_7000=${DEST_7000:-arc:projects/mauve/products/v3tk_v7.6.8_7000}
CADC_USER=${CADC_USER:-RongjunHuang}

BASE_OVERLAY=${BASE_OVERLAY:-/software/projects/pawsey1308/containers/cadc_overlay.img}
DEFAULT_OVERLAY_DIR="${MYSCRATCH:-/scratch/pawsey1308/$USER}/cadc_upload_overlays"
OVERLAY_DIR=${OVERLAY_DIR:-$DEFAULT_OVERLAY_DIR}
CHECKSUM_PYTHON=${CHECKSUM_PYTHON:-python3}
CHECKSUM_MANIFEST_NAME=${CHECKSUM_MANIFEST_NAME:-.checksum_manifest.json}
CHECKSUM_RECEIPT_DIR=${CHECKSUM_RECEIPT_DIR:-$OVERLAY_DIR/checksum_receipts}
AUTO_CHECKSUM=${AUTO_CHECKSUM:-1}
NO_CHECKSUM=0
FORCE_OVERWRITE=0
CHECKSUM_MODE_OPTION=0
CHECKSUM_WAIT_SECONDS=${CHECKSUM_WAIT_SECONDS:-1800}
CHECKSUM_POLL_SECONDS=${CHECKSUM_POLL_SECONDS:-10}
CANFAR_API=${CANFAR_API:-https://ws-uv.canfar.net/skaha/v1}
CANFAR_CHECKSUM_IMAGE=${CANFAR_CHECKSUM_IMAGE:-images.canfar.net/skaha/astroml:latest}
UPLOADER_PATH=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")

# Process up to five galaxies concurrently.
# Set JOBS explicitly to use fewer workers.
MAX_JOBS=5
JOBS=${JOBS:-5}
# Shared across all galaxies, rather than multiplied by JOBS.
TRANSFER_JOBS=${TRANSFER_JOBS:-10}
MAX_TRANSFER_JOBS=10
CADC_READ_ONLY=${CADC_READ_ONLY:-auto}

# Retry policy for each individual file.
FILE_RETRIES=${FILE_RETRIES:-5}
RETRY_BASE_SLEEP=${RETRY_BASE_SLEEP:-30}
RETRY_MAX_SLEEP=${RETRY_MAX_SLEEP:-240}


RUN_ID=$$

find_command() {
    local command_name=$1
    local override_var
    override_var=$(printf '%s_CMD' "$command_name" | tr '[:lower:]' '[:upper:]')

    if [ -n "${!override_var:-}" ]; then
        if [ -x "${!override_var}" ]; then
            printf '%s\n' "${!override_var}"
            return 0
        fi
        echo "ERROR: $override_var is set but not executable: ${!override_var}" >&2
        return 1
    fi

    local candidate
    if command -v "$command_name" >/dev/null 2>&1; then
        command -v "$command_name"
        return 0
    fi

    for candidate in \
        "/home/$USER/bin/$command_name" \
        "/work/$USER/bin/$command_name" \
        "/cadcenv/bin/$command_name"; do
        if [ -x "$candidate" ]; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done

    return 1
}

make_cadc_python_wrapper() {
    if [ -n "${CADC_PYTHON_CMD:-}" ]; then
        [ -x "$CADC_PYTHON_CMD" ] || return 1
        return 0
    fi
    CADC_PYTHON_CMD="$worker_dir/cadc_python_wrapper"
    "$CHECKSUM_PYTHON" - "$VCP_CMD" "$CADC_PYTHON_CMD" <<'CADC_PYTHON_WRAPPER'
from pathlib import Path
import sys
source = Path(sys.argv[1]).read_text()
command = '/cadcenv/bin/vcp "$@"'
if not source.startswith("#!/") or source.count(command) != 1:
    sys.exit("ERROR: unsupported vcp wrapper; set CADC_PYTHON_CMD to a CADC-container Python wrapper")
target = Path(sys.argv[2])
target.write_text(source.replace(command, '/cadcenv/bin/python - "$@"'))
target.chmod(0o700)
CADC_PYTHON_WRAPPER
}

ensure_cadc_certificate() {
    local result
    if CADC_OVERLAY="$BASE_OVERLAY:ro" "$CADC_PYTHON_CMD" <<'CADC_CERT_CHECK'
# CADC_CERT_CHECK
from pathlib import Path
import ssl
import sys
import time
path = Path.home() / ".ssl/cadcproxy.pem"
if not path.exists():
    sys.exit(10)
try:
    cert = ssl._ssl._test_decode_cert(str(path))
    starts = ssl.cert_time_to_seconds(cert["notBefore"])
    expires = ssl.cert_time_to_seconds(cert["notAfter"])
except Exception:
    sys.exit(12)
if expires <= time.time():
    sys.exit(11)
if starts > time.time():
    sys.exit(12)
print("Reusing valid CADC certificate (expires " + cert["notAfter"] + ").")
CADC_CERT_CHECK
    then
        return 0
    else
        result=$?
    fi
    if [ "$result" -ne 10 ] && [ "$result" -ne 11 ]; then
        echo "ERROR: could not validate CADC certificate; refusing an unexpected login prompt." >&2
        return 1
    fi
    echo "CADC certificate missing or expired; renewing now."
    cadc-get-cert -u "$CADC_USER" || return 1
    # Verify renewal without prompting a second time.
    validate_cadc_certificate
}

validate_cadc_certificate() {
    CADC_OVERLAY="$BASE_OVERLAY:ro" "$CADC_PYTHON_CMD" <<'CADC_CERT_RECHECK'
from pathlib import Path
import ssl
import sys
import time
try:
    cert = ssl._ssl._test_decode_cert(str(Path.home()/".ssl/cadcproxy.pem"))
    now = time.time()
    sys.exit(0 if ssl.cert_time_to_seconds(cert["notBefore"]) <= now < ssl.cert_time_to_seconds(cert["notAfter"]) else 1)
except Exception:
    sys.exit(1)
CADC_CERT_RECHECK
}

refresh_cadc_certificate() {
    # All galaxy/file subprocesses share this invocation's lock and markers.
    # Keep prompts off stdout, which may contain a session ID/JSON result.
    (
        flock 9 || return 1
        if [ -f "$worker_dir/cert_refresh.attempted" ]; then
            [ -f "$worker_dir/cert_refresh.ok" ]
            return $?
        fi
        : > "$worker_dir/cert_refresh.attempted" || return 1
        echo "WARNING: CANFAR rejected the current certificate or reset authentication; refreshing CADC certificate once." >&2
        # Background workers may inherit /dev/null as stdin. Use the terminal
        # for the normal password prompt when available; never store a password.
        if ( : </dev/tty ) 2>/dev/null; then
            cadc-get-cert -u "$CADC_USER" </dev/tty >&2 || return 1
        else
            cadc-get-cert -u "$CADC_USER" >&2 || return 1
        fi
        validate_cadc_certificate || return 1
        : > "$worker_dir/cert_refresh.ok" || return 1
        echo "CADC certificate refreshed and validity checked." >&2
    ) 9>"$worker_dir/cert_refresh.lock"
}

cadc_with_recovery() {
    local mode=$1 result output error already_refreshed=0
    shift
    [ ! -f "$worker_dir/cert_refresh.ok" ] || already_refreshed=1
    output=$(mktemp "$worker_dir/cadc_stdout.XXXXXX") || return 1
    error=$(mktemp "$worker_dir/cadc_stderr.XXXXXX") || { rm -f "$output"; return 1; }
    if "$@" >"$output" 2>"$error"; then
        cat "$output"
        cat "$error" >&2
        rm -f "$output" "$error"
        return 0
    else
        result=$?
    fi
    cat "$output" "$error" >&2
    # Match the supplied outage signatures, not every HTTP 500 or permission
    # failure. A service outage can persist even after certificate renewal.
    if [ "$already_refreshed" -eq 0 ] && LC_ALL=C grep -Eiq \
        'failed to check membership with group service|SSLHandshakeException|decrypt_error|Connection reset by peer|HTTP[ /]+401([^0-9]|$)|certificate (has )?expired' "$output" "$error"; then
        if refresh_cadc_certificate; then
            rm -f "$output" "$error"
            if [ "$mode" = retry ]; then
                echo "Retrying CADC operation once with the refreshed certificate." >&2
                "$@"
                return $?
            fi
            echo "Certificate refreshed; recovery did not replay the failed transfer/mutation. Existing file retry rules still apply." >&2
            return "$result"
        fi
        echo "ERROR: CADC certificate refresh failed; no further automatic refresh in this run." >&2
    fi
    rm -f "$output" "$error"
    return "$result"
}

checksum_job_api() {
    local overlay=$1 action=$2
    if [ "$action" = submit ]; then
        # Recover authentication with an idempotent GET before the POST.
        cadc_with_recovery retry checksum_job_api_raw "$overlay" probe || return 1
        cadc_with_recovery no-retry checksum_job_api_raw "$@"
    else
        cadc_with_recovery retry checksum_job_api_raw "$@"
    fi
}

checksum_job_api_raw() {
    local overlay=$1
    shift
    CADC_OVERLAY="${overlay%:ro}:ro" "$CADC_PYTHON_CMD" "$CANFAR_API" "$@" <<'CADC_JOB_API'
# CADC_JOB_API
from pathlib import Path
import re
import requests
import shlex
import sys
base, action = sys.argv[1:3]
if not base.startswith("https://"):
    sys.exit("ERROR: CANFAR_API must use HTTPS")
session = requests.Session()
cert = str(Path.home()/".ssl/cadcproxy.pem")
session.cert = (cert, cert)
try:
    if action == "probe":
        response = session.get(base.rstrip('/')+'/session', timeout=(10,30))
        response.raise_for_status()
    elif action == "submit":
        image, worker, request, name = sys.argv[3:7]
        parameters = {"name":name, "image":image, "type":"headless",
            "cores":1, "ram":1, "cmd":"/bin/bash",
            "args":shlex.quote(worker) + " --canfar-request " + shlex.quote(request)}
        # POST is deliberately not retried: a timed-out submission can have
        # created a job, so repeating it could launch duplicate jobs.
        response = session.post(base.rstrip('/')+'/session', data=parameters, timeout=(15,60))
        if response.status_code >= 400:
            detail = re.sub(r'https://[^\s]+', '[URL redacted]', response.text[:2000])
            raise ValueError("Session creation HTTP " + str(response.status_code) + ": " + detail)
        response.raise_for_status()
        job_id = response.text.strip()
        if not re.fullmatch(r"[A-Za-z0-9_-]+", job_id):
            raise ValueError("Unexpected session ID response")
        print(job_id)
    else:
        job_id = sys.argv[3]
        if not re.fullmatch(r"[A-Za-z0-9_-]+", job_id):
            raise ValueError("Invalid session ID")
        url = base.rstrip('/')+'/session/'+job_id
        if action == "delete":
            response = session.delete(url, timeout=(10,30))
            if response.status_code != 404:
                response.raise_for_status()
        elif action == "status":
            response = session.get(url, timeout=(10,30))
            if response.status_code >= 400:
                raise ValueError("Session status HTTP " + str(response.status_code) + ": " + response.text[:2000])
            response.raise_for_status()
            print(response.json().get("status", "Unknown"))
        elif action in ("logs", "events"):
            response = session.get(url, params={"view":action}, timeout=(10,30))
            response.raise_for_status()
            print(re.sub(r'https://[^\s]+', '[URL redacted]', response.text[-4000:]))
        else:
            raise ValueError("Unknown API operation")
except Exception as error:
    # requests.HTTPError normally omits the response body containing the
    # server-side group-service/TLS diagnostics needed for recovery.
    response = getattr(error, "response", None)
    detail = "" if response is None else ": " + response.text[:2000]
    message = re.sub(r'https://[^\s]+', '[URL redacted]', str(error) + detail)
    print("ERROR: CANFAR checksum job API: " + message, file=sys.stderr)
    sys.exit(1)
CADC_JOB_API
}

resolve_vmv() {
    local wrapper=$1 command_path
    if command_path=$(find_command vmv); then
        printf '%s\n' "$command_path"
        return 0
    fi
    # An explicit override must not be silently ignored.
    if [ -n "${VMV_CMD:-}" ]; then return 1; fi
    # Clone only the known CADC shell-wrapper command, preserving all binds,
    # HOME, container paths and CADC_OVERLAY handling. The clone is temporary.
    "$CHECKSUM_PYTHON" - "$VCP_CMD" "$wrapper" <<'CADC_VMV_WRAPPER'
from pathlib import Path
import sys
try:
    source = Path(sys.argv[1]).read_text()
    command = '/cadcenv/bin/vcp "$@"'
    if not source.startswith("#!/") or source.count(command) != 1:
        raise ValueError("vcp wrapper does not contain the expected CADC command")
    target = Path(sys.argv[2])
    target.write_text(source.replace(command, '/cadcenv/bin/vmv "$@"'))
    target.chmod(0o700)
except Exception as error:
    print("ERROR: cannot derive vmv container wrapper: " + str(error), file=sys.stderr)
    sys.exit(1)
CADC_VMV_WRAPPER
    if [ "$?" -ne 0 ]; then return 1; fi
    printf '%s\n' "$wrapper"
}

ALL_GALAXIES=(
    IC3392
    NGC4064
    NGC4189
    NGC4192
    NGC4216
    NGC4222
    NGC4254
    NGC4293
    NGC4294
    NGC4298
    NGC4302
    NGC4321
    NGC4330
    NGC4351
    NGC4380
    NGC4383
    NGC4388
    NGC4394
    NGC4396
    NGC4405
    NGC4402
    NGC4419
    NGC4424
    NGC4450
    NGC4457
    NGC4501
    NGC4522
    NGC4535
    NGC4548
    NGC4567_8
    NGC4569
    NGC4579
    NGC4580
    NGC4606
    NGC4607
    NGC4654
    NGC4689
    NGC4694
    NGC4698
)

PRODUCT_SUFFIXES=(
    CONFIG
    _sfh_maps.fits
    _gas_bin_maps.fits
    _sfh_weights.fits
    _gas_spaxel_maps.fits
    _spatial_binning_maps.fits
    _kin_maps.fits
    _cont_cube.fits
    LOGFILE
    _mask.fits
)


format_runtime() {
    local total=$1
    printf '%02d:%02d:%02d' \
        $((total / 3600)) \
        $(((total % 3600) / 60)) \
        $((total % 60))
}

retry_sleep_seconds() {
    local failed_attempt=$1
    local delay=$RETRY_BASE_SLEEP
    local i

    # attempt 1 -> base, attempt 2 -> 2*base, etc.
    for ((i = 1; i < failed_attempt; i++)); do
        delay=$((delay * 2))
        if [ "$delay" -ge "$RETRY_MAX_SLEEP" ]; then
            delay=$RETRY_MAX_SLEEP
            break
        fi
    done

    printf '%s\n' "$delay"
}

is_positive_integer() {
    [[ "$1" =~ ^[0-9]+$ ]] && [ "$1" -ge 1 ]
}

is_known_galaxy() {
    local candidate=$1
    local galaxy
    for galaxy in "${ALL_GALAXIES[@]}"; do
        if [ "$candidate" = "$galaxy" ]; then
            return 0
        fi
    done
    return 1
}

source_path() {
    local source_base=$1
    local galaxy=$2
    local suffix=$3
    case "$suffix" in
        CONFIG|LOGFILE)
            printf '%s/%s/%s\n' "$source_base" "$galaxy" "$suffix"
            ;;
        *)
            printf '%s/%s/%s%s\n' "$source_base" "$galaxy" "$galaxy" "$suffix"
            ;;
    esac
}

run_source_base() {
    case "$1" in
        normal) printf '%s\n' "$SOURCE_NORMAL" ;;
        7000) printf '%s\n' "$SOURCE_7000" ;;
    esac
}

run_dest_base() {
    case "$1" in
        normal) printf '%s\n' "$DEST_NORMAL" ;;
        7000) printf '%s\n' "$DEST_7000" ;;
    esac
}

available_galaxies_for_run() {
    local source_base=$1
    local galaxy
    for galaxy in "${ALL_GALAXIES[@]}"; do
        if [ -d "$source_base/$galaxy" ]; then
            printf '%s\n' "$galaxy"
        fi
    done
}

wait_for_overlay() {
    local overlay=$1
    local tries=0
    while ! flock -n "$overlay" -c true 2>/dev/null; do
        tries=$((tries + 1))
        if [ "$tries" -ge 30 ]; then
            echo "ERROR: overlay still busy after 30 seconds: $overlay" >&2
            return 1
        fi
        sleep 1
    done
}

prepare_worker_overlays() {
    local use_read_only=0 i worker_overlay
    if [ "$CADC_READ_ONLY" = 1 ]; then
        use_read_only=1
    elif [ "$CADC_READ_ONLY" = auto ]; then
        if "$CHECKSUM_PYTHON" - "$VCP_CMD" <<'CADC_OVERLAY_MODE'
from pathlib import Path
import sys
try:
    source = Path(sys.argv[1]).read_text()
except (OSError, UnicodeError):
    sys.exit(1)
sys.exit(0 if source.startswith("#!/") and source.count('/cadcenv/bin/vcp "$@"') == 1 else 1)
CADC_OVERLAY_MODE
        then use_read_only=1; fi
    fi
    if [ "$use_read_only" -eq 1 ]; then
        echo "CADC overlay mode: shared read-only (no worker image copies)."
        for ((i = 0; i < TRANSFER_JOBS; i++)); do
            worker_overlays+=("$BASE_OVERLAY:ro")
        done
    else
        echo "CADC overlay mode: separate writable worker copies."
        for ((i = 0; i < TRANSFER_JOBS; i++)); do
            worker_overlay="$OVERLAY_DIR/cadc_overlay_v3tk_upload_${RUN_ID}_${i}.img"
            rm -f "$worker_overlay"
            if ! cp --reflink=auto "$BASE_OVERLAY" "$worker_overlay"; then
                rm -f "$worker_overlay"
                echo "ERROR: failed to copy CADC overlay for slot $i: $worker_overlay" >&2
                return 1
            fi
            wait_for_overlay "$worker_overlay" || return 1
            worker_overlays+=("$worker_overlay")
        done
    fi
}

DRY_RUN=0
CONT_ONLY=0
REQUESTED_RUNS=()
REQUESTED_GALAXIES=()

while [ "$#" -gt 0 ]; do
    case "$1" in
        -n|--dry-run)
            DRY_RUN=1
            ;;
        --cont-only)
            CONT_ONLY=1
            ;;
        --auto-checksum)
            AUTO_CHECKSUM=1
            CHECKSUM_MODE_OPTION=1
            ;;
        --manual-checksum)
            AUTO_CHECKSUM=0
            CHECKSUM_MODE_OPTION=1
            ;;
        --force-overwrite)
            FORCE_OVERWRITE=1
            ;;
        --no-checksum)
            NO_CHECKSUM=1
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --)
            shift
            while [ "$#" -gt 0 ]; do
                REQUESTED_GALAXIES+=("$1")
                shift
            done
            break
            ;;
        -*)
            echo "ERROR: unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
        normal|7000)
            REQUESTED_RUNS+=("$1")
            ;;
        *)
            REQUESTED_GALAXIES+=("$1")
            ;;
    esac
    shift
done

if [ "$FORCE_OVERWRITE" -ne "$NO_CHECKSUM" ]; then
    echo "ERROR: --force-overwrite requires --no-checksum, and --no-checksum requires --force-overwrite." >&2
    exit 2
fi
if [ "$NO_CHECKSUM" -eq 1 ]; then
    if [ "$CHECKSUM_MODE_OPTION" -eq 1 ]; then
        echo "ERROR: cannot combine --no-checksum with --auto-checksum or --manual-checksum." >&2
        exit 2
    fi
    AUTO_CHECKSUM=0
fi

if [ "$CONT_ONLY" -eq 1 ]; then
    PRODUCT_SUFFIXES=(_cont_cube.fits)
fi

if ! is_positive_integer "$CHECKSUM_WAIT_SECONDS" || ! is_positive_integer "$CHECKSUM_POLL_SECONDS"; then
    echo "ERROR: checksum wait/poll intervals must be positive integers." >&2
    exit 2
fi
if [[ "$AUTO_CHECKSUM" != 0 && "$AUTO_CHECKSUM" != 1 ]]; then
    echo "ERROR: AUTO_CHECKSUM must be 0 or 1." >&2
    exit 2
fi

if ! is_positive_integer "$JOBS"; then
    echo "ERROR: JOBS must be a positive integer, got: $JOBS" >&2
    exit 2
fi
if [ "$JOBS" -gt "$MAX_JOBS" ]; then
    echo "Capping JOBS from $JOBS to $MAX_JOBS for transfer reliability."
    JOBS=$MAX_JOBS
fi
if ! is_positive_integer "$TRANSFER_JOBS"; then
    echo "ERROR: TRANSFER_JOBS must be a positive integer, got: $TRANSFER_JOBS" >&2
    exit 2
fi
if [ "$TRANSFER_JOBS" -gt "$MAX_TRANSFER_JOBS" ]; then
    echo "Capping TRANSFER_JOBS from $TRANSFER_JOBS to $MAX_TRANSFER_JOBS."
    TRANSFER_JOBS=$MAX_TRANSFER_JOBS
fi
case "$CADC_READ_ONLY" in
    auto|0|1) ;;
    *) echo "ERROR: CADC_READ_ONLY must be auto, 0, or 1." >&2; exit 2 ;;
esac

if ! is_positive_integer "$FILE_RETRIES"; then
    echo "ERROR: FILE_RETRIES must be a positive integer, got: $FILE_RETRIES" >&2
    exit 2
fi
if ! is_positive_integer "$RETRY_BASE_SLEEP"; then
    echo "ERROR: RETRY_BASE_SLEEP must be a positive integer, got: $RETRY_BASE_SLEEP" >&2
    exit 2
fi
if ! is_positive_integer "$RETRY_MAX_SLEEP"; then
    echo "ERROR: RETRY_MAX_SLEEP must be a positive integer, got: $RETRY_MAX_SLEEP" >&2
    exit 2
fi

if [ "${#REQUESTED_RUNS[@]}" -gt 0 ]; then
    RUNS=("${REQUESTED_RUNS[@]}")
else
    RUNS=(normal 7000)
fi

WORK_ITEMS=()
for run in "${RUNS[@]}"; do
    source_base=$(run_source_base "$run")
    if [ "${#REQUESTED_GALAXIES[@]}" -gt 0 ]; then
        GALAXIES=("${REQUESTED_GALAXIES[@]}")
    else
        mapfile -t GALAXIES < <(available_galaxies_for_run "$source_base")
    fi

    for galaxy in "${GALAXIES[@]}"; do
        if ! is_known_galaxy "$galaxy"; then
            echo "ERROR: unknown galaxy ID: $galaxy" >&2
            echo "Known galaxy IDs: ${ALL_GALAXIES[*]}" >&2
            exit 2
        fi
        WORK_ITEMS+=("${run}|${galaxy}")
    done
done

if [ "${#WORK_ITEMS[@]}" -eq 0 ]; then
    echo "No available galaxies selected."
    exit 0
fi

EFFECTIVE_JOBS=$JOBS
if [ "${#WORK_ITEMS[@]}" -lt "$EFFECTIVE_JOBS" ]; then
    EFFECTIVE_JOBS=${#WORK_ITEMS[@]}
fi
# Avoid unused overlays for small selections such as one continuum cube.
selected_file_limit=$((${#WORK_ITEMS[@]} * ${#PRODUCT_SUFFIXES[@]}))
if [ "$TRANSFER_JOBS" -gt "$selected_file_limit" ]; then
    TRANSFER_JOBS=$selected_file_limit
fi

if [ "$DRY_RUN" -eq 1 ]; then
    if [ "$NO_CHECKSUM" -eq 1 ]; then echo "Mode: FORCE OVERWRITE; destination bytes UNVERIFIED (CANFAR checksums disabled)."; fi
    echo "Requested runs: ${RUNS[*]}"
    if [ "${#REQUESTED_GALAXIES[@]}" -gt 0 ]; then
        echo "Requested galaxies: ${REQUESTED_GALAXIES[*]}"
    else
        echo "Requested galaxies: all available"
    fi
    echo "Work item count: ${#WORK_ITEMS[@]}"
    echo "Effective workers: $EFFECTIVE_JOBS"
    echo "Shared parallel file slots: $TRANSFER_JOBS"
    echo "File retries: $FILE_RETRIES"

    for work_item in "${WORK_ITEMS[@]}"; do
        run=${work_item%%|*}
        galaxy=${work_item#*|}
        source_base=$(run_source_base "$run")
        dest_base=$(run_dest_base "$run")
        echo "$run $galaxy:"
        echo "  Destination: ${dest_base}/${galaxy}/"
        for suffix in "${PRODUCT_SUFFIXES[@]}"; do
            path=$(source_path "$source_base" "$galaxy" "$suffix")
            if [ -f "$path" ]; then
                echo "  [present] $path"
            else
                echo "  [missing] $path"
            fi
        done
    done
    exit 0
fi

CURRENT_STAGE=preflight
command -v cadc-get-cert >/dev/null 2>&1 || {
    echo "ERROR: cadc-get-cert is not available in PATH." >&2
    exit 1
}
VCP_CMD=$(find_command vcp) || {
    echo "ERROR: vcp is not available in PATH." >&2
    exit 1
}
command -v flock >/dev/null 2>&1 || {
    echo "ERROR: flock is not available in PATH." >&2
    exit 1
}
for required_command in sha256sum stat head pgrep "$CHECKSUM_PYTHON"; do
    command -v "$required_command" >/dev/null 2>&1 || {
        echo "ERROR: $required_command is not available in PATH." >&2
        exit 1
    }
done

if [ ! -f "$BASE_OVERLAY" ]; then
    echo "ERROR: CADC overlay not found: $BASE_OVERLAY" >&2
    exit 1
fi

if ! mkdir -p "$OVERLAY_DIR"; then
    echo "ERROR: could not create worker overlay directory: $OVERLAY_DIR" >&2
    exit 1
fi
if [ ! -w "$OVERLAY_DIR" ]; then
    echo "ERROR: worker overlay directory is not writable: $OVERLAY_DIR" >&2
    exit 1
fi

echo "Base CADC overlay: $BASE_OVERLAY"
echo "Worker overlay directory: $OVERLAY_DIR"
echo "Worker count: $EFFECTIVE_JOBS"
echo "Shared parallel file slots: $TRANSFER_JOBS"
echo "Per-file retries: $FILE_RETRIES"
echo "Retry backoff: ${RETRY_BASE_SLEEP}s -> max ${RETRY_MAX_SLEEP}s"
if [ "$NO_CHECKSUM" -eq 1 ]; then
    echo "Mode: FORCE OVERWRITE; destination bytes UNVERIFIED (CANFAR checksums disabled)."
fi

worker_dir=$(mktemp -d)
ISSUE_DIR="$worker_dir/issues"
worker_overlays=()

cleanup() {
    local overlay
    local job_file job_id
    for job_file in "$worker_dir"/checksum_job_*.id; do
        [ -f "$job_file" ] || continue
        read -r job_id < "$job_file" || continue
        if [ -n "${CADC_PYTHON_CMD:-}" ] && [ -x "$CADC_PYTHON_CMD" ]; then
            checksum_job_api "$BASE_OVERLAY" delete "$job_id" || echo "WARNING: could not remove checksum job $job_id." >&2
        fi
    done
    for overlay in "${worker_overlays[@]:-}"; do
        [[ "$overlay" == *:ro ]] && continue
        rm -f "$overlay"
    done
    rm -rf "$worker_dir"
}

stop_workers() {
    echo "Interrupted. Stopping active transfers..." >&2
    stop_process_tree $$
    wait 2>/dev/null || true
    exit 130
}
stop_process_tree() {
    local parent=$1 child
    while read -r child; do
        [ -n "$child" ] || continue
        # Freeze the parent before terminating descendants: otherwise a killed
        # upload command can return and let its shell start the next operation.
        kill -STOP "$child" 2>/dev/null || true
        stop_process_tree "$child"
        kill -TERM "$child" 2>/dev/null || true
        kill -CONT "$child" 2>/dev/null || true
    done < <(pgrep -P "$parent" 2>/dev/null || true)
}
trap stop_workers INT TERM

CURRENT_STAGE=authentication
wait_for_overlay "$BASE_OVERLAY"
make_cadc_python_wrapper || exit 1
ensure_cadc_certificate || exit 1

CURRENT_STAGE=overlay-preparation
prepare_worker_overlays || exit 1
CURRENT_STAGE=galaxy-processing

for ((i = 0; i < EFFECTIVE_JOBS; i++)); do
    : > "$worker_dir/part_${i}.txt"
done

for ((i = 0; i < ${#WORK_ITEMS[@]}; i++)); do
    printf '%s\n' "${WORK_ITEMS[$i]}" >> "$worker_dir/part_$((i % EFFECTIVE_JOBS)).txt"
done

remote_file_exists() {
    cadc_with_recovery retry remote_file_exists_raw "$@"
}

remote_file_exists_raw() {
    local overlay=$1 remote_path=$2
    CADC_OVERLAY="${overlay%:ro}:ro" "$CADC_PYTHON_CMD" "$remote_path" <<'CADC_NODE_EXISTS'
import re
import sys
import vos
from cadcutils.exceptions import NotFoundException
try:
    # Read only the node definition; never download product bytes or use
    # cached checksums as evidence of integrity.
    node = vos.Client().get_node(sys.argv[1], limit=0, force=True)
    if node.isdir():
        raise ValueError("Destination is a directory, not a product file")
except NotFoundException:
    sys.exit(3)
except Exception as error:
    message = re.sub(r'https://[^\s]+', '[URL redacted]', str(error))
    print("ERROR: destination metadata lookup failed: " + message, file=sys.stderr)
    sys.exit(1)
CADC_NODE_EXISTS
}

backup_forced_destination() {
    local run=$1 galaxy=$2 filename=$3 overlay=$4 remote_path=$5 wrapper=$6
    local exists_status vmv_cmd backup_path
    if remote_file_exists "$overlay" "$remote_path"; then
        exists_status=0
    else
        exists_status=$?
    fi
    if [ "$exists_status" -eq 3 ]; then return 0; fi
    if [ "$exists_status" -ne 0 ]; then
        record_issue ERROR "$run" "$galaxy" "$filename" destination-lookup "Cannot establish whether destination exists; forced replacement refused."
        return 1
    fi
    vmv_cmd=$(resolve_vmv "$wrapper") || {
        record_issue ERROR "$run" "$galaxy" "$filename" backup-wrapper "vmv unavailable or wrapper unsupported; forced replacement refused."
        return 1
    }
    backup_path="${remote_path}.overwrite_backup_$(date -u +%Y%m%dT%H%M%SZ)_${RUN_ID}_${RANDOM}_${RANDOM}"
    echo "[$run $galaxy] OVERWRITE BACKUP: $remote_path -> $backup_path"
    if ! CADC_OVERLAY="$overlay" cadc_with_recovery no-retry "$vmv_cmd" -v "$remote_path" "$backup_path"; then
        record_issue ERROR "$run" "$galaxy" "$filename" backup "Remote backup move failed; replacement upload refused."
        return 1
    fi
}

upload_one_file() {
    local run=$1 galaxy=$2 path=$3 dest_base=$4 overlay=$5 stage_galaxy_dir=$6
    overlay=${CADC_SLOT_OVERLAY:-$overlay}
    local filename staged_path remote_path source_bytes source_hash current_hash entry
    local remote_bytes remote_hash manifest_id receipt_key receipt pending_manifest pending_hash
    local vmv_cmd quarantine_path attempt delay
    local manifest_path="${stage_galaxy_dir%/*}/checksum_manifest.json"
    filename=$(basename "$path")
    staged_path="${stage_galaxy_dir}/${filename}"
    remote_path="${dest_base}/${galaxy}/${filename}"

    if [[ "$filename" == *.fits ]] && [ "$(head -c 9 "$path")" != 'SIMPLE  =' ]; then
        echo "ERROR [$run $galaxy]: $filename has no primary FITS SIMPLE card; refusing upload." >&2
        record_issue ERROR "$run" "$galaxy" "$filename" source-header "Missing primary FITS SIMPLE card; upload refused."
        return 1
    fi
    source_bytes=$(stat -c %s "$path") || { record_issue ERROR "$run" "$galaxy" "$filename" source-stat "Could not read source byte count."; return 1; }
    source_hash=$(sha256sum "$path") || { record_issue ERROR "$run" "$galaxy" "$filename" source-hash "Could not hash source bytes."; return 1; }
    source_hash=${source_hash%% *}
    echo "[$run $galaxy] Source: $path ($source_bytes bytes; SHA256 $source_hash)"
    echo "[$run $galaxy] Destination: $remote_path"
    if [ "${NO_CHECKSUM:-0}" -eq 1 ]; then
        remote_bytes=MISSING
        remote_hash=UNVERIFIED
        manifest_id="forced_${RUN_ID}_${RANDOM}_${RANDOM}"
    else
        entry=$(checksum_manifest lookup --manifest "$manifest_path" --run "$run" --galaxy "$galaxy" --filename "$filename") || { record_issue ERROR "$run" "$galaxy" "$filename" initial-manifest-lookup "Manifest entry missing, invalid or unreadable; see lookup diagnostic."; return 1; }
        read -r remote_bytes remote_hash manifest_id <<<"$entry"
    fi
    receipt_key=$(printf '%s' "$remote_path" | sha256sum) || { record_issue ERROR "$run" "$galaxy" "$filename" receipt-key "Could not compute receipt key."; return 1; }
    receipt="${CHECKSUM_RECEIPT_DIR}/${receipt_key%% *}.receipt"
    if [ "${NO_CHECKSUM:-0}" -eq 0 ] && [ -f "$receipt" ]; then
        read -r pending_manifest pending_hash < "$receipt" || { record_issue ERROR "$run" "$galaxy" "$filename" receipt-read "Could not read pending receipt."; return 1; }
        if [[ "$pending_manifest" == forced_* ]] && [ "${AUTO_CHECKSUM:-0}" -eq 0 ]; then
            echo "[$run $galaxy] AWAITING VERIFICATION: $filename was force-uploaded without a snapshot; rerun in default automatic checksum mode to obtain a fresh response."
            return 2
        fi
        if [ "$pending_manifest" = "$manifest_id" ]; then
            echo "[$run $galaxy] AWAITING VERIFICATION: $filename was already attempted against this manifest. Regenerate the CANFAR manifest before proceeding."
            return 2
        fi
    fi
    if [ "$remote_bytes" = "$source_bytes" ] && [ "$remote_hash" = "$source_hash" ]; then
        echo "[$run $galaxy] VERIFIED BY CANFAR MANIFEST: $filename ($source_bytes bytes; SHA256 $source_hash)"
        rm -f "$receipt" || { record_issue ERROR "$run" "$galaxy" "$filename" receipt-cleanup "Could not clear receipt for matching file."; return 1; }
        return 0
    fi

    current_hash=$(sha256sum "$path") || { record_issue ERROR "$run" "$galaxy" "$filename" pre-upload-hash "Could not rehash source before repair."; return 1; }
    current_hash=${current_hash%% *}
    if [ "$current_hash" != "$source_hash" ]; then
        echo "ERROR [$run $galaxy]: source changed while checking $filename; refusing remote changes." >&2
        record_issue ERROR "$run" "$galaxy" "$filename" pre-upload-hash "Source changed while checking; remote changes refused."
        return 1
    fi
    mkdir -p "$CHECKSUM_RECEIPT_DIR" || { record_issue ERROR "$run" "$galaxy" "$filename" receipt-directory "Could not create receipt directory."; return 1; }
    rm -f "$stage_galaxy_dir"/*
    if ! ln "$path" "$staged_path" 2>/dev/null; then
        cp -p "$path" "$staged_path" || { record_issue ERROR "$run" "$galaxy" "$filename" file-staging "Could not link or copy source into staging."; return 1; }
    fi
    if [ "$remote_bytes" != MISSING ]; then
        echo "[$run $galaxy] CHECKSUM MISMATCH: CANFAR=$remote_bytes/$remote_hash Setonix=$source_bytes/$source_hash"
        vmv_cmd=$(resolve_vmv "${stage_galaxy_dir%/*}/vmv_container_wrapper") || {
            echo "ERROR [$run $galaxy]: vmv unavailable and vcp wrapper unsupported; set VMV_CMD." >&2
            record_issue ERROR "$run" "$galaxy" "$filename" quarantine-wrapper "vmv unavailable and vcp wrapper unsupported; set VMV_CMD."
            rm -f "$staged_path"
            return 1
        }
        quarantine_path="${remote_path}.corrupt_$(date -u +%Y%m%dT%H%M%SZ)_${RUN_ID}_${RANDOM}_${RANDOM}"
        echo "[$run $galaxy] QUARANTINE: $remote_path -> $quarantine_path"
        if ! CADC_OVERLAY="$overlay" cadc_with_recovery no-retry "$vmv_cmd" -v "$remote_path" "$quarantine_path"; then
            echo "ERROR [$run $galaxy]: quarantine failed; refusing upload." >&2
            record_issue ERROR "$run" "$galaxy" "$filename" quarantine "Remote move failed; replacement upload refused."
            rm -f "$staged_path"
            return 1
        fi
    fi
    # This snapshot becomes stale as soon as any upload is attempted. Keep a
    # persistent receipt so rerunning it cannot quarantine our new replacement.
    printf '%s %s\n' "$manifest_id" "$source_hash" >"${receipt}.tmp_${RUN_ID}" || { record_issue ERROR "$run" "$galaxy" "$filename" receipt-write "Could not write pending receipt."; return 1; }
    mv "${receipt}.tmp_${RUN_ID}" "$receipt" || { record_issue ERROR "$run" "$galaxy" "$filename" receipt-write "Could not install pending receipt."; return 1; }
    for ((attempt = 1; attempt <= FILE_RETRIES; attempt++)); do
        if [ "${NO_CHECKSUM:-0}" -eq 1 ]; then
            # Recheck before every attempt: a failed transfer may leave a
            # partial destination that vcp could otherwise skip on retry.
            if ! backup_forced_destination "$run" "$galaxy" "$filename" "$overlay" "$remote_path" "${stage_galaxy_dir%/*}/vmv_container_wrapper"; then
                rm -f "$staged_path"
                return 1
            fi
        fi
        echo "[$run $galaxy] UPLOAD: $filename -- attempt $attempt/$FILE_RETRIES"
        if CADC_OVERLAY="$overlay" cadc_with_recovery no-retry "$VCP_CMD" -v "$stage_galaxy_dir" "${dest_base}/"; then
            current_hash=$(sha256sum "$path") || current_hash=unknown
            current_hash=${current_hash%% *}
            rm -f "$staged_path"
            if [ "$current_hash" != "$source_hash" ]; then
                echo "ERROR [$run $galaxy]: source changed during upload; regenerate the manifest and inspect $filename." >&2
                record_issue ERROR "$run" "$galaxy" "$filename" post-upload-source-hash "Source changed or could not be rehashed after upload."
                return 1
            fi
            if [ "${NO_CHECKSUM:-0}" -eq 1 ]; then
                echo "[$run $galaxy] UPLOADED; UNVERIFIED: $filename (CANFAR checksum verification skipped)."
                return 0
            elif [ "${AUTO_CHECKSUM:-0}" -eq 1 ]; then
                echo "[$run $galaxy] UPLOADED; queued for automatic CANFAR verification after this galaxy's uploads."
            else
                echo "[$run $galaxy] UPLOADED; AWAITING VERIFICATION: regenerate the manifest inside CANFAR, then rerun this command. No FITS read-back was performed."
            fi
            return 2
        fi
        if [ "$attempt" -lt "$FILE_RETRIES" ]; then
            delay=$(retry_sleep_seconds "$attempt")
            echo "WARNING [$run $galaxy]: upload failed; retrying in ${delay}s..." >&2
            sleep "$delay"
        fi
    done
    rm -f "$staged_path"
    if [ "${NO_CHECKSUM:-0}" -eq 1 ]; then
        echo "ERROR [$run $galaxy]: forced upload failed; inspect destination and preserved backups before retrying." >&2
        record_issue ERROR "$run" "$galaxy" "$filename" upload "Forced upload failed after $FILE_RETRIES attempt(s); destination unverified; inspect preserved backups."
    else
        echo "ERROR [$run $galaxy]: upload failed; regenerate the CANFAR manifest before another run." >&2
        record_issue ERROR "$run" "$galaxy" "$filename" upload "Upload failed after $FILE_RETRIES attempt(s); obtain fresh checksums before retrying."
    fi
    return 1
}

with_cadc_slot() {
    # mkdir is an atomic lock shared by all galaxy/file subprocesses. Each
    # slot owns its writable overlay, or shares the base mounted read-only.
    local slot slot_lock result
    if [ -z "${TRANSFER_JOBS:-}" ]; then
        "$@"
        return $?
    fi
    while true; do
        for ((slot = 0; slot < TRANSFER_JOBS; slot++)); do
            slot_lock="$worker_dir/slot_${slot}.lock"
            if mkdir "$slot_lock" 2>/dev/null; then
                if (
                    trap 'rmdir "$slot_lock" 2>/dev/null || true' EXIT
                    trap 'exit 130' INT TERM
                    CADC_SLOT_OVERLAY=${worker_overlays[$slot]}
                    "$@"
                ); then result=0; else result=$?; fi
                return "$result"
            fi
        done
        sleep 0.2
    done
}

slot_vcp() {
    local overlay=$1
    shift
    CADC_OVERLAY="${CADC_SLOT_OVERLAY:-$overlay}" cadc_with_recovery retry "$VCP_CMD" -v "$@"
}

slot_job_api() {
    local overlay=$1
    shift
    checksum_job_api "${CADC_SLOT_OVERLAY:-$overlay}" "$@"
}

request_fresh_manifest() {
    local run=$1 galaxy=$2 dest_base=$3 overlay=$4 manifest_path=$5
    shift 5
    local request_path="${manifest_path}.request" request_id started=$SECONDS path
    local remote_worker remote_request job_id job_file job_status=Unknown finished_at=0
    local previous_status= events last_events=-30
    local arguments=()
    for path in "$@"; do arguments+=(--filename "$(basename "$path")"); done
    request_error_stage=request-generation
    request_error_detail="Could not generate checksum request JSON."
    request_id=$(checksum_manifest request --output "$request_path" --run "$run" --galaxy "$galaxy" "${arguments[@]}") || return 1
    if [[ "$dest_base" != arc:* ]]; then
        request_error_stage=destination
        request_error_detail="Automatic checksums require an arc: destination."
        echo "ERROR: automatic checksum jobs currently require arc: destinations." >&2
        return 1
    fi
    remote_worker="$dest_base/.checksum_worker_${request_id}.sh"
    remote_request="$dest_base/.checksum_request_${request_id}.json"
    echo "[$run $galaxy] Requesting fresh checksums INSIDE CANFAR ($request_id)."
    request_error_stage=request-upload
    request_error_detail="Could not upload checksum worker or request to CANFAR."
    if ! with_cadc_slot slot_vcp "$overlay" "$UPLOADER_PATH" "$remote_worker" || ! with_cadc_slot slot_vcp "$overlay" "$request_path" "$remote_request"; then
        echo "ERROR [$run $galaxy]: checksum request upload failed." >&2
        return 1
    fi
    request_error_stage=job-submission
    request_error_detail="CANFAR checksum job submission failed; see API diagnostic."
    job_id=$(with_cadc_slot slot_job_api "$overlay" submit "$CANFAR_CHECKSUM_IMAGE" "/arc/${remote_worker#arc:}" "/arc/${remote_request#arc:}" "checksum-${request_id:0:12}") || return 1
    job_file="$worker_dir/checksum_job_${request_id}.id"
    request_error_stage=job-record
    request_error_detail="Could not save cleanup record for checksum job $job_id."
    printf '%s\n' "$job_id" > "$job_file" || return 1
    echo "[$run $galaxy] Started CANFAR checksum job $job_id."
    request_error_stage=response-timeout
    request_error_detail="Checksum job $job_id produced no valid response within the polling timeout."
    while [ "$((SECONDS - started))" -lt "$CHECKSUM_WAIT_SECONDS" ]; do
        rm -f "$manifest_path"
        if with_cadc_slot slot_vcp "$overlay" "$dest_base/.checksum_response_${request_id}.json" "$manifest_path" >"${manifest_path}.poll_log" 2>&1; then
            if checksum_manifest response-matches --manifest "$manifest_path" --request-id "$request_id"; then
                echo "[$run $galaxy] Received fresh CANFAR-computed checksum response."
                if with_cadc_slot slot_job_api "$overlay" delete "$job_id"; then
                    rm -f "$job_file"
                else
                    echo "WARNING: job $job_id will be cleaned up again at exit." >&2
                fi
                return 0
            fi
        fi
        request_error_stage=job-status
        request_error_detail="Could not read status of checksum job $job_id."
        job_status=$(with_cadc_slot slot_job_api "$overlay" status "$job_id") || break
        if [[ "$job_status" != "$previous_status" ]]; then
            echo "[$run $galaxy] CANFAR checksum job $job_id status: $job_status."
            previous_status=$job_status
        fi
        if [[ "$job_status" == Pending || "$job_status" == Queued ]]; then
            if [ "$((SECONDS - last_events))" -ge 30 ]; then
                last_events=$SECONDS
                events=$(with_cadc_slot slot_job_api "$overlay" events "$job_id" 2>&1) || events=
                if [[ "$events" == *"Back-off pulling image"* || "$events" == *ImagePullBackOff* || "$events" == *ErrImagePull* || "$events" == *"Failed to pull image"* ]]; then
                    request_error_stage=job-image-pull
                    request_error_detail="CANFAR cannot start checksum job $job_id: image pull failed for $CANFAR_CHECKSUM_IMAGE."
                    if [[ "$events" == *"no route to host"* ]]; then
                        request_error_detail+=" CANFAR compute node cannot reach the image registry (no route to host)."
                    fi
                    break
                fi
            fi
        fi
        case "$job_status" in
            Failed|Terminated|Error|failed|terminated)
                # A finished job without a readable response has not verified
                # anything. Preserve diagnostics and fail instead of uploading.
                echo "ERROR [$run $galaxy]: checksum job $job_id ended ($job_status) without a valid response." >&2
                request_error_stage=job-execution
                request_error_detail="Checksum job $job_id ended with status $job_status without a valid response."
                break
                ;;
            Completed|Succeeded|completed|succeeded)
                if [ "$finished_at" -eq 0 ]; then finished_at=$SECONDS; fi
                if [ "$((SECONDS - finished_at))" -ge 30 ]; then
                    echo "ERROR [$run $galaxy]: completed job $job_id has no readable response." >&2
                    request_error_stage=response
                    request_error_detail="Completed checksum job $job_id has no valid readable response."
                    break
                fi
                ;;
        esac
        sleep "$CHECKSUM_POLL_SECONDS"
        request_error_stage=response-timeout
        request_error_detail="Checksum job $job_id produced no valid response within ${CHECKSUM_WAIT_SECONDS}s; last status=$job_status; image=$CANFAR_CHECKSUM_IMAGE."
    done
    # Capture evidence before deleting the job, including failures before its
    # container starts (which have events but no application logs).
    echo "ERROR [$run $galaxy]: $request_error_detail" >&2
    echo "[$run $galaxy] CANFAR job $job_id events:" >&2
    with_cadc_slot slot_job_api "$overlay" events "$job_id" >&2 || true
    echo "[$run $galaxy] CANFAR job $job_id logs:" >&2
    with_cadc_slot slot_job_api "$overlay" logs "$job_id" >&2 || true
    with_cadc_slot slot_job_api "$overlay" delete "$job_id" && rm -f "$job_file"
    echo "ERROR [$run $galaxy]: checksum job produced no verified response before timeout/failure. No verification claimed." >&2
    return 1
}

verify_uploaded_file() {
    local run=$1 galaxy=$2 path=$3 dest_base=$4 manifest_path=$5
    local filename entry remote_bytes remote_hash manifest_id source_bytes source_hash receipt_key
    filename=$(basename "$path")
    verification_error_stage=manifest-lookup
    verification_error_detail="Fresh verification manifest entry is missing, invalid or unreadable."
    entry=$(checksum_manifest lookup --manifest "$manifest_path" --run "$run" --galaxy "$galaxy" --filename "$filename") || return 1
    read -r remote_bytes remote_hash manifest_id <<<"$entry"
    verification_error_stage=source-stat
    verification_error_detail="Could not read source byte count during final verification."
    source_bytes=$(stat -c %s "$path") || return 1
    verification_error_stage=source-hash
    verification_error_detail="Could not hash source during final verification."
    source_hash=$(sha256sum "$path") || return 1
    source_hash=${source_hash%% *}
    if [ "$remote_bytes" != "$source_bytes" ] || [ "$remote_hash" != "$source_hash" ]; then
        echo "ERROR [$run $galaxy]: fresh CANFAR checksum still mismatches $filename; stopping without another repair." >&2
        verification_error_stage=checksum-mismatch
        verification_error_detail="Fresh CANFAR byte count or SHA256 still differs from Setonix; no further repair attempted."
        return 1
    fi
    verification_error_stage=receipt-key
    verification_error_detail="Could not compute receipt key after verification."
    receipt_key=$(printf '%s' "$dest_base/$galaxy/$filename" | sha256sum) || return 1
    rm -f "$CHECKSUM_RECEIPT_DIR/${receipt_key%% *}.receipt" || { verification_error_stage=receipt-cleanup; verification_error_detail="Could not clear verified receipt."; return 1; }
    echo "[$run $galaxy] VERIFIED BY FRESH CANFAR CHECKSUM: $filename ($source_bytes bytes; SHA256 $source_hash)"
    return 0
}

upload_galaxy() {
    local run=$1
    local galaxy=$2
    local overlay=$3
    local source_base dest_base
    local suffix path
    local stage_root stage_galaxy_dir manifest_path file_status
    local sources=()
    local pending_paths=()
    local file_pids=() file_stages=()
    local file_index=0 index file_stage result_file directory_ready=0
    local request_error_stage request_error_detail verification_error_stage verification_error_detail
    local galaxy_status=0
    local succeeded=0
    local failed=0
    local pending=0

    source_base=$(run_source_base "$run")
    dest_base=$(run_dest_base "$run")

    for suffix in "${PRODUCT_SUFFIXES[@]}"; do
        path=$(source_path "$source_base" "$galaxy" "$suffix")
        if [ -f "$path" ]; then
            sources+=("$path")
        else
            echo "WARNING [$run $galaxy]: missing source, skipping: $path" >&2
            if [ "$suffix" = _cont_cube.fits ]; then
                echo "ERROR [$run $galaxy]: required continuum cube is missing." >&2
                record_issue ERROR "$run" "$galaxy" "$(basename "$path")" source-discovery "Required continuum cube is missing on Setonix."
                galaxy_status=1
            fi
        fi
    done

    if [ "${#sources[@]}" -eq 0 ]; then
        echo "ERROR [$run $galaxy]: no requested products found." >&2
        record_issue ERROR "$run" "$galaxy" '(all selected products)' source-discovery "No selected source products were found."
        return 1
    fi

    if [ "${NO_CHECKSUM:-0}" -eq 1 ]; then
        echo "Force uploading $run $galaxy (${#sources[@]} files); CANFAR bytes UNVERIFIED..."
    else
        echo "Comparing $run $galaxy (${#sources[@]} files) against CANFAR-computed checksums..."
    fi

    stage_root=$(mktemp -d "${source_base}/.vcp_upload_stage_${RUN_ID}_${run}_${galaxy}.XXXXXX") || { record_issue ERROR "$run" "$galaxy" '(all selected products)' galaxy-staging "Could not create galaxy staging directory."; return 1; }
    stage_galaxy_dir="${stage_root}/${galaxy}"
    if ! mkdir -p "$stage_galaxy_dir"; then
        record_issue ERROR "$run" "$galaxy" '(all selected products)' galaxy-staging "Could not initialize galaxy staging directory."
        rm -rf "$stage_root"
        return 1
    fi
    manifest_path="$stage_root/checksum_manifest.json"
    if [ "${NO_CHECKSUM:-0}" -eq 1 ]; then
        : # Explicitly skip both automatic jobs and static manifest downloads.
    elif [ "$AUTO_CHECKSUM" -eq 1 ]; then
        if ! request_fresh_manifest "$run" "$galaxy" "$dest_base" "$overlay" "$manifest_path" "${sources[@]}"; then
            for path in "${sources[@]}"; do
                record_issue ERROR "$run" "$galaxy" "$(basename "$path")" "initial-checksum/${request_error_stage:-request-job}" "${request_error_detail:-Fresh CANFAR checksum request failed; no upload attempted.}"
            done
            rm -rf "$stage_root"
            return 1
        fi
    elif ! with_cadc_slot slot_vcp "$overlay" "$dest_base/$CHECKSUM_MANIFEST_NAME" "$manifest_path"; then
        echo "ERROR [$run $galaxy]: manifest unavailable; generate it inside CANFAR before uploading. No FITS download will be attempted." >&2
        for path in "${sources[@]}"; do
            record_issue ERROR "$run" "$galaxy" "$(basename "$path")" initial-checksum/manifest-download "Manual checksum manifest unavailable; no upload attempted."
        done
        rm -rf "$stage_root"
        return 1
    fi

    # Each staging directory contains exactly one product under the original
    # galaxy name. Its adjacent manifest is an immutable copy of the snapshot.
    for path in "${sources[@]}"; do
        file_stage="$stage_root/file_${file_index}"
        if ! mkdir -p "$file_stage/$galaxy" || { [ "${NO_CHECKSUM:-0}" -eq 0 ] && ! cp "$manifest_path" "$file_stage/checksum_manifest.json"; }; then
            echo "ERROR [$run $galaxy]: could not prepare isolated staging for $(basename "$path")." >&2
            for ((index = file_index; index < ${#sources[@]}; index++)); do
                record_issue ERROR "$run" "$galaxy" "$(basename "${sources[$index]}")" file-staging "Staging failed; this file was not launched."
            done
            galaxy_status=1
            break
        fi
        result_file="$file_stage/result"
        (
            if with_cadc_slot upload_one_file "$run" "$galaxy" "$path" "$dest_base" "$overlay" "$file_stage/$galaxy"; then
                printf '0\n' > "$result_file"
            else
                printf '%s\n' "$?" > "$result_file"
            fi
        ) &
        file_pids+=("$!")
        file_stages+=("$file_stage")
        # vcp creates a missing galaxy directory. Establish it with one file
        # before concurrent uploads, avoiding CADC DuplicateNode races. If
        # that file fails, try the next sequentially until one succeeds.
        if [ "$directory_ready" -eq 0 ]; then
            if wait "${file_pids[$file_index]}" && read -r file_status < "$result_file"; then
                if [ "$file_status" -eq 0 ] || [ "$file_status" -eq 2 ]; then
                    directory_ready=1
                fi
            fi
        fi
        file_index=$((file_index + 1))
    done

    # Wait for every file before refreshing the galaxy's checksum manifest.
    for ((index = 0; index < ${#file_pids[@]}; index++)); do
        path=${sources[$index]}
        file_status=1
        if wait "${file_pids[$index]}" && [ -f "${file_stages[$index]}/result" ]; then
            read -r file_status < "${file_stages[$index]}/result" || file_status=1
        else
            record_issue ERROR "$run" "$galaxy" "$(basename "$path")" file-worker "File task terminated or produced no result."
        fi
        if [ "$file_status" -eq 0 ]; then
            succeeded=$((succeeded + 1))
        else
            if [ "$file_status" -eq 2 ] && [ "$AUTO_CHECKSUM" -eq 1 ]; then
                pending_paths+=("$path")
            fi
            if [ "$file_status" -eq 2 ]; then
                pending=$((pending + 1))
                if [ "$galaxy_status" -eq 0 ]; then galaxy_status=2; fi
            else
                failed=$((failed + 1))
                galaxy_status=1
            fi
        fi
    done

    # Verify all newly uploaded products in one job, after uploads finish.
    if [ "$AUTO_CHECKSUM" -eq 1 ] && [ "${#pending_paths[@]}" -gt 0 ]; then
        if request_fresh_manifest "$run" "$galaxy" "$dest_base" "$overlay" "$manifest_path" "${pending_paths[@]}"; then
            for path in "${pending_paths[@]}"; do
                pending=$((pending - 1))
                if verify_uploaded_file "$run" "$galaxy" "$path" "$dest_base" "$manifest_path"; then
                    succeeded=$((succeeded + 1))
                else
                    record_issue ERROR "$run" "$galaxy" "$(basename "$path")" "post-upload-verification/${verification_error_stage:-verification}" "${verification_error_detail:-Final checksum verification failed.}"
                    failed=$((failed + 1))
                    galaxy_status=1
                fi
            done
            if [ "$galaxy_status" -eq 2 ]; then galaxy_status=0; fi
        else
            for path in "${pending_paths[@]}"; do
                record_issue ERROR "$run" "$galaxy" "$(basename "$path")" "post-upload-checksum/${request_error_stage:-request-job}" "${request_error_detail:-Fresh post-upload checksum request failed; uploaded file remains unverified.}"
            done
            failed=$((failed + pending))
            pending=0
            galaxy_status=1
        fi
    fi
    if [ "$AUTO_CHECKSUM" -eq 0 ] && [ "$pending" -gt 0 ]; then
        for ((index = 0; index < ${#file_stages[@]}; index++)); do
            read -r file_status < "${file_stages[$index]}/result" || continue
            if [ "$file_status" -eq 2 ]; then
                record_issue PENDING "$run" "$galaxy" "$(basename "${sources[$index]}")" verification "Obtain fresh checksums and rerun; previously force-uploaded files require default automatic checksum mode."
            fi
        done
    fi
    rm -rf "$stage_root"

    if [ "${NO_CHECKSUM:-0}" -eq 1 ]; then
        echo "Finished $run $galaxy: $succeeded uploaded UNVERIFIED, $failed failed; CANFAR checksum verification skipped."
    elif [ "$galaxy_status" -eq 0 ]; then
        echo "Finished $run $galaxy: ${succeeded}/${#sources[@]} files verified by CANFAR manifest."
    else
        echo "Finished $run $galaxy: $succeeded verified, $pending awaiting fresh manifest, $failed failed." >&2
    fi

    return "$galaxy_status"
}

pids=()
for ((i = 0; i < EFFECTIVE_JOBS; i++)); do
    part_file="$worker_dir/part_${i}.txt"
    worker_overlay="${worker_overlays[$((i % TRANSFER_JOBS))]}"

    (
        worker_status=0
        while IFS= read -r work_item; do
            run=${work_item%%|*}
            galaxy=${work_item#*|}
            if upload_galaxy "$run" "$galaxy" "$worker_overlay"; then
                :
            else
                galaxy_result=$?
                if [ "$galaxy_result" -eq 1 ] || [ "$worker_status" -eq 0 ]; then
                    worker_status=$galaxy_result
                fi
            fi
        done < "$part_file"
        exit "$worker_status"
    ) &

    pids+=("$!")
done

status=0
for pid in "${pids[@]}"; do
    if wait "$pid"; then
        :
    else
        worker_result=$?
        if [ "$worker_result" -eq 1 ] || [ "$status" -eq 0 ]; then status=$worker_result; fi
    fi
done

runtime_secs=$((SECONDS - start_secs))
echo "Total runtime: $(format_runtime "$runtime_secs") (${runtime_secs}s)"

if [ "$status" -eq 2 ]; then
    echo "Uploads attempted; awaiting a fresh CANFAR-generated manifest. Regenerate it and rerun to verify. No FITS products were downloaded."
    exit 2
elif [ "$status" -ne 0 ]; then
    echo "One or more files/galaxies failed. Check manifest coverage, source files, and transfer errors." >&2
    exit "$status"
fi

if [ "$NO_CHECKSUM" -eq 1 ]; then
    echo "All selected available file transfers completed; destination bytes UNVERIFIED. CANFAR checksum verification skipped; no FITS read-back downloads."
else
    echo "All selected available files verified by CANFAR-computed manifest; no FITS read-back downloads."
fi
