#!/usr/bin/env bash
set -euo pipefail

log_file="$(pwd)/v3tk_to_VRI.log"
: > "$log_file"
exec > >(tee -a "$log_file") 2>&1
echo "Logging to: $log_file"

echo "Start time: $(date -Is)"
echo "Host: $(hostname)"
echo "User: $(whoami)"
echo "CWD: $(pwd)"
echo "Uname: $(uname -a)"
echo ""
echo "Resource snapshot (best effort):"
ulimit -a || true
command -v free >/dev/null 2>&1 && free -h || true
df -h . || true
command -v quota >/dev/null 2>&1 && quota -s || true
echo ""

target_dir="/arc/projects/mauve/cubes/v3tk"

if [[ ! -d "$target_dir" ]]; then
	echo "ERROR: target directory does not exist: $target_dir" >&2
	exit 1
fi

if ! command -v conda >/dev/null 2>&1; then
	echo "ERROR: 'conda' not found in PATH. Load conda first, then re-run." >&2
	exit 1
fi

conda_base="$(conda info --base)"
if [[ ! -f "$conda_base/etc/profile.d/conda.sh" ]]; then
	echo "ERROR: conda init script not found: $conda_base/etc/profile.d/conda.sh" >&2
	exit 1
fi

# Needed for 'conda activate' in non-interactive shells
source "$conda_base/etc/profile.d/conda.sh"

shopt -s nullglob
files=("$target_dir"/*_v3tk.fits.gz)

time_cmd=""
time_supports_verbose=0
if command -v /usr/bin/time >/dev/null 2>&1; then
	time_cmd="/usr/bin/time"
	if /usr/bin/time -v true >/dev/null 2>&1; then
		time_supports_verbose=1
	fi
fi

echo "Target directory: $target_dir"
echo "Working directory: $(pwd)"

if (( ${#files[@]} == 0 )); then
	echo "No matches found for: $target_dir/*_v3tk.fits.gz"
	exit 0
fi

echo "Found ${#files[@]} file(s):"
for f in "${files[@]}"; do
	echo "- $f"
done

total_start_epoch="$(date +%s)"

for src_gz in "${files[@]}"; do
	per_start_epoch="$(date +%s)"
	base_gz="$(basename "$src_gz")"
	dest_gz="$(pwd)/$base_gz"
	dest_fits="${dest_gz%.gz}"

	echo ""
	echo "=== Processing: $src_gz ==="

	# 0) conda activate ICRAR (called per-file; no-op if already active)
	conda activate ICRAR

	# 1) copy to pwd (prefer rsync when available)
	rm -f "$dest_gz" "$dest_fits"
	if command -v rsync >/dev/null 2>&1; then
		rsync -a "$src_gz" "$dest_gz"
	else
		cp -f "$src_gz" "$dest_gz"
	fi

	# 2) unzip to .fits and remove the .gz in pwd
	gunzip -f "$dest_gz"

	echo "Local input: $dest_fits"
	ls -lh "$dest_fits" || true

	# 3) run conversion
	echo "Running: python v3tk_to_VRI.py $dest_fits"
	set +e
	if [[ -n "$time_cmd" && "$time_supports_verbose" -eq 1 ]]; then
		"$time_cmd" -v python v3tk_to_VRI.py "$dest_fits"
	elif [[ -n "$time_cmd" ]]; then
		"$time_cmd" -p python v3tk_to_VRI.py "$dest_fits"
	else
		python v3tk_to_VRI.py "$dest_fits"
	fi
	py_status=$?
	set -e
	if [[ $py_status -ne 0 ]]; then
		echo "ERROR: python failed for $dest_fits (exit status: $py_status)" >&2
		if [[ $py_status -eq 137 || $py_status -eq 9 ]]; then
			echo "HINT: Exit status $py_status usually means the process was SIGKILL'ed." >&2
			echo "      Most common cause is out-of-memory (OOM) or a memory/cgroup limit from the system/scheduler." >&2
			echo "      Check: available memory (free -h), job mem limits, and whether this run is inside a batch job." >&2
		fi
		echo "Resource snapshot after failure (best effort):" >&2
		command -v free >/dev/null 2>&1 && free -h || true
		df -h . || true
		command -v quota >/dev/null 2>&1 && quota -s || true
		exit $py_status
	fi

	# 4) cleanup input .fits before moving on
	rm -f "$dest_fits"

	per_end_epoch="$(date +%s)"
	per_runtime="$((per_end_epoch - per_start_epoch))"
	echo "Done: $src_gz"
	echo "Runtime (this file): ${per_runtime}s"
done

total_end_epoch="$(date +%s)"
total_runtime="$((total_end_epoch - total_start_epoch))"
echo ""
echo "All done. Total runtime: ${total_runtime}s"
echo "End time: $(date -Is)"
