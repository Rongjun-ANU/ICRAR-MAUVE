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

target_dir="${TARGET_DIR:-/arc/projects/mauve/cubes/v3tk}"
phangs_native_vos_dir="vos:phangs/RELEASES/PHANGS-MUSE/DR1.0/DATACUBES"
phangs_native_files=(
	NGC4254_PHANGS_DATACUBE_native.fits
	NGC4321_PHANGS_DATACUBE_native.fits
	NGC4535_PHANGS_DATACUBE_native.fits
)

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

is_phangs_native_galid() {
	local candidate="$1"
	local phangs_file phangs_galid
	for phangs_file in "${phangs_native_files[@]}"; do
		phangs_galid="${phangs_file%%_PHANGS_DATACUBE_native.fits}"
		if [[ "$candidate" == "$phangs_galid" ]]; then
			return 0
		fi
	done
	return 1
}

galid_from_cube_name() {
	local name="$1"
	name="${name%.gz}"
	if [[ "$name" == *_PHANGS_DATACUBE_native.fits ]]; then
		printf '%s\n' "${name%%_PHANGS_DATACUBE_native.fits}"
	elif [[ "$name" == *_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk.fits ]]; then
		printf '%s\n' "${name%%_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk.fits}"
	else
		printf '%s\n' "${name%%_*}"
	fi
}

seen_input_stems=()
files=()

append_source_once() {
	local source_path="$1"
	local source_name="$2"
	local stem_key already_seen seen_stem
	stem_key="${source_name%.gz}"
	already_seen=0
	for seen_stem in "${seen_input_stems[@]}"; do
		if [[ "$seen_stem" == "$stem_key" ]]; then
			already_seen=1
			break
		fi
	done
	if [[ "$already_seen" -eq 1 ]]; then
		return
	fi
	seen_input_stems+=("$stem_key")
	files+=("$source_path")
}

shopt -s nullglob
if [[ -d "$target_dir" ]]; then
	raw_files=("$target_dir"/*_v3tk.fits "$target_dir"/*_v3tk.fits.gz)
else
	echo "WARNING: target directory does not exist, skipping local v3tk discovery: $target_dir" >&2
	raw_files=()
fi

for candidate in "${raw_files[@]}"; do
	base_candidate="$(basename "$candidate")"
	if is_phangs_native_galid "$(galid_from_cube_name "$base_candidate")"; then
		continue
	fi
	append_source_once "$candidate" "$base_candidate"
done

for phangs_file in "${phangs_native_files[@]}"; do
	local_phangs="$target_dir/$phangs_file"
	if [[ -f "$local_phangs" ]]; then
		append_source_once "$local_phangs" "$phangs_file"
	else
		append_source_once "$phangs_native_vos_dir/$phangs_file" "$phangs_file"
	fi
done

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
	echo "No input cubes found from local v3tk patterns or PHANGS native public sources"
	exit 0
fi

echo "Found ${#files[@]} file(s):"
for f in "${files[@]}"; do
	echo "- $f"
done

total_start_epoch="$(date +%s)"

for src_input in "${files[@]}"; do
	per_start_epoch="$(date +%s)"
	base_input="$(basename "$src_input")"
	dest_input="$(pwd)/$base_input"

	echo ""
	echo "=== Processing: $src_input ==="

	# 0) conda activate ICRAR (called per-file; no-op if already active)
	conda activate ICRAR

	# 1) copy to pwd (prefer rsync when available)
	rm -f "$dest_input"
	if [[ "$src_input" == vos:* ]]; then
		if ! command -v vcp >/dev/null 2>&1; then
			echo "ERROR: vcp is required to stage public VOSpace input: $src_input" >&2
			exit 1
		fi
		vcp "$src_input" "$dest_input"
	elif command -v rsync >/dev/null 2>&1; then
		rsync -a "$src_input" "$dest_input"
	else
		cp -f "$src_input" "$dest_input"
	fi

	echo "Local input: $dest_input"
	ls -lh "$dest_input" || true

	# 2) run conversion; Astropy reads both .fits and .fits.gz directly.
	echo "Running: python v3tk_to_VRI.py $dest_input"
	set +e
	if [[ -n "$time_cmd" && "$time_supports_verbose" -eq 1 ]]; then
		"$time_cmd" -v python v3tk_to_VRI.py "$dest_input"
	elif [[ -n "$time_cmd" ]]; then
		"$time_cmd" -p python v3tk_to_VRI.py "$dest_input"
	else
		python v3tk_to_VRI.py "$dest_input"
	fi
	py_status=$?
	set -e
	if [[ $py_status -ne 0 ]]; then
		echo "ERROR: python failed for $dest_input (exit status: $py_status)" >&2
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

	# 3) cleanup copied input before moving on
	rm -f "$dest_input"

	per_end_epoch="$(date +%s)"
	per_runtime="$((per_end_epoch - per_start_epoch))"
	echo "Done: $src_input"
	echo "Runtime (this file): ${per_runtime}s"
done

total_end_epoch="$(date +%s)"
total_runtime="$((total_end_epoch - total_start_epoch))"
echo ""
echo "All done. Total runtime: ${total_runtime}s"
echo "End time: $(date -Is)"
