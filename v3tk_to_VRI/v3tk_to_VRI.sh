#!/usr/bin/env bash
set -euo pipefail

log_file="$(pwd)/v3tk_to_VRI.log"
if [[ "${V3TK_TO_VRI_LOG_CAPTURED:-0}" != "1" ]]; then
	: > "$log_file"
	export V3TK_TO_VRI_LOG_CAPTURED=1
	"${BASH:-bash}" "$0" "$@" 2>&1 | tee -a "$log_file"
	exit "${PIPESTATUS[0]}"
fi
echo "Logging to: $log_file"

timestamp_now() {
	date -u +"%Y-%m-%dT%H:%M:%SZ"
}

echo "Start time: $(timestamp_now)"
echo "Host: $(hostname)"
echo "User: $(whoami)"
echo "CWD: $(pwd)"
echo "Uname: $(uname -a)"
echo ""
echo "Resource snapshot (best effort):"
ulimit -a || true
command -v free >/dev/null 2>&1 && free -h || true
df -h . || true
command -v quota >/dev/null 2>&1 && quota -s 2>/dev/null || true
echo ""

target_dir="${TARGET_DIR:-/arc/projects/mauve/cubes/v3tk}"
pwd_dir="$(pwd)"
phangs_native_vos_dir="vos:phangs/RELEASES/PHANGS-MUSE/DR1.0/DATACUBES"
phangs_native_files=(
	NGC4254_PHANGS_DATACUBE_native.fits
	NGC4321_PHANGS_DATACUBE_native.fits
	NGC4535_PHANGS_DATACUBE_native.fits
)

usage() {
	cat <<'USAGE'
Usage:
  ./v3tk_to_VRI.sh [--dry-run] [GALID ...]

Examples:
  ./v3tk_to_VRI.sh
  ./v3tk_to_VRI.sh NGC4254 NGC4321 NGC4535
  ./v3tk_to_VRI.sh --dry-run NGC4254

Without GALID arguments, all local v3tk cubes plus the supported PHANGS-native
public cubes are selected. With GALID arguments, only those galaxies are staged
and converted.
USAGE
}

dry_run=0
while [[ $# -gt 0 ]]; do
	case "$1" in
		--dry-run|-n)
			dry_run=1
			shift
			;;
		-h|--help)
			usage
			exit 0
			;;
		--)
			shift
			break
			;;
		-*)
			echo "ERROR: unknown option: $1" >&2
			usage >&2
			exit 2
			;;
		*)
			break
			;;
	esac
done
requested_galids=("$@")

normalize_galid() {
	printf '%s' "$1" | tr '[:lower:]' '[:upper:]' | tr -d ' '
}

normalized_requested_galids=()
for requested_galid in "${requested_galids[@]}"; do
	normalized_requested_galids+=("$(normalize_galid "$requested_galid")")
done

if [[ "$dry_run" -eq 1 ]]; then
	echo "Dry run: no files will be copied, converted, or removed."
fi

# Needed for 'conda activate' in non-interactive shells during real runs.
activate_conda() {
	if ! command -v conda >/dev/null 2>&1; then
		echo "ERROR: 'conda' not found in PATH. Load conda first, then re-run." >&2
		exit 1
	fi

	conda_base="$(conda info --base)"
	if [[ ! -f "$conda_base/etc/profile.d/conda.sh" ]]; then
		echo "ERROR: conda init script not found: $conda_base/etc/profile.d/conda.sh" >&2
		exit 1
	fi

	source "$conda_base/etc/profile.d/conda.sh"
}

phangs_native_file_for_galid() {
	local candidate="$1"
	local phangs_file phangs_galid
	for phangs_file in "${phangs_native_files[@]}"; do
		phangs_galid="${phangs_file%%_PHANGS_DATACUBE_native.fits}"
		if [[ "$candidate" == "$phangs_galid" ]]; then
			printf '%s\n' "$phangs_file"
			return 0
		fi
	done
	return 1
}

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

local_v3tk_source_for_galid_in_dir() {
	local galid="$1"
	local search_dir="$2"
	local candidate
	local matches=()

	for candidate in \
		"$search_dir/${galid}_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk.fits" \
		"$search_dir/${galid}_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk.fits.gz"
	do
		if [[ -f "$candidate" ]]; then
			printf '%s\n' "$candidate"
			return 0
		fi
	done

	if [[ -d "$search_dir" ]]; then
		matches=(
			"$search_dir/${galid}_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk"*.fits
			"$search_dir/${galid}_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk"*.fits.gz
		)
		for candidate in "${matches[@]}"; do
			if [[ -f "$candidate" ]]; then
				printf '%s\n' "$candidate"
				return 0
			fi
		done
	fi

	return 1
}

local_v3tk_source_for_galid() {
	local galid="$1"
	local local_source

	if local_source="$(local_v3tk_source_for_galid_in_dir "$galid" "$pwd_dir")"; then
		printf '%s\n' "$local_source"
		return 0
	fi

	if local_source="$(local_v3tk_source_for_galid_in_dir "$galid" "$target_dir")"; then
		printf '%s\n' "$local_source"
		return 0
	fi

	return 1
}

source_for_galid() {
	local galid="$1"
	local phangs_file local_phangs local_source

	if phangs_file="$(phangs_native_file_for_galid "$galid")"; then
		local_phangs="$target_dir/$phangs_file"
		if [[ -f "$local_phangs" ]]; then
			printf '%s\n' "$local_phangs"
		else
			printf '%s\n' "$phangs_native_vos_dir/$phangs_file"
		fi
		return 0
	fi

	if local_source="$(local_v3tk_source_for_galid "$galid")"; then
		printf '%s\n' "$local_source"
		return 0
	fi

	echo "ERROR: no local v3tk cube found for selected GALID $galid under $pwd_dir or $target_dir" >&2
	return 1
}

seen_input_stems=()
files=()

append_source_once() {
	local source_path="$1"
	local source_name="$2"
	local stem_key already_seen seen_stem
	stem_key="${source_name%.gz}"
	already_seen=0
	if (( ${#seen_input_stems[@]} > 0 )); then
		for seen_stem in "${seen_input_stems[@]}"; do
			if [[ "$seen_stem" == "$stem_key" ]]; then
				already_seen=1
				break
			fi
		done
	fi
	if [[ "$already_seen" -eq 1 ]]; then
		return
	fi
	seen_input_stems+=("$stem_key")
	files+=("$source_path")
}

shopt -s nullglob
if (( ${#normalized_requested_galids[@]} > 0 )); then
	for requested_galid in "${normalized_requested_galids[@]}"; do
		source_path="$(source_for_galid "$requested_galid")"
		append_source_once "$source_path" "$(basename "$source_path")"
	done
else
	raw_files=("$pwd_dir"/*_v3tk.fits "$pwd_dir"/*_v3tk.fits.gz)

	for candidate in "${raw_files[@]}"; do
		base_candidate="$(basename "$candidate")"
		if is_phangs_native_galid "$(galid_from_cube_name "$base_candidate")"; then
			continue
		fi
		append_source_once "$candidate" "$base_candidate"
	done

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
fi

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

if [[ "$dry_run" -eq 1 ]]; then
	echo "Dry run complete."
	exit 0
fi

activate_conda

total_start_epoch="$(date +%s)"

for src_input in "${files[@]}"; do
	per_start_epoch="$(date +%s)"
	base_input="$(basename "$src_input")"
	output_base="${base_input%.gz}"
	output_path="$(pwd)/${output_base%.fits}_VRI.fits"
	run_input="$src_input"

	echo ""
	echo "=== Processing: $src_input ==="

	# 0) conda activate ICRAR (called per-file; no-op if already active)
	conda activate ICRAR

	# 1) use local filesystem inputs in place; Astropy reads both .fits and .fits.gz directly.
	if [[ "$src_input" == vos:* ]]; then
		run_input="$(pwd)/$base_input"
		if ! command -v vcp >/dev/null 2>&1; then
			echo "ERROR: vcp is required to stage public VOSpace input: $src_input" >&2
			exit 1
		fi
		vcp "$src_input" "$run_input"
	fi

	echo "Input: $run_input"
	ls -lh "$run_input" || true

	# 2) run conversion; Astropy reads both .fits and .fits.gz directly.
	echo "Running: python v3tk_to_VRI.py $run_input --output $output_path"
	set +e
	if [[ -n "$time_cmd" && "$time_supports_verbose" -eq 1 ]]; then
		"$time_cmd" -v python v3tk_to_VRI.py "$run_input" --output "$output_path"
	elif [[ -n "$time_cmd" ]]; then
		"$time_cmd" -p python v3tk_to_VRI.py "$run_input" --output "$output_path"
	else
		python v3tk_to_VRI.py "$run_input" --output "$output_path"
	fi
	py_status=$?
	set -e
	if [[ $py_status -ne 0 ]]; then
		echo "ERROR: python failed for $run_input (exit status: $py_status)" >&2
		if [[ $py_status -eq 137 || $py_status -eq 9 ]]; then
			echo "HINT: Exit status $py_status usually means the process was SIGKILL'ed." >&2
			echo "      Most common cause is out-of-memory (OOM) or a memory/cgroup limit from the system/scheduler." >&2
			echo "      Check: available memory (free -h), job mem limits, and whether this run is inside a batch job." >&2
		fi
		echo "Resource snapshot after failure (best effort):" >&2
		command -v free >/dev/null 2>&1 && free -h || true
		df -h . || true
		command -v quota >/dev/null 2>&1 && quota -s 2>/dev/null || true
		exit $py_status
	fi

	per_end_epoch="$(date +%s)"
	per_runtime="$((per_end_epoch - per_start_epoch))"
	echo "Done: $src_input"
	echo "Runtime (this file): ${per_runtime}s"
done

total_end_epoch="$(date +%s)"
total_runtime="$((total_end_epoch - total_start_epoch))"
echo ""
echo "All done. Total runtime: ${total_runtime}s"
echo "End time: $(timestamp_now)"
