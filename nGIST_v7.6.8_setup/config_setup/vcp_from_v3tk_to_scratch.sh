#!/usr/bin/env bash
set -euo pipefail

start_secs=$SECONDS
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

format_runtime() {
	local total=$1
	printf '%02d:%02d:%02d' $((total / 3600)) $(((total % 3600) / 60)) $((total % 60))
}

usage() {
	cat >&2 <<'EOF'
Usage: ./vcp_from_v3tk_to_scratch.sh [--dry-run] [GALID ...]

With no GALID arguments, transfer the full MAUVE v3tk manifest plus the three
supported PHANGS native cubes. With GALID arguments, transfer only those
galaxies. Use --dry-run to print the selected sources and worker count without
contacting CADC or copying files.
EOF
}

BASE_OVERLAY=/software/projects/pawsey1308/containers/cadc_overlay.img
OVERLAY_DIR=/software/projects/pawsey1308/containers
SRC_GLOB='arc:projects/mauve/cubes/v3tk/*v3tk.fits.gz'
SRC_PREFIX='arc:projects/mauve/cubes/v3tk/'
V3TK_CUBE_SUFFIX='_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk.fits.gz'
PHANGS_NATIVE_VOS_DIR='vos:phangs/RELEASES/PHANGS-MUSE/DR1.0/DATACUBES'
PHANGS_NATIVE_GALIDS=(
	NGC4254
	NGC4321
	NGC4535
)
PHANGS_NATIVE_SOURCES=(
	'vos:phangs/RELEASES/PHANGS-MUSE/DR1.0/DATACUBES/NGC4254_PHANGS_DATACUBE_native.fits'
	'vos:phangs/RELEASES/PHANGS-MUSE/DR1.0/DATACUBES/NGC4321_PHANGS_DATACUBE_native.fits'
	'vos:phangs/RELEASES/PHANGS-MUSE/DR1.0/DATACUBES/NGC4535_PHANGS_DATACUBE_native.fits'
)
DEST_DIR=/scratch/pawsey1308/mauve/cubes/v3tk/
RUN_ID=$$

# Setonix project quota cannot tolerate many overlay copies. Keep at most five
# concurrent workers, and later shrink further when fewer sources are selected.
MAX_JOBS=5
JOBS=${JOBS:-5}
if [[ ! "$JOBS" =~ ^[0-9]+$ ]]; then
	echo "ERROR: JOBS must be a positive integer, got: $JOBS" >&2
	exit 2
fi
if [ "$JOBS" -lt 1 ]; then
	JOBS=1
fi
if [ "$JOBS" -gt "$MAX_JOBS" ]; then
	echo "Capping JOBS from $JOBS to $MAX_JOBS for Setonix overlay quota safety."
	JOBS=$MAX_JOBS
fi

DRY_RUN=0
REQUESTED_GALIDS=()
while [ "$#" -gt 0 ]; do
	case "$1" in
		-n|--dry-run)
			DRY_RUN=1
			;;
		-h|--help)
			usage
			exit 0
			;;
		--)
			shift
			while [ "$#" -gt 0 ]; do
				REQUESTED_GALIDS+=("$1")
				shift
			done
			break
			;;
		-*)
			echo "ERROR: unknown option: $1" >&2
			usage
			exit 2
			;;
		*)
			REQUESTED_GALIDS+=("$1")
			;;
	esac
	shift
done

if [ -f "$SCRIPT_DIR/27_galaxies.sh" ]; then
	# Reuse the canonical local galaxy allowlist when this script is run from config_setup.
	# If copied alone to Setonix home, the regex fallback below still supports selected IDs.
	# shellcheck disable=SC1091
	. "$SCRIPT_DIR/27_galaxies.sh"
else
	looks_like_galid() {
		local candidate="$1"
		[[ "$candidate" =~ ^(IC|NGC)[0-9][0-9_]*$ ]]
	}

	is_known_galid() {
		looks_like_galid "$1"
	}
fi

is_phangs_native_galid() {
	local candidate="$1"
	local phangs_galid
	for phangs_galid in "${PHANGS_NATIVE_GALIDS[@]}"; do
		if [ "$candidate" = "$phangs_galid" ]; then
			return 0
		fi
	done
	return 1
}

source_for_galid() {
	local galid="$1"
	if is_phangs_native_galid "$galid"; then
		printf '%s/%s_PHANGS_DATACUBE_native.fits\n' "$PHANGS_NATIVE_VOS_DIR" "$galid"
	else
		printf '%s%s%s\n' "$SRC_PREFIX" "$galid" "$V3TK_CUBE_SUFFIX"
	fi
}

manifest=""
part_dir=""

# On Ctrl+C, stop child processes and wait for lock holders to exit cleanly.
cleanup_on_interrupt() {
	echo "Interrupted. Stopping active CADC transfers..."
	pkill -TERM -P $$ 2>/dev/null || true
	sleep 1
	pkill -KILL -P $$ 2>/dev/null || true
	wait 2>/dev/null || true
	if [ -n "${part_dir:-}" ] && [ -d "$part_dir" ]; then
		rm -rf "$part_dir"
	fi
	if [ -n "${manifest:-}" ] && [ -f "$manifest" ]; then
		rm -f "$manifest"
	fi
	rm -f "$OVERLAY_DIR"/cadc_overlay_worker_"$RUN_ID"_*.img 2>/dev/null || true
	exit 130
}
trap cleanup_on_interrupt INT TERM

# Cleanup should also run on normal exit.
cleanup_on_exit() {
	if [ -n "${part_dir:-}" ] && [ -d "$part_dir" ]; then
		rm -rf "$part_dir"
	fi
	if [ -n "${manifest:-}" ] && [ -f "$manifest" ]; then
		rm -f "$manifest"
	fi
	rm -f "$OVERLAY_DIR"/cadc_overlay_worker_"$RUN_ID"_*.img 2>/dev/null || true
}
trap cleanup_on_exit EXIT

# If previous runs were interrupted, stop stale transfer processes first.
pkill -TERM -u "$USER" -f '/cadcenv/bin/vcp|/home/rhuang/bin/vcp|/cadcenv/bin/vls' 2>/dev/null || true
sleep 1
pkill -KILL -u "$USER" -f '/cadcenv/bin/vcp|/home/rhuang/bin/vcp|/cadcenv/bin/vls' 2>/dev/null || true

# Wait briefly if the overlay is still locked by a previous process.
wait_for_overlay() {
	local overlay=$1
	local tries=0
	while ! flock -n "$overlay" -c true 2>/dev/null; do
		tries=$((tries + 1))
		if [ "$tries" -ge 30 ]; then
			echo "Overlay still busy after 30s: $overlay"
			return 1
		fi
		sleep 1
	done
}

manifest=$(mktemp)

if [ "${#REQUESTED_GALIDS[@]}" -gt 0 ]; then
	for galid in "${REQUESTED_GALIDS[@]}"; do
		if ! is_known_galid "$galid"; then
			echo "ERROR: unknown galaxy ID: $galid" >&2
			exit 2
		fi
		source_for_galid "$galid" >> "$manifest"
	done
else
	if [ "$DRY_RUN" -eq 0 ]; then
		cadc-get-cert -u RongjunHuang
	fi
	# Build full source manifest once.
	vls "$SRC_GLOB" | sed "s#^#$SRC_PREFIX#" > "$manifest"
	printf '%s\n' "${PHANGS_NATIVE_SOURCES[@]}" >> "$manifest"
fi

source_count=$(wc -l < "$manifest" | tr -d ' ')
if [ "$source_count" -eq 0 ]; then
	echo "No sources selected."
	exit 0
fi

EFFECTIVE_JOBS="$JOBS"
if [ "$source_count" -lt "$EFFECTIVE_JOBS" ]; then
	EFFECTIVE_JOBS="$source_count"
fi

if [ "$DRY_RUN" -eq 1 ]; then
	if [ "${#REQUESTED_GALIDS[@]}" -gt 0 ]; then
		echo "Requested galaxies: ${REQUESTED_GALIDS[*]}"
	else
		echo "Requested galaxies: ALL"
	fi
	echo "Sources selected: $source_count"
	echo "Effective workers: $EFFECTIVE_JOBS"
	echo "Sources:"
	cat "$manifest"
	exit 0
fi

if [ "${#REQUESTED_GALIDS[@]}" -gt 0 ]; then
	cadc-get-cert -u RongjunHuang
fi

wait_for_overlay "$BASE_OVERLAY" || exit 1

# Prepare per-worker overlays so each concurrent vcp process has its own lock.
for i in $(seq 0 $((EFFECTIVE_JOBS - 1))); do
	worker_overlay="$OVERLAY_DIR/cadc_overlay_worker_${RUN_ID}_${i}.img"
	cp --reflink=auto "$BASE_OVERLAY" "$worker_overlay"
	wait_for_overlay "$worker_overlay" || exit 1
done

part_dir=$(mktemp -d)

# Split work round-robin across workers.
awk -v jobs="$EFFECTIVE_JOBS" -v dir="$part_dir" '{print > (dir "/part_" ((NR-1)%jobs) ".txt")}' "$manifest"

# Start one sequential worker per overlay.
for i in $(seq 0 $((EFFECTIVE_JOBS - 1))); do
	part_file="$part_dir/part_${i}.txt"
	worker_overlay="$OVERLAY_DIR/cadc_overlay_worker_${RUN_ID}_${i}.img"
	(
		if [ -f "$part_file" ]; then
			while IFS= read -r src; do
				CADC_OVERLAY="$worker_overlay" vcp "$src" "$DEST_DIR"
			done < "$part_file"
		fi
	) &
done

wait

cd "$DEST_DIR"
shopt -s nullglob
for f in *v3tk.fits.gz; do
	gunzip -kf "$f"
done

runtime_secs=$((SECONDS - start_secs))
echo "Total runtime: $(format_runtime "$runtime_secs") (${runtime_secs}s)"
