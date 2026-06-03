#!/usr/bin/env bash
set -u

start_secs=$SECONDS

format_runtime() {
	local total=$1
	printf '%02d:%02d:%02d' $((total / 3600)) $(((total % 3600) / 60)) $((total % 60))
}

BASE_OVERLAY=/software/projects/pawsey1308/containers/cadc_overlay.img
OVERLAY_DIR=/software/projects/pawsey1308/containers
SRC_GLOB='arc:projects/mauve/cubes/v3tk/*v3tk.fits.gz'
SRC_PREFIX='arc:projects/mauve/cubes/v3tk/'
PHANGS_NATIVE_SOURCES=(
	'vos:phangs/RELEASES/PHANGS-MUSE/DR1.0/DATACUBES/NGC4254_PHANGS_DATACUBE_native.fits'
	'vos:phangs/RELEASES/PHANGS-MUSE/DR1.0/DATACUBES/NGC4321_PHANGS_DATACUBE_native.fits'
	'vos:phangs/RELEASES/PHANGS-MUSE/DR1.0/DATACUBES/NGC4535_PHANGS_DATACUBE_native.fits'
)
DEST_DIR=/scratch/pawsey1308/mauve/cubes/v3tk/
RUN_ID=$$

# Tune this (or export JOBS before running). 4 is a safe starting point.
JOBS=${JOBS:-25}

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

wait_for_overlay "$BASE_OVERLAY" || exit 1
cadc-get-cert -u RongjunHuang

# Prepare per-worker overlays so each concurrent vcp process has its own lock.
for i in $(seq 0 $((JOBS - 1))); do
	worker_overlay="$OVERLAY_DIR/cadc_overlay_worker_${RUN_ID}_${i}.img"
	cp --reflink=auto "$BASE_OVERLAY" "$worker_overlay"
	wait_for_overlay "$worker_overlay" || exit 1
done

manifest=$(mktemp)
part_dir=$(mktemp -d)

# Build source manifest once.
vls "$SRC_GLOB" | sed "s#^#$SRC_PREFIX#" > "$manifest"
printf '%s\n' "${PHANGS_NATIVE_SOURCES[@]}" >> "$manifest"

# Split work round-robin across workers.
awk -v jobs="$JOBS" -v dir="$part_dir" '{print > (dir "/part_" ((NR-1)%jobs) ".txt")}' "$manifest"

# Start one sequential worker per overlay.
for i in $(seq 0 $((JOBS - 1))); do
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
