#!/usr/bin/env bash
# 2026-10-30: Stage missing PHANGS-native cubes as ordinary Mass.sh does.
# 2026-10-30: Reuse neighbouring 7000 binning/KIN/gas with new SFH16 inputs.
# 2026-10-29: Mass.py followed by QC for 40 galaxies in 39 SFH16 cubes.
# See QC_ngist_v3tk_v768_7000_SFH16.md for inputs, outputs and examples.
set -euo pipefail
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
version=v3tk_v7.6.8_7000
product_root="/arc/projects/mauve/products/$version"
[[ -d "$product_root" ]] || product_root="$PWD/$version"
cube_root="${CUBE_ROOT:-/arc/projects/mauve/cubes/v3tk}"
output_dir="$PWD/QC"
dry_run=false
requested=()
all_galids=(
  IC3392 NGC4064 NGC4189 NGC4192 NGC4216 NGC4222 NGC4254 NGC4293
  NGC4294 NGC4298 NGC4302 NGC4321 NGC4330 NGC4351 NGC4380 NGC4383
  NGC4388 NGC4394 NGC4396 NGC4405 NGC4402 NGC4419 NGC4424 NGC4450
  NGC4457 NGC4501 NGC4522 NGC4535 NGC4548 NGC4567_8 NGC4569 NGC4579
  NGC4580 NGC4606 NGC4607 NGC4654 NGC4689 NGC4694 NGC4698
)
while [[ $# -gt 0 ]]; do
  case "$1" in
    --root|--cube-root|--output-dir)
      [[ $# -ge 2 && -n "$2" ]] || { echo "ERROR: $1 needs a value" >&2; exit 2; }
      case "$1" in
        --root) product_root="$2" ;;
        --cube-root) cube_root="$2" ;;
        --output-dir) output_dir="$2" ;;
      esac
      shift 2 ;;
    --dry-run) dry_run=true; shift ;;
    -h|--help)
      echo "Usage: bash ${0##*/} [--root DIR] [--cube-root DIR] [--output-dir DIR] [--dry-run] [GALID ...]"
      echo "Defaults to all 39 cube IDs (40 galaxies); missing inputs are skipped."
      echo "Explicit missing targets fail. GALID may include _SFH16. Runs sequentially."
      echo "Reads SFH16 CONFIG/LOGFILE/maps/weights; reuses binning/KIN/gas from neighbouring GALID."
      exit 0 ;;
    --*) echo "ERROR: unknown option $1" >&2; exit 2 ;;
    *) requested+=("${1%_SFH16}"); shift ;;
  esac
done
explicit=false
if [[ ${#requested[@]} -gt 0 ]]; then
  explicit=true
else
  requested=("${all_galids[@]}")
fi
# Validate the entire selection before writing any products.
selected=()
for galid in "${requested[@]}"; do
  known=false
  for candidate in "${all_galids[@]}"; do
    [[ "$galid" != "$candidate" ]] || known=true
  done
  [[ "$known" == true ]] || { echo "ERROR: unknown cube ID $galid" >&2; exit 2; }
  duplicate=false
  for candidate in "${selected[@]+${selected[@]}}"; do
    [[ "$galid" != "$candidate" ]] || duplicate=true
  done
  [[ "$duplicate" == true ]] || selected+=("$galid")
done
echo "SFH16 QC: ${#selected[@]} cube IDs selected; full sample = 39 cube IDs / 40 galaxies."
if [[ "$dry_run" == true ]]; then
  for galid in "${selected[@]}"; do
    echo "PLAN ${galid}_SFH16: CONFIG, LOGFILE, SFH maps/weights from $product_root/${galid}_SFH16"
    echo "REUSE $galid: spatial binning, KIN, gas BIN/SPAXEL from $product_root/$galid"
    echo "PLAN Mass.py -g $galid via private input view -> $product_root/${galid}_SFH16/${galid}_SFH16_spatial_binning_maps_further.fits"
    echo "PLAN QC: existing QC plotting code via private input view; output $output_dir/${galid}_SFH16_${version}_QC.pdf"
  done
  exit 0
fi
[[ -d "$product_root" ]] || { echo "ERROR: missing product root $product_root" >&2; exit 1; }
product_root="$(cd "$product_root" && pwd -P)"
python_bin="${PYTHON_BIN:-}"
if [[ -z "$python_bin" ]]; then
  if [[ -n "${CONDA_PREFIX:-}" && -x "$CONDA_PREFIX/bin/python" ]]; then
    python_bin="$CONDA_PREFIX/bin/python"
  elif [[ -x /opt/miniconda3/envs/ICRAR/bin/python ]]; then
    python_bin=/opt/miniconda3/envs/ICRAR/bin/python
  else
    python_bin=python
  fi
fi
command -v "$python_bin" >/dev/null || { echo "ERROR: Python not found: $python_bin" >&2; exit 1; }
mkdir -p "$output_dir/logs_SFH16"
output_dir="$(cd "$output_dir" && pwd -P)"
work_view="$(mktemp -d "${TMPDIR:-/tmp}/qc_SFH16.XXXXXX")"
work_view="$(cd "$work_view" && pwd -P)"
phangs_stage_dir=""
release_phangs_cube() {
  if [[ -n "$phangs_stage_dir" ]]; then
    rm -rf -- "$phangs_stage_dir"
    phangs_stage_dir=""
  fi
}
cleanup() {
  release_phangs_cube
  rm -rf -- "$work_view"
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
ln -s "$output_dir" "$work_view/QC"
# Match established FITS spellings, without selecting old SFH16 duplicates.
find_fits() {
  local folder=$1 prefix=$2 suffix extension path
  shift 2
  for suffix in "$@"; do
    for extension in fits fits.gz; do
      path="$folder/${prefix}_${suffix}.$extension"
      if [[ -f "$path" ]]; then printf '%s\n' "$path"; return 0; fi
    done
  done
  return 1
}
mass_args=()
[[ -z "${TEMPLATE_LIBRARY_DIR:-}" ]] || mass_args+=(--template-library-dir "$TEMPLATE_LIBRARY_DIR")
phangs_args=()
prepare_phangs_cube() {
  phangs_args=()
  case "$galid" in NGC4254|NGC4321|NGC4535) ;; *) return 0 ;; esac
  local filename="${galid}_PHANGS_DATACUBE_native.fits"
  local folder stage_root source rc
  stage_root="${PHANGS_STAGING_ROOT:-/scratch}"
  for folder in "${PHANGS_CUBE_ROOT:-}" "$cube_root" "$stage_root"; do
    [[ -n "$folder" ]] || continue
    if [[ -s "$folder/$filename" || -s "$folder/$filename.gz" ]]; then
      echo "Using local PHANGS cube: $folder/$filename[.gz]"
      phangs_args=(--phangs-cube-root "$folder")
      return 0
    fi
  done
  command -v vcp >/dev/null 2>&1 || {
    echo "ERROR: no local PHANGS cube for $galid and vcp is unavailable." >&2
    return 1
  }
  mkdir -p "$stage_root" || return $?
  phangs_stage_dir=$(mktemp -d "${stage_root%/}/qc_SFH16_phangs.XXXXXX") || return $?
  source="${PHANGS_NATIVE_VOS_DIR:-vos:phangs/RELEASES/PHANGS-MUSE/DR1.0/DATACUBES}"
  echo "Staging PHANGS cube: ${source%/}/$filename -> $phangs_stage_dir/$filename"
  if vcp "${source%/}/$filename" "$phangs_stage_dir/$filename"; then
    [[ -s "$phangs_stage_dir/$filename" ]] || {
      echo "ERROR: vcp returned success without a nonempty cube." >&2
      return 1
    }
    phangs_args=(--phangs-cube-root "$phangs_stage_dir")
  else
    rc=$?
    echo "ERROR: PHANGS cube download failed for $galid (exit $rc)." >&2
    return "$rc"
  fi
}
passed=0; skipped=0; failed=0
for galid in "${selected[@]}"; do
  run_id="${galid}_SFH16"
  run_dir="$product_root/$run_id"
  original_dir="$product_root/$galid"
  missing=()
  for name in CONFIG LOGFILE; do
    [[ -f "$run_dir/$name" ]] || missing+=("$run_dir/$name")
  done
  sfh_map=$(find_fits "$run_dir" "$run_id" sfh_maps SFH_maps) || missing+=("$run_dir/${run_id}_sfh_maps.fits[.gz]")
  sfh_weights=$(find_fits "$run_dir" "$run_id" sfh_weights sfh-weights) || missing+=("$run_dir/${run_id}_sfh_weights.fits[.gz]")
  binning=$(find_fits "$original_dir" "$galid" spatial_binning_maps SPATIAL_BINNING_maps) || missing+=("$original_dir/${galid}_spatial_binning_maps.fits[.gz]")
  kin=$(find_fits "$original_dir" "$galid" kin_maps KIN_maps) || missing+=("$original_dir/${galid}_kin_maps.fits[.gz]")
  gas_bin=$(find_fits "$original_dir" "$galid" gas_bin_maps gas_BIN_maps GAS_BIN_maps) || missing+=("$original_dir/${galid}_gas_bin_maps.fits[.gz]")
  gas_spaxel=$(find_fits "$original_dir" "$galid" gas_spaxel_maps gas_SPAXEL_maps GAS_SPAXEL_maps) || missing+=("$original_dir/${galid}_gas_spaxel_maps.fits[.gz]")
  if [[ ${#missing[@]} -gt 0 ]]; then
    echo "SKIP $run_id: missing ${missing[*]}"
    if [[ "$explicit" == true ]]; then failed=$((failed+1)); else skipped=$((skipped+1)); fi
    continue
  fi
  # Mass resolves symlinks before matching the shared SFH16 filename keyword.
  # A hardlink keeps the renamed binning basename; across filesystems, only
  # this binning map needs a private temporary copy. All other inputs are links.
  run_view="$work_view/$galid"
  mkdir "$run_view"
  ln -s "$galid" "$work_view/$run_id"
  binning_alias="$run_view/${run_id}_spatial_binning_maps.fits"
  [[ "$binning" != *.gz ]] || binning_alias+=.gz
  if ! ln -L "$binning" "$binning_alias" 2>/dev/null; then
    cp "$binning" "$binning_alias"
  fi
  ln -s "$sfh_map" "$run_view/${run_id}_sfh_maps.fits"
  ln -s "$sfh_weights" "$run_view/${run_id}_sfh_weights.fits"
  ln -s "$kin" "$run_view/${run_id}_kin_maps.fits"
  ln -s "$gas_bin" "$run_view/${run_id}_gas_bin_maps.fits"
  ln -s "$gas_spaxel" "$run_view/${run_id}_gas_spaxel_maps.fits"
  mass_log="$output_dir/logs_SFH16/${run_id}_mass.log"
  qc_log="$output_dir/logs_SFH16/${run_id}_qc.log"
  echo "MASS $run_id -> $mass_log"
  echo "REUSE $galid: $binning; $kin; $gas_bin; $gas_spaxel"
  if {
    printf 'SFH16 source: %s\nOriginal 7000 binning: %s\nSFH maps: %s\nSFH weights: %s\n' \
      "$run_dir" "$binning" "$sfh_map" "$sfh_weights"
    prepare_phangs_cube && "$python_bin" -u "$script_dir/Mass.py" -g "$galid" \
      --root "$work_view" --product-subdir . --fallback-root "$work_view" \
      --ngist-config "$run_dir/CONFIG" \
      --cube-root "$cube_root" ${phangs_args[@]+"${phangs_args[@]}"} \
      ${mass_args[@]+"${mass_args[@]}"}
  } >"$mass_log" 2>&1; then
    release_phangs_cube
    fresh_mass="$run_view/${run_id}_spatial_binning_maps_further.fits"
    if [[ ! -s "$fresh_mass" ]]; then
      echo "FAIL $run_id: Mass exited 0 but did not produce $fresh_mass"
      failed=$((failed+1)); continue
    fi
  else
    rc=$?
    release_phangs_cube
    echo "FAIL $run_id: Mass exit $rc; see $mass_log"
    failed=$((failed+1)); continue
  fi
  mass_product="$run_dir/${run_id}_spatial_binning_maps_further.fits"
  mv -f -- "$fresh_mass" "$mass_product"
  ln -s "$mass_product" "$fresh_mass"
  echo "QC $run_id -> $qc_log"
  # Override only QC's top-level data_folder assignment in memory. The actual
  # source file, plotting code and scientific calculations remain unchanged.
  # This makes the private view work even when /arc's priority root exists.
  if (cd "$work_view" && MPLBACKEND=Agg "$python_bin" - \
      "$script_dir/QC_ngist_v3tk_v768_7000.py" "$run_id" "$work_view" <<'SFH16_QC_PYTHON'
import ast
import sys
from pathlib import Path

source, run_id, product_view = sys.argv[1:]
tree = ast.parse(Path(source).read_text(), filename=source)
assignments = [node for node in tree.body if isinstance(node, ast.Assign)
               and any(isinstance(target, ast.Name) and target.id == 'data_folder'
                       for target in node.targets)]
if len(assignments) != 1:
    raise RuntimeError('QC data_folder assignment changed; review the SFH16 launcher.')
assignments[0].value = ast.copy_location(ast.Constant(product_view), assignments[0].value)
ast.fix_missing_locations(tree)
sys.argv = [source, run_id]
sys.path.insert(0, str(Path(source).parent))
exec(compile(tree, source, 'exec'), {'__name__': '__main__', '__file__': source})
SFH16_QC_PYTHON
  ) >"$qc_log" 2>&1; then
    echo "PASS $run_id"
    passed=$((passed+1))
  else
    rc=$?
    echo "FAIL $run_id: QC exit $rc; see $qc_log"
    failed=$((failed+1))
  fi
done
echo "SFH16 7000: $passed passed, $skipped skipped, $failed failed."
[[ "$failed" -eq 0 && "$passed" -gt 0 ]]
