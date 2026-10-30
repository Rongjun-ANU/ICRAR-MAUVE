#!/usr/bin/env bash
# Usage: bash SFH16_send.sh [USER[@HOST]] [GALID ...]
set -euo pipefail
cd "$(dirname "$0")"
. ./SFH16_galaxies.sh
remote_login=rhuang@setonix.pawsey.org.au
if [[ $# -gt 0 ]] && ! looks_like_galid "$1"; then
  remote_login="$1"
  shift
  [[ "$remote_login" == *@* ]] || remote_login="${remote_login}@setonix.pawsey.org.au"
fi
[[ "$remote_login" =~ ^[A-Za-z0-9._-]+@[A-Za-z0-9.-]+$ ]] || {
  echo "ERROR: expected USER or USER@HOST" >&2; exit 2;
}
select_sfh16_galids "$@"
run_dir=/software/projects/pawsey1308/ngist_supplementary_public/ngistTutorial
yaml_files=()
slurm_files=()
for galid in "${SELECTED_GALIDS[@]}"; do
  yaml_files+=("${galid}_MAUVE_MasterConfig_v7.6.8_7000_SFH16_setonix.yaml")
  slurm_files+=("${galid}_v3tk_v7.6.8_7000_SFH16_setonix.slurm")
done
control_files=(27_galaxies.sh SFH16_ready.sh SFH16_galaxies.sh SFH16_workflow.py)
for file in "${yaml_files[@]}" "${slurm_files[@]}" "${control_files[@]}"; do
  [[ -s "$file" ]] || { echo "ERROR: missing $file; run SFH16_creation.sh first" >&2; exit 1; }
done
echo "Sending ${#yaml_files[@]} SFH16 configs and jobs to ${remote_login}"
ssh "$remote_login" "mkdir -p '$run_dir/configFiles'"
COPYFILE_DISABLE=1 tar -cf - "${yaml_files[@]}" | ssh "$remote_login" \
  "tar -xf - -C '$run_dir/configFiles'"
COPYFILE_DISABLE=1 tar -cf - "${slurm_files[@]}" "${control_files[@]}" | ssh "$remote_login" \
  "tar -xf - -C '$run_dir' && chmod +x '$run_dir/SFH16_galaxies.sh' '$run_dir/SFH16_ready.sh'"
echo "Transfer complete. On Setonix, run: bash SFH16_ready.sh ${SELECTED_GALIDS[*]}"
echo "Then submit with: ./SFH16_galaxies.sh ${SELECTED_GALIDS[*]}"
