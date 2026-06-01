#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG_DIR="${REPO_ROOT}/config_setup"
TMP_DIRS=()

fail() {
  echo "FAIL: $*" >&2
  exit 1
}

cleanup() {
  local tmp_dir
  if [[ ${#TMP_DIRS[@]} -eq 0 ]]; then
    return
  fi
  for tmp_dir in "${TMP_DIRS[@]}"; do
    rm -rf "$tmp_dir"
  done
}
trap cleanup EXIT

new_workdir() {
  local tmp_dir
  tmp_dir="$(mktemp -d)"
  TMP_DIRS+=("$tmp_dir")

  cp "${CONFIG_DIR}/27_creation.sh" "$tmp_dir/"
  cp "${CONFIG_DIR}/27_send.sh" "$tmp_dir/"
  cp "${CONFIG_DIR}/27_setonix.sh" "$tmp_dir/"
  cp "${CONFIG_DIR}/27_status.sh" "$tmp_dir/"
  cp "${CONFIG_DIR}/v3tk_v7.6.8_setonix.slurm" "$tmp_dir/"
  if [[ -f "${CONFIG_DIR}/27_galaxies.sh" ]]; then
    cp "${CONFIG_DIR}/27_galaxies.sh" "$tmp_dir/"
  fi

  printf "%s\n" "$tmp_dir"
}

install_make_gist_stub() {
  local workdir="$1"
  cat > "${workdir}/make_gist_config_try.py" <<'STUB'
#!/usr/bin/env bash
set -euo pipefail
galid="$1"
printf "%s\n" "$galid" >> created_galaxies.txt
touch "${galid}_MAUVE_MasterConfig_v7.6.8_setonix.yaml"
STUB
  chmod +x "${workdir}/make_gist_config_try.py"
}

install_sbatch_stub() {
  local workdir="$1"
  mkdir -p "${workdir}/bin"
  cat > "${workdir}/bin/sbatch" <<'STUB'
#!/usr/bin/env bash
set -euo pipefail
printf "%s\n" "$*" >> "${SBATCH_LOG}"
STUB
  chmod +x "${workdir}/bin/sbatch"
}

install_ssh_stub() {
  local workdir="$1"
  mkdir -p "${workdir}/bin"
cat > "${workdir}/bin/ssh" <<'STUB'
#!/usr/bin/env bash
set -euo pipefail
printf "%s\n" "$*" >> "${SSH_LOG}"
case "$*" in
  *"tar -xf"*)
    count="$(find "${SSH_TAR_DIR}" -type f -name 'payload_*.tar' 2>/dev/null | wc -l | tr -d ' ')"
    cat > "${SSH_TAR_DIR}/payload_${count}.tar"
    ;;
esac
STUB
  chmod +x "${workdir}/bin/ssh"
}

assert_file() {
  [[ -f "$1" ]] || fail "expected file to exist: $1"
}

assert_no_file() {
  [[ ! -f "$1" ]] || fail "expected file not to exist: $1"
}

assert_contains() {
  local file="$1"
  local text="$2"
  grep -Fq "$text" "$file" || fail "expected ${file} to contain: ${text}"
}

assert_not_contains() {
  local file="$1"
  local text="$2"
  if grep -Fq "$text" "$file"; then
    fail "expected ${file} not to contain: ${text}"
  fi
}

test_creation_uses_requested_galaxies_only() {
  local workdir
  workdir="$(new_workdir)"
  install_make_gist_stub "$workdir"

  (cd "$workdir" && ./27_creation.sh IC3392 NGC4383 > creation.out)

  diff -u <(printf "IC3392\nNGC4383\n") "${workdir}/created_galaxies.txt"
  assert_file "${workdir}/IC3392_v3tk_v7.6.8_setonix.slurm"
  assert_file "${workdir}/NGC4383_v3tk_v7.6.8_setonix.slurm"
  assert_no_file "${workdir}/NGC4698_v3tk_v7.6.8_setonix.slurm"
  assert_contains "${workdir}/creation.out" "Created 2 YAML files and 2 slurm scripts."
}

test_setonix_submits_requested_galaxies_only() {
  local workdir
  workdir="$(new_workdir)"
  touch "${workdir}/IC3392_v3tk_v7.6.8_setonix.slurm"
  touch "${workdir}/NGC4383_v3tk_v7.6.8_setonix.slurm"
  install_sbatch_stub "$workdir"

  (
    export SBATCH_LOG="${workdir}/sbatch.log"
    export PATH="${workdir}/bin:${PATH}"
    cd "$workdir"
    ./27_setonix.sh NGC4383 IC3392 > setonix.out
  )

  diff -u <(printf "NGC4383_v3tk_v7.6.8_setonix.slurm\nIC3392_v3tk_v7.6.8_setonix.slurm\n") "${workdir}/sbatch.log"
  assert_contains "${workdir}/setonix.out" "Submitted 2 Setonix jobs."
}

test_status_reports_requested_galaxies_only() {
  local workdir
  workdir="$(new_workdir)"

  (cd "$workdir" && ./27_status.sh IC3392 NGC4383 > status.out)

  assert_contains "${workdir}/status.out" "IC3392"
  assert_contains "${workdir}/status.out" "NGC4383"
  assert_not_contains "${workdir}/status.out" "NGC4698"
}

test_send_copies_requested_galaxies_only() {
  local workdir
  workdir="$(new_workdir)"
  install_ssh_stub "$workdir"
  touch "${workdir}/NGC4383_MAUVE_MasterConfig_v7.6.8_setonix.yaml"
  touch "${workdir}/NGC4419_MAUVE_MasterConfig_v7.6.8_setonix.yaml"
  touch "${workdir}/NGC4383_v3tk_v7.6.8_setonix.slurm"
  touch "${workdir}/NGC4419_v3tk_v7.6.8_setonix.slurm"
  mkdir -p "${workdir}/ssh_payloads"

  (
    export SSH_LOG="${workdir}/ssh.log"
    export SSH_TAR_DIR="${workdir}/ssh_payloads"
    export PATH="${workdir}/bin:${PATH}"
    cd "$workdir"
    ./27_send.sh rhuang NGC4383 NGC4419 > send.out
  )

  assert_contains "${workdir}/send.out" "Sending 2 YAML files"
  assert_contains "${workdir}/send.out" "Sending 2 slurm scripts plus control scripts"
  assert_contains "${workdir}/send.out" "Remote login: rhuang@setonix.pawsey.org.au"
  tar -tf "${workdir}/ssh_payloads/payload_1.tar" > "${workdir}/control_payload.txt"
  assert_contains "${workdir}/control_payload.txt" "27_galaxies.sh"
}

test_unknown_galaxy_is_rejected() {
  local workdir
  workdir="$(new_workdir)"

  if (cd "$workdir" && ./27_setonix.sh NGC0000 > unknown.out 2>&1); then
    fail "unknown galaxy should have failed"
  fi

  assert_contains "${workdir}/unknown.out" "unknown galaxy ID: NGC0000"
}

test_creation_uses_requested_galaxies_only
test_setonix_submits_requested_galaxies_only
test_status_reports_requested_galaxies_only
test_send_copies_requested_galaxies_only
test_unknown_galaxy_is_rejected

echo "All 27-galaxy selection tests passed."
