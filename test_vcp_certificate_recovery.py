"""Offline certificate recovery checks; never contact CADC or read credentials."""
import os
from pathlib import Path
import re
import subprocess
import sys

SCRIPT = Path(__file__).with_name('vcp_scratch_v3tk_v768_to_canfar.sh')


def run_case(tmp_path, body):
    functions = '\n'.join(re.findall(r'(?ms)^\w+\(\) \{.*?^\}', SCRIPT.read_text()))
    harness = functions + r'''
set -euo pipefail
worker_dir="$TEST_ROOT"
CADC_USER=fixture_user
cadc-get-cert() { echo "$*" >> "$worker_dir/renewals"; return "${RENEW_RESULT:-0}"; }
validate_cadc_certificate() { return 0; }
# Portable offline flock replacement; the production script requires Linux flock.
flock() {
    python3 -c 'import fcntl,sys; fcntl.flock(int(sys.argv[1]), fcntl.LOCK_EX)' "$1"
}
operation() {
    echo call >> "$worker_dir/calls"
    if [ ! -f "$worker_dir/cert_refresh.ok" ] || [ "${PERSIST:-0}" = 1 ]; then
        echo "${ERROR_TEXT:-HTTP 500: failed to check membership with group service}" >&2
        return 7
    fi
    echo recovered
}
''' + body
    return subprocess.run(['/bin/bash', '-c', harness],
                          env=dict(os.environ, TEST_ROOT=str(tmp_path)),
                          capture_output=True, text=True)


def test_refresh_valid_but_rejected_certificate_and_retry(tmp_path):
    result = run_case(tmp_path, 'cadc_with_recovery retry operation\n')
    assert result.returncode == 0, result.stderr
    assert result.stdout == 'recovered\n'
    assert (tmp_path/'renewals').read_text() == '-u fixture_user\n'
    assert len((tmp_path/'calls').read_text().splitlines()) == 2


def test_concurrent_failures_share_one_refresh(tmp_path):
    result = run_case(tmp_path, '''
cadc_with_recovery retry operation &
one=$!
cadc_with_recovery retry operation &
two=$!
wait "$one"
wait "$two"
''')
    assert result.returncode == 0, result.stderr
    assert len((tmp_path/'renewals').read_text().splitlines()) == 1


def test_persistent_service_error_stops_after_one_retry(tmp_path):
    result = run_case(tmp_path, 'PERSIST=1; cadc_with_recovery retry operation\n')
    assert result.returncode == 7, result.stderr
    assert len((tmp_path/'calls').read_text().splitlines()) == 2
    assert len((tmp_path/'renewals').read_text().splitlines()) == 1


def test_ambiguous_mutation_refreshes_without_replay(tmp_path):
    result = run_case(tmp_path, 'cadc_with_recovery no-retry operation\n')
    assert result.returncode == 7, result.stderr
    assert len((tmp_path/'calls').read_text().splitlines()) == 1
    assert (tmp_path/'renewals').exists()


def test_unrelated_failure_never_prompts(tmp_path):
    result = run_case(tmp_path, 'ERROR_TEXT="HTTP 500: database unavailable"; cadc_with_recovery retry operation\n')
    assert result.returncode == 7, result.stderr
    assert not (tmp_path/'renewals').exists()


def test_failed_refresh_is_not_repeated(tmp_path):
    result = run_case(tmp_path, '''
RENEW_RESULT=1
cadc_with_recovery retry operation || true
cadc_with_recovery retry operation
''')
    assert result.returncode == 7, result.stderr
    assert len((tmp_path/'renewals').read_text().splitlines()) == 1
    assert len((tmp_path/'calls').read_text().splitlines()) == 2


def test_tls_error_refreshes(tmp_path):
    result = run_case(tmp_path, 'ERROR_TEXT="SSLHandshakeException: (decrypt_error) Received fatal alert: decrypt_error"; cadc_with_recovery retry operation\n')
    assert result.returncode == 0, result.stderr
    assert (tmp_path/'renewals').exists()


def test_connection_reset_refreshes(tmp_path):
    result = run_case(tmp_path, 'ERROR_TEXT="Connection reset by peer"; cadc_with_recovery retry operation\n')
    assert result.returncode == 0, result.stderr
    assert (tmp_path/'renewals').exists()


def test_api_probe_recovers_before_single_post(tmp_path):
    (tmp_path/'requests.py').write_text('''
import os
from pathlib import Path
root = Path(os.environ['TEST_ROOT'])
class Response:
    def __init__(self, status, text):
        self.status_code, self.text = status, text
    def raise_for_status(self):
        if self.status_code >= 400:
            error = RuntimeError('500 Server Error')
            error.response = self
            raise error
class Session:
    def get(self, url, timeout):
        with (root/'api_calls').open('a') as stream: stream.write('GET\\n')
        if not (root/'cert_refresh.ok').exists():
            return Response(500, 'failed to check membership with group service')
        return Response(200, '[]')
    def post(self, url, data, timeout):
        with (root/'api_calls').open('a') as stream: stream.write('POST\\n')
        return Response(200, 'fixture_job')
''')
    wrapper = tmp_path/'python-wrapper'
    wrapper.write_text('#!/bin/bash\nexec "'+sys.executable+'" - "$@"\n')
    wrapper.chmod(0o700)
    result = run_case(tmp_path, '''
export PYTHONPATH="$TEST_ROOT"
CADC_PYTHON_CMD="$TEST_ROOT/python-wrapper"
CANFAR_API=https://example.invalid/skaha/v1
BASE_OVERLAY=fixture
checksum_job_api fixture submit image worker request name
''')
    assert result.returncode == 0, result.stderr
    assert result.stdout == 'fixture_job\n'
    assert (tmp_path/'api_calls').read_text().splitlines() == ['GET', 'GET', 'POST']


def test_failed_submit_never_reposts(tmp_path):
    result = run_case(tmp_path, '''
checksum_job_api_raw() {
    echo "$2" >> "$worker_dir/api_calls"
    [ "$2" != submit ] || operation
}
checksum_job_api fixture submit image worker request name
''')
    assert result.returncode == 7, result.stderr
    assert (tmp_path/'api_calls').read_text().splitlines() == ['probe', 'submit']
    assert len((tmp_path/'renewals').read_text().splitlines()) == 1


def test_invalid_refreshed_certificate_prevents_retry(tmp_path):
    result = run_case(tmp_path, '''
validate_cadc_certificate() { return 1; }
cadc_with_recovery retry operation
''')
    assert result.returncode == 7, result.stderr
    assert len((tmp_path/'calls').read_text().splitlines()) == 1
    assert not (tmp_path/'cert_refresh.ok').exists()
