import os
from pathlib import Path
import re
import subprocess
import pytest

SCRIPT=Path(os.environ.get('VCP_SCRIPT_UNDER_TEST',Path(__file__).with_name('vcp_scratch_v3tk_v768_to_canfar.sh')))

def test_events_api_uses_canfar_events_view(tmp_path):
    (tmp_path/'requests.py').write_text('''
class Response:
    text='Warning Failed ErrImagePull'
    def raise_for_status(self): pass
class Session:
    def get(self,url,params,timeout):
        assert params == {'view':'events'}
        assert url == 'https://example.test/skaha/v1/session/testjob'
        return Response()
''')
    wrapper=tmp_path/'python-wrapper'
    wrapper.write_text('#!/bin/bash\nexec "'+__import__('sys').executable+'" - "$@"\n')
    wrapper.chmod(0o700)
    env=dict(os.environ,PYTHONPATH=str(tmp_path),CADC_PYTHON_CMD=str(wrapper),worker_dir=str(tmp_path))
    code='\n'.join(re.findall(r'(?ms)^\w+\(\) \{.*?^\}',SCRIPT.read_text()))
    result=subprocess.run(['/bin/bash','-c',code+'\nCANFAR_API=https://example.test/skaha/v1\nchecksum_job_api overlay events testjob'],env=env,capture_output=True,text=True)
    assert result.returncode==0,result.stderr
    assert 'ErrImagePull' in result.stdout

@pytest.mark.parametrize('status,event,stage', [('Pending','Back-off pulling image "test:image"; no route to host','job-image-pull'), ('Running','Scheduled successfully','response-timeout'), ('Pending','Waiting for resources','response-timeout')])
def test_job_failure_reports_diagnostics_before_cleanup(tmp_path,status,event,stage):
    result=run_functions(tmp_path,'''\nworker_dir="$TEST_ROOT"
UPLOADER_PATH=unused
CANFAR_CHECKSUM_IMAGE=test:image
CHECKSUM_WAIT_SECONDS=2
CHECKSUM_POLL_SECONDS=1
checksum_manifest() { echo fixture_uuid; }
with_cadc_slot() { "$@"; }
slot_vcp() { case "$2" in *.checksum_response_*) echo "Response not available" >&2; return 1;; esac; return 0; }
slot_job_api() {
    echo "$2" >> "$TEST_ROOT/actions"
    case "$2" in
        submit) echo testjob;;
        status) echo "'''+status+'''";;
        events) echo '''+"'"+event+"'"+''';;
        logs) echo CHECKSUM_JOB_LOG;;
    esac
}
request_fresh_manifest 7000 NGC4698 arc:projects/test overlay "$TEST_ROOT/manifest" source.fits
status=$?
echo "STAGE=$request_error_stage DETAIL=$request_error_detail"
exit "$status"
''')
    assert result.returncode==1,result.stdout+result.stderr
    assert 'STAGE='+stage in result.stdout
    if 'no route to host' in event:
        assert 'no route to host' in result.stdout
    assert 'testjob' in result.stdout
    assert status in result.stdout
    assert event in result.stderr
    assert 'CHECKSUM_JOB_LOG' in result.stderr
    actions=(tmp_path/'actions').read_text().splitlines()
    assert actions.index('events')<actions.index('delete')
    assert actions.index('logs')<actions.index('delete')

def run_functions(tmp_path,body):
    code='\n'.join(re.findall(r'(?ms)^\w+\(\) \{.*?^\}',SCRIPT.read_text()))
    env=dict(os.environ,ISSUE_DIR=str(tmp_path/'issues'),TEST_ROOT=str(tmp_path))
    return subprocess.run(['/bin/bash','-c',code+'\n'+body],env=env,capture_output=True,text=True,timeout=15)

def test_parallel_errors_are_complete_and_deduplicated(tmp_path):
    result=run_functions(tmp_path,'''
for i in {1..20}; do
    record_issue ERROR 7000 NGC4698 NGC4698_cont_cube.fits upload "Attempts exhausted" &
    record_issue ERROR normal NGC4351 CONFIG quarantine "Remote move failed" &
done
wait
print_error_summary
''')
    assert result.returncode==0,result.stderr
    assert 'FINAL ERROR SUMMARY' in result.stderr
    assert result.stderr.count('run=7000 | galaxy=NGC4698 | file=NGC4698_cont_cube.fits | stage=upload')==1
    assert result.stderr.count('run=normal | galaxy=NGC4351 | file=CONFIG | stage=quarantine')==1

@pytest.mark.parametrize('phase', ['initial', 'post-upload'])
def test_checksum_failure_identifies_each_affected_file(tmp_path,phase):
    galaxy=tmp_path/'source/NGC4698'
    galaxy.mkdir(parents=True)
    for name in ('CONFIG','NGC4698_cont_cube.fits'):
        (galaxy/name).write_bytes(b'SIMPLE  =fixture')
    body='''
PRODUCT_SUFFIXES=(CONFIG _cont_cube.fits)
SOURCE_7000="$TEST_ROOT/source"
DEST_7000=arc:test
AUTO_CHECKSUM=1
RUN_ID=test
request_fresh_manifest() {
    if [ "$PHASE" = initial ] || [ -f "$TEST_ROOT/requested" ]; then
        request_error_stage=job-submission
        request_error_detail="API refused checksum job"
        return 1
    fi
    touch "$TEST_ROOT/requested"
    printf '{}' > "$5"
}
upload_one_file() { return 2; }
upload_galaxy 7000 NGC4698 unused
status=$?
print_error_summary
exit "$status"
'''
    result=run_functions(tmp_path,'PHASE='+phase+'\n'+body)
    assert result.returncode==1,result.stdout+result.stderr
    for name in ('CONFIG','NGC4698_cont_cube.fits'):
        assert 'run=7000 | galaxy=NGC4698 | file='+name+' | stage='+phase+'-checksum/job-submission' in result.stderr
    assert 'API refused checksum job' in result.stderr

def test_invalid_fits_source_is_in_summary(tmp_path):
    (tmp_path/'bad.fits').write_bytes(b'zero header')
    result=run_functions(tmp_path,'''
upload_one_file 7000 NGC4698 "$TEST_ROOT/bad.fits" arc:test overlay "$TEST_ROOT/stage"
status=$?
print_error_summary
exit "$status"
''')
    assert result.returncode==1
    assert 'file=bad.fits | stage=source-header' in result.stderr

def test_pending_is_separate_and_empty_summary_is_silent(tmp_path):
    result=run_functions(tmp_path,'print_error_summary')
    assert result.returncode==0 and result.stderr==''
    result=run_functions(tmp_path,'record_issue PENDING 7000 NGC4698 CONFIG verification "Fresh manual manifest required"\nprint_error_summary')
    assert result.returncode==0
    assert 'AWAITING VERIFICATION SUMMARY' in result.stderr
    assert 'FINAL ERROR SUMMARY' not in result.stderr

def test_exit_summary_survives_cleanup_and_preserves_status(tmp_path):
    result=run_functions(tmp_path,'''
worker_dir="$TEST_ROOT"
cleanup() { rm -rf "$ISSUE_DIR"; echo CLEANED >&2; }
trap finish_reporting EXIT
record_issue ERROR 7000 NGC4698 CONFIG quarantine "Move failed"
exit 1
''')
    assert result.returncode==1
    assert result.stderr.index('CLEANED')<result.stderr.index('FINAL ERROR SUMMARY')
    assert 'file=CONFIG | stage=quarantine' in result.stderr
    assert not (tmp_path/'issues').exists()

def test_cli_argument_error_has_final_summary_and_keeps_exit_code():
    result=subprocess.run(['/bin/bash',str(SCRIPT),'--dry-run','7000','NGC4698'],env=dict(os.environ,JOBS='0'),capture_output=True,text=True)
    assert result.returncode==2
    assert 'FINAL ERROR SUMMARY' in result.stderr
    assert 'run=7000 | galaxy=NGC4698 | file=(not identified) | stage=argument-validation' in result.stderr

def test_missing_required_source_is_in_summary(tmp_path):
    (tmp_path/'source/NGC4698').mkdir(parents=True)
    result=run_functions(tmp_path,'''
PRODUCT_SUFFIXES=(_cont_cube.fits)
SOURCE_7000="$TEST_ROOT/source"
DEST_7000=arc:test
AUTO_CHECKSUM=1
upload_galaxy 7000 NGC4698 overlay
status=$?
print_error_summary
exit "$status"
''')
    assert result.returncode==1
    assert 'file=NGC4698_cont_cube.fits | stage=source-discovery' in result.stderr
