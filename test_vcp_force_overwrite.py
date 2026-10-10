import os
from pathlib import Path
import re
import subprocess
import sys
import pytest

SCRIPT=Path(os.environ.get('VCP_SCRIPT_UNDER_TEST',Path(__file__).with_name('vcp_scratch_v3tk_v768_to_canfar.sh')))

def test_usage_is_at_beginning_and_cli_accepts_flags(tmp_path):
    result=subprocess.run(['/bin/bash',str(SCRIPT),'--dry-run','7000','NGC4216','NGC4689','--force-overwrite','--no-checksum'],capture_output=True,text=True)
    assert result.returncode==0,result.stderr
    assert 'FORCE OVERWRITE' in result.stdout and 'UNVERIFIED' in result.stdout
    assert SCRIPT.read_text().index('usage() {')<SCRIPT.read_text().index('checksum_manifest() {')

@pytest.mark.parametrize('flags',[['--force-overwrite'],['--no-checksum'],['--no-checksum','--force-overwrite','--auto-checksum'],['--manual-checksum','--force-overwrite','--no-checksum']])
def test_unsafe_or_conflicting_flags_rejected(flags):
    result=subprocess.run(['/bin/bash',str(SCRIPT),'--dry-run']+flags+['7000','NGC4698'],capture_output=True,text=True)
    assert result.returncode==2
    assert 'requires' in result.stderr or 'cannot combine' in result.stderr

@pytest.mark.parametrize('existing,fail_first,metadata_failure,move_failure',[(True,False,False,False),(False,False,False,False),(True,True,False,False),(True,False,True,False),(True,False,False,True)])
def test_force_upload_preserves_old_bytes_and_never_uses_checksums(tmp_path,existing,fail_first,metadata_failure,move_failure):
    source=tmp_path/'source/NGC4698'
    source.mkdir(parents=True)
    original=b'SIMPLE  =old-data'
    new=b'SIMPLE  =new-data'
    (source/'NGC4698_cont_cube.fits').write_bytes(new)
    target=tmp_path/'remote/NGC4698/NGC4698_cont_cube.fits'
    if existing:
        target.parent.mkdir(parents=True)
        target.write_bytes(original)
    tool=tmp_path/'tool.py'
    tool.write_text('''#!'''+sys.executable+'''
import os,pathlib,sys,shutil
root=pathlib.Path(os.environ['TEST_ROOT'])
args=[a for a in sys.argv[1:] if a!='-v']
mode=pathlib.Path(sys.argv[0]).name
with (root/'actions').open('a') as f:f.write(mode+' '+args[0]+'\\n')
def remote(uri):return root/'remote'/uri.removeprefix('arc:test/')
if mode=='vmv':
 if os.environ['MOVE_FAILURE']=='1':sys.exit(1)
 remote(args[0]).rename(remote(args[1]))
else:
 assert not args[0].startswith('arc:'), 'downloads forbidden'
 src=next(pathlib.Path(args[0]).iterdir());dst=remote(args[1])/pathlib.Path(args[0]).name/src.name
 dst.parent.mkdir(parents=True,exist_ok=True)
 assert not dst.exists(), 'vcp could skip existing destination'
 if os.environ['FAIL_FIRST']=='1' and not (root/'first_failed').exists():
  dst.write_bytes(b'SIMPLE  =partial');(root/'first_failed').touch();sys.exit(1)
 shutil.copyfile(src,dst)
''')
    tool.chmod(0o700)
    for name in ('vcp','vmv'):(tmp_path/name).symlink_to(tool)
    functions='\n'.join(re.findall(r'(?ms)^\w+\(\) \{.*?^\}',SCRIPT.read_text()))
    harness=functions+'''
set -eu
PRODUCT_SUFFIXES=(_cont_cube.fits)
SOURCE_7000="$TEST_ROOT/source"
DEST_7000=arc:test
AUTO_CHECKSUM=0
NO_CHECKSUM=1
FORCE_OVERWRITE=1
RUN_ID=fixture
FILE_RETRIES=2
RETRY_BASE_SLEEP=1
RETRY_MAX_SLEEP=1
CHECKSUM_RECEIPT_DIR="$TEST_ROOT/receipts"
stat() { python3 -c 'import os,sys; print(os.stat(sys.argv[-1]).st_size)' "$@"; }
sha256sum() { shasum -a 256 "$@"; }
sleep() { :; }
remote_file_exists() {
    [ "$METADATA_FAILURE" = 0 ] || return 1
    if [ -f "$TEST_ROOT/remote/${2#arc:test/}" ]; then return 0; else return 3; fi
}
checksum_manifest() { echo FORBIDDEN_MANIFEST >&2; return 1; }
request_fresh_manifest() { echo FORBIDDEN_JOB >&2; return 1; }
verify_uploaded_file() { echo FORBIDDEN_VERIFY >&2; return 1; }
if upload_galaxy 7000 NGC4698 overlay; then status=0; else status=$?; fi
print_error_summary
exit "$status"
'''
    env=dict(os.environ,TEST_ROOT=str(tmp_path),worker_dir=str(tmp_path),ISSUE_DIR=str(tmp_path/'issues'),VCP_CMD=str(tmp_path/'vcp'),VMV_CMD=str(tmp_path/'vmv'),FAIL_FIRST=str(int(fail_first)),METADATA_FAILURE=str(int(metadata_failure)),MOVE_FAILURE=str(int(move_failure)))
    result=subprocess.run(['/bin/bash','-c',harness],env=env,capture_output=True,text=True,timeout=15)
    assert result.returncode==int(metadata_failure or move_failure),result.stdout+result.stderr
    assert 'FORBIDDEN' not in result.stdout+result.stderr
    if metadata_failure or move_failure:
        assert target.read_bytes()==original
        assert 'vcp ' not in (tmp_path/'actions').read_text() if (tmp_path/'actions').exists() else True
        assert 'FINAL ERROR SUMMARY' in result.stderr
        assert 'stage='+('destination-lookup' if metadata_failure else 'backup') in result.stderr
    else:
        assert target.read_bytes()==new
        backups=list(target.parent.glob('*.overwrite_backup_*'))
        assert len(backups)==int(existing)+int(fail_first)
        if existing:assert any(p.read_bytes()==original for p in backups)
        if fail_first:assert any(p.read_bytes()==b'SIMPLE  =partial' for p in backups)
        assert 'UNVERIFIED' in result.stdout
        assert 'files verified by CANFAR manifest' not in result.stdout

def test_forced_receipt_blocks_old_manual_manifest(tmp_path):
    source=tmp_path/'source.fits';source.write_bytes(b'SIMPLE  =source')
    stage=tmp_path/'stage/NGC4698';stage.mkdir(parents=True)
    receipts=tmp_path/'receipts';receipts.mkdir()
    import hashlib
    key=hashlib.sha256(b'arc:test/NGC4698/source.fits').hexdigest()
    digest=hashlib.sha256(source.read_bytes()).hexdigest()
    (receipts/(key+'.receipt')).write_text('forced_fixture '+digest+'\n')
    functions='\n'.join(re.findall(r'(?ms)^\w+\(\) \{.*?^\}',SCRIPT.read_text()))
    harness=functions+'''
AUTO_CHECKSUM=0
NO_CHECKSUM=0
stat() { echo 15; }
sha256sum() { shasum -a 256 "$@"; }
checksum_manifest() { echo "15 '''+digest+''' aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"; }
upload_one_file 7000 NGC4698 "$TEST_ROOT/source.fits" arc:test overlay "$TEST_ROOT/stage/NGC4698"
'''
    result=subprocess.run(['/bin/bash','-c',harness],env=dict(os.environ,TEST_ROOT=str(tmp_path),CHECKSUM_RECEIPT_DIR=str(receipts)),capture_output=True,text=True)
    assert result.returncode==2,result.stdout+result.stderr
    assert 'automatic' in result.stdout
    assert (receipts/(key+'.receipt')).exists()
