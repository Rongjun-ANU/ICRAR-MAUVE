import os
from pathlib import Path
import re
import subprocess
import json
import pytest

SCRIPT = Path(os.environ.get('VCP_SCRIPT_UNDER_TEST', Path(__file__).with_name('vcp_scratch_v3tk_v768_to_canfar.sh')))

@pytest.mark.parametrize('galaxies', [1, 2])
@pytest.mark.parametrize('fail_file', [False, True])
@pytest.mark.parametrize('new_directory', [False, True])
def test_shared_limit_isolation_and_verification_barrier(tmp_path, galaxies, fail_file, new_directory):
    tracker = tmp_path / 'tracker.py'
    tracker.write_text('''import fcntl,json,os,pathlib,sys,time
root=pathlib.Path(os.environ['TEST_ROOT'])
def change(start):
 with (root/'guard').open('a') as guard:
  fcntl.flock(guard,fcntl.LOCK_EX)
  file=root/'state.json'
  state=json.loads(file.read_text()) if file.exists() else dict(active=[],maximum=0,done=[])
  key=sys.argv[1]
  overlay=sys.argv[2]
  if start:
   assert overlay not in [x[1] for x in state['active']], 'overlay shared concurrently'
   state['active'].append([key,overlay])
   state['maximum']=max(state['maximum'],len(state['active']))
  else:
   state['active'].remove([key,overlay]); state['done'].append(key)
  file.write_text(json.dumps(state))
change(True)
time.sleep(.25)
change(False)
''')
    names=['CONFIG', '_cont_cube.fits', '_mask.fits', '_kin_maps.fits', '_sfh_maps.fits', '_sfh_weights.fits']
    for galaxy in ['NGC4698', 'NGC4607'][:galaxies]:
        folder=tmp_path/'source'/galaxy
        folder.mkdir(parents=True)
        for suffix in names:
            (folder/(suffix if suffix=='CONFIG' else galaxy+suffix)).write_bytes(b'SIMPLE  =fixture')
    code='\n'.join(re.findall(r'(?ms)^\w+\(\) \{.*?^\}', SCRIPT.read_text()))
    harness=code+'''
set -eu
PRODUCT_SUFFIXES=(CONFIG _cont_cube.fits _mask.fits _kin_maps.fits _sfh_maps.fits _sfh_weights.fits)
SOURCE_7000="$TEST_ROOT/source"
DEST_7000=arc:test
AUTO_CHECKSUM=1
RUN_ID=test
TRANSFER_JOBS=3
worker_dir="$TEST_ROOT/workers"
mkdir -p "$worker_dir"
worker_overlays=(overlay0 overlay1 overlay2)
request_fresh_manifest() { echo "$2 $#" >> "$TEST_ROOT/requests"; printf '{}' > "$5"; }
upload_one_file() {
    [ "$(basename "$6")" = "$2" ] || return 1
    [ -f "${6%/*}/checksum_manifest.json" ] || return 1
    if [ "$NEW_DIRECTORY" = 1 ] && [ "$(basename "$3")" != CONFIG ]; then
        [ -f "$TEST_ROOT/directory_$2" ] || return 1
    fi
    printf '%s\n' "$6" >> "$TEST_ROOT/stages"
    python3 "$TEST_ROOT/tracker.py" "$2/$(basename "$3")" "${CADC_SLOT_OVERLAY:-$5}" || return 1
    if [ "$(basename "$3")" = CONFIG ]; then touch "$TEST_ROOT/directory_$2"; fi
    if [ "$FAIL_FILE" = 1 ] && [[ "$3" == *_mask.fits ]]; then return 1; fi
    return 2
}
verify_uploaded_file() {
    python3 -c 'import json,os,sys; s=json.load(open(os.environ["TEST_ROOT"]+"/state.json")); assert len([x for x in s["done"] if x.startswith(sys.argv[1]+"/")])==6' "$2"
}
upload_galaxy 7000 NGC4698 overlay0 &
p1=$!
if [ "$GALAXIES" = 2 ]; then upload_galaxy 7000 NGC4607 overlay1 & p2=$!; fi
status=0
wait "$p1" || status=1
if [ "$GALAXIES" = 2 ]; then wait "$p2" || status=1; fi
exit "$status"
'''
    env=dict(os.environ, TEST_ROOT=str(tmp_path), GALAXIES=str(galaxies), FAIL_FILE=str(int(fail_file)), NEW_DIRECTORY=str(int(new_directory)))
    result=subprocess.run(['/bin/bash','-c',harness],env=env,capture_output=True,text=True,timeout=25)
    assert result.returncode == int(fail_file),result.stdout+result.stderr
    state=json.loads((tmp_path/'state.json').read_text())
    assert state['maximum']==3, state
    assert len(state['done'])==6*galaxies
    stages=(tmp_path/'stages').read_text().splitlines()
    assert len(set(stages))==6*galaxies, 'files shared staging directories'
    assert not list((tmp_path/'workers').glob('slot_*.lock'))
    requests=(tmp_path/'requests').read_text().splitlines()
    assert sorted(requests)==sorted([g+' 11' for g in ['NGC4698','NGC4607'][:galaxies]]+[g+(' 10' if fail_file else ' 11') for g in ['NGC4698','NGC4607'][:galaxies]])

def test_interrupt_releases_slot_and_stops_descendants(tmp_path):
    import signal
    import time
    code='\n'.join(re.findall(r'(?ms)^\w+\(\) \{.*?^\}',SCRIPT.read_text()))
    harness=code+'''
worker_dir="$TEST_ROOT"
worker_overlays=(overlay0)
TRANSFER_JOBS=1
trap stop_workers TERM INT
task() { echo ready > "$TEST_ROOT/ready"; sleep 30; echo leaked > "$TEST_ROOT/leaked"; }
with_cadc_slot task &
wait
'''
    env=dict(os.environ,TEST_ROOT=str(tmp_path))
    process=subprocess.Popen(['/bin/bash','-c',harness],env=env,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,start_new_session=True)
    try:
        deadline=time.monotonic()+5
        while not (tmp_path/'ready').exists() and time.monotonic()<deadline:
            time.sleep(.05)
        assert (tmp_path/'ready').exists()
        process.send_signal(signal.SIGTERM)
        stdout,stderr=process.communicate(timeout=5)
        assert process.returncode==130,stdout+stderr
        assert not list(tmp_path.glob('slot_*.lock'))
        assert not (tmp_path/'leaked').exists()
    finally:
        try: os.killpg(process.pid,signal.SIGKILL)
        except ProcessLookupError: pass

@pytest.mark.parametrize('known_wrapper,mode,shared', [(True,'auto',True),(False,'auto',False),(True,'0',False)])
def test_overlay_preparation_shares_read_only_or_copies(tmp_path,known_wrapper,mode,shared):
    code='\n'.join(re.findall(r'(?ms)^\w+\(\) \{.*?^\}',SCRIPT.read_text()))
    wrapper=tmp_path/'vcp'
    wrapper.write_text('#!/bin/bash\n'+('/cadcenv/bin/vcp "$@"\n' if known_wrapper else 'echo custom\n'))
    (tmp_path/'base.img').write_bytes(b'base')
    harness=code+'''
set -eu
TRANSFER_JOBS=3
RUN_ID=test
BASE_OVERLAY="$TEST_ROOT/base.img"
OVERLAY_DIR="$TEST_ROOT/overlays"
VCP_CMD="$TEST_ROOT/vcp"
CHECKSUM_PYTHON=python3
worker_overlays=()
mkdir -p "$OVERLAY_DIR"
wait_for_overlay() { :; }
cp() { command cp "$2" "$3"; }
prepare_worker_overlays
printf '%s\n' "${worker_overlays[@]}" > "$TEST_ROOT/result"
'''
    result=subprocess.run(['/bin/bash','-c',harness],env=dict(os.environ,TEST_ROOT=str(tmp_path),CADC_READ_ONLY=mode),capture_output=True,text=True)
    assert result.returncode==0,result.stdout+result.stderr
    overlays=(tmp_path/'result').read_text().splitlines()
    assert len(overlays)==3
    if shared:
        assert overlays==[str(tmp_path/'base.img')+':ro']*3
        assert not list((tmp_path/'overlays').iterdir())
    else:
        assert len(set(overlays))==3
        assert all(Path(x).read_bytes()==b'base' for x in overlays)

def test_api_read_only_suffix_is_not_duplicated(tmp_path):
    code='\n'.join(re.findall(r'(?ms)^\w+\(\) \{.*?^\}',SCRIPT.read_text()))
    wrapper=tmp_path/'python'
    wrapper.write_text('#!/bin/bash\ncat >/dev/null\nprintf "%s" "$CADC_OVERLAY"\n')
    wrapper.chmod(0o755)
    result=subprocess.run(['/bin/bash','-c',code+'\nCANFAR_API=https://example.invalid\nchecksum_job_api base.img:ro status job'],env=dict(os.environ,CADC_PYTHON_CMD=str(wrapper),worker_dir=str(tmp_path)),capture_output=True,text=True)
    assert result.returncode==0,result.stderr
    assert result.stdout=='base.img:ro'
