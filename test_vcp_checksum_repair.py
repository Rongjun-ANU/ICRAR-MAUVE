"""Offline checks: CANFAR byte-hash manifests, no FITS read-back, scoped repair."""

import importlib.util
import os
from pathlib import Path
import re
import subprocess
import sys
from types import SimpleNamespace

import pytest

SCRIPT = Path(os.environ.get("VCP_SCRIPT_UNDER_TEST", Path(__file__).with_name("vcp_scratch_v3tk_v768_to_canfar.sh")))
HELPER = SCRIPT.with_name("canfar_checksum_manifest.py")
spec = importlib.util.spec_from_file_location("manifest_helper", HELPER)
manifest_helper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(manifest_helper)


@pytest.fixture
def transfer(tmp_path):
    source = tmp_path / "source" / "NGC4698" / "NGC4698_cont_cube.fits"
    source.parent.mkdir(parents=True)
    source.write_bytes(b"SIMPLE  =" + b"valid-source-content" * 30)
    remote = tmp_path / "remote"
    target = remote / "v3tk_v7.6.8_7000" / "NGC4698" / source.name
    target.parent.mkdir(parents=True)
    target.write_bytes(b"damaged-existing-file")
    tools = tmp_path / "tools"
    tools.mkdir()
    stage = tmp_path / "stage" / "NGC4698"
    stage.mkdir(parents=True)
    fake = tools / "fake.py"
    fake.write_text('''#!/usr/bin/env python3
import hashlib, os, pathlib, shutil, sys
mode = pathlib.Path(sys.argv[0]).name
args = [arg for arg in sys.argv[1:] if arg != '-v']
remote = pathlib.Path(os.environ['FAKE_REMOTE'])
events = pathlib.Path(os.environ['FAKE_EVENTS'])
def event(text):
    with events.open('a') as f: f.write(text + '\\n')
def node(uri): return remote / uri.removeprefix('arc:products/')
if mode in ('flock', 'cadc-get-cert'):
    pass
elif mode == 'cp':
    shutil.copyfile(args[-2], args[-1])
elif mode == 'stat':
    print(pathlib.Path(args[-1]).stat().st_size)
elif mode == 'sha256sum':
    if not args: print(hashlib.sha256(sys.stdin.buffer.read()).hexdigest(), '-'); sys.exit(0)
    p = pathlib.Path(args[0])
    counter = events.with_name('hash_count')
    count = int(counter.read_text()) + 1 if counter.exists() else 1
    counter.write_text(str(count))
    if count >= 2 and os.environ.get('SOURCE_CHANGES') == '1': p.write_bytes(b'SIMPLE  =changed')
    print(hashlib.sha256(p.read_bytes()).hexdigest(), p)
elif mode == 'vmv':
    event('MOVE ' + args[0] + ' ' + args[1])
    if os.environ.get('MOVE_FAIL') == '1': sys.exit(1)
    node(args[0]).rename(node(args[1]))
elif mode == 'vcp':
    if args[0].startswith('arc:'):
        event('GET ' + args[0])
        if args[0].endswith('.fits'): raise RuntimeError('FITS downloads are forbidden')
        pathlib.Path(args[1]).write_bytes(node(args[0]).read_bytes())
    else:
        src = next(pathlib.Path(args[0]).iterdir()); dest = node(args[1]) / pathlib.Path(args[0]).name / src.name
        if dest.exists(): event('SKIP')
        else:
            event('UPLOAD')
            dest.write_bytes(src.read_bytes())
''')
    fake.chmod(0o755)
    for name in ("vcp", "vmv", "stat", "sha256sum", "cp", "flock", "cadc-get-cert"):
        (tools / name).symlink_to(fake)
    events = tmp_path / "events"
    environment = os.environ.copy()
    environment.update(PATH=str(tools) + os.pathsep + environment["PATH"],
        FAKE_REMOTE=str(remote), FAKE_SOURCE=str(source), FAKE_EVENTS=str(events),
        VCP_CMD=str(tools / "vcp"), VMV_CMD=str(tools / "vmv"),
        CHECKSUM_PYTHON=sys.executable, CHECKSUM_HELPER=str(HELPER),
        CHECKSUM_RECEIPT_DIR=str(tmp_path / "receipts"), worker_dir=str(tmp_path))

    def regenerate():
        manifest_helper.generate(SimpleNamespace(root=str(remote), run="7000",
            galaxies=["NGC4698"], cont_only=True, output=None))

    regenerate()
    base_overlay = tmp_path / 'base.img'
    base_overlay.write_bytes(b'offline overlay')
    environment.update(BASE_OVERLAY=str(base_overlay), OVERLAY_DIR=str(tmp_path / 'overlays'),
        SOURCE_7000=str(source.parent.parent), DEST_7000='arc:products/v3tk_v7.6.8_7000')

    def run(**settings):
        functions = "\n".join(re.findall(r"(?ms)^\w+\(\) \{.*?^\}", SCRIPT.read_text()))
        harness = functions + '''
FILE_RETRIES=1
RETRY_BASE_SLEEP=1
RETRY_MAX_SLEEP=1
RUN_ID=offline_test
"$VCP_CMD" arc:products/v3tk_v7.6.8_7000/.checksum_manifest.json "$FAKE_STAGE/../checksum_manifest.json" || exit 1
upload_one_file 7000 NGC4698 "$FAKE_SOURCE" arc:products/v3tk_v7.6.8_7000 overlay "$FAKE_STAGE"
'''
        env = environment | {"FAKE_STAGE": str(stage)} | settings
        result = subprocess.run(["/bin/bash", "-c", harness], env=env, capture_output=True, text=True)
        actions = events.read_text().splitlines() if events.exists() else []
        assert not any(action.startswith("GET ") and action.endswith(".fits") for action in actions), result.stdout + result.stderr
        return result, actions

    return run, source, target, regenerate


@pytest.fixture
def full_transfer(transfer):
    # Reuse the real fixture bytes/tools, but execute the entire shell script.
    run, source, target, regenerate = transfer
    tmp_path = source.parent.parent.parent
    tools = tmp_path / 'tools'
    # Authentication is mocked; this test must never inspect a real certificate.
    python_wrapper = tools / 'cadc-python'
    python_wrapper.write_text('#!/bin/bash\ncat >/dev/null\nexit 0\n')
    python_wrapper.chmod(0o700)
    environment = os.environ.copy()
    environment.update(PATH=str(tools) + os.pathsep + environment['PATH'],
        FAKE_REMOTE=str(tmp_path / 'remote'), FAKE_SOURCE=str(source),
        FAKE_EVENTS=str(tmp_path / 'events'), VCP_CMD=str(tools / 'vcp'), VMV_CMD=str(tools / 'vmv'),
        CHECKSUM_PYTHON=sys.executable, CHECKSUM_HELPER=str(HELPER),
        CHECKSUM_RECEIPT_DIR=str(tmp_path / 'receipts'), BASE_OVERLAY=str(tmp_path / 'base.img'),
        OVERLAY_DIR=str(tmp_path / 'overlays'), SOURCE_7000=str(source.parent.parent),
        DEST_7000='arc:products/v3tk_v7.6.8_7000', FILE_RETRIES='1',
        CADC_PYTHON_CMD=str(python_wrapper))
    def execute():
        result = subprocess.run(['/bin/bash', str(SCRIPT), '--manual-checksum', '--cont-only', '7000', 'NGC4698'],
            env=environment, capture_output=True, text=True)
        actions = (tmp_path / 'events').read_text().splitlines()
        assert not any(a.startswith('GET ') and a.endswith('.fits') for a in actions)
        return result
    return execute, regenerate


def test_corrupt_manifest_entry_is_replaced_without_fits_download(transfer):
    run, source, target, _ = transfer
    result, actions = run()
    assert result.returncode == 2, result.stdout + result.stderr
    assert target.read_bytes() == source.read_bytes()
    backups = list(target.parent.glob("*.corrupt_*"))
    assert len(backups) == 1 and backups[0].read_bytes() == b"damaged-existing-file"
    assert actions.count("UPLOAD") == 1
    assert 'AWAITING VERIFICATION' in result.stdout


def test_matching_manifest_skips_transfer_and_does_not_need_vmv(transfer):
    run, source, target, regenerate = transfer
    target.write_bytes(source.read_bytes())
    regenerate()
    result, actions = run(VMV_CMD="/nonexistent/vmv")
    assert result.returncode == 0, result.stdout + result.stderr
    assert len(actions) == 1 and actions[0].startswith('GET ')


def test_missing_remote_file_is_uploaded_without_quarantine(transfer):
    run, source, target, regenerate = transfer
    target.unlink()
    regenerate()
    result, actions = run(VMV_CMD="/nonexistent/vmv")
    assert result.returncode == 2, result.stdout + result.stderr
    assert target.read_bytes() == source.read_bytes()
    assert not any(action.startswith('MOVE ') for action in actions)


def test_reusing_preupload_manifest_does_not_repeat_repair(transfer):
    run, _, _, _ = transfer
    assert run()[0].returncode == 2
    result, actions = run()
    assert result.returncode == 2, result.stdout + result.stderr
    assert actions.count('UPLOAD') == 1
    assert sum(action.startswith('MOVE ') for action in actions) == 1


def test_fresh_manifest_verifies_replacement_without_reupload(transfer):
    run, _, _, regenerate = transfer
    assert run()[0].returncode == 2
    regenerate()
    result, actions = run()
    assert result.returncode == 0, result.stdout + result.stderr
    assert actions.count('UPLOAD') == 1
    assert 'VERIFIED BY CANFAR MANIFEST' in result.stdout


def test_source_change_prevents_quarantine(transfer):
    run, _, target, _ = transfer
    result, actions = run(SOURCE_CHANGES='1')
    assert result.returncode == 1
    assert target.read_bytes() == b'damaged-existing-file'
    assert not any(action.startswith('MOVE ') for action in actions)


def test_quarantine_failure_prevents_upload(transfer):
    run, _, target, _ = transfer
    result, actions = run(MOVE_FAIL='1')
    assert result.returncode == 1
    assert target.read_bytes() == b'damaged-existing-file'
    assert 'UPLOAD' not in actions


def test_manifest_lookup_rejects_wrong_run_or_uncovered_file(transfer):
    _, _, target, _ = transfer
    manifest = target.parent.parent / '.checksum_manifest.json'
    for run, filename in [('normal', target.name), ('7000', 'NGC4698_sfh_maps.fits')]:
        with pytest.raises(ValueError):
            manifest_helper.lookup(SimpleNamespace(manifest=str(manifest), run=run,
                galaxy='NGC4698', filename=filename))


def test_entire_script_propagates_pending_then_verified_exit_codes(full_transfer):
    execute, regenerate = full_transfer
    result = execute()
    assert result.returncode == 2, result.stdout + result.stderr
    assert 'awaiting a fresh CANFAR-generated manifest' in result.stdout
    regenerate()
    result = execute()
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'All selected available files verified' in result.stdout
