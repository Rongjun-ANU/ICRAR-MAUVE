import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

import pytest


ROOT = Path(__file__).parent
SCRIPT = ROOT / 'vcp_scratch_v3tk_v768_7000_SFH16_to_canfar.sh'
SAMPLE = ROOT / 'nGIST_v7.6.8_setup/config_setup/27_galaxies.sh'
GALIDS = re.search(r'ALL_GALIDS=\((.*?)\)', SAMPLE.read_text(), re.S).group(1).split()
SUFFIXES = ('CONFIG', 'LOGFILE', '_sfh_maps.fits', '_sfh_weights.fits')
REUSED_SUFFIXES = ('_mask.fits', '_kin_maps.fits', '_spatial_binning_maps.fits',
                   '_gas_bin_maps.fits', '_gas_spaxel_maps.fits')


def run_script(tmp_path, *args):
    env = dict(os.environ, SOURCE_7000=str(tmp_path / 'source'),
               SOURCE_ORIGINAL_7000=str(tmp_path / 'original'))
    return subprocess.run(['bash', str(SCRIPT), *args], env=env,
                          capture_output=True, text=True, timeout=30)


def test_default_dry_run_covers_exact_full_sample(tmp_path):
    result = run_script(tmp_path, '--dry-run')
    assert result.returncode == 0, result.stderr
    selected = re.findall(r'^7000 (\w+)_SFH16:$', result.stdout, re.M)
    assert selected == GALIDS
    assert 'Work item count: 39' in result.stdout
    assert result.stdout.count('  Target file: ') == 39 * 4


@pytest.mark.parametrize('galid', GALIDS)
def test_every_cube_and_run_id_selects_only_four_sfh_products(tmp_path, galid):
    plain = run_script(tmp_path, '--dry-run', galid)
    suffixed = run_script(tmp_path, '--dry-run', galid + '_SFH16')
    assert plain.returncode == suffixed.returncode == 0, plain.stderr + suffixed.stderr
    assert plain.stdout == suffixed.stdout
    assert 'Work item count: 1' in plain.stdout
    for suffix in SUFFIXES:
        filename = suffix if suffix in ('CONFIG', 'LOGFILE') else galid + '_SFH16' + suffix
        assert '  Target file: ' + filename + '\n' in plain.stdout
        prefix = tmp_path / 'source'
        source_id = galid + '_SFH16'
        source_name = suffix if suffix in ('CONFIG', 'LOGFILE') else source_id + suffix
        assert str(prefix / source_id / source_name) in plain.stdout
    for suffix in REUSED_SUFFIXES:
        assert suffix not in plain.stdout


@pytest.mark.parametrize('args', [('NGC9999',), ('../IC3392',), ('IC3392', 'IC3392_SFH16')])
def test_unsupported_and_duplicate_ids_rejected(tmp_path, args):
    result = run_script(tmp_path, '--dry-run', *args)
    assert result.returncode == 2


def test_unfinished_runs_are_skipped_or_explicitly_rejected(tmp_path):
    result = run_script(tmp_path)
    assert result.returncode == 0, result.stderr
    assert result.stdout.count('Skipping ') == 39
    explicit = run_script(tmp_path, 'IC3392')
    assert explicit.returncode == 1
    assert 'not ready' in explicit.stderr


def test_readiness_needs_only_the_four_sfh_run_products(tmp_path):
    folder = tmp_path / 'source/IC3392_SFH16'
    folder.mkdir(parents=True)
    for suffix in SUFFIXES:
        name = suffix if suffix in ('CONFIG', 'LOGFILE') else 'IC3392_SFH16' + suffix
        (folder / name).write_text('nGIST completed successfully.\n')
    result = run_script(tmp_path, '--dry-run', 'IC3392')
    assert result.returncode == 0, result.stderr
    assert '  Readiness: ready\n' in result.stdout
    assert result.stdout.count('  Target file: ') == 4


def test_checksum_worker_hashes_all_39_runs(tmp_path):
    root = tmp_path / 'canfar'
    expected = {}
    for galid in GALIDS:
        run_id = galid + '_SFH16'
        folder = root / 'v3tk_v7.6.8_7000' / run_id
        folder.mkdir(parents=True)
        for suffix in SUFFIXES:
            filename = suffix if suffix in ('CONFIG', 'LOGFILE') else run_id + suffix
            data = (run_id + '/' + filename).encode()
            (folder / filename).write_bytes(data)
            expected[run_id + '/' + filename] = hashlib.sha256(data).hexdigest()
    output = tmp_path / 'manifest.json'
    result = run_script(tmp_path, '--canfar-manifest', '--root', str(root),
                        '--output', str(output), '7000')
    assert result.returncode == 0, result.stderr
    entries = json.loads(output.read_text())['entries']
    assert set(entries) == set(expected)
    assert all(entries[key]['sha256'] == value for key, value in expected.items())


def test_automatic_checksum_request_accepts_new_cube(tmp_path):
    folder = tmp_path / 'v3tk_v7.6.8_7000'
    source = folder / 'IC3392_SFH16'
    source.mkdir(parents=True)
    (source / 'CONFIG').write_text('new sample checksum fixture')
    request_id = 'a' * 32
    request = folder / ('.checksum_request_' + request_id + '.json')
    request.write_text(json.dumps(dict(request_id=request_id, run='7000',
                                      galaxy='IC3392_SFH16', filenames=['CONFIG'])))
    result = run_script(tmp_path, '--canfar-request', str(request))
    assert result.returncode == 0, result.stderr
    response = json.loads((folder / ('.checksum_response_' + request_id + '.json')).read_text())
    assert response['request_id'] == request_id
    assert response['entries']['IC3392_SFH16/CONFIG']['status'] == 'ok'


@pytest.mark.parametrize('suffix', REUSED_SUFFIXES)
def test_checksum_worker_rejects_reused_products(tmp_path, suffix):
    folder = tmp_path / 'v3tk_v7.6.8_7000'
    folder.mkdir()
    request_id = 'b' * 32
    request = folder / ('.checksum_request_' + request_id + '.json')
    request.write_text(json.dumps(dict(request_id=request_id, run='7000',
                                      galaxy='IC3392_SFH16',
                                      filenames=['IC3392_SFH16' + suffix])))
    result = run_script(tmp_path, '--canfar-request', str(request))
    assert result.returncode == 1
    assert 'Unexpected product filename' in result.stderr
    assert not (folder / ('.checksum_response_' + request_id + '.json')).exists()
