"""Check split-run routing with actual Mass path helpers, without science runs."""
import ast
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

SOURCE = Path(__file__).resolve().parents[1]
SCRIPT = SOURCE / 'QC_ngist_v3tk_v768_7000_SFH16.sh'
MASS_STANDIN = r'''
import ast, os, sys
from pathlib import Path
a = sys.argv[1:]
g = a[a.index('-g')+1]
root = Path(a[a.index('--root')+1])
with open(os.environ['TRACE'], 'a') as f:
    f.write(repr(('mass', a)) + '\n')
assert Path(a[a.index('--fallback-root')+1]) == root
if g in {'NGC4254', 'NGC4321', 'NGC4535'}:
    cube_root = Path(a[a.index('--phangs-cube-root')+1])
    name = g+'_PHANGS_DATACUBE_native.fits'
    assert (cube_root/name).is_file() or (cube_root/(name+'.gz')).is_file()
assert Path(a[a.index('--ngist-config')+1]).resolve() == (Path(os.environ['PRODUCTS'])/(g+'_SFH16')/'CONFIG').resolve()
names = {'_unique_paths', '_has_glob_chars', '_fits_path_alternates',
         '_glob_candidate_patterns', '_candidate_paths', '_input_search_dirs',
         '_keyword_sort_key', '_format_keyword', 'extract_shared_keyword',
         'collect_keyworded_input_paths', 'resolve_keyworded_input_group',
         'build_further_output_path'}
tree = ast.parse(Path(os.environ['REAL_MASS']).read_text())
tree.body = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in names]
exec(compile(tree, os.environ['REAL_MASS'], 'exec'))
keyword, paths = resolve_keyworded_input_group('test products', root, root, Path(g), g, {
    'bin': ['SPATIAL_BINNING_maps.fits', 'spatial_binning_maps.fits'],
    'sfh': ['SFH_maps.fits', 'sfh_maps.fits'],
    'weights': ['sfh-weights.fits', 'sfh_weights.fits']})
assert keyword == 'SFH16', (keyword, paths)
assert paths['bin'].read_bytes() == b'ordinary:spatial_binning_maps'
assert paths['sfh'].read_bytes() == b'new:sfh_maps'
assert paths['weights'].read_bytes() == b'new:sfh_weights'
if os.environ.get('FAIL_MASS') == g: sys.exit(9)
if os.environ.get('NO_MASS_OUTPUT') == g: sys.exit(0)
out = root/g/build_further_output_path(paths['bin']).name
assert out.name == g+'_SFH16_spatial_binning_maps_further.fits', out
out.write_bytes(b'new:mass')
'''
QC_STANDIN = r'''
import os, sys
from pathlib import Path
version = 'v3tk_v7.6.8_7000'
priority_data_folder = os.environ['PRODUCTS']
fallback_data_folder = os.path.join(os.getcwd(), version)
data_folder = priority_data_folder if os.path.isdir(priority_data_folder) else fallback_data_folder
galaxyid = sys.argv[1]
if __name__ == '__main__':
    with open(os.environ['TRACE'], 'a') as f:
        f.write(repr(('qc', sys.argv[1:], data_folder)) + '\n')
    g = galaxyid.removesuffix('_SFH16')
    folder = Path(data_folder)/galaxyid
    assert folder.resolve() != (Path(os.environ['PRODUCTS'])/galaxyid).resolve()
    for suffix in ('kin_maps', 'gas_bin_maps', 'gas_spaxel_maps'):
        assert (folder/(galaxyid+'_'+suffix+'.fits')).read_bytes() == ('ordinary:'+suffix).encode()
    assert (folder/(galaxyid+'_sfh_maps.fits')).read_bytes() == b'new:sfh_maps'
    assert (folder/(galaxyid+'_spatial_binning_maps_further.fits')).read_bytes() == b'new:mass'
    if os.environ.get('FAIL_QC') == g: sys.exit(7)
    (Path.cwd()/'QC'/(galaxyid+'_'+version+'_QC.pdf')).write_bytes(b'new:pdf')
'''


class SFH16WrapperTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.products = self.root/'products with spaces'
        self.products.mkdir()
        self.trace = self.root/'trace'
        self.bin = self.root/'scripts'
        self.bin.mkdir()
        self.script = self.bin/SCRIPT.name
        shutil.copy2(SCRIPT, self.script)
        (self.bin/'Mass.py').write_text(MASS_STANDIN)
        (self.bin/'QC_ngist_v3tk_v768_7000.py').write_text(QC_STANDIN)
        self.scratch = self.root/'scratch'
        self.scratch.mkdir()
        self.env = dict(os.environ, PYTHON_BIN=sys.executable, TRACE=str(self.trace),
                        PRODUCTS=str(self.products), REAL_MASS=str(SOURCE/'Mass.py'),
                        TMPDIR=str(self.scratch))

    def seed(self, g, variants=False):
        run = g+'_SFH16'
        folder = self.products/run
        folder.mkdir()
        (folder/'CONFIG').write_text('SFH16 configuration')
        (folder/'LOGFILE').write_text('nGIST completed successfully.\n')
        for suffix in ('sfh_maps', 'sfh_weights'):
            (folder/f'{run}_{suffix}.fits').write_bytes(('new:'+suffix).encode())
        original = self.products/g
        original.mkdir()
        for suffix in ('spatial_binning_maps', 'kin_maps', 'gas_bin_maps', 'gas_spaxel_maps'):
            name = f'{g}_{suffix}.fits'
            if variants:
                name = name.replace('gas_bin', 'gas_BIN').replace('gas_spaxel', 'gas_SPAXEL')+'.gz'
            (original/name).write_bytes(('ordinary:'+suffix).encode())
        (original/f'{g}_spatial_binning_maps_further.fits').write_bytes(b'original:mass')

    def run_wrapper(self, *args):
        return subprocess.run(['bash', str(self.script), '--root', str(self.products),
                               '--output-dir', str(self.root/'QC'), *args],
                              env=self.env, text=True, capture_output=True)

    def calls(self):
        if not self.trace.exists():
            return []
        return [ast.literal_eval(s) for s in self.trace.read_text().splitlines()]

    def assert_scratch_clean(self):
        self.assertEqual(list(self.scratch.iterdir()), [])

    def test_combined_id_reads_split_sources_and_writes_sfh16_outputs(self):
        g = 'NGC4567_8'
        self.seed(g)
        before = {p: p.read_bytes() for p in self.products.rglob('*') if p.is_file()}
        r = self.run_wrapper(g+'_SFH16')
        self.assertEqual(r.returncode, 0, r.stdout+r.stderr)
        calls = self.calls()
        self.assertEqual([c[0] for c in calls], ['mass', 'qc'])
        self.assertEqual(calls[0][1][calls[0][1].index('-g')+1], g)
        self.assertEqual(calls[1][1], [g+'_SFH16'])
        for p, content in before.items():
            self.assertEqual(p.read_bytes(), content, p)
        run = self.products/(g+'_SFH16')
        self.assertEqual((run/(g+'_SFH16_spatial_binning_maps_further.fits')).read_bytes(), b'new:mass')
        self.assertEqual(len(list(run.iterdir())), 5)
        self.assertTrue((self.root/'QC'/(g+'_SFH16_v3tk_v7.6.8_7000_QC.pdf')).is_file())
        self.assert_scratch_clean()

    def test_old_sfh16_duplicates_are_ignored_and_preserved(self):
        g = 'NGC4189'
        self.seed(g)
        for suffix in ('spatial_binning_maps', 'kin_maps', 'gas_bin_maps', 'gas_spaxel_maps'):
            (self.products/(g+'_SFH16')/f'{g}_SFH16_{suffix}.fits').write_bytes(b'stale duplicate')
        r = self.run_wrapper(g)
        self.assertEqual(r.returncode, 0, r.stdout+r.stderr)
        self.assertEqual((self.products/(g+'_SFH16')/f'{g}_SFH16_kin_maps.fits').read_bytes(), b'stale duplicate')

    def test_compressed_and_uppercase_gas_source_names(self):
        self.seed('NGC4189', variants=True)
        r = self.run_wrapper('NGC4189')
        self.assertEqual(r.returncode, 0, r.stdout+r.stderr)

    def test_temporary_binning_copy_when_hardlinks_unavailable(self):
        self.seed('NGC4189')
        fakebin = self.root/'fakebin'
        fakebin.mkdir()
        (fakebin/'ln').write_text('#!/usr/bin/env bash\n[[ "$1" == -s ]] || exit 1\nexec /bin/ln "$@"\n')
        (fakebin/'ln').chmod(0o755)
        self.env['PATH'] = str(fakebin)+os.pathsep+os.environ['PATH']
        r = self.run_wrapper('NGC4189')
        self.assertEqual(r.returncode, 0, r.stdout+r.stderr)
        self.assert_scratch_clean()

    def test_original_binning_may_itself_be_a_symlink(self):
        g = 'NGC4189'
        self.seed(g)
        binning = self.products/g/f'{g}_spatial_binning_maps.fits'
        target = self.products/g/'inherited_binning.fits'
        binning.rename(target)
        binning.symlink_to(target)
        r = self.run_wrapper(g)
        self.assertEqual(r.returncode, 0, r.stdout+r.stderr)
        self.assert_scratch_clean()

    def test_mass_failure_prevents_qc_and_preserves_existing_output(self):
        g = 'NGC4189'
        self.seed(g)
        output = self.products/(g+'_SFH16')/f'{g}_SFH16_spatial_binning_maps_further.fits'
        output.write_bytes(b'previous mass')
        self.env['FAIL_MASS'] = g
        r = self.run_wrapper(g)
        self.assertNotEqual(r.returncode, 0)
        self.assertEqual([c[0] for c in self.calls()], ['mass'])
        self.assertEqual(output.read_bytes(), b'previous mass')
        self.assert_scratch_clean()

    def test_success_without_fresh_mass_output_does_not_use_stale_product(self):
        g = 'NGC4189'
        self.seed(g)
        output = self.products/(g+'_SFH16')/f'{g}_SFH16_spatial_binning_maps_further.fits'
        output.write_bytes(b'previous mass')
        self.env['NO_MASS_OUTPUT'] = g
        r = self.run_wrapper(g)
        self.assertNotEqual(r.returncode, 0)
        self.assertEqual([c[0] for c in self.calls()], ['mass'])
        self.assertEqual(output.read_bytes(), b'previous mass')
        self.assert_scratch_clean()

    def test_qc_failure_continues_to_next_run_and_cleans_view(self):
        self.seed('NGC4189')
        self.seed('NGC4380')
        self.env['FAIL_QC'] = 'NGC4189'
        r = self.run_wrapper('NGC4189', 'NGC4380')
        self.assertNotEqual(r.returncode, 0)
        self.assertIn('1 passed, 0 skipped, 1 failed', r.stdout)
        self.assertEqual(len(self.calls()), 4)
        self.assert_scratch_clean()

    def test_missing_neighbour_cannot_fall_back_to_duplicate(self):
        g = 'NGC4189'
        self.seed(g)
        (self.products/g/f'{g}_kin_maps.fits').rename(self.products/(g+'_SFH16')/f'{g}_SFH16_kin_maps.fits')
        r = self.run_wrapper(g)
        self.assertNotEqual(r.returncode, 0)
        self.assertFalse(self.trace.exists())
        self.assertIn(str(self.products/g), r.stdout+r.stderr)

    def test_missing_sfh16_config_does_not_borrow_original(self):
        g = 'NGC4189'
        self.seed(g)
        (self.products/(g+'_SFH16')/'CONFIG').unlink()
        (self.products/g/'CONFIG').write_text('ordinary configuration')
        r = self.run_wrapper(g)
        self.assertNotEqual(r.returncode, 0)
        self.assertFalse(self.trace.exists())
        self.assertIn('CONFIG', r.stdout+r.stderr)

    def test_explicit_missing_inputs_fail_without_invocation(self):
        r = self.run_wrapper('NGC4579')
        self.assertNotEqual(r.returncode, 0)
        self.assertFalse(self.trace.exists())

    def test_default_sample_skips_missing_and_processes_ready(self):
        self.seed('NGC4189')
        r = self.run_wrapper()
        self.assertEqual(r.returncode, 0, r.stdout+r.stderr)
        self.assertIn('38 skipped', r.stdout)
        self.assertEqual(len(self.calls()), 2)

    def test_dry_run_is_read_only_and_lists_full_sample_and_both_sources(self):
        r = self.run_wrapper('--dry-run')
        self.assertEqual(r.returncode, 0, r.stdout+r.stderr)
        self.assertIn('NGC4579_SFH16', r.stdout)
        self.assertIn('39 cube IDs', r.stdout)
        self.assertIn('REUSE', r.stdout)
        self.assertFalse(self.trace.exists())
        self.assertFalse((self.root/'QC').exists())
        self.assert_scratch_clean()

    def test_unknown_id_rejected_before_any_work(self):
        self.seed('NGC4189')
        r = self.run_wrapper('NGC4189', 'NGC9999')
        self.assertNotEqual(r.returncode, 0)
        self.assertFalse(self.trace.exists())

    def test_duplicate_plain_and_suffixed_id_runs_once(self):
        self.seed('NGC4189')
        r = self.run_wrapper('NGC4189', 'NGC4189_SFH16')
        self.assertEqual(r.returncode, 0, r.stdout+r.stderr)
        self.assertEqual(len(self.calls()), 2)

    def setup_phangs(self):
        self.seed('NGC4254')
        stage = self.root/'stage'
        stage.mkdir()
        self.env['PHANGS_STAGING_ROOT'] = str(stage)
        fakebin = self.root/'vcpbin'
        fakebin.mkdir()
        vcp = fakebin/'vcp'
        vcp.write_text('#!/usr/bin/env bash\nprintf "%s\\n" "$1" >> "$VCP_TRACE"\nprintf cube > "$2"\nexit "${VCP_EXIT:-0}"\n')
        vcp.chmod(0o755)
        self.env['VCP_TRACE'] = str(self.root/'vcp_trace')
        self.env['PATH'] = str(fakebin)+os.pathsep+os.environ['PATH']
        return stage

    def test_phangs_download_passed_to_mass_and_removed(self):
        stage = self.setup_phangs()
        r = self.run_wrapper('--cube-root', str(self.root/'missing'), 'NGC4254')
        self.assertEqual(r.returncode, 0, r.stdout+r.stderr)
        self.assertIn('vos:phangs/RELEASES/PHANGS-MUSE/DR1.0/DATACUBES/NGC4254_PHANGS_DATACUBE_native.fits', (self.root/'vcp_trace').read_text())
        self.assertEqual(list(stage.iterdir()), [])

    def test_phangs_failed_download_blocks_mass_and_cleans_partial(self):
        stage = self.setup_phangs()
        self.env['VCP_EXIT'] = '8'
        r = self.run_wrapper('--cube-root', str(self.root/'missing'), 'NGC4254')
        self.assertNotEqual(r.returncode, 0)
        self.assertEqual(self.calls(), [])
        self.assertEqual(list(stage.iterdir()), [])

    def test_phangs_local_compressed_cube_is_preserved_without_download(self):
        stage = self.setup_phangs()
        cubes = self.root/'native cubes'
        cubes.mkdir()
        cube = cubes/'NGC4254_PHANGS_DATACUBE_native.fits.gz'
        cube.write_bytes(b'existing cube')
        self.env['PHANGS_CUBE_ROOT'] = str(cubes)
        r = self.run_wrapper('NGC4254')
        self.assertEqual(r.returncode, 0, r.stdout+r.stderr)
        self.assertEqual(cube.read_bytes(), b'existing cube')
        self.assertFalse((self.root/'vcp_trace').exists())
        self.assertEqual(list(stage.iterdir()), [])

    def test_phangs_staged_cube_removed_after_mass_failure(self):
        stage = self.setup_phangs()
        self.env['FAIL_MASS'] = 'NGC4254'
        r = self.run_wrapper('--cube-root', str(self.root/'missing'), 'NGC4254')
        self.assertNotEqual(r.returncode, 0)
        self.assertEqual([c[0] for c in self.calls()], ['mass'])
        self.assertEqual(list(stage.iterdir()), [])
        self.assert_scratch_clean()


if __name__ == '__main__':
    unittest.main()
