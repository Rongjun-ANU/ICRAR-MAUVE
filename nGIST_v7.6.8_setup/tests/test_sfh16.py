"""Safe staging regression tests; no real products or remote connections."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from astropy.io import fits
import h5py
import numpy as np
import yaml

CONFIG = Path(__file__).resolve().parents[1] / "config_setup"
GALIDS = ["NGC4189", "NGC4380", "NGC4450", "NGC4567_8", "NGC4689"]
SUFFIXES = ["_mask.fits", "_table.fits", "_bin_spectra.hdf5", "_kin.fits"]


class SFH16Tests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.work = Path(self.tmp.name)
        self.cfg = self.work / "configFiles"
        self.cfg.mkdir()
        self.products = self.work / "products"
        self.products.mkdir()
        for path in CONFIG.glob("*"):
            if path.is_file() and (path.name.startswith("SFH16_") or
                                  any(path.name.startswith(g + "_") for g in GALIDS)):
                shutil.copy2(path, self.cfg / path.name)
        for galid in GALIDS:
            path = self.cfg / f"{galid}_MAUVE_MasterConfig_v7.6.8_7000_setonix.yaml"
            config = yaml.safe_load(path.read_text())
            config["GENERAL"]["OUTPUT"] = str(self.products)
            path.write_text(yaml.safe_dump(config))
        self.env = dict(os.environ, PYTHON=os.sys.executable,
                        SFH16_PRODUCTS_DIR=str(self.products),
                        SFH16_CONFIG_DIR=str(self.cfg))

    def run_script(self, name, *args, ok=True):
        result = subprocess.run(["bash", str(self.cfg / name), *args],
                                env=self.env, text=True, capture_output=True)
        self.assertEqual(result.returncode == 0, ok, result.stdout + result.stderr)
        return result

    def generate(self):
        self.run_script("SFH16_creation.sh")

    def source(self, galid="NGC4689"):
        source = self.products / galid
        source.mkdir()
        base = yaml.safe_load((self.cfg / f"{galid}_MAUVE_MasterConfig_v7.6.8_7000_setonix.yaml").read_text())
        base["GENERAL"]["OUTPUT"] = str(source)
        (source / "CONFIG").write_text(yaml.safe_dump(base))
        for suffix, cols in [
            ("_mask.fits", [fits.Column(name="MASK", format="J", array=[0, 0, 0])]),
            ("_table.fits", [fits.Column(name="BIN_ID", format="J", array=[0, 1, 1])]),
            ("_kin.fits", [fits.Column(name=name, format="D", array=[0., 1.])
                           for name in ["V", "SIGMA", "H3", "H4"]]),
        ]:
            fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU.from_columns(cols)]).writeto(source / (galid + suffix))
        with h5py.File(source / (galid + "_bin_spectra.hdf5"), "w") as spectra:
            spectra["SPEC"] = np.ones((3, 2))
            spectra["ESPEC"] = np.ones((3, 2))
            spectra["LOGLAM"] = np.log([4801., 5000., 6999.])
            spectra.attrs["VELSCALE"] = 41
        # Old SFH must never get staged, or nGIST would silently skip the fit.
        (source / (galid + "_sfh.fits")).write_bytes(b"old fit")
        return source

    def test_generated_science_changes_only_mdeg(self):
        self.generate()
        disabled = {"SPATIAL_MASKING", "SPATIAL_BINNING", "PREPARE_SPECTRA", "KIN", "CONT", "GAS", "LS", "UMOD"}
        for galid in GALIDS:
            base = yaml.safe_load((self.cfg / f"{galid}_MAUVE_MasterConfig_v7.6.8_7000_setonix.yaml").read_text())
            new = yaml.safe_load((self.cfg / f"{galid}_MAUVE_MasterConfig_v7.6.8_7000_SFH16_setonix.yaml").read_text())
            base["GENERAL"]["RUN_ID"] = galid + "_SFH16"
            base["SFH"]["MDEG"] = 16
            for section in disabled:
                base[section]["METHOD"] = False
            self.assertEqual(new, base)

    def test_all_sfh16_jobs_use_work_resources(self):
        self.generate()
        for galid in GALIDS:
            text = (self.cfg / f"{galid}_v3tk_v7.6.8_7000_SFH16_setonix.slurm").read_text()
            for line in ["#SBATCH --partition=work", "#SBATCH --mem=230G",
                         "#SBATCH --time=24:00:00", "#SBATCH --cpus-per-task=128"]:
                self.assertIn(line + "\n", text)
            self.assertNotIn("#SBATCH --partition=highmem", text)
            self.assertNotIn("#SBATCH --partition=long", text)

    def test_galaxies_command_submits_all_or_selected_only_when_executed(self):
        self.generate()
        bindir = self.work / "bin"
        bindir.mkdir()
        stub = bindir / "sbatch"
        stub.write_text('#!/usr/bin/env bash\nprintf "%s\\n" "$*" >> "$SBATCH_TEST_LOG"\n')
        stub.chmod(0o755)
        log = self.work / "submitted"
        log.touch()
        self.env.update(PATH=str(bindir) + os.pathsep + self.env["PATH"],
                        SBATCH_TEST_LOG=str(log))
        self.run_script("SFH16_galaxies.sh")
        self.assertEqual(log.read_text().splitlines(),
                         [f"{g}_v3tk_v7.6.8_7000_SFH16_setonix.slurm" for g in GALIDS])
        log.write_text("")
        self.run_script("SFH16_galaxies.sh", "NGC4689", "NGC4567_8")
        self.assertEqual(log.read_text().splitlines(),
                         [f"{g}_v3tk_v7.6.8_7000_SFH16_setonix.slurm" for g in ["NGC4689", "NGC4567_8"]])
        log.write_text("")
        subprocess.run(["bash", "-c", 'source "$1"', "bash",
                        str(self.cfg / "SFH16_galaxies.sh")], env=self.env, check=True)
        self.assertEqual(log.read_text(), "")
        (self.cfg / "NGC4689_v3tk_v7.6.8_7000_SFH16_setonix.slurm").unlink()
        self.run_script("SFH16_galaxies.sh", ok=False)
        self.assertEqual(log.read_text(), "")

    def test_stage_exact_inputs_and_verify_resume(self):
        self.generate()
        source = self.source()
        self.run_script("SFH16_ready.sh", "NGC4689")
        target = self.products / "NGC4689_SFH16"
        for suffix in SUFFIXES:
            self.assertEqual((target / ("NGC4689_SFH16" + suffix)).read_bytes(),
                             (source / ("NGC4689" + suffix)).read_bytes())
        self.assertFalse((target / "NGC4689_SFH16_sfh.fits").exists())
        (target / "NGC4689_SFH16_sfh.fits").write_bytes(b"new fit")
        self.run_script("SFH16_ready.sh", "--check", "NGC4689")
        self.run_script("SFH16_ready.sh", "NGC4689")
        self.assertEqual((target / "NGC4689_SFH16_sfh.fits").read_bytes(), b"new fit")
        (target / "NGC4689_SFH16_kin.fits").write_bytes(b"changed")
        self.run_script("SFH16_ready.sh", "--check", "NGC4689", ok=False)

    def test_existing_folders_are_skipped_without_loading_sources(self):
        for galid in GALIDS:
            target = self.products / (galid + "_SFH16")
            target.mkdir()
            (target / "keep.txt").write_text("existing result")
        result = self.run_script("SFH16_ready.sh")
        self.assertEqual(result.stdout.count("Skipped "), 5)
        for galid in GALIDS:
            target = self.products / (galid + "_SFH16")
            self.assertEqual([p.name for p in target.iterdir()], ["keep.txt"])
            self.assertEqual((target / "keep.txt").read_text(), "existing result")

    def test_missing_source_rejected_and_existing_target_skipped(self):
        self.generate()
        source = self.source()
        kin_bytes = (source / "NGC4689_kin.fits").read_bytes()
        (source / "NGC4689_kin.fits").unlink()
        self.run_script("SFH16_ready.sh", "NGC4689", ok=False)
        self.assertFalse((self.products / "NGC4689_SFH16").exists())
        (source / "NGC4689_kin.fits").write_bytes(kin_bytes)
        target = self.products / "NGC4689_SFH16"
        target.mkdir()
        result = self.run_script("SFH16_ready.sh", "NGC4689")
        self.assertIn("Skipped ", result.stdout)
        self.assertEqual(list(target.iterdir()), [])
        self.run_script("SFH16_ready.sh", "--check", "NGC4689", ok=False)

    def test_wrong_bin_indexing_rejected(self):
        self.generate()
        source = self.source()
        with fits.open(source / "NGC4689_table.fits", mode="update") as hdus:
            hdus[1].data["BIN_ID"] = [0, 2, 2]
        result = self.run_script("SFH16_ready.sh", "NGC4689", ok=False)
        self.assertIn("Table bin IDs", result.stderr)
        self.assertFalse((self.products / "NGC4689_SFH16").exists())

    def test_send_transfers_selected_jobs_and_preparation(self):
        self.generate()
        bindir = self.work / "bin"
        bindir.mkdir()
        stub = bindir / "ssh"
        stub.write_text("#!/usr/bin/env bash\nset -eu\n"
                        'echo "$*" >> "$SSH_TEST_DIR/calls"\n'
                        'if [[ "$*" == *"tar -xf"* ]]; then\n'
                        '  payload="$SSH_TEST_DIR/payload_$(wc -l < "$SSH_TEST_DIR/calls").tar"\n'
                        '  cat > "$payload"\nfi\n')
        stub.chmod(0o755)
        self.env.update(PATH=str(bindir) + os.pathsep + self.env["PATH"], SSH_TEST_DIR=str(self.work))
        self.run_script("SFH16_send.sh", "NGC4567_8")
        import tarfile
        contents = []
        for archive in sorted(self.work.glob("payload_*.tar")):
            with tarfile.open(archive) as stream:
                contents.extend(stream.getnames())
        self.assertEqual(set(contents), {
            "NGC4567_8_MAUVE_MasterConfig_v7.6.8_7000_SFH16_setonix.yaml",
            "NGC4567_8_v3tk_v7.6.8_7000_SFH16_setonix.slurm",
            "SFH16_ready.sh", "SFH16_galaxies.sh", "SFH16_workflow.py"})
        self.assertIn("rhuang@setonix.pawsey.org.au", (self.work / "calls").read_text())

    def test_source_config_mismatch_rejected(self):
        self.generate()
        source = self.source()
        config = yaml.safe_load((source / "CONFIG").read_text())
        config["SFH"]["SPEC_MASK"] = "different_mask"
        (source / "CONFIG").write_text(yaml.safe_dump(config))
        self.run_script("SFH16_ready.sh", "NGC4689", ok=False)
        self.assertFalse((self.products / "NGC4689_SFH16").exists())


if __name__ == "__main__":
    unittest.main()
