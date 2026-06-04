import os
import shutil
import subprocess
import sys
import tempfile
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import yaml
from astropy.io import fits
from astropy.table import Table


ROOT = Path(__file__).resolve().parents[1]
CONFIG_DIR = ROOT / "config_setup"
MASTER_CONFIG = "MAUVE_MasterConfig_v7.6.8_setonix.yaml"
PHANGS_VOS_DIR = "vos:phangs/RELEASES/PHANGS-MUSE/DR1.0/DATACUBES"


def write_minimal_master_config(path: Path):
    path.write_text(
        "\n".join(
            [
                "GENERAL:",
                "  RUN_ID: ''",
                "  INPUT: ''",
                "  OUTPUT: ''",
                "  REDSHIFT: 0",
                "  NCPU: 0",
                "READ_DATA:",
                "  ORIGIN: ''",
                "  EBmV: 0",
                "SPATIAL_MASKING:",
                "  MASK: ''",
                "KIN:",
                "  SIGMA: 0",
                "CONT:",
                "  SIGMA: 0",
                "",
            ]
        ),
        encoding="utf-8",
    )


def write_setup_table(path: Path):
    table = Table(
        {
            "Galaxy": ["NGC4254", "NGC4321", "NGC4535"],
            "z": [0.008026, 0.005240, 0.006551],
            "EBV": [0.0334, 0.0225, 0.0166],
            "sigma": [128.0, 154.0, 130.0],
        }
    )
    table.write(path, format="fits", overwrite=True)


def populate_workdir(workdir: Path):
    shutil.copy2(CONFIG_DIR / "make_gist_config_try.py", workdir / "make_gist_config_try.py")
    write_minimal_master_config(workdir / MASTER_CONFIG)
    write_setup_table(workdir / "GIST_setupinput_v1.fits")


@contextmanager
def make_workdir():
    with tempfile.TemporaryDirectory() as tmpname:
        workdir = Path(tmpname)
        populate_workdir(workdir)
        yield workdir


def run_make_config(workdir: Path, cube_dir: Path, galid: str):
    env = dict(os.environ)
    env["MAUVE_CUBE_DIR"] = str(cube_dir)
    env["MAUVE_PRODUCTS_DIR"] = "/scratch/pawsey1308/mauve/products/v3tk_v7.6.8"
    result = subprocess.run(
        [sys.executable, str(workdir / "make_gist_config_try.py"), galid],
        cwd=workdir,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    return yaml.safe_load((workdir / f"{galid}_{MASTER_CONFIG}").read_text(encoding="utf-8"))


def test_phangs_native_galaxies_use_native_cube_paths():
    with make_workdir() as workdir:
        (workdir / "cube_centers_v3tk.csv").write_text(
            "ID,file,shape,nz,ny,nx,center_y,center_x,status\n"
            "NGC4254,/scratch/cube.fits,\"(1, 20, 40)\",1,20,40,10,20,ok\n",
            encoding="utf-8",
        )
        cube_dir = workdir / "cubes"
        cube_dir.mkdir()

        config = run_make_config(workdir, cube_dir, "NGC4254")

        assert config["GENERAL"]["INPUT"] == str(cube_dir / "NGC4254_PHANGS_DATACUBE_native.fits")
        assert config["READ_DATA"]["ORIGIN"] == "20,10"


def test_missing_center_csv_row_falls_back_to_cube_midpoint():
    with make_workdir() as workdir:
        (workdir / "cube_centers_v3tk.csv").write_text(
            "ID,file,shape,nz,ny,nx,center_y,center_x,status\n",
            encoding="utf-8",
        )
        cube_dir = workdir / "cubes"
        cube_dir.mkdir()
        fits.HDUList(
            [
                fits.PrimaryHDU(),
                fits.ImageHDU(
                    data=np.zeros((2, 5, 6), dtype=np.float32),
                    name="DATA",
                ),
            ]
        ).writeto(cube_dir / "NGC4321_PHANGS_DATACUBE_native.fits")

        config = run_make_config(workdir, cube_dir, "NGC4321")

        assert config["GENERAL"]["INPUT"] == str(cube_dir / "NGC4321_PHANGS_DATACUBE_native.fits")
        assert config["READ_DATA"]["ORIGIN"] == "3,3"


def test_vcp_transfer_script_lists_public_native_phangs_cubes_only():
    script = (CONFIG_DIR / "vcp_from_v3tk_to_scratch.sh").read_text(encoding="utf-8")

    for galid in ("NGC4254", "NGC4321", "NGC4535"):
        assert f"{PHANGS_VOS_DIR}/{galid}_PHANGS_DATACUBE_native.fits" in script
    assert "copt" not in script.lower()


def test_vcp_transfer_script_accepts_selected_galaxies():
    script = (CONFIG_DIR / "vcp_from_v3tk_to_scratch.sh").read_text(encoding="utf-8")

    assert "source_for_galid()" in script
    assert "REQUESTED_GALIDS" in script
    assert "V3TK_CUBE_SUFFIX='_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk.fits.gz'" in script
    assert 'printf \'%s%s%s\\n\' "$SRC_PREFIX" "$galid" "$V3TK_CUBE_SUFFIX"' in script
    assert "--dry-run" in script


def test_vcp_transfer_script_caps_workers_at_five():
    script = (CONFIG_DIR / "vcp_from_v3tk_to_scratch.sh").read_text(encoding="utf-8")

    assert "MAX_JOBS=5" in script
    assert "JOBS=${JOBS:-5}" in script
    assert "Capping JOBS" in script
    assert "EFFECTIVE_JOBS" in script


def test_cube_check_scripts_default_to_v3tk_and_phangs_patterns():
    for script_name in ("check_v3tk_cube_centers.sh", "check_v3tk_cube_sizes.sh"):
        script = (CONFIG_DIR / script_name).read_text(encoding="utf-8")

        assert "DEFAULT_GLOB_PATTERNS=\"*_DATACUBE*.fits,*_DATACUBE*.fits.gz\"" in script
        assert "GLOB_PATTERNS=\"${2:-$DEFAULT_GLOB_PATTERNS}\"" in script
        assert "for part in pattern_arg.split(',')" in script


def test_cube_check_scripts_strip_phangs_from_cube_id():
    for script_name in ("check_v3tk_cube_centers.sh", "check_v3tk_cube_sizes.sh"):
        script = (CONFIG_DIR / script_name).read_text(encoding="utf-8")

        assert "def cube_id_from_filename(file_name: str) -> str:" in script
        assert 'cube_id = file_name.split("_DATACUBE", 1)[0]' in script
        assert 'if cube_id.endswith("_PHANGS"):' in script
        assert 'cube_id = cube_id[: -len("_PHANGS")]' in script


def test_cube_check_scripts_use_run_local_ngist_overlays():
    for script_name in ("check_v3tk_cube_centers.sh", "check_v3tk_cube_sizes.sh"):
        script = (CONFIG_DIR / script_name).read_text(encoding="utf-8")

        assert "BASE_OVERLAY=" in script
        assert "RUN_OVERLAY=" in script
        assert 'cp --reflink=auto "$BASE_OVERLAY" "$RUN_OVERLAY"' in script
        assert 'wait_for_overlay "$BASE_OVERLAY" || exit 1' in script
        assert 'wait_for_overlay "$RUN_OVERLAY" || exit 1' in script
        assert 'NGIST_OVERLAY="$RUN_OVERLAY"' in script
        assert 'rm -f "$RUN_OVERLAY"' in script


if __name__ == "__main__":
    for name, func in sorted(globals().items()):
        if name.startswith("test_") and callable(func):
            func()
            print(f"PASS {name}")
