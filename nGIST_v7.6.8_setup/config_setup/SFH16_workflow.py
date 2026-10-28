"""Generate SFH16 configs and stage immutable original SFH inputs.

Checked against upstream nGIST 7.6.8, commit 038e9fa4839758e9dae53baeaf93a82bd75d8026.
No SFH outputs are copied: their presence would cause the fit to be skipped.
"""
import copy
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile

import yaml

PRODUCTS = Path(os.environ.get("SFH16_PRODUCTS_DIR", "/scratch/pawsey1308/mauve/products/v3tk_v7.6.8_7000"))
HERE = Path(__file__).resolve().parent
CONFIG_DIR = Path(os.environ.get("SFH16_CONFIG_DIR", HERE / "configFiles" if (HERE / "configFiles").is_dir() else HERE))
DISABLED = ("SPATIAL_MASKING", "SPATIAL_BINNING", "PREPARE_SPECTRA", "KIN", "CONT", "GAS", "LS", "UMOD")
SUFFIXES = ("_mask.fits", "_table.fits", "_bin_spectra.hdf5", "_kin.fits")
MANIFEST = "SFH16_INPUTS.json"


def load(path):
    with path.open() as stream:
        return yaml.safe_load(stream)


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def config_path(galid, special=True):
    suffix = "_SFH16" if special else ""
    return CONFIG_DIR / f"{galid}_MAUVE_MasterConfig_v7.6.8_7000{suffix}_setonix.yaml"


def validate_sfh(config, galid, special):
    require(config["GENERAL"]["RUN_ID"] == galid + ("_SFH16" if special else ""), "Unexpected RUN_ID")
    require(config["SFH"]["MDEG"] == (16 if special else 12), "Expected SFH MDEG 12 -> 16")
    require(config["SFH"]["FIXED"] is True, "SFH.FIXED must be true to reuse KIN")
    require(config["SFH"]["SPEC_EMICLEAN"] is False, "This workflow requires SFH.SPEC_EMICLEAN=false")
    require(config["SFH"]["METHOD"] == "ppxf", "Expected SFH ppxf")
    require(config["SFH"]["MC_PPXF"] == 0, "Expected MC_PPXF=0")


def create(galid):
    base = load(config_path(galid, False))
    validate_sfh(base, galid, False)
    require(base["GENERAL"]["OW_CONFIG"] is True and base["GENERAL"]["OW_OUTPUT"] is False,
            "Require OW_CONFIG=true and OW_OUTPUT=false")
    config = copy.deepcopy(base)
    config["GENERAL"]["RUN_ID"] = galid + "_SFH16"
    config["SFH"]["MDEG"] = 16
    for section in DISABLED:
        config[section]["METHOD"] = False
    slurm_base = HERE / f"{galid}_v3tk_v7.6.8_7000_setonix.slurm"
    text = slurm_base.read_text()
    # These jobs reuse spectra/KIN and run SFH only. Do not inherit the full
    # pipeline's GAS/CONT high-memory or long-queue resource requests.
    for option, value in [("partition", "work"), ("mem", "230G"), ("time", "24:00:00")]:
        text, count = re.subn(rf"^#SBATCH --{option}=.*$", f"#SBATCH --{option}={value}", text, flags=re.MULTILINE)
        require(count == 1, f"Expected one SBATCH --{option} in {slurm_base}")
    # Keep the original CPU request and change all job/config/log prefixes.
    text = text.replace(f"{galid}_v3tk_v7.6.8_7000", f"{galid}_v3tk_v7.6.8_7000_SFH16")
    text = text.replace(f"{galid}_MAUVE_MasterConfig_v7.6.8_7000_setonix.yaml",
                        config_path(galid).name)
    marker = "# run + stream output"
    require(marker in text, f"Missing pipeline marker in {slurm_base}")
    text = text.replace(marker, f"# Refuse missing, changed or unverified reused inputs.\nbash ./SFH16_ready.sh --check {galid}\n\n{marker}", 1)
    config_path(galid).write_text(yaml.safe_dump(config, sort_keys=False))
    output = HERE / f"{galid}_v3tk_v7.6.8_7000_SFH16_setonix.slurm"
    output.write_text(text)
    output.chmod(0o755)
    print(f"Created {config_path(galid).name} and {output.name}")


def check_config(galid, original, new):
    validate_sfh(original, galid, False)
    validate_sfh(new, galid, True)
    require(new["GENERAL"]["OW_CONFIG"] is True and new["GENERAL"]["OW_OUTPUT"] is False,
            "Require OW_CONFIG=true, OW_OUTPUT=false")
    require(new["GENERAL"].get("CREATE_HTML", False) is False, "CREATE_HTML must be false")
    require(Path(new["GENERAL"]["OUTPUT"]).resolve() == PRODUCTS.resolve(),
            "Config OUTPUT differs from preparation product root")
    for section in DISABLED:
        require(new[section]["METHOD"] is False, f"{section} must be disabled")
    # Saved CONFIG has expanded OUTPUT and extra runtime paths. All supplied
    # science settings must agree, including otherwise inactive module settings.
    for section, fields in new.items():
        for key, value in fields.items():
            if section == "GENERAL" and key in {"RUN_ID", "OUTPUT", "OW_CONFIG", "OW_OUTPUT", "NCPU", "PARALLEL", "CREATE_HTML"}:
                continue
            if section in DISABLED and key == "METHOD":
                continue
            expected = 12 if (section, key) == ("SFH", "MDEG") else value
            require(key in original.get(section, {}) and original[section][key] == expected,
                    f"Original CONFIG differs at {section}.{key}; regenerate from the actual baseline before running")


def check_dimensions(source, galid, config):
    # Check bin indexing before fixing each SFH spectrum to the matching KIN row.
    # No science fit is executed and large spectral arrays are not loaded.
    import h5py
    import numpy as np
    from astropy.io import fits

    with h5py.File(source / (galid + "_bin_spectra.hdf5"), "r") as spectra:
        shape = spectra["SPEC"].shape
        require(len(shape) == 2 and min(shape) > 0, "Invalid binned spectral dimensions")
        require(spectra["ESPEC"].shape == shape, "SPEC/ESPEC dimensions differ")
        require(spectra["LOGLAM"].shape == (shape[0],), "LOGLAM dimension differs")
        require(np.isclose(spectra.attrs["VELSCALE"], config["PREPARE_SPECTRA"]["VELSCALE"]),
                "Binned VELSCALE differs from CONFIG")
        nbins = shape[1]
    with fits.open(source / (galid + "_kin.fits"), memmap=True) as kin:
        require(len(kin[1].data) == nbins, "KIN rows differ from binned spectra")
        columns = kin[1].columns.names[:config["KIN"]["MOM"]]
        require([name.upper() for name in columns] == ["V", "SIGMA", "H3", "H4", "H5", "H6"][:config["KIN"]["MOM"]],
                "Unexpected KIN moment column order")
    with fits.open(source / (galid + "_table.fits"), memmap=True) as table:
        ids = table[1].data["BIN_ID"]
        require(np.array_equal(np.unique(ids[ids >= 0]), np.arange(nbins)),
                "Table bin IDs differ from spectral/KIN row indexing")
        nspaxels = len(ids)
    with fits.open(source / (galid + "_mask.fits"), memmap=True) as mask:
        require(len(mask[1].data) == nspaxels, "Mask/table spaxel counts differ")


def ready(galid, check_only=False):
    source = PRODUCTS / galid
    target = PRODUCTS / (galid + "_SFH16")
    if not check_only and target.exists():
        require(target.is_dir() and not target.is_symlink(), f"Expected a regular product directory: {target}")
        print(f"Skipped {target}: folder already exists; left unchanged (use --check to verify inputs)")
        return
    original_path = source / "CONFIG"
    new_path = config_path(galid)
    original, new = load(original_path), load(new_path)
    check_config(galid, original, new)
    files = {}
    for suffix in SUFFIXES:
        path = source / (galid + suffix)
        require(path.is_file() and path.stat().st_size > 0, f"Missing or empty source: {path}")
        files[suffix] = digest(path)
    check_dimensions(source, galid, new)
    expected = {"galid": galid, "source": str(source.resolve()),
                "source_config_sha256": digest(original_path),
                "run_config_sha256": digest(new_path), "inputs_sha256": files}
    if target.exists():
        require(not target.is_symlink(), f"Refusing symlink target: {target}")
        manifest = target / MANIFEST
        require(manifest.is_file(), f"Unverified existing target: {target}. Preserve/move it aside manually before preparing a fresh run.")
        require(json.loads(manifest.read_text()) == expected, f"Changed config or source inputs for {target}")
        for suffix, sha in files.items():
            path = target / (galid + "_SFH16" + suffix)
            require(not path.is_symlink() and digest(path) == sha, f"Changed staged input: {path}")
        require(digest(target / "SFH16_SOURCE_CONFIG.yaml") == expected["source_config_sha256"],
                "Changed saved source CONFIG")
        print(f"Verified {target}; existing SFH results retained")
        return
    require(not check_only, f"Run bash SFH16_ready.sh {galid} before submission")
    # Publish only a complete verified directory. Independent copies prevent
    # later changes in the baseline from altering a running SFH16 experiment.
    temp = Path(tempfile.mkdtemp(prefix=f".{galid}_SFH16-", dir=PRODUCTS))
    try:
        for suffix, sha in files.items():
            path = temp / (galid + "_SFH16" + suffix)
            shutil.copy2(source / (galid + suffix), path)
            require(digest(path) == sha, f"Source changed while copying {suffix}")
        shutil.copy2(original_path, temp / "SFH16_SOURCE_CONFIG.yaml")
        require(digest(temp / "SFH16_SOURCE_CONFIG.yaml") == expected["source_config_sha256"],
                "Source CONFIG changed during preparation")
        (temp / MANIFEST).write_text(json.dumps(expected, indent=2) + "\n")
        require(not target.exists(), f"Target appeared during preparation: {target}")
        temp.rename(target)
    finally:
        if temp.exists():
            shutil.rmtree(temp)
    print(f"Prepared {target} (four copied inputs; SFH outputs absent)")


if __name__ == "__main__":
    try:
        mode, *galids = sys.argv[1:]
        require(mode in {"create", "ready", "check"} and galids, "Use the SFH16_*.sh entry points")
        for galid in galids:
            require(galid in {"NGC4189", "NGC4380", "NGC4450", "NGC4567_8", "NGC4689"}, "Unknown SFH16 galaxy")
            if mode == "create":
                create(galid)
            else:
                ready(galid, mode == "check")
    except (ValueError, OSError, KeyError) as error:
        sys.exit(f"ERROR: {error}")
