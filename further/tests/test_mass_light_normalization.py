"""Regression checks for the 2026-10-29 nGIST LIGHT normalization correction."""
import sys
import ast
import hashlib
import re
import scipy
import astropy
import ppxf as ppxf_package
from pathlib import Path

import numpy as np
import pytest
import yaml
from astropy.io import fits
from ppxf.ppxf_util import gaussian_filter1d, log_rebin
from scipy.interpolate import interp1d

ROOT = Path(__file__).resolve().parents[1]
MASS_PATH = ROOT / "Mass.py"
# Mass.py retains its existing executable layout. Compile its actual function
# definitions, rather than importing it and starting a science pipeline.
_function_names = {"sha256_file", "parse_template_filename_key", "_resource_path",
                   "load_run_setup", "load_lsf", "build_template_norm_flux_lookup",
                   "r_light_mass_to_light"}
_constant_names = {"NGIST_REVISION", "TEMPLATE_VELSCALE_RATIO", "_TEMPLATE_RE"}
_tree = ast.parse(MASS_PATH.read_text())
_nodes = [node for node in _tree.body
          if isinstance(node, ast.FunctionDef) and node.name in _function_names
          or isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id in _constant_names for t in node.targets)]
_namespace = dict(globals(), __file__=str(MASS_PATH))
exec(compile(ast.Module(body=_nodes, type_ignores=[]), str(MASS_PATH), "exec"), _namespace)
globals().update({name: _namespace[name] for name in _function_names})


def test_normalization_is_implemented_inside_mass_py():
    root = Path(__file__).resolve().parents[1]
    tree = ast.parse((root / "Mass.py").read_text())
    functions = {node.name for node in tree.body if isinstance(node, ast.FunctionDef)}
    assert {"load_run_setup", "load_lsf", "build_template_norm_flux_lookup",
            "r_light_mass_to_light"} <= functions
    assert not any(isinstance(node, ast.ImportFrom) and node.module == "mass_light_normalization"
                   for node in tree.body)
    assert not (root / "mass_light_normalization.py").exists()


def write_template(folder, name="Ech1.30Zm0.40T01.0000_iTp0.00_baseFe.fits", offset=0):
    wave = np.arange(4400., 7501.)
    flux = (wave / 5000.) ** 3 + .4 * np.exp(-((wave - 4740) / 20) ** 2)
    fits.writeto(folder / name, flux, fits.Header({"CRVAL1": 4400.+offset, "CDELT1": 1.}))
    return wave, flux


def test_prepared_mean_replaces_original_window(tmp_path):
    wave, flux = write_template(tmp_path)
    selected = (wave > 4650) & (wave < 7150)
    prepared = log_rebin(
        [wave[selected][0], wave[selected][-1]],
        gaussian_filter1d(flux[selected], np.ones(selected.sum())), velscale=20.5,
    )[0]
    result, metadata = build_template_norm_flux_lookup(
        tmp_path, 4800., 7000., lsf_data=lambda x: np.full_like(x, np.sqrt(2.51**2+2.355**2)),
        lsf_templates=lambda x: np.full_like(x, 2.51), template_velscale=20.5,
    )
    assert result[(1., -.4, 0.)] == pytest.approx(prepared.mean(), rel=1e-13)
    assert not np.isclose(prepared.mean(), flux[(wave >= 4800)&(wave <= 7000)].mean())
    assert metadata["prepared_range"] == [4651., 7149.]


def test_duplicate_ssp_id_is_rejected(tmp_path):
    write_template(tmp_path)
    write_template(tmp_path, "Mch1.30Zm0.40T01.0000_iTp0.00_baseFe.fits")
    with pytest.raises(ValueError, match="Duplicate"):
        build_template_norm_flux_lookup(tmp_path, 4800., 7000.,
            lsf_data=lambda x: np.full_like(x, 3.), lsf_templates=lambda x: np.full_like(x, 2.51),
            template_velscale=20.5)


def test_nonfinite_prepared_input_is_rejected(tmp_path):
    write_template(tmp_path)
    p = next(tmp_path.glob("*.fits"))
    with fits.open(p, mode="update") as h:
        h[0].data[700] = np.nan
    with pytest.raises(ValueError, match="Non-finite"):
        build_template_norm_flux_lookup(tmp_path, 4800., 7000.,
            lsf_data=lambda x: np.full_like(x, 3.), lsf_templates=lambda x: np.full_like(x, 2.51),
            template_velscale=20.5)


def test_r_fractions_and_rescaling_invariance():
    w = np.array([[.5, .5], [0., 1.]])
    ml, fractions = r_light_mass_to_light(w, np.array([1., 4.]), np.array([1., 6.]), np.array([1., 3.]))
    np.testing.assert_allclose(fractions[0], [1/3, 2/3])
    np.testing.assert_allclose(fractions.sum(axis=1), 1.)
    np.testing.assert_allclose(ml, [3., 4.])
    scaled, _ = r_light_mass_to_light(w*np.array([[7.], [2.]]), np.array([1., 4.]), np.array([1., 6.]), np.array([1., 3.]))
    np.testing.assert_allclose(scaled, ml)
    with pytest.raises(ValueError, match="nonnegative"):
        r_light_mass_to_light(np.array([[1., -.1]]), np.ones(2), np.ones(2), np.ones(2))


def test_light_r_only_agrees_when_band_conversion_does_not_change_the_mixture():
    weights = np.array([[.5, .5], [1., 0.]])
    ml = np.array([1., 4.])
    corrected, _ = r_light_mass_to_light(weights, ml, np.array([1., 6.]), np.array([1., 3.]))
    light_r = np.sum(weights * ml[None, :], axis=1)
    np.testing.assert_allclose(corrected, [3., 1.])
    np.testing.assert_allclose(light_r, [2.5, 1.])
    assert np.log10(light_r[0]/corrected[0]) < 0
    equal_conversion, _ = r_light_mass_to_light(weights, ml, np.array([2., 6.]), np.array([1., 3.]))
    np.testing.assert_allclose(equal_conversion, light_r)


def test_invalid_lsf_and_zero_weight_sum_are_rejected(tmp_path):
    write_template(tmp_path)
    with pytest.raises(ValueError, match="resolution"):
        build_template_norm_flux_lookup(tmp_path,4800.,7000.,
            lsf_data=lambda x: np.full_like(x,2.), lsf_templates=lambda x:np.full_like(x,2.51),
            template_velscale=20.5)
    with pytest.raises(ValueError, match="contribution sums"):
        r_light_mass_to_light(np.zeros((1,2)),np.ones(2),np.ones(2),np.ones(2))


def make_setup(tmp_path):
    config = {"GENERAL": {"RUN_ID": "TEST", "REDSHIFT": .005566, "LSF_DATA": "data_lsf", "CONFIG_DIR": "."},
              "PREPARE_SPECTRA": {"VELSCALE": 41.},
              "SFH": {"LMIN": 4800., "LMAX": 7000., "NORM_TEMP": "LIGHT", "TEMPLATE_SET": "miles", "LIBRARY": "EMILES_baseFe_safe/", "LSF_TEMP": "template_lsf"}}
    (tmp_path / "CONFIG").write_text(yaml.safe_dump(config))
    np.savetxt(tmp_path/"data_lsf", [[4000, 3.], [10000, 2.54]])
    np.savetxt(tmp_path/"template_lsf", [[4000, 2.51], [9000, 2.51]])
    return fits.Header(config["SFH"])


def test_config_crosscheck_and_config_velscale(tmp_path):
    header = make_setup(tmp_path)
    setup = load_run_setup(tmp_path/"TEST_sfh_weights.fits", header, tmp_path)
    assert setup["template_velscale"] == 20.5
    assert setup["velscale_source"] == "CONFIG.PREPARE_SPECTRA.VELSCALE"
    header["LMAX"] = 8900.
    with pytest.raises(ValueError, match="LMAX"):
        load_run_setup(tmp_path/"TEST_sfh_weights.fits", header, tmp_path)


def test_config_from_another_galaxy_is_rejected(tmp_path):
    header = make_setup(tmp_path)
    with pytest.raises(ValueError, match="RUN_ID"):
        load_run_setup(tmp_path/"OTHER_sfh_weights.fits", header, tmp_path)


def test_missing_config_cannot_be_replaced_by_approximation(tmp_path):
    with pytest.raises(FileNotFoundError, match="CONFIG"):
        load_run_setup(tmp_path/"TEST_sfh_weights.fits", {}, tmp_path)


def test_missing_grid_identity_and_inconsistent_wavelength_grid(tmp_path):
    write_template(tmp_path)
    kwargs = dict(lsf_data=lambda x: np.full_like(x, 3.),
                  lsf_templates=lambda x: np.full_like(x, 2.51), template_velscale=20.5)
    with pytest.raises(ValueError, match="Missing"):
        build_template_norm_flux_lookup(tmp_path,4800.,7000.,grid_keys=[(99.,-.4,0.)],**kwargs)
    write_template(tmp_path,"Ech1.30Zm0.40T02.0000_iTp0.00_baseFe.fits",offset=.1)
    with pytest.raises(ValueError, match="grid differs"):
        build_template_norm_flux_lookup(tmp_path,4800.,7000.,**kwargs)


def test_config_velscale_does_not_open_intermediate_hdf5(tmp_path):
    header = make_setup(tmp_path)
    # Final CANFAR products omit intermediates. Even an unrelated/unreadable
    # leftover must not become a post-processing input or override CONFIG.
    (tmp_path/"TEST_bin_spectra.hdf5").write_bytes(b"not an HDF5 file")
    setup = load_run_setup(tmp_path/"TEST_sfh_weights.fits", header, tmp_path)
    assert setup["template_velscale"] == 20.5
    assert setup["velscale_source"] == "CONFIG.PREPARE_SPECTRA.VELSCALE"


def test_lsf_preserves_floor_and_extrapolation(tmp_path):
    make_setup(tmp_path)
    np.savetxt(tmp_path/"data_lsf", [[4000, 3.], [10000, 2.52]])
    z = .005566
    data, template = load_lsf(tmp_path/"data_lsf", tmp_path/"template_lsf", z)
    np.testing.assert_allclose(template(np.array([7000., 9100., 9250.])), 2.51)
    # The historical double-redshift floor check is deliberately preserved.
    np.testing.assert_allclose(data(np.array([10000./(1+z)])), 2.54*(1+z))
