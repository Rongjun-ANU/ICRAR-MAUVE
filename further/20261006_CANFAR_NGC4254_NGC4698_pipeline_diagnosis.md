# CANFAR pipeline diagnosis: NGC4254 and NGC4698

Date: 2026-10-06 (Australia/Perth).

## Findings

The two galaxies failed at different stages. These statuses come from the
downloaded CANFAR logs, not from the presence of output filenames.

| Galaxy | Mass | SFR+Z | proxy EWHa |
| --- | --- | --- | --- |
| NGC4254 | Success, 7m35s | Failed, 158m28s: pyqz integer overflow | Success, 35m22s |
| NGC4698 | Success, 6m54s | Success, 17m27s | Failed, 1m27s: invalid continuum FITS |

Inputs inspected:

- `/Volumes/Untitled/v3tk_v7.6.8_7000/product/NGC4254/`
- `/Volumes/Untitled/v3tk_v7.6.8_7000/product/NGC4698/`
- Corresponding logs under `mass_logs`, `sfr_logs`, and `proxy_ewha_logs`.
- Local `Mass.py`, `Mass.sh`, `SFR+Z.py`, `SFR.sh`, `proxy_EWHa.py`,
  `proxy_EWHa.sh`, and the model-grid helper modules and tests.

## NGC4254: pyqz QC flag storage

The SFR log reports all 40,755 eligible spatial-bin pyqz calculations complete,
then fails while converting result records into arrays:

```text
model_grid_diagnostics.py: _records_to_arrays
OverflowError: Python integer 91234 out of bounds for int16
```

The installed pyqz source constructs its raw flag by concatenating decimal
QC digits. A five-digit flag such as `91234` is therefore possible. It cannot
fit in signed 16-bit storage, whose maximum is 32767. This is a storage bug;
the repair does not change the model grid, fitting parameters, fluxes, masks,
or treatment of non-finite model results.

Four conversions could lose this value: record conversion, bin broadcasting,
HII/SF selection, and FITS serialization. The FITS schema validator also
required 16-bit storage. All these boundaries were addressed:

- `model_grid_diagnostics.py`: pyqz `flag` records use `np.int32`; bin
  broadcasting preserves wider supplied integer types.
- `SFR+Z.py`: region selection preserves the integer type; `PYQZ_FLAG_HII`
  and `PYQZ_FLAG_SF` are written and validated as `BITPIX=32`.
- `README.md`: documents the changed flag storage contract.
- `tests/test_pyqz_large_flags.py`: four targeted regression tests.

The other model QC arrays keep their existing 16-bit production types.
Changing only the initial conversion would be insufficient: the old
broadcast/masking code silently changed `91234` to `25698`.

Deploy BOTH corrected `SFR+Z.py` and `model_grid_diagnostics.py` to the existing
CANFAR script folder. Then, from the existing ICRAR science environment:

```bash
cd /arc/home/RongjunHuang/ICRAR/further
bash ./SFR.sh 7000 NGC4254
```

This production command overwrites the corresponding SFR+Z output and log.
Preserve prior versions if needed before running it. Mass does not need to be
repeated to resolve this failure. EWHa already succeeded and uses observed
Halpha flux and velocity, which this flag-storage fix does not alter.

The full NGC4254 SFR/model-grid pipeline was NOT rerun locally. The fix is
verified at the failing storage boundary and subsequent FITS boundaries.
Production completion still needs to be confirmed from the new CANFAR log.

## NGC4698: damaged continuum file

Both CANFAR and the local reproduction fail when Astropy opens:

```text
NGC4698_cont_cube.fits
OSError: No SIMPLE card found, this file does not appear to be a valid FITS file.
```

The downloaded file is 5,106,657,600 bytes, but its first 5,760 bytes are all
zero. Additional sampled 64-KiB blocks at offsets 100,000, 1,000,000 and
10,000,000 are also all zero. Samples farther into the file contain nonzero
data. This is a damaged file, not just a missing FITS keyword.

`ignore_missing_simple=True` cannot recover the missing header and data. The
Mass and SFR scripts do not read this continuum cube, which explains why they
can finish while EWHa fails. Their inputs and the gas/binning map spatial
dimensions were readable and compatible in the local checks.

Required replacement: an intact `NGC4698_cont_cube.fits` for the SAME 7000-A
run, with spatial shape `(1189, 608)`. The downloaded CONFIG identifies the
original Setonix output directory as:

```text
/scratch/pawsey1308/mauve/products/v3tk_v7.6.8_7000/NGC4698
```

Check that original copy before downloading it. The same error was already
present on CANFAR, so downloading the unchanged CANFAR file again may simply
copy the same damage. If the original is intact, replace the damaged CANFAR
and local continuum copies with it. Compare checksums across transfers. If
the original is also damaged, regenerate the matching nGIST continuum product.

The evidence does not establish whether the damage happened during nGIST
writing, storage, or a later transfer. No header reconstruction or source-file
modification was attempted.

After replacing and validating the continuum file, rerun only EWHa:

```bash
cd /arc/home/RongjunHuang/ICRAR/further
bash ./proxy_EWHa.sh 7000 NGC4698
```

## Old outputs and additional log failures

NGC4698 already has a readable EWHa output, but its preserved modification
time is 2026-07-18, while its current SFR+Z gas product is dated 2026-09-13.
The failed EWHa run writes no new output; the older file remains. Its presence
does not validate the failed run.

NGC4254's gas further product is dated 2026-09-10, while its EWHa product is
dated 2026-09-13. A failed SFR run likewise leaves an earlier gas further
product in place. Both separate `*_SFR_maps_further.fits` files are dated
2026-07-27; the current local `SFR+Z.py` writes its SFR layers into
`*_gas_bin_maps_further.fits` instead. Preserved file times are supporting
evidence; the logged exit status determines whether a particular run completed.

A read-only scan of all supplied logs found:

| Stage | Logs | Successes | Failures |
| --- | --- | --- | --- |
| Mass | 35 | 35 | 0 |
| SFR | 35 | 34 | 1: NGC4254 |
| EWHa | 35 | 28 | 7 |

The seven EWHa failures are NGC4351, NGC4402, NGC4405, NGC4450, NGC4457,
NGC4535, and NGC4698. All seven logs report the same `No SIMPLE card` error
while opening the continuum file. Only NGC4698's damaged bytes were inspected;
the other six continuum files were not supplied here.

## Executed validation

Python: `/opt/miniconda3/envs/ICRAR/bin/python` (3.13.3).

1. Before the patch, all four new regression tests failed for the expected
   overflow, truncation, or schema mismatch.
2. After the patch:

   ```bash
   /opt/miniconda3/envs/ICRAR/bin/python -m pytest -q \
     tests/test_pyqz_large_flags.py \
     tests/test_model_grid_diagnostics.py \
     tests/test_model_grid_compat.py
   ```

   Result: **50 passed in 49.25 s**. The new tests pass a controlled raw pyqz
   result through the real adapter, broadcast and masking functions, and the
   actual FITS writer call. They verify a FITS round trip of `91234`, the
   `-99` sentinel, 32-bit flag storage, and preservation of a spatial WCS card.
3. AST syntax checks passed for the changed Python code and new tests.
4. `git diff --check` passed for the tracked changes; the exact patch was
   reviewed. Repository root: `/Users/Igniz/Desktop/ICRAR`.
5. Both galaxies were run through the unchanged `proxy_EWHa.py` with temporary
   output destinations. The local `new_redshifts` file was absent, so a
   temporary two-row table used the exact redshifts printed in the CANFAR logs
   (`0.008026` and `0.003366`). No substitute redshift was inferred.

   ```bash
   /opt/miniconda3/envs/ICRAR/bin/python -u proxy_EWHa.py \
     -g NGC4254 \
     --root /Volumes/Untitled/v3tk_v7.6.8_7000 \
     --product-subdir product \
     --fallback-root /Users/Igniz/Desktop/ICRAR/further \
     --redshift-file /private/tmp/20261006_CANFAR_logged_redshifts \
     --out /private/tmp/20261006_NGC4254_proxy_EW_reproduction.fits
   ```

   NGC4254: exit 0; 1,045,283 valid proxy-EWHa pixels. `CONT_HA_MEAN` and
   `PROXY_EWHA` are exactly equal to the downloaded output. The legacy
   `PSEUDO_EWHA` agrees within floating-point roundoff (maximum relative
   discrepancy about 3.7e-15; checked with `rtol=1e-13`, zero absolute tolerance).

   The equivalent NGC4698 command: exit 1; reproduced the CANFAR `No SIMPLE
   card` failure at continuum opening. No replacement EWHa output was produced.

No full Mass or SFR pipeline was executed; no downloaded FITS file or CANFAR
file was overwritten. No memory or external literature was needed.
