# 20260517 ICRAR local Conda environment inventory

Generated in `/Users/Igniz/Desktop/ICRAR` on 2026-05-17.

This documents the local environment named `ICRAR`, currently located at:

```text
/opt/miniconda3/envs/ICRAR
```

## Checks performed

```bash
which conda
conda --version
which mamba
mamba --version
conda env list
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate ICRAR
python -V
which mamba
mamba --version
conda env export -n ICRAR --from-history
conda list -n ICRAR --json
```

## Conda/Mamba status

- Base Conda executable: `/opt/miniconda3/bin/conda`
- Conda version: `26.1.1`
- Base shell `mamba`: not found on `PATH`
- After `conda activate ICRAR`, Mamba is available at `/opt/miniconda3/envs/ICRAR/bin/mamba`
- Activated-env Mamba version: `2.1.1`
- Activated-env Python: `Python 3.13.3`
- Local platform: `macOS-26.4.1-arm64-arm-64bit-Mach-O`

Important interpretation note: Conda does not reliably preserve a perfect "the user typed this exact install command" audit trail for every dependency. The cleanest split is:

- `conda env export --from-history`: explicit Conda packages requested in this environment history.
- `conda list`: source split by channel. Packages whose channel is `pypi` were installed by Pip; packages whose channel is not `pypi` were installed/resolved by Conda/Mamba.

The current `conda list` source split is:

- Total packages: `404`
- Conda/Mamba-sourced packages: `388`
- Pip-sourced packages: `16`

The full package list is in:

```text
/Users/Igniz/Desktop/ICRAR/20260517_ICRAR_local_env_full_package_list.md
```

The machine-readable reference used for Setonix verification is:

```text
/Users/Igniz/Desktop/ICRAR/20260517_ICRAR_local_conda_list_reference.json
```

## Explicit Conda packages from environment history

These are the packages recorded by Conda as explicit requested specs:

- `astroquery`
- `ipykernel`
- `ipympl`
- `ipywidgets`
- `jupyter`
- `jupyter_client`
- `jupyter_core`
- `mamba`
- `matplotlib`
- `numpy`
- `opencv`
- `pandas`
- `pillow`
- `plotly`
- `reproject`
- `scikit-learn`
- `speclite`
- `spectral-cube`
- `synphot`
- `tqdm`

## Pip-sourced packages

These are the packages whose `conda list` channel is `pypi`:

- `asdf==4.1.0`
- `asdf-astropy==0.7.1`
- `asdf-coordinates-schemas==0.4.0`
- `asdf-standard==1.2.0`
- `asdf-transform-schemas==0.6.0`
- `asdf-wcs-schemas==0.4.0`
- `build==1.2.2.post1`
- `capfit==2.6.6`
- `gwcs==0.24.0`
- `montage-wrapper==0.9.9`
- `mpdaf==3.6`
- `ndcube==2.3.1`
- `ppxf==9.4.2`
- `pyproject-hooks==1.2.0`
- `semantic-version==2.10.0`
- `specutils==1.20.2`

## Setonix rebuild recipe produced from this audit

For Setonix, use the stricter Linux-target recipe:

```text
/Users/Igniz/Desktop/ICRAR/20260517_ICRAR_Rongjun_environment.yml
```

That file now pins all Python-facing Conda packages and all Pip packages with exact versions. It intentionally does not pin every low-level binary dependency, because the local env is macOS arm64 and Setonix is Linux. A literal full export contains macOS-only packages such as `libcxx`, `pyobjc-core`, and `pyobjc-framework-cocoa`, and a full exact macOS export cannot be installed on Setonix.

The strict recipe was dry-run solved for `linux-64` with Conda using Linux virtual-package overrides. A stricter candidate that pinned every Conda package failed because macOS/Setonix low-level runtime packages are not interchangeable.

After installing on Setonix, verify the result with:

```text
/Users/Igniz/Desktop/ICRAR/20260517_compare_ICRAR_env.py
```

## Commands to regenerate this audit

```bash
conda env export -n ICRAR --from-history
conda list -n ICRAR
conda list -n ICRAR --json
conda env export -n ICRAR --no-builds
```
