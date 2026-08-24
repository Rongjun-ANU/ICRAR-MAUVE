This repository contains the photoionization model grids used in Peng et al. (2026). The grids are based on the photoionization models of Ji & Yan (2020) and generated using CLOUDY C17.03.

photoionization_grid_native.fits
    Native CLOUDY model grid.

photoionization_grid_interpolated.fits
    40 × 40 interpolated grid in the (log Z/Z☉, log U) parameter space used in the Bayesian inference.

Each row corresponds to one photoionization model with:

- log(Z/Zsun)
- log(U)
- log(N/O)

Predicted emission-line intensities and commonly used line ratios are provided.

The CLOUDY input files used to generate the models are also included.

These grids were used throughout the metallicity inference presented in the paper.