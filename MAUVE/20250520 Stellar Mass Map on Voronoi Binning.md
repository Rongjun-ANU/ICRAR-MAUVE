# 20250520 Stellar Mass Map on Voronoi Binning

## 1. Load data

I have two files. One is the voronoi binng map:

```python
Filename: IC3392_SPATIAL_BINNING_maps.fits
No.    Name      Ver    Type      Cards   Dimensions   Format
  0  PRIMARY       1 PrimaryHDU       4   ()      
  1  BINID         1 ImageHDU        26   (437, 438)   float64   
  2  FLUX          1 ImageHDU        26   (437, 438)   float64   
  3  SNR           1 ImageHDU        26   (437, 438)   float64   
  4  SNRBIN        1 ImageHDU        26   (437, 438)   float64   
  5  XBIN          1 ImageHDU        26   (437, 438)   float64   
  6  YBIN          1 ImageHDU        26   (437, 438)   float64   
None
```

The other is about the SFH, including `GIRD` for ages and metallicities and `WEIGHTS` for each bin:

```python
Filename: IC3392_sfh-weights.fits
No.    Name      Ver    Type      Cards   Dimensions   Format
  0  PRIMARY       1 PrimaryHDU      23   ()      
  1  WEIGHTS       1 BinTableHDU     27   4077R x 1C   [477D]   
  2  GRID          1 BinTableHDU     31   477R x 3C   [D, D, D]   
None
```

The data size makes sense for me. The IC3392 datacube is stacking into a $437\times438$  map with 4077 voronoi bins. By further examination, I know that each bin has corresponding weights on a 9 (rows of metallicities) $\times$ 53 (columns of ages) SPS grid. The header of `WEIGHTS` also tells me that wavelength range is $4800-7000\AA$ and the SPS template is `MILES`. I also check that in each bin, the `WEIGHTS` are summed up to be 1. 

## 2. Getting the Mass-to-Light ratio

I go the the [`MILES` website](https://research.iac.es/proyecto/miles/pages/predicted-masses-and-photometric-observables-based-on-photometric-libraries.php) and find that they provide tables for model predictions: 

![image-20250520160120079](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250520 Stellar Mass Map on Voronoi Binning.assets/image-20250520160120079.png)

By looking at the [Fraser-McKelvie et al. (2024)](https://arxiv.org/abs/2411.03430) and checking the exact value in `GRID` data, I think **BaSTI+Chabrier** is the correct one that I should use, so I download the header and table into `BaSTI+Chabrier.dat` and extract the R-band Mass-to-Light ratio ($M/L_R$). 

Now for each bin, I can match the $9\times53$ METAL-LOGAGE `GRID` with  **BaSTI+Chabrier** table to get their $M/L_R$. Then I dot product the `WEIGHTS` with $M/L_R$ to get effective $M/L_R$ at each bin. And therefore I get $M/L_R$ map and stored it as a new `ImageHDU` called `ML_R` in `IC3392_SPATIAL_BINNING_maps_extended.fits`. 

Here I perform a check for a comparison with . Since I also have access to `FLUX` data of each bin, I can compute the flux-weighted mean R-band Mass-to-Light ratio:
$$
\overline{M/L_R} = \frac{\Sigma_i^{4077} (M/L_R)_i\cdot F_i}{\Sigma_i^{4077}F_i}.
$$
And this gives $\overline{M/L_R} = 1.442 = \log(0.159)$. Recall that the one I get from legacy data with `BC03` template is $1.40 = \log(0.15)$ and the other one I get from `pPXF` with `E-MILES` template is $2.105 = \log(0.323)$. Now here is a question, how can `MILES` one become so close to `BC03` one rather than `E-MILES`.

## 3. Stellar mass map

To get stellar mass map, I make such assumptions: 

1. The distance to IC3392 is 11.5 Mpc
2. $4800-7000\AA$ is SDSS r-band coverage with effective wavelegth at $6231\AA$
3. Solar R-band Magnitude is 4.64

Then I construce the stellar mass map and also save it into an extra array `LOGMSTAR` of `IC3392_SPATIAL_BINNING_maps_extended.fits`.  

In this way, the total r band magnitude is 12.323 and total stellar mass is $\log(9.207)$.  As a comparison, these values are 11.958 & $\log(9.32)$ in legacy data and 12.369 & $\log(9.355)$ from previous `pPXF`. 

Below is the maps I store in `IC3392_SPATIAL_BINNING_maps_extended.fits`. Top left is the stellar mass map, bottom left is the $M/L_R$ map, and the one on the RHS is just the flux. 

![IC3392_SPATIAL_BINNING_maps_extended.fits.HDU_8_LOGMSTAR-IC3392_SPATIAL_BINNING_maps_extended.fits.HDU_2_FLUX-IC3392_SPATIAL_BINNING_maps_extended.fits.HDU_7_ML_R-image-2025-05-20-15-33-02](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250520 Stellar Mass Map on Voronoi Binning.assets/IC3392_SPATIAL_BINNING_maps_extended.fits.HDU_8_LOGMSTAR-IC3392_SPATIAL_BINNING_maps_extended.fits.HDU_2_FLUX-IC3392_SPATIAL_BINNING_maps_extended.fits.HDU_7_ML_R-image-2025-05-20-15-33-02.png)

