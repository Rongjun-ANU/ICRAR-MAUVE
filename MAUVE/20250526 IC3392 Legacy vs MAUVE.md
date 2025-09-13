# 20250526 IC3392 Legacy vs MAUVE

From now on, all the distance set to 16.5 Mpc. 

## 1. Reprojected Legacy data

In order to have a fair comparison between Legacy and MAUVE data, I need to make sure they are in the same size, same pixel scale and same voronoi binning. Here I use `reproject_exact` for this transformation.

Then I use the same procedure as before to extract the corrected r-band magnitude, stellar mass and $M/L_r$ map:

1. Galactic extinction. E(B-V) = 0.0369 from Corteste et al. 2012 and A/E(B-V) = 2.165 from Legacy data description. 
2. Isophotal mask. Though it does not change anything here. 
3. Galacy contour. Though it does not change anything here.
4. r-band magnitude map. I convert nanomaggies to r-band magnitude.
5. Stellar mass map. I use g-i map and calibration from Taylor et al. 2011 to get stellar mass. 
6. $M/L_r$ map. With r-band luminosity and stellar mass map, I can calculate the $M/L$ in r-band. 

## 2. MAUVE data

The way to deal with MAUVE data is different:

1. $M/L_R$ map. As described last week, I use the predicted table provided on [MILES website](https://research.iac.es/proyecto/miles/pages/predicted-masses-and-photometric-observables-based-on-photometric-libraries.php) to extract the $M/L$ in R-band. **I need to point out that the R-band here is Cousins filters (SDSS filters are provided in E-MILES table).** 
2. Galactic extinction. E(B-V) map is provided by `IC3392_SFH_maps.fits`. According to Table 6 of [Schlafly & Finkbeiner 2011](https://iopscience.iop.org/article/10.1088/0004-637X/737/2/103), I use A/E(B-V) = 2.285 for SDSS r-band at $R_V=3.1$. 
3. r-band magnitude map. I use `speclite.filters` and specify `sdss2010-r` filter to get r-band magnitude of MAUVE data. 
4. Stellar mass map. With r-band luminosity and $M/L_R$ map, I can calculate the stellar mass. 

## 3. Comparion

Here I show the corrected r-band magnitude, $M/L_r$ and stellar mass map from Legacy and MAUVE data.

![image-20250526153617880](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250526153617880.png)

Then I further make a r-band M/L ratio map. 

![image-20250526154246280](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250526154246280.png)

or in log scale:

![image-20250526154945572](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250526154945572.png)

It seems that Mass-to-Light Color Relation (MLCR) tends to predict older population in star-forming regions (the inner elliptical structure) than `pPXF`, and vice versa for outer regions. 