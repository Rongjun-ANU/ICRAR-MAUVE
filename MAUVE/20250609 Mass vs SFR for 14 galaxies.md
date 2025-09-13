## O'Donnell+1994 to Calzetti+2000

When calculating the gas $E(B-V)$, I previously adopt the same approach as [Belfiore et al. (2023)](https://ui.adsabs.harvard.edu/#abs/2023A%26A...670A..67B), i.e., $k_{H_\beta}=3.609$ and $k_{H_\alpha}=2.535$, with $R_V = 3.1$. However, in `nGIST`, we adopt Calzetti-like extinction curve, so now I change to $k_{H_\beta}=3.609$ and $k_{H_\alpha}=2.535$, with $R_V = 4.05$.

## Kroupa to Chabrier

I look at the Figure 4 of Madau & Dickson (2014) and it turns out that converting SFRs from Kroupa to Chabrier IMF is 0.63/0.67. `nGIST` adopts Chabrier IMF in SPS templates and SFR coefficient $C_{H\alpha} = 5.3\times10^{-42}$ from Calzetti et al. (2007) adopts Kroupa, so I times that constant to keep consistency. 

## Pipeline 

I have done converting my previous `.ipynb` to an automatic and universal `.py` to create stellar mass and SFR map. 

Stellar mass map and relevant properties are stored in `*_SPATIAL_BINNING_maps_extended.fits` (an extended version of original `_SPATIAL_BINNING_maps.fits`). 

SFR map and relevant properties are stored in `*_gas_BIN_maps_extended.fits` (an extended version of original `_gas_BIN_maps.fits`). 

All are located at `/home/RongjunHuang/ICRAR/extended` on CANFAR. 

## Spatially-resolved $\Sigma_*\,\text{vs.}\,\Sigma_{\mathrm{SFR}}$ in SF regions

Now I can create Spatially-resolved $\Sigma_*\,\text{vs.}\,\Sigma_{\mathrm{SFR}}$ in SF regions

First, for IC3392 (hide in the bottom) and NGC4501:

![image-20250609220114566](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250609220114566.png)

Then, for all 14 galaxies

![image-20250609214356829](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250609214356829.png)