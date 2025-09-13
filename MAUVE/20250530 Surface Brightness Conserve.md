# 20250530 Surface Brightness Conserve

I nail it. 

## Misleading 1

At first by eyes, it seems that even from the raw data (i.e., without projection, binning, and galactic extinction), the flux is descending faster in MAUVE than Legacy data. Particularly, the SDSS r-band magnitudes are 11.72 and 12.31 for Legacy data after projection and MAUVE data, respectively. And I have checked that the images from Legacy look similar among raw, projection and binning. Also I try different reproject function (i.e., `_exact`, `_interp` and `adaptive`) and download Legacy data at different pixel scales, but the results all look the same. So that is why I say I think Legacy stuff is fine, but there could be some issues with getting the r-band flux/magnitude from MAUVE data on Wednesday meeting.  Here I use `get_ab_magnitude` or `get_ab_maggies` from `speclite.filters`. 

However, is it truely the issue in MAUVE data? I also find something supscious. The data from [NED](https://ned.ipac.caltech.edu/byname?objname=ic3392&hconst=67.8&omegam=0.308&omegav=0.692&wmap=4&corr_z=1) suggests that the IC3392 magnitude in SDSS r-band is 12.3.  Also my previous calculation (i.e. isophotal+contour) and and the data from Cotese et al. 2011 are both 11.96. Since the area from MAUVE seems not completely cover the entire area of IC3392, I think there is no reason that after projection (shrinking) into MAUVE area, the magnitude from Legacy data is even lower than 12. 

A further investigation confirms my idea. Yesterday I try using another package call `mpdaf` to see what r-band magnitude I can get for MAUVE data in this parallel way. And it turns out it is still $\sim12.3$, so I believe must be something wrong with Legacy side, because this is an independent check. 

## Misleading 2

Ok now I turn into examine my procedure with Legacy data. Since I have found that all the flux maps look the same in Legacy side, I feel like there must be something wrong at the very beginning, i.e. the projection. 

As a test, I try drawing an fix 20 arcsec aperture for both Legacy and MAUVE raw data and see the difference of the total flux inside the aperture. 

![image-20250530124108966](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250530124108966.png)

`Total r-band flux inside ±20″ aperture: 8035.229 nanomaggies`

![image-20250530124403073](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250530124403073.png)

`Total r-band flux inside ±20″ aperture: 7968.000 nanomaggies`

It turns out the flux is almost identical between them, even though the color graident looks different. I will say that is probably due to the different PSFs and pixel scales. 

The key is here. When I look at the value that just after projection for Legacy data, it goes to $\sim 13000$, so I need to carefully check if the function I use to do projection is truely flux-conserved. 

Luckily, I notice that the ratio of 13000 and 8000 is $\sim1.7$, which is coincidently the same as the scaling factor due to pixel scales (*Attention is all you need!*). Recall that the pixel scales of  Legacy and MAUVE data are 0.262 and 0.2 arcsec/pixel, respectively, so the scaling factor is $(0.262/0.2)^2=1.72$. Hence, I suspect that the flux does not resacle properly during projection. 

The **most straightfoward way** to fix this is just convert the flux from something per pixel square to surface brightness unit (i.e., something per arcsec square), which properly accounts for different pixel scales. By doing this, the results satify me (r-band magnitude now is 12.3). Below is the Legacy data after projection now.

![image-20250530135158149](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250530135158149.png)

`Total r-band flux inside ±20″ aperture: 7966.787 nanomaggies`

Now the flux or surface brightness is really conserved (check with first figure), so the total flux variation is within $1\%$. 

And it turns out on [Github](https://github.com/astropy/reproject) people are still debating if [reproject_exact() is not conserving flux](https://github.com/astropy/reproject/issues/199). 

Actually, the **simplest** way to do is replacing `reproject_exact` with `reproject_adaptive(conserve_flux=True)`, **even though they claim that `reproject_exact` is flux-conserving!** I guess it is conserved in flux per pixel square, but not conserved in flux per arcsec square (that is why i say it should be surface brightness conserve). 

## Missleading 3

As mentioned, I do try using Legacy data at different pixel scales and it does not work. Based on these results, I am doubting if Legacy website truely provide the data at different pixel scales or do they resacle the data properly. And indeed, the total flux in same aperture from Legacy data differs at different pixel scales. On Legacy website, they say 0.262 arcsec/pixel is their native resolution when running the reduction pipeline. My suggestion is that for safety, just **download and use the Legacy data at 0.262 arcsec/pixel and do not trust the data at different pixel scales!** 

## Updated results

![image-20250530135510968](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250530135510968.png)

Again, I show the 1) r-band magnitude map after projection+binning; 2) r-band magnitude map after projection, binning and galactic extinction correction; 3) r-band mass-to-light ratio map; and 4) stellar mass map. 