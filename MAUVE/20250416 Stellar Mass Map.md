# 20250416 Stellar Mass Map

## IC3392

coordiate query by NOIRlab https://astroarchive.noirlab.edu/portal/search/

https://www.legacysurvey.org/viewer/jpeg-cutout?ra=187.18031864729997&dec=14.99943425062&layer=ls-dr10&pixscale=0.262&size=512&bands=g

![image-20250416235624186](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250416 Stellar Mass Map.assets/image-20250416235624186.png)

wget "https://www.legacysurvey.org/viewer/fits-cutout?ra=187.18031864729997&dec=14.99943425062&layer=ls-dr10&pixscale=0.262&bands=gi&size=512" -O IC3392_DESI_gi.fits

I extract `.fits` file and check the g band image first:

![image-20250416235924749](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250416 Stellar Mass Map.assets/image-20250416235924749.png)

According to the DESI DR10 document, https://www.legacysurvey.org/dr10/description/, the brightnesses of objects are all stored as linear fluxes in units of nanomaggies. The conversion from linear fluxes to magnitudes is $m=22.5−2.5\log_{10}(flux)$. 

Then i convert nanomaggies to AB magnitude:

![image-20250417000110464](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250416 Stellar Mass Map.assets/image-20250417000110464.png)

Similarly, I can do for i band:

![image-20250417000156761](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250416 Stellar Mass Map.assets/image-20250417000156761.png)

The blank rectangle area is `nan` data (negative flux). 

In Taylor+2011 (https://academic.oup.com/mnras/article/418/3/1587/1060932#91710255), they purpose an empirical relation between (*g*−*i*) colour, *i*-band luminosity and stellar mass as:

log(M_*/M_sun) = 1.15+0.7(g-i)-0.4M_i, where M_i is the absolute magnitude in the restframe i-band, expressed in the AB system. 

Here is the g-i color map:

![image-20250417112212265](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250416 Stellar Mass Map.assets/image-20250417112212265.png)

The distance I use to calculate absolute i band magnitude is 13.63Mpc from Soria+2021: https://academic.oup.com/mnras/article/512/3/3284/6517474. It is a median redshift-independent distance from NED, and standard deviation of the published values.

Now i can derive the stellar mass map, and then sum up to get total stellar mass in solar masses: 4.41e+11, or log10(11.64).

![image-20250417112223655](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250416 Stellar Mass Map.assets/image-20250417112223655.png)

Alternatively, I estimate the stellar mass by intergrating g and i band flux first. In this case, I have total stellar mass (integrated flux): 3.10e+09, or log10(9.49) solar masses, which is quite different from the previous method. 

My understanding is that, 

1. mathematically, the convertion formula is non-linear, which means $\sum_{i}M(Flux_i) \neq M(\sum_{i}Flux_i)$;
2. this empirical relation is calibrated by the total magnitude from the galaxy
3. when we calcualte the flux from a single pixel, we actually include the flux emitted from adjacent pixels, and so integrating pixel-by-pixel leads to overestimated total magnitude, and thus much higher stellar mass in this case. It is a binning effect. 

In fact, this should shows a convergence if i choose large bin sizes. Here i choose the bin sizes in $1\times1$ (single pixel per bin), $4\times4$, $16\times16$, $64\times64$, $256\times256$, and $512\times512$ (total galaxy):

![image-20250417080415552](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250416 Stellar Mass Map.assets/image-20250417080415552.png)

And it indeed converges to log(9.49) with increading bin size. As a comparison, the stellar mass in Soria+2021 is log(9.75 ± 0.09), which is derived from WISE W1 and W2 fluxes, as calibrated by Cluver et al. (2014) (see also Jarrett et al. 2019). The median NED distances were adopted.  

## NGC4383

![image-20250417080835266](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250416 Stellar Mass Map.assets/image-20250417080835266.png)

![image-20250417112254155](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250416 Stellar Mass Map.assets/image-20250417112254155.png)

In Adam's paper (https://academic.oup.com/mnras/article/530/2/1968/7642869),
he assume that NGC 4383 is at the distance of Virgo, 16.5 Mpc (Mei et al. 2007)
and NGC 4383 is an intermediate-mass galaxy (log(M⋆/M⊙) = 9.44; Leroy et al. 2019). 

I use the same distance and get log10(9.41) solar masses in total. Below is stellar mass map. 

![image-20250417112308265](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250416 Stellar Mass Map.assets/image-20250417112308265.png)

![image-20250417112317487](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250416 Stellar Mass Map.assets/image-20250417112317487.png)

## NGC4501

![image-20250417081525153](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250416 Stellar Mass Map.assets/image-20250417081525153.png)

![image-20250417112417740](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250416 Stellar Mass Map.assets/image-20250417112417740.png)

In Soria+2021, they have 17.10 ± 4.58Mpc and log(11.03 ± 0.09) solar masses, while I get log(10.75). 

![image-20250417112428956](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250416 Stellar Mass Map.assets/image-20250417112428956.png)

![image-20250417112440076](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250416 Stellar Mass Map.assets/image-20250417112440076.png)

## NGC4192

![image-20250417081721711](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250416 Stellar Mass Map.assets/image-20250417081721711.png)

![image-20250417112454471](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250416 Stellar Mass Map.assets/image-20250417112454471.png)

In Soria+2021, they have 13.55 ± 2.93Myr and log(10.63 ± 0.09) solar masses, while I get log(10.25). 

![image-20250417112506104](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250416 Stellar Mass Map.assets/image-20250417112506104.png)

![image-20250417112515643](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250416 Stellar Mass Map.assets/image-20250417112515643.png)