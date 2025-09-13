## IC3392 HRS172

coordiate query by NOIRlab https://astroarchive.noirlab.edu/portal/search/

https://www.legacysurvey.org/viewer/jpeg-cutout?ra=187.18031864729997&dec=14.99943425062&layer=ls-dr10&pixscale=0.262&size=512&bands=g

wget "https://www.legacysurvey.org/viewer/jpeg-cutout?ra=187.18031864729997&dec=14.99943425062&layer=ls-dr10&pixscale=0.262&size=512&bands=g" -O IC3392_g.jpeg

https://www.legacysurvey.org/viewer/fits-cutout?ra=187.18031864729997&dec=14.99943425062&layer=ls-dr10&pixscale=0.262&bands=gi&width=512&height=512

wget "https://www.legacysurvey.org/viewer/fits-cutout?ra=187.18031864729997&dec=14.99943425062&layer=ls-dr10&pixscale=0.262&bands=gi&size=3000" -O IC3392_DESI_gi.fits

wget "https://www.legacysurvey.org/viewer/fits-cutout?ra=187.18031864729997&dec=14.99943425062&layer=ls-dr10-model&pixscale=0.262&bands=gi&size=512" -O IC3392_DESI_gi_model.fits

wget "https://www.legacysurvey.org/viewer/fits-cutout?ra=187.18031864729997&dec=14.99943425062&layer=ls-dr10&pixscale=0.262&bands=gri&size=3000" -O IC3392_DESI_gri.fits

wget "https://www.legacysurvey.org/viewer/fits-cutout?ra=187.18031864729997&dec=14.99943425062&layer=ls-dr10&pixscale=0.262&bands=grih&size=3000" -O IC3392_DESI_grih.fits

https://www.legacysurvey.org/dr10/description/
The brightnesses of objects are all stored as linear fluxes in units of nanomaggies. The conversion from linear fluxes to magnitudes is $m=22.5−2.5\log_{10}(flux)$. These linear fluxes are well-defined even at the faint end, and the errors on the linear fluxes should be very close to a normal distribution. The fluxes can be negative for faint objects, and indeed we expect many such cases for the faintest objects.

log(M_*/M_sun) = 1.15+0.7(g-i)-0.4M_i, where M_i is the absolute magnitude in the restframe i-band, expressed in the AB system. Here we assume the distance to the IC3392 is 11.1Mpc from Tully+2016.

13.65 ± 3.86Myr, 9.75 ± 0.09 in Soria+2021, I got 9.49
Note from Soria+2021: Median redshift-independent distance from NED, and standard deviation of the published values. When a galaxy has Cepheid distance measurements, we used the median and standard deviation of those values only. Stellar mass, derived from WISE W1 and W2 fluxes, as calibrated by Cluver et al. (2014) (see also Jarrett et al. 2019). The median NED distances were adopted. 

HRS172: 12.739, 11.660, 9.77 

Mass to light ratio in g-band: log(0.23) or 1.72
Mass to light ratio in r-band: log(0.13) or 1.35
Mass to light ratio in i-band: log(0.06) or 1.14

Total g-band luminosity: 1.31e+09 L_sun = log(L/L_sun) = 9.12
Total r-band luminosity: 1.66e+09 L_sun = log(L/L_sun) = 9.22
Total i-band luminosity: 1.97e+09 L_sun = log(L/L_sun) = 9.30

Stellar mass sum up in IC3392: 2.24e+09 M_sun or log (M_*/M_sun) = 9.35

## Doing ppxf stellar fitting
ppxf_stellar_fitting.ipynb

## Parallel check with Legacy data 
reproject.ipynb 
DESI_reproject.ipynb
map_MILES.ipynb
Legacy_vs_MUSE.ipynb
aperture_check.ipynb


## NGC4383 HRS142

wget "https://www.legacysurvey.org/viewer/fits-cutout?ra=186.3561&dec=16.4699&layer=ls-dr10&pixscale=0.262&bands=gi&size=3000" -O NGC4383_DESI_gi.fits

20.1Mpc from Theureau+2007
or in Adam's paper: https://academic.oup.com/mnras/article/530/2/1968/7642869
Throughout this work, we assume that NGC 4383 is at the distance of Virgo, 16.5 Mpc (Mei et al. 2007)
NGC 4383 is an intermediate-mass galaxy (log(M⋆/M⊙) = 9.44; Leroy et al. 2019)

HRS142: 12.392, 11.654, 9.42

## NGC4501 HRS190

wget "https://www.legacysurvey.org/viewer/fits-cutout?ra=187.99646878&dec=14.42031916&layer=ls-dr10&pixscale=0.262&bands=gi&size=3000" -O NGC4501_DESI_gi.fits

wget "https://www.legacysurvey.org/viewer/fits-cutout?ra=187.99646878&dec=14.42031916&layer=ls-dr10&pixscale=0.262&bands=grih&size=3000" -O NGC4501_DESI_grih.fits

15.3Mpc from Mandel+2011

17.10 ± 4.58Mpc, 11.03 ± 0.09 in Soria+2021, I got 10.75

HRS142: 10.011, 8.832, 10.98

## NGC4192 HRS91

wget "https://www.legacysurvey.org/viewer/fits-cutout?ra=183.45121278623998&dec=14.900542632&layer=ls-dr10&pixscale=0.262&bands=gi&size=3000" -O NGC4192_DESI_gi.fits

17.6Mpc from Mandel+2011

13.55 ± 2.93Myr, 10.63 ± 0.09 in Soria+2021, I got 10.25

HRS91: 10.560, 9.461, 10.65

https://www.legacysurvey.org/viewer/ls-dr10/1/14/7963/6632.jpg
https://www.legacysurvey.org/viewer/ls-dr10-model/1/14/7963/6632.jpg
https://www.legacysurvey.org/viewer/ls-dr10-resid/1/14/7963/6632.jpg

## NGC4380

wget "https://www.legacysurvey.org/viewer/fits-cutout?ra=186.3424083&dec=10.0167056&layer=ls-dr10&pixscale=0.262&bands=gi&size=3000" -O NGC4380_DESI_gi.fits

## NGC4302

wget "https://www.legacysurvey.org/viewer/fits-cutout?ra=185.4269875&dec=14.5977611&layer=ls-dr10&pixscale=0.262&bands=gi&size=3000" -O NGC4302_DESI_gi.fits

## NGC4254

wget "https://www.legacysurvey.org/viewer/fits-cutout?ra=184.7067708&dec=14.4164889&layer=ls-dr10&pixscale=0.262&bands=gi&size=3000" -O NGC4254_DESI_gi.fits

## NGC4064

wget "https://www.legacysurvey.org/viewer/fits-cutout?ra=181.0466292&dec=18.44345&layer=ls-dr10&pixscale=0.262&bands=gi&size=3000" -O NGC4064_DESI_gi.fits

wget --http-user=RongjunHuang \
     --http-password=@Rj19990726 \
     "https://cloud.datacentral.org.au/remote.php/webdav/\
mauve/teamdata/cubes/red_v3/\
NGC4064_DATACUBE_FINAL_WCS_Pall_mad_red_v3.fits"


## NGC4536

wget "https://www.legacysurvey.org/viewer/fits-cutout?ra=188.6128&dec=2.1883&layer=ls-dr10&pixscale=0.262&bands=gi&size=3000" -O NGC4536_DESI_gi.fits

## NGC4535

wget "https://www.legacysurvey.org/viewer/fits-cutout?ra=188.58476813196&dec=8.19775235781&layer=ls-dr10&pixscale=0.262&bands=gi&size=3000" -O NGC4535_DESI_gi.fits

## NGC4698

wget "https://www.legacysurvey.org/viewer/fits-cutout?ra=192.09545117940996&dec=8.48740765341&layer=ls-dr10&pixscale=0.262&bands=gi&size=3000" -O NGC4698_DESI_gi.fits

## Leo

wget "https://www.legacysurvey.org/viewer/fits-cutout?ra=169.95&dec=13.27&layer=ls-dr10&pixscale=1&bands=griz&size=3000" -O Leo_DESI_griz.fits

## M97 Owl Nebura
wget "https://www.legacysurvey.org/viewer/fits-cutout?ra=168.69880123&dec=55.01902301&layer=ls-dr10&pixscale=0.262&bands=griz&size=1000" -O M97_DESI_griz.fits

## Cap
delta (Deneb Algedi): ra=326.76018433125&dec=-16.12728708528
alpha2 (Algedi): ra=304.5135649821425&dec=-12.54485212019527

wget -O Capricornus_W1W2_FULLfield.fits \
  "https://www.legacysurvey.org/viewer/fits-cutout?\
ra=318.5&dec=-18.5&layer=unwise-neo7&pixscale=45&size=3000&bands=W1W2"

wget "https://www.legacysurvey.org/viewer/fits-cutout?ra=318.5&dec=-18.5&layer=ls-dr10&pixscale=45&bands=griz&size=3000" -O Capricornus_DESI_griz.fits

