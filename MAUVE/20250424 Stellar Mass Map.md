# 20250424 Stellar Mass Map

## NGC4383 (HRS142)

I choose NGC4383 as a baseline because it is more complete in both g and i-bands. 

### Data Acquisition

First, I cutout data by 

```bash
wget "https://www.legacysurvey.org/viewer/fits-cutout?ra=186.3561&dec=16.4699&layer=ls-dr10&pixscale=0.262&bands=gi&size=3000" -O NGC4383_DESI_gi.fits
```

RA and Dec are acquired by `https://astroarchive.noirlab.edu/portal/search`. 

Pixel scale is set to 0.262''/pix, which is  (approximately) the native pixels used by the [Tractor](https://github.com/dstndstn/tractor). 

Resolution is fixed at $3000\times3000$ (the maximum i can get) to make sure large enough range to cover the full size of galaxies. 

Brackground subtraction is done by pipelines. The only correction I need to do is galactic extinction. 

![image-20250424111158020](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250424 Stellar Mass Map.assets/image-20250424111158020.png)

### Galactic Extinction

The $E(B-V)$ is provided by [Table 6](https://www.aanda.org/articles/aa/full_html/2012/08/aa19312-12/T3.html) of [Cortese et al. 2012](https://www.aanda.org/articles/aa/full_html/2012/08/aa19312-12/aa19312-12.html): 

```python
    "IC3392": 0.0369,
    "NGC4383": 0.0237, 
    "NGC4501": 0.0380, 
    "NGC4192": 0.0350, 
```

The [galactic extinction coefficients](https://www.legacysurvey.org/dr9/catalogs/#galactic-extinction-coefficients) provided by legacy survey are  $𝐴/𝐸(𝐵−𝑉) = 3.214, 1.592$ for the DECam g and i filters. 

Hence, we can get corrected flux by

$$flux_{correction} = flux\times10^{0.4\times𝐴/𝐸(𝐵−𝑉)\times(𝐵−𝑉)}$$

![image-20250424113259734](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250424 Stellar Mass Map.assets/image-20250424113259734.png)

Then we can convert flux to AB magnitude by

$$m = 22.5-2.5\times\log(flux_{correction})$$

![image-20250424114008055](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250424 Stellar Mass Map.assets/image-20250424114008055.png)

```python
g-band magnitude range: 44.52 to 18.87
i-band magnitude range: 44.96 to 19.45
```

### Isophotal Threshold

To determine the size of the galaxy, I first apply the isophotal threshold suggested by [Cortese et al. 2012](https://www.aanda.org/articles/aa/full_html/2012/08/aa19312-12/aa19312-12.html): 

```python
g_band_limit = 24.5
i_band_limit = 23.5
```

Note that my pixel scale is set to 0.262 $arcsec/pixcel$, so I need to convert this limit to magnitude at each pixel:

$$mag/pixel^2 = mag/arcsec^2-2.5\times\log((0.262arcsec/pixcel)^2)$$

```python
g-band isophotal threshold in pixel scale: 27.41 mag/pixel^2
i-band isophotal threshold in pixel scale: 26.41 mag/pixel^2
```

Then I apply this mask to the previous image and get

![image-20250424115252128](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250424 Stellar Mass Map.assets/image-20250424115252128.png)

![image-20250424115710461](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250424 Stellar Mass Map.assets/image-20250424115710461.png)

### Galaxy Contour

After applying isophotal mask, I still need to determine the "actual" size of the target galaxy. 

I can search for the best-fitting ellipse (or other models) to get the size of the galaxy, but my idea is that the morphology of a galaxy is not a perfect shape. Instead of fitting, I use a CV-based segmentation to isolate the traget galaxy from other sources. 

First I convert the isophotal mask to a binary map (i.e., any magnitude below this threshold set to white, otherwise set to black). Then I draw all the contours on this binary map. And I pick the biggest coutour as the interested region. I do this for both g and i bands and create a final galaxy mask: **the region of galaxy consists of any pixel that sit in biggest contour in either g or i bands but must pass both isophotal thresholds**.  

Note. The contours can be tweaked by `erosion` and `closing` to handle some complex cases like the merger. But here i turn these two options off because I think searching contours directly is good enough for our samples. I may underestimate the actual size of the galaxy, because some surrounding pixels that are not inside the biggest contour but may belong to the target galaxy are not counted. However, I believe that these pixels should contribute negligible to the total size of galaxy, so I am happy with the current contour. 

![image-20250424122341852](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250424 Stellar Mass Map.assets/image-20250424122341852.png)

![image-20250424122348411](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250424 Stellar Mass Map.assets/image-20250424122348411.png)

![image-20250424122719103](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250424 Stellar Mass Map.assets/image-20250424122719103.png)

![image-20250424122726713](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250424 Stellar Mass Map.assets/image-20250424122726713.png)

In this case, I can derive the total magnitude of NGC4383:

```python
Total g-band magnitude: 12.389
Total i-band magnitude: 11.723
```

As a reference, these values in Cortese et al. 2012 is 12.392 and 11.654. 

I think g-band magnitude is close enough (0.003 mag underestimated), while i-band is 0.07 mag overestimated. 

### Stellar Mass

Now assuming the distance is $17Mpc$, we can derive the i-band absolute magnitude. In [Taylor et al. 2011](https://academic.oup.com/mnras/article/418/3/1587/1060932#91710255), they purpose an empirical relation between (*g*−*i*) colour, *i*-band luminosity and stellar mass as:

$$\log(M_*/M_\odot) = 1.15+0.7(g-i)-0.4M_i$$

where M_i is the absolute magnitude in the restframe i-band, expressed in the AB system. 

Below is g-i and stellar mass map:

![image-20250424131718373](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250424 Stellar Mass Map.assets/image-20250424131718373.png)

![image-20250424131723553](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250424 Stellar Mass Map.assets/image-20250424131723553.png)

Here is two ways to get total stellar mass, one is summing up pixel-by-pixel, another is derived from total magnitude: 

```python
Stellar mass sum up in NGC4383: 3.31e+09 M_sun or log (M_*/M_sun) = 9.52
Stellar mass total in NGC4383: 2.44e+09 M_sun or log (M_*/M_sun) = 9.39
```

I take the second one as the true value, because: 

1. mathematically, the convertion formula is non-linear, which means $\sum_{i}M(Flux_i) \neq M(\sum_{i}Flux_i)$;
2. this empirical relation is calibrated by the total magnitude from the galaxy
3. when we calcualte the flux from a single pixel, we actually include the flux emitted from adjacent pixels, and so integrating pixel-by-pixel leads to overestimated total magnitude, and thus much higher stellar mass in this case. It is a binning effect. 

But overall, they dont differ a lot. 

As a reference, the stellar mass in [Table 6](https://www.aanda.org/articles/aa/full_html/2012/08/aa19312-12/T3.html) of [Cortese et al. 2012](https://www.aanda.org/articles/aa/full_html/2012/08/aa19312-12/aa19312-12.html) is $\log(M_*/M_\odot)=9.42$. For a fair comparison, I also estimate it by the formula adpoted in [Cortese et al. 2012](https://www.aanda.org/articles/aa/full_html/2012/08/aa19312-12/aa19312-12.html): 

$$\log(M_*/L_i) = -0.963+1.032(g-i) = 0.869+1.032(g-i)-0.4M_i$$

and it returns $\log(M_*/M_\odot)=9.33$. 

I think $\sim0.1dex$ difference is primary due to underestimation of i-band magnitude by $0.07mag$. 

Similarly, if I take the magnitude from [Cortese et al. 2012](https://www.aanda.org/articles/aa/full_html/2012/08/aa19312-12/aa19312-12.html) and apply [Taylor et al. 2011](https://academic.oup.com/mnras/article/418/3/1587/1060932#91710255), this will give $\log(M_*/M_\odot)=9.47$, which is $0.05 dex$ higher than its origianl value. In fact, if we can rule out the differences in magnitude, difference between using these two models to estimate total stellar mass is within $0.08dex$ for total $g-i$ in the range of 0.6 to 1.1 ([Taylor et al. 2011](https://academic.oup.com/mnras/article/418/3/1587/1060932#91710255) is less steeper and these two models intersect at $g-i=0.85$). 

In summary, my results:

```python
Total g-band magnitude: 12.389
Total i-band magnitude: 11.723
Total absolute g-band Magnitude: -18.76
Total absolute i-band Magnitude: -19.43
Stellar mass sum up in NGC4383: 3.31e+09 M_sun or log (M_*/M_sun) = 9.52
Stellar mass total in NGC4383: 2.44e+09 M_sun or log (M_*/M_sun) = 9.39
Stellar mass total in NGC4383 (Corteste et al. 2012): 2.13e+09 M_sun or log (M_*/M_sun) = 9.33
```

HRS142 in [Cortese et al. 2012](https://www.aanda.org/articles/aa/full_html/2012/08/aa19312-12/aa19312-12.html) (g & i-band magnitude and $\log(M_*/M_\odot)$): 

```python
12.392, 11.654, 9.42
```

## IC3392 (HRS172)

![image-20250424140135481](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250424 Stellar Mass Map.assets/image-20250424140135481.png)

My results: 

```python
Total g-band magnitude: 12.631
Total i-band magnitude: 11.623
Total absolute g-band Magnitude: -18.52
Total absolute i-band Magnitude: -19.53
Stellar mass sum up in IC3392: 4.90e+09 M_sun or log (M_*/M_sun) = 9.69
Stellar mass total in IC3392: 4.65e+09 M_sun or log (M_*/M_sun) = 9.67
Stellar mass total in IC3392 (Corteste et al. 2012): 5.26e+09 M_sun or log (M_*/M_sun) = 9.72
```

HRS172 in [Cortese et al. 2012](https://www.aanda.org/articles/aa/full_html/2012/08/aa19312-12/aa19312-12.html) (g & i-band magnitude and $\log(M_*/M_\odot)$): 

```python
12.739, 11.660, 9.77
```

## NGC4501 (HRS190)

![image-20250424140015184](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250424 Stellar Mass Map.assets/image-20250424140015184.png)

My results: 

```python
Total g-band magnitude: 9.809
Total i-band magnitude: 8.774
Total absolute g-band Magnitude: -21.34
Total absolute i-band Magnitude: -22.38
Stellar mass sum up in NGC4501: 7.60e+10 M_sun or log (M_*/M_sun) = 10.88
Stellar mass total in NGC4501: 6.69e+10 M_sun or log (M_*/M_sun) = 10.83
Stellar mass total in NGC4501 (Corteste et al. 2012): 7.72e+10 M_sun or log (M_*/M_sun) = 10.89
```

HRS190 in [Cortese et al. 2012](https://www.aanda.org/articles/aa/full_html/2012/08/aa19312-12/aa19312-12.html) (g & i-band magnitude and $\log(M_*/M_\odot)$): 

```python
10.011, 8.832, 10.98
```

## NGC4192 (HRS91)

![image-20250424140405887](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250424 Stellar Mass Map.assets/image-20250424140405887.png)

My results: 

```python
Total g-band magnitude: 10.471
Total i-band magnitude: 9.436
Total absolute g-band Magnitude: -20.68
Total absolute i-band Magnitude: -21.72
Stellar mass sum up in NGC4192: 4.40e+10 M_sun or log (M_*/M_sun) = 10.64
Stellar mass total in NGC4192: 3.64e+10 M_sun or log (M_*/M_sun) = 10.56
Stellar mass total in NGC4192 (Corteste et al. 2012): 4.20e+10 M_sun or log (M_*/M_sun) = 10.62
```

HRS91 in [Cortese et al. 2012](https://www.aanda.org/articles/aa/full_html/2012/08/aa19312-12/aa19312-12.html) (g & i-band magnitude and $\log(M_*/M_\odot)$): 

```python
10.560, 9.461, 10.65
```

## Summary

In short, the cutoff taken from [Legacy Survey Viewer](https://www.legacysurvey.org/viewer) is already calibrated, but I still need to correct for the galactic extinction. Then I apply isophotal threshold and CV-based contours to determine the size of the galaxy. I adopt the relation from [Taylor et al. 2011](https://academic.oup.com/mnras/article/418/3/1587/1060932#91710255) to estiamte the stellar mass and finally compare my results with [Cortese et al. 2012](https://www.aanda.org/articles/aa/full_html/2012/08/aa19312-12/aa19312-12.html). I found:

1. My total magnitude measurement differs from [Cortese et al. 2012](https://www.aanda.org/articles/aa/full_html/2012/08/aa19312-12/aa19312-12.html) up to $0.2dex$, but in most cases within $0.1dex$. I will say they are consistent.
2. My total stellar mass is underestimated by $\sim0.15dex$. If I apply the same formula in [Cortese et al. 2012](https://www.aanda.org/articles/aa/full_html/2012/08/aa19312-12/aa19312-12.html), it will reduce to within $0.1dex$; while switching model in these 4 samples ($0.7<g-i<1.0$, and my $g-i$ values are $0.1\sim0.2mag$ smaller than [Cortese et al. 2012](https://www.aanda.org/articles/aa/full_html/2012/08/aa19312-12/aa19312-12.html)) will lead to $0.05\sim0.06dex$ differences in stellar mass. It seems to me that correct measurement of magnitude matters in my samples. However, in the ideal case that the magnitudes are correct, changing model can give very different stellar mass if $g-i$ values vary a lot. For example, in GAMA data, $g-i$ values range from -0.2 to 1.6, and that results in up to $0.35dex$ differences for estimating by these two models, especially in lower end. Since [Taylor et al. 2011](https://academic.oup.com/mnras/article/418/3/1587/1060932#91710255) is less steep, it will overestimate the low mass end but underestimate the high mass end. Also, my feeling is that, if $g-i$ is near 0.85 (e.g., 0.7 to 1.0), then accurate measurement of magnitudes is important; if $g-i$ differs a lot from 0.85, then models matter. For me, I prefer [Taylor et al. 2011](https://academic.oup.com/mnras/article/418/3/1587/1060932#91710255) because this calibration is less sensitive to the magnitude changes. 
3. Now the binning effect only change the stellar mass within a factor of 1.2. 