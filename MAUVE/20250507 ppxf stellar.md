# 20250507 ppxf stellar

## 1. Data

I am still using IC3392

```python
Cube dimensions  →  nz = 3761,  ny = 438,  nx = 437
```

## 2. Wavelength cutoff and velocity scale

''Sky subtraction is clearly not perfect, but the best that we can do for the moment. Below 7000 Å it is generally acceptable, at longer wavelengths the situation is worse.'' 

So I make a cutoff at $7000\AA$. Actually, now it remains $4750-7000\AA$. 

```python
After λ‑cut (<=7000 Å)  →  nz = 1800
```

Then I compute the velocity scale:
```python
c_kms   = c.c.to(u.km/u.s).value        # 299 792.458
dlnλ    = np.diff(np.log(lam_ang))  # dlnλ in Å
velscale = np.min(c_kms * dlnλ) # km/s per pixel
```

This gives $53.55 km/s$. Ok, but I still have questions, why choose the minmum one, why not using the one return by `log_rebin` below, and what exactly is the meaning of velocity scale?

Here is the native spectrum:

![image-20250507135903876](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250507 ppxf stellar.assets/image-20250507135903876.png)

## 3. Mask emission lines

Since we are only interested in the continuum, I remove the emission lines from galaxy and air:

![image-20250507140112273](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250507 ppxf stellar.assets/image-20250507140112273.png)

## 4. `log_rebin`

Now I need to do `log_rebin`. It seems that this is one of the requirements in `pPXF`. But the question is that, why in natual log rather than $\log_{10}$ (same question for `velscale`)? No need to multiply an extra constant when taking derivative?

I still do `log_rebin` anyway for both flux and noise, and I force `velscale=velscale`, so I got: 

```python
Log‑grid length : 2171 pixels
velscale        : 53.549 km/s
```

![image-20250507140816111](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250507 ppxf stellar.assets/image-20250507140816111.png)

## 5. SPS templates: E-MILES

Then I load the SPS templates. Here I choose `spectra_emiles_9.0.npz` because it seems to be more suitable for IFS data. 

[E-MILES](http://miles.iac.es/) SPS model templates:  [Vazdekis et al. (2016)](https://ui.adsabs.harvard.edu/abs/2016MNRAS.463.3409V).

## 6. FWHM and MUSE LSF

I am still confused starting from here. 

`pPXF` requires thestellar templates and the galaxy spectrum to have the same instrumental resolution before it adds any extra broadening for the LOSVD. And I think this resolution is related to FWHM of the instrument LSF? 

[Emsellem+2022](https://www.aanda.org/articles/aa/full_html/2022/03/aa41727-21/aa41727-21.html) use this equation for MUSE LSF:

$$ \begin{aligned} {FWHM} \, (\lambda \, [{\AA }]) = 5.866\times 10^{-8} \lambda ^2 - 9.187\times 10^{-4}\lambda + 6.040. \end{aligned} $$

But then? 

Convolve the SPS templates to MUSE resolution? To boarden them to avoid sharpening the spectrum in fitting? 

Technically, I am not sure what to do with this.

After that I can `log_rebin` the template and ready for `pPXF` fitting. 