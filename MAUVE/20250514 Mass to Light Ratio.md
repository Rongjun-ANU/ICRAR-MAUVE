# 20250514 Mass to Light Ratio

## Including Dust Attentuation in `pPXF` fitting

Note that for now, regularization parameter `regul` is turning off. I will discuss it later. 

Previously, I end up with r band $M_*/L_* = 2.139 = \log(0.330)$. 

```python
 Best Fit:       Vel     sigma
 comp.  0:      1675        59
chi2/DOF: 888.9; DOF: 1833; degree = 12; mdegree = 0
method = capfit; Jac calls: 4; Func calls: 14; Status: 4
linear_method = lsq_box; Nonzero Templates (>0.1%): 6/150
(M*/L)=2.139 (SDSS/r at z=0.0056)
nodust | regul=    0 | χ²/DOF= 888.87 | M/L_r=  2.14 = log(  0.330) | 
```

![image-20250513191650122](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250514 Mass to Light Ratio.assets/image-20250513191650122.png)

![image-20250513191658715](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250514 Mass to Light Ratio.assets/image-20250513191658715.png)

Considering that dust attentuation will largely obscure young stellar population, I turn on the `dust` parameter in `pPXF` fitting to account for this effect. And this gives r band $M_*/L_* = 2.068 = \log(0.316)$. 

```python
 Best Fit:       Vel     sigma
 comp.  0:      1676        65
Attenuation Parameters 0: 0.935 -1.000
chi2/DOF: 831.2; DOF: 1831; degree = 12; mdegree = 0
method = capfit; Jac calls: 4; Func calls: 22; Status: 4
linear_method = lsq_box; Nonzero Templates (>0.1%): 10/150
(M*/L)=2.068 (SDSS/r at z=0.0056)
dust   | regul=    0 | χ²/DOF= 831.18 | M/L_r=  2.07 = log( 0.316) | 
```

![image-20250513192439477](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250514 Mass to Light Ratio.assets/image-20250513192439477.png)

![image-20250513192444898](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250514 Mass to Light Ratio.assets/image-20250513192444898.png)

So the dust attenuation parameters here are: $A_V = 0.935$ and $\delta = -1$. 

![image-20250513193656364](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250514 Mass to Light Ratio.assets/image-20250513193656364.png)

So $A_V = 0.935$ looks reasonable for me, but $\delta = -1$ is even steaper than SMC-like galaxies. That means IC3392 is very metal-poor ([Shivaei et al. 2020](https://iopscience.iop.org/article/10.3847/2041-8213/abc1ef))? 

## Different $M_*/L_*$

Recall that using Legacy data and applying [Taylor et al. 2011](https://doi.org/10.1111/j.1365-2966.2011.19536.x)'s approach, I get r band $M_*/L_* = 1.35 = \log(0.13)$. Also if take $g-i = 1.008$ and use [Zibetti et al. 2009](https://doi.org/10.1111/j.1365-2966.2009.15528.x)'s calibration, it will be r band $M_*/L_* = 1.54 = \log(0.19)$. The former use `BC03` templates while the latter use `CB07` (2007 version of BC03). In contract, in `pPXF` I adopt `E-MILES`, so you can see in general they differs by $0.13\sim0.2$ dex. 

In Figure E.15 of [Pessa et al. 2023](https://www.aanda.org/articles/aa/full_html/2023/05/aa45673-22/aa45673-22.html), they find stellar mass surface density $\log(\Sigma_*)$ of the star-forming region derived using `E-MILES` is roughly 0.2 dex higher than `CB07`, regardless of adding nebular correction and removing most metal-poor templates. But in general, all templates produce unexpected very metal-poor artefacts ($[Z/H] \lesssim −1.3$) in the star-forming ring. 

However, in [Lee et al. 2025](https://iopscience.iop.org/article/10.3847/1538-3881/adb285), they show different resluts by compareing `E-MILES`, `BC03`, `CB19` (2019 version of `CB07`), and `FSPS`. In Figure 1, they show that the $M_*/L_*$ curves are nearly coincident in SDSS r band across different templates; while in $3.6\mu m$, `BC03` and `CB19` lie $\approx$ 0.25 dex below `E-MILES`/`FSPS` from $\log(8.7-10.2)M_\odot$. 

I notice that they only compare in the range of $\log(8.7-10.2)M_\odot$, but in fact `E-MILES` models lack spectral templates for stellar ages younger than 63 Myr. I think this may `pPXF` tends to estimate the weights (and therefore $M_*/L_*$) in a higher parameter space of $\log(Age/yr)$. Indeed, if swtich to `GALAXEV` (updated 2016 version of the `BC03` templates), the weight pattern will be lefter than `E-MILES`: 

![image-20250514104535450](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250514 Mass to Light Ratio.assets/image-20250514104535450.png)

And this gives $M_*/L_* = 1.91 = \log(0.290)$, so it is slightly decrease but still higher than color-color relation. And it seems kind of breaking the age-metallicity degeneracy because of no more young but very metal poor stars. 

Switching different templates still not accounts for the gap between them, so I guess the smoking gun is the different approaches I adopted, i.e. colour-color relation tends to underestimate the $M_*/L_*$. [Ge et al. 2021](https://academic.oup.com/mnras/article/507/2/2488/6350571#298387008) concluded that we should be careful with using $M_*/L_*$ color relation for those galaxies with young luminosity-weighted stellar ages. My understanding is that color relation will be intrinsically biased by the luminosity and misrepresents older stars with younger populations. Therefore, it should be safer to use `pPXF` if data are eligible. 

One more thing, should we customize SPS templates and adopt mass-weighted properties? 