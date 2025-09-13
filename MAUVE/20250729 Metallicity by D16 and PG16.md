# 20250729 Metallicity by D16 and PG16

Before everything, I guess the [N II]$\lambda6484$ line I seen in some papers is actually the historical typo of [N II]$\lambda6584$ in D16, right?

## D16 & PG16 to get [O/H] and comparison with Adam's paper

Below is copy from Adam's paper: "We adopt two different metallicity calibrations to ensure the robustness of our results. The first is the calibration of Dopita et al. ([2016](https://doi.org/10.1007/s10509-016-2657-8)), which uses the [N II]$_{6583}$⁠, H$\alpha$, and [S II] lines. The second is the ‘Scal’ calibration from Pilyugin & Grebel ([2016](https://doi.org/10.1093/mnras/stw238)), which uses the ⁠[O III], H$\beta$, ⁠[N II], and [S II] lines, and we use [O III] = 1.33[O III]$_{5007}$ and ⁠[N II] = 1.34[N II]$_{6583}$."

I apply exactly the same way to NGC4383 to see if I can recreate the metallicity gradient map shown in Figure 7 of Adam's paper. For D16, I use 5th order polunomial calibration.  

Note that I dont see why it is 1.33 and 1.34 here, is it a common sense here? Also it seems that Adam use a customized colormap. 

![image-20250729192936630](assets/image-20250729192936630.png)

I also show the maximum value inside a 2-arcsec square at the galaxy’s centre for each calibration: 8.661 and 8.442 (Adam says 8.65 and 8.44). Simiarly, I see D16 reports more prominent gradient. Also, [O/H] by D16 seems to span in broaded range:  

```python
O_H_D16_SF - Min: 8.0765, Max: 8.6614
O_H_PG16_SF - Min: 8.1901, Max: 8.4423
```

## All 14 galaxies

So NGC4383 looks fine for me and I am sure my code is correct, but thing is getting crazy when applying to other galaxies. 

![image-20250729195031269](assets/image-20250729195031269.png)

![image-20250729195038799](assets/image-20250729195038799.png)

![image-20250729195045568](assets/image-20250729195045568.png)

![image-20250729195052881](assets/image-20250729195052881.png)

![image-20250729195112446](assets/image-20250729195112446.png)

![image-20250729195120069](assets/image-20250729195120069.png)

![image-20250729195131098](assets/image-20250729195131098.png)

![image-20250729195140260](assets/image-20250729195140260.png)

![image-20250729195146882](assets/image-20250729195146882.png)

![image-20250729195154136](assets/image-20250729195154136.png)

![image-20250729195203063](assets/image-20250729195203063.png)

![image-20250729195211631](assets/image-20250729195211631.png)

![image-20250729195219640](assets/image-20250729195219640.png)

![image-20250729195227346](assets/image-20250729195227346.png)

Ok, first I realise that they must be only valid in proper [O/H] range. 

It seems that D16 only say they are linear from 7.6 to 9.05, while PG16 are very confident: “The oxygen and nitrogen abundances estimated through the suggested calibrations agree with the $\rm T_e$‑based abundances within $\sim0.1$ dex over the whole metallicity range of H II regions.” But I think as suggested by Kewley+2019, we may limit both between 7.63 and 9.23. But here is the question, for those outside this range, should we just simply drop them or maybe still assign some values to them and call them upper/lower limits? 

## Redo Calculation

My idea is remove those points do rerun my code and indeed looks better. 

![image-20250729212217010](assets/image-20250729212217010.png)

![image-20250729212229888](assets/image-20250729212229888.png)

![image-20250729212239904](assets/image-20250729212239904.png)

![image-20250729212252402](assets/image-20250729212252402.png)

![image-20250729212305046](assets/image-20250729212305046.png)

![image-20250729212323603](assets/image-20250729212323603.png)

![image-20250729212336438](assets/image-20250729212336438.png)

![image-20250729212348761](assets/image-20250729212348761.png)

![image-20250729212401800](assets/image-20250729212401800.png)

![image-20250729212415910](assets/image-20250729212415910.png)

![image-20250729212429091](assets/image-20250729212429091.png)

![image-20250729212505526](assets/image-20250729212505526.png)

![image-20250729212513935](assets/image-20250729212513935.png)

![image-20250729212526357](assets/image-20250729212526357.png)

## How about other calibration?

It looks MUSE data is still sufficient to use direct $\rm T_e$ method instead of those strong line calibrations? 