# 20250702 Consistency in [NII] and [SII] BPT

## 1. 4 catalogs (updated)

I have the QC for each line: $\rm FLUX/FLUX_{ERR} \geq 3$ and $\rm FLUX \geq 20\times10\ erg/s/cm^2$. 

Now I have divided all non-nan spaxels into 4 catalogs (here I use [NII] BPT as an example):

1. `Upper` (in color blue): not detected in Balmer lines, i.e., either H$\alpha$ or H$\beta$ or both does not pass the QC. 
2. `Unclassified` (in color red): detected in Balmer lines, but error bars are too large to locate the spaxel's position on the [NII] BPT diagram. 
3. `nonSF` (in color orange): detected in Balmer lines, constrainted on [NII] BPT diagram, and above the Kewley+2001 curve (red curve). 
4. `SF` (in color green): detected in Balmer lines, constrainted on [NII] BPT diagram, and below the Kewley+2001 curve (red curve). 

## 2. Consistency on [NII] and [SII] BPT

Now I want to see if `SF` and `nonSF` are consistent on [NII] and [SII] BPT diagrams, so in the `classified` spaxels in both BPT diagrams, I further have: 

```python
mask_SF_N2_SF_S2 = mask_classified_both & mask_SF_N2 & mask_SF_S2
mask_nonSF_N2_nonSF_S2 = mask_classified_both & mask_nonSF_N2 & mask_nonSF_S2
mask_nonSF_N2_SF_S2 = mask_classified_both & mask_nonSF_N2 & mask_SF_S2
mask_SF_N2_nonSF_S2 = mask_classified_both & mask_SF_N2 & mask_nonSF_S2
```

As expected, most of them are `SF` and they are consistent on both [NII] and [SII] BPT (in green), while there is 0 spaxels `SF` on [SII] but `nonSF` on [NII] BPT. However, the case that is `SF` in [NII] but `nonSF` on [SII] BPT become sifnificant in NGC4383, NGC4396, NGC4522. I can see that they are still not far from the boundary between HII and Comp, so I think they may be due to ambigus classification of SF region and harder ionisation.  

### 2.1 IC3392

![image-20250702194804528](assets/image-20250702194804528.png)

### 2.2 NGC4064

![image-20250702194812289](assets/image-20250702194812289.png)

### 2.3 NGC4192

![image-20250702195555012](assets/image-20250702195555012.png)

### 2.4 NGC4293

![image-20250702195012179](assets/image-20250702195012179.png)

### 2.5 NGC4298

![image-20250702195413655](assets/image-20250702195413655.png)

### 2.6 NGC4330

![image-20250702195036918](assets/image-20250702195036918.png)

### 2.7 NGC4383

![image-20250702195046747](assets/image-20250702195046747.png)

### 2.8 NGC4396

![image-20250702195056421](assets/image-20250702195056421.png)

### 2.9 NGC4419

![image-20250702195119399](assets/image-20250702195119399.png)

### 2.10 NGC4457

![image-20250702195242003](assets/image-20250702195242003.png)

### 2.11 NGC4501

![image-20250702195845810](assets/image-20250702195845810.png)

### 2.12 NGC4522

![image-20250702195315856](assets/image-20250702195315856.png)

### 2.13 NGC4694

![image-20250702195332337](assets/image-20250702195332337.png)

### 2.14 NGC4698

![image-20250702195347091](assets/image-20250702195347091.png)
