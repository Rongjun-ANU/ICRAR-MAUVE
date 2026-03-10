# 20260310 Skyline Spectral Masking

## 1. Creating Skyline Spectral Masking

We confirmed that we will use `DEG` = 23 for the upcoming $4800-9100\AA$ runs.

![image-20260310173748588](assets/image-20260310173748588.png)

Below is my way to create the Skyline Spectral Masking for $\geq 6800\AA$: 

1. detect all the peaks in the sky spectrum with flux >= 50 at $6800-9200\AA$

2. then we distinguish two groups of peaks: green is already masked in `specMask_KIN`, while red is not masked yet.

3. For green peaks, then use the existing `specMask_KIN` to determine the masking, while for red peaks, use the detected peak wavelength as the center. As for the width, we widen all the $5\AA$ linewidth to $6\AA$ or $7\AA$ (including the exising ones in `specMask_KIN` and the newly detected ones). 

4. As for the special `7628.00             84             sky` (telluric feature) in `specMask_KIN`, update its range from $84\AA$ to $96\AA$.

5. We further mask the noisy $7200-7240\AA$ region, which seems to minimise the bump right before [Ar III]$\lambda7135$. 

6. The original `specMask_KIN` only contains [Ar III]$\lambda7135$ for those emission lines $\geq 6800\AA$. So I checkout the `emissionLinesPHANGS_full.config` (my modified version that turn on the fitting for all emissioin lines $\geq 6800\AA$ in original `emissionLinesPHANGS.config`), and add those gas lines, with their widths all set to $20\AA$ as existing lines:

   ```python
       ("OII", 7318.92),
       ("OII", 7319.99),
       ("OII", 7329.67),
       ("OII", 7330.73),
       ("HPa11", 8862.89),
       ("HPa10", 9015.3),
       ("SIII", 9068.6),
   ```

So now I have `specMask_KIN_narrow6` and `specMask_KIN_narrow7`. 

![image-20260310180249341](assets/image-20260310180249341.png)

![image-20260310175325579](assets/image-20260310175325579.png)

## 2. Checking the fitting results

Here i show the original `MDEG` =23, and the cases with `test3`,  `specMask_KIN_narrow6` and `specMask_KIN_narrow7` masks of means from our spatial mask regions. I can see slight improve of reducing two bumps from original mask. It seems that `test3` provides the worst residual, so probably due to too wide masking. 

![image-20260310201852476](assets/image-20260310201852476.png)

Then I show the same plots from both central (high S/N) and faint (low S/N) regions. 

![image-20260310205418446](assets/image-20260310205418446.png)

![image-20260310205434138](assets/image-20260310205434138.png)

## 3. Checking in `Mapviewer`

Note that now the spectra are in rest frame, so the skyline masking is blue shifted here. 

I use `specMask_KIN_narrow7` as I feel like $6\AA$ still not doesn't cover some edges (mismatch) of skylines, especially $>7500\AA$.  

### 3.1 Bright `Mapviewer`

First, we still look at the central spaxel and it looks fine. The red dot roughly at (0, 0).

![image-20260310210548780](assets/image-20260310210548780.png)

Zoom in to see $6800\sim7000\AA$:

![image-20260310210809525](assets/image-20260310210809525.png)

$7000\sim7150\AA$ (Argon line the bump before it): 

![image-20260310211114677](assets/image-20260310211114677.png)

$7150\sim7500\AA$ (added noisy sky $7200-7240\AA$ region and [O II] lines): 

![image-20260310211628503](assets/image-20260310211628503.png)

$7500\sim7785\AA$ (telluric feature):

![image-20260310211835878](assets/image-20260310211835878.png)

$7785\sim8120\AA$:

![image-20260310212246690](assets/image-20260310212246690.png)

$8120\sim8420\AA$:

![image-20260310212456482](assets/image-20260310212456482.png)

$8420\sim8700\AA$ (CaT absorption features, but dominated by sky emission for IC3392 unfortunately): 

![image-20260310212728229](assets/image-20260310212728229.png)

$8700\sim9100\AA$ (Paschen and [S III] lines): 

![image-20260310213629529](assets/image-20260310213629529.png)

### 3.2 Faint `Mapviewer`

Then, we focus on the isolated upper right bin (still donted by the red dot). Clearly, the spectrum is very noisy but those noisy spikes are masked.

![image-20260310214246029](assets/image-20260310214246029.png)

Zoom in to see $6800\sim7000\AA$:

![image-20260310214350384](assets/image-20260310214350384.png)

$7000\sim7150\AA$ (Argon line the bump before it): 

![image-20260310214628327](assets/image-20260310214628327.png)

$7150\sim7500\AA$ (added noisy sky $7200-7240\AA$ region and [O II] lines): 

![image-20260310214717491](assets/image-20260310214717491.png)

$7500\sim7785\AA$ (telluric feature):

![image-20260310214839506](assets/image-20260310214839506.png)

$7785\sim8120\AA$:

![image-20260310214954367](assets/image-20260310214954367.png)

$8120\sim8420\AA$:

![image-20260310215059357](assets/image-20260310215059357.png)

$8420\sim8700\AA$ (CaT absorption features, but dominated by sky emission for IC3392 unfortunately): 

![image-20260310215146784](assets/image-20260310215146784.png)

$8700\sim9100\AA$ (Paschen and [S III] lines): 

![image-20260310215209718](assets/image-20260310215209718.png)
