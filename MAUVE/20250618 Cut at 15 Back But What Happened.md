# 20250618 Cut at 15 Back But What Happened

I now add $\rm Flux/Flux_{err} \geq 15$ back and see how does it affect my resolved SFMS. 

This time I use NGC4298 as the example.

## No Foreground Stars

I use a simple mask (based on median$\pm$5$\sigma$ outliters in velocities or velocity dispersions in `*KIN_map.fits`) to remove the foreground stars in NGC4064, NGC4298 and NGC4694. A note here, even after removing the stars, the surrounding areas still seem to be a bit problematic. At least those surrounding areas are still quite small in Balmer Decrement (BD), i.e., <2 or even <1. See the figures next, but luckily they are removed when applying $\rm Flux/Flux_{err} \geq 15$ cut to other lines, i.e., O[III], N[II], S[II] doublets. 

##   $\rm Flux/Flux_{err} \geq 15$ and check the BD map

After applying the quality control cut, I check the BD map. Also, for those spaxels removed by the cut, I assume their BD is 2.86 so that I can have a upper limit map.

![image-20250617224810576](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250617224810576.png)

So there are still some spaxels having BD<2, especially near the postion of the foreground star. This is same case for NGC4064 and NGC4694:

![image-20250617225600523](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250617225600523.png)

![image-20250617225606519](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250617225606519.png)

Also, even though some spaxels have passed the quality control, they are still low in BD. Again, something wrong with `pPXF` and maybe also simply set these values to be 2.86? 

Here I show the cut as green contours in all lines:

![image-20250618111905160](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250618111905160.png)

So clearly, the "SNR" of H$\beta$ dominates $\rm Flux/Flux_{err} \geq 15$ when determing the BD; while "SNR" of O[III] dominates $\rm Flux/Flux_{err} \geq 15$ when determing the SF in BPT diagrams later. 

Here I define two "SNR" selection:

1. `combined`: $\rm Flux/Flux_{err} \geq 15$ for both H$\alpha$ and H$\beta$. 
2. `combined`: $\rm Flux/Flux_{err} \geq 15$ for all 6 lines: H$\alpha$, H$\beta$, O[III], N[II], S[II] doublets. 

## BD and $E(B-V)_{BD}$

![image-20250618114110344](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250618114110344.png)

```python
Lowest Balmer Decrement: 0.3662206123463743
Highest Balmer Decrement: 10.399312980553152
Lowest 5 unique non-NaN Balmer Decrement values: [0.36622061 0.38014533 0.48627324 0.51733327 0.51902884]
```

For those spaxels that are removed by `combined`, I set their BD to be 2.86 (and thus 0 in $E(B-V)_{BD}$). However, some spaxels selected by `combined` are still quite smaller than 2.86. Although these spaxels will be further removed by `combinedd` (not SF selection in BPT diagrams), I am thinking do we also need to set them as 2.86 so that they will not cause further confusion.

Then I use $E(B-V)_{BD}$ to correct the all the line fluxes. 

## BPT diagrams

![image-20250618121733432](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250618121733432.png)

Note, the points here on BPT diagrams are only those spaxels selected by `combined`, but when I determing the `SF` selection, I use all the spaxels so that I can combine with `combinedd` selection to see which SF spaxels are truely SF under `SF` & `combinedd` selection (left panel in the next figure); while others are "potentially" SF, i.e., selected by `SF` but excluded by `combinedd` selection (right panel in the next figure). 

![image-20250618122501686](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250618122501686.png)

Below i show the $E(B-V)_{BD}$ maps of SF spaxels and all spaxels. 

![image-20250618122741771](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250618122741771.png)

```python
Lowest Balmer Decrement: 3.1286756389787187
Highest Balmer Decrement: 10.399312980553152
Lowest 5 unique non-NaN Balmer Decrement values: [3.12867564 3.12929884 3.14896565 3.17144341 3.17667899]
```

So clearly no more negative $E(B-V)_{BD}$ in SF spaxels compared to those under `combined` selection, and I can confirm that these are primarily removed by `combinedd` not `SF` selection. 

## Flux -> Luminosity -> SFR -> $\Sigma_{\rm SFR}$

Now we have three different data:

1. SF spaxels (`_SF`): extinction corrected, under `SF` & `combinedd` selection. 
2. Potential SF spaxels (`_out_SF`): extinction corrected, under `SF` & `~combinedd` selection. They are shown as `SF` in BPT diagrams but removed due to `combinedd` selection. 
3. All spaxels (`_upper`): only extinction corrected. `_upper` - (`_SF`+`_out_SF`) to be all other spaxels. 

I perform the same procedure to get SFR surface density:

![image-20250618124516945](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250618124516945.png)

## Same selection for $\Sigma_*$

![image-20250618124610734](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250618124610734.png)

## Resolved SFMS

![image-20250618140402802](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250618140402802.png)

Again, blue dots are "true" SF spaxels (`SF` & `combinedd` selection), green dots are "potential" SF spaxels (under `SF` & `~combinedd` selection), and red dots are other spaxels that are definitely not SF, respectively. 

The purpose of this plot is to show what happened when applying the $\rm Flux/Flux_{err} \geq 15$ cut and see if this cut biases our resolved SFMS.  

It seems that the green dots eroses the blue dots, but the answer is i dont know, because i really dont know to what extend I can say I am happy with this $\rm Flux/Flux_{err}$ value. So that means I need to play a bit with this value.  

Below I show the cases that the cut = 25, 20, 15, 10, and 5. 

![image-20250618145725055](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250618145725055.png)

## ![image-20250618141044086](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250618141044086.png)

![image-20250618145804951](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250618145804951.png)

![image-20250618141237815](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250618141237815.png)

![image-20250618145645843](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250618145645843.png)

I can see that as the cut decreases, more green spaxels are scaled up in resolved SFMS diagrams compared to previous high cut case. So these spaxels are actually having BD>2.86, but due to conservative cut, they are set to BD=2.86. When we loose the restrictions, more spaxels pass the selection and preserve their BD>2.86, so after extinction correction, they are shifted to higher SFR. 