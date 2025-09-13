# 20250615 No More Cut at 15

**No more $\rm Flux/Flux_{err}\geq15$ and no more foreground stars**: 

1. $\rm Flux/Flux_{err}$ does matter in terms of quality control, but it doesnt matter in terms of the approach to do quality control. In short, we are doing this by not doing this. 

2. And I have confirmed that those galaxies with the "two-line" feature (i.e., NGC4298 and NGC4694) are just due to the foreground stars. 

## Go to BPT directly and do SF selection

The quality control previously, $\rm Flux/Flux_{err}>15$, is important for removing bad pixels (i.e., not fitting well by `pPXF`), but by doing this we may: 

1. manually introduce an artificial selection effect in resolved SFMS if this cut is too conservative
2. lost some potential star-forming spaxels if this cut is too conservative
3. still remain some bad spaxels if this cut is too low. 

Instead, to still have quality control but not having these issues, we can just simply use BPT diagram and Balmer Decrement (BD) to do such selection. Below is the description of my current apporach: Here I take the example of NGC4298. 

![image-20250615222704533](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250615222704533.png)

Notice that there is a foreground star, and it also behaves as very low BD:

![image-20250615223658656](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250615223658656.png)

First, after reading the data, we directly go to BPT diagram:

![image-20250615202831304](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250615202831304.png)

```python
Galaxy flux-weighted representative points:
[N II] BPT: logN2 = 0.003, logO3 = 1.583
[S II] BPT: logS2 = -0.189, logO3 = 1.192
Number of spaxels in [N II] BPT regions:
HII: 182075, Comp: 112546, AGN: 61438
Number of spaxels in [S II] BPT regions:
HII: 235646, Seyfert: 26197, LINER: 85916
```

And we still select the SF region (i.e., HII+Comp in [N II] BPT and HII in [S II] BPT)

![image-20250615223803505](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250615223803505.png)

## BD selection

![image-20250615225342885](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250615225342885.png)

![image-20250615225403110](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250615225403110.png)

Unfortunately the foreground star is till there after SF selection, so now I need to do the cut in BD. 

**Here, I set a threshold at 2.5. ** I will drop all the spaxels that $BD<2.5$, which are "bad" spaxels. Any spaxels higher than 2.86, I will say they are "good". As for those "between", I will say I will consider them as the spaxels that they maybe actually $BD=2.86$, but just unluckily get $BD<2.86$ because `pPXF` does not fit well. That means there are some spaxels are actually "good" with $BD=2.86$, and I torelate these variations down to $BD=2.5$. 

**I know this 2.5 threshold is still mannualy setting a value. I will porpose another cut based on statistcs later, but I still use this 2.5 to show the ability of removing the "bad" spaxels, including the foreground stars.**

Below i show the histogram of "SNR" (actually the $\rm Flux/Flux_{err}$) for all spaxels and SF region only. 

![image-20250615231343491](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250615231343491.png)

![image-20250615231349419](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250615231349419.png)

So these histograms here provide some information:

1. The previous $\rm Flux/Flux_{err}$ is not a good approach. Some "bad" spaxels can have very high SNR; while some "good" spaxels can have very low SNR. 
2. The SF selection indeed performs a quallity control, but I notice that it seems not that effectively for high SNR in H$\beta$. Again, I suspect that the fitting at H$\beta$ in `pPXF` does not perform well. 
3. SF selection effectively removes "bad" spaxels with high SNR in [S II] doublets. So if in the futrure we only want to use one BPT diagram instead of both, I will say keep using [S II] BPT. 

 Anyway, now I apply the 2.5 threshold and optain the updated BD map:

![image-20250615234413007](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250615234413007.png)

The green spaxels are the "between" spaxels that I still trust them and I will set them to $BD=2.86$; while the blue spaxels are the "bad" spaxels, so the foreground star now is removed. 

As for the rest, everything remains exactly the same as before: getting $E(B-V)_{BD}$, getting corrected H$\alpha$ flux, converting to SFR and finally $\Sigma_{SFR}$. 

## The resolved SFMS again

First, as a comparison, I show the previous plot:

![image-20250615235108336](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250615235108336.png)

Now this is the plot using current approach

![image-20250615235045184](/Users/Igniz/Desktop/ICRAR/MAUVE/assets/image-20250615235045184.png)



Clearly, no more foreground stars in NGC4298 and NGC4694 after BD selection and therefore no more "two-line" feature. Also, the foregound star in NGC4064 is removed, though it was already disappeared by conservative $\rm Flux/Flux_{err}$ cut in the previous plot. On the other hand, without $\rm Flux/Flux_{err}$ cut, I have more SF spaxels which are dropped due to previous conservative $\rm Flux/Flux_{err}$ cut.  

However, there are still a few things I need to investigate:

1. IC3392. I notice that now many spaxels in the outer region (especially the edge of galaxy, or in other words, very low H$\alpha$ flux and SFR) are classified as SF, but are they really photonized by OB stars? I guess they are DIG? It seems to me that in my new approach, I include some DIG or maybe something else that are really SF area. In fact, by simpliy excluding the "Comp" part in [N II] BPT diagram, I can remove most of them, but again, this falls back to too conservative cut and will lost some spaxels in the inner region. I am thinking of alternative way to do with them. In fact, some other galaxies also have this phenomenon. 
2. NGC4192. It shows two components in the resolved SFMS plot. However, I have my preliminary guess: they are real, so the lower and major component is the outer region; while the upper and minor component is the central region. That says, unlike NGC4698 that almost zero SF area in the inner region, NGC4192 remains some SF area at the galactic centre. And I guess that is also why NGC4698 is classified as Seyfert galaxy, while NGC4192 is not. 

In fact, I will look at each galaxy, but will start by these two. 

## A more meaningful BD selection

As mentioned before, the threshold at 2.5 for BD selection falls back to the question that why 2.5 is chosen? Here, I come out with a idea that we can adopt this threshold by the data itself. 

The reason why I set this threshold at a value lower than the intrinsic 2.86 is that I allow the uncertainty from `pPXF`, i.e., the fitting uncertainties of both H$\alpha$ and H$\beta$ lines.  Thus, we can use uncertainty of this 2.86, to quantify the question: to what extend we can allow this uncertainty? And this leads to the error propagation. 

For the function that $f = \frac{A}{B}$, the standard deviation of $f$ is defined by
$$
\sigma_f = |f|\sqrt{(\frac{\sigma_A}{A})^2+(\frac{\sigma_B}{B})^2-(\frac{\sigma_{AB}}{AB})^2}
$$
Since H$\alpha$ and H$\beta$ lines are fitted independently so the covariance term $(\frac{\sigma_{AB}}{AB})^2$ is dropped. 

Hence, we have
$$
\begin{align}
\sigma_f &= |f|\sqrt{(\frac{\sigma_A}{A})^2+(\frac{\sigma_B}{B})^2} \\
&= 2.86\sqrt{(\frac{\sigma_{A}}{A})^2+(\frac{\sigma_B}{B})^2} \\
&= 2.86\sqrt{(\frac{<\sigma_{H\alpha}>}{F_{H\alpha}})^2+(\frac{<\sigma_{H\beta}>}{F_{H\beta}})^2} \\
&= 2.86\sqrt{(\frac{\frac{\Sigma Err_{H\alpha}F_{H\alpha}}{\Sigma F_{H\alpha}}}{\frac{\Sigma F_{H\alpha}^2}{\Sigma F_{H\alpha}}})^2+(\frac{\frac{\Sigma Err_{H\beta}F_{H\beta}}{\Sigma F_{H\alpha}}}{\frac{\Sigma F_{H\beta}^2}{\Sigma F_{H\beta}}})^2} \\
&= 2.86\sqrt{(\frac{\Sigma Err_{H\alpha}F_{H\alpha}}{\Sigma F_{H\alpha}^2})^2+(\frac{\Sigma Err_{H\beta}F_{H\beta}}{\Sigma F_{H\beta}^2})^2}. \\
\end{align}
$$
And finally, we allow down to $3\sigma_f$ variation from the intrinsic value. Normally, this $\sigma_f$ is bwtween 0.05 to 0.1. I will switch to this approach if it makes sense. 
