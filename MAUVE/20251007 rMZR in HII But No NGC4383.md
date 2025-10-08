# 20251007 rMZR in HII But No NGC4383

Note.

1. Try removing the abnnormal NGC4383.  The reason is that it almost occupies the entire independent parameter space in rMZR.  So if we want to study the SFR dependence in rMZR, it will be biased by this single galaxy.  Therefore, we can try removing it and see how the results change.
2. Only spaxels from HII regions are included.
3. **Update to 6 $\Sigma_\mathrm{SFR}$ bins, which requires same spaxel counts in each $\Sigma_\mathrm{SFR}$ bin (so the binwidth will change in different panels).**
4. Stick with these 4 indicators: D16, PG16, O3N2-M13 and O3N2-C20. 
5. Update to require at least 20 unique datapoints in each median bin of each $\Sigma_\mathrm{SFR}$ bin. 

Copy table from previous notetable (* are those "abnormalies" in rMZR, together with the removed NGC4383):

|    ID    |   Position   |    correlation     | Expectation |
| :------: | :----------: | :----------------: | :---------: |
|  IC3392  |     High     |      Positive      |     Yes     |
| NGC4064  |     High     |      Positive      |     Yes     |
| NGC4192  |    Cross     |      Inverse       |     Yes     |
| NGC4293  | Few (High?)  |         -          |      -      |
| NGC4298  |    Cross     | Positive (mostly?) |     No?     |
| NGC4330* |     Low      |      Negative      |     Yes     |
| NGC4396* |     Low      |      Negative      |     Yes     |
| NGC4419  |     High     |      Positive      |     Yes     |
| NGC4457  | Few (High?)  |         -          |      -      |
| NGC4501  |    Cross     |      Positive      |     No      |
| NGC4522* |     Low      |   - (Positive?)    |     No?     |
| NGC4694* |     High     |      Positive      |     Yes     |
| NGC4698  | Few (Cross?) |         -          |      -      |

## Combined

![image-20251007173533392](assets/image-20251007173533392.png)

## Individual

![image-20251007173600490](assets/image-20251007173600490.png)

![image-20251007173611363](assets/image-20251007173611363.png)

![image-20251007173639613](assets/image-20251007173639613.png)

![image-20251007173716755](assets/image-20251007173716755.png)

![image-20251007173725182](assets/image-20251007173725182.png)

![image-20251007173732391](assets/image-20251007173732391.png)

![image-20251007173741024](assets/image-20251007173741024.png)

![image-20251007173750120](assets/image-20251007173750120.png)

![image-20251007173802274](assets/image-20251007173802274.png)

![image-20251007173809995](assets/image-20251007173809995.png)

![image-20251007173818696](assets/image-20251007173818696.png)

![image-20251007173826584](assets/image-20251007173826584.png)

![image-20251007173834147](assets/image-20251007173834147.png)

## What if we try seperate the combinded rMZR plots into two branches: outliers and majority
Remember, except the NGC4383 that has been excluded, the other abnormal galaxies are: NGC4330, NGC4396, NGC4522 and NGC4694.

1. NGC4330: low mass end, metallicity anti-correlate with SFR, match expected trend
2. NGC4396: low mass end, metallicity anti-correlate with SFR, match expected trend
3. NGC4522: low mass end, can't see clear trend
4. NGC4694: high mass end, metallicity anti-correlate with SFR, match expected trend

![image-20251007174003110](assets/image-20251007174003110.png)

![image-20251007174013265](assets/image-20251007174013265.png)

So now we can may conclude that the anti-correlation between Z and $\Sigma_\mathrm{SFR}$ in low mass end are driven by 3 outliers: NGC4330, NGC4396, NGC4522. Below I combine two above plots together

![image-20251007175708131](assets/image-20251007175708131.png)