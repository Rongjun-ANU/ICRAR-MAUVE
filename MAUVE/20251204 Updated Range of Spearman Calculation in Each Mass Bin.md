# 20251204 Updated Range of Spearman Calculation in Each Mass Bin

So previously, we have only selected data within 95% contour to do Spearman test and hope that can remove some bias. 

![image-20251204174113442](assets/image-20251204174113442.png)

So yes the last mass bin still looks problematic.

![image-20251204173742815](assets/image-20251204173742815.png)

And I check that we still get 0.518 just because the extreme outliers. Particacularly, the top right red dots: that bin represents 36 spaxels. 

Hence, instead of select by 95% contour, I clip the data out of 3$\sigma$ in both SFR and O/H, and the result looks as expected.

![image-20251204174410261](assets/image-20251204174410261.png) 

And when I update to other mass bins, results are still similar.

![image-20251204174619183](assets/image-20251204174619183.png)

And so i think we can update to all the galaxies:

![image-20251204174853563](assets/image-20251204174853563.png)

And all indicators:

![image-20251204174915920](assets/image-20251204174915920.png)