# nGIST tasks

Here I am listing some instructions and to-do-list to start getting up to speed with nGIST while I am away and gradually help with improving the way we create data products. 

## Go-to places 

- nGIST code on github - [ngist/CHANGELOG at dev-branch · geckos-survey/ngist](https://github.com/geckos-survey/ngist/blob/dev-branch/CHANGELOG) note that you should use the "dev-branch" not the main branch, as the main branch is now quite old (don't ask 😭 )
- Check the dwg-gist-development channel for all the most updated info on changes. In case of doubt, feel free to contact Jesse van de Sande. Also the #GIST set-up channel on the MAUVE wiki includes some early conversations with a CANVAS listing the original choices made.
- A good starting point to understand most of the changes implemented since our last data release is Amelia's paper [The GECKOS survey: Identifying kinematic sub-structures in edge-on galaxies](https://ui.adsabs.harvard.edu/abs/2025A%26A...700A.237F/abstract)

## Things you could help with (sorted in order of complexity - and the order I would personally follow)

- Create nGIST compatible masks for artefacts/background galaxies/foreground stars for the new cubes. As for some cubes the pixel grid has changed, the ones you already put together may not be suitable. My suggestion is to create a nice collage where - for each cube - you show the color image before and after applying the mask, so that you are sure you are masking most of the clear bck/foreground objects. Also, make sure the masks are in the correct format and try to run nGIST on at least one of them to check that it runs.
- Check that the default mask of skylines in nGIST properly mask the sky (i.e., the masks are not too narrow) for MAUVE galaxies. The idea is to run simply the kinematic module of nGIST and then use the "mapviewer" tool to check that the sky lines are properly masked. Note that below 7000 this is already the case, but I have not tested the default (yet) between 7000 and 9000A. Since it would make sense to create the next iterations of products using the wider wavelength range it is worth checking. Make sure you check spaxels with relatively low S/N as the high S/N spaxels may not show any sky residuals. 
- Improve the quality of the continuum subtracted cubes. Continuum substracted cubes are a very useful data products. However, we are not entirely happy (yet) with the final result. This could be due to masking (see point above), choice of polynomial fit or stellar population templates. So, I would suggest the following approach.
  - Test different MDEG polynomial on the 4800-7000 A wavelength range. We used MDEG 8 in the past and a general rule of thumb is 1 order per 300 A, but I would compare these results with higher orders and see how things go. Elements to check are a) if the continuum where there are no lines is around zero or there is residual flux; b) see if there are artefacts close to the lines. c) compare the flux in the various versions. d) if things to do not dramatically change, you could test moving from the "EMLINES" version of the templates ot the "SAFE" one. The code will take much longer to run!
  - Once you have tested the 4800-7000 A, I would expand the test to the full wavelength range and see if things improve or not (of course, you will have to increase the MDEG order. 

## How to actually get things running with nGIST

Since our last internal DR in December 2024, nGIST has gone through a major upgrade. As our focus has been on the data reduction, we haven't yet updated the entire nGIST workflow. So, setting things up for you may take some time (sorry for this!) as things are available in different places. 

\* Important - please do not edit files already present either on CANFAR or OzSTAR but create your own copy *

Step 1: Install your own nGIST version on CANFAR. For the testing phase, I do suggest to work on CANFAR. It is not as fast as OzSTAR but it is "free".

Step 2: Work on the set-up of the nGIST configuration file. I wrote a simple script to create the nGIST configuration file for every galaxy run in the most automatic way to avoid mistakes. You find everything under /arc/projects/mauve/products/configFiles/create_config 

There is a README file but, briefly, the "GIST_setupinput" has all the input parameters needed (e.g., velocity, Ebv, etc.); the "make_gist_config.py" is the actual script; and "MAUVE_Master_Config..." is the template .yam config file for nGIST. This last file will need to be changed and the overall script to be updated as the new nGIST has a different config file (you will also have to update the script accordingly). Please remember to create copies and not to touch the original folder/file.

Step 3: Set-up the stellar population models. Before running nGIST, you need to set-up your stellar population models for the various modules. I think you have already MILES and in the "ngist_supplementary_public" github repository you should find some libraries. However, probably the easiest way to setup things is to copy the libraries I have assembled on OzStar for some recent tests I made. Under /fred/oz084/mauve/ngist_supplementary_public/gistTutorial/spectralTemplates/ you will find a lot of directories. The ones to start with are the "EMLINES" and "safe" versions for both EMILES and MILES. Do not copy the gz but copy the folders directly and check that the range of metallicities and ages is the same for EMILES and MILES versions. Going through the GECKOS nGIST Slack channel you may notice that they have tested other stellar population models that may be available on github. Tests them only as last resource. 

I think this is all you need to run nGIST on CANFAR. It is possible that you may hit memory limits. If so (and only if you have issues), an alternative is to use ozSTAR. Note that, in this case, you will use the allocation we have been provided (which is not limitless), so make sure you use this wisely. 😀

In theory, on ozSTAR everything should be already set-up for you under /fred/oz084/mauve 

You would simply need to copy of the relevant cube and create the configuration file you plan to use. Then, you simply need to create a slurm file to launch you batch job. Here an example (with highlighted in blue parts to change:

\-----------------

 			 #!/bin/bash

 	#

  \#SBATCH --account=oz084

​      \#SBATCH --job-name= name

​     \#

  \#SBATCH --ntask= number of CPUs to use  

​     \#SBATCH --time=00:30:00

​     \#SBATCH --mem=60GB              

  module load python-scientific/3.13.1-foss-2025a

​     . /fred/oz084/mauve/py_ven/bin/activate

  ngistPipeline --config /fred/oz084/mauve/configFiles/NGC4388_MAUVE_emiles_short_sn100_constant_dust.yaml --default-dir /home/lcortese/oz084/mauve/configFiles/NGC4383.defaultDir 

 \------------------

Then, you simply run the slurm file as sbatch slurmfile 

This should be your last resort 😉