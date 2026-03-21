import argparse
import sys

import astropy.io.fits as fits
import numpy as np
import yaml

#################################################
#### Version 1 - Created by L. Cortese - Jan 12, 2024
####
#### last update: Sep 30, 2024
#################################################


###########################################################
#### SET-UP INFO ####
#########################################################
### Name of gist Master config file
### This is the only file where modifications of set-up should be made.
### If you change the file, please make sure you change its name to reflect any changes you made
master_config = "MAUVE_MasterConfig_v7.2_sn40.yaml"
##
### MAUVE sample info file containing values of redshift, ebv and sigma to use
mauve_info_file = "GIST_setupinput_v1.fits"
##
# MUSE cube subscript
cube_sub = "_DATACUBE_FINAL_WCS_Pall_mad_red_v2.fits"
##
########################################################
########################################################


parser = argparse.ArgumentParser(
    description="Creates GIST yaml configuration file for MAUVE",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument("galid", help="MAUVE GALAXY ID")
parser.add_argument("-cpu", default=16, help="cpu used")  ## default CPU=16
parser.add_argument(
    "-center", default="219,219", help="coordinates of the center of cube in pixels"
)  ## default cube centre 219,219
args = parser.parse_args()
setup = vars(args)
print(setup)

########################################################
########################################################


###
### gets galaxy name when lanching
galid = sys.argv[1]

###################################
### read MasterConfig YAML file
###################################
with open(master_config, "r") as f:
    configfile = yaml.safe_load(f)


######################################
### read MAUVE GIST set-up input file
######################################
mauve_sample = fits.open(mauve_info_file)

### extract z and (if id is correct) also ebv and sigma
z = mauve_sample[1].data["z"][np.where(mauve_sample[1].data["Galaxy"] == galid)]

### check id is correct
if z.size == 0:
    print("--ERROR--")
    print(
        "No MAUVE galaxy found with id",
        galid,
        ". Please check that you have entered the correct galaxy id",
    )
    sys.exit()

ebv = mauve_sample[1].data["EBV"][np.where(mauve_sample[1].data["Galaxy"] == galid)]


sigma = mauve_sample[1].data["sigma"][np.where(mauve_sample[1].data["Galaxy"] == galid)]

######################################
### Modify the YAML config file
######################################

configfile["GENERAL"]["RUN_ID"] = galid
configfile["GENERAL"]["INPUT"] = galid + cube_sub
configfile["GENERAL"]["REDSHIFT"] = round(z.item(), 6)
configfile["GENERAL"]["NCPU"] = int(setup["cpu"])

configfile["READ_DATA"]["ORIGIN"] = setup["center"]
configfile["READ_DATA"]["EBmV"] = round(ebv.item(), 6)

configfile["SPATIAL_MASKING"]["MASK"] = galid + "_mask.fits"

configfile["KIN"]["SIGMA"] = round(sigma.item(), 0)

configfile["CONT"]["SIGMA"] = round(sigma.item(), 0)


### write output
outname = galid + "_" + master_config
with open(outname, "w") as outfile:
    yaml.dump(configfile, outfile, default_flow_style=False, sort_keys=False)
