#!/usr/bin/env python
import argparse
import csv
import os
import sys

import astropy.io.fits as fits
import numpy as np
import yaml

#################################################
#### Version 1.1 - Created by L. Cortese - Jan 12, 2024
####
#### try version for local validation - Mar 20, 2026
#################################################


###########################################################
#### SET-UP INFO ####
#########################################################
### Name of gist Master config file
### This is the only file where modifications of set-up should be made.
### If you change the file, please make sure you change its name to reflect any changes you made
master_config = "MAUVE_MasterConfig_v7.6.8_setonix.yaml"
##
### MAUVE sample info file containing values of redshift, ebv and sigma to use
### Supports either a real FITS file or a text pointer to a FITS file path.
mauve_info_file = "GIST_setupinput_v1.fits"
##
### CSV file containing cube centers (ID, center_y, center_x)
cube_centers_file = "cube_centers_v3tk.csv"
##
### Combined cube IDs that are represented by separate galaxies in the setup table.
### The generator keeps the combined cube ID for RUN_ID/INPUT/MASK, but averages
### z, EBV, and sigma from these setup-table rows.
combined_setup_ids = {
    "NGC4567_8": ("NGC4567", "NGC4568"),
}
##
# MUSE cube subscript
cube_sub = "_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk.fits"
# path
pathcube = '/scratch/pawsey1308/mauve/cubes/v3tk/'
pathproducts = '/scratch/pawsey1308/mauve/products/v3tk_v7.6.8/'

##
########################################################
########################################################


parser = argparse.ArgumentParser(
    description="Creates GIST yaml configuration file for MAUVE",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument("galid", help="MAUVE GALAXY ID")
parser.add_argument("-cpu", default=128, help="cpu used")
parser.add_argument(
    "-center",
    default=None,
    help="optional manual center override as 'x,y' (otherwise read from cube_centers_v3tk.csv)",
)
args = parser.parse_args()
setup = vars(args)
print(setup)


def resolve_mauve_info_file(path):
    # If this is a text stub with a FITS path inside, resolve it relative to this script.
    with open(path, "rb") as f:
        header = f.read(6)

    if header == b"SIMPLE":
        return path

    with open(path, "r", encoding="utf-8") as f:
        target = f.read().strip()

    if not target:
        raise ValueError(f"{path} is not a valid FITS file and does not contain a target path")

    resolved = os.path.join(os.path.dirname(path), target)
    if not os.path.exists(resolved):
        raise FileNotFoundError(
            f"Resolved MAUVE setup file does not exist: {resolved} (from {path})"
        )

    return resolved


def resolve_cube_centers_file(path):
    # Resolve the centers CSV relative to this script.
    resolved = os.path.join(os.path.dirname(__file__), path)
    if not os.path.exists(resolved):
        raise FileNotFoundError(f"Cube centers CSV not found: {resolved}")
    return resolved


def get_center_from_csv(galaxy_id, csv_path):
    with open(csv_path, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("ID", "").strip() == galaxy_id:
                center_y = row.get("center_y")
                center_x = row.get("center_x")
                if center_y is None or center_x is None:
                    raise ValueError(
                        f"Missing center_y/center_x for galaxy {galaxy_id} in {csv_path}"
                    )
                return f"{int(center_x)},{int(center_y)}"

    raise ValueError(f"Galaxy ID {galaxy_id} not found in cube centers CSV: {csv_path}")


def get_setup_mask(sample_data, galaxy_id):
    galaxy_names = np.char.strip(sample_data["Galaxy"].astype(str))
    mask = galaxy_names == galaxy_id
    if np.any(mask):
        return mask, (galaxy_id,)

    setup_ids = combined_setup_ids.get(galaxy_id)
    if setup_ids is None:
        return mask, ()

    mask = np.isin(galaxy_names, setup_ids)
    found = set(galaxy_names[mask])
    missing = [setup_id for setup_id in setup_ids if setup_id not in found]
    if missing:
        raise ValueError(
            f"Combined cube {galaxy_id} is missing setup rows: {', '.join(missing)}"
        )

    print(
        f"Using combined setup rows for {galaxy_id}: "
        f"{', '.join(setup_ids)}"
    )
    return mask, setup_ids

########################################################
########################################################


###
### gets galaxy name when lanching
galid = setup["galid"]

###################################
### read MasterConfig YAML file
###################################
with open(master_config, "r") as f:
    configfile = yaml.safe_load(f)


######################################
### read MAUVE GIST set-up input file
######################################
mauve_info_resolved = resolve_mauve_info_file(mauve_info_file)
print(f"Using MAUVE setup file: {mauve_info_resolved}")
mauve_sample = fits.open(mauve_info_resolved)

if setup["center"] is None:
    centers_csv_resolved = resolve_cube_centers_file(cube_centers_file)
    setup["center"] = get_center_from_csv(galid, centers_csv_resolved)
    print(f"Using center from CSV: {setup['center']} ({centers_csv_resolved})")
else:
    print(f"Using manually provided center: {setup['center']}")

### extract z and (if id is correct) also ebv and sigma
setup_mask, setup_ids = get_setup_mask(mauve_sample[1].data, galid)

### check id is correct
if not np.any(setup_mask):
    print("--ERROR--")
    print(
        "No MAUVE galaxy found with id",
        galid,
        ". Please check that you have entered the correct galaxy id",
    )
    sys.exit()

z = mauve_sample[1].data["z"][setup_mask]
ebv = mauve_sample[1].data["EBV"][setup_mask]
sigma = mauve_sample[1].data["sigma"][setup_mask]
galaxies_matched = np.char.strip(mauve_sample[1].data["Galaxy"].astype(str))[setup_mask]

print("Recorded values in FITS file:")
for g_id, g_z, g_ebv, g_sigma in zip(galaxies_matched, z, ebv, sigma):
    print(f"  {g_id}: z={g_z:.6f}, EBV={g_ebv:.6f}, sigma={g_sigma:.1f}")

z_value = float(np.mean(z))
ebv_value = float(np.mean(ebv))
sigma_value = float(np.mean(sigma))

print("Final values used in YAML:")
print(f"  z={z_value:.6f}, EBV={ebv_value:.6f}, sigma={sigma_value:.1f}")

######################################
### Modify the YAML config file
######################################

configfile["GENERAL"]["RUN_ID"] = galid
configfile["GENERAL"]["INPUT"] = pathcube + galid + cube_sub
configfile["GENERAL"]["OUTPUT"] = pathproducts
configfile["GENERAL"]["REDSHIFT"] = round(z_value, 6)
configfile["GENERAL"]["NCPU"] = int(setup["cpu"])

configfile["READ_DATA"]["ORIGIN"] = setup["center"]
configfile["READ_DATA"]["EBmV"] = round(ebv_value, 6)

configfile["SPATIAL_MASKING"]["MASK"] = galid + "_mask.fits"

configfile["KIN"]["SIGMA"] = int(round(sigma_value))

configfile["CONT"]["SIGMA"] = int(round(sigma_value))


### write output
outname = galid + "_" + master_config
with open(outname, "w") as outfile:
    yaml.dump(configfile, outfile, default_flow_style=False, sort_keys=False)
