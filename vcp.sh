# download canfar to local
# from arc:home/RongjunHuang/ICRAR/extended
# to local: /Users/Igniz/Desktop/ICRAR/extended
# transfer files: *gas*extended.fits, *BINNING*extended.fits, 
# transfer folders: mass_logs, sfr_logs

# vcp -v -i --overwrite \
#     --include '*gas*extended.fits' \
#     --include '*BINNING*extended.fits' \
#     --include 'mass_logs/***' \
#     --include 'sfr_logs/***' \
#     --exclude '*' \
#     arc:home/RongjunHuang/ICRAR/extended /Users/Igniz/Desktop/ICRAR/extended

# vcp arc:home/RongjunHuang/ICRAR/extended/sfr_logs /Users/Igniz/Desktop/ICRAR/extended/
# vcp arc:home/RongjunHuang/ICRAR/extended/mass_logs /Users/Igniz/Desktop/ICRAR/extended/
# vcp arc:home/RongjunHuang/ICRAR/extended/*gas*extended.fits /Users/Igniz/Desktop/ICRAR/extended/
# vcp arc:home/RongjunHuang/ICRAR/extended/*BINNING*extended.fits /Users/Igniz/Desktop/ICRAR/extended/

cadc-get-cert -u RongjunHuang

vcp -v arc:projects/mauve/extended /Users/Igniz/Desktop/ICRAR/
vcp -v arc:projects/mauve/alma/products/v0.3/ /Users/Igniz/Desktop/ICRAR/CO/products_v0.3/
vcp -v arc:projects/mauve/nGIST_dev_7.6.4/ /Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4/
vcp -v arc:projects/mauve/nGIST_dev_7.6.4/IC3392_9300_KIN /Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4/IC3392_9300_KIN
vcp -v arc:projects/mauve/v3tk_to_R/ /Users/Igniz/Desktop/ICRAR/v3tk_to_R/
vcp -v arc:projects/mauve/v3tk_to_VRI/ /Users/Igniz/Desktop/ICRAR/v3tk_to_VRI/
vcp -v arc:projects/mauve/nGIST_dev_7.6.4/IC3392_v3tk_mask_uptoCONT/IC3392_v3tk_mask_uptoCONT_line_cube.fits /Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4/IC3392_v3tk_mask_uptoCONT/IC3392_v3tk_mask_uptoCONT_line_cube.fits
vcp -v arc:projects/mauve/nGIST_dev_7.6.4/IC3392_v3tk_mask_uptoCONT_MDEG12/IC3392_v3tk_mask_uptoCONT_MDEG12_line_cube.fits /Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4/IC3392_v3tk_mask_uptoCONT_MDEG12/IC3392_v3tk_mask_uptoCONT_MDEG12_line_cube.fits
vcp -v arc:projects/mauve/nGIST_dev_7.6.4/IC3392_v3tk_mask_uptoCONT_MDEG16/IC3392_v3tk_mask_uptoCONT_MDEG16_line_cube.fits /Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4/IC3392_v3tk_mask_uptoCONT_MDEG16/IC3392_v3tk_mask_uptoCONT_MDEG16_line_cube.fits
vcp -v arc:projects/mauve/nGIST_dev_7.6.4/IC3392_v3tk_mask_uptoCONT_MDEG20/IC3392_v3tk_mask_uptoCONT_MDEG20_line_cube.fits /Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4/IC3392_v3tk_mask_uptoCONT_MDEG20/IC3392_v3tk_mask_uptoCONT_MDEG20_line_cube.fits
vcp -v arc:projects/mauve/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG8/IC3392_v3tk_9100_mask_uptoCONT_MDEG8_line_cube.fits /Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG8/IC3392_v3tk_9100_mask_uptoCONT_MDEG8_line_cube.fits
vcp -v arc:projects/mauve/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG12/IC3392_v3tk_9100_mask_uptoCONT_MDEG12_line_cube.fits /Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG12/IC3392_v3tk_9100_mask_uptoCONT_MDEG12_line_cube.fits
vcp -v arc:projects/mauve/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG16/IC3392_v3tk_9100_mask_uptoCONT_MDEG16_line_cube.fits /Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG16/IC3392_v3tk_9100_mask_uptoCONT_MDEG16_line_cube.fits
vcp -v arc:projects/mauve/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG20/IC3392_v3tk_9100_mask_uptoCONT_MDEG20_line_cube.fits /Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG20/IC3392_v3tk_9100_mask_uptoCONT_MDEG20_line_cube.fits
vcp -v arc:projects/mauve/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG14/IC3392_v3tk_9100_mask_uptoCONT_MDEG14_line_cube.fits /Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG14/IC3392_v3tk_9100_mask_uptoCONT_MDEG14_line_cube.fits
vcp -v arc:projects/mauve/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG23/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_line_cube.fits /Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG23/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_line_cube.fits
vcp -v arc:projects/mauve/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG8/IC3392_v3tk_9100_mask_uptoCONT_MDEG8_kin_maps.fits /Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG8/IC3392_v3tk_9100_mask_uptoCONT_MDEG8_kin_maps.fits
vcp -v arc:projects/mauve/nGIST_dev_7.6.4/IC3392_v3tk_mask_uptoCONT/IC3392_v3tk_mask_uptoCONT_kin_maps.fits /Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4/IC3392_v3tk_mask_uptoCONT/IC3392_v3tk_mask_uptoCONT_kin_maps.fits
vcp -v arc:projects/mauve/nGIST_dev_7.6.4/IC3392_v3tk_mask_EMILES_uptoCONT/IC3392_v3tk_mask_EMILES_uptoCONT_line_cube.fits /Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4/IC3392_v3tk_mask_EMILES_uptoCONT/IC3392_v3tk_mask_EMILES_uptoCONT_line_cube.fits
vcp -v arc:projects/mauve/nGIST_dev_7.6.4/IC3392_v3tk_mask_EMILES_uptoCONT/IC3392_v3tk_mask_EMILES_uptoCONT_kin_maps.fits /Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4/IC3392_v3tk_mask_EMILES_uptoCONT/IC3392_v3tk_mask_EMILES_uptoCONT_kin_maps.fits
vcp -v arc:projects/mauve/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG50/IC3392_v3tk_9100_mask_uptoCONT_MDEG50_line_cube.fits /Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG50/IC3392_v3tk_9100_mask_uptoCONT_MDEG50_line_cube.fits
vcp -v arc:projects/mauve/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG100/IC3392_v3tk_9100_mask_uptoCONT_MDEG100_line_cube.fits /Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG100/IC3392_v3tk_9100_mask_uptoCONT_MDEG100_line_cube.fits
vcp -v arc:projects/mauve/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG30/IC3392_v3tk_9100_mask_uptoCONT_MDEG30_line_cube.fits /Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG30/IC3392_v3tk_9100_mask_uptoCONT_MDEG30_line_cube.fits
vcp -v arc:projects/mauve/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG40/IC3392_v3tk_9100_mask_uptoCONT_MDEG40_line_cube.fits /Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG40/IC3392_v3tk_9100_mask_uptoCONT_MDEG40_line_cube.fits
vcp -v arc:projects/mauve/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG0/IC3392_v3tk_9100_mask_uptoCONT_MDEG0_line_cube.fits /Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG0/IC3392_v3tk_9100_mask_uptoCONT_MDEG0_line_cube.fits
vcp -v arc:projects/mauve/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_sky20/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_sky20_line_cube.fits /Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_sky20/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_sky20_line_cube.fits
vcp -v arc:projects/mauve/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_sky20_no7000/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_sky20_no7000_line_cube.fits /Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_sky20_no7000/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_sky20_no7000_line_cube.fits
vcp -v arc:projects/mauve/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG12_sky20_no7000/IC3392_v3tk_9100_mask_uptoCONT_MDEG12_sky20_no7000_line_cube.fits /Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG12_sky20_no7000/IC3392_v3tk_9100_mask_uptoCONT_MDEG12_sky20_no7000_line_cube.fits
vcp -v arc:projects/mauve/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_sky20_no7000_test1/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_sky20_no7000_test1_line_cube.fits /Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_sky20_no7000_test1/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_sky20_no7000_test1_line_cube.fits
vcp -v arc:projects/mauve/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_sky20_no7000_test2/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_sky20_no7000_test2_line_cube.fits /Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_sky20_no7000_test2/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_sky20_no7000_test2_line_cube.fits
vcp -v arc:projects/mauve/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_sky20_no7000_test3/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_sky20_no7000_test3_line_cube.fits /Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_sky20_no7000_test3/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_sky20_no7000_test3_line_cube.fits
vcp -v arc:projects/mauve/nGIST_dev_7.6.4/IC3392_v3tk_mask_uptoCONT_allMILES/IC3392_v3tk_mask_uptoCONT_allMILES_line_cube.fits /Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4/IC3392_v3tk_mask_uptoCONT_allMILES/IC3392_v3tk_mask_uptoCONT_allMILES_line_cube.fits
vcp -v arc:projects/mauve/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_narrow6/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_narrow6_line_cube.fits /Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_narrow6/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_narrow6_line_cube.fits
vcp -v arc:projects/mauve/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_narrow7/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_narrow7_line_cube.fits /Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_narrow7/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_narrow7_line_cube.fits
vcp -v arc:projects/mauve/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_narrow8/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_narrow8_line_cube.fits /Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_narrow8/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_narrow8_line_cube.fits
vcp -v arc:projects/mauve/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_narrow10/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_narrow10_line_cube.fits /Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_narrow10/IC3392_v3tk_9100_mask_uptoCONT_MDEG23_narrow10_line_cube.fits
vcp -v arc:projects/mauve/products/v0.6/IC3392/ /Users/Igniz/Desktop/ICRAR/further/IC3392_v3_v0.6

# Download all v0.7dev_gas_BIN_maps.fits from XXX_v0.7dev folders
GALAXIES=("IC3392" "NGC4064" "NGC4192" "NGC4293" "NGC4298" "NGC4330" "NGC4383" "NGC4396" "NGC4419" "NGC4457" "NGC4501" "NGC4522" "NGC4694" "NGC4698")
for gal in "${GALAXIES[@]}"; do
    vcp -v "arc:projects/mauve/toby_sandbox/products/${gal}_v0.7dev/${gal}_v0.7dev_gas_BIN_maps.fits" /Users/Igniz/Desktop/ICRAR/extended/ &
done
wait

vcp -v arc:/projects/mauve/products/v0.6/IC3392/IC3392_CONTcube.fits /Users/Igniz/Desktop/ICRAR/extended/IC3392_CONTcube.fits
