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

vcp -v arc:projects/mauve/extended /Users/Igniz/Desktop/ICRAR/