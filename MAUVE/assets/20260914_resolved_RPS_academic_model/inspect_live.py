import json, pathlib, os, contextlib, hashlib, ast, time
import numpy as np
root=pathlib.Path("/Users/Igniz/Desktop/ICRAR/further")
out=pathlib.Path("/Users/Igniz/Desktop/ICRAR/MAUVE/assets/20260914_resolved_RPS_academic_model")
os.chdir(root)
p=root/"20260909_check_Sigma_SFR_intensity_by_stage.ipynb"
nb=json.loads(p.read_text()); ns={}
with (out/"live_sfr_execution.log").open("w") as log, contextlib.redirect_stdout(log):
 for cell in [3,5,7,9,11]:
  exec(compile("".join(nb["cells"][cell-1]["source"]),f"{p.name}:cell{cell}","exec"),ns)
 for key in ["STAGE_SFR_PROFILES","SFR_BY_GALAXY_BIN","SAMPLE_QC","SNR_RETAINED_SUMMARY"]:
  if key in ns: ns[key].to_csv(out/(key.lower()+".csv"),index=False)
print("SFR_SELECTED_CELLS_PASS", ns["sample"].groupby("stage").size().to_dict(), flush=True)
# Inspect Balmer-decrement invariance on the same retained SF/NSF pixels, using live maps.
records=[]
from astropy.io import fits
for row in ns["sample"].itertuples(index=False):
 maps=ns["load_galaxy_maps"](row,ns["geometry"])
 keep=maps["valid_disc_mask"] & np.isfinite(maps["snr_postfit"]) & (maps["snr_postfit"]>25)
 with fits.open(row.sfr_path,memmap=True) as hd:
  raw=hd["HA6562_FLUX"].data/hd["HB4861_FLUX"].data
  corr=hd["HA6562_FLUX_CORR"].data/hd["HB4861_FLUX_CORR"].data
  for cat in ["sf","nsf"]:
   use=keep & maps[cat+"_mask"] & np.isfinite(corr)
   records.append({"GALID":row.GALID,"stage":row.stage,"category":cat,"n":int(use.sum()),"n_raw_below_286":int((use & (raw<2.86)).sum()),"corr_min":float(np.min(corr[use])) if use.any() else None,"corr_max":float(np.max(corr[use])) if use.any() else None,"CUT_SN":hd[0].header["CUT_SN"],"NOISE":hd[0].header["NOISE"],"SFRCOEF":hd[0].header["SFRCOEF"]})
 del maps
(out/"balmer_live_check.json").write_text(json.dumps(records,indent=2))
print("BALMER_LIVE_CHECK_PASS",flush=True)
