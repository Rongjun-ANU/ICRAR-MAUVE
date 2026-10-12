from pathlib import Path
import json, hashlib, shutil, re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path('/Users/Igniz/Desktop/ICRAR/further')
MAUVE=ROOT.parent/'MAUVE'
TMP=Path('/private/tmp/mauve_second_paper_20261011')
ASSET=MAUVE/'assets/20261011_paper_II_synthesis'
ASSET.mkdir(exist_ok=True)
REPORT=MAUVE/'20261011_MAUVE_Paper_II_Findings_Synthesis_and_Paper_Plan.md'
rows=json.loads((TMP/'notebook_inventory.json').read_text())
families={
 'Reduction, fitting and calibration provenance':list(range(1,29))+[30,33,34,38,43,50,51],
 'Metallicity branch (science excluded from this report)':[29,31,32,35,36,37,39,41,42,44,46,47,48,49,52,65,79,80,81,82,83,84,100],
 'Early structural/WHAN diagnostics':[40,45],
 'Spiral/KTZ branch (separate exploration)':[53,54,55,56,127],
 'Stage rSFMS and structural controls':[57,58,59,60,61,62,63,64,66,73,74,75,76,77,78],
 'HI pressure and stage-context experiments':list(range(67,73)),
 'SF occupancy and survivor distributions':[85,86,87,91,92,94,95,96,97,98,99,102,103,104,105],
 'Spatial asymmetry and directional comparisons':[88,89,90,93,106,126],
 'Line fading, source mixtures and excitation':list(range(108,126)),
 'Presentation, computing and scratch':[101,107,128],
}
flat=[n for ns in families.values() for n in ns]
assert sorted(flat)==list(range(1,129)) and len(set(flat))==128
parts=[]
for title,nums in families.items():
 parts.append('### '+title+'\n')
 for num in nums:
  r=rows[num-1];r['family']=title
  path=ROOT/r['name']; assert hashlib.sha256(path.read_bytes()).hexdigest()==r['sha256']
  status=f'E/C {r["executed"]}/{r["code"]}; {r["images"]} saved PNG outputs'
  if r['errors']:status+='; stored exception: '+str(r['errors'][0][1])
  parts.append(f'- **{r["id"]}** [{r["name"]}]({path}) - {status}.')
 parts.append('')
reports=[
('R01','20260731 rSFMS and rMZR by infall stages.md','First stage-rSFMS impression; line 5.'),
('R02','20260803 Pre peak rSFMS.md','Selection dependence of the PHANGS comparison; lines 3-15.'),
('R03','20260827 Resolved Star Formation Population through RPS.md','Occupancy and survivor-intensity distinction; lines 3-35.'),
('R04','20260828 SF fraction in Pre-peak.md','Individual-galaxy diversity; lines 3-15.'),
('R05','20260906 Where Star Formation Remains along Infall Stages.md','Strict-S/N spatial category pattern; lines 3-21.'),
('R06','20260831 Pre-Peak rSFMS with SF, ND and NSF.md','ND proxy and category interpretation; lines 3-9; compare current N092.'),
('R07','20260901 Close-to-peak and Post-peak rSFMS with SF, NSF and ND.md','Historical galaxy-by-galaxy interpretation; later quantitative work qualifies it.'),
('R08','20260909 Remained SF and Increasing NSF along Infall Stages.md','Early synthesis, including wording later corrected; lines 21-59.'),
('R09','20260911_Physical_Scenario_for_Remaining_SF_and_Increasing_NSF.md','Source audit and two-stage resampling; sections near lines 115-177 and 568-578.'),
('R10','20260914_Resolved_Gas_Loss_Star_Formation_and_Ionization_Model.md','Key empirical extraction and attenuation fit; sections 1.3-1.5 and 6.3-6.4; lines 87-165 and 1043-1076.'),
('R11','20260930_HI_Stripping_and_HII_DIG_Line_Ratio_Model.md','Why proportionally fading young components cannot change their mixture; lines 472-498.'),
('R12','20261002_Model0_Gas_Fading_and_HOLMES_Derivation.md','Photon normalization and quantitative model tensions; section 8, especially lines 624-853.'),
('R13','20261003 The Origin of Fading Emission Lines and Increasing BPT Line Ratios along Infall Stages.md','Simplified analytical narrative; model time is not the observed stage clock.'),
('R14','20261008 Distribution of SF NSF TNSF BPT Distances Across Infall Stages.md','Additive PDF interpretation; complement shorthand must be checked against the current dual-BPT mask.'),
('R15','20260812 Comprehensive Literature Atlas and Environmental Histories of the MAUVE Virgo Galaxy Sample.md','Environmental-history context; 40-galaxy scope, not the current optical sample.'),
('R16','20260623 IC3392 7000 vs 8900.md','Wavelength-range sensitivity; numerical table at lines 16-74.'),
('R17','20260624 NGC4064 7000 vs 8900.md','Second range-comparison example; numerical table at lines 16-74.'),
('R18','20261006 Halpha Dispersion along Infall Stages.md','Support-sensitive linewidth interpretation; historical TNSF shorthand superseded by current source.'),
('R19','20261001_Connected_HI_Stripping_and_HOLMES_Line_Ratio_Model.md','Intermediate population-normalized model and retained limitations.'),
('R20','20260830_RPS_Local_Star_Formation_Theory_and_Toy_Models.md','Local enhancement theoretical possibilities; not an observed triggering measurement.'),
('R21','20260831_Molecular_Gas_Supply_Definitions_and_Observational_Inference.md','Distinguishes reservoir supply, HI-to-H2 conversion, and cloud accretion.'),
('R22','20260518 SFH and Stellar Mass Check.md','Mass sensitivity in tests that also vary range and velocity scale; not a pure 7000/8900 correction.'),
]
rp=[];report_hashes=[]
for rid,name,note in reports:
 path=MAUVE/name
 assert path.exists(),str(path)
 rp.append(f'- **{rid}** [{name}]({str(path).replace(" ","%20")}) - {note}')
 report_hashes.append(dict(id=rid,path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
text=REPORT.read_text()
text=text.replace('<!-- NOTEBOOK_INVENTORY -->','\n'.join(parts))
text=text.replace('<!-- REPORT_REGISTER -->','\n\n'.join(rp))
REPORT.write_text(text)
(ASSET/'notebook_inventory.json').write_text(json.dumps(rows,indent=2))
(ASSET/'report_sources.json').write_text(json.dumps(report_hashes,indent=2))
for file in ['notebook_inventory.txt','report_index.md','report_inventory.tsv','check_central_results.py','live_check.log','central_result_check.json']:
 shutil.copy2(TMP/file,ASSET/file)
for file in TMP.glob('live_*.csv'):shutil.copy2(file,ASSET/file.name)
selected=['N087_c19_o1_t1.csv','N087_c19_o2_t1.csv','N103_c13_o2_t1.csv','N104_c15_o4_t1.csv','N105_c31_o1_t1.csv','N109_c43_o3_t1.csv','N109_c37_o1_t1.csv','N120_c21_o1_t1.csv','N121_c21_o1_t1.csv','N124_c13_o2_t1.csv','N126_c9_o2_t1.csv','N126_c16_o2_t1.csv','N126_c14_o2_t1.csv','N126_c19_o1_t1.csv','N126_c19_o2_t1.csv']
for file in selected:shutil.copy2(TMP/file,ASSET/file)

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':200})
stages=['pre-peak','close-to-peak','post-peak']
colors=['#2166ac','#1a9850','#d73027']
cats=['SF','ND','NSF'];catcolors=['#2166ac','#e69f00','#b2182b']
fra=pd.read_csv(ASSET/'live_category_fractions_by_stage.csv').set_index('stage').loc[stages]
d=pd.read_csv(ASSET/'live_sfr_decomposition.csv')
fig,axes=plt.subplots(1,3,figsize=(13.5,4.8),gridspec_kw={'width_ratios':[1,1.15,1.05]})
fig.subplots_adjust(left=.065,right=.99,bottom=.30,top=.89,wspace=.30)
left=np.zeros(3)
for cat,color in zip(cats,catcolors):
 vals=fra['f_'+cat].to_numpy()
 axes[0].barh(range(3),vals,left=left,color=color,label=cat,height=.55)
 for i,v in enumerate(vals):axes[0].text(left[i]+v/2,i,f'{v:.2f}',ha='center',va='center',fontsize=9,color='white' if cat!='ND' else 'black')
 left+=vals
axes[0].set(yticks=range(3),yticklabels=['Pre (8)','Close (5)','Post (13)'],xlim=(0,1),xlabel='Mean usable-footprint fraction',title='A. Full-footprint category fractions')
axes[0].invert_yaxis()
fig.legend(*axes[0].get_legend_handles_labels(),loc='lower left',bbox_to_anchor=(.065,.025),ncol=3,frameon=False)
terms=['occupancy','survivor_intensity','total_hii']; labels=['Occupancy','Survivor intensity','HII-traced total'];cs=['#7b3294','#e69f00','#222222']
for ax,variable,centres in [(axes[1],'log_sigma_star',[7.625,8.125,8.625]),(axes[2],'radius_re',[.625,1.625])]:
 sel=d[d.stage.eq('post-peak')&d.variable.eq(variable)].copy();sel['center']=(sel.bin_left+sel.bin_right)/2;sel=sel.set_index('center').loc[centres]
 xx=np.arange(len(centres));width=.23
 for k,(term,label,color) in enumerate(zip(terms,labels,cs)):
  ax.bar(xx+(k-1)*width,sel['delta_log_'+term],width,color=color,label=label)
 ax.axhline(0,color='.3',lw=.8);ax.set(xticks=xx,xticklabels=[str(x) for x in centres],ylim=(-1.7,.3),ylabel='Post / pre contrast [dex]')
 ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
axes[1].set(xlabel=r'$\log_{10}\Sigma_\star$ bin centre [$M_\odot\,\mathrm{kpc}^{-2}$]',title='B. At fixed stellar density')
axes[2].set(xlabel=r'$R/R_e$ bin centre',title='C. At fixed radius')
fig.legend(*axes[1].get_legend_handles_labels(),loc='lower left',bbox_to_anchor=(.40,.015),ncol=1,frameon=False,fontsize=9)
fig.savefig(ASSET/'figure_01_central_findings.png',bbox_inches='tight');plt.close(fig)

fig,axes=plt.subplots(1,2,figsize=(11,4.7),layout='constrained')
q=pd.read_csv(ASSET/'N087_c19_o1_t1.csv');q=q[q.metric.eq('M')].set_index('stage').loc[stages]
for i,(s,c) in enumerate(zip(stages,colors)):
 row=q.loc[s];axes[0].errorbar(i,row['median'],yerr=[[row['median']-row.boot_p2p5],[row.boot_p97p5-row['median']]],fmt='o',color=c,capsize=5,ms=7)
axes[0].axhline(0,color='.5',ls='--',lw=1);axes[0].set(xticks=range(3),xticklabels=['Pre','Close','Post'],ylabel='Product median survivor residual [dex]',title='A. August Q3: saved 95% intervals')
ng=pd.read_csv(ASSET/'N126_c16_o2_t1.csv')
for i,row in ng.iterrows():
 axes[1].errorbar(row.estimate_dex,i,xerr=[[row.estimate_dex-row.CI95_low],[row.CI95_high-row.estimate_dex]],fmt='o',color='#344e64',capsize=5,ms=7)
axes[1].axvline(0,color='.5',ls='--',lw=1);axes[1].set(yticks=range(3),yticklabels=['Positive side','Negative side','Side difference'],xlabel='Held-out external residual [dex]',title='B. NGC4654: saved 95% intervals');axes[1].invert_yaxis()
for ax in axes:ax.grid(alpha=.15);ax.set_axisbelow(True)
fig.savefig(ASSET/'figure_02_survivors_and_asymmetry.png',bbox_inches='tight');plt.close(fig)

ledger=[
('occupancy','N104 cells 15,19,24; N103 cells 3,5,7,9,11','live central reproduction','live_category_fractions_by_stage.csv'),
('SFR decomposition','N103 cells 9,11,13','live central reproduction; stored uncertainty','live_sfr_decomposition.csv'),
('survivor distribution','N087 cells 19,34,35','saved notebook output','N087_c19_o1_t1.csv'),
('NSF intensity','N105 cell 31; R10 section 1.3','saved notebook/report result','N105_c31_o1_t1.csv'),
('HII invariance','N109 cell 43','saved notebook output','N109_c43_o3_t1.csv'),
('BPT distance','N120 and N121 cell 21','saved notebook output','N120_c21_o1_t1.csv'),
('fixed Halpha ratios','N124 cells 11,13','saved notebook output','N124_c13_o2_t1.csv'),
('NGC4654','N126 cells 9,14,16,19,21,23','saved notebook output','N126_c16_o2_t1.csv'),
('gas/HOLMES models','R12 section 8; N110 cells 1,26','conditional model demonstration',''),
]
pd.DataFrame(ledger,columns=['claim','source','status','evidence_file']).to_csv(ASSET/'claim_evidence_ledger.csv',index=False)
shutil.copy2(__file__,ASSET/'prepare_report_assets.py')
print('ASSETS_READY',len(rows),'notebooks',len(reports),'report sources; words',len(REPORT.read_text().split()))
