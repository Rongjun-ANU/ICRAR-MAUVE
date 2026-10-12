"""Rebuild the report figures from the companion evidence tables only."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ASSET = Path(__file__).resolve().parent
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
                     'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':200})
stages=['pre-peak','close-to-peak','post-peak']
stage_colors=['#2166ac','#1a9850','#d73027']
fra=pd.read_csv(ASSET/'live_category_fractions_by_stage.csv').set_index('stage').loc[stages]
d=pd.read_csv(ASSET/'live_sfr_decomposition.csv')
fig=plt.figure(figsize=(7.2,8.0))
gs=fig.add_gridspec(2,2,height_ratios=[.8,1])
axes=[fig.add_subplot(gs[0,:]),fig.add_subplot(gs[1,0]),fig.add_subplot(gs[1,1])]
fig.subplots_adjust(left=.12,right=.975,top=.94,bottom=.21,hspace=.80,wspace=.48)
left=np.zeros(3)
for cat,color in zip(['SF','ND','NSF'],['#2166ac','#e69f00','#b2182b']):
    vals=fra['f_'+cat].to_numpy()
    axes[0].barh(range(3),vals,left=left,color=color,label=cat,height=.55)
    for i,v in enumerate(vals):
        axes[0].text(left[i]+v/2,i,f'{v:.2f}',ha='center',va='center',fontsize=10,
                     color='white' if cat!='ND' else 'black')
    left+=vals
axes[0].set(yticks=range(3),yticklabels=['Pre (8)','Close (5)','Post (13)'],xlim=(0,1),
            xlabel='Mean usable-footprint fraction',title='A. Full-footprint category fractions')
axes[0].invert_yaxis()
axes[0].legend(loc='upper center',bbox_to_anchor=(.5,-.38),ncol=3,frameon=False)
terms=['occupancy','survivor_intensity','total_hii']
labels=['Occupancy','Survivor intensity','HII-traced total']
term_colors=['#7b3294','#e69f00','#222222']
for ax,variable,centres in [(axes[1],'log_sigma_star',[7.625,8.125,8.625]),
                            (axes[2],'radius_re',[.625,1.625])]:
    sel=d[d.stage.eq('post-peak')&d.variable.eq(variable)].copy()
    sel['center']=(sel.bin_left+sel.bin_right)/2
    sel=sel.set_index('center').loc[centres]
    xx=np.arange(len(centres));width=.23
    for k,(term,label,color) in enumerate(zip(terms,labels,term_colors)):
        ax.bar(xx+(k-1)*width,sel['delta_log_'+term],width,color=color,label=label)
    ax.axhline(0,color='.3',lw=.8)
    ax.set(xticks=xx,xticklabels=[str(x) for x in centres],ylim=(-1.7,.3),
           ylabel='Post / pre contrast [dex]')
    ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
axes[1].set(xlabel=r'$\log_{10}[\Sigma_\star/(M_\odot\,\mathrm{kpc}^{-2})]$'+'\nbin centre',
            title='B. Fixed stellar density')
axes[2].set(xlabel=r'$R/R_e$ bin centre',title='C. Fixed radius')
fig.legend(*axes[1].get_legend_handles_labels(),loc='lower center',
           bbox_to_anchor=(.54,.005),ncol=1,frameon=False,fontsize=10)
fig.savefig(ASSET/'figure_01_central_findings.png',bbox_inches='tight');plt.close(fig)

fig,axes=plt.subplots(2,1,figsize=(7.2,7.2),layout='constrained')
q=pd.read_csv(ASSET/'N087_c19_o1_t1.csv')
q=q[q.metric.eq('M')].set_index('stage').loc[stages]
for i,(s,c) in enumerate(zip(stages,stage_colors)):
    row=q.loc[s]
    axes[0].errorbar(i,row['median'],
                     yerr=[[row['median']-row.boot_p2p5],[row.boot_p97p5-row['median']]],
                     fmt='o',color=c,capsize=5,ms=7)
axes[0].axhline(0,color='.5',ls='--',lw=1)
axes[0].set(xticks=range(3),xticklabels=['Pre','Close','Post'],
            ylabel='Product median survivor\nresidual [dex]',
            title='A. August Q3: saved 95% intervals',xlim=(-.4,2.4))
ng=pd.read_csv(ASSET/'N126_c16_o2_t1.csv')
for i,row in ng.iterrows():
    axes[1].errorbar(row.estimate_dex,i,
                     xerr=[[row.estimate_dex-row.CI95_low],[row.CI95_high-row.estimate_dex]],
                     fmt='o',color='#344e64',capsize=5,ms=7)
axes[1].axvline(0,color='.5',ls='--',lw=1)
axes[1].set(yticks=range(3),yticklabels=['Positive side','Negative side','Side difference'],
            xlabel='Held-out external residual [dex]',
            title='B. NGC4654: saved 95% intervals',ylim=(-.4,2.4))
axes[1].invert_yaxis()
for ax in axes:
    ax.grid(alpha=.15);ax.set_axisbelow(True)
fig.savefig(ASSET/'figure_02_survivors_and_asymmetry.png',bbox_inches='tight');plt.close(fig)
print('REPORT_FIGURES_WRITTEN: live central values and saved intervals only')
