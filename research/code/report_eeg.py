from pathlib import Path
import json,numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'results'/'eeg';MAN=ROOT/'manuscript';FIG=ROOT.parent/'output'/'QSEEG'/'frame'
FIG.mkdir(parents=True,exist_ok=True)
(ROOT/'qa').mkdir(parents=True,exist_ok=True)
f=np.arange(1,61)*.5;t=np.arange(1,10)/10
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42,'axes.titlesize':10})
# Verify every derived result archive before using it.
for p in D.glob('*.npz'):
    with np.load(p) as z:
        for k in z.files:
            a=z[k]
            if np.issubdtype(a.dtype,np.number):assert np.isfinite(a).all(),(p,k)

def read(name,method='quantile'):
    with np.load(D/(name+'_'+method+'.npz')) as z:
        b=z['effect'];se=z['se'];c=float(z['critical'])
    return b,se,c

def add(ax,b,se,c,color):
    ax.axhline(0,color='#444444',ls='--',lw=.7)
    ax.fill_between(f,b-c*se,b+c*se,color=color,alpha=.16,lw=0)
    ax.plot(f,b,color=color,lw=1.5)
    ax.set_xlim(.5,30);ax.set_xticks([.5,5,10,15,20,25,30]);ax.xaxis.set_major_formatter(FormatStrFormatter('%g'));ax.grid(axis='y',alpha=.18)
fig,axs=plt.subplots(2,3,figsize=(7.1,4.4),sharex=True,sharey=True,layout='constrained')
for row,ch in enumerate([0,1]):
    b,se,c=read(f'SC_2_{ch}_main');b=b.reshape(9,60);se=se.reshape(9,60)
    for col,qi in enumerate([0,4,8]):
        add(axs[row,col],b[qi],se[qi],c,'#165B83' if ch==0 else '#A45B20')
        axs[row,col].set_title(f'({"abcdef"[row*3+col]}) '+('Fpz-Cz' if ch==0 else 'Pz-Oz')+fr', $\tau = {t[qi]:.1f}$',loc='left',fontsize=9)
        if row==1:axs[row,col].set_xlabel('Frequency (Hz)')
    axs[row,0].set_ylabel('Age coefficient\n(per decade)')
fig.savefig(FIG/'age_associations.pdf');fig.savefig(ROOT/'qa'/'age_associations.png',dpi=180);plt.close(fig)
fig,axs=plt.subplots(2,2,figsize=(7,4.7),layout='constrained')
b,se,c=read('ST_2_0_main');b=b.reshape(9,60);se=se.reshape(9,60)
for ax,qi,letter in zip(axs.ravel()[:3],[0,4,8],'abc'):
    add(ax,b[qi],se[qi],c,'#165B83');ax.set_ylim(-.23,.56);ax.set_title(fr'({letter}) Quantile spectrum, $\tau = {t[qi]:.1f}$',loc='left')
bo,so,co=read('ST_2_0_main','ordinary');add(axs[1,1],bo,so,co,'#A45B20');axs[1,1].set_title('(d) Ordinary log power',loc='left')
for ax in axs.ravel():ax.set_xlabel('Frequency (Hz)');ax.set_ylabel('Adjusted paired contrast')
fig.savefig(FIG/'paired_medication.pdf');fig.savefig(ROOT/'qa'/'paired_medication.png',dpi=180);plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(7,2.8),sharey=True,layout='constrained')
b,se,c=read('SC_2_1_main','heterogeneity');b=b.reshape(8,60);se=se.reshape(8,60)
for ax,k,label in zip(axs,[0,7],[r'(a) $\tau=0.1$ minus median',r'(b) $\tau=0.9$ minus median']):
    add(ax,b[k],se[k],c,'#165B83');ax.set_xlabel('Frequency (Hz)');ax.set_title(label,loc='left')
axs[0].set_ylabel('Difference in age coefficients')
fig.savefig(FIG/'age_quantile_contrasts.pdf');fig.savefig(ROOT/'qa'/'age_quantile_contrasts.png',dpi=180);plt.close(fig)
# All quantile levels are displayed in the supplementary surface figure.
fig,axs=plt.subplots(1,3,figsize=(7.1,2.9),layout='constrained')
for ax,name,title in zip(axs,['SC_2_0_main','SC_2_1_main','ST_2_0_main'],['(a) Age: Fpz-Cz','(b) Age: Pz-Oz','(c) Paired: Fpz-Cz']):
    b,se,c=read(name);b=b.reshape(9,60);se=se.reshape(9,60);lim=np.max(abs(b))
    mesh=ax.pcolormesh(f,t,b,shading='nearest',cmap='RdBu_r',vmin=-lim,vmax=lim,rasterized=True,antialiased=False,edgecolors='none')
    significant=(abs(b)>c*se)
    ax.contour(f,t,significant.astype(float),levels=[.5],colors='black',linewidths=.6)
    ax.set_title(title,loc='left',fontsize=9);ax.set_xlabel('Frequency (Hz)');ax.set_yticks([.1,.3,.5,.7,.9]);ax.set_xticks([5,15,25]);fig.colorbar(mesh,ax=ax,orientation='horizontal',pad=.15,shrink=.9)
axs[0].set_ylabel('Quantile level')
fig.savefig(FIG/'full_surfaces.pdf',dpi=400);fig.savefig(ROOT/'qa'/'full_surfaces.png',dpi=180);plt.close(fig)
fig,axs=plt.subplots(2,3,figsize=(7.1,4.6),sharex=True,layout='constrained')
variants=[('main','60 s, 1.5 Hz','#165B83','-'),('narrow','60 s, 1 Hz','#A45B20','--'),('wide','60 s, 2 Hz','#42806E','-.'),('short','30 s, 1.5 Hz','#777777',':')]
for row,study in enumerate(['SC','ST']):
    for var,label,col,ls in variants:
        b,se,c=read(f'{study}_2_0_{var}');b=b.reshape(9,60)
        for colidx,qi in enumerate([0,4,8]):axs[row,colidx].plot(f,b[qi],color=col,ls=ls,lw=1.2,label=label)
    for colidx,qi in enumerate([0,4,8]):
        ax=axs[row,colidx];ax.axhline(0,color='#444444',lw=.6);ax.set_title(f'({"abcdef"[row*3+colidx]}) '+('Age' if row==0 else 'Paired')+fr', $\tau = {t[qi]:.1f}$',loc='left',fontsize=9);ax.grid(axis='y',alpha=.18)
        ax.set_xlim(.5,30);ax.set_xticks([.5,10,20,30]);ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
        if row==1:ax.set_xlabel('Frequency (Hz)')
    axs[row,0].set_ylabel('Coefficient')
handles,labels=axs[0,1].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncol=4,frameon=False,fontsize=8)
fig.savefig(FIG/'sensitivity_eeg.pdf');fig.savefig(ROOT/'qa'/'sensitivity_eeg.png',dpi=180);plt.close(fig)

def line(v):return ' & '.join(v)+r' \\'+'\n'
def pformat(p):return '$<0.001$' if p<.001 else f'{p:.3f}'
def holm(vals):
    vals=np.asarray(vals);order=np.argsort(vals);r=np.empty(len(vals));r[order]=np.minimum(1,np.maximum.accumulate((len(vals)-np.arange(len(vals)))*vals[order]));return r
j={p.stem:json.loads(p.read_text()) for p in D.glob('*.json')}
text=r'''\begin{table}[htbp]
\centering
\caption{Numbers of included subjects and segments and global $p$-values. Quantile tests use all frequency--quantile grid points. The heterogeneity test compares every nonmedian quantile curve with the median curve. Adjusted heterogeneity $p$-values use Holm's procedure across the four stage--channel analyses within each study.}\label{tab:eeg}
\small
\setlength{\tabcolsep}{4pt}
\begin{tabular}{lllrrrr}
\toprule
Study / stage & Channel & $N$ & Segments & Quantile & Ordinary & Heterogeneity\\
 & & & & $p$ & $p$ & adjusted $p$\\
\midrule
'''
summary=[]
for study in ['SC','ST']:
    keys=[f'{study}_{stage}_{ch}_main' for stage in [2,5] for ch in [0,1]]
    adj=holm([j[k]['heterogeneity']['p'] for k in keys]);gadj=holm([j[k]['quantile']['p'] for k in keys])
    for k,hp,gp in zip(keys,adj,gadj):
        r=j[k];text+=line([('Age' if study=='SC' else 'Paired')+' / '+('stage 2' if r['stage']==2 else 'REM'),'Fpz--Cz' if r['channel']==0 else 'Pz--Oz',str(r['N']),str(r['segments']),pformat(r['quantile']['p']),pformat(r['ordinary']['p']),pformat(hp)])
        summary.append({'name':k,'heterogeneity_holm':float(hp),'quantile_holm':float(gp),**r})
    if study=='SC':text+='\\addlinespace\n'
text+='\\bottomrule\n\\end{tabular}\n\\end{table}\n';(MAN/'table_eeg.tex').write_text(text)
(D/'adjusted_summary.json').write_text(json.dumps(summary,indent=2))
text=r'''\begin{table}[htbp]
\centering
\caption{Sensitivity of the primary Fpz--Cz stage-2 analyses to window length and smoothing bandwidth. The same candidate intervals, quality exclusions, covariates, and subject weighting are used.}\label{tab:eeg-sensitivity}
\small
\begin{tabular}{lrrr r}
\toprule
Analysis & Window (s) & Bandwidth (Hz) & Quantile $p$ & Ordinary $p$\\
\midrule
'''
for study in ['SC','ST']:
    for variant,length,h in [('main',60,1.5),('narrow',60,1),('wide',60,2),('short',30,1.5)]:
        r=j[f'{study}_2_0_{variant}'];text+=line(['Age' if study=='SC' else 'Paired',str(length),str(h),pformat(r['quantile']['p']),pformat(r['ordinary']['p'])])
text+='\\bottomrule\n\\end{tabular}\n\\end{table}\n';(MAN/'table_eeg_sensitivity.tex').write_text(text)
print('Verified EEG figures and tables generated')
