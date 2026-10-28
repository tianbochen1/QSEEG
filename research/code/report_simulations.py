from pathlib import Path
import json,numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results';MAN=ROOT/'manuscript';FIG=ROOT.parent/'output'/'QSEEG'/'frame'
FIG.mkdir(parents=True,exist_ok=True)
(ROOT/'qa').mkdir(parents=True,exist_ok=True)
rows=[];unique={};counts={}
for filename in ['simulation_final.jsonl','extra_final.jsonl']:
    for l in (OUT/filename).read_text().splitlines():
        r=json.loads(l);key=(r['config']['id'],r['rep'])
        if key in unique:
            for m in r['metrics']:
                for k,v in r['metrics'][m].items():assert np.isclose(v,unique[key]['metrics'][m][k],atol=1e-10)
        else:unique[key]=r
for (cid,rep),r in unique.items():
    counts[cid]=counts.get(cid,0)+1
    for method,v in r['metrics'].items():rows.append({**r['config'],'method':method,**v})
assert len(counts)==25 and set(counts.values())=={1000},counts
(OUT/'simulation_verified.jsonl').write_text(''.join(json.dumps(unique[k])+'\n' for k in sorted(unique)))
df=pd.DataFrame(rows)
summ=df.groupby(['id','kind','N','m','effect','method'],dropna=False).agg(reject=('reject','mean'),cover=('cover','mean'),width=('width','mean'),rmse=('rmse',lambda x:np.sqrt(x.pow(2).mean()))).reset_index()
summ.to_csv(OUT/'simulation_summary_final.csv',index=False)
by={(int(r.id),r.method):r for _,r in summ.iterrows()}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42,'axes.titlesize':10})
colors={'Cluster':'#165B83','Average':'#A45B20','Independent':'#727272','Quantile':'#165B83','Ordinary':'#A45B20','Median':'#42806E'}
def mc_interval_errors(probability,n=1000):
    z=1.959963984540054
    p=np.asarray(probability)
    den=1+z*z/n
    center=(p+z*z/(2*n))/den
    half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return np.maximum(np.stack([p-np.maximum(0,center-half),np.minimum(1,center+half)-p]),0)
fig,axs=plt.subplots(1,2,figsize=(7,3.2),layout='constrained')
x=np.arange(4);labels=['80 / 4','80 / 12','160 / 4','160 / 12']
for method,marker in [('Cluster','o'),('Average','s'),('Independent','^')]:
    yy=np.array([by[(i,method)]['cover'] for i in range(4)])
    axs[0].errorbar(x,yy,yerr=mc_interval_errors(yy),marker=marker,color=colors[method],lw=1.2,capsize=2,label=method)
    axs[1].plot(x,[by[(i,method)]['width'] for i in range(4)],marker=marker,color=colors[method],lw=1.2,label=method)
axs[0].axhline(.95,color='#333333',ls='--',lw=.8);axs[0].set_ylim(0,1.03);axs[0].set_ylabel('Simultaneous coverage')
axs[1].set_ylabel('Mean band width')
for a,title in zip(axs,['(a) Coverage under the null','(b) Width under the null']):
    a.set_xticks(x,labels);a.set_xlabel('Subjects / segments per subject');a.set_title(title,loc='left');a.grid(axis='y',alpha=.18)
handles,labels=axs[1].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncol=3,frameon=False,fontsize=8)
fig.savefig(FIG/'simulation_calibration.pdf');fig.savefig(ROOT/'qa'/'simulation_calibration.png',dpi=180);plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(7,3.1),layout='constrained')
for ax,ids,title in [(axs[0],[9,100,101,10],'(a) Gaussian autoregression'),(axs[1],[12,103,104,13],'(b) Stochastic volatility')]:
    for method,marker in [('Quantile','o'),('Ordinary','s'),('Median','^')]:
        yy=np.array([by[(i,method)]['reject'] for i in ids]);ax.errorbar([0,.5,1,2],yy,yerr=mc_interval_errors(yy),marker=marker,lw=1.2,capsize=2,color=colors[method],label=method)
    ax.set_title(title,loc='left');ax.set_xlabel(r'Exposure effect in persistence ($\delta$)');ax.set_ylim(-.02,1.04);ax.set_xticks([0,.5,1,2]);ax.grid(axis='y',alpha=.18)
axs[0].set_ylabel('Global rejection probability');handles,labels=axs[1].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncol=3,frameon=False,fontsize=8)
fig.savefig(FIG/'simulation_sensitivity.pdf');fig.savefig(ROOT/'qa'/'simulation_sensitivity.png',dpi=180);plt.close(fig)

def line(vals):return ' & '.join(vals)+r' \\'+'\n'
def table_start(caption,label,cols,head):return '\\begin{table}[htbp]\n\\centering\n\\caption{'+caption+'}\\label{'+label+'}\n\\small\n\\begin{tabular}{'+cols+'}\n\\toprule\n'+line(head)+'\\midrule\n'
def table_end():return '\\bottomrule\n\\end{tabular}\n\\end{table}\n'
text=table_start('Global rejection probability, simultaneous coverage, and estimation error in the regression-surface experiment. Width is averaged over the grid. The nominal rejection probability under $\\delta=0$ is 0.05, and nominal simultaneous coverage is 0.95.','tab:surface','rrrlrrrr',['$N$','$m$','$\\delta$','Method','Reject','Cover','Width','RMSE'])
for cid in [0,1,2,3,4,5,6,7]:
    for method in ['Cluster','Average','Independent']:
        r=by[cid,method];text+=line([str(int(r.N)),str(int(r.m)),f'{r.effect:.2f}',method]+[f'{r[k]:.3f}' for k in ['reject','cover','width','rmse']])
    if cid!=7:text+='\\addlinespace[2pt]\n'
text+=table_end();(MAN/'table_surface.tex').write_text(text)
text=table_start('Informative segment counts with $N=160$ and $\\delta=0$. Coverage and error are evaluated for the subject-weighted exposure coefficient. Segment-weighted regression retains a subject-level sandwich covariance but targets a different population projection.','tab:counts','lrrrr',['Method','Reject','Cover','Width','RMSE'])
for method in ['Cluster','Average','Independent','Segment-weighted']:
    r=by[8,method];text+=line([method]+[f'{r[k]:.3f}' for k in ['reject','cover','width','rmse']])
text+=table_end();(MAN/'table_counts.tex').write_text(text)
text=table_start('Global rejection probabilities for raw time-series experiments. All designs have six segments per subject. The full quantile method tests all five quantiles jointly; median and ordinary methods test frequency curves. The independent analysis uses the quantile surface.','tab:raw','lrrrrrrr',['Process','$N$','$n$','$\\delta$','Quantile','Median','Ordinary','Independent'])
configs={r['config']['id']:r['config'] for r in unique.values()}
for cid in [9,102,100,101,10,11,12,105,103,104,13,14,15]:
    c=configs[cid];lab='AR' if c['kind']=='ar' else 'SV'
    if cid==16:lab='SV + contamination'
    text+=line([lab,str(c['N']),str(c['n']),f"{c['effect']:.1f}"]+[f'{by[cid,m]["reject"]:.3f}' for m in ['Quantile','Median','Ordinary','Independent']])
text+=table_end();(MAN/'table_raw.tex').write_text(text)
text=table_start('Sensitivity analysis with $N=24$ subjects and $\\delta=0$. The surface model has four segments per subject, and its coverage refers to simultaneous bands for the known exposure coefficient. The stochastic-volatility model has six segments and 512 observations per segment.','tab:small','llrr',['Design','Method','Reject','Cover'])
for cid,methods in [(106,['Cluster','Average','Independent']),(107,['Quantile','Median','Ordinary','Independent'])]:
    for m in methods:
        r=by[cid,m];text+=line(['Surface' if cid==106 else 'Volatility',m,f'{r["reject"]:.3f}',f'{r["cover"]:.3f}' if cid==106 else '--'])
text+=table_end();(MAN/'table_small.tex').write_text(text)
# Gaussian multipliers with fixed standard errors.
old=[]
for l in (OUT/'gaussian_comparison_final.jsonl').read_text().splitlines():
    r=json.loads(l)
    if r['config']['id'] in [0,1,2,3]:old.append({'id':r['config']['id'],**r['metrics']['Cluster']})
odf=pd.DataFrame(old).groupby('id')[['reject','cover','width']].mean()
text=table_start('Comparison of bootstrap constructions in the null surface experiment. Both use the same subject-weighted estimator and observed cluster sandwich standard errors. They differ in the multiplier distribution and whether standard errors are recomputed after refitting.','tab:studentization','rrlrrr',['$N$','$m$','Bootstrap','Reject','Cover','Width'])
for cid in range(4):
    c=configs[cid]
    for label,r in [('Gaussian, fixed standard error',odf.loc[cid]),('Wild, recomputed standard error',by[cid,'Cluster'])]:text+=line([str(c['N']),str(c['m']),label]+[f'{r[k]:.3f}' for k in ['reject','cover','width']])
text+=table_end();(MAN/'table_studentization.tex').write_text(text)
print(summ.round(3).to_string(index=False))
