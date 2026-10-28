import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import numpy as np,pandas as pd,json
from core import fit,TAUS,FREQS
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'/'eeg';OUT.mkdir(exist_ok=True)
FEATURES=ROOT/'data'/'features'

def load(study,stage,channel,variant):
    qs=[];ops=[];shapes=[];frames=[];ncandidate=0
    for path in sorted(FEATURES.glob(study+'*.npz')):
        with np.load(path) as a:
            h=json.loads(str(a['metadata']));rows=a['rows'];qual=a['quality']
            select=(rows[:,0]==stage);ncandidate+=int(select.sum());select &= qual[:,channel]
            if not select.any():continue
            qs.append(a[variant+'_q'][select,channel]);ops.append(a[variant+'_o'][select,channel]);shapes.append(a[variant+'_s'][select,channel])
            frame=pd.DataFrame({k:[h[k]]*int(select.sum()) for k in ['subject','night','age','female','treatment','study','record']})
            frame['elapsed']=rows[select,2];frames.append(frame)
    return pd.concat(frames,ignore_index=True),np.concatenate(qs),np.concatenate(ops),np.concatenate(shapes),ncandidate

def run(study,stage,channel,variant='main'):
    name=f'{study}_{stage}_{channel}_{variant}'
    if (OUT/(name+'.json')).exists():return
    df,q,o,sh,ncand=load(study,stage,channel,variant)
    originalM=len(df);originalrecords=df.record.nunique()
    if study=='SC':
        g=df.subject.to_numpy();elapsed=df.elapsed.to_numpy()-3
        X=np.c_[np.ones(len(df)),(df.age.to_numpy()-60)/10,df.female.to_numpy(),(df.night.to_numpy()==2).astype(int),elapsed,elapsed**2]
        cnt=df.groupby(['subject','night']).subject.transform('size').to_numpy()
        nights=df.groupby('subject').night.transform('nunique').to_numpy()
        w=1/(cnt*nights);coef=1
    else:
        complete=df.groupby('subject').treatment.nunique();ids=complete[complete==2].index.to_numpy()
        qq=[];oo=[];ss=[];xx=[]
        for i in ids:
            a=(df.subject.to_numpy()==i)&(df.treatment.to_numpy()==1);b=(df.subject.to_numpy()==i)&(df.treatment.to_numpy()==0)
            qq.append(q[a].mean(0)-q[b].mean(0));oo.append(o[a].mean(0)-o[b].mean(0));ss.append(sh[a].mean(0)-sh[b].mean(0))
            xx.append([1,df.loc[a,'night'].iloc[0]-df.loc[b,'night'].iloc[0],df.loc[a,'elapsed'].mean()-df.loc[b,'elapsed'].mean()])
        q=np.asarray(qq);o=np.asarray(oo);sh=np.asarray(ss);X=np.asarray(xx);g=ids;coef=0;w=None
    rng=np.random.default_rng(284316+stage+channel*101+sum(map(ord,variant)))
    B=4999 if variant=='main' else 1999
    results={}
    methods={'quantile':q.reshape(len(g),-1),'ordinary':o,'shape':sh,'median':q[:,4,:]}
    if variant=='main':
        # Direct contrasts against the median, excluding the identically zero row.
        methods['heterogeneity']=(q[:,np.arange(9)!=4,:]-q[:,4:5,:]).reshape(len(g),-1)
    for method,Y in methods.items():
        if variant!='main' and method not in ['quantile','ordinary']:continue
        v=fit(X,Y,g,coef=coef,bootstrap=B,rng=rng,weights=w)
        np.savez_compressed(OUT/(name+'_'+method+'.npz'),**v)
        lower=v['effect']-v['critical']*v['se'];upper=v['effect']+v['critical']*v['se']
        results[method]={'p':v['p'],'maxstat':v['maxstat'],'critical':v['critical'],'median_se':float(np.median(v['se'])),'min_effect':float(v['effect'].min()),'max_effect':float(v['effect'].max()),'positive_points':int((lower>0).sum()),'negative_points':int((upper<0).sum()),'gridpoints':len(Y[0])}
    if study=='SC' and variant=='main':
        naive=fit(X,q.reshape(len(g),-1),g,coef=coef,cluster=False,bootstrap=199,rng=rng,studentized=False,weights=w)
        with np.load(OUT/(name+'_quantile.npz')) as z:
            results['se_ratio_cluster_independent']=float(np.median(z['se']/naive['se']))
    subjects=df.drop_duplicates('subject')
    results.update(study=study,stage=int(stage),channel=int(channel),variant=variant,N=len(np.unique(g)),segments=originalM,recordings=int(originalrecords),candidate_segments=ncand,excluded_segments=ncand-originalM,age_min=float(subjects.age.min()),age_max=float(subjects.age.max()),age_median=float(subjects.age.median()),female=int(subjects.female.sum()),paired_subjects=g.tolist() if study=='ST' else None)
    (OUT/(name+'.json')).write_text(json.dumps(results,indent=2))
    df.to_csv(OUT/(name+'_included.csv'),index=False)
    print(name,results['N'],originalM,{k:round(v['p'],4) for k,v in results.items() if isinstance(v,dict)},flush=True)

if __name__=='__main__':
    assert len(list(FEATURES.glob('*.npz')))==197,'All recordings must be processed before analysis'
    for study in ['SC','ST']:
        for stage in [2,5]:
            for channel in [0,1]:run(study,stage,channel)
        for variant in ['narrow','wide','short']:run(study,2,0,variant)
