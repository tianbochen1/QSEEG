import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
from hashlib import sha256
from contextlib import contextmanager
from tempfile import NamedTemporaryFile
import numpy as np,pandas as pd,json
from core import fit,SPECTRAL_FLOOR
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'/'eeg';OUT.mkdir(exist_ok=True)
FEATURES=ROOT/'data'/'features'

@contextmanager
def result_file(path, mode):
    """Replace each result only after it has been written in full."""
    temporary=NamedTemporaryFile(dir=path.parent,mode=mode,delete=False)
    pending=Path(temporary.name)
    try:
        with temporary:
            yield temporary
        pending.replace(path)
    finally:
        pending.unlink(missing_ok=True)

def load(study,stage,channel,variant):
    qs=[];ops=[];shapes=[];frames=[];ncandidate=0
    for path in sorted(FEATURES.glob(study+'*.npz')):
        with np.load(path) as a:
            h=json.loads(str(a['metadata']));rows=a['rows'];qual=a['quality']
            select=(rows[:,0]==stage);ncandidate+=int(select.sum());select &= qual[:,channel]
            if not select.any():continue
            q=a[variant+'_q'][:,channel]
            finite=np.isfinite(q).all(axis=(-2,-1))
            if variant+'_valid' in a:
                select &= a[variant+'_valid'][:,channel]
            else:
                # A floor-only curve in an older archive can also be a valid spectral zero.
                # The segment observations are needed to distinguish these two cases.
                floor_only=(q==q.dtype.type(np.log(SPECTRAL_FLOOR))).all(axis=-1).any(axis=-1)
                if np.any(select & finite & floor_only):
                    raise ValueError(f'{path.name}: rebuild features from its segments to verify the rank grid.')
            select &= finite
            if not select.any():continue
            qs.append(q[select]);ops.append(a[variant+'_o'][select,channel]);shapes.append(a[variant+'_s'][select,channel])
            frame=pd.DataFrame({k:[h[k]]*int(select.sum()) for k in ['subject','night','age','female','treatment','study','record']})
            frame['elapsed']=rows[select,2];frames.append(frame)
    if not frames:
        raise ValueError(f'No eligible segments for {study}, stage {stage}, channel {channel}, {variant}.')
    return pd.concat(frames,ignore_index=True),np.concatenate(qs),np.concatenate(ops),np.concatenate(shapes),ncandidate

def run(study,stage,channel,variant='main'):
    name=f'{study}_{stage}_{channel}_{variant}'
    df,q,o,sh,ncand=load(study,stage,channel,variant)
    digest=sha256(str(ncand).encode())
    digest.update(pd.util.hash_pandas_object(df,index=False).to_numpy().tobytes())
    for values in (q,o,sh):
        digest.update(np.ascontiguousarray(values).view(np.uint8))
    signature=digest.hexdigest()
    metadata=OUT/(name+'.json')
    outputs=['quantile','ordinary','shape','median','heterogeneity'] if variant=='main' else ['quantile','ordinary']
    if metadata.is_file() and (OUT/(name+'_included.csv')).is_file():
        try:
            saved_signature=json.loads(metadata.read_text()).get('input_sha256')
        except json.JSONDecodeError:
            saved_signature=None
        if saved_signature==signature and all(
                (OUT/(name+'_'+method+'.npz')).is_file() for method in outputs):
            return
    if study=='ST':
        paired=df.groupby('subject').treatment.transform('nunique').eq(2).to_numpy()
        if not paired.any():
            raise ValueError('The paired analysis requires retained segments in both conditions.')
        df=df.loc[paired].reset_index(drop=True);q,o,sh=q[paired],o[paired],sh[paired]
    originalM=len(df);originalrecords=df.record.nunique()
    if study=='SC':
        g=df.subject.to_numpy();elapsed=df.elapsed.to_numpy()-3
        X=np.c_[np.ones(len(df)),(df.age.to_numpy()-60)/10,df.female.to_numpy(),(df.night.to_numpy()==2).astype(int),elapsed,elapsed**2]
        cnt=df.groupby(['subject','night']).subject.transform('size').to_numpy()
        nights=df.groupby('subject').night.transform('nunique').to_numpy()
        w=1/(cnt*nights);coef=1
    else:
        ids=np.sort(df.subject.unique())
        qq=[];oo=[];ss=[];xx=[]
        for i in ids:
            a=(df.subject.to_numpy()==i)&(df.treatment.to_numpy()==1);b=(df.subject.to_numpy()==i)&(df.treatment.to_numpy()==0)
            qq.append(q[a].mean(0)-q[b].mean(0));oo.append(o[a].mean(0)-o[b].mean(0));ss.append(sh[a].mean(0)-sh[b].mean(0))
            xx.append([1,df.loc[a,'night'].iloc[0]-df.loc[b,'night'].iloc[0],df.loc[a,'elapsed'].mean()-df.loc[b,'elapsed'].mean()])
        q=np.asarray(qq);o=np.asarray(oo);sh=np.asarray(ss);X=np.asarray(xx);g=ids;coef=0;w=None
    rng=np.random.default_rng(284316+stage+channel*101+sum(map(ord,variant)))
    B=4999 if variant=='main' else 1999
    results={'input_sha256':signature}
    methods={'quantile':q.reshape(len(g),-1),'ordinary':o,'shape':sh,'median':q[:,4,:]}
    if variant=='main':
        # Direct contrasts against the median, excluding the identically zero row.
        methods['heterogeneity']=(q[:,np.arange(9)!=4,:]-q[:,4:5,:]).reshape(len(g),-1)
    metadata.unlink(missing_ok=True)
    for method,Y in methods.items():
        if variant!='main' and method not in ['quantile','ordinary']:continue
        v=fit(X,Y,g,coef=coef,bootstrap=B,rng=rng,weights=w)
        with result_file(OUT/(name+'_'+method+'.npz'),'wb') as stream:
            np.savez_compressed(stream,**v)
        lower=v['effect']-v['critical']*v['se'];upper=v['effect']+v['critical']*v['se']
        results[method]={'p':v['p'],'maxstat':v['maxstat'],'critical':v['critical'],'median_se':float(np.median(v['se'])),'min_effect':float(v['effect'].min()),'max_effect':float(v['effect'].max()),'positive_points':int((lower>0).sum()),'negative_points':int((upper<0).sum()),'gridpoints':len(Y[0])}
    if study=='SC' and variant=='main':
        naive=fit(X,q.reshape(len(g),-1),g,coef=coef,cluster=False,bootstrap=199,rng=rng,studentized=False,weights=w)
        with np.load(OUT/(name+'_quantile.npz')) as z:
            results['se_ratio_cluster_independent']=float(np.median(z['se']/naive['se']))
    subjects=df.drop_duplicates('subject')
    results.update(study=study,stage=int(stage),channel=int(channel),variant=variant,N=len(np.unique(g)),segments=originalM,recordings=int(originalrecords),candidate_segments=ncand,excluded_segments=ncand-originalM,age_min=float(subjects.age.min()),age_max=float(subjects.age.max()),age_median=float(subjects.age.median()),female=int(subjects.female.sum()),paired_subjects=g.tolist() if study=='ST' else None)
    with result_file(OUT/(name+'_included.csv'),'w') as stream:
        df.to_csv(stream,index=False)
    with result_file(metadata,'w') as stream:
        json.dump(results,stream,indent=2)
    print(name,results['N'],originalM,{k:round(v['p'],4) for k,v in results.items() if isinstance(v,dict)},flush=True)

if __name__=='__main__':
    assert len(list(FEATURES.glob('*.npz')))==197,'All recordings must be processed before analysis'
    for study in ['SC','ST']:
        for stage in [2,5]:
            for channel in [0,1]:run(study,stage,channel)
        for variant in ['narrow','wide','short']:run(study,2,0,variant)
