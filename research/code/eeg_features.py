import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import numpy as np
import json,time,concurrent.futures as cf
from core import spectral_features,TAUS,FREQS
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data'/'features';OUT.mkdir(exist_ok=True)

def one(p):
    target=OUT/p.name
    if target.exists():
        try:
            with np.load(target) as cached:
                cached['main_q']
            return p.stem,'cached'
        except Exception:pass
    with np.load(p) as a:
        sig=a['signals'];rows=a['rows'];qual=a['quality'];meta=a['metadata']
    n,ch,t=sig.shape
    vals={}
    for label,h,l in [('main',.015,6000),('narrow',.01,6000),('wide',.02,6000),('short',.015,3000)]:
        q,o,s=spectral_features(sig[:,:,:l].reshape(n*ch,l),TAUS,FREQS,bandwidth=h)
        vals[label+'_q']=q.reshape(n,ch,len(TAUS),len(FREQS)).astype('float32')
        vals[label+'_o']=o.reshape(n,ch,len(FREQS)).astype('float32')
        vals[label+'_s']=s.reshape(n,ch,len(FREQS)).astype('float32')
    temporary=target.with_suffix('.tmp')
    with open(temporary,'wb') as dest:
        np.savez_compressed(dest,rows=rows,quality=qual,metadata=meta,**vals)
    os.replace(temporary,target)
    return p.stem,len(rows)

if __name__=='__main__':
    files=sorted((ROOT/'data'/'segments').glob('*.npz'));t=time.time()
    with cf.ProcessPoolExecutor(max_workers=3) as ex:
        for k,res in enumerate(ex.map(one,files),1):
            if k%10==0:print(k,len(files),res,round(time.time()-t),flush=True)
