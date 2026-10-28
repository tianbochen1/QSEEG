import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import concurrent.futures as cf
import json,time
from simulate import one

configs=[]
for kind in ['ar','sv']:
    for effect in [.5,1.]:configs.append(dict(kind=kind,N=80,m=6,n=512,effect=effect))
    configs.append(dict(kind=kind,N=160,m=6,n=512,effect=0))
configs.append(dict(kind='surface',N=24,m=4,effect=0))
configs.append(dict(kind='sv',N=24,m=6,n=512,effect=0))
for k,c in enumerate(configs):c['id']=100+k

if __name__=='__main__':
    dest=Path(__file__).resolve().parents[1]/'results'/'extra_final.jsonl'
    done=set()
    if dest.exists():
        for l in dest.read_text().splitlines():
            r=json.loads(l);done.add((r['config']['id'],r['rep']))
    jobs=[(c,k) for c in configs for k in range(1000) if (c['id'],k) not in done]
    start=time.time()
    with cf.ProcessPoolExecutor(max_workers=3) as ex,open(dest,'a') as out:
        for k,r in enumerate(ex.map(one,jobs,chunksize=2),1):
            out.write(json.dumps(r)+'\n');out.flush()
            if k%100==0: print(k,'elapsed',round(time.time()-start),flush=True)
