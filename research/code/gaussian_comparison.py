import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import concurrent.futures as cf,json
import core,simulate
original=core.fit

def fixed_gaussian(*args,**kwargs):
    kwargs['studentized']=False
    return original(*args,**kwargs)
core.fit=fixed_gaussian
simulate.fit=fixed_gaussian

if __name__=='__main__':
    configs=simulate.configurations()[:4]
    dest=Path(__file__).resolve().parents[1]/'results'/'gaussian_comparison_final.jsonl'
    with cf.ProcessPoolExecutor(max_workers=2) as ex,open(dest,'w') as out:
        for r in ex.map(simulate.one,[(c,k) for c in configs for k in range(1000)],chunksize=4):out.write(json.dumps(r)+'\n')
    print('Complete Gaussian comparison')
