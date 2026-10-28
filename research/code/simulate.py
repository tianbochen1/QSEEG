import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
from pathlib import Path
import concurrent.futures as cf
import argparse,json,time
import numpy as np
from scipy.special import expit
from core import fit,aggregate_fit,spectral_features,ar_samples

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
F=np.linspace(.015,.30,30);T=np.array([.1,.3,.5,.7,.9])

def one(job):
    cfg,rep=job;rng=np.random.default_rng(736281+100000*cfg['id']+rep)
    N=cfg['N'];m=cfg['m'];D=len(F)*len(T)
    z=rng.standard_normal(N);a=.5*z+rng.standard_normal(N)
    u=rng.standard_normal(N)
    if cfg['kind']=='surface':
        counts=np.full(N,m)
        if cfg.get('unbalanced'):counts=np.where(a*u>0,12,4)
        g=np.repeat(np.arange(N),counts);M=len(g)
        v=rng.standard_normal(N)
        w=.6*a[g]+.2*z[g]+.5*v[g]+rng.standard_normal(M)
        X=np.c_[np.ones(M),a[g],z[g],w]
        s1=np.exp(-.5*((F-.08)/.035)**2)
        s2=np.exp(-.5*((F-.21)/.045)**2)
        bA=cfg['effect']*((2*T[:,None]-1)**2)*(s1[None,:]-.6*s2[None,:])
        bZ=.2*np.ones((len(T),1))*np.cos(2*np.pi*F)[None,:]
        bW=.4*(1+.5*(T[:,None]-.5))*s2[None,:]
        basis=np.array([np.outer(np.cos(np.pi*k*T),np.cos(2*np.pi*l*F/.6)).ravel() for k,l in [(0,0),(0,1),(1,1),(2,1),(0,2),(1,2)]])
        basis/=np.sqrt(6)
        bc=rng.standard_normal((N,6));bc[:,0]=u
        b=bc@basis*.55
        eps=rng.standard_normal((M,6))@basis*.35
        # Additional within-subject dependence beyond the random spectral curve.
        for i in range(N):
            ix=np.flatnonzero(g==i)
            for j in range(1,len(ix)):eps[ix[j]]=.4*eps[ix[j-1]]+np.sqrt(.84)*eps[ix[j]]
        Y=a[g,None]*bA.ravel()+z[g,None]*bZ.ravel()+w[:,None]*bW.ravel()+b[g]+eps
        truth=bA.ravel()
        methods={
          'Cluster':fit(X,Y,g,bootstrap=399,rng=rng),
          'Independent':fit(X,Y,g,cluster=False,bootstrap=399,rng=rng),
          'Average':aggregate_fit(X,Y,g,bootstrap=399,rng=rng)}
        if cfg.get('unbalanced'):
            methods['Segment-weighted']=fit(X,Y,g,equal_subject=False,bootstrap=399,rng=rng)
        result={}
        for name,v in methods.items():
            result[name]={'reject':int(v['p']<.05),'cover':int(np.all(np.abs(v['effect']-truth)<=v['critical']*v['se'])),'width':float(np.mean(2*v['critical']*v['se'])),'rmse':float(np.sqrt(np.mean((v['effect']-truth)**2)))}
        return {'config':cfg,'rep':rep,'metrics':result}
    n=cfg['n'];a=np.r_[np.full(N//2,-.5),np.full(N-N//2,.5)];rng.shuffle(a)
    phi=.08+.84*expit(cfg['effect']*a+.4*z+.9*u)
    if cfg['kind']=='ar':x=ar_samples(phi,m,n,rng)
    else:
        h=1.5*ar_samples(phi,m,n,rng)
        x=np.exp((h-1.5**2/2)/2)*rng.standard_normal(h.shape)
    if cfg.get('contamination',0):
        mask=rng.random(x.shape)<cfg['contamination'];x=x+mask*10*rng.standard_normal(x.shape)
    q,pg,sh=spectral_features(x,T,F,bandwidth=.02)
    g=np.repeat(np.arange(N),m);X=np.c_[np.ones(N*m),a[g],z[g]]
    methods={
      'Quantile':fit(X,q.reshape(N*m,-1),g,bootstrap=399,rng=rng),
      'Ordinary':fit(X,pg,g,bootstrap=399,rng=rng),
      'Median':fit(X,q[:,2,:],g,bootstrap=399,rng=rng),
      'Independent':fit(X,q.reshape(N*m,-1),g,cluster=False,bootstrap=399,rng=rng)}
    return {'config':cfg,'rep':rep,'metrics':{k:{'reject':int(v['p']<.05),'width':float(np.mean(2*v['critical']*v['se']))} for k,v in methods.items()}}

def configurations():
    c=[]
    for effect in [0,.35]:
        for N in [80,160]:
            for m in [4,12]:c.append(dict(kind='surface',N=N,m=m,effect=effect))
    c.append(dict(kind='surface',N=160,m=8,effect=0,unbalanced=True))
    for kind in ['ar','sv']:
        c.append(dict(kind=kind,N=80,m=6,n=512,effect=0))
        for N in [80,160]:
            c.append(dict(kind=kind,N=N,m=6,n=512,effect=2.0))
    c.append(dict(kind='sv',N=80,m=6,n=1024,effect=2.0))
    c.append(dict(kind='sv',N=80,m=6,n=512,effect=2.0,contamination=.01))
    for k,v in enumerate(c):v['id']=k
    return c

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--reps',type=int,default=500);p.add_argument('--workers',type=int,default=6);p.add_argument('--pilot',action='store_true');a=p.parse_args()
    configs=configurations();dest=OUT/('pilot.jsonl' if a.pilot else 'simulation_final.jsonl')
    done=set()
    if dest.exists():
        for line in dest.read_text().splitlines():
            r=json.loads(line);done.add((r['config']['id'],r['rep']))
    jobs=[(c,k) for c in configs for k in range(a.reps) if (c['id'],k) not in done]
    start=time.time();print('jobs',len(jobs),flush=True)
    with cf.ProcessPoolExecutor(max_workers=a.workers) as ex,open(dest,'a') as out:
        for k,r in enumerate(ex.map(one,jobs,chunksize=2),1):
            out.write(json.dumps(r)+'\n');out.flush()
            if k%25==0:print('done',k,'seconds',round(time.time()-start),flush=True)
    print('finished',round(time.time()-start),flush=True)
