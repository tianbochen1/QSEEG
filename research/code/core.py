import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.stats import rankdata

TAUS=np.arange(1,10)/10
FREQS=np.arange(1,61)*.005

def spectral_features(x, taus=TAUS, freqs=FREQS, bandwidth=.015, floor=1e-6):
    x=np.asarray(x,dtype=float)
    n=x.shape[-1]
    r=rankdata(x,axis=-1,method='average')/(n+1)
    score=(r[...,None,:]<=np.asarray(taus)[:,None]).astype(float)
    score-=score.mean(axis=-1,keepdims=True)
    den=score.var(axis=-1,keepdims=True)
    ps=np.abs(np.fft.fft(score,axis=-1))**2/n/np.maximum(den,1e-12)
    ps=gaussian_filter1d(ps,bandwidth*n,axis=-1,mode='wrap',truncate=4)
    idx=np.asarray(freqs)*n;lo=np.floor(idx).astype(int);w=idx-lo
    s=ps[...,lo]*(1-w)+ps[...,lo+1]*w
    q=np.log(s+floor)
    xx=x-x.mean(axis=-1,keepdims=True)
    pg=np.abs(np.fft.fft(xx,axis=-1))**2/n
    pg=gaussian_filter1d(pg,bandwidth*n,axis=-1,mode='wrap',truncate=4)
    # Ordinary log power retains physical amplitude; standardized shape is also available.
    ss=pg[...,lo]*(1-w)+pg[...,lo+1]*w
    ordinary=np.log(ss+floor)
    shape=np.log(ss/np.maximum(xx.var(axis=-1,keepdims=True),1e-12)+floor)
    return q,ordinary,shape

def fit(X,Y,groups,coef=1,cluster=True,equal_subject=True,bootstrap=499,rng=None,studentized=True,weights=None):
    X=np.asarray(X,float);Y=np.asarray(Y,float)
    _,g=np.unique(groups,return_inverse=True)
    N=g.max()+1;M=len(g);p=X.shape[1]
    count=np.bincount(g)
    w=(1/count[g] if equal_subject else np.ones(M)) if weights is None else np.asarray(weights,float)
    inv=np.linalg.inv(X.T@(w[:,None]*X))
    bhat=inv@(X.T@(w[:,None]*Y))
    resid=Y-X@bhat
    row=(X@inv[coef])*w
    ri=row[:,None]*resid
    if cluster:
        inf=np.zeros((N,Y.shape[1]));np.add.at(inf,g,ri)
        # Cluster HC1, tending to one as the number of subjects increases.
        factor=N/(N-1)*(M-1)/(M-p)
    else:
        inf=ri;factor=M/(M-p)
    inf*=np.sqrt(factor)
    se=np.sqrt(np.sum(inf**2,axis=0))
    se=np.maximum(se,1e-12)
    rng=np.random.default_rng() if rng is None else rng
    multipliers=rng.choice([-1.,1.],size=(bootstrap,len(inf))) if studentized else rng.standard_normal((bootstrap,len(inf)))
    if studentized:
        # Recompute the sandwich standard error under the subject-level wild bootstrap.
        # This is algebraically identical to refitting the linear model in every draw.
        gg=g if cluster else np.arange(M)
        nn=len(inf)
        score=np.zeros((nn,p,Y.shape[1])); gram=np.zeros((nn,p,p))
        np.add.at(score,gg,w[:,None,None]*X[:,:,None]*resid[:,None,:])
        np.add.at(gram,gg,w[:,None,None]*X[:,:,None]*X[:,None,:])
        direction=np.einsum('p,npq->nq',inv[coef],gram)
        rawinf=inf/np.sqrt(factor)
        total_score=(multipliers@score.reshape(nn,-1)).reshape(bootstrap,p,-1)
        delta=np.einsum('pq,bqd->bpd',inv,total_score,optimize=True)
        cross=(multipliers@(direction[:,:,None]*rawinf[:,None,:]).reshape(nn,-1)).reshape(bootstrap,p,-1)
        base=np.sum(rawinf**2,axis=0)
        gram_direction=direction.T@direction
        varstar=base[None,:]-2*np.sum(delta*cross,axis=1)+np.einsum('bpd,pq,bqd->bd',delta,gram_direction,delta,optimize=True)
        sestar=np.sqrt(factor*np.maximum(varstar,1e-24))
        maxdist=np.max(np.abs(delta[:,coef,:]/sestar),axis=1)
    else:
        maxdist=np.max(np.abs(multipliers@(inf/se)),axis=1)
    # Match the strict Monte Carlo tail-probability test to its inversion.
    critical_index=min(bootstrap-1,int(np.floor(.95*(bootstrap+1))))
    crit=float(np.sort(maxdist)[critical_index])
    observed=float(np.max(np.abs(bhat[coef]/se)))
    pvalue=float((1+np.sum(maxdist>=observed))/(bootstrap+1))
    return {'beta':bhat,'effect':bhat[coef],'se':se,'influence':inf,'critical':crit,'p':pvalue,'maxstat':observed,'N':N,'M':M,'maxdist':maxdist}

def aggregate_fit(X,Y,g,coef=1,bootstrap=499,rng=None):
    u,gi=np.unique(g,return_inverse=True);n=len(u)
    xx=np.zeros((n,X.shape[1]));yy=np.zeros((n,Y.shape[1]));counts=np.bincount(gi)
    np.add.at(xx,gi,X);np.add.at(yy,gi,Y)
    return fit(xx/counts[:,None],yy/counts[:,None],np.arange(n),coef,True,True,bootstrap,rng)

def ar_samples(phi,m,n,rng):
    phi=np.repeat(np.asarray(phi),m)
    out=np.empty((len(phi),n));out[:,0]=rng.standard_normal(len(phi))
    sd=np.sqrt(1-phi**2)
    e=rng.standard_normal((len(phi),n))
    for t in range(1,n):out[:,t]=phi*out[:,t-1]+sd*e[:,t]
    return out
