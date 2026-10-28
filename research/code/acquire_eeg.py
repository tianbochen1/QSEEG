from pathlib import Path
import concurrent.futures as cf
import hashlib, json, time, urllib.request, traceback
import numpy as np
import pandas as pd
import pyedflib
from read_annotations import read_hypnogram

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'data'
BASE = 'https://physionet-open.s3.amazonaws.com/sleep-edfx/1.0.0/'
META = pd.read_excel(DATA/'SC-subjects.xls')
ST = pd.read_excel(DATA/'ST-subjects.xls', header=1)
HASHES = dict((s.split()[1], s.split()[0]) for s in (DATA/'SHA256SUMS.txt').read_text().splitlines() if len(s.split()) == 2)
RECORDS = (DATA/'RECORDS').read_text().splitlines()
OUT = DATA/'segments'
OUT.mkdir(exist_ok=True)

def fetch(rel, path):
    if path.exists() and hashlib.sha256(path.read_bytes()).hexdigest() == HASHES.get(rel):
        return
    last = None
    for attempt in range(3):
        try:
            req = urllib.request.Request(BASE+rel, headers={'User-Agent':'Research-data-analysis/1.0'})
            with urllib.request.urlopen(req, timeout=90) as r, open(str(path)+'.part','wb') as f:
                if r.status != 200:
                    raise RuntimeError(f'Unexpected HTTP status {r.status}')
                while True:
                    b=r.read(1024*1024)
                    if not b: break
                    f.write(b)
            pp=Path(str(path)+'.part')
            digest=hashlib.sha256(pp.read_bytes()).hexdigest()
            if rel in HASHES and digest != HASHES[rel]:
                raise RuntimeError('SHA256 mismatch')
            pp.replace(path)
            return
        except Exception as e:
            last=e
    raise last

def process(rel):
    name=Path(rel).stem.replace('-PSG','')
    target=OUT/(name+'.npz')
    if target.exists(): return name,'cached'
    subject=int(name[3:5]); night=int(name[5]); study=name[:2]
    hyprel=[x for x in HASHES if Path(x).name.startswith(name[:6]) and x.endswith('Hypnogram.edf')]
    if len(hyprel)!=1: raise RuntimeError((name,hyprel))
    hyp=DATA/Path(hyprel[0]).name
    raw=DATA/(name+'-PSG.edf')
    fetch(hyprel[0],hyp); fetch(rel,raw)
    onset,duration,label,hypstart=read_hypnogram(hyp)
    intervals=[(float(o),float(d),str(l)) for o,d,l in zip(onset,duration,label)]
    with pyedflib.EdfReader(str(raw)) as r:
        start=r.getStartdatetime(); offset=(hypstart-start).total_seconds()
        labels=r.getSignalLabels(); channels=[labels.index(x) for x in ['EEG Fpz-Cz','EEG Pz-Oz']]
        fs=np.array([r.getSampleFrequency(c) for c in channels])
        if not np.all(fs==100): raise RuntimeError(('sampling',name,fs))
        total=r.getFileDuration()
        sleep_onsets=[o+offset for o,d,l in intervals if l in ['Sleep stage 1','Sleep stage 2','Sleep stage 3','Sleep stage 4','Sleep stage R']]
        sleep_start=min(sleep_onsets)
        rows=[]; signals=[]; quality=[]
        for stage in ['Sleep stage 2','Sleep stage R']:
            candidates=[]
            for o,d,l in intervals:
                if l!=stage:continue
                # Use 60-second windows wholly contained in one scored stage.
                for a in np.arange(o+offset, o+offset+d-60+0.1, 60):
                    if a>=0 and a+60<=total:candidates.append(float(a))
            if not candidates:continue
            # Fixed uniform coverage of the eligible stage, independent of EEG values.
            inds=np.unique(np.round(np.linspace(0,len(candidates)-1,min(20,len(candidates)))).astype(int))
            for ix in inds:
                a=candidates[ix]
                arr=np.stack([r.readSignal(c,start=int(round(a*100)),n=6000) for c in channels])
                dig=np.stack([r.readSignal(c,start=int(round(a*100)),n=6000,digital=True) for c in channels])
                amp=np.max(np.abs(arr-np.median(arr,axis=1,keepdims=True)),axis=1)
                sd=np.std(arr,axis=1)
                clipped=np.array([np.mean((dig[k]<=r.getDigitalMinimum(c))|(dig[k]>=r.getDigitalMaximum(c))) for k,c in enumerate(channels)])
                good=(amp<=500)&(sd>=0.5)&(clipped<0.01)&np.all(np.isfinite(arr),axis=1)
                signals.append(arr.astype('float32')); quality.append(good)
                rows.append([2 if stage=='Sleep stage 2' else 5,a,(a-sleep_start)/3600])
        header={'record':name,'subject':subject,'night':night,'study':study,'fs':100,'channels':['Fpz-Cz','Pz-Oz'],'psg_sha256':HASHES.get(rel),'hypnogram':hyprel[0],'hypnogram_sha256':HASHES.get(hyprel[0]),'header_patient':r.getPatientCode(),'psg_start':str(start),'hypnogram_start':str(hypstart),'sampling_rule':'20 uniformly spaced eligible nonoverlapping 60-second windows per stage per night'}
    if study=='SC':
        m=META[(META.subject==subject)&(META.night==night)].iloc[0]
        header.update(age=float(m['age']),female=int(m['sex (F=1)']==1),treatment=0)
    else:
        m=ST[ST['Nr']==subject].iloc[0]
        header.update(age=float(m['Age']),female=int(m['M1/F2']==2),treatment=int(night==m['night nr.1']))
    np.savez_compressed(target,signals=np.asarray(signals),rows=np.asarray(rows),quality=np.asarray(quality),metadata=json.dumps(header))
    raw.unlink()
    return name,f'{len(signals)} windows'

if __name__=='__main__':
    errors=[]
    t=time.time()
    with cf.ThreadPoolExecutor(max_workers=10) as ex:
        futures={ex.submit(process,r):r for r in RECORDS}
        for k,f in enumerate(cf.as_completed(futures),1):
            try:
                name,status=f.result();print(k,len(RECORDS),name,status,'elapsed',round(time.time()-t),flush=True)
            except Exception as e:
                errors.append({'record':futures[f],'error':repr(e),'traceback':traceback.format_exc()})
                print('ERROR',futures[f],repr(e),flush=True)
    (DATA/'acquisition_errors.json').write_text(json.dumps(errors,indent=2))
    print('DONE',len(RECORDS)-len(errors),'errors',len(errors),flush=True)
