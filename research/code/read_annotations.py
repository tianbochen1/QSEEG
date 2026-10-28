from pathlib import Path
from datetime import datetime
import numpy as np

def read_hypnogram(path):
    raw=Path(path).read_bytes()
    header=int(raw[184:192]);ns=int(raw[252:256])
    if ns!=1 or raw[256:272].decode().strip()!='EDF Annotations':raise ValueError('Unexpected annotation layout')
    date=raw[168:176].decode();time=raw[176:184].decode()
    day,month,year=map(int,date.split('.'));year+=1900 if year>=85 else 2000
    hour,minute,second=map(int,time.split('.'))
    start=datetime(year,month,day,hour,minute,second)
    onset=[];duration=[];labels=[]
    for chunk in raw[header:].split(b'\x00'):
        if not chunk:continue
        parts=chunk.split(b'\x14');timing=parts[0].split(b'\x15')
        if len(timing)!=2:continue
        o=float(timing[0]);d=float(timing[1])
        for text in parts[1:]:
            if text:onset.append(o);duration.append(d);labels.append(text.decode('utf-8'))
    return np.array(onset),np.array(duration),np.array(labels),start

if __name__=='__main__':
    import pyedflib
    n=0
    for p in Path(__file__).resolve().parents[1].joinpath('data').glob('*Hypnogram.edf'):
        a,b,c,d=read_hypnogram(p)
        try:
            with pyedflib.EdfReader(str(p)) as r:
                x,y,z=r.readAnnotations();dt=r.getStartdatetime()
            assert np.array_equal(a,x) and np.array_equal(b,y) and np.array_equal(c,z) and d==dt,p.name
            n+=1
        except OSError:pass
    print('Verified against pyedflib:',n)
