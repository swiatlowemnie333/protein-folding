import sys,glob,math,hashlib,pickle
import numpy as np
exec(open('/tmp/gly/gly_scan.py').read().split("seen=set()")[0].replace("D=sys.argv[1]","D=''"))
files=sorted(glob.glob('/tmp/all20/set/*.pdb'))
kept=[];rows=[];pid=0
def ident(a,b):
    if abs(len(a)-len(b))>3:return 0
    m=min(len(a),len(b));return sum(x==y for x,y in zip(a,b))/m
for f in files:
    try:r=parse(f)
    except Exception:continue
    if len(r)<40:continue
    seq=''.join(x['aa'] for x in r)
    if any(ident(seq,s)>0.8 for s in kept):continue
    ca=np.array([x['CA'] for x in r]);n=len(r)
    if np.isnan(ca).any():continue
    kept.append(seq);pid+=1
    ok=[np.linalg.norm(ca[i+1]-ca[i])<4.2 for i in range(n-1)]
    dm=np.linalg.norm(ca[:,None]-ca[None],axis=2);bur=(dm<10).sum(1)-1
    ss=['C']*n
    for i in range(n-4):
        if all(ok[i:i+4]) and 4.6<dm[i,i+3]<6.1 and 5.6<dm[i,i+4]<7.0:
            for j in range(i,i+5):ss[j]='H'
    for i in range(1,n-1):
        if ss[i]=='C' and ok[i-1] and ok[i] and dm[i-1,i+1]>6.2:ss[i]='E'
    for i in range(2,n-2):
        if not all(ok[i-2:i+2]):continue
        try:
            phi=dih(r[i-1]['C'],r[i]['N'],r[i]['CA'],r[i]['C']);psi=dih(r[i]['N'],r[i]['CA'],r[i]['C'],r[i+1]['N'])
        except Exception:continue
        rows.append(dict(aa=r[i]['aa'],phi=phi,psi=psi,d2=dm[i-1,i+1],d3=dm[i,i+3] if i+3<n else np.nan,d4=dm[i,i+4] if i+4<n else np.nan,
            bur=bur[i],ss=ss[i],pss=ss[i-1],nss=ss[i+1],p1=r[i-1]['aa'],n1=r[i+1]['aa'],p2=r[i-2]['aa'],n2=r[i+2]['aa'],rel=i/n,pid=pid,L=n))
print('proteins',pid,'residues',len(rows))
pickle.dump(rows,open('/tmp/all20/rows.pkl','wb'))
