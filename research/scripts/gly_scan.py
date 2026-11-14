import sys, glob, os, math, hashlib
import numpy as np
D=sys.argv[1]
three={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
def parse(f):
    res={};order=[];chain0=None
    for l in open(f):
        if l.startswith('ENDMDL'):break
        if not l.startswith('ATOM'):continue
        if l[16] not in ' A':continue
        ch=l[21]
        if chain0 is None:chain0=ch
        if ch!=chain0:continue
        rn=l[17:20]
        if rn not in three:continue
        k=(l[22:27])
        if k not in res:
            res[k]={'aa':three[rn]};order.append(k)
        res[k][l[12:16].strip()]=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
    return [res[k] for k in order if all(a in res[k] for a in ('N','CA','C'))]
def dih(p0,p1,p2,p3):
    b0=p0-p1;b1=p2-p1;b2=p3-p2
    b1n=b1/np.linalg.norm(b1)
    v=b0-np.dot(b0,b1n)*b1n;w=b2-np.dot(b2,b1n)*b1n
    x=np.dot(v,w);y=np.dot(np.cross(b1n,v),w)
    return math.degrees(math.atan2(y,x))
seen=set();rows=[];nprot=0
for f in sorted(glob.glob(D+'/*.pdb')):
    r=parse(f)
    if len(r)<30:continue
    seq=''.join(x['aa'] for x in r)
    h=hashlib.md5(seq.encode()).hexdigest()
    if h in seen:continue
    seen.add(h);nprot+=1
    ca=np.array([x['CA'] for x in r]);n=len(r)
    # contiguity
    ok=[np.linalg.norm(ca[i+1]-ca[i])<4.2 for i in range(n-1)]
    # burial: CA neighbours within 10A
    dm=np.linalg.norm(ca[:,None]-ca[None],axis=2)
    bur=(dm<10).sum(1)-1
    # SS from CA geometry
    ss=['C']*n
    for i in range(n-4):
        if all(ok[i:i+4]) and 4.6<dm[i,i+3]<6.1 and 5.6<dm[i,i+4]<7.0:
            for j in range(i,i+5):ss[j]='H'
    for i in range(1,n-1):
        if ss[i]=='C' and ok[i-1] and ok[i] and dm[i-1,i+1]>6.2: ss[i]='E'
    for i in range(2,n-2):
        if all(ok[i-2:i+2]):
            phi=psi=None
            if True:
                phi=dih(r[i-1]['C'],r[i]['N'],r[i]['CA'],r[i]['C'])
                psi=dih(r[i]['N'],r[i]['CA'],r[i]['C'],r[i+1]['N'])
            d2=dm[i-1,i+1]
            ang=math.degrees(math.acos(np.clip(np.dot(ca[i-1]-ca[i],ca[i+1]-ca[i])/(3.8*3.8),-1,1)))
            tors=dih(ca[i-2],ca[i-1],ca[i],ca[i+1])
            tors2=dih(ca[i-1],ca[i],ca[i+1],ca[i+2]) if i+2<n else np.nan
            rows.append((r[i]['aa'],phi,psi,d2,ang,tors,tors2,bur[i],ss[i],r[i-1]['aa'],r[i+1]['aa'],ss[i-1],ss[i+1],i/n,h[:6],i))
print('unique proteins',nprot,'residues',len(rows))
import pickle;pickle.dump(rows,open('/tmp/gly/rows.pkl','wb'))
aa=np.array([x[0] for x in rows]);phi=np.array([x[1] for x in rows]);psi=np.array([x[2] for x in rows])
d2=np.array([x[3] for x in rows]);ang=np.array([x[4] for x in rows]);bur=np.array([x[7] for x in rows]);ss=np.array([x[8] for x in rows])
G=aa=='G'
print('\nGly n=%d (%.1f%%)'%(G.sum(),100*G.mean()))
def st(name,a):print('%-34s Gly %.2f±%.2f | rest %.2f±%.2f'%(name,a[G].mean(),a[G].std(),a[~G].mean(),a[~G].std()))
st('d(i-1,i+1) CA-CA',d2);st('virtual CA angle',ang);st('burial (CA<10A)',bur.astype(float))
print('phi>0 fraction: Gly %.3f | rest %.3f'%((phi[G]>0).mean(),(phi[~G]>0).mean()))
for s in 'HEC':print('SS %s: Gly %.3f | rest %.3f'%(s,(ss[G]==s).mean(),(ss[~G]==s).mean()))
# Gly phi>0 per SS
for s in 'HEC':
    m=G&(ss==s);print('  Gly in %s: n=%d phi>0 %.2f  d2 mean %.2f sd %.2f'%(s,m.sum(),(phi[m]>0).mean(),d2[m].mean(),d2[m].std()))
# neighbours
from collections import Counter
prev=np.array([x[9] for x in rows]);nxt=np.array([x[10] for x in rows])
base=Counter(aa)
tot=len(aa)
print('\nWhat precedes/follows Gly (enrichment vs background):')
for nm,arr in (('prev',prev),('next',nxt)):
    c=Counter(arr[G]);out=[]
    for k,v in c.items():
        e=(v/G.sum())/(base[k]/tot);out.append((e,k,v))
    out.sort(reverse=True);print(nm,' '.join('%s%.2f(%d)'%(k,e,v) for e,k,v in out[:6]),'| lowest',' '.join('%s%.2f'%(k,e) for e,k,v in out[-4:]))
# Gly-Gly clusters
gg=G&((prev=='G')|(nxt=='G'));print('Gly with Gly neighbour: %.1f%% of Gly (expected %.1f%%)'%(100*gg.sum()/G.sum(),100*(1-(1-G.mean())**2)))
# phi>0 for Gly: what drives it?
m=G
X_phi=phi[m]>0
print('\nGly phi>0 rate by prev aa class and by neighbour-is-Gly:')
for cls,name in (('GNDST','small/polar'),('P','Pro'),('VILMFWY','hydrophobic/arom'),('AEQKRH','other')):
    pm=np.isin(prev[m],list(cls));nm=np.isin(nxt[m],list(cls))
    print('  prev %-17s n=%3d phi>0 %.2f | next n=%3d phi>0 %.2f'%(name,pm.sum(),X_phi[pm].mean() if pm.sum() else float('nan'),nm.sum(),X_phi[nm].mean() if nm.sum() else float('nan')))
for lo,hi in ((0,20),(20,26),(26,32),(32,99)):
    mm=(bur[m]>=lo)&(bur[m]<hi);print('  burial %d-%d n=%d phi>0 %.2f'%(lo,hi,mm.sum(),X_phi[mm].mean() if mm.sum() else float('nan')))
# Decision tree for Gly phi>0
from sklearn.tree import DecisionTreeClassifier,export_text
from sklearn.model_selection import cross_val_score
aas='ACDEFGHIKLMNPQRSTVWY'
def feat(i):
    x=rows[i];f=[aas.index(x[9]),aas.index(x[10]),x[7],x[13],'HEC'.index(x[11]),'HEC'.index(x[12]),x[3]]
    return f
idx=[i for i in range(len(rows)) if rows[i][0]=='G']
X=np.array([feat(i)[:6] for i in idx]);y=np.array([rows[i][1]>0 for i in idx])
names=['prev_aa','next_aa','burial','relpos','prevSS','nextSS']
for dpt in (2,3,4):
    t=DecisionTreeClassifier(max_depth=dpt,min_samples_leaf=15,random_state=0)
    print('tree depth',dpt,'CV acc %.3f'%cross_val_score(t,X,y,cv=5).mean(),'baseline %.3f'%max(y.mean(),1-y.mean()))
t=DecisionTreeClassifier(max_depth=3,min_samples_leaf=15,random_state=0).fit(X,y)
print(export_text(t,feature_names=names))
# Gly: predicting d(i-1,i+1) variance explained, vs other aa
from sklearn.tree import DecisionTreeRegressor
from sklearn.model_selection import cross_val_predict
for lab,mm in (('Gly',G),('non-Gly',~G)):
    ii=[i for i in range(len(rows)) if mm[i]][:6000]
    Xr=np.array([feat(i)[:6] for i in ii]);yr=np.array([rows[i][3] for i in ii])
    p=cross_val_predict(DecisionTreeRegressor(max_depth=5,min_samples_leaf=20,random_state=0),Xr,yr,cv=5)
    print('%s: SD d(i-1,i+1) %.2f -> residual after tree %.2f (R2 %.2f)'%(lab,yr.std(),(yr-p).std(),1-((yr-p)**2).sum()/((yr-yr.mean())**2).sum()))
