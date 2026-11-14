# comparator_loop.py - petla komparatora: krzywa uczenia + aktywny dobor bialek + polowanie na nieuwzglednione zmienne.
# Uzycie: python3 -I comparator_loop.py <katalog_z_pdb> [max_bialek]
# Cel NIE jest "blad = 0" (to byloby zapamietanie), tylko: blad na NIEWIDZIANYCH bialkach vs liczba bialek uczacych.
import sys,glob,math,hashlib,random,time
import numpy as np
from collections import Counter
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import log_loss
D=sys.argv[1];MAXP=int(sys.argv[2]) if len(sys.argv)>2 else 10**9
three={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
AAS='ACDEFGHIKLMNPQRSTVWY';IDX={a:i for i,a in enumerate(AAS)}
KD=dict(zip('ARNDCQEGHILKMFPSTWYV',[1.8,-4.5,-3.5,-3.5,2.5,-3.5,-3.5,-0.4,-3.2,4.5,3.8,-3.9,1.9,2.8,-1.6,-0.8,-0.7,-0.9,-1.3,4.2]))
def parse(f):
    res={};order=[];ch0=None
    for l in open(f):
        if l.startswith('ENDMDL'):break
        if not l.startswith('ATOM') or l[16] not in ' A':continue
        if ch0 is None:ch0=l[21]
        if l[21]!=ch0:continue
        rn=l[17:20]
        if rn not in three:continue
        k=l[22:27]
        if k not in res:res[k]={'aa':three[rn]};order.append(k)
        res[k][l[12:16].strip()]=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
    return [res[k] for k in order if all(a in res[k] for a in('N','CA','C'))]
def dih(p0,p1,p2,p3):
    b0=p0-p1;b1=p2-p1;b2=p3-p2;b1n=b1/np.linalg.norm(b1)
    v=b0-np.dot(b0,b1n)*b1n;w=b2-np.dot(b2,b1n)*b1n
    return math.degrees(math.atan2(np.dot(np.cross(b1n,v),w),np.dot(v,w)))
def reg(p,s):
    if p>0:return 3 if False else 2
    if -160<p<-20 and -120<s<50:return 0
    if p<-45 and (s>=90 or s<=-150):return 1
    return 3
def ident(a,b):
    if abs(len(a)-len(b))>3:return 0
    m=min(len(a),len(b));return sum(x==y for x,y in zip(a,b))/m
t0=time.time();prots=[];kept=[]
for f in sorted(glob.glob(D+'/**/*.pdb',recursive=True)):
    if len(prots)>=MAXP:break
    try:r=parse(f)
    except Exception:continue
    if len(r)<40:continue
    seq=''.join(x['aa'] for x in r)
    if any(ident(seq,s)>0.8 for s in kept):continue
    ca=np.array([x['CA'] for x in r]);n=len(r)
    ok=[np.linalg.norm(ca[i+1]-ca[i])<4.2 for i in range(n-1)]
    y=np.full(n,-1)
    for i in range(1,n-1):
        if ok[i-1] and ok[i]:
            try:y[i]=reg(dih(r[i-1]['C'],r[i]['N'],r[i]['CA'],r[i]['C']),dih(r[i]['N'],r[i]['CA'],r[i]['C'],r[i+1]['N']))
            except Exception:pass
    kept.append(seq);prots.append((seq,y))
print('bialek uzytych: %d (po usunieciu >80%% identycznych); wczytanie %.0fs'%(len(prots),time.time()-t0))
W=7
def feats(seq,y,groups):
    n=len(seq);rows=[];ys=[]
    hf=sum(KD[a]>1.5 for a in seq)/n;cf=sum(a in 'KRDE' for a in seq)/n
    for i in range(n):
        if y[i]<0:continue
        f=[]
        for k in range(-W,W+1):
            j=i+k;f.append(IDX[seq[j]] if 0<=j<n else 20)
        g=[]
        g+=[f[W]]                                   # self
        if 'nb1' in groups:g+=[f[W-1],f[W+1]]
        if 'nb3' in groups:g+=[f[W-2],f[W+2],f[W-3],f[W+3]]
        if 'nb7' in groups:g+=[f[W-4],f[W+4],f[W-5],f[W+5],f[W-6],f[W+6],f[W-7],f[W+7]]
        if 'pos' in groups:g+=[i/n,i,n-i]
        if 'glob' in groups:g+=[n,hf,cf]
        rows.append(g);ys.append(y[i])
    return rows,ys
def build(idx,groups):
    X=[];Y=[];P=[]
    for p in idx:
        r,y=feats(prots[p][0],prots[p][1],groups);X+=r;Y+=y;P+=[p]*len(y)
    return np.array(X,float),np.array(Y),np.array(P)
def ncat(groups):
    c=1+(2 if 'nb1' in groups else 0)+(4 if 'nb3' in groups else 0)+(8 if 'nb7' in groups else 0);return c
def fit(Xtr,Ytr,groups):
    m=np.zeros(Xtr.shape[1],bool);m[:ncat(groups)]=True
    return HistGradientBoostingClassifier(max_iter=120,learning_rate=0.08,max_leaf_nodes=15,categorical_features=m,random_state=0).fit(Xtr,Ytr)
def score(m,X,Y):
    p=m.predict_proba(X);return (m.predict(X)==Y).mean(),log_loss(Y,p,labels=[0,1,2,3])
random.seed(0);allp=list(range(len(prots)));random.shuffle(allp)
nt=max(10,len(allp)//5);test=allp[:nt];pool=allp[nt:]
print('kontrolny (niewidziany) zbior: %d bialek; pula uczaca: %d'%(len(test),len(pool)))
G0=['self'];G1=['nb1','nb3']
Xte,Yte,_=build(test,G1)
base=Counter(Yte).most_common(1)[0][1]/len(Yte)
print('baza (najczestsza klasa) acc=%.3f\n'%base)
sizes=[s for s in (5,10,20,40,80,160,320,640,1280) if s<=len(pool)]
if sizes[-1]!=len(pool):sizes.append(len(pool))
def curve(order_name):
    out=[];cur=[]
    for s in sizes:
        if order_name=='random':
            cur=pool[:s]
        else:
            if not cur:cur=pool[:sizes[0]]
            else:
                Xc,Yc,Pc=build(cur,G1);m=fit(Xc,Yc,G1)
                rest=[p for p in pool if p not in cur]
                ls=[]
                for p in rest:
                    Xp,Yp,_=build([p],G1)
                    if len(Yp)==0:ls.append((-1,p));continue
                    ls.append((score(m,Xp,Yp)[1],p))
                ls.sort(reverse=True);cur=cur+[p for _,p in ls[:s-len(cur)]]
        Xtr,Ytr,_=build(cur,G1);m=fit(Xtr,Ytr,G1);a,l=score(m,Xte,Yte);out.append((s,len(Ytr),a,l))
    return out
print('=== KRZYWA UCZENIA: blad na niewidzianych bialkach vs liczba bialek uczacych ===')
print('bialek  reszt   |  losowy dobor: acc  logloss  |  aktywny dobor (najtrudniejsze): acc  logloss')
import os
cr=curve('random');ca_=curve('active') if not os.environ.get('NOACTIVE') else None
for k,a in enumerate(cr):
    b=ca_[k] if ca_ else (0,0,float('nan'),float('nan'))
    print('%5d %7d   |   %.3f  %.3f               |   %.3f  %.3f'%(a[0],a[1],a[2],a[3],b[2],b[3]))
print('\n=== POLOWANIE NA NIEUWZGLEDNIONE ZMIENNE (przyrost na niewidzianych, pelna pula) ===')
order=[('self',['self']),('+ sasiedzi +-1',['nb1']),('+ sasiedzi +-2,3',['nb1','nb3']),('+ sasiedzi +-4..7 (periodycznosc helisy)',['nb1','nb3','nb7']),('+ pozycja w lancuchu',['nb1','nb3','nb7','pos']),('+ cechy globalne bialka (dlugosc, hydrofobowosc, ladunek)',['nb1','nb3','nb7','pos','glob'])]
prev=None
for name,gr in order:
    Xtr,Ytr,_=build(pool,gr);Xt,Yt,_=build(test,gr);m=fit(Xtr,Ytr,gr);a,l=score(m,Xt,Yt)
    print('%-62s acc=%.3f logloss=%.3f %s'%(name,a,l,'' if prev is None else '(zmiana logloss %+.3f)'%(l-prev)));prev=l
print('\nczas: %.0fs'%(time.time()-t0))
