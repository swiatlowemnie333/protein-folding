import pickle,numpy as np,math
from collections import Counter,defaultdict
rows=pickle.load(open('/tmp/all20/rows.pkl','rb'))
N=len(rows)
def reg(x):
    p,s=x['phi'],x['psi']
    if p>0:return 'L'
    if -160<p<-20 and -120<s<50:return 'A'
    if p<-45 and (s>=90 or s<=-150):return 'B'
    return 'O'
for x in rows:x['reg']=reg(x)
aas='ACDEFGHIKLMNPQRSTVWY'
KD=dict(zip('ARNDCQEGHILKMFPSTWYV',[1.8,-4.5,-3.5,-3.5,2.5,-3.5,-3.5,-0.4,-3.2,4.5,3.8,-3.9,1.9,2.8,-1.6,-0.8,-0.7,-0.9,-1.3,4.2]))
VOL=dict(zip('ARNDCQEGHILKMFPSTWYV',[88.6,173.4,114.1,111.1,108.5,143.8,138.4,60.1,153.2,166.7,166.7,168.6,162.9,189.9,112.7,89.0,116.1,227.8,193.6,140.0]))
CHG=defaultdict(float,{'K':1,'R':1,'D':-1,'E':-1})
print('residues',N,'proteins',len(set(x['pid'] for x in rows)))
print('\n=== 1. Ramachandran-region mix per amino acid (A=alpha B=beta L=left(phi>0) O=other) ===')
print('aa    n   A    B    L    O   | mean d(i-1,i+1)  sd')
base={}
for a in aas:
    s=[x for x in rows if x['aa']==a];c=Counter(x['reg'] for x in s);n=len(s)
    d=np.array([x['d2'] for x in s]);base[a]=c
    print('%s %5d %.2f %.2f %.2f %.2f | %.2f  %.2f'%(a,n,c['A']/n,c['B']/n,c['L']/n,c['O']/n,d.mean(),d.std()))
# 2. group CV: how much does context explain, per aa
from sklearn.ensemble import RandomForestClassifier,RandomForestRegressor
from sklearn.model_selection import GroupKFold,cross_val_predict
def fe(x,mode):
    f=[x['bur'],x['rel'],'HEC'.index(x['pss']),'HEC'.index(x['nss'])]
    if mode>=1:
        for k in ('p1','n1','p2','n2'):f+= [KD[x[k]],VOL[x[k]],CHG[x[k]],x[k]=='G',x[k]=='P']
    return f
print('\n=== 2. predictability of Ramachandran region (leave-protein-out, random forest) ===')
print('aa  majority | context-acc | gain   | top variable')
names=['burial','relpos','prevSS','nextSS']+[k+'_'+q for k in ('p1','n1','p2','n2') for q in ('hyd','vol','chg','isG','isP')]
res={}
for a in aas:
    s=[x for x in rows if x['aa']==a]
    if len(s)<200:continue
    y=np.array([x['reg'] for x in s]);g=np.array([x['pid'] for x in s])
    X=np.array([fe(x,1) for x in s],float)
    maj=Counter(y).most_common(1)[0][1]/len(y)
    rf=RandomForestClassifier(150,min_samples_leaf=5,random_state=0,n_jobs=-1)
    p=cross_val_predict(rf,X,y,cv=GroupKFold(5),groups=g)
    acc=(p==y).mean();rf.fit(X,y);imp=rf.feature_importances_;top=np.argsort(-imp)[:2]
    res[a]=(maj,acc)
    print('%s   %.3f  |  %.3f      | %+.3f | %s(%.2f) %s(%.2f)'%(a,maj,acc,acc-maj,names[top[0]],imp[top[0]],names[top[1]],imp[top[1]]))
# 3. hunt strongest single conditions per aa (z-score of region rate vs aa baseline)
print('\n=== 3. strongest single conditions per amino acid (|z|>=4, n>=40) ===')
def conds(x):
    c=[('prevSS=%s'%x['pss']),('nextSS=%s'%x['nss']),('SS=%s'%x['ss']),
       ('bur<%d'%(10 if x['bur']<12 else 16 if x['bur']<16 else 20 if x['bur']<20 else 99)),
       ('prev=%s'%x['p1']),('next=%s'%x['n1']),('prev2=%s'%x['p2']),('next2=%s'%x['n2']),
       ('pos=%s'%('N-term' if x['rel']<0.1 else 'C-term' if x['rel']>0.9 else 'mid'))]
    c.append('prev_isHydrophobic=%s'%(KD[x['p1']]>1.5));c.append('next_isHydrophobic=%s'%(KD[x['n1']]>1.5))
    c.append('prev_charged=%s'%(CHG[x['p1']]!=0));c.append('next_charged=%s'%(CHG[x['n1']]!=0))
    return c
found=[]
for a in aas:
    s=[x for x in rows if x['aa']==a]
    for target in ('A','B','L'):
        pb=sum(1 for x in s if x['reg']==target)/len(s)
        cnt=defaultdict(lambda:[0,0])
        for x in s:
            for c in conds(x):
                cnt[c][0]+=1;cnt[c][1]+= x['reg']==target
        for c,(n,k) in cnt.items():
            if n<40:continue
            se=math.sqrt(max(pb*(1-pb),1e-4)/n);z=(k/n-pb)/se
            if abs(z)>=4:found.append((abs(z),a,target,c,n,k/n,pb))
found.sort(reverse=True)
byaa=defaultdict(list)
for f in found:byaa[f[1]].append(f)
for a in aas:
    for f in byaa[a][:3]:
        print('%s  region %s | %-22s n=%4d rate %.2f vs baseline %.2f (z=%.1f)'%(f[1],f[2],f[3],f[4],f[5],f[6],f[0]))
print('total strong conditions:',len(found),'(multiple tests -> expect some false positives; confirm on new data)')
# 4. split-half replication of conditions
import random
random.seed(1);pids=sorted(set(x['pid'] for x in rows));random.shuffle(pids);half=set(pids[:len(pids)//2])
def cond_rates(sub):
    out={}
    for a in aas:
        s=[x for x in sub if x['aa']==a]
        for target in 'ABL':
            cnt=defaultdict(lambda:[0,0])
            for x in s:
                for c in conds(x):cnt[c][0]+=1;cnt[c][1]+=x['reg']==target
            for c,(n,k) in cnt.items():out[(a,target,c)]=(n,k/n if n else 0)
    return out
h1=cond_rates([x for x in rows if x['pid'] in half]);h2=cond_rates([x for x in rows if x['pid'] not in half])
rep=0;tot=0
for z,a,t,c,n,r,pb in found:
    k=(a,t,c)
    if k in h1 and k in h2 and h1[k][0]>=15 and h2[k][0]>=15:
        tot+=1;d1=h1[k][1]-pb;d2_=h2[k][1]-pb
        if d1*d2_>0:rep+=1
print('replication: %d of %d strong conditions keep the SAME direction in both independent halves'%(rep,tot))
