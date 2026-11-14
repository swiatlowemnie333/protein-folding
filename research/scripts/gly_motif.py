import pickle,numpy as np
from collections import defaultdict
rows=pickle.load(open('/tmp/gly/rows.pkl','rb'))
# (aa,phi,psi,d2,ang,tors,tors2,bur,ss,prev,next,pss,nss,relpos,hash,i)
def stats(sel,label,allrows=rows):
    s=[x for x in allrows if sel(x)]
    if len(s)<8:return
    ph=np.array([x[1] for x in s]);d2=np.array([x[3] for x in s]);an=np.array([x[4] for x in s]);t=np.array([x[5] for x in s])
    print('%-44s n=%4d phi>0 %.2f | d2 %.2f±%.2f | angle %.0f±%.0f | vtors %+.0f±%.0f'%(label,len(s),(ph>0).mean(),d2.mean(),d2.std(),an.mean(),an.std(),t.mean(),t.std()))
print('== Gly by (prev SS, next SS) ==')
for p in 'HEC':
    for n in 'HEC':
        stats(lambda x:x[0]=='G' and x[11]==p and x[12]==n,'Gly prevSS=%s nextSS=%s'%(p,n))
print('\n== same for ALL non-Gly (reference) ==')
for p,n in (('H','C'),('E','E'),('C','C'),('E','C'),('C','E')):
    stats(lambda x:x[0]!='G' and x[11]==p and x[12]==n,'nonGly prevSS=%s nextSS=%s'%(p,n))
print('\n== motifs ==')
stats(lambda x:x[0]=='G' and x[9] in 'ND','Gly after N/D (turn motif N-G, D-G)')
stats(lambda x:x[0]=='G' and x[9]=='N','Gly after N')
stats(lambda x:x[0]=='G' and x[9]=='D','Gly after D')
stats(lambda x:x[0]=='G' and x[9] not in 'ND','Gly NOT after N/D')
stats(lambda x:x[0]=='G' and x[11]=='H' and x[12]!='H','Gly = C-cap of helix (prev H, next not H)')
stats(lambda x:x[0]=='G' and x[11]=='E' and x[12]=='E','Gly between two strand residues')
stats(lambda x:x[0]=='G' and x[11]!='H' and x[12]=='H','Gly = N-cap of helix')
stats(lambda x:x[0]=='G' and x[9]=='G' or x[0]=='G' and x[10]=='G','Gly next to Gly')
stats(lambda x:x[0]=='G' and x[9]=='P' or x[0]=='G' and x[10]=='P','Gly next to Pro')
# Does phi>0 Gly change neighbours? Compare d2 of residue i-1 and i+1 around Gly(phi>0) vs Gly(phi<0)
byprot=defaultdict(dict)
for x in rows:byprot[x[14]][x[15]]=x
def nb(cond,label):
    v=[]
    for h,d in byprot.items():
        for i,x in d.items():
            if x[0]=='G' and cond(x):
                for k in (i-1,i+1):
                    if k in d and d[k][0]!='G':v.append(d[k][3])
    if v:print('%-44s neighbours: d2 %.2f±%.2f n=%d'%(label,np.mean(v),np.std(v),len(v)))
nb(lambda x:x[1]>0,'neighbours of Gly with phi>0')
nb(lambda x:x[1]<0,'neighbours of Gly with phi<0')
# Within coil: is "Gly" alone predictive? information gain on phi>0 using aa identity vs SS
print('\n== how much of phi>0 is explained ==')
allphi=np.array([x[1]>0 for x in rows]);isG=np.array([x[0]=='G' for x in rows])
print('phi>0 overall %.3f ; among Gly %.3f ; among non-Gly %.3f'%(allphi.mean(),allphi[isG].mean(),allphi[~isG].mean()))
print('share of all phi>0 residues that are Gly: %.1f%% (Gly is %.1f%% of residues)'%(100*(allphi&isG).sum()/allphi.sum(),100*isG.mean()))
# non-Gly with phi>0: which aa?
from collections import Counter
c=Counter(x[0] for x in rows if x[0]!='G' and x[1]>0);t=Counter(x[0] for x in rows if x[0]!='G')
print('non-Gly phi>0 rate by aa:',' '.join('%s%.2f'%(k,c[k]/t[k]) for k in sorted(t,key=lambda k:-c[k]/t[k])[:7]))
# 1PGA actual
for x in rows:
    pass
