"""Opiekun reguł, warstwa 3: geometria szkieletu (stałe silnika vs pomiar na PDB). python3 rules_keeper_geom.py DIR"""
import sys, os, glob, numpy as np
d = sys.argv[1]
def ang(a,b,c):
    u=a-b; v=c-b; return np.degrees(np.arccos(np.dot(u,v)/np.linalg.norm(u)/np.linalg.norm(v)))
def dih(p0,p1,p2,p3):
    b0=p0-p1; b1=p2-p1; b2=p3-p2; b1/=np.linalg.norm(b1)
    v=b0-np.dot(b0,b1)*b1; w=b2-np.dot(b2,b1)*b1
    return np.degrees(np.arctan2(np.dot(np.cross(b1,v),w),np.dot(v,w)))
S={k:[] for k in ['N-CA','CA-C','C-N','CA-CA','NCAC','NCAC_preP','CACN','CNCA','CNCA_P','omega','omega_cis','NCAC_G']}
for f in glob.glob(d+'/*.pdb'):
    res={}; order=[]
    for l in open(f):
        if l.startswith('ATOM') and l[12:16].strip() in ('N','CA','C') and l[16] in ' A':
            k=(l[21],l[22:27]); 
            if k not in res: res[k]={'aa':l[17:20]}; order.append(k)
            res[k][l[12:16].strip()]=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
    for a,b in zip(order,order[1:]):
        A,B=res[a],res[b]
        if not all(x in A for x in ('N','CA','C')) or not all(x in B for x in ('N','CA','C')): continue
        if np.linalg.norm(A['C']-B['N'])>1.8: continue
        S['N-CA'].append(np.linalg.norm(A['N']-A['CA'])); S['CA-C'].append(np.linalg.norm(A['CA']-A['C']))
        S['C-N'].append(np.linalg.norm(A['C']-B['N'])); S['CA-CA'].append(np.linalg.norm(A['CA']-B['CA']))
        t=ang(A['N'],A['CA'],A['C']); S['NCAC'].append(t)
        if A['aa']=='GLY': S['NCAC_G'].append(t)
        if B['aa']=='PRO': S['NCAC_preP'].append(t)
        S['CACN'].append(ang(A['CA'],A['C'],B['N']))
        c=ang(A['C'],B['N'],B['CA']); S['CNCA'].append(c)
        if B['aa']=='PRO': S['CNCA_P'].append(c)
        om=dih(A['CA'],A['C'],B['N'],B['CA']); S['omega'].append(abs(om) if abs(om)>90 else 180-abs(om))
        if abs(om)<30: S['omega_cis'].append(1)
ref={'N-CA':1.458,'CA-C':1.525,'C-N':1.329,'CA-CA':3.80,'NCAC':111.2,'NCAC_preP':110.2,'CACN':116.2,'CNCA':121.7}
print(f"{'wielkość':10} {'n':>7} {'średnia':>8} {'mediana':>8} {'SD':>6}  silnik  rozbieżność")
for k,v in S.items():
    if k=='omega_cis': print(f"cis-peptydy: {len(v)} z {len(S['omega'])} = {100*len(v)/len(S['omega']):.2f}%"); continue
    v=np.array(v); r=ref.get(k)
    flag='' if r is None else ('  <== PODEJRZANE' if abs(np.median(v)-r)>2*v.std() or abs(np.median(v)-r)>(0.02 if k in('N-CA','CA-C','C-N') else 1.5) else '  ok')
    print(f"{k:10} {len(v):7d} {v.mean():8.3f} {np.median(v):8.3f} {v.std():6.3f}  {r if r else '-':>6}{flag}")
