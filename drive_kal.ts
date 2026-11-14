import { readFileSync } from 'fs';
import { runMultiAgentFoldingToConvergence, predictSecondaryStructureGOR , setAngleCal} from './src/snaplockMultiAgentEngine';

const ONE1: Record<string,string> = {ALA:'A',ARG:'R',ASN:'N',ASP:'D',CYS:'C',GLN:'Q',GLU:'E',GLY:'G',HIS:'H',ILE:'I',LEU:'L',LYS:'K',MET:'M',PHE:'F',PRO:'P',SER:'S',THR:'T',TRP:'W',TYR:'Y',VAL:'V'};

function getSpecLocal(pdbPath: string, chain: string) {
  const txt = readFileSync(pdbPath,'utf8'); const lines:string[]=[]; const seq3:string[]=[]; const seen=new Set<string>();
  for (const l of txt.split('\n')) { if (l.startsWith('ENDMDL')) break;
    if (!l.startsWith('ATOM')||l[21]!==chain||!(l[16]===' '||l[16]==='A')) continue;
    const rn=l.slice(17,20), id=l.slice(22,27); if (!ONE1[rn]) continue; lines.push(l);
    if (l.slice(12,16).trim()==='CA'&&!seen.has(id)){seen.add(id);seq3.push(rn);} }
  return { sequence3Letter: seq3, referencePdbText: lines.join('\n')+'\nEND' };
}
function caCoords(pdb:string){ const xyz:number[][]=[]; const seen=new Set<string>(); let chain='';
  for (const l of (pdb||'').split('\n')){ if(!l.startsWith('ATOM')||l.slice(12,16).trim()!=='CA'||!(l[16]===' '||l[16]==='A'))continue;
    if(!chain)chain=l[21]; if(l[21]!==chain)continue; const id=l.slice(22,27); if(seen.has(id))continue; seen.add(id);
    xyz.push([+l.slice(30,38),+l.slice(38,46),+l.slice(46,54)]); } return xyz; }
function nativeContacts(pdb:string,n:number,minSep=6){ const ca=caCoords(pdb); const m=Math.min(n,ca.length); const out:{i:number;j:number}[]=[];
  for(let i=0;i<m;i++)for(let j=i+minSep;j<m;j++) if(Math.hypot(ca[i][0]-ca[j][0],ca[i][1]-ca[j][1],ca[i][2]-ca[j][2])<8) out.push({i,j}); return out; }
const gorSS=(seq3:string[])=>predictSecondaryStructureGOR(seq3).map((x:string)=>x[0]==='H'?'H':x[0]==='S'?'E':'C').join('');

function optimalRmsd(P:number[][],Q:number[][]):number{ const n=Math.min(P.length,Q.length); if(n<3)return NaN;
  const c=(A:number[][])=>[0,1,2].map(k=>A.slice(0,n).reduce((s,p)=>s+p[k],0)/n); const cp=c(P),cq=c(Q);
  const X=P.slice(0,n).map(p=>p.map((x,k)=>x-cp[k])),Y=Q.slice(0,n).map(p=>p.map((x,k)=>x-cq[k]));
  const S=[[0,0,0],[0,0,0],[0,0,0]]; for(let i=0;i<n;i++)for(let a=0;a<3;a++)for(let b=0;b<3;b++)S[a][b]+=X[i][a]*Y[i][b];
  const [[xx,xy,xz],[yx,yy,yz],[zx,zy,zz]]=S;
  const N=[[xx+yy+zz,yz-zy,zx-xz,xy-yx],[yz-zy,xx-yy-zz,xy+yx,zx+xz],[zx-xz,xy+yx,-xx+yy-zz,yz+zy],[xy-yx,zx+xz,yz+zy,-xx-yy+zz]];
  const e0=X.reduce((s,p)=>s+p[0]**2+p[1]**2+p[2]**2,0)+Y.reduce((s,p)=>s+p[0]**2+p[1]**2+p[2]**2,0);
  let q=[0.5,0.5,0.5,0.5]; for(let it=0;it<20000;it++){const w=N.map((r,i)=>r.reduce((s,x,j)=>s+x*q[j],0)+e0*q[i]);const no=Math.hypot(...w);q=w.map(x=>x/no);}
  const lam=N.map(r=>r.reduce((s,x,j)=>s+x*q[j],0)).reduce((s,x,j)=>s+x*q[j],0); return Math.sqrt(Math.max(0,(e0-2*lam)/n)); }

// odzysk natywnych kontaktow dalekiego zasiegu (arkusz beta = kontakty miedzy odleglymi odcinkami)
function lrContactRecovery(model:number[][], native:number[][], minSep:number){
  const n=Math.min(model.length,native.length); let tot=0,hit=0;
  for(let i=0;i<n;i++)for(let j=i+minSep;j<n;j++){
    const dn=Math.hypot(native[i][0]-native[j][0],native[i][1]-native[j][1],native[i][2]-native[j][2]);
    if(dn>=8)continue; tot++;
    const dm=Math.hypot(model[i][0]-model[j][0],model[i][1]-model[j][1],model[i][2]-model[j][2]);
    if(dm<8)hit++; }
  return {tot, recovery: tot?hit/tot:0}; }

function runOne(spec:any,key:string,cfg:string,seed:number,frames:number,mcSteps:number){
  setAngleCal(cfg.includes('K'));
  const run=runMultiAgentFoldingToConvergence({ sequenceCodes:spec.sequence3Letter, cdAlphaTargetPct:0, cdBetaTargetPct:0,
    explicitDisulfidePairs:[], referencePdbText:spec.referencePdbText, referenceId:key, enableTargetGuidedRefinement:false,
    randomSeed:seed, initialNoiseDeg:15, mjContactScale:cfg.includes('B')?1.0:0, hBondScale:cfg.includes('C')?1.5:0,
    secondaryStructurePredictor:cfg.includes('D')?'GOR':'LEGACY',
    monteCarloSteps:(cfg.includes('E')||cfg.includes('F')||cfg.includes('H')||cfg.includes('P')||cfg.includes('Z')||cfg.includes('R'))?mcSteps:0,
    zipperFraction:cfg.includes('Z')?0.05:0, ramaScale:cfg.includes('R')?Number(process.env.RS||1):0, betaScale:cfg.includes('P')?3:0, betaSS:cfg.includes('P')?gorSS(spec.sequence3Letter):undefined,
    burialScale:cfg.includes('H')?20:0, contactRestraints:cfg.includes('F')?nativeContacts(spec.referencePdbText,spec.sequence3Letter.length):undefined,
    contactWeight:2.0 } as any, frames, 3);
  const g:any=run.groundTruthValidation;
  const model=(run.caPositions as any[]).map((v:any)=>[v[0],v[1],v[2]]);
  return { rmsd:g?.caRmsdAngstrom??null, tm:g?.tmScore??null, model };
}
function d2err(model:number[][],nat:number[][],seq:string[]){ const d=(P:number[][],i:number,j:number)=>Math.hypot(P[i][0]-P[j][0],P[i][1]-P[j][1],P[i][2]-P[j][2]);
  let all=0,na=0,gl=0,ng=0; for(let i=1;i<nat.length-1&&i<model.length-1;i++){ const e=Math.abs(d(model,i-1,i+1)-d(nat,i-1,i+1)); all+=e;na++; if(seq[i]==='GLY'){gl+=e;ng++;} } return {all:all/Math.max(1,na),gly:ng?gl/ng:NaN}; }
const mean=(a:number[])=>a.reduce((x,y)=>x+y,0)/Math.max(1,a.length);
const sd=(a:number[])=>{const m=mean(a);return Math.sqrt(mean(a.map(x=>(x-m)**2)));};

const TARGETS=[ {key:'PROTEIN_G_56', path:'tests/fixtures/1PGA.pdb', chain:'A', label:'Bialko G B1 (1PGA, bb-a-bb)'},
                {key:'ENGRAILED_54', path:'tests/fixtures/1ENH.pdb', chain:'A', label:'Engrailed (1ENH, 3x alfa, kontrola)'} ];
const CONFIGS=(process.env.CFGS||'ACPZ,ACPZR').split(',');
const SEEDS=Number(process.env.SEEDS||5), FRAMES=150, MC=5000;

for (const tgt of TARGETS){
  const spec=getSpecLocal(tgt.path,tgt.chain);
  const nat=caCoords(spec.referencePdbText);
  const ss=gorSS(spec.sequence3Letter);
  const nStrand=(ss.match(/E/g)||[]).length;
  const lrNat=lrContactRecovery(nat,nat,20).tot; // ile natywnych kontaktow dalekich istnieje
  console.log(`\n==== ${tgt.label} | n=${spec.sequence3Letter.length} | GOR E-reszt=${nStrand} | nat. kontaktow dalekich(|i-j|>=20)=${lrNat} ====`);
  console.log('SS:', ss);
  console.log('cfg   RMSD(sr+/-sd)   best    TM      odzysk kontaktow dalekich |i-j|>=20');
  for (const cfg of CONFIGS){
    const rs:number[]=[],ts:number[]=[],recs:number[]=[],la:number[]=[],lg:number[]=[]; let best=1e9,bestModel:number[][]|null=null;
    for (let s=1;s<=SEEDS;s++){ const r=runOne(spec,tgt.key,cfg,s,FRAMES,MC);
      if(r.rmsd!=null){rs.push(r.rmsd);ts.push(r.tm??0); const rec=lrContactRecovery(r.model,nat,20); recs.push(rec.recovery); const de=d2err(r.model,nat,spec.sequence3Letter); la.push(de.all); if(!isNaN(de.gly))lg.push(de.gly);
        if(r.rmsd<best){best=r.rmsd;bestModel=r.model;}} }
    console.log(`${cfg.padEnd(5)} ${mean(rs).toFixed(2)}+/-${sd(rs).toFixed(2)}    ${best.toFixed(2)}   ${mean(ts).toFixed(3)}   ${(mean(recs)*100).toFixed(0)}%  | lokalny blad d(i-1,i+1): wszystkie ${mean(la).toFixed(3)} A, Gly ${mean(lg).toFixed(3)} A`);
  }
}
console.log('\n[gotowe]');
