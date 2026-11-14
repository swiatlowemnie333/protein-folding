// Skaner korekcji bledu: nasza struktura vs natywna (1PGA). Gdzie siedzi blad,
// ile dalaby korekcja danego fragmentu, i z jaka struktura drugorzedowa sie wiaze.
import { readFileSync } from 'fs';
import { runMultiAgentFoldingToConvergence, predictSecondaryStructureGOR } from './src/snaplockMultiAgentEngine';
const ONE1:Record<string,string>={ALA:'A',ARG:'R',ASN:'N',ASP:'D',CYS:'C',GLN:'Q',GLU:'E',GLY:'G',HIS:'H',ILE:'I',LEU:'L',LYS:'K',MET:'M',PHE:'F',PRO:'P',SER:'S',THR:'T',TRP:'W',TYR:'Y',VAL:'V'};
function getSpec(p:string,ch:string){const txt=readFileSync(p,'utf8');const lines:string[]=[];const seq3:string[]=[];const seen=new Set<string>();for(const l of txt.split('\n')){if(l.startsWith('ENDMDL'))break;if(!l.startsWith('ATOM')||l[21]!==ch||!(l[16]===' '||l[16]==='A'))continue;const rn=l.slice(17,20),id=l.slice(22,27);if(!ONE1[rn])continue;lines.push(l);if(l.slice(12,16).trim()==='CA'&&!seen.has(id)){seen.add(id);seq3.push(rn);}}return{seq3,ref:lines.join('\n')+'\nEND'};}
function caCoords(pdb:string){const xyz:number[][]=[];const seen=new Set<string>();let chain='';for(const l of pdb.split('\n')){if(!l.startsWith('ATOM')||l.slice(12,16).trim()!=='CA'||!(l[16]===' '||l[16]==='A'))continue;if(!chain)chain=l[21];if(l[21]!==chain)continue;const id=l.slice(22,27);if(seen.has(id))continue;seen.add(id);xyz.push([+l.slice(30,38),+l.slice(38,46),+l.slice(46,54)]);}return xyz;}
const nativeContacts=(pdb:string,n:number)=>{const ca=caCoords(pdb);const m=Math.min(n,ca.length);const o:any[]=[];for(let i=0;i<m;i++)for(let j=i+6;j<m;j++)if(Math.hypot(ca[i][0]-ca[j][0],ca[i][1]-ca[j][1],ca[i][2]-ca[j][2])<8)o.push({i,j});return o;};
// Kabsch zwracajacy nalozony model
function kabsch(P:number[][],Q:number[][]){const n=Math.min(P.length,Q.length);const c=(A:number[][])=>[0,1,2].map(k=>A.slice(0,n).reduce((s,p)=>s+p[k],0)/n);const cp=c(P),cq=c(Q);const X=P.slice(0,n).map(p=>p.map((x,k)=>x-cp[k])),Y=Q.slice(0,n).map(p=>p.map((x,k)=>x-cq[k]));const S=[[0,0,0],[0,0,0],[0,0,0]];for(let i=0;i<n;i++)for(let a=0;a<3;a++)for(let b=0;b<3;b++)S[a][b]+=X[i][a]*Y[i][b];const[[xx,xy,xz],[yx,yy,yz],[zx,zy,zz]]=S;const N=[[xx+yy+zz,yz-zy,zx-xz,xy-yx],[yz-zy,xx-yy-zz,xy+yx,zx+xz],[zx-xz,xy+yx,-xx+yy-zz,yz+zy],[xy-yx,zx+xz,yz+zy,-xx-yy+zz]];let q=[1,0,0,0];for(let it=0;it<30000;it++){const w=N.map((r,i)=>r.reduce((s,x,j)=>s+x*q[j],0)+50*q[i]);const no=Math.hypot(...w);q=w.map(x=>x/no);}const[w0,x,y,z]=q;const R=[[1-2*(y*y+z*z),2*(x*y-w0*z),2*(x*z+w0*y)],[2*(x*y+w0*z),1-2*(x*x+z*z),2*(y*z-w0*x)],[2*(x*z-w0*y),2*(y*z+w0*x),1-2*(x*x+y*y)]];const A=X.map(p=>[R[0][0]*p[0]+R[0][1]*p[1]+R[0][2]*p[2]+cq[0],R[1][0]*p[0]+R[1][1]*p[1]+R[1][2]*p[2]+cq[1],R[2][0]*p[0]+R[2][1]*p[1]+R[2][2]*p[2]+cq[2]]);return A;} // model nalozony na Q(natyw centr. wraca do cq)
const gorSS=(s:string[])=>predictSecondaryStructureGOR(s).map((x:string)=>x[0]==='H'?'H':x[0]==='S'?'E':'C').join('');
const spec=getSpec('tests/fixtures/1PGA.pdb','A');const nat=caCoords(spec.ref);const n=nat.length;const ss=gorSS(spec.seq3);
// nasz najlepszy model (AF = z kontaktami, kilka startow, bierzemy najlepszy)
let best:number[][]|null=null,bestR=1e9;
for(let s=1;s<=6;s++){const run=runMultiAgentFoldingToConvergence({sequenceCodes:spec.seq3,cdAlphaTargetPct:0,cdBetaTargetPct:0,explicitDisulfidePairs:[],referencePdbText:spec.ref,referenceId:'1PGA',enableTargetGuidedRefinement:false,randomSeed:s,initialNoiseDeg:15,hBondScale:1.5,secondaryStructurePredictor:'GOR',monteCarloSteps:5000,betaScale:3,betaSS:ss,contactRestraints:nativeContacts(spec.ref,n),contactWeight:2.0}as any,150,3);const g:any=run.groundTruthValidation;const m=(run.caPositions as any[]).map(v=>[v[0],v[1],v[2]]);if(g.caRmsdAngstrom<bestR){bestR=g.caRmsdAngstrom;best=m;}}
const model=kabsch(best!,nat);
const dev=model.map((p,i)=>Math.hypot(p[0]-nat[i][0],p[1]-nat[i][1],p[2]-nat[i][2]));
const totSq=dev.reduce((a,d)=>a+d*d,0);
console.log(`1PGA — nasz najlepszy model (z kontaktami): RMSD ${Math.sqrt(totSq/n).toFixed(2)} Å, n=${n}`);
console.log(`\nBLAD WG STRUKTURY DRUGORZEDOWEJ (gdzie siedzi):`);
for(const t of['H','E','C']){const idx=[...Array(n).keys()].filter(i=>ss[i]===t);if(!idx.length)continue;const rms=Math.sqrt(idx.reduce((a,i)=>a+dev[i]*dev[i],0)/idx.length);const share=100*idx.reduce((a,i)=>a+dev[i]*dev[i],0)/totSq;console.log(`  ${t==='H'?'helisa':t==='E'?'nic   ':'zwoj  '} reszt ${String(idx.length).padStart(2)}  blad RMS ${rms.toFixed(1).padStart(4)} Å  udzial w bledzie ${share.toFixed(0)}%`);}
console.log(`\nPRZELOT PO FRAGMENTACH (okno 8): „gdyby ten fragment poprawic do natywnego -> RMSD spadnie do"`);
const W=8;const rows:any[]=[];
for(let st=0;st+W<=n;st+=2){let sq=0;for(let i=0;i<n;i++){if(i>=st&&i<st+W)continue;sq+=dev[i]*dev[i];}const rmsdFixed=Math.sqrt(sq/n);rows.push({st,end:st+W-1,rmsdFixed,ssFrag:ss.slice(st,st+W)});}
rows.sort((a,b)=>a.rmsdFixed-b.rmsdFixed);
console.log(`  (pelny RMSD ${Math.sqrt(totSq/n).toFixed(2)} Å; im nizej po poprawie, tym wiekszy blad tam siedzi)`);
rows.slice(0,6).forEach(r=>console.log(`  reszty ${String(r.st+1).padStart(2)}-${String(r.end+1).padStart(2)} [${r.ssFrag}] -> ${r.rmsdFixed.toFixed(2)} Å`));
console.log(`\nNAJGORSZE POJEDYNCZE RESZTY (odchylenie Å):`);
[...dev.keys()].sort((a,b)=>dev[b]-dev[a]).slice(0,8).forEach(i=>console.log(`  reszta ${String(i+1).padStart(2)} (${spec.seq3[i]}, ${ss[i]}): ${dev[i].toFixed(1)} Å`));
