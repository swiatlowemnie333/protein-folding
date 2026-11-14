// Moduł dalekiego zasięgu: globalne osadzanie struktury z więzów odległości
// (majoryzacja stresu / zejście gradientowe), zamiast lokalnego MC po jednym kącie.
// Test na 1PGA: pełna macierz (sanity) / same kontakty / kontakty+SS+pary dalekie.
import { readFileSync } from 'fs';
function caCoords(pdb:string){const xyz:number[][]=[];const seen=new Set<string>();let chain='';
  for(const l of pdb.split('\n')){if(!l.startsWith('ATOM')||l.slice(12,16).trim()!=='CA'||!(l[16]===' '||l[16]==='A'))continue;if(!chain)chain=l[21];if(l[21]!==chain)continue;const id=l.slice(22,27);if(seen.has(id))continue;seen.add(id);xyz.push([+l.slice(30,38),+l.slice(38,46),+l.slice(46,54)]);}return xyz;}
const dist=(a:number[],b:number[])=>Math.hypot(a[0]-b[0],a[1]-b[1],a[2]-b[2]);
// RMSD z dopuszczeniem odbicia (geometria odległości nie ustala chiralności)
function optimalRmsd(P:number[][],Q:number[][]):number{const n=Math.min(P.length,Q.length);if(n<3)return NaN;const c=(A:number[][])=>[0,1,2].map(k=>A.slice(0,n).reduce((s,p)=>s+p[k],0)/n);const cp=c(P),cq=c(Q);const X=P.slice(0,n).map(p=>p.map((x,k)=>x-cp[k])),Y=Q.slice(0,n).map(p=>p.map((x,k)=>x-cq[k]));const S=[[0,0,0],[0,0,0],[0,0,0]];for(let i=0;i<n;i++)for(let a=0;a<3;a++)for(let b=0;b<3;b++)S[a][b]+=X[i][a]*Y[i][b];const[[xx,xy,xz],[yx,yy,yz],[zx,zy,zz]]=S;const N=[[xx+yy+zz,yz-zy,zx-xz,xy-yx],[yz-zy,xx-yy-zz,xy+yx,zx+xz],[zx-xz,xy+yx,-xx+yy-zz,yz+zy],[xy-yx,zx+xz,yz+zy,-xx-yy+zz]];const e0=X.reduce((s,p)=>s+p[0]**2+p[1]**2+p[2]**2,0)+Y.reduce((s,p)=>s+p[0]**2+p[1]**2+p[2]**2,0);let q=[0.5,0.5,0.5,0.5];for(let it=0;it<20000;it++){const w=N.map((r,i)=>r.reduce((s,x,j)=>s+x*q[j],0)+e0*q[i]);const no=Math.hypot(...w);q=w.map(x=>x/no);}const lam=N.map(r=>r.reduce((s,x,j)=>s+x*q[j],0)).reduce((s,x,j)=>s+x*q[j],0);return Math.sqrt(Math.max(0,(e0-2*lam)/n));}
function rmsdChiral(P:number[][],Q:number[][]){const m=P.map(p=>[-p[0],p[1],p[2]]);return Math.min(optimalRmsd(P,Q),optimalRmsd(m,Q));}
function lrRecovery(model:number[][],nat:number[][],sep:number){const n=Math.min(model.length,nat.length);let tot=0,hit=0;for(let i=0;i<n;i++)for(let j=i+sep;j<n;j++){if(dist(nat[i],nat[j])>=8)continue;tot++;if(dist(model[i],model[j])<8)hit++;}return tot?hit/tot:0;}
// RNG deterministyczny
function mulberry32(a:number){return function(){a|=0;a=a+0x6D2B79F5|0;let t=Math.imul(a^a>>>15,1|a);t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296;};}
// Zejście gradientowe na stresie: S = sum w_ij (|xi-xj| - d_ij)^2
function embed(n:number, targets:{i:number;j:number;d:number;w:number}[], seed:number, steps=300){
  const rnd=mulberry32(seed);const X:number[][]=Array.from({length:n},()=>[rnd()*20-10,rnd()*20-10,rnd()*20-10]);
  const Wi=new Array(n).fill(0);for(const t of targets){Wi[t.i]+=t.w;Wi[t.j]+=t.w;}
  for(let it=0;it<steps;it++){const Xn:number[][]=X.map(p=>[0,0,0]);
    for(const t of targets){const a=X[t.i],b=X[t.j];const dx=a[0]-b[0],dy=a[1]-b[1],dz=a[2]-b[2];const r=Math.hypot(dx,dy,dz)||1e-6;const s=t.d/r;
      // wklad i: x_j + d*(x_i-x_j)/r ; wklad j: x_i + d*(x_j-x_i)/r
      Xn[t.i][0]+=t.w*(b[0]+s*dx);Xn[t.i][1]+=t.w*(b[1]+s*dy);Xn[t.i][2]+=t.w*(b[2]+s*dz);
      Xn[t.j][0]+=t.w*(a[0]-s*dx);Xn[t.j][1]+=t.w*(a[1]-s*dy);Xn[t.j][2]+=t.w*(a[2]-s*dz);}
    for(let i=0;i<n;i++){if(Wi[i]>0){X[i][0]=Xn[i][0]/Wi[i];X[i][1]=Xn[i][1]/Wi[i];X[i][2]=Xn[i][2]/Wi[i];}}}
  return X;}
function stress(X:number[][],targets:{i:number;j:number;d:number;w:number}[]){let s=0;for(const t of targets){const r=dist(X[t.i],X[t.j]);s+=t.w*(r-t.d)**2;}return s;}
function bestEmbed(n:number,targets:any[],restarts=6,steps=4000){let bestX:number[][]|null=null,bestS=1e18;for(let s=1;s<=restarts;s++){const X=embed(n,targets,s*7+1,steps);const st=stress(X,targets);if(st<bestS){bestS=st;bestX=X;}}return bestX!;}

const nat=caCoords(readFileSync('tests/fixtures/1PGA.pdb','utf8'));const n=nat.length;
// natywne kontakty i SS (helisa 18-31 wg GOR -> pary i,i+3/i,i+4; nici -> lokalne)
const chain:any[]=[];for(let i=0;i<n-1;i++)chain.push({i,j:i+1,d:dist(nat[i],nat[i+1]),w:5});
const contacts6:any[]=[];for(let i=0;i<n;i++)for(let j=i+6;j<n;j++)if(dist(nat[i],nat[j])<8)contacts6.push({i,j,d:dist(nat[i],nat[j]),w:1});
const full:any[]=[];for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)full.push({i,j,d:dist(nat[i],nat[j]),w:1});
// SS lokalna: helisa z GOR (indeksy ~17-30), pary i,i+3 oraz i,i+4
const ssLocal:any[]=[];for(let i=17;i<=29;i++){if(i+3<n)ssLocal.push({i,j:i+3,d:dist(nat[i],nat[i+3]),w:3});if(i+4<n)ssLocal.push({i,j:i+4,d:dist(nat[i],nat[i+4]),w:3});}
// pary DALEKIE: tylko kontakty |i-j|>=20 (krzyżowe wiazania arkusza)
const farOnly:any[]=[];for(let i=0;i<n;i++)for(let j=i+20;j<n;j++)if(dist(nat[i],nat[j])<8)farOnly.push({i,j,d:dist(nat[i],nat[j]),w:2});

function report(name:string,X:number[][]){const r=rmsdChiral(X,nat);console.log(`${name.padEnd(42)} RMSD(z odbiciem)=${r.toFixed(2)} Å  odzysk|i-j|>=20=${(lrRecovery(X,nat,20)*100).toFixed(0)}%  |i-j|>=6=${(lrRecovery(X,nat,6)*100).toFixed(0)}%`);}
console.log(`1PGA n=${n}, kontaktow(|i-j|>=6)=${contacts6.length}, dalekich(|i-j|>=20)=${farOnly.length}\n`);
console.log('--- globalne osadzanie (zejscie gradientowe, 6 restartow) ---');
report('pelna macierz odleglosci (sanity)', bestEmbed(n,[...full],8,400));
report('lancuch + same kontakty (<8A, |i-j|>=6)', bestEmbed(n,[...chain,...contacts6],8,400));
report('lancuch + SS + pary DALEKIE(|i-j|>=20)', bestEmbed(n,[...chain,...ssLocal,...farOnly],8,400));
console.log('\n[dla porownania] silnik AF (lokalne MC, natywne kontakty): RMSD 7,36  odzysk|i-j|>=20 64%');
