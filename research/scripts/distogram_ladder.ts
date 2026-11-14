// Drabina informacji o odleglosciach na 1PGA: ile "gestosci" potrzeba, by zlozyc fald.
// binarny kontakt (tylko "<8A") -> distogram (odleglosci w koszykach 2A) -> pelna macierz.
import { readFileSync } from 'fs';
function caCoords(pdb:string){const xyz:number[][]=[];const seen=new Set<string>();let chain='';
  for(const l of pdb.split('\n')){if(!l.startsWith('ATOM')||l.slice(12,16).trim()!=='CA'||!(l[16]===' '||l[16]==='A'))continue;if(!chain)chain=l[21];if(l[21]!==chain)continue;const id=l.slice(22,27);if(seen.has(id))continue;seen.add(id);xyz.push([+l.slice(30,38),+l.slice(38,46),+l.slice(46,54)]);}return xyz;}
const dist=(a:number[],b:number[])=>Math.hypot(a[0]-b[0],a[1]-b[1],a[2]-b[2]);
function optimalRmsd(P:number[][],Q:number[][]):number{const n=Math.min(P.length,Q.length);if(n<3)return NaN;const c=(A:number[][])=>[0,1,2].map(k=>A.slice(0,n).reduce((s,p)=>s+p[k],0)/n);const cp=c(P),cq=c(Q);const X=P.slice(0,n).map(p=>p.map((x,k)=>x-cp[k])),Y=Q.slice(0,n).map(p=>p.map((x,k)=>x-cq[k]));const S=[[0,0,0],[0,0,0],[0,0,0]];for(let i=0;i<n;i++)for(let a=0;a<3;a++)for(let b=0;b<3;b++)S[a][b]+=X[i][a]*Y[i][b];const[[xx,xy,xz],[yx,yy,yz],[zx,zy,zz]]=S;const N=[[xx+yy+zz,yz-zy,zx-xz,xy-yx],[yz-zy,xx-yy-zz,xy+yx,zx+xz],[zx-xz,xy+yx,-xx+yy-zz,yz+zy],[xy-yx,zx+xz,yz+zy,-xx-yy+zz]];const e0=X.reduce((s,p)=>s+p[0]**2+p[1]**2+p[2]**2,0)+Y.reduce((s,p)=>s+p[0]**2+p[1]**2+p[2]**2,0);let q=[0.5,0.5,0.5,0.5];for(let it=0;it<20000;it++){const w=N.map((r,i)=>r.reduce((s,x,j)=>s+x*q[j],0)+e0*q[i]);const no=Math.hypot(...w);q=w.map(x=>x/no);}const lam=N.map(r=>r.reduce((s,x,j)=>s+x*q[j],0)).reduce((s,x,j)=>s+x*q[j],0);return Math.sqrt(Math.max(0,(e0-2*lam)/n));}
function rmsdChiral(P:number[][],Q:number[][]){const m=P.map(p=>[-p[0],p[1],p[2]]);return Math.min(optimalRmsd(P,Q),optimalRmsd(m,Q));}
function lrRec(model:number[][],nat:number[][],sep:number){const n=Math.min(model.length,nat.length);let tot=0,hit=0;for(let i=0;i<n;i++)for(let j=i+sep;j<n;j++){if(dist(nat[i],nat[j])>=8)continue;tot++;if(dist(model[i],model[j])<8)hit++;}return tot?hit/tot:0;}
function mulberry32(a:number){return function(){a|=0;a=a+0x6D2B79F5|0;let t=Math.imul(a^a>>>15,1|a);t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296;};}
function embed(n:number,targets:{i:number;j:number;d:number;w:number}[],seed:number,steps=400){
  const rnd=mulberry32(seed);const X:number[][]=Array.from({length:n},()=>[rnd()*20-10,rnd()*20-10,rnd()*20-10]);
  const Wi=new Array(n).fill(0);for(const t of targets){Wi[t.i]+=t.w;Wi[t.j]+=t.w;}
  for(let it=0;it<steps;it++){const Xn:number[][]=X.map(()=>[0,0,0]);
    for(const t of targets){const a=X[t.i],b=X[t.j];const dx=a[0]-b[0],dy=a[1]-b[1],dz=a[2]-b[2];const r=Math.hypot(dx,dy,dz)||1e-6;const s=t.d/r;
      Xn[t.i][0]+=t.w*(b[0]+s*dx);Xn[t.i][1]+=t.w*(b[1]+s*dy);Xn[t.i][2]+=t.w*(b[2]+s*dz);
      Xn[t.j][0]+=t.w*(a[0]-s*dx);Xn[t.j][1]+=t.w*(a[1]-s*dy);Xn[t.j][2]+=t.w*(a[2]-s*dz);}
    for(let i=0;i<n;i++)if(Wi[i]>0){X[i][0]=Xn[i][0]/Wi[i];X[i][1]=Xn[i][1]/Wi[i];X[i][2]=Xn[i][2]/Wi[i];}}
  return X;}
function stress(X:number[][],t:any[]){let s=0;for(const x of t){const r=dist(X[x.i],X[x.j]);s+=x.w*(r-x.d)**2;}return s;}
function best(n:number,t:any[],restarts=8){let bX:number[][]|null=null,bS=1e18;for(let s=1;s<=restarts;s++){const X=embed(n,t,s*7+1);const st=stress(X,t);if(st<bS){bS=st;bX=X;}}return bX!;}

const nat=caCoords(readFileSync('tests/fixtures/1PGA.pdb','utf8'));const n=nat.length;
const chain:any[]=[];for(let i=0;i<n-1;i++)chain.push({i,j:i+1,d:dist(nat[i],nat[i+1]),w:5});
// 1. binarny kontakt: wiemy TYLKO ze <8A (nie znamy odleglosci) -> cel 6.5A, luzno
const binary:any[]=[...chain];for(let i=0;i<n;i++)for(let j=i+6;j<n;j++)if(dist(nat[i],nat[j])<8)binary.push({i,j,d:6.5,w:1});
// 2. distogram: odleglosci w koszykach 2A dla par <16A (jak distogram AlphaFold, ale z wyroczni)
const bin2=(d:number)=>Math.round(d/2)*2;
const disto16:any[]=[...chain];for(let i=0;i<n;i++)for(let j=i+2;j<n;j++){const d=dist(nat[i],nat[j]);if(d<16)disto16.push({i,j,d:bin2(d),w:1});}
// 3. distogram szerszy <20A
const disto20:any[]=[...chain];for(let i=0;i<n;i++)for(let j=i+2;j<n;j++){const d=dist(nat[i],nat[j]);if(d<20)disto20.push({i,j,d:bin2(d),w:1});}
// 4. pelna macierz (ciagle)
const full:any[]=[];for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)full.push({i,j,d:dist(nat[i],nat[j]),w:1});

function rep(name:string,t:any[]){const X=best(n,t);console.log(`${name.padEnd(46)} par=${String(t.length).padStart(4)}  RMSD=${rmsdChiral(X,nat).toFixed(2)} Å  odzysk|i-j|>=20=${(lrRec(X,nat,20)*100).toFixed(0)}%`);}
console.log(`1PGA n=${n}\n`);
rep('1. binarny kontakt (tylko <8A, cel 6,5)', binary);
rep('2. DISTOGRAM koszyki 2A, pary <16A', disto16);
rep('3. DISTOGRAM koszyki 2A, pary <20A', disto20);
rep('4. pelna macierz ciagla (sanity)', full);
console.log('\n[porownanie] silnik AF lokalne MC: 7,36 Å / 64% ; embedding exact-kontakty: 9,6 Å / 100%');
