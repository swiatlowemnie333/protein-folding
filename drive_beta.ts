// Headless driver ablacji SNAPLOCK dla 1PGA (β+α) i 1ENH (α, kontrola).
// Odtwarza logikę AblationLabView.runGrid, ale bez przeglądarki. Silnik liczy RMSD sam (groundTruthValidation).
import { readFileSync } from 'fs';
import { runMultiAgentFoldingToConvergence, predictSecondaryStructureGOR } from './src/snaplockMultiAgentEngine';

const ONE1: Record<string,string> = {ALA:'A',ARG:'R',ASN:'N',ASP:'D',CYS:'C',GLN:'Q',GLU:'E',GLY:'G',HIS:'H',ILE:'I',LEU:'L',LYS:'K',MET:'M',PHE:'F',PRO:'P',SER:'S',THR:'T',TRP:'W',TYR:'Y',VAL:'V'};

function getSpecLocal(pdbPath: string, chain: string) {
  const txt = readFileSync(pdbPath, 'utf8');
  const lines: string[] = []; const seq3: string[] = []; const seen = new Set<string>();
  for (const l of txt.split('\n')) {
    if (l.startsWith('ENDMDL')) break;
    if (!l.startsWith('ATOM') || l[21] !== chain || !(l[16] === ' ' || l[16] === 'A')) continue;
    const rn = l.slice(17,20), id = l.slice(22,27);
    if (!ONE1[rn]) continue;
    lines.push(l);
    if (l.slice(12,16).trim() === 'CA' && !seen.has(id)) { seen.add(id); seq3.push(rn); }
  }
  return { sequence3Letter: seq3, cdAlphaTargetPct: 0, cdBetaTargetPct: 0, explicitDisulfidePairs0Based: [] as [number,number][], referencePdbText: lines.join('\n') + '\nEND' };
}
function nativeContacts(pdb: string, n: number) {
  const ca: number[][] = []; const seen = new Set<number>();
  for (const l of (pdb||'').split('\n')) {
    if (!l.startsWith('ATOM') || l.slice(12,16).trim() !== 'CA' || !(l[16] === ' ' || l[16] === 'A')) continue;
    const r = +l.slice(22,26); if (seen.has(r)) continue; seen.add(r);
    ca.push([+l.slice(30,38), +l.slice(38,46), +l.slice(46,54)]);
  }
  const m = Math.min(n, ca.length), out: {i:number;j:number}[] = [];
  for (let i=0;i<m;i++) for (let j=i+6;j<m;j++) if (Math.hypot(ca[i][0]-ca[j][0],ca[i][1]-ca[j][1],ca[i][2]-ca[j][2])<8) out.push({i,j});
  return out;
}
const gorSS = (seq3: string[]): string => predictSecondaryStructureGOR(seq3).map((x:string)=> x[0]==='H'?'H':x[0]==='S'?'E':'C').join('');

function runOne(spec: any, key: string, cfg: string, seed: number, frames: number, mcSteps: number) {
  const run = runMultiAgentFoldingToConvergence({
    sequenceCodes: spec.sequence3Letter,
    cdAlphaTargetPct: spec.cdAlphaTargetPct, cdBetaTargetPct: spec.cdBetaTargetPct,
    explicitDisulfidePairs: spec.explicitDisulfidePairs0Based ?? [],
    referencePdbText: spec.referencePdbText, referenceId: key,
    enableTargetGuidedRefinement: false, randomSeed: seed, initialNoiseDeg: 15,
    mjContactScale: cfg.includes('B') ? 1.0 : 0,
    hBondScale: cfg.includes('C') ? 1.5 : 0,
    secondaryStructurePredictor: cfg.includes('D') ? 'GOR' : 'LEGACY',
    monteCarloSteps: (cfg.includes('E')||cfg.includes('F')||cfg.includes('H')||cfg.includes('P')||cfg.includes('Z')) ? mcSteps : 0,
    zipperFraction: cfg.includes('Z') ? 0.05 : 0,
    betaScale: cfg.includes('P') ? 3 : 0,
    betaSS: cfg.includes('P') ? gorSS(spec.sequence3Letter) : undefined,
    burialScale: cfg.includes('H') ? 20 : 0,
    contactRestraints: cfg.includes('F') ? nativeContacts(spec.referencePdbText, spec.sequence3Letter.length) : undefined,
    contactWeight: 2.0,
  } as any, frames, 3);
  const g: any = run.groundTruthValidation;
  return { rmsd: g?.caRmsdAngstrom ?? null, tm: g?.tmScore ?? null, caPositions: run.caPositions };
}

// --- smoke test ---
const pga = getSpecLocal('tests/fixtures/1PGA.pdb', 'A');
console.log('1PGA dlugosc sekwencji:', pga.sequence3Letter.length);
console.log('1PGA GOR SS:', gorSS(pga.sequence3Letter));
const t0 = Date.now();
const r = runOne(pga, 'PROTEIN_G_56', 'A', 1, 60, 1000);
console.log('smoke A seed1: RMSD', r.rmsd?.toFixed(2), 'TM', r.tm?.toFixed(3), '| czas', ((Date.now()-t0)/1000).toFixed(1), 's');
