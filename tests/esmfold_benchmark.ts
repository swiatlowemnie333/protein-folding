// ESMFold vs SNAPLOCK - uczciwy benchmark na tych samych celach
import fs from 'node:fs';
import path from 'node:path';
import { validateAgainstGroundTruth } from '../src/ground_truth_validation';
import { getCanonicalBenchmarkTargets } from '../src/snaplockAutocalibrationLoop';

const URL_ESM = 'https://api.esmatlas.com/foldSequence/v1/pdb/';
const OUT = path.join(process.cwd(), 'results', 'esmfold');
const T: Record<string, string> = { ALA:'A',ARG:'R',ASN:'N',ASP:'D',CYS:'C',GLN:'Q',GLU:'E',GLY:'G',HIS:'H',ILE:'I',
  LEU:'L',LYS:'K',MET:'M',PHE:'F',PRO:'P',SER:'S',THR:'T',TRP:'W',TYR:'Y',VAL:'V' };
const one = (c: string[]) => c.map((x) => (x.length === 1 ? x.toUpperCase() : T[x.toUpperCase()] || 'X')).join('');

function plddt(pdb: string): number | null {
  const v: number[] = [];
  for (const l of pdb.split('\n')) if (l.startsWith('ATOM') && l.substring(12, 16).trim() === 'CA') {
    const b = parseFloat(l.substring(60, 66)); if (!Number.isNaN(b)) v.push(b <= 1 ? b * 100 : b);
  }
  return v.length ? v.reduce((a, b) => a + b, 0) / v.length : null;
}

async function fold(seq: string): Promise<string> {
  for (let i = 1; i <= 3; i++) {
    try {
      const r = await fetch(URL_ESM, { method: 'POST', headers: { 'Content-Type': 'text/plain' }, body: seq });
      const t = await r.text();
      if (r.ok && t.includes('ATOM')) return t;
      console.log(`  HTTP ${r.status}, proba ${i}/3`);
    } catch (e) { console.log(`  blad sieci, proba ${i}/3: ${e}`); }
    await new Promise((r) => setTimeout(r, 5000 * i));
  }
  throw new Error('ESMFold niedostepny');
}

async function main() {
  fs.mkdirSync(OUT, { recursive: true });
  const rows: any[] = [];
  for (const t of getCanonicalBenchmarkTargets()) {
    const seq = one(t.sequenceCodes);
    console.log(`>> ${t.pdbId}:${t.chainId} [${t.split}] ${seq.length} aa`);
    if (seq.length < 10 || seq.length > 400 || seq.includes('X')) continue;
    try {
      const pdb = await fold(seq);
      fs.writeFileSync(path.join(OUT, `${t.targetId}.pdb`), pdb);
      const g = validateAgainstGroundTruth(pdb, t.referencePdbText, { modelChainId: 'A', referenceChainId: t.chainId,
        referenceId: t.pdbId, referenceSource: t.referenceSource, allowSubsequenceAlignment: false });
      rows.push({ cel: `${t.pdbId}:${t.chainId} ${t.split} ${seq.length}aa`, rmsd: g.caRmsdAngstrom?.toFixed(2),
        tm: g.tmScore?.toFixed(3), plddt: plddt(pdb)?.toFixed(1) });
    } catch (e) { rows.push({ cel: `${t.pdbId}:${t.chainId}`, rmsd: String(e) }); }
    await new Promise((r) => setTimeout(r, 2000));
  }
  fs.writeFileSync(path.join(OUT, 'esmfold_benchmark.json'), JSON.stringify(rows, null, 2));
  console.table(rows);
}
main();
