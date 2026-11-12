import React, { useEffect, useRef, useState } from 'react';
import { runMultiAgentFoldingToConvergence, predictSecondaryStructureGOR, pickFragments, annealFragmentsMonteCarlo } from './src/snaplockMultiAgentEngine';
import { GLOBAL_PROTEIN_TARGETS } from './src/snaplockGlobalTargets';

// =====================================================================================
// KARTA J: Laboratorium wariantów A/B + Dziennik wniosków
// A = ścieżka bazowa. Każdy moduł (B, C, D, E) dokładany osobno i w zestawach.
// Każdy wariant liczony z kilku startów (średnia ± odchylenie), bo zwijanie jest chaotyczne.
// Wszystko na ślepo: silnik NIE widzi struktury docelowej (enableTargetGuidedRefinement: false).
// =====================================================================================

type ModuleKey = 'B' | 'C' | 'D' | 'E' | 'F' | 'H' | 'I' | 'P';
const MODULES: { key: ModuleKey; name: string; desc: string }[] = [
  { key: 'B', name: 'Kontakty MJ', desc: 'Energie kontaktu Miyazawa-Jernigan 1996 (210 par aminokwasów)' },
  { key: 'C', name: 'Wiązania H', desc: 'Wiązania wodorowe szkieletu O…N (2,9 Å)' },
  { key: 'D', name: 'GOR', desc: 'Struktura drugorzędowa wyuczona z 691 białek PDB (Q3 62,6%)' },
  { key: 'E', name: 'Monte Carlo', desc: 'Wyżarzanie Monte Carlo w przestrzeni kątów (φ, ψ)' },
  { key: 'F', name: 'Kontakty (test kontrolny)', desc: 'UWAGA: NIE na ślepo. Kontakty z prawdziwej struktury = górna granica tego, co da koewolucja/DCA. Włącza też Monte Carlo.' },
  { key: 'H', name: 'Zakopywanie (SASA)', desc: 'Efekt hydrofobowy: reszty o dużej tłustej powierzchni (SASA × hydropatia Kyte) nagradzane za bycie w środku białka. Na ślepo. Włącza też Monte Carlo.' },
  { key: 'P', name: 'Parowanie nici β', desc: 'Nagroda za „drabinki” nici β: pary reszt (i, j) blisko siebie razem z sąsiednimi parami, równolegle lub antyrównolegle. Waga wyższa tam, gdzie GOR przewiduje nić. Na ślepo. Włącza też Monte Carlo.' },
  { key: 'I', name: 'Hybryda: fragmenty PDB', desc: 'Zamiast kręcić jednym kątem, wstawia kawałki 9 i 3 aminokwasów z 255 prawdziwych białek (jak Rosetta), dobrane po sekwencji + GOR. Na ślepo. Biblioteka pobierana z RCSB przy pierwszym użyciu (~20 s).' }
];

// Kontakty CA–CA < 8 Å, |i-j| >= 6, z pliku PDB (tylko do testu kontrolnego F)
function nativeContacts(pdb: string, n: number): { i: number; j: number }[] {
  const ca: number[][] = []; const seen = new Set<number>();
  for (const l of (pdb || '').split('\n')) {
    if (!l.startsWith('ATOM') || l.slice(12, 16).trim() !== 'CA' || !(l[16] === ' ' || l[16] === 'A')) continue;
    const r = +l.slice(22, 26); if (seen.has(r)) continue; seen.add(r);
    ca.push([+l.slice(30, 38), +l.slice(38, 46), +l.slice(46, 54)]);
  }
  const m = Math.min(n, ca.length), out: { i: number; j: number }[] = [];
  for (let i = 0; i < m; i++) for (let j = i + 6; j < m; j++) if (Math.hypot(ca[i][0] - ca[j][0], ca[i][1] - ca[j][1], ca[i][2] - ca[j][2]) < 8) out.push({ i, j });
  return out;
}
const HISTORY_KEY = 'snaplock_historia_przebiegow_v1';

// ---- Moduł I: biblioteka fragmentów (kąty phi/psi z 255 łańcuchów PDB, pula treningowa GOR bez białek testowych) ----
const FRAG_IDS = '1O4W,2Z6F,1I2K,1WOL,1M1S,2RIL,2P8G,2QFL,3ES1,2HJN,1U53,1K8U,2RH3,1DIX,2HUH,2OI0,3CBN,12IS,2P14,1YN4,1VJK,3B84,1KP6,2CKK,1DVO,1BYI,1CXC,3EN8,2XPW,1BKB,2HDO,1UFY,2O71,3F6Y,2WND,1UV4,2OPC,3FOV,2YCI,1V77,2EJ9,1VI4,1UXO,2PND,2AJ6,3ESM,2NLR,2YF9,1JNI,1RC9,2AV4,1MG4,2DXA,2SGA,3EO5,2VDF,1RTQ,2BL1,1ZD0,3AV3,2UUR,2VCL,2HSH,3B7H,1EDT,2QSW,1WHO,2IUG,2NPN,1QQF,2GKP,2Y2Z,2HC1,1JBZ,2OSO,2VJW,1KNG,1QCS,2W6K,3CGX,2GYZ,2RBK,1D2T,2PTH,1ES9,2BZG,3BER,2I9I,1QOI,1QCZ,2V2P,1VBW,1BOL,1SFS,2XNQ,3FGH,1VK1,2IW0,2VSD,2ERF,2FMA,3E9V,1Y9L,3AS8,1A1X,2HNG,2VM9,1MOF,2WFO,2Z08,1TUK,2OD5,1G5T,2R16,2GFH,2YCD,3AGN,2V9V,2JEN,2R2Y,2G3A,3CC8,2EWH,1A58,2XDH,2XMF,2DP9,2PQ3,1DG6,1VIN,1R26,1R9W,1JLN,3G85,1DUW,1IFR,1J27,1CI3,2P5K,1B2L,1ZI8,1WLU,2ZZJ,1M1H,1RFE,1IGD,1QAD,1PZW,2Y1B,3ED5,3CHV,1YZF,2AMH,2P84,2CUL,2RKN,2QML,1LMI,1CTF,2C4W,2HSB,2EO4,2NYC,1RUT,1XPH,1KTH,1X3S,2VEB,3FTJ,2OA2,1BRT,1UXX,1BX7,2I9W,1OCY,1DY2,1Q42,3FLJ,2Z98,1YOY,3DHN,2XGR,2RFA,3FDJ,2C07,2IHD,2FJ9,2VXN,2VG9,1M1Q,1XTP,2CWY,1CFB,1SZ7,2DPM,2A4V,2ISB,2V3G,2CYJ,2CBP,1UCS,3BPU,2GKE,1OXJ,3G8Z,1KEA,3EBT,2ZK9,1OSH,3B7C,2O9U,1EP0,1OKS,2YBY,2QHP,2HAQ,1IWL,1U2H,1ZX6,1MVL,1HUW,2YZT,3BYR,2O9S,1BXE,2AR5,2JAY,2ZPM,1ZUH,2Y7B,1J6O,2YYO,2OLG,2OKM,3FGE,2FSR,2C8B,1WQ8,2P25,1S2W,3DEO,3CEC,1FXD,1KQ6,1ZCE,2R0J,3ERS,3E8T,1K4N,1WLJ,1FX2,3D6J,1S7I,1OOT,1WM3'.split(',');
const FRAG_KEY = 'snaplock_biblioteka_fragmentow_v1';
const A3: Record<string, string> = { ALA:'A',ARG:'R',ASN:'N',ASP:'D',CYS:'C',GLN:'Q',GLU:'E',GLY:'G',HIS:'H',ILE:'I',LEU:'L',LYS:'K',MET:'M',PHE:'F',PRO:'P',SER:'S',THR:'T',TRP:'W',TYR:'Y',VAL:'V',MSE:'M' };
async function loadFragmentLibrary(onProgress: (s: string) => void): Promise<{ aa: string; phi: number[]; psi: number[] }> {
  try { const c = localStorage.getItem(FRAG_KEY); if (c) return JSON.parse(c); } catch { /* ignore */ }
  const sub = (a: number[], b: number[]) => [a[0]-b[0], a[1]-b[1], a[2]-b[2]], cr = (a: number[], b: number[]) => [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]], dot = (a: number[], b: number[]) => a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
  const dih = (p0: number[], p1: number[], p2: number[], p3: number[]) => { const b0 = sub(p0, p1), b1 = sub(p2, p1), b2 = sub(p3, p2); const n = Math.hypot(b1[0], b1[1], b1[2]); const u = b1.map(x => x / n); const v = sub(b0, u.map(x => x * dot(b0, u))), w = sub(b2, u.map(x => x * dot(b2, u))); return Math.atan2(dot(cr(u, v), w), dot(v, w)) * 180 / Math.PI; };
  let aa = ''; const phi: number[] = [], psi: number[] = []; let k = 0, done = 0;
  const one = async (id: string) => {
    const t = await (await fetch(`https://files.rcsb.org/download/${id}.pdb`)).text(); const res: any[] = []; let chain: string | null = null;
    for (const l of t.split('\n')) { if (l.startsWith('ENDMDL')) break; if (!(l.startsWith('ATOM') || l.startsWith('HETATM'))) continue; const rn = l.slice(17, 20); if (!A3[rn]) continue; const ch = l[21]; if (chain === null) chain = ch; if (ch !== chain) break; const alt = l[16]; if (alt !== ' ' && alt !== 'A') continue; const key = l.slice(22, 27); let r = res[res.length - 1]; if (!r || r.key !== key) { r = { key, aa: A3[rn] }; res.push(r); } const nm = l.slice(12, 16).trim(); if (nm === 'N' || nm === 'CA' || nm === 'C') r[nm] = [+l.slice(30, 38), +l.slice(38, 46), +l.slice(46, 54)]; }
    const ok = res.filter(r => r.N && r.CA && r.C); let s = ''; const p: number[] = [], q: number[] = [];
    for (let i = 0; i < ok.length; i++) { const r = ok[i], a = ok[i - 1], b = ok[i + 1]; const brkP = !a || Math.hypot(...sub(r.N, a.C)) > 2, brkN = !b || Math.hypot(...sub(b.N, r.C)) > 2; s += r.aa; p.push(brkP ? 999 : Math.round(dih(a.C, r.N, r.CA, r.C))); q.push(brkN ? 999 : Math.round(dih(r.N, r.CA, r.C, b.N))); }
    aa += s + '/'; phi.push(...p, 999); psi.push(...q, 999);
  };
  const work = async () => { while (k < FRAG_IDS.length) { const id = FRAG_IDS[k++]; try { await one(id); } catch { /* pomiń */ } done++; onProgress(`Biblioteka fragmentów: ${done}/${FRAG_IDS.length} struktur z PDB`); } };
  await Promise.all(Array.from({ length: 8 }, work));
  const lib = { aa, phi, psi }; try { localStorage.setItem(FRAG_KEY, JSON.stringify(lib)); } catch { /* za duże – trzymamy w pamięci */ }
  return lib;
}
const ONE1: Record<string, string> = { ALA:'A',ARG:'R',ASN:'N',ASP:'D',CYS:'C',GLN:'Q',GLU:'E',GLY:'G',HIS:'H',ILE:'I',LEU:'L',LYS:'K',MET:'M',PHE:'F',PRO:'P',SER:'S',THR:'T',TRP:'W',TYR:'Y',VAL:'V' };

const PROTEINS: { key: string; label: string }[] = [
  { key: 'INSULIN_A_21', label: 'Insulina A (21)' },
  { key: 'CRAMBIN_46', label: 'Krambina (46)' },
  { key: 'UBIQUITIN_76', label: 'Ubikwityna (76)' },
  { key: 'LYSOZYME_129', label: 'Lizozym (129)' }
];

interface CellStats { rmsdMean: number; rmsdSd: number; tmMean: number; tmSd: number; best: number; n: number; }
interface GridRow { cfg: string; cells: Record<string, CellStats>; }
interface JournalEntry { date: string; text: string; }

const JOURNAL_KEY = 'snaplock_dziennik_wnioskow_v1';
const DEFAULT_JOURNAL: JournalEntry[] = [
  { date: '2026-11-12', text: 'WYNIK modułu P (waga 3, 5 startów, MC 5000): HIPOTEZA NIEPOTWIERDZONA dla ubikwityny – E 11,37±1,01 vs P 11,80±0,87 Å, EC 11,36±1,62 vs PC 11,68±1,33 Å (szum, bez poprawy). Krambina: PC 7,37±0,57 vs EC 8,04±0,43 Å (ok. 0,7 Å lepiej, t≈2,1 – sugestia, nie dowód). Moduł P zostaje jako opcja, NIE domyślny. Następny krok: sprawdzić, czy człon β w ogóle obniża energię natywnej ubikwityny względem wabików (problem energii czy problem przeszukiwania?).' },
  { date: '2026-11-12', text: 'HIPOTEZA (przed testem): ubikwityna ma fałszywe minimum, bo energia nie nagradza parowania nici β (38% β). Dodany moduł P – „drabinki” CA–CA < 5,5 Å z sąsiednią parą, waga wg GOR. Test A/B: E vs EP na ubikwitynie + kontrola na krambinie (α/β) i insulinie A (α), ≥ 5 startów. Wdrażamy jako domyślny tylko, jeśli poprawa > 2× szum i bez pogorszenia kontroli.' },
  { date: '2026-11-10', text: 'ESMFold (API ESM Atlas, darmowe, bez klucza) zmierzony na prawdziwych strukturach PDB: ubikwityna 0,83 Å (pLDDT 91), krambina 7,86 Å (pLDDT 50). Nasz najlepszy na krambinie: FRAG+K r1 7,53 ± 1,15 Å (12 startów) – na równi z ESMFold. Ubikwitynę ESMFold prawie na pewno widział w treningu, więc to nie jest uczciwe porównanie de novo.' },
  { date: '2026-11-10', text: 'Replikacja 12 startów: krambina MC 8,71 ± 0,51 → FRAG+K r1 7,53 ± 1,15 Å (p≈0,006), r2phys 7,85 ± 0,93 (p≈0,01). Ubikwityna bez poprawy (MC 11,31 → K 12,3): fałszywe minimum energii (natywna E −107,0 vs wabiki −110,9). Hipoteza: brakuje członu parowania nici β (ubikwityna 38% β).' },
  { date: '2026-11-10', text: 'Energia wyuczona r2phys z ograniczeniami znaków (VdW ≥ 0,3; mostki, wiązania H, pary PDB ≥ 0): natywna w top-3 w 42/45 białek CV, na 1. miejscu na obu białkach testowych. Runda 3 bez ograniczeń się wysypała (ujemna waga mostków) – odrzucona.' },
  { date: '2026-11-08', text: 'Tabela 24 właściwości: większość wartości zgodna z literaturą. Poprawione: LYS objętość 135 Å³ (Richards), SASA wg Millera 1987, HIS ładunek 0,04, chiralność GLY. Wpływ na RMSD: w granicach szumu.' },
  { date: '2026-11-08', text: 'BŁĄD KRYTYCZNY naprawiony: budowa łańcucha (NeRF) liczyła ψ względem C(i-2) zamiast N(i-1). Odbudowa prawdziwej krambiny z jej kątów: 7,84 Å -> 1,33 Å.' },
  { date: '2026-11-08', text: 'Test energii: prawdziwa struktura ma NIŻSZĄ energię niż wynik silnika (krambina, ubikwityna, lizozym). Wzory energii idą w dobrą stronę; wąskie gardło = przeszukiwanie.' },
  { date: '2026-11-09', text: 'GOR (struktura drugorzędowa z 691 białek PDB): Q3 62,6% vs 53,5% stara reguła vs 46,1% Chou-Fasman. Sufit metod z jednej sekwencji ~60-65%. Do 3D samo GOR nie wystarcza.' },
  { date: '2026-11-09', text: 'Siatka A/AB/AC/.../ABCD, 5 startów: różnice między wariantami mniejsze niż rozrzut (1-4 Å). Najlepszy trend: lizozym ABCD 14,5±1,3 Å vs A 15,8±1,4 Å; ubikwityna ACD 11,2±1,6 vs A 14,1±2,3. Żaden moduł nie jest jeszcze przełomem. Następny krok: koewolucja (kontakty z rodzin białek).' },
  { date: '2026-11-09', text: 'Moduł F, test kontrolny: nawet PRAWDZIWE kontakty w Monte Carlo dają mało (krambina 8,7 -> 7,2 Å; ubikwityna 10,6 -> 9,4 Å; realistyczne DCA 70%: 7,9 i 10,2 Å). Monte Carlo po jednym kącie nie umie spełnić kontaktów daleko w łańcuchu.' },
  { date: '2026-11-09', text: 'Moduł G (geometria odległości, poza aplikacją): pełna macierz odległości -> 0,3 Å (optymalizator działa). Same kontakty ~7 Å. Kontakty + pary DALEKO + prawdziwa struktura drugorzędowa -> 4,4 Å. WNIOSEK: kolejność = 1) lepsza struktura drugorzędowa (profil z rodziny białek), 2) przewidywanie odległości, 3) dopiero potem DCA.' },
  { date: '2026-11-09', text: 'Tabela 24 właściwości: w fizyce działa 13, martwych 6 (masa, SASA, polarność, rozpuszczalność, pI, rotamery). Masa nie wpływa na kształt równowagowy. Moduł H (SASA × hydropatia, zakopywanie), 3 starty, MC 10000: krambina 8,56±0,55 -> 8,33±0,52 Å, ubikwityna 11,28±1,45 -> 11,07±0,94 Å = w granicach szumu. Jak przy F: energia ok, przeszukiwanie nie nadąża.' },
  { date: '2026-11-09', text: 'Moduł I – HYBRYDA (fragmenty 9/3 z 255 białek PDB, jak Rosetta), 3 starty, 15000 kroków: krambina 8,17 vs zwykłe MC 8,13 Å; ubikwityna 10,93 vs 10,90 Å = bez zmian.' },
  { date: '2026-11-09', text: 'TEST KLUCZOWY (12 prób, wybór najniższej energii, krambina): najlepsza próba 6,5 Å, ale energia wybrała 8,3 Å. Korelacja energia–RMSD = −0,2 (zero). PRAWDZIWE WĄSKIE GARDŁO: funkcja energii nie odróżnia struktur 6 Å od 9 Å. Następny krok: energia statystyczna wyuczona z PDB (odległości par aminokwasów) i wybór przez klastrowanie prób.' },
  { date: '2026-11-09', text: 'Moduł J – energia statystyczna z PDB (3,08 mln par CB–CB z 255 łańcuchów, 210 par × 1 Å). Sama nauczyła się fizyki: Leu–Leu przyciąganie, Lys–Lys odpychanie, Lys–Glu mostek solny, Cys–Cys bardzo silnie. Prawdziwa krambina −9,5, ta sama struktura z potasowaną sekwencją +5,8; wszystkie próby silnika od −0,5 do +15 = potencjał wyraźnie wskazuje naturę.' },
  { date: '2026-11-09', text: 'Moduł J w przeszukiwaniu (waga 3–10): silnik znajduje struktury o energii NIŻSZEJ niż prawdziwa (−14,7 vs −9,5), ale błędne (krambina 8,8 Å, ubikwityna 12,1 Å). WNIOSEK: każdy składnik energii osobno da się oszukać (fałszywe minima). Następny krok: złożyć energię z kilku składników (pary J + wiązania H + zakopywanie H + zwartość) i dobrać wagi tak, by prawdziwa struktura była najniżej względem prób (uczenie wag na decoyach).' },
  { date: '2026-11-09', text: 'BŁĄD POMIARU naprawiony: MJ i wiązania H liczą się w składniku soap, a wcześniejsze testy decoyów miały je wyłączone. Poprawnie policzona energia silnika koreluje z RMSD: krambina +0,59, ubikwityna +0,35 (nie −0,2).' },
  { date: '2026-11-09', text: 'UCZENIE WAG (9 składników energii, 47 białek PDB, ~20 struktur każde: natywna + bliskie + hybryda; ranking parami). Walidacja krzyżowa na 45 białkach: prawdziwa struktura w top-3 w 24/45 przypadków (silnik: 9/45) = 2,7× lepiej. Test ubikwityna: korelacja +0,35 -> +0,58; krambina: natywna z 4. na 1. miejsce. ALE: wśród samych złych prób (>8 Å) żadna energia nie wybiera lepiej niż los. Wniosek: energia już rozpoznaje naturę; teraz przeszukiwanie musi dojść w jej pobliże (lejek).' }
];

// ---- Poprzeczka ESMFold: CA-RMSD po optymalnym nałożeniu (kwaterniony Horna) + średnie pLDDT ----
function caCoords(pdb: string): { xyz: number[][]; b: number[] } {
  const xyz: number[][] = []; const b: number[] = []; const seen = new Set<string>(); let chain = '';
  for (const l of (pdb || '').split('\n')) {
    if (!l.startsWith('ATOM') || l.slice(12, 16).trim() !== 'CA' || !(l[16] === ' ' || l[16] === 'A')) continue;
    if (!chain) chain = l[21]; if (l[21] !== chain) continue;
    const id = l.slice(22, 27); if (seen.has(id)) continue; seen.add(id);
    xyz.push([+l.slice(30, 38), +l.slice(38, 46), +l.slice(46, 54)]); b.push(+l.slice(60, 66));
  }
  return { xyz, b };
}
function optimalRmsd(P: number[][], Q: number[][]): number {
  const n = Math.min(P.length, Q.length); if (n < 3) return NaN;
  const c = (A: number[][]) => [0, 1, 2].map(k => A.slice(0, n).reduce((s, p) => s + p[k], 0) / n);
  const cp = c(P), cq = c(Q);
  const X = P.slice(0, n).map(p => p.map((x, k) => x - cp[k])), Y = Q.slice(0, n).map(p => p.map((x, k) => x - cq[k]));
  const S = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let i = 0; i < n; i++) for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) S[a][b] += X[i][a] * Y[i][b];
  const [[xx, xy, xz], [yx, yy, yz], [zx, zy, zz]] = S;
  const N = [[xx + yy + zz, yz - zy, zx - xz, xy - yx], [yz - zy, xx - yy - zz, xy + yx, zx + xz], [zx - xz, xy + yx, -xx + yy - zz, yz + zy], [xy - yx, zx + xz, yz + zy, -xx - yy + zz]];
  const e0 = X.reduce((s, p) => s + p[0] ** 2 + p[1] ** 2 + p[2] ** 2, 0) + Y.reduce((s, p) => s + p[0] ** 2 + p[1] ** 2 + p[2] ** 2, 0);
  let q = [0.5, 0.5, 0.5, 0.5];
  for (let it = 0; it < 20000; it++) { const w = N.map((r, i) => r.reduce((s, x, j) => s + x * q[j], 0) + e0 * q[i]); const no = Math.hypot(...w); q = w.map(x => x / no); }
  const lam = N.map(r => r.reduce((s, x, j) => s + x * q[j], 0)).reduce((s, x, j) => s + x * q[j], 0);
  return Math.sqrt(Math.max(0, (e0 - 2 * lam) / n));
}
async function esmfoldPdb(seq1: string): Promise<string> {
  try {
    const r = await fetch('/api/esmfold', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ sequence: seq1 }) });
    const d = await r.json(); if (r.ok && d.pdb) return d.pdb;
  } catch { /* spróbuj bezpośrednio */ }
  const r2 = await fetch('https://api.esmatlas.com/foldSequence/v1/pdb/', { method: 'POST', body: seq1 });
  if (!r2.ok) throw new Error('ESMFold HTTP ' + r2.status);
  return r2.text();
}
interface EsmRow { key: string; n: number; esmRmsd: number; plddt: number; ours: number | null; oursCfg: string; }

// Struktura drugorzędowa z GOR jako napis H/E/C (dla modułów I i P)
const gorSS = (seq3: string[]): string => predictSecondaryStructureGOR(seq3).map((x: string) => x[0] === 'H' ? 'H' : x[0] === 'S' ? 'E' : 'C').join('');

function allConfigs(selected: ModuleKey[]): string[] {
  const out: string[] = [];
  const n = selected.length;
  for (let mask = 0; mask < 1 << n; mask++) {
    out.push('A' + selected.filter((_, i) => mask & (1 << i)).join(''));
  }
  return out;
}

const mean = (a: number[]) => a.reduce((x, y) => x + y, 0) / Math.max(1, a.length);
const sd = (a: number[]) => { const m = mean(a); return Math.sqrt(mean(a.map(x => (x - m) ** 2))); };

export default function AblationLabView() {
  const [selectedModules, setSelectedModules] = useState<ModuleKey[]>(['B', 'C', 'D']);
  const [selectedProteins, setSelectedProteins] = useState<string[]>(['CRAMBIN_46', 'UBIQUITIN_76', 'LYSOZYME_129']);
  const [seeds, setSeeds] = useState(5);
  const [mcSteps, setMcSteps] = useState(5000);
  const [frames, setFrames] = useState(150);
  const [rows, setRows] = useState<GridRow[]>([]);
  const [running, setRunning] = useState(false);
  const [progress, setProgress] = useState('');
  const stopRef = useRef(false);
  const [esmRows, setEsmRows] = useState<EsmRow[]>([]);
  const [esmBusy, setEsmBusy] = useState(false);
  const runEsmBaseline = async () => {
    setEsmBusy(true); setEsmRows([]);
    const out: EsmRow[] = [];
    for (const key of selectedProteins) {
      const spec: any = (GLOBAL_PROTEIN_TARGETS as any)[key]; if (!spec) continue;
      const seq1 = spec.sequence3Letter.map((c: string) => ONE1[c] ?? 'X').join('');
      setProgress('ESMFold · ' + key);
      try {
        const pdb = await esmfoldPdb(seq1);
        const P = caCoords(pdb), Q = caCoords(spec.referencePdbText);
        const pl = P.b.length ? P.b.reduce((a, x) => a + x, 0) / P.b.length : 0;
        let ours: number | null = null, oursCfg = '';
        for (const r of rows) { const c = r.cells[key]; if (c && (ours === null || c.rmsdMean < ours)) { ours = c.rmsdMean; oursCfg = r.cfg; } }
        out.push({ key, n: Math.min(P.xyz.length, Q.xyz.length), esmRmsd: optimalRmsd(P.xyz, Q.xyz), plddt: pl <= 1 ? pl * 100 : pl, ours, oursCfg });
        setEsmRows([...out]);
      } catch (e: any) { setProgress('ESMFold błąd: ' + (e?.message || e)); }
    }
    if (out.length) {
      const date = new Date().toISOString().slice(0, 10);
      setJournal(j => [{ date, text: '[auto] Poprzeczka ESMFold (AI Mety, nie na ślepo: model mógł widzieć te białka w treningu): ' + out.map(r => r.key + ' ' + r.esmRmsd.toFixed(2) + ' Å, pLDDT ' + r.plddt.toFixed(0) + (r.ours !== null ? ' vs nasz ' + r.oursCfg + ' ' + r.ours.toFixed(2) + ' Å' : '')).join('; ') }, ...j]);
    }
    setProgress(''); setEsmBusy(false);
  };

  const [journal, setJournal] = useState<JournalEntry[]>(() => {
    try { const raw = localStorage.getItem(JOURNAL_KEY); if (raw) return JSON.parse(raw); } catch { /* ignore */ }
    return DEFAULT_JOURNAL;
  });
  // dopisz nowe wnioski domyślne, których jeszcze nie ma w zapisanym dzienniku
  useEffect(() => { setJournal(j => { const missing = DEFAULT_JOURNAL.filter(d => !j.some(x => x.text === d.text)); return missing.length ? [...missing.reverse(), ...j] : j; }); }, []);
  const [newNote, setNewNote] = useState('');
  useEffect(() => { try { localStorage.setItem(JOURNAL_KEY, JSON.stringify(journal)); } catch { /* ignore */ } }, [journal]);

  const toggle = <T,>(list: T[], v: T) => (list.includes(v) ? list.filter(x => x !== v) : [...list, v]);

  const runGrid = async () => {
    stopRef.current = false;
    setRunning(true);
    setRows([]);
    const configs = allConfigs(selectedModules);
    const total = configs.length * selectedProteins.length * seeds;
    let done = 0;
    const result: GridRow[] = [];
    const fragLib = configs.some(c => c.includes('I')) ? await loadFragmentLibrary(setProgress) : null;
    const fragCache: Record<string, any> = {};
    for (const cfg of configs) {
      const row: GridRow = { cfg, cells: {} };
      for (const key of selectedProteins) {
        const spec: any = (GLOBAL_PROTEIN_TARGETS as any)[key];
        if (!spec) continue;
        const r: number[] = [], t: number[] = [];
        for (let seed = 1; seed <= seeds; seed++) {
          if (stopRef.current) break;
          setProgress(`${cfg} · ${key} · start ${seed}/${seeds}  (${done}/${total})`);
          await new Promise(res => setTimeout(res, 0)); // oddaj sterowanie przeglądarce
          const run = runMultiAgentFoldingToConvergence(
            {
              sequenceCodes: spec.sequence3Letter,
              cdAlphaTargetPct: spec.cdAlphaTargetPct,
              cdBetaTargetPct: spec.cdBetaTargetPct,
              explicitDisulfidePairs: spec.explicitDisulfidePairs0Based ?? [],
              referencePdbText: spec.referencePdbText,
              referenceId: key,
              enableTargetGuidedRefinement: false,
              randomSeed: seed,
              initialNoiseDeg: 15,
              mjContactScale: cfg.includes('B') ? 1.0 : 0,
              hBondScale: cfg.includes('C') ? 1.5 : 0,
              secondaryStructurePredictor: cfg.includes('D') ? 'GOR' : 'LEGACY',
              monteCarloSteps: (!cfg.includes('I') && (cfg.includes('E') || cfg.includes('F') || cfg.includes('H') || cfg.includes('P'))) ? mcSteps : 0,
              betaScale: cfg.includes('P') ? 3 : 0,
              betaSS: cfg.includes('P') ? gorSS(spec.sequence3Letter) : undefined,
              burialScale: cfg.includes('H') ? 20 : 0,
              contactRestraints: cfg.includes('F') ? nativeContacts(spec.referencePdbText, spec.sequence3Letter.length) : undefined,
              contactWeight: 2.0
            },
            frames,
            3
          );
          let g: any = run.groundTruthValidation;
          if (cfg.includes('I') && fragLib) {
            if (!fragCache[key]) {
              const seq1 = spec.sequence3Letter.map((c: string) => ONE1[c] ?? 'X').join('');
              const ss = predictSecondaryStructureGOR(spec.sequence3Letter).map((x: string) => x[0] === 'H' ? 'H' : x[0] === 'S' ? 'E' : 'C').join('');
              fragCache[key] = { f9: pickFragments(fragLib, seq1, ss, 9), f3: pickFragments(fragLib, seq1, ss, 3) };
            }
            const sim: any = (run as any).simulation;
            annealFragmentsMonteCarlo(sim, fragLib, fragCache[key].f9, fragCache[key].f3, { steps: mcSteps, seed,
              contacts: cfg.includes('F') ? nativeContacts(spec.referencePdbText, spec.sequence3Letter.length) : undefined, burialScale: cfg.includes('H') ? 20 : 0, betaScale: cfg.includes('P') ? 3 : 0, betaSS: cfg.includes('P') ? gorSS(spec.sequence3Letter) : undefined });
            g = sim.exportFinalAllAtomAndValidate().groundTruthValidation;
          }
          if (g && g.caRmsdAngstrom != null) { r.push(g.caRmsdAngstrom); t.push(g.tmScore ?? 0); }
          done++;
        }
        row.cells[key] = { rmsdMean: mean(r), rmsdSd: sd(r), tmMean: mean(t), tmSd: sd(t), best: r.length ? Math.min(...r) : NaN, n: r.length };
      }
      result.push(row);
      setRows([...result]);
      if (stopRef.current) break;
    }
    setProgress(stopRef.current ? 'Zatrzymano.' : `Gotowe: ${result.length} wariantów × ${selectedProteins.length} białek × ${seeds} startów.`);
    // UCZENIE: każdy pełny przebieg zostaje w historii i dopisuje automatyczny wniosek do dziennika
    if (!stopRef.current && result.length) {
      const date = new Date().toISOString().slice(0, 10);
      try { const h = JSON.parse(localStorage.getItem(HISTORY_KEY) || '[]'); h.push({ date, seeds, mcSteps, frames, rows: result }); localStorage.setItem(HISTORY_KEY, JSON.stringify(h.slice(-200))); } catch { /* ignore */ }
      const a = result.find(r => r.cfg === 'A');
      const parts = selectedProteins.map(key => {
        const best = [...result].filter(r => r.cells[key]?.n).sort((x, y) => x.cells[key].rmsdMean - y.cells[key].rmsdMean)[0];
        if (!best) return '';
        const b = best.cells[key], base = a?.cells[key];
        let tag = '';
        if (base && best.cfg !== 'A') { const noise = Math.sqrt((base.rmsdSd ** 2 + b.rmsdSd ** 2) / Math.max(1, Math.min(base.n, b.n))); const diff = base.rmsdMean - b.rmsdMean; tag = diff > 2 * noise ? ' (LEPIEJ niż A)' : ' (w granicach szumu)'; }
        return `${key}: najlepszy ${best.cfg} ${b.rmsdMean.toFixed(2)}±${b.rmsdSd.toFixed(2)} Å` + (base ? ` vs A ${base.rmsdMean.toFixed(2)}` : '') + tag;
      }).filter(Boolean);
      setJournal(j => [{ date, text: `[auto] Siatka ${result.map(r => r.cfg).join('/')}, ${seeds} startów, MC ${mcSteps}: ` + parts.join('; ') }, ...j]);
    }
    setRunning(false);
  };

  const base = rows.find(r => r.cfg === 'A');
  const verdict = (row: GridRow, key: string) => {
    const a = base?.cells[key], c = row.cells[key];
    if (!a || !c || row.cfg === 'A') return '';
    const diff = a.rmsdMean - c.rmsdMean; // dodatnie = lepiej niż A
    const noise = Math.sqrt((a.rmsdSd ** 2 + c.rmsdSd ** 2) / Math.max(1, Math.min(a.n, c.n)));
    if (diff > 2 * noise) return 'lepiej';
    if (diff < -2 * noise) return 'gorzej';
    return 'szum';
  };

  // Synergia: efekt(XY) - [efekt(X) + efekt(Y)], gdzie efekt = poprawa RMSD względem A (Å). Dodatnia = moduły się wzmacniają.
  const synergy = (row: GridRow, key: string): number | null => {
    const mods = row.cfg.slice(1).split('');
    if (mods.length < 2 || !base?.cells[key] || !row.cells[key]) return null;
    const a = base.cells[key].rmsdMean;
    let sumSingles = 0;
    for (const m of mods) { const single = rows.find(r => r.cfg === 'A' + m)?.cells[key]; if (!single) return null; sumSingles += a - single.rmsdMean; }
    return (a - row.cells[key].rmsdMean) - sumSingles;
  };

  const exportCsv = () => {
    const head = ['wariant', ...selectedProteins.flatMap(k => [`${k}_rmsd`, `${k}_sd`, `${k}_tm`])].join(',');
    const lines = rows.map(r => [r.cfg, ...selectedProteins.flatMap(k => { const c = r.cells[k]; return c ? [c.rmsdMean.toFixed(3), c.rmsdSd.toFixed(3), c.tmMean.toFixed(4)] : ['', '', '']; })].join(','));
    const blob = new Blob([[head, ...lines].join('\n')], { type: 'text/csv' });
    const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = 'siatka_wariantow.csv'; a.click();
  };

  return (
    <div className="flex-1 w-full h-full overflow-auto p-4 text-zinc-200 font-mono text-[12px] space-y-4">
      <div>
        <h2 className="text-lg font-bold text-white">J: Laboratorium wariantów A/B</h2>
        <p className="text-zinc-400">A = ścieżka bazowa. Każdy moduł dokładany osobno i w zestawach. Każdy wariant z kilku startów — „lepiej/gorzej” tylko gdy różnica &gt; 2× szum. „syn” = synergia: czy moduły razem dają więcej (+) czy mniej (−) niż suma osobno. Wszystko na ślepo.</p>
      </div>

      <div className="grid md:grid-cols-3 gap-3">
        <div className="border border-zinc-800 rounded p-3 space-y-1">
          <div className="font-bold text-white mb-1">Moduły</div>
          {MODULES.map(m => (
            <label key={m.key} className="flex gap-2 items-start cursor-pointer">
              <input type="checkbox" checked={selectedModules.includes(m.key)} onChange={() => setSelectedModules(toggle(selectedModules, m.key))} disabled={running} />
              <span><b>{m.key}</b>: {m.name} <span className="text-zinc-500">— {m.desc}</span></span>
            </label>
          ))}
        </div>
        <div className="border border-zinc-800 rounded p-3 space-y-1">
          <div className="font-bold text-white mb-1">Białka</div>
          {PROTEINS.map(p => (
            <label key={p.key} className="flex gap-2 items-center cursor-pointer">
              <input type="checkbox" checked={selectedProteins.includes(p.key)} onChange={() => setSelectedProteins(toggle(selectedProteins, p.key))} disabled={running} />
              <span>{p.label}</span>
            </label>
          ))}
        </div>
        <div className="border border-zinc-800 rounded p-3 space-y-2">
          <div className="font-bold text-white">Ustawienia</div>
          <label className="flex justify-between gap-2">Starty na wariant <input className="w-16 bg-zinc-900 px-1" type="number" min={1} max={20} value={seeds} onChange={e => setSeeds(+e.target.value)} disabled={running} /></label>
          <label className="flex justify-between gap-2">Klatki symulacji <input className="w-16 bg-zinc-900 px-1" type="number" min={30} max={2000} value={frames} onChange={e => setFrames(+e.target.value)} disabled={running} /></label>
          <label className="flex justify-between gap-2">Kroki Monte Carlo (E) <input className="w-20 bg-zinc-900 px-1" type="number" min={0} max={200000} step={1000} value={mcSteps} onChange={e => setMcSteps(+e.target.value)} disabled={running} /></label>
          <div className="flex gap-2 pt-1">
            {!running ? (
              <button className="px-3 py-1 rounded bg-emerald-700 hover:bg-emerald-600 text-white font-bold" onClick={runGrid}>Uruchom siatkę ({allConfigs(selectedModules).length} wariantów)</button>
            ) : (
              <button className="px-3 py-1 rounded bg-rose-700 hover:bg-rose-600 text-white font-bold" onClick={() => (stopRef.current = true)}>Stop</button>
            )}
            {rows.length > 0 && !running && <button className="px-3 py-1 rounded bg-zinc-800 hover:bg-zinc-700" onClick={exportCsv}>CSV</button>}
          </div>
          <div className="text-zinc-400 min-h-[1.2em]">{progress}</div>
        </div>
      </div>

      {rows.length > 0 && (
        <div className="overflow-auto border border-zinc-800 rounded">
          <table className="w-full text-left">
            <thead className="bg-zinc-900 text-zinc-300">
              <tr>
                <th className="px-2 py-1">Wariant</th>
                {selectedProteins.map(k => <th key={k} className="px-2 py-1">{PROTEINS.find(p => p.key === k)?.label ?? k}<br /><span className="text-zinc-500 font-normal">RMSD Å (śr ± odch.) · TM · najlepszy</span></th>)}
              </tr>
            </thead>
            <tbody>
              {rows.map(r => (
                <tr key={r.cfg} className="border-t border-zinc-800">
                  <td className="px-2 py-1 font-bold text-white">{r.cfg}</td>
                  {selectedProteins.map(k => {
                    const c = r.cells[k]; const v = verdict(r, k);
                    return (
                      <td key={k} className="px-2 py-1">
                        {c ? <>
                          {c.rmsdMean.toFixed(2)} ± {c.rmsdSd.toFixed(2)} · {c.tmMean.toFixed(3)} · {c.best.toFixed(2)}
                          {v && <span className={`ml-2 px-1 rounded ${v === 'lepiej' ? 'bg-emerald-800 text-emerald-200' : v === 'gorzej' ? 'bg-rose-800 text-rose-200' : 'bg-zinc-800 text-zinc-400'}`}>{v}</span>}
                          {(() => { const syn = synergy(r, k); return syn === null ? null : <span className={`ml-2 ${syn > 0 ? 'text-emerald-300' : 'text-rose-300'}`} title="Synergia: o ile Å kombinacja jest lepsza (+) lub gorsza (−) niż suma pojedynczych modułów">syn {syn > 0 ? '+' : ''}{syn.toFixed(2)}</span>; })()}
                        </> : '—'}
                      </td>
                    );
                  })}
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}

      <div className="border border-sky-900 rounded p-3 space-y-2">
        <div className="flex justify-between items-center gap-2">
          <div className="font-bold text-white">Poprzeczka: ESMFold (AI Mety) dla zaznaczonych białek</div>
          <button className="px-3 py-1 rounded bg-sky-700 hover:bg-sky-600 disabled:opacity-50 text-white" disabled={esmBusy || running} onClick={runEsmBaseline}>{esmBusy ? 'Liczę…' : 'Policz ESMFold'}</button>
        </div>
        <p className="text-zinc-400">To NIE jest wynik na ślepo: ESMFold uczył się na PDB i mógł widzieć te białka. Służy jako uczciwa poprzeczka. pLDDT = pewność modelu (0–100); poniżej 70 model zgaduje. „Nasz” = najlepszy wariant z siatki powyżej (uruchom siatkę, żeby porównać).</p>
        {esmRows.length > 0 && (
          <table className="w-full text-left">
            <thead><tr className="text-zinc-400"><th>Białko</th><th>reszt</th><th>ESMFold RMSD</th><th>pLDDT</th><th>Nasz najlepszy</th></tr></thead>
            <tbody>{esmRows.map(r => (
              <tr key={r.key} className="border-t border-zinc-800">
                <td>{r.key}</td><td>{r.n}</td><td>{r.esmRmsd.toFixed(2)} Å</td>
                <td className={r.plddt >= 70 ? 'text-emerald-300' : 'text-amber-300'}>{r.plddt.toFixed(0)}</td>
                <td>{r.ours !== null ? <span className={r.ours <= r.esmRmsd ? 'text-emerald-300 font-bold' : ''}>{r.ours.toFixed(2)} Å ({r.oursCfg}){r.ours <= r.esmRmsd ? ' ✓ lepiej lub równo' : ''}</span> : '—'}</td>
              </tr>))}</tbody>
          </table>
        )}
      </div>

      <div className="border border-zinc-800 rounded p-3 space-y-2">
        <div className="font-bold text-white">Dziennik wniosków (zapisujemy, nie wdrażamy bez dowodu)</div>
        <div className="flex gap-2">
          <input className="flex-1 bg-zinc-900 px-2 py-1 rounded" placeholder="Nowy wniosek z badania…" value={newNote} onChange={e => setNewNote(e.target.value)} />
          <button className="px-3 py-1 rounded bg-indigo-700 hover:bg-indigo-600 text-white" onClick={() => { if (newNote.trim()) { setJournal([{ date: new Date().toISOString().slice(0, 10), text: newNote.trim() }, ...journal]); setNewNote(''); } }}>Dodaj</button>
        </div>
        <ul className="space-y-1">
          {journal.map((j, i) => (
            <li key={i} className="border-t border-zinc-800 pt-1"><span className="text-zinc-500">{j.date}</span> — {j.text}</li>
          ))}
        </ul>
      </div>
    </div>
  );
}
