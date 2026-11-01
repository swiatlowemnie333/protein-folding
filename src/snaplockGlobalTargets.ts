// Global Protein Target Registry — Single Source of Truth across ALL SNAPLOCK Tabs (A..H)
// Enforces 1:1 sequence & reference correspondence (e.g. Ubiquitin = exact 76-aa 1UBQ/P0CG48[1-76],
// Insulin A = exact 21-aa mature Chain A 3I40:A, NEVER 110-aa preproinsulin P01308 or 685-aa polyubiquitin).

import { EXPERIMENTAL_GROUND_TRUTH_PDBS } from './snaplockExperimentalReferences';
import { buildBackboneFromTorsionsNeRF, ResidueTorsionState } from './snaplockPhysicalFoldingEngine';
import { buildAllAtomStructure, exportStandardAllAtomPdb } from './snaplockAllAtomBuilder';

export type GlobalProteinTargetId =
  | 'UBIQUITIN_76'
  | 'INSULIN_A_21'
  | 'CRAMBIN_46'
  | 'HELIX_BASELINE_18'
  | 'CONOTOXIN_25'
  | 'LYSOZYME_129'
  | 'CUSTOM';

export interface TargetMutationVariant {
  id: string;
  name: string;
  gene: string;
  clinicalImpact: string;
  healthyResidue: string;
  mutatedResidue: string;
  position: number; // 1-based residue index
  brokenBridges: string[];
  track5Status: 'LOCKED' | 'BROKEN' | 'NONE';
  coreStatus: 'COMPACT' | 'DESTABILIZED';
}

export interface GlobalProteinTargetSpec {
  id: Exclude<GlobalProteinTargetId, 'CUSTOM'>;
  shortLabel: string;
  titlePl: string;
  fastaHeader: string;
  sequence1Letter: string;
  sequence3Letter: string[];
  residueCount: number;
  uniprotLabel: string;
  pdbReferenceId: string;
  pdbChainId: string;
  referenceSource: 'RCSB_PDB' | 'CUSTOM_PDB';
  splitRole: 'TRAIN' | 'VALIDATION' | 'TEST' | 'HELIX_BASELINE';
  topologyClass: 'PURE_ALPHA_HELIX' | 'ALPHA_DISULFIDE' | 'ALPHA_BETA_SS' | 'ALPHA_BETA_GRASP' | 'ICK_KNOT';
  cdAlphaTargetPct: number;
  cdBetaTargetPct: number;
  explicitDisulfidePairs0Based: [number, number][]; // 0-based indices
  subdomainId: string;
  autocalibrationTargetId: string;
  referencePdbText: string;
  mutations: TargetMutationVariant[];
}

const ONE_TO_THREE: Record<string, string> = {
  A: 'ALA', R: 'ARG', N: 'ASN', D: 'ASP', C: 'CYS',
  Q: 'GLN', E: 'GLU', G: 'GLY', H: 'HIS', I: 'ILE',
  L: 'LEU', K: 'LYS', M: 'MET', F: 'PHE', P: 'PRO',
  S: 'SER', T: 'THR', W: 'TRP', Y: 'TYR', V: 'VAL',
  U: 'SEC'
};

function seq1ToSeq3(seq1: string): string[] {
  return seq1
    .trim()
    .toUpperCase()
    .split('')
    .map(ch => ONE_TO_THREE[ch] || 'ALA');
}

/**
 * Extracts strictly Chain A (21 residues: GLY A 1 .. ASN A 21 + intra-chain SSBOND CYS A 6 - CYS A 11)
 * from RCSB PDB 3I40 so mature 21-aa Insulin A is NEVER mixed with Chain B (30 aa) or P01308 preproinsulin (110 aa).
 */
function extractInsulinChainA21Pdb(): string {
  const raw3I40 = EXPERIMENTAL_GROUND_TRUTH_PDBS['3I40']?.pdbText || '';
  const lines = raw3I40.split(/\r?\n/);
  const out: string[] = [
    'HEADER    HORMONE (MATURE CHAIN A 21 AA)          01-JUL-09   3I40              ',
    'TITLE     HUMAN INSULIN MATURE A-CHAIN (RESIDUES 1-21, UNIPROT P01308[90-110])  ',
    'SSBOND   1 CYS A    6    CYS A   11                          1555   1555  2.03  '
  ];
  for (const line of lines) {
    if (line.startsWith('ATOM  ') || line.startsWith('TER   ')) {
      const chainChar = line.charAt(21);
      if (chainChar === 'A') {
        out.push(line);
      }
    }
  }
  out.push('END');
  return out.join('\n');
}

/**
 * Builds canonical Marqusee-Baldwin ideal α-helix reference PDB (18 aa: AEAAAKEAAAKEAAAKAA)
 * with Engh & Huber covalent bond geometry and canonical Pauling-Corey (φ=-57.8°, ψ=-47.0°, ω=180°)
 * stabilized by i -> i+4 GLU-LYS salt bridges.
 */
function buildCanonicalHelix18ReferencePdb(): string {
  const seq1 = 'AEAAAKEAAAKEAAAKAA';
  const seq3 = seq1ToSeq3(seq1);
  const torsions: ResidueTorsionState[] = seq3.map((code3, idx) => ({
    resIndex: idx,
    resSeq: idx + 1,
    resCode: code3,
    phiDeg: idx === 0 ? -62.0 : idx === seq3.length - 1 ? -61.0 : -57.8,
    psiDeg: idx === seq3.length - 1 ? -41.0 : -47.0,
    omegaDeg: 180.0,
    ssType: 'HELIX',
    segmentId: 0,
    isFlexibleHinge: false
  }));
  const bb = buildBackboneFromTorsionsNeRF(torsions, true);
  const caCoords: [number, number, number][] = bb.map(r => r.CA);
  const secTypes: ('HELIX' | 'SHEET' | 'TURN' | 'COIL')[] = torsions.map(() => 'HELIX');
  const allAtom = buildAllAtomStructure({
    sequenceCodes: seq3,
    caPositions: caCoords,
    ssBondPairs: [],
    secondaryTypes: secTypes,
    explicitBackboneByRes: bb
  });
  return exportStandardAllAtomPdb(
    'CANONICAL ALPHA HELIX BASELINE 18AA (MARQUSEE-BALDWIN)',
    seq3,
    allAtom,
    secTypes,
    { pipelineMode: 'CANONICAL_PAULING_COREY_ALPHA_HELIX' }
  );
}

const INSULIN_A_21_PDB = extractInsulinChainA21Pdb();
const HELIX_18_REF_PDB = buildCanonicalHelix18ReferencePdb();

export const GLOBAL_PROTEIN_TARGETS: Record<Exclude<GlobalProteinTargetId, 'CUSTOM'>, GlobalProteinTargetSpec> = {
  UBIQUITIN_76: {
    id: 'UBIQUITIN_76',
    shortLabel: 'Ubikwityna (76 aa)',
    titlePl: 'Ubikwityna dojrzała (UniProt P0CG48 [1–76] / PDB 1UBQ:A, 76 aa)',
    fastaHeader: '1UBQ_A|UniProt_P0CG48_1_76 Human Mature Ubiquitin [76 aa]',
    sequence1Letter: 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG',
    sequence3Letter: seq1ToSeq3('MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'),
    residueCount: 76,
    uniprotLabel: 'P0CG48 [1–76]',
    pdbReferenceId: '1UBQ',
    pdbChainId: 'A',
    referenceSource: 'RCSB_PDB',
    splitRole: 'VALIDATION',
    topologyClass: 'ALPHA_BETA_GRASP',
    cdAlphaTargetPct: 24,
    cdBetaTargetPct: 38,
    explicitDisulfidePairs0Based: [],
    subdomainId: 'ubiquitin_beta_hairpin',
    autocalibrationTargetId: 'ubiquitin_1ubq_A',
    referencePdbText: EXPERIMENTAL_GROUND_TRUTH_PDBS['1UBQ'].pdbText,
    mutations: [
      {
        id: 'UBQ_I44A',
        name: 'Mutacja ILE44 ➔ ALA (Łatka Hydrofobowa I44)',
        gene: 'UBB / P0CG48 (p.Ile44Ala)',
        clinicalImpact: 'Zniszczenie kluczowej łatki hydrofobowej Leu8-Ile44-Val70 rozpoznawanej przez proteasom 26S. Utrata degradacji białek ubikwitynowanych.',
        healthyResidue: 'ILE',
        mutatedResidue: 'ALA',
        position: 44,
        brokenBridges: ['Utrata kontaktu hydrofobowego ILE44 ↔ VAL70 / LEU8'],
        track5Status: 'NONE',
        coreStatus: 'DESTABILIZED'
      },
      {
        id: 'UBQ_L30E',
        name: 'Mutacja ILE30 ➔ GLU (Rdzeń Helisy α1)',
        gene: 'UBB / P0CG48 (p.Ile30Glu)',
        clinicalImpact: 'Wprowadzenie ładunku ujemnego (GLU, q=-1) do suchego rdzenia hydrofobowego między helisą α1 a 5-niciowym arkuszem β-grasp. Rozpad fałdu.',
        healthyResidue: 'ILE',
        mutatedResidue: 'GLU',
        position: 30,
        brokenBridges: ['Rozbicie rdzenia hydrofobowego α1–β1/β2 (ILE30 ↔ VAL17/ILE23)'],
        track5Status: 'NONE',
        coreStatus: 'DESTABILIZED'
      },
      {
        id: 'UBQ_K48R',
        name: 'Mutacja LYS48 ➔ ARG (Blokada Łańcuchów K48)',
        gene: 'UBB / P0CG48 (p.Lys48Arg)',
        clinicalImpact: 'Zablokowanie tworzenia poliubikwitynowych łańcuchów K48 kierujących substraty do proteasomu przy zachowaniu częściowej zwartości rdzenia.',
        healthyResidue: 'LYS',
        mutatedResidue: 'ARG',
        position: 48,
        brokenBridges: ['Zmiana objętości VdW na powierzchni pętli β3–β4 (135 ➔ 148 Å³)'],
        track5Status: 'NONE',
        coreStatus: 'COMPACT'
      }
    ]
  },

  INSULIN_A_21: {
    id: 'INSULIN_A_21',
    shortLabel: 'Insulina A (21 aa)',
    titlePl: 'Insulina ludzka — dojrzały łańcuch A (PDB 3I40:A [1–21], 21 aa, CYS6–CYS11)',
    fastaHeader: '3I40_A|Mature_Insulin_Chain_A [Homo sapiens, 21 aa, NOT 110-aa precursor]',
    sequence1Letter: 'GIVEQCCTSICSLYQLENYCN',
    sequence3Letter: seq1ToSeq3('GIVEQCCTSICSLYQLENYCN'),
    residueCount: 21,
    uniprotLabel: 'P01308 [90–110]',
    pdbReferenceId: '3I40_A_21',
    pdbChainId: 'A',
    referenceSource: 'RCSB_PDB',
    splitRole: 'TRAIN',
    topologyClass: 'ALPHA_DISULFIDE',
    cdAlphaTargetPct: 52,
    cdBetaTargetPct: 0,
    explicitDisulfidePairs0Based: [[5, 10]],
    subdomainId: 'insulin_a_core',
    autocalibrationTargetId: 'insulin_3i40_A',
    referencePdbText: INSULIN_A_21_PDB,
    mutations: [
      {
        id: 'INS_C6A',
        name: 'Mutacja CYS6 ➔ ALA (Klamra Wewnątrzłańcuchowa A6–A11)',
        gene: 'INS (p.Cys6Ala / Akita-like)',
        clinicalImpact: 'Zerwanie wewnątrzłańcuchowego mostka dwusiarczkowego CYS6–CYS11 spajającego helisę N-końcową A2–A8 z pętlą rdzenia. Cukrzyca neonatalna MODY.',
        healthyResidue: 'CYS',
        mutatedResidue: 'ALA',
        position: 6,
        brokenBridges: ['Mostek S-S CYS6 ↔ CYS11 (Zerwany, d(Sγ,Sγ) > 2.05 Å)'],
        track5Status: 'BROKEN',
        coreStatus: 'DESTABILIZED'
      },
      {
        id: 'INS_I2E',
        name: 'Mutacja ILE2 ➔ GLU (Katastrofa Rdzenia N-Helisy)',
        gene: 'INS (p.Ile2Glu)',
        clinicalImpact: 'Zastąpienie silnie hydrofobowej izoleucyny (Kyte +4.5) naładowanym glutaminianem (-3.5) niszczy amfipatyczny styk helis A1–A8 i A13–A19.',
        healthyResidue: 'ILE',
        mutatedResidue: 'GLU',
        position: 2,
        brokenBridges: ['Mostek CYS6–CYS11 zachowany, lecz rdzeń VdW rozepchnięty elektrostatycznie'],
        track5Status: 'LOCKED',
        coreStatus: 'DESTABILIZED'
      },
      {
        id: 'INS_L16P',
        name: 'Mutacja LEU16 ➔ PRO (Złamanie Helisy C-Końcowej)',
        gene: 'INS (p.Leu16Pro)',
        clinicalImpact: 'Wstawienie cyklicznej proliny (zablokowany kąt φ ≈ -63°) w środek helisy A13–A19 rozrywa wiązanie wodorowe i → i+4.',
        healthyResidue: 'LEU',
        mutatedResidue: 'PRO',
        position: 16,
        brokenBridges: ['Zerwanie wiązania wodorowego helisy α2 (SER12 O ⋯ LEU16 N)'],
        track5Status: 'LOCKED',
        coreStatus: 'DESTABILIZED'
      }
    ]
  },

  CRAMBIN_46: {
    id: 'CRAMBIN_46',
    shortLabel: 'Krambina (46 aa)',
    titlePl: 'Krambina roślinna 3x S-S (PDB 1CRN:A, X-ray 0.945 Å, 46 aa)',
    fastaHeader: '1CRN_A|Crambin [Crambe hispanica, 46 aa, 3x S-S]',
    sequence1Letter: 'TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN',
    sequence3Letter: seq1ToSeq3('TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN'),
    residueCount: 46,
    uniprotLabel: 'P01542 [1–46]',
    pdbReferenceId: '1CRN',
    pdbChainId: 'A',
    referenceSource: 'RCSB_PDB',
    splitRole: 'TRAIN',
    topologyClass: 'ALPHA_BETA_SS',
    cdAlphaTargetPct: 41,
    cdBetaTargetPct: 18,
    explicitDisulfidePairs0Based: [
      [2, 39],
      [3, 31],
      [15, 25]
    ],
    subdomainId: 'crambin_helix1',
    autocalibrationTargetId: 'crambin_1crn_A',
    referencePdbText: EXPERIMENTAL_GROUND_TRUTH_PDBS['1CRN'].pdbText,
    mutations: [
      {
        id: 'CRN_C16A',
        name: 'Mutacja CYS16 ➔ ALA (Mostek Helisa α1 ↔ α2)',
        gene: 'CRN (p.Cys16Ala)',
        clinicalImpact: 'Zerwanie wewnętrznej klamry CYS16–CYS26 spinającej dwie antyrównoległe helisy amfipatyczne krambiny.',
        healthyResidue: 'CYS',
        mutatedResidue: 'ALA',
        position: 16,
        brokenBridges: ['Mostek S-S CYS16 ↔ CYS26 (Zerwany)'],
        track5Status: 'BROKEN',
        coreStatus: 'DESTABILIZED'
      },
      {
        id: 'CRN_C3A',
        name: 'Mutacja CYS3 ➔ ALA (Makro-Klamra N–C CYS3–CYS40)',
        gene: 'CRN (p.Cys3Ala)',
        clinicalImpact: 'Otwarcie dalekozasięgowej klamry CYS3–CYS40 domykającej koniec N z pętlą C-końcową.',
        healthyResidue: 'CYS',
        mutatedResidue: 'ALA',
        position: 3,
        brokenBridges: ['Mostek S-S CYS3 ↔ CYS40 (Zerwany)'],
        track5Status: 'BROKEN',
        coreStatus: 'DESTABILIZED'
      }
    ]
  },

  HELIX_BASELINE_18: {
    id: 'HELIX_BASELINE_18',
    shortLabel: 'Bazowa α-Helisa (18 aa)',
    titlePl: 'Kanoniczna α-Helisa Bazowa Marqusee–Baldwina (18 aa, mostki solne E–K i→i+4)',
    fastaHeader: 'HELIX_18|Marqusee_Baldwin_Canonical_Alpha_Helix [18 aa]',
    sequence1Letter: 'AEAAAKEAAAKEAAAKAA',
    sequence3Letter: seq1ToSeq3('AEAAAKEAAAKEAAAKAA'),
    residueCount: 18,
    uniprotLabel: 'DESIGNED-HELIX-18',
    pdbReferenceId: 'HELIX18_REF',
    pdbChainId: 'A',
    referenceSource: 'CUSTOM_PDB',
    splitRole: 'HELIX_BASELINE',
    topologyClass: 'PURE_ALPHA_HELIX',
    cdAlphaTargetPct: 88,
    cdBetaTargetPct: 0,
    explicitDisulfidePairs0Based: [],
    subdomainId: 'crambin_helix1',
    autocalibrationTargetId: 'crambin_1crn_A',
    referencePdbText: HELIX_18_REF_PDB,
    mutations: [
      {
        id: 'HELIX_A9P',
        name: 'Mutacja ALA9 ➔ PRO (Helix-Breaker Kink w środku helisy)',
        gene: 'DESIGNED (p.Ala9Pro)',
        clinicalImpact: 'Klasyczny test zgięcia helisy α: prolina w pozycji 9 nie posiada protonu amidowego N-H i wymusza złamanie osi helisy o ~26°.',
        healthyResidue: 'ALA',
        mutatedResidue: 'PRO',
        position: 9,
        brokenBridges: ['Zerwanie wiązania wodorowego C=O(5) ⋯ H-N(9) w rdzeniu α-helisy'],
        track5Status: 'NONE',
        coreStatus: 'DESTABILIZED'
      },
      {
        id: 'HELIX_K6E',
        name: 'Mutacja LYS6 ➔ GLU (Odwrócenie Mostka Solnego i→i+4)',
        gene: 'DESIGNED (p.Lys6Glu)',
        clinicalImpact: 'Zamiana dodatniej lizyny (+1) na ujemny glutaminian (-1) zamienia stabilizujący mostek solny GLU2–LYS6 w odpychanie kulombowskie GLU2–GLU6.',
        healthyResidue: 'LYS',
        mutatedResidue: 'GLU',
        position: 6,
        brokenBridges: ['Zastąpienie mostka solnego GLU2 ↔ LYS6 odpychaniem elektrostatycznym (-1/-1)'],
        track5Status: 'NONE',
        coreStatus: 'DESTABILIZED'
      }
    ]
  },

  CONOTOXIN_25: {
    id: 'CONOTOXIN_25',
    shortLabel: 'ω-Konotoksyna (25 aa)',
    titlePl: 'ω-Konotoksyna MVIIA — węzeł cysteinowy ICK 3x S-S (PDB 1OMG:A, 25 aa)',
    fastaHeader: '1OMG_A|Omega_Conotoxin_MVIIA [Conus magus, 25 aa, 3x S-S ICK]',
    sequence1Letter: 'CKGKGAKCSRLMYDCCTGSCRSGKC',
    sequence3Letter: seq1ToSeq3('CKGKGAKCSRLMYDCCTGSCRSGKC'),
    residueCount: 25,
    uniprotLabel: 'P05484 [Mature 1–25]',
    pdbReferenceId: '1OMG',
    pdbChainId: 'A',
    referenceSource: 'RCSB_PDB',
    splitRole: 'VALIDATION',
    topologyClass: 'ICK_KNOT',
    cdAlphaTargetPct: 10,
    cdBetaTargetPct: 32,
    explicitDisulfidePairs0Based: [
      [0, 15],
      [7, 19],
      [14, 24]
    ],
    subdomainId: 'conotoxin_ick_knot',
    autocalibrationTargetId: 'ubiquitin_1ubq_A',
    referencePdbText: EXPERIMENTAL_GROUND_TRUTH_PDBS['1OMG'].pdbText,
    mutations: [
      {
        id: 'OMG_C1A',
        name: 'Mutacja CYS1 ➔ ALA (Rozwiązanie Węzła Cysteinowego ICK)',
        gene: 'MVIIA (p.Cys1Ala)',
        clinicalImpact: 'Zerwanie pierwszej klamry ICK (CYS1–CYS16) znosi blokadę kanału wapniowego CaV2.2.',
        healthyResidue: 'CYS',
        mutatedResidue: 'ALA',
        position: 1,
        brokenBridges: ['Mostek S-S CYS1 ↔ CYS16 (Zerwany)'],
        track5Status: 'BROKEN',
        coreStatus: 'DESTABILIZED'
      }
    ]
  },

  LYSOZYME_129: {
    id: 'LYSOZYME_129',
    shortLabel: 'Lizozym HEWL (129 aa)',
    titlePl: 'Lizozym HEWL 4x S-S (PDB 1LYZ:A, 129 aa — Held-Out TEST)',
    fastaHeader: '1LYZ_A|Hen_Egg_White_Lysozyme [Gallus gallus, 129 aa, 4x S-S]',
    sequence1Letter:
      'KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL',
    sequence3Letter: seq1ToSeq3(
      'KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL'
    ),
    residueCount: 129,
    uniprotLabel: 'P00698 [19–147]',
    pdbReferenceId: '1LYZ',
    pdbChainId: 'A',
    referenceSource: 'RCSB_PDB',
    splitRole: 'TEST',
    topologyClass: 'ALPHA_BETA_SS',
    cdAlphaTargetPct: 40,
    cdBetaTargetPct: 12,
    explicitDisulfidePairs0Based: [
      [5, 126],
      [29, 114],
      [63, 79],
      [75, 93]
    ],
    subdomainId: 'crambin_helix1',
    autocalibrationTargetId: 'lysozyme_1lyz_A',
    referencePdbText: EXPERIMENTAL_GROUND_TRUTH_PDBS['1LYZ'].pdbText,
    mutations: [
      {
        id: 'LYZ_C6A',
        name: 'Mutacja CYS6 ➔ ALA (Klamra N–C CYS6–CYS127)',
        gene: 'LYZ (p.Cys6Ala)',
        clinicalImpact: 'Utrata głównej makro-klamry spinającej koniec N z końcem C lizozymu (121 reszt dystansu sekwencyjnego); amyloidoza lizozymowa.',
        healthyResidue: 'CYS',
        mutatedResidue: 'ALA',
        position: 6,
        brokenBridges: ['Mostek S-S CYS6 ↔ CYS127 (Zerwany)'],
        track5Status: 'BROKEN',
        coreStatus: 'DESTABILIZED'
      }
    ]
  }
};

export const GLOBAL_TARGET_ORDER: Exclude<GlobalProteinTargetId, 'CUSTOM'>[] = [
  'UBIQUITIN_76',
  'INSULIN_A_21',
  'CRAMBIN_46',
  'HELIX_BASELINE_18',
  'CONOTOXIN_25',
  'LYSOZYME_129'
];

/**
 * Resolves any PDB / UniProt / Target alias strictly to the exact matching single-chain reference structure
 * when it corresponds to one of the canonical targets (preventing 110-aa P01308 or 685-aa P0CG48 precursor mismatches!).
 */
export function resolveStrictTargetReference(queryOrTargetId: string): {
  targetSpec: GlobalProteinTargetSpec;
  pdbText: string;
  resolvedId: string;
  chainId: string;
  source: 'RCSB_PDB' | 'CUSTOM_PDB';
  description: string;
} | null {
  const q = (queryOrTargetId || '').trim().toUpperCase();
  if (!q) return null;

  if (
    q === 'UBIQUITIN_76' ||
    q === '1UBQ' ||
    q === 'P0CG48' ||
    q === 'P0CG48_76' ||
    q === 'AF-P0CG48-F1' ||
    q === 'UBQ'
  ) {
    const t = GLOBAL_PROTEIN_TARGETS.UBIQUITIN_76;
    return {
      targetSpec: t,
      pdbText: t.referencePdbText,
      resolvedId: '1UBQ',
      chainId: 'A',
      source: 'RCSB_PDB',
      description: t.titlePl
    };
  }

  if (
    q === 'INSULIN_A_21' ||
    q === '3I40_A_21' ||
    q === '3I40:A' ||
    q === '3I40_A' ||
    q === 'INSULIN_A' ||
    q === 'P01308' ||
    q === 'AF-P01308-F1'
  ) {
    const t = GLOBAL_PROTEIN_TARGETS.INSULIN_A_21;
    return {
      targetSpec: t,
      pdbText: t.referencePdbText,
      resolvedId: '3I40_A_21',
      chainId: 'A',
      source: 'RCSB_PDB',
      description: t.titlePl
    };
  }

  if (q === 'CRAMBIN_46' || q === '1CRN' || q === 'P01542' || q === 'CRAMBIN') {
    const t = GLOBAL_PROTEIN_TARGETS.CRAMBIN_46;
    return {
      targetSpec: t,
      pdbText: t.referencePdbText,
      resolvedId: '1CRN',
      chainId: 'A',
      source: 'RCSB_PDB',
      description: t.titlePl
    };
  }

  if (q === 'HELIX_BASELINE_18' || q === 'HELIX18_REF' || q === 'HELIX_18') {
    const t = GLOBAL_PROTEIN_TARGETS.HELIX_BASELINE_18;
    return {
      targetSpec: t,
      pdbText: t.referencePdbText,
      resolvedId: 'HELIX18_REF',
      chainId: 'A',
      source: 'CUSTOM_PDB',
      description: t.titlePl
    };
  }

  if (q === 'CONOTOXIN_25' || q === '1OMG' || q === 'P05484' || q === 'AF-P05484-F1') {
    const t = GLOBAL_PROTEIN_TARGETS.CONOTOXIN_25;
    return {
      targetSpec: t,
      pdbText: t.referencePdbText,
      resolvedId: '1OMG',
      chainId: 'A',
      source: 'RCSB_PDB',
      description: t.titlePl
    };
  }

  if (q === 'LYSOZYME_129' || q === '1LYZ' || q === 'P00698' || q === 'AF-P00698-F1') {
    const t = GLOBAL_PROTEIN_TARGETS.LYSOZYME_129;
    return {
      targetSpec: t,
      pdbText: t.referencePdbText,
      resolvedId: '1LYZ',
      chainId: 'A',
      source: 'RCSB_PDB',
      description: t.titlePl
    };
  }

  return null;
}

/**
 * Detects if a given 3-letter sequence array matches one of the canonical global targets 1:1.
 */
export function matchSequenceToGlobalTarget(sequence3Letter: string[]): GlobalProteinTargetId {
  const joined = sequence3Letter.map(s => s.trim().toUpperCase()).filter(Boolean).join(' ');
  for (const tid of GLOBAL_TARGET_ORDER) {
    if (GLOBAL_PROTEIN_TARGETS[tid].sequence3Letter.join(' ') === joined) {
      return tid;
    }
  }
  return 'CUSTOM';
}
