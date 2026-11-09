// Independent Ground-Truth Structural Validation Module for SNAPLOCK
// Completely decoupled from internal E8/L2 heuristic scores and internal quantization metrics.
// Implements:
// 1. Strict wwPDB coordinate & sequence parser (explicit chain selection, multi-MODEL, gaps, UNK, NaN checks)
// 2. Sequence-based residue mapping (Needleman-Wunsch with strict amino-acid identity verification)
// 3. Optimal rigid-body superposition via Kabsch algorithm with explicit 3x3 SVD (det(R) = +1)
// 4. Exact XYZ Cα RMSD: sqrt(mean(sum((X_ref - X_aligned)^2, axis=1))) + Backbone & All-Atom RMSD
// 5. Independent secondary RMSD re-computation from raw coordinates + Davenport quaternion cross-check
// 6. Physical protein geometry verification (peptide Cα-Cα bounds, steric clash exclusion)
// 7. Iteration & artifact hash binding + configurable pass criteria

export type Vec3Tuple = [number, number, number];
export type Mat3Tuple = [
  [number, number, number],
  [number, number, number],
  [number, number, number]
];

export const STANDARD_21_AMINO_ACIDS = new Set<string>([
  'ALA', 'ARG', 'ASN', 'ASP', 'CYS',
  'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
  'LEU', 'LYS', 'MET', 'PHE', 'PRO',
  'SER', 'THR', 'TRP', 'TYR', 'VAL',
  'SEC', 'PYL'
]);

export const MODIFIED_RESIDUE_CANONICAL_MAP: Record<string, string> = {
  MSE: 'MET',
  CSE: 'SEC',
  HSD: 'HIS',
  HSE: 'HIS',
  HSP: 'HIS',
  HID: 'HIS',
  HIE: 'HIS',
  HIP: 'HIS',
  CYX: 'CYS',
  CYM: 'CYS'
};

export const THREE_TO_ONE_MAP: Record<string, string> = {
  ALA: 'A', ARG: 'R', ASN: 'N', ASP: 'D', CYS: 'C',
  GLN: 'Q', GLU: 'E', GLY: 'G', HIS: 'H', ILE: 'I',
  LEU: 'L', LYS: 'K', MET: 'M', PHE: 'F', PRO: 'P',
  SER: 'S', THR: 'T', TRP: 'W', TYR: 'Y', VAL: 'V',
  SEC: 'U', PYL: 'O'
};

const SOLVENT_AND_ION_NAMES = new Set<string>([
  'HOH', 'WAT', 'DOD', 'SO4', 'PO4', 'GOL', 'EDO', 'PEG', 'ACT', 'DMS',
  'NA', 'CL', 'MG', 'ZN', 'CA', 'K', 'FE', 'MN', 'CU', 'CO', 'NI', 'IOD'
]);

export interface ParsedPdbResidue {
  indexInChain: number;
  chainId: string;
  resSeq: number;
  iCode: string;
  rawResName: string;
  canonicalResName: string; // Standard 3-letter AA or 'UNK'
  isKnownAminoAcid: boolean;
  hasCaAtom: boolean;
  caCoord: Vec3Tuple | null;
  atoms: Map<string, Vec3Tuple>;
  bFactorMean: number;
}

export interface SequenceNumberingGap {
  afterResSeq: number;
  beforeResSeq: number;
  missingCountEstimate: number;
}

export interface ParsedPdbStructure {
  headerId: string;
  title: string;
  totalModelsDetected: number;
  selectedModelNumber: number;
  availableChains: string[];
  selectedChainId: string;
  residues: ParsedPdbResidue[];
  residuesWithCa: ParsedPdbResidue[];
  residuesMissingCa: { resSeq: number; resName: string }[];
  unknownResidues: { resSeq: number; rawResName: string }[];
  modifiedResiduesNormalized: { resSeq: number; rawResName: string; mappedTo: string }[];
  sequenceNumberingGaps: SequenceNumberingGap[];
  invalidCoordinateLines: { lineNumber: number; reason: string; rawLine: string }[];
  altLocFilteredCount: number;
  duplicateAtomIgnoredCount: number;
  rawPdbHash: string;
}

export interface ResiduePairMapping {
  modelChainIndex: number;
  refChainIndex: number;
  modelResSeq: number;
  refResSeq: number;
  modelResName: string;
  refResName: string;
  aminoAcidIdentical: boolean;
}

export interface SequenceMappingSummary {
  pairedResidues: ResiduePairMapping[];
  mismatchedPairs: ResiduePairMapping[];
  unpairedModelResidues: { resSeq: number; resName: string }[];
  unpairedRefResidues: { resSeq: number; resName: string }[];
  pairedCount: number;
  exactMatchCount: number;
  modelTotalCaResidues: number;
  refTotalCaResidues: number;
  modelCoveragePct: number;
  referenceCoveragePct: number;
  sequenceIdentityPct: number;
  alignmentMode: 'GLOBAL_NEEDLEMAN_WUNSCH' | 'SEMIGLOBAL_SUBSEQUENCE';
}

export interface Svd3x3Result {
  U: Mat3Tuple;
  S: Vec3Tuple;
  V: Mat3Tuple;
  reconstructionMaxError: number;
  orthogonalityMaxError: number;
}

export interface KabschSuperpositionResult {
  rotationMatrix: Mat3Tuple;       // R such that x_aligned = R * x_model + translationVec
  translationVec: Vec3Tuple;       // t = mu_ref - R * mu_model
  modelCentroid: Vec3Tuple;
  referenceCentroid: Vec3Tuple;
  determinantR: number;
  improperReflectionCorrected: boolean;
  svdSingularValues: Vec3Tuple;
  svdReconstructionError: number;
  alignedModelCaCoords: Vec3Tuple[];
  referenceCaCoords: Vec3Tuple[];
  caRmsdAngstrom: number;
}

export interface PhysicalGeometryCheckResult {
  geometryValid: boolean;
  checkedResidueCount: number;
  caBondCount: number;
  meanCaCaDistAngstrom: number;
  minCaCaDistAngstrom: number;
  maxCaCaDistAngstrom: number;
  caBondViolations: { resSeqA: number; resSeqB: number; distanceAngstrom: number; issue: 'TOO_SHORT' | 'TOO_LONG' }[];
  nonBondedCaClashes: { resSeqA: number; resSeqB: number; distanceAngstrom: number }[];
  heavyAtomStericClashesCount: number;
  collapsedChainDetected: boolean;
  reasons: string[];
}

export interface GroundTruthResidueDeviation {
  modelResSeq: number;
  refResSeq: number;
  modelResName: string;
  refResName: string;
  code1: string;
  caDistAngstrom: number;
  backboneRmsdAngstrom: number;
  sidechainRmsdAngstrom: number;
  allAtomRmsdAngstrom: number;
  maxAtomDistAngstrom: number;
  maxAtomName: string;
  matchedAtomCount: number;
  matchedSidechainAtomCount: number;
  refBfactorOrPlddt: number;
}

export interface GroundTruthPassCriteria {
  testId?: string;
  testName?: string;
  targetCaRmsdAngstrom: number;          // Default: 1.0 Å
  intermediateCaRmsdAngstrom?: number;   // Optional intermediate/fallback threshold (e.g. 1.5 Å if explicitly set)
  minModelCoveragePct: number;           // Default: 95.0%
  minReferenceCoveragePct: number;       // Default: 95.0% for full-chain tests, 0% when validating mature chain vs precursor
  requireExactSequenceIdentity: boolean; // Default: true
  requireValidPhysicalGeometry: boolean; // Default: true
}

export const DEFAULT_HIGH_ACCURACY_CRITERIA: GroundTruthPassCriteria = {
  targetCaRmsdAngstrom: 1.0,
  minModelCoveragePct: 95.0,
  minReferenceCoveragePct: 95.0,
  requireExactSequenceIdentity: true,
  requireValidPhysicalGeometry: true
};

export type GroundTruthIterationStatus =
  | 'PASSED_TARGET'
  | 'PASSED_INTERMEDIATE'
  | 'BELOW_THRESHOLD'
  | 'GEOMETRY_VIOLATION'
  | 'VALIDATION_ERROR'
  | 'NO_REFERENCE';

export interface GroundTruthValidationOptions {
  modelChainId?: string;
  referenceChainId?: string;
  modelNumber?: number;
  referenceModelNumber?: number;
  referenceId?: string;
  referenceName?: string;
  referenceSource?: 'RCSB_PDB' | 'ALPHAFOLD_DB' | 'CUSTOM_PDB' | 'TEST_FIXTURE';
  iterationId?: string | number;
  parameterVersion?: string;
  modelArtifactPath?: string;
  allowSequenceMismatches?: boolean;
  allowSubsequenceAlignment?: boolean; // True when aligning e.g. 21-aa mature Insulin A to 110-aa P01308 precursor
  numericalToleranceAngstrom?: number; // Default 1e-7 Å
  criteria?: Partial<GroundTruthPassCriteria>;
}

export interface GroundTruthValidationReport {
  valid: boolean;
  hasGroundTruthReference: boolean;
  status: GroundTruthIterationStatus;
  passedTargetCriterion: boolean;
  passedIntermediateCriterion: boolean;
  failureReason: string | null;
  warnings: string[];

  // Provenance & iteration binding (prevents stale/cached result reuse)
  iterationId: string;
  parameterVersion: string;
  timestampIso: string;
  modelArtifactPath: string | null;
  modelPdbHash: string;
  referencePdbHash: string | null;

  // Reference & chain metadata
  referenceId: string | null;
  referenceName: string | null;
  referenceSource: 'RCSB_PDB' | 'ALPHAFOLD_DB' | 'CUSTOM_PDB' | 'TEST_FIXTURE' | 'NONE';
  selectedModelChain: string;
  availableModelChains: string[];
  selectedReferenceChain: string | null;
  availableReferenceChains: string[];
  selectedModelNumber: number;
  totalModelStructures: number;
  selectedReferenceModelNumber: number | null;
  totalReferenceStructures: number;

  // True Ground-Truth Metrics (NULL when no reference or validation error — NEVER 0.0!)
  caRmsdAngstrom: number | null;
  independentCaRmsdAngstrom: number | null;
  quaternionCrossCheckCaRmsdAngstrom: number | null;
  numericalAgreementDeltaAngstrom: number | null;
  independentVerificationPassed: boolean;

  backboneRmsdAngstrom: number | null;
  allHeavyAtomRmsdAngstrom: number | null;
  tmScore: number | null;
  gdtTsPct: number | null;
  meanReferencePlddtOrBfactor: number | null;

  // Residue pairing & coverage
  pairedResidueCount: number;
  modelResidueCount: number;
  referenceResidueCount: number;
  modelCoveragePct: number;
  referenceCoveragePct: number;
  sequenceIdentityPct: number;
  mismatchedResidueCount: number;

  // Superposition transform (model -> reference)
  superposition: {
    rotationMatrix: Mat3Tuple;
    translationVec: Vec3Tuple;
    modelCentroid: Vec3Tuple;
    referenceCentroid: Vec3Tuple;
    determinantR: number;
    svdSingularValues: Vec3Tuple;
  } | null;

  // Detailed diagnostics
  sequenceMapping: SequenceMappingSummary | null;
  modelGeometry: PhysicalGeometryCheckResult;
  perResidueDeviations: GroundTruthResidueDeviation[];
  alignedModelCaCoords: Vec3Tuple[];
  pairedReferenceCaCoords: Vec3Tuple[];
  criteriaApplied: GroundTruthPassCriteria;
}

// -----------------------------------------------------------------------------
// 1. DETERMINISTIC 64-BIT FNV-1a ARTIFACT HASH
// -----------------------------------------------------------------------------

export function computeDeterministicArtifactHash(content: string): string {
  let h1 = 0x811c9dc5;
  let h2 = 0x811c9dc5 ^ 0x9e3779b9;
  for (let i = 0; i < content.length; i++) {
    const ch = content.charCodeAt(i);
    h1 ^= ch;
    h1 = Math.imul(h1, 0x01000193);
    h2 ^= ch ^ (i & 0xff);
    h2 = Math.imul(h2, 0x01000193);
  }
  const hex1 = (h1 >>> 0).toString(16).padStart(8, '0');
  const hex2 = (h2 >>> 0).toString(16).padStart(8, '0');
  return `fnv64:${hex1}${hex2}`;
}

// -----------------------------------------------------------------------------
// 2. STRICT PDB PARSER WITH CHAIN, MULTI-MODEL, GAP & COORDINATE SANITY CHECKS
// -----------------------------------------------------------------------------

export function parsePdbForGroundTruth(
  pdbText: string,
  options?: { chainId?: string; modelNumber?: number }
): ParsedPdbStructure {
  const raw = typeof pdbText === 'string' ? pdbText : '';
  const rawPdbHash = computeDeterministicArtifactHash(raw);
  const lines = raw.split(/\r?\n/);

  let headerId = 'UNKNOWN';
  let titleParts: string[] = [];

  // First pass: discover models & chains, validate coordinate formatting
  let currentModel = 1;
  let explicitModelRecordsCount = 0;
  const targetModel = options?.modelNumber ?? 1;

  interface RawAtomEntry {
    modelNum: number;
    lineNumber: number;
    recordType: 'ATOM' | 'HETATM';
    atomName: string;
    altLoc: string;
    rawResName: string;
    chainId: string;
    resSeq: number;
    iCode: string;
    x: number;
    y: number;
    z: number;
    bFactor: number;
  }

  const invalidCoordinateLines: { lineNumber: number; reason: string; rawLine: string }[] = [];
  const atomsByModel = new Map<number, RawAtomEntry[]>();
  let altLocFilteredCount = 0;

  for (let idx = 0; idx < lines.length; idx++) {
    const line = lines[idx];
    const lineNo = idx + 1;
    if (!line) continue;

    if (line.startsWith('HEADER')) {
      const codeCandidate = line.substring(62, 66).trim();
      if (codeCandidate) headerId = codeCandidate;
      else {
        const rest = line.substring(6).trim();
        if (rest) headerId = rest.split(/\s+/)[0].substring(0, 24);
      }
      continue;
    }

    if (line.startsWith('TITLE ')) {
      titleParts.push(line.substring(10).trim());
      continue;
    }

    if (line.startsWith('MODEL ')) {
      explicitModelRecordsCount++;
      const parsedModelNum = parseInt(line.substring(6).trim(), 10);
      currentModel = Number.isFinite(parsedModelNum) && parsedModelNum > 0
        ? parsedModelNum
        : explicitModelRecordsCount;
      continue;
    }

    if (line.startsWith('ENDMDL')) {
      continue;
    }

    const isAtom = line.startsWith('ATOM  ');
    const isHetAtm = line.startsWith('HETATM');
    if (!isAtom && !isHetAtm) continue;

    const recordType: 'ATOM' | 'HETATM' = isAtom ? 'ATOM' : 'HETATM';
    const rawResName = line.substring(17, 20).trim().toUpperCase();

    // Ignore bulk water and inorganic ions
    if (SOLVENT_AND_ION_NAMES.has(rawResName)) continue;
    // For HETATM, only include modified polymer amino acids (MSE, SEC, CSE, etc.) or UNK
    if (
      isHetAtm &&
      !STANDARD_21_AMINO_ACIDS.has(rawResName) &&
      !MODIFIED_RESIDUE_CANONICAL_MAP[rawResName] &&
      rawResName !== 'UNK'
    ) {
      continue;
    }

    const altLoc = line.length > 16 ? line[16] : ' ';
    if (altLoc !== ' ' && altLoc !== 'A' && altLoc !== '1') {
      altLocFilteredCount++;
      continue;
    }

    const atomName = line.substring(12, 16).trim().toUpperCase();
    const chainChar = line.length > 21 ? line[21].trim() : '';
    const chainId = chainChar || 'A';
    const resSeqRaw = line.substring(22, 26).trim();
    const resSeq = parseInt(resSeqRaw, 10);
    const iCode = line.length > 26 ? line[26].trim() : '';

    if (!Number.isFinite(resSeq)) {
      invalidCoordinateLines.push({
        lineNumber: lineNo,
        reason: `Nieprawidłowy numer reszty (resSeq="${resSeqRaw}")`,
        rawLine: line.substring(0, 80)
      });
      continue;
    }

    const xStr = line.substring(30, 38).trim();
    const yStr = line.substring(38, 46).trim();
    const zStr = line.substring(46, 54).trim();

    const x = Number(xStr);
    const y = Number(yStr);
    const z = Number(zStr);

    if (
      !xStr || !yStr || !zStr ||
      !Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)
    ) {
      invalidCoordinateLines.push({
        lineNumber: lineNo,
        reason: `Nieprawidłowe lub nieliczbowe współrzędne XYZ ("${xStr}", "${yStr}", "${zStr}")`,
        rawLine: line.substring(0, 80)
      });
      continue;
    }

    if (Math.abs(x) > 1e5 || Math.abs(y) > 1e5 || Math.abs(z) > 1e5) {
      invalidCoordinateLines.push({
        lineNumber: lineNo,
        reason: `Współrzędne poza fizycznym zakresem białkowym (>100000 Å): (${x}, ${y}, ${z})`,
        rawLine: line.substring(0, 80)
      });
      continue;
    }

    const bStr = line.length >= 66 ? line.substring(60, 66).trim() : '';
    const bFactor = bStr && Number.isFinite(Number(bStr)) ? Number(bStr) : 0.0;

    const modelKey = explicitModelRecordsCount > 0 ? currentModel : 1;
    let bucket = atomsByModel.get(modelKey);
    if (!bucket) {
      bucket = [];
      atomsByModel.set(modelKey, bucket);
    }
    bucket.push({
      modelNum: modelKey,
      lineNumber: lineNo,
      recordType,
      atomName,
      altLoc,
      rawResName,
      chainId,
      resSeq,
      iCode,
      x,
      y,
      z,
      bFactor
    });
  }

  const modelKeys = Array.from(atomsByModel.keys()).sort((a, b) => a - b);
  const totalModelsDetected = explicitModelRecordsCount > 0 ? Math.max(explicitModelRecordsCount, modelKeys.length) : (modelKeys.length > 0 ? 1 : 0);
  const selectedModelNumber = atomsByModel.has(targetModel)
    ? targetModel
    : (modelKeys[0] ?? 1);

  const modelAtoms = atomsByModel.get(selectedModelNumber) ?? [];

  // Discover available chains in selected model
  const chainOrder: string[] = [];
  const chainSeen = new Set<string>();
  for (const a of modelAtoms) {
    if (!chainSeen.has(a.chainId)) {
      chainSeen.add(a.chainId);
      chainOrder.push(a.chainId);
    }
  }

  let selectedChainId = options?.chainId?.trim().toUpperCase() || '';
  if (!selectedChainId) {
    // Pick first chain that has at least one CA atom, or first chain overall
    let bestChain = chainOrder[0] || 'A';
    for (const cid of chainOrder) {
      if (modelAtoms.some(a => a.chainId === cid && a.atomName === 'CA')) {
        bestChain = cid;
        break;
      }
    }
    selectedChainId = bestChain;
  }

  const chainAtoms = modelAtoms.filter(a => a.chainId.toUpperCase() === selectedChainId);

  // Group into ordered residues along the selected chain
  const residueMap = new Map<string, {
    chainId: string;
    resSeq: number;
    iCode: string;
    rawResName: string;
    atoms: Map<string, Vec3Tuple>;
    bSum: number;
    bCnt: number;
  }>();
  const residueOrder: string[] = [];
  let duplicateAtomIgnoredCount = 0;

  for (const a of chainAtoms) {
    const key = `${a.chainId}:${a.resSeq}:${a.iCode}`;
    let entry = residueMap.get(key);
    if (!entry) {
      entry = {
        chainId: a.chainId,
        resSeq: a.resSeq,
        iCode: a.iCode,
        rawResName: a.rawResName,
        atoms: new Map<string, Vec3Tuple>(),
        bSum: 0,
        bCnt: 0
      };
      residueMap.set(key, entry);
      residueOrder.push(key);
    }
    if (entry.atoms.has(a.atomName)) {
      duplicateAtomIgnoredCount++;
      continue;
    }
    entry.atoms.set(a.atomName, [a.x, a.y, a.z]);
    entry.bSum += a.bFactor;
    entry.bCnt += 1;
  }

  const residues: ParsedPdbResidue[] = [];
  const residuesWithCa: ParsedPdbResidue[] = [];
  const residuesMissingCa: { resSeq: number; resName: string }[] = [];
  const unknownResidues: { resSeq: number; rawResName: string }[] = [];
  const modifiedResiduesNormalized: { resSeq: number; rawResName: string; mappedTo: string }[] = [];

  residueOrder.forEach((key, idx) => {
    const item = residueMap.get(key)!;
    let canonical = item.rawResName;
    let isKnown = STANDARD_21_AMINO_ACIDS.has(canonical);

    if (!isKnown && MODIFIED_RESIDUE_CANONICAL_MAP[canonical]) {
      const mapped = MODIFIED_RESIDUE_CANONICAL_MAP[canonical];
      modifiedResiduesNormalized.push({
        resSeq: item.resSeq,
        rawResName: canonical,
        mappedTo: mapped
      });
      canonical = mapped;
      isKnown = true;
    }

    if (!isKnown) {
      unknownResidues.push({ resSeq: item.resSeq, rawResName: item.rawResName });
      canonical = 'UNK';
    }

    const ca = item.atoms.get('CA') ?? null;
    const parsedRes: ParsedPdbResidue = {
      indexInChain: idx,
      chainId: item.chainId,
      resSeq: item.resSeq,
      iCode: item.iCode,
      rawResName: item.rawResName,
      canonicalResName: canonical,
      isKnownAminoAcid: isKnown,
      hasCaAtom: ca !== null,
      caCoord: ca,
      atoms: item.atoms,
      bFactorMean: item.bCnt > 0 ? item.bSum / item.bCnt : 0
    };

    residues.push(parsedRes);
    if (ca !== null) {
      residuesWithCa.push(parsedRes);
    } else {
      residuesMissingCa.push({ resSeq: item.resSeq, resName: item.rawResName });
    }
  });

  // Detect numbering gaps in residuesWithCa
  const sequenceNumberingGaps: SequenceNumberingGap[] = [];
  for (let i = 0; i < residuesWithCa.length - 1; i++) {
    const cur = residuesWithCa[i];
    const nxt = residuesWithCa[i + 1];
    const step = nxt.resSeq - cur.resSeq;
    if (step > 1) {
      sequenceNumberingGaps.push({
        afterResSeq: cur.resSeq,
        beforeResSeq: nxt.resSeq,
        missingCountEstimate: step - 1
      });
    }
  }

  return {
    headerId,
    title: titleParts.join(' ').trim() || headerId,
    totalModelsDetected,
    selectedModelNumber,
    availableChains: chainOrder,
    selectedChainId,
    residues,
    residuesWithCa,
    residuesMissingCa,
    unknownResidues,
    modifiedResiduesNormalized,
    sequenceNumberingGaps,
    invalidCoordinateLines,
    altLocFilteredCount,
    duplicateAtomIgnoredCount,
    rawPdbHash
  };
}

// -----------------------------------------------------------------------------
// 3. SEQUENCE-BASED RESIDUE MAPPING (NEEDLEMAN-WUNSCH + STRICT IDENTITY CHECK)
// -----------------------------------------------------------------------------

export function mapResiduesBySequence(
  modelResiduesWithCa: ParsedPdbResidue[],
  refResiduesWithCa: ParsedPdbResidue[],
  allowSubsequenceAlignment: boolean = false
): SequenceMappingSummary {
  const mLen = modelResiduesWithCa.length;
  const rLen = refResiduesWithCa.length;

  if (mLen === 0 || rLen === 0) {
    return {
      pairedResidues: [],
      mismatchedPairs: [],
      unpairedModelResidues: modelResiduesWithCa.map(r => ({ resSeq: r.resSeq, resName: r.canonicalResName })),
      unpairedRefResidues: refResiduesWithCa.map(r => ({ resSeq: r.resSeq, resName: r.canonicalResName })),
      pairedCount: 0,
      exactMatchCount: 0,
      modelTotalCaResidues: mLen,
      refTotalCaResidues: rLen,
      modelCoveragePct: 0,
      referenceCoveragePct: 0,
      sequenceIdentityPct: 0,
      alignmentMode: allowSubsequenceAlignment ? 'SEMIGLOBAL_SUBSEQUENCE' : 'GLOBAL_NEEDLEMAN_WUNSCH'
    };
  }

  // Decide semi-global (free terminal gaps on reference) ONLY if explicitly requested
  // or if reference is a precursor > 1.5x longer than model (e.g. P01308 110aa vs Insulin Chain A 21aa)
  const useSemiGlobalOnRef = allowSubsequenceAlignment;

  const MATCH_SCORE = 5;
  const MISMATCH_SCORE = -4;
  const GAP_PENALTY = -3;

  const cols = rLen + 1;
  const dp = new Int32Array((mLen + 1) * cols);
  const trace = new Uint8Array((mLen + 1) * cols); // 1 = diag, 2 = up (gap in ref), 3 = left (gap in model)

  for (let i = 1; i <= mLen; i++) {
    dp[i * cols + 0] = i * GAP_PENALTY;
    trace[i * cols + 0] = 2;
  }
  for (let j = 1; j <= rLen; j++) {
    dp[0 * cols + j] = useSemiGlobalOnRef ? 0 : j * GAP_PENALTY;
    trace[0 * cols + j] = 3;
  }

  for (let i = 1; i <= mLen; i++) {
    const mAA = modelResiduesWithCa[i - 1].canonicalResName;
    for (let j = 1; j <= rLen; j++) {
      const rAA = refResiduesWithCa[j - 1].canonicalResName;
      const isMatch = mAA !== 'UNK' && rAA !== 'UNK' && mAA === rAA;
      const diag = dp[(i - 1) * cols + (j - 1)] + (isMatch ? MATCH_SCORE : MISMATCH_SCORE);
      const up = dp[(i - 1) * cols + j] + GAP_PENALTY;
      const left = dp[i * cols + (j - 1)] + GAP_PENALTY;

      if (diag >= up && diag >= left) {
        dp[i * cols + j] = diag;
        trace[i * cols + j] = 1;
      } else if (up >= left) {
        dp[i * cols + j] = up;
        trace[i * cols + j] = 2;
      } else {
        dp[i * cols + j] = left;
        trace[i * cols + j] = 3;
      }
    }
  }

  let endI = mLen;
  let endJ = rLen;
  if (useSemiGlobalOnRef) {
    let bestScore = -1e9;
    for (let j = 0; j <= rLen; j++) {
      const s = dp[mLen * cols + j];
      if (s > bestScore) {
        bestScore = s;
        endJ = j;
      }
    }
  }

  const pairedReversed: ResiduePairMapping[] = [];
  const pairedModelIdxSet = new Set<number>();
  const pairedRefIdxSet = new Set<number>();

  let ci = endI;
  let cj = endJ;
  while (ci > 0 || cj > 0) {
    if (ci === 0) {
      cj--;
      continue;
    }
    if (cj === 0) {
      ci--;
      continue;
    }
    const dir = trace[ci * cols + cj];
    if (dir === 1) {
      const mRes = modelResiduesWithCa[ci - 1];
      const rRes = refResiduesWithCa[cj - 1];
      const identical =
        mRes.canonicalResName !== 'UNK' &&
        rRes.canonicalResName !== 'UNK' &&
        mRes.canonicalResName === rRes.canonicalResName;
      pairedReversed.push({
        modelChainIndex: ci - 1,
        refChainIndex: cj - 1,
        modelResSeq: mRes.resSeq,
        refResSeq: rRes.resSeq,
        modelResName: mRes.canonicalResName,
        refResName: rRes.canonicalResName,
        aminoAcidIdentical: identical
      });
      pairedModelIdxSet.add(ci - 1);
      pairedRefIdxSet.add(cj - 1);
      ci--;
      cj--;
    } else if (dir === 2) {
      ci--;
    } else {
      cj--;
    }
  }

  const pairedResidues = pairedReversed.reverse();
  const mismatchedPairs = pairedResidues.filter(p => !p.aminoAcidIdentical);
  const exactMatchCount = pairedResidues.length - mismatchedPairs.length;

  const unpairedModelResidues = modelResiduesWithCa
    .filter((_, idx) => !pairedModelIdxSet.has(idx))
    .map(r => ({ resSeq: r.resSeq, resName: r.canonicalResName }));

  const unpairedRefResidues = refResiduesWithCa
    .filter((_, idx) => !pairedRefIdxSet.has(idx))
    .map(r => ({ resSeq: r.resSeq, resName: r.canonicalResName }));

  return {
    pairedResidues,
    mismatchedPairs,
    unpairedModelResidues,
    unpairedRefResidues,
    pairedCount: pairedResidues.length,
    exactMatchCount,
    modelTotalCaResidues: mLen,
    refTotalCaResidues: rLen,
    modelCoveragePct: parseFloat(((pairedResidues.length / mLen) * 100).toFixed(2)),
    referenceCoveragePct: parseFloat(((pairedResidues.length / rLen) * 100).toFixed(2)),
    sequenceIdentityPct:
      pairedResidues.length > 0
        ? parseFloat(((exactMatchCount / pairedResidues.length) * 100).toFixed(2))
        : 0,
    alignmentMode: useSemiGlobalOnRef ? 'SEMIGLOBAL_SUBSEQUENCE' : 'GLOBAL_NEEDLEMAN_WUNSCH'
  };
}

// -----------------------------------------------------------------------------
// 4. EXACT 3x3 SINGULAR VALUE DECOMPOSITION (SVD) & KABSCH SUPERPOSITION
// -----------------------------------------------------------------------------

function mat3Mul(A: Mat3Tuple, B: Mat3Tuple): Mat3Tuple {
  const C: Mat3Tuple = [
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0]
  ];
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      C[i][j] = A[i][0] * B[0][j] + A[i][1] * B[1][j] + A[i][2] * B[2][j];
    }
  }
  return C;
}

function mat3Transpose(A: Mat3Tuple): Mat3Tuple {
  return [
    [A[0][0], A[1][0], A[2][0]],
    [A[0][1], A[1][1], A[2][1]],
    [A[0][2], A[1][2], A[2][2]]
  ];
}

export function mat3Determinant(M: Mat3Tuple): number {
  return (
    M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) -
    M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) +
    M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0])
  );
}

function vec3Dot(a: Vec3Tuple, b: Vec3Tuple): number {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

function vec3Cross(a: Vec3Tuple, b: Vec3Tuple): Vec3Tuple {
  return [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0]
  ];
}

function vec3Norm(a: Vec3Tuple): number {
  return Math.hypot(a[0], a[1], a[2]);
}

function vec3Normalize(a: Vec3Tuple): Vec3Tuple {
  const n = vec3Norm(a);
  if (n < 1e-15) return [0, 0, 0];
  return [a[0] / n, a[1] / n, a[2] / n];
}

function completeOrthonormalBasis(u0: Vec3Tuple, u1Candidate?: Vec3Tuple): [Vec3Tuple, Vec3Tuple, Vec3Tuple] {
  let e0 = vec3Normalize(u0);
  if (vec3Norm(e0) < 0.5) e0 = [1, 0, 0];

  let e1: Vec3Tuple = [0, 0, 0];
  if (u1Candidate) {
    const proj = vec3Dot(u1Candidate, e0);
    const orth: Vec3Tuple = [
      u1Candidate[0] - proj * e0[0],
      u1Candidate[1] - proj * e0[1],
      u1Candidate[2] - proj * e0[2]
    ];
    if (vec3Norm(orth) > 1e-9) {
      e1 = vec3Normalize(orth);
    }
  }

  if (vec3Norm(e1) < 0.5) {
    // Pick coordinate axis least aligned with e0
    const ax: Vec3Tuple =
      Math.abs(e0[0]) <= Math.abs(e0[1]) && Math.abs(e0[0]) <= Math.abs(e0[2])
        ? [1, 0, 0]
        : Math.abs(e0[1]) <= Math.abs(e0[2])
        ? [0, 1, 0]
        : [0, 0, 1];
    const proj = vec3Dot(ax, e0);
    e1 = vec3Normalize([ax[0] - proj * e0[0], ax[1] - proj * e0[1], ax[2] - proj * e0[2]]);
  }

  const e2 = vec3Normalize(vec3Cross(e0, e1));
  return [e0, e1, e2];
}

/**
 * Exact 3x3 Singular Value Decomposition H = U * diag(S) * V^T
 * via symmetric Jacobi eigendecomposition of H^T H with full rank-deficiency completion.
 */
export function computeSvd3x3(H: Mat3Tuple): Svd3x3Result {
  const Ht = mat3Transpose(H);
  const A = mat3Mul(Ht, H); // Symmetric positive semi-definite 3x3

  const V: Mat3Tuple = [
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]
  ];

  // Symmetric Jacobi sweeps on A
  const pairs: [number, number][] = [
    [0, 1],
    [0, 2],
    [1, 2]
  ];
  for (let sweep = 0; sweep < 60; sweep++) {
    const offDiag = Math.abs(A[0][1]) + Math.abs(A[0][2]) + Math.abs(A[1][2]);
    if (offDiag < 1e-15) break;

    for (const [p, q] of pairs) {
      const apq = A[p][q];
      if (Math.abs(apq) < 1e-16) continue;
      const app = A[p][p];
      const aqq = A[q][q];
      const tau = (aqq - app) / (2.0 * apq);
      const t = tau >= 0
        ? 1.0 / (tau + Math.hypot(1.0, tau))
        : -1.0 / (-tau + Math.hypot(1.0, tau));
      const c = 1.0 / Math.hypot(1.0, t);
      const s = t * c;

      A[p][p] = app - t * apq;
      A[q][q] = aqq + t * apq;
      A[p][q] = 0.0;
      A[q][p] = 0.0;

      for (let r = 0; r < 3; r++) {
        if (r !== p && r !== q) {
          const arp = A[r][p];
          const arq = A[r][q];
          A[r][p] = c * arp - s * arq;
          A[p][r] = A[r][p];
          A[r][q] = s * arp + c * arq;
          A[q][r] = A[r][q];
        }
      }

      for (let r = 0; r < 3; r++) {
        const vrp = V[r][p];
        const vrq = V[r][q];
        V[r][p] = c * vrp - s * vrq;
        V[r][q] = s * vrp + c * vrq;
      }
    }
  }

  // Sort eigenvalues descending: sigma_0 >= sigma_1 >= sigma_2 >= 0
  const order = [0, 1, 2].sort((a, b) => A[b][b] - A[a][a]);
  const sVals: Vec3Tuple = [
    Math.sqrt(Math.max(0, A[order[0]][order[0]])),
    Math.sqrt(Math.max(0, A[order[1]][order[1]])),
    Math.sqrt(Math.max(0, A[order[2]][order[2]]))
  ];

  let vCols: [Vec3Tuple, Vec3Tuple, Vec3Tuple] = [
    [V[0][order[0]], V[1][order[0]], V[2][order[0]]],
    [V[0][order[1]], V[1][order[1]], V[2][order[1]]],
    [V[0][order[2]], V[1][order[2]], V[2][order[2]]]
  ];

  // Ensure V is strictly right-handed orthonormal (det(V) = +1)
  vCols = completeOrthonormalBasis(vCols[0], vCols[1]);
  const Vsorted: Mat3Tuple = [
    [vCols[0][0], vCols[1][0], vCols[2][0]],
    [vCols[0][1], vCols[1][1], vCols[2][1]],
    [vCols[0][2], vCols[1][2], vCols[2][2]]
  ];

  // Construct U columns from H * v_k / sigma_k
  const hv = (v: Vec3Tuple): Vec3Tuple => [
    H[0][0] * v[0] + H[0][1] * v[1] + H[0][2] * v[2],
    H[1][0] * v[0] + H[1][1] * v[1] + H[1][2] * v[2],
    H[2][0] * v[0] + H[2][1] * v[1] + H[2][2] * v[2]
  ];

  const epsRank = Math.max(1e-12, sVals[0] * 1e-11);
  let u0: Vec3Tuple = sVals[0] > epsRank ? vec3Normalize(hv(vCols[0])) : [1, 0, 0];
  let u1Raw: Vec3Tuple | undefined = sVals[1] > epsRank ? vec3Normalize(hv(vCols[1])) : undefined;
  let [uCol0, uCol1, uCol2] = completeOrthonormalBasis(u0, u1Raw);

  if (sVals[2] > epsRank) {
    const hv2 = hv(vCols[2]);
    // Align sign of uCol2 with H * vCols[2] so sigma_2 >= 0 satisfies H v_2 = sigma_2 u_2
    if (vec3Dot(uCol2, hv2) < 0) {
      uCol2 = [-uCol2[0], -uCol2[1], -uCol2[2]];
    }
  }

  const U: Mat3Tuple = [
    [uCol0[0], uCol1[0], uCol2[0]],
    [uCol0[1], uCol1[1], uCol2[1]],
    [uCol0[2], uCol1[2], uCol2[2]]
  ];

  // Check reconstruction error ||U * diag(S) * V^T - H||_max and orthogonality
  const US: Mat3Tuple = [
    [U[0][0] * sVals[0], U[0][1] * sVals[1], U[0][2] * sVals[2]],
    [U[1][0] * sVals[0], U[1][1] * sVals[1], U[1][2] * sVals[2]],
    [U[2][0] * sVals[0], U[2][1] * sVals[1], U[2][2] * sVals[2]]
  ];
  const recon = mat3Mul(US, mat3Transpose(Vsorted));
  let maxReconErr = 0;
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      maxReconErr = Math.max(maxReconErr, Math.abs(recon[i][j] - H[i][j]));
    }
  }

  const UtU = mat3Mul(mat3Transpose(U), U);
  const VtV = mat3Mul(mat3Transpose(Vsorted), Vsorted);
  let maxOrthErr = 0;
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      const target = i === j ? 1 : 0;
      maxOrthErr = Math.max(maxOrthErr, Math.abs(UtU[i][j] - target), Math.abs(VtV[i][j] - target));
    }
  }

  return {
    U,
    S: sVals,
    V: Vsorted,
    reconstructionMaxError: maxReconErr,
    orthogonalityMaxError: maxOrthErr
  };
}

/**
 * Optimal rigid-body superposition of modelCoords onto referenceCoords via Kabsch SVD.
 * Strictly enforces proper rotation det(R) = +1 (no mirror reflections) and exact order:
 *   x_aligned = R * (x_model - mu_model) + mu_ref = R * x_model + t
 * Computes Cα RMSD = sqrt(mean(sum((X_ref - X_aligned)^2, axis=1))).
 */
export function kabschSuperposeSvd(
  modelCoords: Vec3Tuple[],
  referenceCoords: Vec3Tuple[]
): KabschSuperpositionResult {
  const n = Math.min(modelCoords.length, referenceCoords.length);
  if (n < 3) {
    throw new Error(`Algorytm Kabscha wymaga co najmniej 3 sparowanych atomów Cα (otrzymano ${n}).`);
  }

  // 1. Centroids in Å
  let mx = 0, my = 0, mz = 0;
  let rx = 0, ry = 0, rz = 0;
  for (let i = 0; i < n; i++) {
    mx += modelCoords[i][0];
    my += modelCoords[i][1];
    mz += modelCoords[i][2];
    rx += referenceCoords[i][0];
    ry += referenceCoords[i][1];
    rz += referenceCoords[i][2];
  }
  const modelCentroid: Vec3Tuple = [mx / n, my / n, mz / n];
  const referenceCentroid: Vec3Tuple = [rx / n, ry / n, rz / n];

  // 2. Covariance matrix H = P_centered^T * Q_centered
  // where P_centered is model (N x 3) and Q_centered is reference (N x 3)
  const H: Mat3Tuple = [
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0]
  ];
  for (let k = 0; k < n; k++) {
    const px = modelCoords[k][0] - modelCentroid[0];
    const py = modelCoords[k][1] - modelCentroid[1];
    const pz = modelCoords[k][2] - modelCentroid[2];
    const qx = referenceCoords[k][0] - referenceCentroid[0];
    const qy = referenceCoords[k][1] - referenceCentroid[1];
    const qz = referenceCoords[k][2] - referenceCentroid[2];

    H[0][0] += px * qx; H[0][1] += px * qy; H[0][2] += px * qz;
    H[1][0] += py * qx; H[1][1] += py * qy; H[1][2] += py * qz;
    H[2][0] += pz * qx; H[2][1] += pz * qy; H[2][2] += pz * qz;
  }

  // 3. SVD: H = U * diag(S) * V^T
  const svd = computeSvd3x3(H);
  const Ut = mat3Transpose(svd.U);
  const VUt = mat3Mul(svd.V, Ut);
  const detVUt = mat3Determinant(VUt);
  const signD = detVUt < 0 ? -1.0 : 1.0;

  // R = V * diag(1, 1, d) * U^T maps column vectors (x_model - mu_model) -> (x_ref - mu_ref)
  const VCorr: Mat3Tuple = [
    [svd.V[0][0], svd.V[0][1], svd.V[0][2] * signD],
    [svd.V[1][0], svd.V[1][1], svd.V[1][2] * signD],
    [svd.V[2][0], svd.V[2][1], svd.V[2][2] * signD]
  ];
  const R = mat3Mul(VCorr, Ut);
  const detR = mat3Determinant(R);

  // Translation vector t = mu_ref - R * mu_model so x_aligned = R * x_model + t
  const rMuModel: Vec3Tuple = [
    R[0][0] * modelCentroid[0] + R[0][1] * modelCentroid[1] + R[0][2] * modelCentroid[2],
    R[1][0] * modelCentroid[0] + R[1][1] * modelCentroid[1] + R[1][2] * modelCentroid[2],
    R[2][0] * modelCentroid[0] + R[2][1] * modelCentroid[1] + R[2][2] * modelCentroid[2]
  ];
  const translationVec: Vec3Tuple = [
    referenceCentroid[0] - rMuModel[0],
    referenceCentroid[1] - rMuModel[1],
    referenceCentroid[2] - rMuModel[2]
  ];

  // 4. Transform coordinates & compute exact RMSD = sqrt(mean(sum((X_ref - X_aligned)^2, axis=1)))
  const alignedModelCaCoords: Vec3Tuple[] = new Array(n);
  let sumSq = 0.0;
  for (let i = 0; i < n; i++) {
    const xm = modelCoords[i][0] - modelCentroid[0];
    const ym = modelCoords[i][1] - modelCentroid[1];
    const zm = modelCoords[i][2] - modelCentroid[2];

    const ax = R[0][0] * xm + R[0][1] * ym + R[0][2] * zm + referenceCentroid[0];
    const ay = R[1][0] * xm + R[1][1] * ym + R[1][2] * zm + referenceCentroid[1];
    const az = R[2][0] * xm + R[2][1] * ym + R[2][2] * zm + referenceCentroid[2];

    alignedModelCaCoords[i] = [ax, ay, az];

    const dx = referenceCoords[i][0] - ax;
    const dy = referenceCoords[i][1] - ay;
    const dz = referenceCoords[i][2] - az;
    sumSq += dx * dx + dy * dy + dz * dz;
  }

  const caRmsdAngstrom = Math.sqrt(sumSq / n);

  return {
    rotationMatrix: R,
    translationVec,
    modelCentroid,
    referenceCentroid,
    determinantR: detR,
    improperReflectionCorrected: signD < 0,
    svdSingularValues: svd.S,
    svdReconstructionError: svd.reconstructionMaxError,
    alignedModelCaCoords,
    referenceCaCoords: referenceCoords.slice(0, n),
    caRmsdAngstrom
  };
}

// -----------------------------------------------------------------------------
// 5. INDEPENDENT SECONDARY RMSD RE-COMPUTATION (ZERO SHARED STATE / ZERO CACHE)
// -----------------------------------------------------------------------------

/**
 * Independently recomputes Cα RMSD from raw uncentered PDB coordinates without reading
 * X_aligned or caRmsdAngstrom from kabschSuperposeSvd:
 *  (1) Direct raw affine evaluation x'_i = R * x_raw[i] + t against raw reference XYZ, AND
 *  (2) Independent 4x4 Davenport quaternion eigenvalue bound sqrt(max(0, (E0 - 2*lambda_max)/N)).
 */
export function recomputeCaRmsdIndependently(
  rawModelCa: Vec3Tuple[],
  rawRefCa: Vec3Tuple[],
  rotationMatrix: Mat3Tuple,
  translationVec: Vec3Tuple
): {
  directTransformedRmsdAngstrom: number;
  quaternionEigenRmsdAngstrom: number;
} {
  const n = Math.min(rawModelCa.length, rawRefCa.length);
  if (n === 0) {
    return { directTransformedRmsdAngstrom: 0, quaternionEigenRmsdAngstrom: 0 };
  }

  // Method A: Direct affine transformation R * x_raw + t on uncentered coordinates
  let sqAcc = 0.0;
  for (let i = 0; i < n; i++) {
    const mx = rawModelCa[i][0];
    const my = rawModelCa[i][1];
    const mz = rawModelCa[i][2];

    const tx = rotationMatrix[0][0] * mx + rotationMatrix[0][1] * my + rotationMatrix[0][2] * mz + translationVec[0];
    const ty = rotationMatrix[1][0] * mx + rotationMatrix[1][1] * my + rotationMatrix[1][2] * mz + translationVec[1];
    const tz = rotationMatrix[2][0] * mx + rotationMatrix[2][1] * my + rotationMatrix[2][2] * mz + translationVec[2];

    const ex = rawRefCa[i][0] - tx;
    const ey = rawRefCa[i][1] - ty;
    const ez = rawRefCa[i][2] - tz;
    sqAcc += ex * ex + ey * ey + ez * ez;
  }
  const directTransformedRmsdAngstrom = Math.sqrt(sqAcc / n);

  // Method B: Independent Horn 4x4 symmetric Davenport quaternion matrix eigenvalue (completely SVD-free!)
  let cmx = 0, cmy = 0, cmz = 0;
  let crx = 0, cry = 0, crz = 0;
  for (let i = 0; i < n; i++) {
    cmx += rawModelCa[i][0]; cmy += rawModelCa[i][1]; cmz += rawModelCa[i][2];
    crx += rawRefCa[i][0]; cry += rawRefCa[i][1]; crz += rawRefCa[i][2];
  }
  cmx /= n; cmy /= n; cmz /= n;
  crx /= n; cry /= n; crz /= n;

  let E0 = 0;
  let Sxx = 0, Sxy = 0, Sxz = 0;
  let Syx = 0, Syy = 0, Syz = 0;
  let Szx = 0, Szy = 0, Szz = 0;
  for (let i = 0; i < n; i++) {
    const px = rawModelCa[i][0] - cmx;
    const py = rawModelCa[i][1] - cmy;
    const pz = rawModelCa[i][2] - cmz;
    const qx = rawRefCa[i][0] - crx;
    const qy = rawRefCa[i][1] - cry;
    const qz = rawRefCa[i][2] - crz;
    E0 += px * px + py * py + pz * pz + qx * qx + qy * qy + qz * qz;
    Sxx += px * qx; Sxy += px * qy; Sxz += px * qz;
    Syx += py * qx; Syy += py * qy; Syz += py * qz;
    Szx += pz * qx; Szy += pz * qy; Szz += pz * qz;
  }

  const K: number[][] = [
    [Sxx + Syy + Szz, Syz - Szy, Szx - Sxz, Sxy - Syx],
    [Syz - Szy, Sxx - Syy - Szz, Sxy + Syx, Szx + Sxz],
    [Szx - Sxz, Sxy + Syx, -Sxx + Syy - Szz, Syz + Szy],
    [Sxy - Syx, Szx + Sxz, Syz + Szy, -Sxx - Syy + Szz]
  ];

  // Symmetric 4x4 Jacobi eigensolver for K
  for (let sweep = 0; sweep < 50; sweep++) {
    let maxOff = 0;
    let p = 0, q = 1;
    for (let i = 0; i < 4; i++) {
      for (let j = i + 1; j < 4; j++) {
        const v = Math.abs(K[i][j]);
        if (v > maxOff) {
          maxOff = v;
          p = i;
          q = j;
        }
      }
    }
    if (maxOff < 1e-14) break;
    const kpq = K[p][q];
    const tau = (K[q][q] - K[p][p]) / (2.0 * kpq);
    const t = tau >= 0 ? 1.0 / (tau + Math.hypot(1.0, tau)) : -1.0 / (-tau + Math.hypot(1.0, tau));
    const c = 1.0 / Math.hypot(1.0, t);
    const s = t * c;
    const kpp = K[p][p];
    const kqq = K[q][q];
    K[p][p] = kpp - t * kpq;
    K[q][q] = kqq + t * kpq;
    K[p][q] = 0;
    K[q][p] = 0;
    for (let r = 0; r < 4; r++) {
      if (r !== p && r !== q) {
        const krp = K[r][p];
        const krq = K[r][q];
        K[r][p] = c * krp - s * krq;
        K[p][r] = K[r][p];
        K[r][q] = s * krp + c * krq;
        K[q][r] = K[r][q];
      }
    }
  }

  const lambdaMax = Math.max(K[0][0], K[1][1], K[2][2], K[3][3]);
  const residual = Math.max(0, (E0 - 2.0 * lambdaMax) / n);
  const quaternionEigenRmsdAngstrom = Math.sqrt(residual);

  return {
    directTransformedRmsdAngstrom,
    quaternionEigenRmsdAngstrom
  };
}

// -----------------------------------------------------------------------------
// 6. PHYSICAL GEOMETRY & STERIC PLAUSIBILITY CHECK
// -----------------------------------------------------------------------------

export function verifyPhysicalProteinGeometry(parsed: ParsedPdbStructure): PhysicalGeometryCheckResult {
  const resList = parsed.residuesWithCa;
  const n = resList.length;
  const reasons: string[] = [];

  if (parsed.invalidCoordinateLines.length > 0) {
    reasons.push(`Wykryto ${parsed.invalidCoordinateLines.length} uszkodzonych rekordów współrzędnych XYZ.`);
  }
  if (n < 3) {
    reasons.push(`Zbyt mało reszt z atomem Cα (${n} < 3).`);
    return {
      geometryValid: false,
      checkedResidueCount: n,
      caBondCount: 0,
      meanCaCaDistAngstrom: 0,
      minCaCaDistAngstrom: 0,
      maxCaCaDistAngstrom: 0,
      caBondViolations: [],
      nonBondedCaClashes: [],
      heavyAtomStericClashesCount: 0,
      collapsedChainDetected: true,
      reasons
    };
  }

  // Check if all coordinates are collapsed into a tiny ball (< 1.0 Å spread)
  let cx = 0, cy = 0, cz = 0;
  for (const r of resList) {
    cx += r.caCoord![0];
    cy += r.caCoord![1];
    cz += r.caCoord![2];
  }
  cx /= n; cy /= n; cz /= n;
  let rgSq = 0;
  for (const r of resList) {
    rgSq +=
      (r.caCoord![0] - cx) ** 2 +
      (r.caCoord![1] - cy) ** 2 +
      (r.caCoord![2] - cz) ** 2;
  }
  const rg = Math.sqrt(rgSq / n);
  const collapsedChainDetected = rg < 1.2;
  if (collapsedChainDetected) {
    reasons.push(`Wykryto zdegenerowany (zapadnięty) łańcuch współrzędnych: Rg = ${rg.toFixed(3)} Å.`);
  }

  // Check sequential Cα-Cα peptide bond lengths (2.60 Å .. 4.60 Å allows cis-Pro 2.9 Å and trans 3.8 Å)
  const caBondViolations: PhysicalGeometryCheckResult['caBondViolations'] = [];
  let sumBond = 0;
  let minBond = Infinity;
  let maxBond = 0;
  let bondCount = 0;

  for (let i = 0; i < n - 1; i++) {
    const a = resList[i];
    const b = resList[i + 1];
    if (b.resSeq - a.resSeq !== 1) continue; // Skip chain numbering gaps
    const d = Math.hypot(
      b.caCoord![0] - a.caCoord![0],
      b.caCoord![1] - a.caCoord![1],
      b.caCoord![2] - a.caCoord![2]
    );
    sumBond += d;
    minBond = Math.min(minBond, d);
    maxBond = Math.max(maxBond, d);
    bondCount++;

    if (d < 2.60) {
      caBondViolations.push({
        resSeqA: a.resSeq,
        resSeqB: b.resSeq,
        distanceAngstrom: parseFloat(d.toFixed(3)),
        issue: 'TOO_SHORT'
      });
    } else if (d > 4.60) {
      caBondViolations.push({
        resSeqA: a.resSeq,
        resSeqB: b.resSeq,
        distanceAngstrom: parseFloat(d.toFixed(3)),
        issue: 'TOO_LONG'
      });
    }
  }

  if (caBondViolations.length > 0) {
    reasons.push(
      `Naruszenie ciągłości wiązań szkieletu Cα-Cα: ${caBondViolations.length} wiązań poza fizycznym zakresem [2.60, 4.60] Å.`
    );
  }

  // Check non-bonded Cα-Cα steric clashes (|i - j| >= 2, d < 2.50 Å)
  const nonBondedCaClashes: PhysicalGeometryCheckResult['nonBondedCaClashes'] = [];
  for (let i = 0; i < n; i++) {
    for (let j = i + 2; j < n; j++) {
      const a = resList[i];
      const b = resList[j];
      const d = Math.hypot(
        b.caCoord![0] - a.caCoord![0],
        b.caCoord![1] - a.caCoord![1],
        b.caCoord![2] - a.caCoord![2]
      );
      if (d < 2.50) {
        nonBondedCaClashes.push({
          resSeqA: a.resSeq,
          resSeqB: b.resSeq,
          distanceAngstrom: parseFloat(d.toFixed(3))
        });
      }
    }
  }

  if (nonBondedCaClashes.length > 0) {
    reasons.push(`Wykryto ${nonBondedCaClashes.length} ciężkich kolizji sterycznych Cα-Cα (< 2.50 Å).`);
  }

  // Check heavy-atom clashes between residues separated by >= 2 sequence positions (< 1.55 Å, excluding S-S/Se-Se)
  let heavyAtomStericClashesCount = 0;
  for (let i = 0; i < n; i++) {
    for (let j = i + 2; j < n; j++) {
      const a = resList[i];
      const b = resList[j];
      // Quick Cα bounding sphere filter (14 Å)
      const caDist = Math.hypot(
        b.caCoord![0] - a.caCoord![0],
        b.caCoord![1] - a.caCoord![1],
        b.caCoord![2] - a.caCoord![2]
      );
      if (caDist > 14.0) continue;

      a.atoms.forEach((coordA, nameA) => {
        b.atoms.forEach((coordB, nameB) => {
          const isDisulfidePair =
            (nameA === 'SG' || nameA === 'SE') && (nameB === 'SG' || nameB === 'SE');
          if (isDisulfidePair) return;
          const d = Math.hypot(
            coordB[0] - coordA[0],
            coordB[1] - coordA[1],
            coordB[2] - coordA[2]
          );
          if (d < 1.55) {
            heavyAtomStericClashesCount++;
          }
        });
      });
    }
  }

  if (heavyAtomStericClashesCount > 0) {
    reasons.push(`Wykryto ${heavyAtomStericClashesCount} nakładających się par ciężkich atomów (< 1.55 Å).`);
  }

  const geometryValid =
    parsed.invalidCoordinateLines.length === 0 &&
    !collapsedChainDetected &&
    caBondViolations.length === 0 &&
    nonBondedCaClashes.length === 0 &&
    heavyAtomStericClashesCount === 0;

  return {
    geometryValid,
    checkedResidueCount: n,
    caBondCount: bondCount,
    meanCaCaDistAngstrom: bondCount > 0 ? parseFloat((sumBond / bondCount).toFixed(3)) : 0,
    minCaCaDistAngstrom: bondCount > 0 && Number.isFinite(minBond) ? parseFloat(minBond.toFixed(3)) : 0,
    maxCaCaDistAngstrom: bondCount > 0 ? parseFloat(maxBond.toFixed(3)) : 0,
    caBondViolations,
    nonBondedCaClashes,
    heavyAtomStericClashesCount,
    collapsedChainDetected,
    reasons
  };
}

// -----------------------------------------------------------------------------
// 7. MAIN ENTRY POINT: INDEPENDENT GROUND-TRUTH VALIDATOR
// -----------------------------------------------------------------------------

export function validateAgainstGroundTruth(
  modelPdbText: string,
  referencePdbText: string | null | undefined,
  options?: GroundTruthValidationOptions
): GroundTruthValidationReport {
  const iterationId = String(options?.iterationId ?? `iter-${Date.now()}`);
  const parameterVersion = options?.parameterVersion ?? 'v1.0.0-default';
  const timestampIso = new Date().toISOString();
  const numericalTol = options?.numericalToleranceAngstrom ?? 1e-7;

  const criteria: GroundTruthPassCriteria = {
    ...DEFAULT_HIGH_ACCURACY_CRITERIA,
    ...(options?.criteria ?? {}),
    minReferenceCoveragePct:
      options?.criteria?.minReferenceCoveragePct !== undefined
        ? options.criteria.minReferenceCoveragePct
        : options?.allowSubsequenceAlignment
        ? 0.0
        : DEFAULT_HIGH_ACCURACY_CRITERIA.minReferenceCoveragePct
  };

  const warnings: string[] = [];

  // Parse model PDB
  const parsedModel = parsePdbForGroundTruth(modelPdbText || '', {
    chainId: options?.modelChainId,
    modelNumber: options?.modelNumber
  });
  const modelGeometry = verifyPhysicalProteinGeometry(parsedModel);

  if (parsedModel.residuesMissingCa.length > 0) {
    warnings.push(
      `Model zawiera ${parsedModel.residuesMissingCa.length} reszt bez atomu Cα.`
    );
  }
  if (parsedModel.unknownResidues.length > 0) {
    warnings.push(
      `Model zawiera ${parsedModel.unknownResidues.length} nieznanych reszt (${parsedModel.unknownResidues.map(u => `${u.rawResName}${u.resSeq}`).join(', ')}).`
    );
  }
  if (parsedModel.sequenceNumberingGaps.length > 0) {
    warnings.push(
      `Model zawiera ${parsedModel.sequenceNumberingGaps.length} luk w numeracji reszt.`
    );
  }

  // Case 0: Model itself is unreadable or empty
  if (parsedModel.invalidCoordinateLines.length > 0 || parsedModel.residuesWithCa.length < 3) {
    const errReason =
      parsedModel.invalidCoordinateLines.length > 0
        ? `Błąd współrzędnych modelu PDB: ${parsedModel.invalidCoordinateLines[0].reason} (linia ${parsedModel.invalidCoordinateLines[0].lineNumber})`
        : `Model PDB nie zawiera wystarczającej liczby poprawnych atomów Cα (znaleziono ${parsedModel.residuesWithCa.length}).`;
    return {
      valid: false,
      hasGroundTruthReference: Boolean(referencePdbText && referencePdbText.trim()),
      status: 'VALIDATION_ERROR',
      passedTargetCriterion: false,
      passedIntermediateCriterion: false,
      failureReason: errReason,
      warnings,
      iterationId,
      parameterVersion,
      timestampIso,
      modelArtifactPath: options?.modelArtifactPath ?? null,
      modelPdbHash: parsedModel.rawPdbHash,
      referencePdbHash: referencePdbText ? computeDeterministicArtifactHash(referencePdbText) : null,
      referenceId: options?.referenceId ?? null,
      referenceName: options?.referenceName ?? null,
      referenceSource: options?.referenceSource ?? 'NONE',
      selectedModelChain: parsedModel.selectedChainId,
      availableModelChains: parsedModel.availableChains,
      selectedReferenceChain: null,
      availableReferenceChains: [],
      selectedModelNumber: parsedModel.selectedModelNumber,
      totalModelStructures: parsedModel.totalModelsDetected,
      selectedReferenceModelNumber: null,
      totalReferenceStructures: 0,
      caRmsdAngstrom: null,
      independentCaRmsdAngstrom: null,
      quaternionCrossCheckCaRmsdAngstrom: null,
      numericalAgreementDeltaAngstrom: null,
      independentVerificationPassed: false,
      backboneRmsdAngstrom: null,
      allHeavyAtomRmsdAngstrom: null,
      tmScore: null,
      gdtTsPct: null,
      meanReferencePlddtOrBfactor: null,
      pairedResidueCount: 0,
      modelResidueCount: parsedModel.residuesWithCa.length,
      referenceResidueCount: 0,
      modelCoveragePct: 0,
      referenceCoveragePct: 0,
      sequenceIdentityPct: 0,
      mismatchedResidueCount: 0,
      superposition: null,
      sequenceMapping: null,
      modelGeometry,
      perResidueDeviations: [],
      alignedModelCaCoords: [],
      pairedReferenceCaCoords: [],
      criteriaApplied: criteria
    };
  }

  // Case 1: No external Ground-Truth reference provided -> NEVER fake a 0.0 Å RMSD!
  if (!referencePdbText || !referencePdbText.trim()) {
    return {
      valid: modelGeometry.geometryValid,
      hasGroundTruthReference: false,
      status: modelGeometry.geometryValid ? 'NO_REFERENCE' : 'GEOMETRY_VIOLATION',
      passedTargetCriterion: false,
      passedIntermediateCriterion: false,
      failureReason: modelGeometry.geometryValid
        ? 'RMSD referencyjny niedostępny — nie podano eksperymentalnej ani testowej struktury referencyjnej PDB.'
        : modelGeometry.reasons.join(' '),
      warnings,
      iterationId,
      parameterVersion,
      timestampIso,
      modelArtifactPath: options?.modelArtifactPath ?? null,
      modelPdbHash: parsedModel.rawPdbHash,
      referencePdbHash: null,
      referenceId: null,
      referenceName: null,
      referenceSource: 'NONE',
      selectedModelChain: parsedModel.selectedChainId,
      availableModelChains: parsedModel.availableChains,
      selectedReferenceChain: null,
      availableReferenceChains: [],
      selectedModelNumber: parsedModel.selectedModelNumber,
      totalModelStructures: parsedModel.totalModelsDetected,
      selectedReferenceModelNumber: null,
      totalReferenceStructures: 0,
      caRmsdAngstrom: null,
      independentCaRmsdAngstrom: null,
      quaternionCrossCheckCaRmsdAngstrom: null,
      numericalAgreementDeltaAngstrom: null,
      independentVerificationPassed: false,
      backboneRmsdAngstrom: null,
      allHeavyAtomRmsdAngstrom: null,
      tmScore: null,
      gdtTsPct: null,
      meanReferencePlddtOrBfactor: null,
      pairedResidueCount: 0,
      modelResidueCount: parsedModel.residuesWithCa.length,
      referenceResidueCount: 0,
      modelCoveragePct: 0,
      referenceCoveragePct: 0,
      sequenceIdentityPct: 0,
      mismatchedResidueCount: 0,
      superposition: null,
      sequenceMapping: null,
      modelGeometry,
      perResidueDeviations: [],
      alignedModelCaCoords: [],
      pairedReferenceCaCoords: [],
      criteriaApplied: criteria
    };
  }

  // Parse reference PDB
  // If referenceChainId is not specified and multiple chains exist in reference (e.g. 4INS / 3I40 Insulin A & B chains),
  // inspect all available chains if default chain doesn't match sequence, and explicitly report selectedReferenceChain!
  let parsedRef = parsePdbForGroundTruth(referencePdbText, {
    chainId: options?.referenceChainId,
    modelNumber: options?.referenceModelNumber
  });

  if (!options?.referenceChainId && parsedRef.availableChains.length > 1) {
    let bestChainParsed = parsedRef;
    let bestMap = mapResiduesBySequence(
      parsedModel.residuesWithCa,
      parsedRef.residuesWithCa,
      Boolean(options?.allowSubsequenceAlignment)
    );
    for (const cid of parsedRef.availableChains) {
      if (cid === parsedRef.selectedChainId) continue;
      const candParsed = parsePdbForGroundTruth(referencePdbText, {
        chainId: cid,
        modelNumber: options?.referenceModelNumber
      });
      const candMap = mapResiduesBySequence(
        parsedModel.residuesWithCa,
        candParsed.residuesWithCa,
        Boolean(options?.allowSubsequenceAlignment)
      );
      if (candMap.exactMatchCount > bestMap.exactMatchCount) {
        bestChainParsed = candParsed;
        bestMap = candMap;
      }
    }
    parsedRef = bestChainParsed;
  }

  if (parsedRef.totalModelsDetected > 1) {
    warnings.push(
      `Struktura referencyjna jest wielomodelowa (${parsedRef.totalModelsDetected} modeli NMR/ensemble); wybrano jawnie MODEL ${parsedRef.selectedModelNumber}.`
    );
  }
  if (parsedRef.residuesMissingCa.length > 0) {
    warnings.push(
      `Struktura referencyjna (łańcuch ${parsedRef.selectedChainId}) ma ${parsedRef.residuesMissingCa.length} reszt bez atomu Cα.`
    );
  }
  if (parsedRef.unknownResidues.length > 0) {
    warnings.push(
      `Struktura referencyjna zawiera nieznane reszty: ${parsedRef.unknownResidues.map(u => `${u.rawResName}${u.resSeq}`).join(', ')}.`
    );
  }
  if (parsedRef.sequenceNumberingGaps.length > 0) {
    warnings.push(
      `Struktura referencyjna zawiera ${parsedRef.sequenceNumberingGaps.length} luk w numeracji reszt (${parsedRef.sequenceNumberingGaps.map(g => `${g.afterResSeq}..${g.beforeResSeq}`).join(', ')}).`
    );
  }

  // Validate reference coordinates
  if (parsedRef.invalidCoordinateLines.length > 0 || parsedRef.residuesWithCa.length < 3) {
    const refErr =
      parsedRef.invalidCoordinateLines.length > 0
        ? `Błąd współrzędnych referencji PDB: ${parsedRef.invalidCoordinateLines[0].reason} (linia ${parsedRef.invalidCoordinateLines[0].lineNumber})`
        : `Referencja PDB (łańcuch ${parsedRef.selectedChainId}) zawiera zbyt mało atomów Cα (${parsedRef.residuesWithCa.length}).`;
    return {
      valid: false,
      hasGroundTruthReference: true,
      status: 'VALIDATION_ERROR',
      passedTargetCriterion: false,
      passedIntermediateCriterion: false,
      failureReason: refErr,
      warnings,
      iterationId,
      parameterVersion,
      timestampIso,
      modelArtifactPath: options?.modelArtifactPath ?? null,
      modelPdbHash: parsedModel.rawPdbHash,
      referencePdbHash: parsedRef.rawPdbHash,
      referenceId: options?.referenceId ?? parsedRef.headerId,
      referenceName: options?.referenceName ?? parsedRef.title,
      referenceSource: options?.referenceSource ?? 'CUSTOM_PDB',
      selectedModelChain: parsedModel.selectedChainId,
      availableModelChains: parsedModel.availableChains,
      selectedReferenceChain: parsedRef.selectedChainId,
      availableReferenceChains: parsedRef.availableChains,
      selectedModelNumber: parsedModel.selectedModelNumber,
      totalModelStructures: parsedModel.totalModelsDetected,
      selectedReferenceModelNumber: parsedRef.selectedModelNumber,
      totalReferenceStructures: parsedRef.totalModelsDetected,
      caRmsdAngstrom: null,
      independentCaRmsdAngstrom: null,
      quaternionCrossCheckCaRmsdAngstrom: null,
      numericalAgreementDeltaAngstrom: null,
      independentVerificationPassed: false,
      backboneRmsdAngstrom: null,
      allHeavyAtomRmsdAngstrom: null,
      tmScore: null,
      gdtTsPct: null,
      meanReferencePlddtOrBfactor: null,
      pairedResidueCount: 0,
      modelResidueCount: parsedModel.residuesWithCa.length,
      referenceResidueCount: parsedRef.residuesWithCa.length,
      modelCoveragePct: 0,
      referenceCoveragePct: 0,
      sequenceIdentityPct: 0,
      mismatchedResidueCount: 0,
      superposition: null,
      sequenceMapping: null,
      modelGeometry,
      perResidueDeviations: [],
      alignedModelCaCoords: [],
      pairedReferenceCaCoords: [],
      criteriaApplied: criteria
    };
  }

  // Guard against accidental self-comparison warning if hashes are identical
  if (parsedModel.rawPdbHash === parsedRef.rawPdbHash) {
    warnings.push(
      'UWAGA DIAGNOSTYCZNA: Plik modelu i plik referencji mają identyczny skrót kryptograficzny (porównanie struktury z samą sobą).'
    );
  }

  // Sequence-based residue mapping
  const seqMapping = mapResiduesBySequence(
    parsedModel.residuesWithCa,
    parsedRef.residuesWithCa,
    Boolean(options?.allowSubsequenceAlignment)
  );

  if (seqMapping.mismatchedPairs.length > 0) {
    const sampleMismatches = seqMapping.mismatchedPairs
      .slice(0, 5)
      .map(m => `${m.modelResName}${m.modelResSeq}≠${m.refResName}${m.refResSeq}`)
      .join(', ');
    warnings.push(
      `Niezgodność sekwencji w ${seqMapping.mismatchedPairs.length} parach reszt (np. ${sampleMismatches}).`
    );
  }

  if (seqMapping.unpairedModelResidues.length > 0) {
    warnings.push(
      `Niepokryte reszty modelu: ${seqMapping.unpairedModelResidues.length} z ${seqMapping.modelTotalCaResidues} (pokrycie modelu ${seqMapping.modelCoveragePct}%).`
    );
  }

  if (seqMapping.unpairedRefResidues.length > 0) {
    warnings.push(
      `Niepokryte reszty referencji: ${seqMapping.unpairedRefResidues.length} z ${seqMapping.refTotalCaResidues} (pokrycie referencji ${seqMapping.referenceCoveragePct}%).`
    );
  }

  const strictSeqMismatchBlocked =
    !options?.allowSequenceMismatches &&
    criteria.requireExactSequenceIdentity &&
    (seqMapping.mismatchedPairs.length > 0 || parsedModel.unknownResidues.length > 0 || parsedRef.unknownResidues.length > 0);

  if (strictSeqMismatchBlocked || seqMapping.pairedCount < 3) {
    const mismatchMsg =
      seqMapping.pairedCount < 3
        ? `Nie można wiarygodnie sparować struktur: tylko ${seqMapping.pairedCount} sparowanych atomów Cα (minimum 3).`
        : seqMapping.mismatchedPairs.length > 0
        ? `Odrzucono dopasowanie niezgodnych sekwencji: ${seqMapping.mismatchedPairs.length} niezgodnych aminokwasów na ${seqMapping.pairedCount} par (identyczność ${seqMapping.sequenceIdentityPct}%).`
        : `Wykryto niezidentyfikowane reszty (UNK) uniemożliwiające jednoznaczną weryfikację tożsamości sekwencji.`;

    return {
      valid: false,
      hasGroundTruthReference: true,
      status: 'VALIDATION_ERROR',
      passedTargetCriterion: false,
      passedIntermediateCriterion: false,
      failureReason: mismatchMsg,
      warnings,
      iterationId,
      parameterVersion,
      timestampIso,
      modelArtifactPath: options?.modelArtifactPath ?? null,
      modelPdbHash: parsedModel.rawPdbHash,
      referencePdbHash: parsedRef.rawPdbHash,
      referenceId: options?.referenceId ?? parsedRef.headerId,
      referenceName: options?.referenceName ?? parsedRef.title,
      referenceSource: options?.referenceSource ?? 'CUSTOM_PDB',
      selectedModelChain: parsedModel.selectedChainId,
      availableModelChains: parsedModel.availableChains,
      selectedReferenceChain: parsedRef.selectedChainId,
      availableReferenceChains: parsedRef.availableChains,
      selectedModelNumber: parsedModel.selectedModelNumber,
      totalModelStructures: parsedModel.totalModelsDetected,
      selectedReferenceModelNumber: parsedRef.selectedModelNumber,
      totalReferenceStructures: parsedRef.totalModelsDetected,
      caRmsdAngstrom: null,
      independentCaRmsdAngstrom: null,
      quaternionCrossCheckCaRmsdAngstrom: null,
      numericalAgreementDeltaAngstrom: null,
      independentVerificationPassed: false,
      backboneRmsdAngstrom: null,
      allHeavyAtomRmsdAngstrom: null,
      tmScore: null,
      gdtTsPct: null,
      meanReferencePlddtOrBfactor: null,
      pairedResidueCount: seqMapping.pairedCount,
      modelResidueCount: parsedModel.residuesWithCa.length,
      referenceResidueCount: parsedRef.residuesWithCa.length,
      modelCoveragePct: seqMapping.modelCoveragePct,
      referenceCoveragePct: seqMapping.referenceCoveragePct,
      sequenceIdentityPct: seqMapping.sequenceIdentityPct,
      mismatchedResidueCount: seqMapping.mismatchedPairs.length,
      superposition: null,
      sequenceMapping: seqMapping,
      modelGeometry,
      perResidueDeviations: [],
      alignedModelCaCoords: [],
      pairedReferenceCaCoords: [],
      criteriaApplied: criteria
    };
  }

  // Extract paired Cα coordinates strictly according to unambiguous sequence mapping
  // If allowSequenceMismatches is true, only use exact or all mapped Cα pairs explicitly
  const rawModelCa: Vec3Tuple[] = [];
  const rawRefCa: Vec3Tuple[] = [];
  for (const pair of seqMapping.pairedResidues) {
    const mRes = parsedModel.residuesWithCa[pair.modelChainIndex];
    const rRes = parsedRef.residuesWithCa[pair.refChainIndex];
    if (mRes?.caCoord && rRes?.caCoord) {
      rawModelCa.push([mRes.caCoord[0], mRes.caCoord[1], mRes.caCoord[2]]);
      rawRefCa.push([rRes.caCoord[0], rRes.caCoord[1], rRes.caCoord[2]]);
    }
  }

  // Primary Kabsch SVD Superposition
  const kabsch = kabschSuperposeSvd(rawModelCa, rawRefCa);

  // Independent Secondary RMSD Re-computation (no shared intermediate arrays or cached RMSD)
  const independent = recomputeCaRmsdIndependently(
    rawModelCa,
    rawRefCa,
    kabsch.rotationMatrix,
    kabsch.translationVec
  );

  const deltaDirect = Math.abs(kabsch.caRmsdAngstrom - independent.directTransformedRmsdAngstrom);
  const deltaQuat = Math.abs(kabsch.caRmsdAngstrom - independent.quaternionEigenRmsdAngstrom);
  const numericalAgreementDelta = Math.max(deltaDirect, deltaQuat);
  const independentVerificationPassed =
    deltaDirect <= numericalTol &&
    deltaQuat <= Math.max(1e-5, numericalTol * 100) &&
    Math.abs(kabsch.determinantR - 1.0) <= 1e-6;

  if (!independentVerificationPassed) {
    return {
      valid: false,
      hasGroundTruthReference: true,
      status: 'VALIDATION_ERROR',
      passedTargetCriterion: false,
      passedIntermediateCriterion: false,
      failureReason: `Niezgodność numeryczna estymatorów RMSD: główny Kabsch SVD = ${kabsch.caRmsdAngstrom.toFixed(8)} Å, niezależny XYZ = ${independent.directTransformedRmsdAngstrom.toFixed(8)} Å, kwaternionowy = ${independent.quaternionEigenRmsdAngstrom.toFixed(8)} Å (det R = ${kabsch.determinantR.toFixed(6)}).`,
      warnings,
      iterationId,
      parameterVersion,
      timestampIso,
      modelArtifactPath: options?.modelArtifactPath ?? null,
      modelPdbHash: parsedModel.rawPdbHash,
      referencePdbHash: parsedRef.rawPdbHash,
      referenceId: options?.referenceId ?? parsedRef.headerId,
      referenceName: options?.referenceName ?? parsedRef.title,
      referenceSource: options?.referenceSource ?? 'CUSTOM_PDB',
      selectedModelChain: parsedModel.selectedChainId,
      availableModelChains: parsedModel.availableChains,
      selectedReferenceChain: parsedRef.selectedChainId,
      availableReferenceChains: parsedRef.availableChains,
      selectedModelNumber: parsedModel.selectedModelNumber,
      totalModelStructures: parsedModel.totalModelsDetected,
      selectedReferenceModelNumber: parsedRef.selectedModelNumber,
      totalReferenceStructures: parsedRef.totalModelsDetected,
      caRmsdAngstrom: kabsch.caRmsdAngstrom,
      independentCaRmsdAngstrom: independent.directTransformedRmsdAngstrom,
      quaternionCrossCheckCaRmsdAngstrom: independent.quaternionEigenRmsdAngstrom,
      numericalAgreementDeltaAngstrom: numericalAgreementDelta,
      independentVerificationPassed: false,
      backboneRmsdAngstrom: null,
      allHeavyAtomRmsdAngstrom: null,
      tmScore: null,
      gdtTsPct: null,
      meanReferencePlddtOrBfactor: null,
      pairedResidueCount: rawModelCa.length,
      modelResidueCount: parsedModel.residuesWithCa.length,
      referenceResidueCount: parsedRef.residuesWithCa.length,
      modelCoveragePct: seqMapping.modelCoveragePct,
      referenceCoveragePct: seqMapping.referenceCoveragePct,
      sequenceIdentityPct: seqMapping.sequenceIdentityPct,
      mismatchedResidueCount: seqMapping.mismatchedPairs.length,
      superposition: {
        rotationMatrix: kabsch.rotationMatrix,
        translationVec: kabsch.translationVec,
        modelCentroid: kabsch.modelCentroid,
        referenceCentroid: kabsch.referenceCentroid,
        determinantR: kabsch.determinantR,
        svdSingularValues: kabsch.svdSingularValues
      },
      sequenceMapping: seqMapping,
      modelGeometry,
      perResidueDeviations: [],
      alignedModelCaCoords: kabsch.alignedModelCaCoords,
      pairedReferenceCaCoords: kabsch.referenceCaCoords,
      criteriaApplied: criteria
    };
  }

  // Compute Backbone (N, CA, C, O), All-Atom RMSD, TM-score, GDT-TS & per-residue profile
  // strictly under the Cα Kabsch transformation x_aligned = R * x_model + t
  const R = kabsch.rotationMatrix;
  const t = kabsch.translationVec;
  const applyTransform = (pt: Vec3Tuple): Vec3Tuple => [
    R[0][0] * pt[0] + R[0][1] * pt[1] + R[0][2] * pt[2] + t[0],
    R[1][0] * pt[0] + R[1][1] * pt[1] + R[1][2] * pt[2] + t[1],
    R[2][0] * pt[0] + R[2][1] * pt[1] + R[2][2] * pt[2] + t[2]
  ];

  const Ltarget = Math.max(15, parsedRef.residuesWithCa.length);
  const d0 = Math.max(0.5, 1.24 * Math.cbrt(Ltarget - 15) - 1.8);

  let bbSqSum = 0, bbCount = 0;
  let allSqSum = 0, allCount = 0;
  let tmSum = 0;
  let gdt1 = 0, gdt2 = 0, gdt4 = 0, gdt8 = 0;
  let bFactorSum = 0;

  const perResidueDeviations: GroundTruthResidueDeviation[] = [];

  for (let idx = 0; idx < seqMapping.pairedResidues.length; idx++) {
    const pair = seqMapping.pairedResidues[idx];
    const mRes = parsedModel.residuesWithCa[pair.modelChainIndex];
    const rRes = parsedRef.residuesWithCa[pair.refChainIndex];

    const alignedCa = kabsch.alignedModelCaCoords[idx];
    const refCa = kabsch.referenceCaCoords[idx];
    const caDist = Math.hypot(
      refCa[0] - alignedCa[0],
      refCa[1] - alignedCa[1],
      refCa[2] - alignedCa[2]
    );

    tmSum += 1.0 / (1.0 + (caDist / d0) * (caDist / d0));
    if (caDist <= 1.0) gdt1++;
    if (caDist <= 2.0) gdt2++;
    if (caDist <= 4.0) gdt4++;
    if (caDist <= 8.0) gdt8++;
    bFactorSum += rRes.bFactorMean;

    let resBbSq = 0, resBbCnt = 0;
    let resScSq = 0, resScCnt = 0;
    let resAllSq = 0, resAllCnt = 0;
    let maxAtomDist = caDist;
    let maxAtomName = 'CA';

    rRes.atoms.forEach((rPt, atomName) => {
      if (atomName === 'OXT') return;
      const mPt = mRes.atoms.get(atomName);
      if (!mPt) return;
      const mAligned = applyTransform(mPt);
      const sq =
        (rPt[0] - mAligned[0]) ** 2 +
        (rPt[1] - mAligned[1]) ** 2 +
        (rPt[2] - mAligned[2]) ** 2;
      const d = Math.sqrt(sq);
      if (d >= maxAtomDist) {
        maxAtomDist = d;
        maxAtomName = atomName;
      }

      resAllSq += sq;
      resAllCnt++;
      allSqSum += sq;
      allCount++;

      if (atomName === 'N' || atomName === 'CA' || atomName === 'C' || atomName === 'O') {
        resBbSq += sq;
        resBbCnt++;
        bbSqSum += sq;
        bbCount++;
      } else {
        resScSq += sq;
        resScCnt++;
      }
    });

    perResidueDeviations.push({
      modelResSeq: mRes.resSeq,
      refResSeq: rRes.resSeq,
      modelResName: mRes.canonicalResName,
      refResName: rRes.canonicalResName,
      code1: THREE_TO_ONE_MAP[mRes.canonicalResName] || 'X',
      caDistAngstrom: parseFloat(caDist.toFixed(4)),
      backboneRmsdAngstrom: resBbCnt > 0 ? parseFloat(Math.sqrt(resBbSq / resBbCnt).toFixed(4)) : parseFloat(caDist.toFixed(4)),
      sidechainRmsdAngstrom: resScCnt > 0 ? parseFloat(Math.sqrt(resScSq / resScCnt).toFixed(4)) : 0,
      allAtomRmsdAngstrom: resAllCnt > 0 ? parseFloat(Math.sqrt(resAllSq / resAllCnt).toFixed(4)) : parseFloat(caDist.toFixed(4)),
      maxAtomDistAngstrom: parseFloat(maxAtomDist.toFixed(4)),
      maxAtomName,
      matchedAtomCount: resAllCnt,
      matchedSidechainAtomCount: resScCnt,
      refBfactorOrPlddt: parseFloat(rRes.bFactorMean.toFixed(2))
    });
  }

  const nPaired = Math.max(1, rawModelCa.length);
  const caRmsdExact = parseFloat(kabsch.caRmsdAngstrom.toFixed(6));
  const indRmsdExact = parseFloat(independent.directTransformedRmsdAngstrom.toFixed(6));
  const quatRmsdExact = parseFloat(independent.quaternionEigenRmsdAngstrom.toFixed(6));
  const bbRmsd = bbCount > 0 ? parseFloat(Math.sqrt(bbSqSum / bbCount).toFixed(4)) : caRmsdExact;
  const allRmsd = allCount > 0 ? parseFloat(Math.sqrt(allSqSum / allCount).toFixed(4)) : bbRmsd;
  const tmScore = parseFloat((tmSum / nPaired).toFixed(4));
  const gdtTsPct = parseFloat((((gdt1 + gdt2 + gdt4 + gdt8) / (4.0 * nPaired)) * 100.0).toFixed(2));
  const meanB = parseFloat((bFactorSum / nPaired).toFixed(2));

  // Evaluate pass/fail criteria strictly
  const coveragePassed =
    seqMapping.modelCoveragePct >= criteria.minModelCoveragePct &&
    seqMapping.referenceCoveragePct >= criteria.minReferenceCoveragePct;
  const seqPassed = !criteria.requireExactSequenceIdentity || seqMapping.mismatchedPairs.length === 0;
  const geomPassed = !criteria.requireValidPhysicalGeometry || modelGeometry.geometryValid;

  let status: GroundTruthIterationStatus;
  let passedTarget = false;
  let passedIntermediate = false;
  let failureReason: string | null = null;

  if (!geomPassed) {
    status = 'GEOMETRY_VIOLATION';
    failureReason = `Naruszenie geometrii fizycznej modelu: ${modelGeometry.reasons.join(' ')}`;
  } else if (!coveragePassed || !seqPassed) {
    status = 'VALIDATION_ERROR';
    failureReason = !seqPassed
      ? `Niezgodność sekwencji (${seqMapping.mismatchedPairs.length} niezgodnych reszt).`
      : `Zbyt niskie pokrycie sekwencji: model ${seqMapping.modelCoveragePct}% (wymagane ≥${criteria.minModelCoveragePct}%), referencja ${seqMapping.referenceCoveragePct}% (wymagane ≥${criteria.minReferenceCoveragePct}%).`;
  } else if (kabsch.caRmsdAngstrom < criteria.targetCaRmsdAngstrom) {
    status = 'PASSED_TARGET';
    passedTarget = true;
    passedIntermediate = true;
  } else if (
    criteria.intermediateCaRmsdAngstrom !== undefined &&
    kabsch.caRmsdAngstrom < criteria.intermediateCaRmsdAngstrom
  ) {
    status = 'PASSED_INTERMEDIATE';
    passedTarget = false;
    passedIntermediate = true;
    failureReason = `Osiągnięto próg pośredni (<${criteria.intermediateCaRmsdAngstrom} Å), lecz Cα RMSD = ${kabsch.caRmsdAngstrom.toFixed(4)} Å przekracza próg docelowy <${criteria.targetCaRmsdAngstrom} Å.`;
  } else {
    status = 'BELOW_THRESHOLD';
    failureReason = `Brak zbieżności do progu docelowego: rzeczywisty Cα RMSD = ${kabsch.caRmsdAngstrom.toFixed(4)} Å ≥ ${criteria.targetCaRmsdAngstrom.toFixed(2)} Å (${rawModelCa.length}/${parsedModel.residuesWithCa.length} sparowanych reszt).`;
  }

  return {
    valid: status !== 'VALIDATION_ERROR',
    hasGroundTruthReference: true,
    status,
    passedTargetCriterion: passedTarget,
    passedIntermediateCriterion: passedIntermediate,
    failureReason,
    warnings,
    iterationId,
    parameterVersion,
    timestampIso,
    modelArtifactPath: options?.modelArtifactPath ?? null,
    modelPdbHash: parsedModel.rawPdbHash,
    referencePdbHash: parsedRef.rawPdbHash,
    referenceId: options?.referenceId ?? parsedRef.headerId,
    referenceName: options?.referenceName ?? parsedRef.title,
    referenceSource: options?.referenceSource ?? 'CUSTOM_PDB',
    selectedModelChain: parsedModel.selectedChainId,
    availableModelChains: parsedModel.availableChains,
    selectedReferenceChain: parsedRef.selectedChainId,
    availableReferenceChains: parsedRef.availableChains,
    selectedModelNumber: parsedModel.selectedModelNumber,
    totalModelStructures: parsedModel.totalModelsDetected,
    selectedReferenceModelNumber: parsedRef.selectedModelNumber,
    totalReferenceStructures: parsedRef.totalModelsDetected,
    caRmsdAngstrom: caRmsdExact,
    independentCaRmsdAngstrom: indRmsdExact,
    quaternionCrossCheckCaRmsdAngstrom: quatRmsdExact,
    numericalAgreementDeltaAngstrom: numericalAgreementDelta,
    independentVerificationPassed: true,
    backboneRmsdAngstrom: bbRmsd,
    allHeavyAtomRmsdAngstrom: allRmsd,
    tmScore,
    gdtTsPct,
    meanReferencePlddtOrBfactor: meanB,
    pairedResidueCount: rawModelCa.length,
    modelResidueCount: parsedModel.residuesWithCa.length,
    referenceResidueCount: parsedRef.residuesWithCa.length,
    modelCoveragePct: seqMapping.modelCoveragePct,
    referenceCoveragePct: seqMapping.referenceCoveragePct,
    sequenceIdentityPct: seqMapping.sequenceIdentityPct,
    mismatchedResidueCount: seqMapping.mismatchedPairs.length,
    superposition: {
      rotationMatrix: kabsch.rotationMatrix,
      translationVec: kabsch.translationVec,
      modelCentroid: kabsch.modelCentroid,
      referenceCentroid: kabsch.referenceCentroid,
      determinantR: kabsch.determinantR,
      svdSingularValues: kabsch.svdSingularValues
    },
    sequenceMapping: seqMapping,
    modelGeometry,
    perResidueDeviations,
    alignedModelCaCoords: kabsch.alignedModelCaCoords,
    pairedReferenceCaCoords: kabsch.referenceCaCoords,
    criteriaApplied: criteria
  };
}
