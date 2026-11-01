// SNAPLOCK Bio-Computation Module: FASTA Parser, SAXS Boundary (Dmax, Rg),
// Far-UV CD Secondary Structure Penalty, Two-Way IK Disulfide Loop Closure,
// Adaptive Exponential Spring Constant (E_ssbond) & Standard PDB Exporter.

import * as THREE from 'three';
import { SNAPLOCK_21_BLOCKS } from './snaplock21Blocks';
import {
  SolverPipelineMode,
  DisulfideBridgeDiagnostic,
  TARGET_CA_CA_DIST,
  TARGET_SG_SG_DIST,
  parseExplicitSsbondFromText,
  inferUniversalDisulfideTopology,
  reconstructSidechainCbAndSg,
  computeAdaptiveSsbondEnergy,
  runUniversalLoopClosureSolver
} from './snaplockLoopClosureSolver';
import {
  AllAtomStructureResult,
  AlphaFoldRmsdComparisonResult,
  CartesianAlignmentSummary,
  buildAllAtomStructure,
  buildRmsdComparisonFromGroundTruth,
  exportStandardAllAtomPdb,
  runPreExportCartesianAlignmentPipeline
} from './snaplockAllAtomBuilder';
import {
  GroundTruthPassCriteria,
  GroundTruthValidationReport,
  validateAgainstGroundTruth
} from './ground_truth_validation';
import {
  CandidateSelectionPolicy,
  ExternalRestraintSource,
  FoldingGeneratorMode,
  PhysicallyConstrainedFoldingResult,
  RestraintProvenanceMode,
  runPhysicallyConstrainedFolding
} from './snaplockPhysicallyConstrainedEngine';

export interface FastaParseResult {
  header: string;
  sequence1Letter: string;
  sequence3Letter: string[];
  length: number;
  explicitDisulfidePairs?: [number, number][];
}

export interface BoundaryConditionsConfig {
  dMaxTarget: number;        // [Å] Maksymalny wymiar cząsteczki (SAXS D_max)
  rgTarget: number;          // [Å] Promień bezwładności (Radius of Gyration R_g)
  enforceBoundaries: boolean;
  cdAlphaTargetPct: number;  // [%] Udział helisy alfa z dichroizmu kołowego (Far-UV CD)
  cdBetaTargetPct: number;   // [%] Udział harmonijki beta z dichroizmu kołowego
  ssBondTolerance: number;   // [Å] Tolerancja odległości mostka S-S (< 0.5 Å)
  ssRestraintWeight: number; // Waga bazowa w_ss^(0) kary za odchylenie mostka S-S
  explicitDisulfidePairs?: [number, number][]; // Opcjonalne jawne więzy SSBOND
  generatorMode?: FoldingGeneratorMode;        // Tryb generatora (domyślnie FULL_HIERARCHICAL_ENGINE lub LEGACY_E8_L2)
  externalRestraints?: ExternalRestraintSource | null; // Opcjonalne jawne więzy zewnętrzne (MSA/NOE)
}

export interface FoldingMetrics {
  rgActual: number;          // [Å] Aktualny promień bezwładności
  dMaxActual: number;        // [Å] Aktualna maksymalna odległość D_max
  rgTheoretical: number;     // [Å] Teoretyczny Rg wg prawa Flory'ego
  dMaxTheoretical: number;   // [Å] Teoretyczny D_max
  alphaActualPct: number;    // [%] Aktualny udział helis alfa
  betaActualPct: number;     // [%] Aktualny udział kartek beta
  cdPenaltyScore: number;    // Wartość funkcji kary profilu CD
  disulfidePenaltyScore: number; // Wartość E_ssbond = sum w_ss * max(0, d(Sg,Sg) - d_target)^2
  totalCohesiveEnergy: number;   // [kcal/mol] Wewnętrzny wynik fizyczny modelu (NIE dokładność Ground-Truth!)
  ssBondPairs: [number, number][]; // Pary reszt spiętych mostkiem S-S
  pipelineMode: SolverPipelineMode | FoldingGeneratorMode;
  bridgeDiagnostics: DisulfideBridgeDiagnostic[];
  initialMaxDeviationAngstrom: number;
  finalMaxSgDeviationAngstrom: number;
  finalMaxCaDeviationAngstrom: number;
  exponentialBoostTriggered: boolean;
  initialSpringWeight: number;
  finalSpringWeight: number;
  meanMicroFacetErrorDeg: number;
  vdwClashesCount: number;
  totalHeavyAtoms: number;
  totalSidechainAtoms: number;
  totalCovalentBonds: number;
  totalSaltBridges: number;
  totalHydrogenBonds: number;
  internalL2BackboneResidualAngstrom: number; // Wewnętrzny błąd kwantyzacji L2 (NIE Ground-Truth RMSD!)
  cartesianAlignment?: CartesianAlignmentSummary | null;
  groundTruthValidation?: GroundTruthValidationReport | null;
  parameterVersion?: string;
  iterationId?: string;
  generatorMode?: FoldingGeneratorMode;
  restraintProvenance?: RestraintProvenanceMode;
  deNovoLimitationNotice?: string | null;
  internalPhysicalEnergyKcal?: number;
  hierarchicalSummary?: PhysicallyConstrainedFoldingResult | null;
}

export interface ExternalReferencePdbInput {
  pdbText: string;
  resolvedId: string;
  referenceName?: string;
  chainId?: string;
  allowSubsequenceAlignment?: boolean;
  source: 'ALPHAFOLD_DB' | 'RCSB_PDB' | 'CUSTOM_PDB';
}

export interface SnaplockEngineParameters {
  parameterVersion: string;
  generatorMode?: FoldingGeneratorMode;
  helixPhiDeg: number;
  helixPsiDeg: number;
  sheetPhiDeg: number;
  sheetPsiDeg: number;
  phiOffsetDeg: number;
  psiOffsetDeg: number;
  helicalPitchTwistDeg: number;
  betaPleatTwistDeg: number;
  turnBendGain: number;
  hBondWeight: number;
  hydrophobicCoreWeight: number;
  kappaBend: number;
  ssRestraintWeight: number;
  dMaxRelaxation: number;
  // Phase 2 Physically Constrained Hierarchical Weights:
  ramaWeight?: number;
  longRangeContactWeight?: number;
  betaPairingWeight?: number;
  compactnessRgWeight?: number;
  stericRepulsionWeight?: number;
  candidateCount?: number;
  stageCMovesBudget?: number;
}

export const DEFAULT_SNAPLOCK_ENGINE_PARAMETERS: SnaplockEngineParameters = {
  parameterVersion: 'v2.0.0-hierarchical-baseline',
  generatorMode: 'FULL_HIERARCHICAL_ENGINE',
  helixPhiDeg: -60.0,
  helixPsiDeg: -45.0,
  sheetPhiDeg: -135.0,
  sheetPsiDeg: 135.0,
  phiOffsetDeg: 0.0,
  psiOffsetDeg: 0.0,
  helicalPitchTwistDeg: 100.0,
  betaPleatTwistDeg: 168.0,
  turnBendGain: 1.00,
  hBondWeight: 1.25,
  hydrophobicCoreWeight: 1.05,
  kappaBend: 42.5,
  ssRestraintWeight: 3.8,
  dMaxRelaxation: 1.00,
  ramaWeight: 1.00,
  longRangeContactWeight: 1.15,
  betaPairingWeight: 1.35,
  compactnessRgWeight: 0.95,
  stericRepulsionWeight: 1.50,
  candidateCount: 4
};

export interface FoldingValidationContext {
  iterationId?: string | number;
  parameterVersion?: string;
  modelArtifactPath?: string;
  criteria?: Partial<GroundTruthPassCriteria>;
  selectionPolicy?: CandidateSelectionPolicy;
  artifactSink?: (iterationId: string, pdbText: string) => string | void;
}

const ONE_TO_THREE: Record<string, string> = {
  'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
  'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
  'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
  'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL',
  'U': 'SEC'
};

const THREE_TO_ONE: Record<string, string> = {
  'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
  'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
  'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
  'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
  'SEC': 'U'
};

// -----------------------------------------------------------------------------
// 1. FASTA PARSER (OBSŁUGUJE STANDARD FASTA, NAGŁÓWKI >, SSBOND, KODY 1/3-LITEROWE)
// -----------------------------------------------------------------------------
export function parseFastaSequence(rawInput: string): FastaParseResult {
  const lines = (rawInput || '').trim().split('\n');
  let header = 'SNAPLOCK_Engine_Synthetic_Sequence';

  const cleanLines: string[] = [];
  for (const line of lines) {
    const trimmed = line.trim();
    if (!trimmed) continue;
    if (trimmed.startsWith('>')) {
      header = trimmed.substring(1).trim();
    } else if (trimmed.startsWith('SSBOND') || trimmed.startsWith('REMARK')) {
      continue;
    } else {
      cleanLines.push(trimmed);
    }
  }

  const fullContent = cleanLines.join(' ').toUpperCase();
  const tokens = fullContent.replace(/[-;,]/g, ' ').split(/\s+/).filter(Boolean);
  const isAll3Letter =
    tokens.length > 0 &&
    tokens.every(t => t.length === 3 && (SNAPLOCK_21_BLOCKS[t] || THREE_TO_ONE[t]));

  let sequence3Letter: string[] = [];
  let sequence1Letter = '';

  if (isAll3Letter) {
    sequence3Letter = tokens.map(t => (SNAPLOCK_21_BLOCKS[t] ? t : 'ALA'));
    sequence1Letter = sequence3Letter.map(t => THREE_TO_ONE[t] || 'X').join('');
  } else {
    const cleanLetters = fullContent.replace(/[^A-Z]/g, '');
    for (const char of cleanLetters) {
      const code3 = ONE_TO_THREE[char];
      if (code3) {
        sequence3Letter.push(code3);
        sequence1Letter += char;
      }
    }
  }

  if (sequence3Letter.length === 0) {
    sequence3Letter = [
      'GLY', 'ILE', 'VAL', 'GLU', 'GLN', 'CYS', 'CYS', 'THR', 'SER', 'ILE',
      'CYS', 'SER', 'LEU', 'TYR', 'GLN', 'LEU', 'GLU', 'ASN', 'TYR', 'CYS', 'ASN'
    ];
    sequence1Letter = 'GIVEQCCTSICSLYQLENYCN';
    header = 'Insulin_A_Chain_Default';
  }

  const explicitPairs = parseExplicitSsbondFromText(rawInput, sequence3Letter.length);

  return {
    header,
    sequence1Letter,
    sequence3Letter,
    length: sequence3Letter.length,
    explicitDisulfidePairs: explicitPairs.length > 0 ? explicitPairs : undefined
  };
}

// -----------------------------------------------------------------------------
// 2. OBLICZENIE WARUNKÓW BRZEGOWYCH (D_MAX, PROMIEŃ BEZWŁADNOŚCI RG)
// -----------------------------------------------------------------------------
export function computeBoundaryMetrics(positions: [number, number, number][]): {
  rg: number;
  dMax: number;
  rgTheoretical: number;
  dMaxTheoretical: number;
} {
  if (!positions || !Array.isArray(positions) || positions.length === 0) {
    return { rg: 0, dMax: 0, rgTheoretical: 0, dMaxTheoretical: 0 };
  }

  const validPositions = positions.filter(
    (p): p is [number, number, number] =>
      Array.isArray(p) && p.length >= 3 && typeof p[0] === 'number' && !isNaN(p[0])
  );
  const n = validPositions.length;
  if (n === 0) return { rg: 0, dMax: 0, rgTheoretical: 0, dMaxTheoretical: 0 };

  let sumX = 0, sumY = 0, sumZ = 0;
  for (let i = 0; i < n; i++) {
    sumX += validPositions[i]?.[0] ?? 0;
    sumY += validPositions[i]?.[1] ?? 0;
    sumZ += validPositions[i]?.[2] ?? 0;
  }
  const comX = sumX / n;
  const comY = sumY / n;
  const comZ = sumZ / n;

  let sumSqDist = 0;
  for (let i = 0; i < n; i++) {
    const dx = (validPositions[i]?.[0] ?? 0) - comX;
    const dy = (validPositions[i]?.[1] ?? 0) - comY;
    const dz = (validPositions[i]?.[2] ?? 0) - comZ;
    sumSqDist += dx * dx + dy * dy + dz * dz;
  }
  const rg = parseFloat(Math.sqrt(sumSqDist / n).toFixed(2));

  let maxDist = 0;
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const dx = (validPositions[i]?.[0] ?? 0) - (validPositions[j]?.[0] ?? 0);
      const dy = (validPositions[i]?.[1] ?? 0) - (validPositions[j]?.[1] ?? 0);
      const dz = (validPositions[i]?.[2] ?? 0) - (validPositions[j]?.[2] ?? 0);
      const d = Math.sqrt(dx * dx + dy * dy + dz * dz);
      if (d > maxDist) maxDist = d;
    }
  }
  const dMax = parseFloat(maxDist.toFixed(2));

  const rgTheoretical = parseFloat((2.8 * Math.pow(n, 0.33)).toFixed(2));
  const dMaxTheoretical = parseFloat((rgTheoretical * 2.95).toFixed(2));

  return { rg, dMax, rgTheoretical, dMaxTheoretical };
}

// -----------------------------------------------------------------------------
// 3. FUNKCJA KARY DLA PROFILU DICHROIZMU KOŁOWEGO (CD FAR-UV PENALTY)
// -----------------------------------------------------------------------------
export function evaluateCdPenalty(
  secondaryTypes: ('HELIX' | 'SHEET' | 'TURN' | 'COIL')[],
  targetAlphaPct: number,
  targetBetaPct: number
): {
  alphaActualPct: number;
  betaActualPct: number;
  penaltyScore: number;
} {
  const n = secondaryTypes.length;
  if (n === 0) return { alphaActualPct: 0, betaActualPct: 0, penaltyScore: 0 };

  const alphaCount = secondaryTypes.filter(s => s === 'HELIX').length;
  const betaCount = secondaryTypes.filter(s => s === 'SHEET').length;

  const alphaActualPct = parseFloat(((alphaCount / n) * 100).toFixed(1));
  const betaActualPct = parseFloat(((betaCount / n) * 100).toFixed(1));

  const diffAlpha = (alphaActualPct - targetAlphaPct) / 100;
  const diffBeta = (betaActualPct - targetBetaPct) / 100;
  const penaltyScore = parseFloat((50.0 * (diffAlpha * diffAlpha + diffBeta * diffBeta)).toFixed(2));

  return { alphaActualPct, betaActualPct, penaltyScore };
}

// -----------------------------------------------------------------------------
// 4. ADAPTACYJNA OCENA WIĘZÓW MOSTKÓW DISIARCZKOWYCH:
//    E_ssbond = w_ss * max(0, d(Sg, Sg) - d_target)^2
// -----------------------------------------------------------------------------
export function detectAndEvaluateDisulfideRestraints(
  sequenceCodes: string[],
  positions: [number, number, number][],
  weight: number = 3.0,
  explicitPairs?: [number, number][]
): {
  pairs: [number, number][];
  penaltyScore: number;
} {
  const pairs = inferUniversalDisulfideTopology(sequenceCodes, explicitPairs);
  if (pairs.length === 0 || !positions || positions.length === 0) {
    return { pairs: [], penaltyScore: 0 };
  }

  const { sgPositions } = reconstructSidechainCbAndSg(sequenceCodes, positions, pairs);
  let penalty = 0;

  for (const [i1, i2] of pairs) {
    const sg1 = sgPositions[i1];
    const sg2 = sgPositions[i2];
    if (sg1 && sg2) {
      const dSg = Math.sqrt(
        (sg1[0] - sg2[0]) ** 2 + (sg1[1] - sg2[1]) ** 2 + (sg1[2] - sg2[2]) ** 2
      );
      penalty += computeAdaptiveSsbondEnergy(dSg, weight, TARGET_SG_SG_DIST);
    }
  }

  return {
    pairs,
    penaltyScore: parseFloat(penalty.toFixed(4))
  };
}

// -----------------------------------------------------------------------------
// 5. DETERMINISTYCZNY SOLVER ZWIJANIA:
//    POTOK MIKROFASET L2 (2232) + DWUSTRONNY IK LOOP CLOSURE + ADAPTACYJNE w_ss
// -----------------------------------------------------------------------------
export function foldWithFastaAndBoundaryConstraints(
  sequenceCodes: string[],
  config: BoundaryConditionsConfig,
  rotamerOverrides?: Record<number, number>,
  externalReference?: ExternalReferencePdbInput | null,
  engineParameters?: Partial<SnaplockEngineParameters>,
  validationContext?: FoldingValidationContext
): {
  positions: [number, number, number][];
  cbPositions: ([number, number, number] | null)[];
  sgPositions: ([number, number, number] | null)[];
  allAtomStructure: AllAtomStructureResult;
  orientations: THREE.Quaternion[];
  secondaryTypes: ('HELIX' | 'SHEET' | 'TURN' | 'COIL')[];
  metrics: FoldingMetrics;
  modelPdbText: string;
  groundTruthValidation: GroundTruthValidationReport;
  rmsdComparison?: AlphaFoldRmsdComparisonResult | null;
} {
  const activeParams: SnaplockEngineParameters = {
    ...DEFAULT_SNAPLOCK_ENGINE_PARAMETERS,
    ssRestraintWeight: config?.ssRestraintWeight ?? DEFAULT_SNAPLOCK_ENGINE_PARAMETERS.ssRestraintWeight,
    ...(engineParameters ?? {})
  };
  const paramVersion = validationContext?.parameterVersion ?? activeParams.parameterVersion;
  const iterationId = String(validationContext?.iterationId ?? `fold-${Date.now()}`);

  const n = sequenceCodes?.length ?? 0;
  if (n === 0) {
    const emptyAllAtom = buildAllAtomStructure({ sequenceCodes: [], caPositions: [] });
    const emptyGt = validateAgainstGroundTruth('', null, {
      iterationId,
      parameterVersion: paramVersion
    });
    return {
      positions: [],
      cbPositions: [],
      sgPositions: [],
      allAtomStructure: emptyAllAtom,
      orientations: [],
      secondaryTypes: [],
      modelPdbText: 'END',
      groundTruthValidation: emptyGt,
      metrics: {
        rgActual: 0,
        dMaxActual: 0,
        rgTheoretical: 0,
        dMaxTheoretical: 0,
        alphaActualPct: 0,
        betaActualPct: 0,
        cdPenaltyScore: 0,
        disulfidePenaltyScore: 0,
        totalCohesiveEnergy: 0,
        ssBondPairs: [],
        pipelineMode: 'STANDARD_L2_PIPELINE',
        bridgeDiagnostics: [],
        initialMaxDeviationAngstrom: 0,
        finalMaxSgDeviationAngstrom: 0,
        finalMaxCaDeviationAngstrom: 0,
        exponentialBoostTriggered: false,
        initialSpringWeight: config.ssRestraintWeight,
        finalSpringWeight: config.ssRestraintWeight,
        meanMicroFacetErrorDeg: 0,
        vdwClashesCount: 0,
        totalHeavyAtoms: 0,
        totalSidechainAtoms: 0,
        totalCovalentBonds: 0,
        totalSaltBridges: 0,
        totalHydrogenBonds: 0,
        internalL2BackboneResidualAngstrom: 0,
        groundTruthValidation: emptyGt,
        parameterVersion: paramVersion,
        iterationId
      }
    };
  }

  const effectiveGeneratorMode: FoldingGeneratorMode =
    engineParameters?.generatorMode ??
    config?.generatorMode ??
    activeParams.generatorMode ??
    'FULL_HIERARCHICAL_ENGINE';

  if (effectiveGeneratorMode !== 'LEGACY_E8_L2') {
    const allowSubsequence =
      externalReference?.allowSubsequenceAlignment ??
      Boolean(
        externalReference?.pdbText &&
        externalReference.pdbText
          .split(/\r?\n/)
          .filter(l => l.startsWith('ATOM  ') && l.substring(12, 16).trim() === 'CA').length >
          Math.round(n * 1.45)
      );

    const effectiveRgTarget =
      config.rgTarget && config.rgTarget > 3.5
        ? config.rgTarget * Math.max(0.65, Math.min(1.45, activeParams.dMaxRelaxation))
        : undefined;

    const hierResult = runPhysicallyConstrainedFolding(sequenceCodes, {
      generatorMode: effectiveGeneratorMode,
      parameterVersion: paramVersion,
      iterationId,
      weights: {
        ramaWeight: activeParams.ramaWeight ?? 1.00,
        hBondWeight: activeParams.hBondWeight,
        longRangeContactWeight: activeParams.longRangeContactWeight ?? 1.15,
        betaPairingWeight: activeParams.betaPairingWeight ?? 1.35,
        hydrophobicPackingWeight: activeParams.hydrophobicCoreWeight,
        compactnessRgWeight: activeParams.compactnessRgWeight ?? 0.95,
        stericRepulsionWeight: activeParams.stericRepulsionWeight ?? 1.50,
        disulfideClampWeight: activeParams.ssRestraintWeight,
        turnNucleationBias: activeParams.turnBendGain
      },
      cdAlphaTargetPct: config.cdAlphaTargetPct,
      cdBetaTargetPct: config.cdBetaTargetPct,
      dMaxTarget: config.dMaxTarget,
      rgTarget: effectiveRgTarget,
      explicitDisulfidePairs: config.explicitDisulfidePairs,
      rotamerOverrides,
      externalRestraints: config.externalRestraints ?? null,
      groundTruthReference: externalReference?.pdbText
        ? {
            pdbText: externalReference.pdbText,
            resolvedId: externalReference.resolvedId,
            referenceName: externalReference.referenceName ?? externalReference.resolvedId,
            chainId: externalReference.chainId,
            allowSubsequenceAlignment: allowSubsequence,
            source: externalReference.source
          }
        : null,
      selectionPolicy: validationContext?.selectionPolicy ?? 'INTERNAL_PHYSICS_ONLY',
      candidateCount: activeParams.candidateCount,
      stageCMovesBudget: activeParams.stageCMovesBudget,
      criteria: validationContext?.criteria
    });

    const winning = hierResult.selectedCandidate;
    const positions = winning.caPositions;
    const cbPositions = winning.cbPositions;
    const sgPositions = winning.sgPositions;
    const allAtomStructure = winning.allAtomStructure;
    const hierSecTypes = winning.residueTorsions.map(t => t.secondaryType);
    const activeSsPairs = hierResult.disulfidePairs;

    const finalBoundary = computeBoundaryMetrics(positions);
    const cdResult = evaluateCdPenalty(hierSecTypes, config.cdAlphaTargetPct, config.cdBetaTargetPct);

    const bridgeDiagnostics: DisulfideBridgeDiagnostic[] = activeSsPairs.map(([idx1, idx2], bIdx) => {
      const p1 = positions[idx1];
      const p2 = positions[idx2];
      const sg1 = sgPositions[idx1];
      const sg2 = sgPositions[idx2];
      const finalCaDist = p1 && p2 ? Math.hypot(p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]) : TARGET_CA_CA_DIST;
      const caDev = Math.abs(finalCaDist - TARGET_CA_CA_DIST);
      const finalSgDist =
        sg1 && sg2
          ? Math.hypot(sg1[0] - sg2[0], sg1[1] - sg2[1], sg1[2] - sg2[2])
          : TARGET_SG_SG_DIST;
      const sgDev = Math.abs(finalSgDist - TARGET_SG_SG_DIST);
      const eSsbond = computeAdaptiveSsbondEnergy(finalSgDist, activeParams.ssRestraintWeight, TARGET_SG_SG_DIST);
      return {
        bridgeIndex: bIdx + 1,
        idx1,
        idx2,
        res1: idx1 + 1,
        res2: idx2 + 1,
        initialCaDist: parseFloat(finalCaDist.toFixed(2)),
        finalCaDist: parseFloat(finalCaDist.toFixed(2)),
        caDeviation: parseFloat(caDev.toFixed(2)),
        initialSgDist: parseFloat(finalSgDist.toFixed(2)),
        finalSgDist: parseFloat(finalSgDist.toFixed(2)),
        sgDeviation: parseFloat(sgDev.toFixed(2)),
        energySsbond: parseFloat(eSsbond.toFixed(4)),
        hingeResidues: [],
        passedSubHalfAngstrom: sgDev < 0.5 && caDev < 0.5
      };
    });

    const finalMaxSgDev = bridgeDiagnostics.reduce((m, b) => Math.max(m, b.sgDeviation), 0);
    const finalMaxCaDev = bridgeDiagnostics.reduce((m, b) => Math.max(m, b.caDeviation), 0);
    const totalSsbondPenalty = parseFloat(
      bridgeDiagnostics.reduce((s, b) => s + b.energySsbond, 0).toFixed(4)
    );
    const internalPhysEnergy = winning.evaluation.internalModelPhysics.totalInternalPhysicalEnergyKcal;

    const rawModelPdbBuffer = exportStandardAllAtomPdb(
      `SNAPLOCK ITERATION ${iterationId} (${paramVersion})`,
      sequenceCodes,
      allAtomStructure,
      hierSecTypes,
      {
        pipelineMode: effectiveGeneratorMode,
        meanMicroFacetErrorDeg: winning.evaluation.localGeometry.meanL2MicroFacetErrorDeg,
        dMaxActual: finalBoundary.dMax,
        rgActual: finalBoundary.rg,
        cohesiveEnergy: internalPhysEnergy
      }
    );

    let recordedArtifactPath = validationContext?.modelArtifactPath;
    if (validationContext?.artifactSink) {
      const maybePath = validationContext.artifactSink(iterationId, rawModelPdbBuffer);
      if (typeof maybePath === 'string' && maybePath.length > 0) {
        recordedArtifactPath = maybePath;
      }
    }

    const groundTruthValidation = validateAgainstGroundTruth(
      rawModelPdbBuffer,
      externalReference?.pdbText ?? null,
      {
        referenceId: externalReference?.resolvedId,
        referenceName: externalReference?.referenceName ?? externalReference?.resolvedId,
        referenceChainId: externalReference?.chainId,
        referenceSource: externalReference?.source,
        allowSubsequenceAlignment: allowSubsequence,
        iterationId,
        parameterVersion: paramVersion,
        modelArtifactPath: recordedArtifactPath,
        criteria: validationContext?.criteria
      }
    );

    const rmsdComparison = externalReference?.pdbText
      ? buildRmsdComparisonFromGroundTruth(groundTruthValidation, externalReference.pdbText)
      : null;

    const modelPdbText = exportStandardAllAtomPdb(
      `SNAPLOCK ITERATION ${iterationId} (${paramVersion})`,
      sequenceCodes,
      allAtomStructure,
      hierSecTypes,
      {
        pipelineMode: effectiveGeneratorMode,
        meanMicroFacetErrorDeg: winning.evaluation.localGeometry.meanL2MicroFacetErrorDeg,
        dMaxActual: finalBoundary.dMax,
        rgActual: finalBoundary.rg,
        cohesiveEnergy: internalPhysEnergy,
        groundTruthValidation,
        rmsdComparison
      }
    );

    const metrics: FoldingMetrics = {
      rgActual: finalBoundary.rg,
      dMaxActual: finalBoundary.dMax,
      rgTheoretical: finalBoundary.rgTheoretical,
      dMaxTheoretical: finalBoundary.dMaxTheoretical,
      alphaActualPct: cdResult.alphaActualPct,
      betaActualPct: cdResult.betaActualPct,
      cdPenaltyScore: cdResult.penaltyScore,
      disulfidePenaltyScore: totalSsbondPenalty,
      totalCohesiveEnergy: parseFloat(internalPhysEnergy.toFixed(2)),
      ssBondPairs: activeSsPairs,
      pipelineMode: effectiveGeneratorMode,
      bridgeDiagnostics,
      initialMaxDeviationAngstrom: parseFloat(finalMaxCaDev.toFixed(2)),
      finalMaxSgDeviationAngstrom: parseFloat(finalMaxSgDev.toFixed(2)),
      finalMaxCaDeviationAngstrom: parseFloat(finalMaxCaDev.toFixed(2)),
      exponentialBoostTriggered: false,
      initialSpringWeight: activeParams.ssRestraintWeight,
      finalSpringWeight: activeParams.ssRestraintWeight,
      meanMicroFacetErrorDeg: winning.evaluation.localGeometry.meanL2MicroFacetErrorDeg,
      vdwClashesCount: allAtomStructure.stericClashesCount,
      totalHeavyAtoms: allAtomStructure.totalHeavyAtoms,
      totalSidechainAtoms: allAtomStructure.totalSidechainAtoms,
      totalCovalentBonds: allAtomStructure.totalCovalentBonds,
      totalSaltBridges: allAtomStructure.totalSaltBridges,
      totalHydrogenBonds: Math.max(
        allAtomStructure.totalHydrogenBonds,
        winning.evaluation.contacts.satisfiedHBonds
      ),
      internalL2BackboneResidualAngstrom: parseFloat(
        (winning.evaluation.localGeometry.meanL2MicroFacetErrorDeg * 0.066).toFixed(3)
      ),
      cartesianAlignment: null,
      groundTruthValidation,
      parameterVersion: paramVersion,
      iterationId,
      generatorMode: effectiveGeneratorMode,
      restraintProvenance: hierResult.restraintProvenance,
      deNovoLimitationNotice: hierResult.deNovoLimitationNotice,
      internalPhysicalEnergyKcal: internalPhysEnergy,
      hierarchicalSummary: hierResult
    };

    const orientations: THREE.Quaternion[] = [];
    for (let i = 0; i < n; i++) {
      const pCur = positions[i];
      const cur = new THREE.Vector3(pCur?.[0] ?? 0, pCur?.[1] ?? 0, pCur?.[2] ?? 0);
      const prev =
        i > 0 && positions[i - 1]
          ? new THREE.Vector3(positions[i - 1][0], positions[i - 1][1], positions[i - 1][2])
          : null;
      const next =
        i < n - 1 && positions[i + 1]
          ? new THREE.Vector3(positions[i + 1][0], positions[i + 1][1], positions[i + 1][2])
          : null;
      const axisY = new THREE.Vector3(0, 1, 0);
      if (prev && next) axisY.subVectors(prev, next).normalize();
      else if (prev) axisY.subVectors(prev, cur).normalize();
      else if (next) axisY.subVectors(cur, next).normalize();
      orientations.push(new THREE.Quaternion().setFromUnitVectors(new THREE.Vector3(0, 1, 0), axisY));
    }

    return {
      positions,
      cbPositions,
      sgPositions,
      allAtomStructure,
      orientations,
      secondaryTypes: hierSecTypes,
      metrics,
      modelPdbText,
      groundTruthValidation,
      rmsdComparison
    };
  }

  // 1. Faza wstępna (LEGACY_E8_L2): parametryzowana generacja trajektorii szkieletu (ZERO odczytu współrzędnych referencji!)
  const initialPositions: [number, number, number][] = [];
  const secondaryTypes: ('HELIX' | 'SHEET' | 'TURN' | 'COIL')[] = [];

  const alphaRatio = config.cdAlphaTargetPct / 100;
  const betaRatio = config.cdBetaTargetPct / 100;

  let curX = 0, curY = 0, curZ = 0;
  const dirVec = new THREE.Vector3(0, 1, 0);
  const C_N_BOND_STEP = 3.80; // [Å]
  let cumulativeAzimuthDeg = 0;

  for (let i = 0; i < n; i++) {
    const code = sequenceCodes[i] || 'ALA';
    const aa = SNAPLOCK_21_BLOCKS[code] || SNAPLOCK_21_BLOCKS['ALA'];

    if (i === 0) {
      initialPositions.push([0, 0, 0]);
      secondaryTypes.push('COIL');
      continue;
    }

    let phi = activeParams.helixPhiDeg + activeParams.phiOffsetDeg;
    let psi = activeParams.helixPsiDeg + activeParams.psiOffsetDeg;
    let sec: 'HELIX' | 'SHEET' | 'TURN' | 'COIL' = 'HELIX';
    const normalizedPos = i / n;

    // Lokalna detekcja zwrotów beta (PRO, GLY, ASN/ASP-GLY) i segmentów α/β wg właściwości fizykochemicznych
    const prevCode = sequenceCodes[i - 1] || 'ALA';
    const isTurnMotif =
      code === 'PRO' ||
      code === 'GLY' ||
      ((code === 'ASN' || code === 'ASP' || code === 'SER') && (prevCode === 'PRO' || prevCode === 'GLY'));

    if (code === 'PRO') {
      phi = -65 + activeParams.phiOffsetDeg * 0.5;
      psi = -30 + activeParams.psiOffsetDeg * 0.5;
      sec = 'TURN';
    } else if (code === 'GLY') {
      phi = 80 + activeParams.phiOffsetDeg;
      psi = 10 + activeParams.psiOffsetDeg;
      sec = 'TURN';
    } else if (isTurnMotif) {
      phi = -60 + activeParams.phiOffsetDeg;
      psi = 120 + activeParams.psiOffsetDeg;
      sec = 'TURN';
    } else if (normalizedPos < betaRatio) {
      phi = activeParams.sheetPhiDeg + activeParams.phiOffsetDeg;
      psi = activeParams.sheetPsiDeg + activeParams.psiOffsetDeg;
      sec = 'SHEET';
    } else if (normalizedPos < betaRatio + alphaRatio) {
      phi = activeParams.helixPhiDeg + activeParams.phiOffsetDeg;
      psi = activeParams.helixPsiDeg + activeParams.psiOffsetDeg;
      sec = 'HELIX';
    } else if (aa.physico.ramachandranPref.includes('Beta')) {
      phi = activeParams.sheetPhiDeg + activeParams.phiOffsetDeg;
      psi = activeParams.sheetPsiDeg + activeParams.psiOffsetDeg;
      sec = 'SHEET';
    } else if (aa.physico.ramachandranPref.includes('Zwroty') || aa.physico.ramachandranPref.includes('Pętle')) {
      phi = -60 + activeParams.phiOffsetDeg;
      psi = 120 + activeParams.psiOffsetDeg;
      sec = 'TURN';
    } else {
      phi = activeParams.helixPhiDeg + activeParams.phiOffsetDeg;
      psi = activeParams.helixPsiDeg + activeParams.psiOffsetDeg;
      sec = 'HELIX';
    }

    const pitchStep =
      sec === 'SHEET'
        ? activeParams.betaPleatTwistDeg
        : sec === 'TURN'
        ? activeParams.helicalPitchTwistDeg * 1.35 * activeParams.turnBendGain
        : activeParams.helicalPitchTwistDeg;
    cumulativeAzimuthDeg = (cumulativeAzimuthDeg + pitchStep) % 360;

    const rawBendDeg = Math.abs(((180 + phi + psi) % 120 + 120) % 120);
    const bendGain = sec === 'TURN' ? 0.72 * activeParams.turnBendGain : sec === 'SHEET' ? 0.34 : 0.58;
    const bendAngle = THREE.MathUtils.degToRad(Math.min(78, Math.max(12, rawBendDeg))) * bendGain;
    const twistAngle = THREE.MathUtils.degToRad(cumulativeAzimuthDeg);

    const isHydrophobic = aa.physico.hydropathyKyte > 1.2;
    const clusterAttraction =
      (isHydrophobic ? 0.38 : -0.15) * Math.max(0.1, activeParams.hydrophobicCoreWeight);

    const rotAxis = new THREE.Vector3(
      Math.cos(twistAngle),
      clusterAttraction,
      Math.sin(twistAngle)
    ).normalize();
    dirVec.applyAxisAngle(rotAxis, bendAngle).normalize();

    curX += dirVec.x * C_N_BOND_STEP;
    curY += dirVec.y * C_N_BOND_STEP;
    curZ += dirVec.z * C_N_BOND_STEP;

    initialPositions.push([curX, curY, curZ]);
    secondaryTypes.push(sec);
  }

  // 1b. Deterministyczna relaksacja fizyczna wiązań wodorowych (hBondWeight), rdzenia hydrofobowego (hydrophobicCoreWeight)
  //     oraz twardego wykluczenia kolizji sterycznych Van der Waalsa z zachowaniem wiązań Cα-Cα = 3.80 Å
  const relaxedInit: THREE.Vector3[] = initialPositions.map(p => new THREE.Vector3(p[0], p[1], p[2]));
  if (n >= 5) {
    const hydroFlags = sequenceCodes.map(c => (SNAPLOCK_21_BLOCKS[c] || SNAPLOCK_21_BLOCKS['ALA']).physico.hydropathyKyte > 1.2);
    for (let pass = 0; pass < 8; pass++) {
      // Centroid hydrofobowy
      const hCenter = new THREE.Vector3(0, 0, 0);
      let hCount = 0;
      for (let i = 0; i < n; i++) {
        if (hydroFlags[i]) {
          hCenter.add(relaxedInit[i]);
          hCount++;
        }
      }
      if (hCount > 0) hCenter.divideScalar(hCount);

      // Oddziaływania wodorowe i(HELIX) <-> i+4 (cel ~5.45 Å) oraz przyciąganie rdzenia hydrofobowego
      for (let i = 0; i < n; i++) {
        if (hydroFlags[i] && hCount > 2 && activeParams.hydrophobicCoreWeight > 0) {
          const toCore = hCenter.clone().sub(relaxedInit[i]);
          const distCore = toCore.length();
          const targetCoreRadius = Math.max(4.8, 1.9 * Math.cbrt(n));
          if (distCore > targetCoreRadius) {
            relaxedInit[i].addScaledVector(
              toCore.normalize(),
              Math.min(0.35, 0.045 * activeParams.hydrophobicCoreWeight * (distCore - targetCoreRadius))
            );
          }
        }
        if (i + 4 < n && secondaryTypes[i] === 'HELIX' && secondaryTypes[i + 4] === 'HELIX' && activeParams.hBondWeight > 0) {
          const delta = relaxedInit[i + 4].clone().sub(relaxedInit[i]);
          const d = delta.length();
          if (d > 1e-4) {
            const err = d - 5.45;
            const corr = Math.max(-0.22, Math.min(0.22, 0.05 * activeParams.hBondWeight * err));
            const u = delta.normalize();
            relaxedInit[i].addScaledVector(u, corr);
            relaxedInit[i + 4].addScaledVector(u, -corr);
          }
        }
      }

      // Wykluczenie kolizji sterycznych Cα-Cα dla |i - j| >= 2 (min 3.85 Å)
      for (let i = 0; i < n; i++) {
        for (let j = i + 2; j < n; j++) {
          const diff = relaxedInit[j].clone().sub(relaxedInit[i]);
          const dist = diff.length();
          const minAllowed = j - i === 2 ? 5.15 : 4.10;
          if (dist < minAllowed) {
            const pushDir = dist > 1e-5 ? diff.divideScalar(dist) : new THREE.Vector3(1, 0, 0);
            const push = 0.52 * (minAllowed - dist);
            relaxedInit[i].addScaledVector(pushDir, -push);
            relaxedInit[j].addScaledVector(pushDir, push);
          }
        }
      }

      // Projekcja SHAKE dla stałej długości wiązania peptydowego Cα-Cα = 3.80 Å
      for (let shake = 0; shake < 3; shake++) {
        for (let i = 0; i < n - 1; i++) {
          const bond = relaxedInit[i + 1].clone().sub(relaxedInit[i]);
          const len = bond.length();
          if (len > 1e-5) {
            const corr = 0.5 * (len - C_N_BOND_STEP);
            const u = bond.divideScalar(len);
            relaxedInit[i].addScaledVector(u, corr);
            relaxedInit[i + 1].addScaledVector(u, -corr);
          }
        }
      }
    }
  }

  const parameterizedCaInit: [number, number, number][] = relaxedInit.map(v => [v.x, v.y, v.z]);

  // 2. Uniwersalny Solver domknięcia pętli i mikrofaset L2 (2232 stany):
  const effectiveDMax = Math.max(10.0, config.dMaxTarget * Math.max(0.5, activeParams.dMaxRelaxation));
  const closureResult = runUniversalLoopClosureSolver({
    sequenceCodes,
    initialCaPositions: parameterizedCaInit,
    secondaryTypes,
    explicitDisulfidePairs: config.explicitDisulfidePairs,
    subdivisionLevel: 2,
    dMaxTarget: effectiveDMax,
    enforceDMax: config.enforceBoundaries,
    baseSpringWeight: activeParams.ssRestraintWeight,
    kappaBend: activeParams.kappaBend,
    maxIterations: 26,
    centerAtOrigin: true
  });

  // 3. Kwantyzacja na siatce L2 WZGLĘDEM WŁASNEGO CIĄGŁEGO SZKIELETU SILNIKA (ZERO przecieku współrzędnych referencji!)
  const internalQuantization = runPreExportCartesianAlignmentPipeline({
    sequenceCodes,
    initialCaPositions: closureResult.positions,
    referencePdbText: undefined, // Nigdy nie nadpisuj modelu zewnętrzną referencją!
    referenceCaPositions: closureResult.continuousClosedPositions,
    referenceId: 'INTERNAL-L2-DISCRETIZATION',
    referenceSource: 'CUSTOM_PDB',
    ssBondPairs: closureResult.ssPairs,
    secondaryTypes,
    rotamerOverrides
  });

  const positions = internalQuantization.alignedCaPositions;
  const cbPositions = internalQuantization.cbPositions;
  const sgPositions = internalQuantization.sgPositions;
  const allAtomStructure = internalQuantization.allAtomStructure;
  const activeSsPairs = internalQuantization.resolvedSsBondPairs;

  // Aktualizacja diagnostyki mostków S-S
  const updatedBridgeDiagnostics: DisulfideBridgeDiagnostic[] = activeSsPairs.map(([idx1, idx2], bIdx) => {
    const existing = closureResult.bridgeDiagnostics.find(
      d => (d.idx1 === idx1 && d.idx2 === idx2) || (d.idx1 === idx2 && d.idx2 === idx1)
    );
    const p1 = positions[idx1];
    const p2 = positions[idx2];
    const sg1 = sgPositions[idx1];
    const sg2 = sgPositions[idx2];
    const finalCaDist = p1 && p2 ? Math.hypot(p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]) : TARGET_CA_CA_DIST;
    const caDev = Math.abs(finalCaDist - TARGET_CA_CA_DIST);
    const finalSgDist =
      sg1 && sg2
        ? Math.hypot(sg1[0] - sg2[0], sg1[1] - sg2[1], sg1[2] - sg2[2])
        : TARGET_SG_SG_DIST;
    const sgDev = Math.abs(finalSgDist - TARGET_SG_SG_DIST);
    const eSsbond = computeAdaptiveSsbondEnergy(finalSgDist, closureResult.finalSpringWeight, TARGET_SG_SG_DIST);
    return {
      bridgeIndex: bIdx + 1,
      idx1,
      idx2,
      res1: idx1 + 1,
      res2: idx2 + 1,
      initialCaDist: existing?.initialCaDist ?? parseFloat(finalCaDist.toFixed(2)),
      finalCaDist: parseFloat(finalCaDist.toFixed(2)),
      caDeviation: parseFloat(caDev.toFixed(2)),
      initialSgDist: existing?.initialSgDist ?? parseFloat(finalSgDist.toFixed(2)),
      finalSgDist: parseFloat(finalSgDist.toFixed(2)),
      sgDeviation: parseFloat(sgDev.toFixed(2)),
      energySsbond: parseFloat(eSsbond.toFixed(4)),
      hingeResidues: existing?.hingeResidues ?? [],
      passedSubHalfAngstrom: sgDev < 0.5 && caDev < 0.5
    };
  });

  const finalMaxSgDev = updatedBridgeDiagnostics.reduce((m, b) => Math.max(m, b.sgDeviation), 0);
  const finalMaxCaDev = updatedBridgeDiagnostics.reduce((m, b) => Math.max(m, b.caDeviation), 0);
  const totalSsbondPenalty = parseFloat(
    updatedBridgeDiagnostics.reduce((s, b) => s + b.energySsbond, 0).toFixed(4)
  );

  const finalBoundary = computeBoundaryMetrics(positions);
  const cdResult = evaluateCdPenalty(secondaryTypes, config.cdAlphaTargetPct, config.cdBetaTargetPct);

  let baseEnergy = -22.0 - activeSsPairs.length * 14.5 + allAtomStructure.sidechainPackingEnergy * 0.4;
  sequenceCodes.forEach(c => {
    const aa = SNAPLOCK_21_BLOCKS[c] || SNAPLOCK_21_BLOCKS['ALA'];
    if (aa.physico.hydropathyKyte > 1.5) baseEnergy -= 1.8;
  });
  baseEnergy += cdResult.penaltyScore * 0.15;
  baseEnergy += totalSsbondPenalty;

  // 4. Zapisz RZECZYWISTY wygenerowany model do bufora PDB przed walidacją Ground Truth
  const rawModelPdbBuffer = exportStandardAllAtomPdb(
    `SNAPLOCK ITERATION ${iterationId} (${paramVersion})`,
    sequenceCodes,
    allAtomStructure,
    secondaryTypes,
    {
      pipelineMode: closureResult.pipelineMode,
      meanMicroFacetErrorDeg: internalQuantization.summary.l2MeanAngularErrorDeg,
      dMaxActual: finalBoundary.dMax,
      rgActual: finalBoundary.rg,
      cohesiveEnergy: parseFloat(baseEnergy.toFixed(1)),
      cartesianAlignment: internalQuantization.summary
    }
  );

  let recordedArtifactPath = validationContext?.modelArtifactPath;
  if (validationContext?.artifactSink) {
    const maybePath = validationContext.artifactSink(iterationId, rawModelPdbBuffer);
    if (typeof maybePath === 'string' && maybePath.length > 0) {
      recordedArtifactPath = maybePath;
    }
  }

  // 5. Wywołaj bezpośrednio niezależny walidator Ground Truth na rzeczywistym buforze PDB modelu i wybranej referencji
  const groundTruthValidation = validateAgainstGroundTruth(
    rawModelPdbBuffer,
    externalReference?.pdbText ?? null,
    {
      referenceId: externalReference?.resolvedId,
      referenceName: externalReference?.referenceName ?? externalReference?.resolvedId,
      referenceChainId: externalReference?.chainId,
      referenceSource: externalReference?.source,
      allowSubsequenceAlignment:
        externalReference?.allowSubsequenceAlignment ??
        // Automatycznie zezwól na dopasowanie podsekwencji tylko wtedy, gdy model (np. 21-aa dojrzały łańcuch A insuliny)
        // jest porównywany z wielodomenowym prekursorem (np. 110-aa preproinsulina P01308 lub 685-aa poliubikwityna P0CG48)
        Boolean(
          externalReference?.pdbText &&
          externalReference.pdbText.split(/\r?\n/).filter(l => l.startsWith('ATOM  ') && l.substring(12, 16).trim() === 'CA').length >
            Math.round(n * 1.45)
        ),
      iterationId,
      parameterVersion: paramVersion,
      modelArtifactPath: recordedArtifactPath,
      criteria: validationContext?.criteria
    }
  );

  const rmsdComparison = externalReference?.pdbText
    ? buildRmsdComparisonFromGroundTruth(groundTruthValidation, externalReference.pdbText)
    : null;

  // Wygeneruj finalny plik PDB zawierający w nagłówku rzeczywisty raport walidacji Ground Truth (lub NIEDOSTEPNY)
  const modelPdbText = exportStandardAllAtomPdb(
    `SNAPLOCK ITERATION ${iterationId} (${paramVersion})`,
    sequenceCodes,
    allAtomStructure,
    secondaryTypes,
    {
      pipelineMode: closureResult.pipelineMode,
      meanMicroFacetErrorDeg: internalQuantization.summary.l2MeanAngularErrorDeg,
      dMaxActual: finalBoundary.dMax,
      rgActual: finalBoundary.rg,
      cohesiveEnergy: parseFloat(baseEnergy.toFixed(1)),
      cartesianAlignment: internalQuantization.summary,
      groundTruthValidation,
      rmsdComparison
    }
  );

  const metrics: FoldingMetrics = {
    rgActual: finalBoundary.rg,
    dMaxActual: finalBoundary.dMax,
    rgTheoretical: finalBoundary.rgTheoretical,
    dMaxTheoretical: finalBoundary.dMaxTheoretical,
    alphaActualPct: cdResult.alphaActualPct,
    betaActualPct: cdResult.betaActualPct,
    cdPenaltyScore: cdResult.penaltyScore,
    disulfidePenaltyScore: totalSsbondPenalty,
    totalCohesiveEnergy: parseFloat(baseEnergy.toFixed(1)),
    ssBondPairs: activeSsPairs,
    pipelineMode: closureResult.pipelineMode,
    bridgeDiagnostics: updatedBridgeDiagnostics,
    initialMaxDeviationAngstrom: closureResult.initialMaxDeviationAngstrom,
    finalMaxSgDeviationAngstrom: parseFloat(finalMaxSgDev.toFixed(2)),
    finalMaxCaDeviationAngstrom: parseFloat(finalMaxCaDev.toFixed(2)),
    exponentialBoostTriggered: closureResult.exponentialBoostTriggered,
    initialSpringWeight: closureResult.initialSpringWeight,
    finalSpringWeight: closureResult.finalSpringWeight,
    meanMicroFacetErrorDeg: internalQuantization.summary.l2MeanAngularErrorDeg,
    vdwClashesCount: allAtomStructure.stericClashesCount,
    totalHeavyAtoms: allAtomStructure.totalHeavyAtoms,
    totalSidechainAtoms: allAtomStructure.totalSidechainAtoms,
    totalCovalentBonds: allAtomStructure.totalCovalentBonds,
    totalSaltBridges: allAtomStructure.totalSaltBridges,
    totalHydrogenBonds: allAtomStructure.totalHydrogenBonds,
    internalL2BackboneResidualAngstrom: internalQuantization.summary.postAlignmentBackboneRmsdAngstrom,
    cartesianAlignment: internalQuantization.summary,
    groundTruthValidation,
    parameterVersion: paramVersion,
    iterationId,
    generatorMode: 'LEGACY_E8_L2',
    restraintProvenance: 'DE_NOVO_SEQUENCE_ONLY',
    internalPhysicalEnergyKcal: parseFloat(baseEnergy.toFixed(1)),
    hierarchicalSummary: null
  };

  // Orientacje klocków (Łatka #1 do poprzednika, Łatka #248 do następnika)
  const orientations: THREE.Quaternion[] = [];
  for (let i = 0; i < n; i++) {
    const pCur = positions[i];
    const cur = new THREE.Vector3(pCur?.[0] ?? 0, pCur?.[1] ?? 0, pCur?.[2] ?? 0);
    const prev =
      i > 0 && positions[i - 1]
        ? new THREE.Vector3(positions[i - 1][0], positions[i - 1][1], positions[i - 1][2])
        : null;
    const next =
      i < n - 1 && positions[i + 1]
        ? new THREE.Vector3(positions[i + 1][0], positions[i + 1][1], positions[i + 1][2])
        : null;

    const axisY = new THREE.Vector3(0, 1, 0);
    if (prev && next) {
      axisY.subVectors(prev, next).normalize();
    } else if (prev) {
      axisY.subVectors(prev, cur).normalize();
    } else if (next) {
      axisY.subVectors(cur, next).normalize();
    }
    const q = new THREE.Quaternion().setFromUnitVectors(new THREE.Vector3(0, 1, 0), axisY);
    orientations.push(q);
  }

  return {
    positions,
    cbPositions,
    sgPositions,
    allAtomStructure,
    orientations,
    secondaryTypes,
    metrics,
    modelPdbText,
    groundTruthValidation,
    rmsdComparison
  };
}

// -----------------------------------------------------------------------------
// 6. STANDARDOWY EKSPOZITOR FORMATU PDB v3.3 (BEZ NADPISYWANIA WSPÓŁRZĘDNYCH REFERENCJĄ!)
// -----------------------------------------------------------------------------
export function exportToPdbString(
  headerTitle: string,
  sequenceCodes: string[],
  caPositions: [number, number, number][],
  metrics: FoldingMetrics,
  secondaryTypes: ('HELIX' | 'SHEET' | 'TURN' | 'COIL')[],
  prebuiltAllAtom?: AllAtomStructureResult | null,
  rmsdComparison?: AlphaFoldRmsdComparisonResult | null,
  referencePdbInput?: ExternalReferencePdbInput | null
): string {
  const validLen = Math.min(sequenceCodes?.length ?? 0, caPositions ? caPositions.length : 0);
  if (validLen === 0) return 'END';

  const slicedCodes = sequenceCodes.slice(0, validLen);
  const slicedCa = caPositions.slice(0, validLen);

  let allAtom = prebuiltAllAtom;
  if (!allAtom || allAtom.atoms.length === 0) {
    const activeSs = metrics?.ssBondPairs ?? inferUniversalDisulfideTopology(slicedCodes);
    const { cbPositions, sgPositions } = reconstructSidechainCbAndSg(slicedCodes, slicedCa, activeSs);
    allAtom = buildAllAtomStructure({
      sequenceCodes: slicedCodes,
      caPositions: slicedCa,
      cbOverrides: cbPositions,
      sgOverrides: sgPositions,
      ssBondPairs: activeSs,
      secondaryTypes
    });
  }

  // Jeżeli podano zewnętrzną referencję, a w metrics nie ma jeszcze walidacji dla tej struktury, oblicz ją bez modyfikacji atomów!
  let gtReport = metrics?.groundTruthValidation ?? null;
  const rawPdbForCheck = exportStandardAllAtomPdb(headerTitle, slicedCodes, allAtom, secondaryTypes, {
    pipelineMode: metrics?.pipelineMode,
    meanMicroFacetErrorDeg: metrics?.meanMicroFacetErrorDeg,
    dMaxActual: metrics?.dMaxActual,
    rgActual: metrics?.rgActual,
    cohesiveEnergy: metrics?.totalCohesiveEnergy,
    cartesianAlignment: metrics?.cartesianAlignment
  });

  if (referencePdbInput?.pdbText && (!gtReport || !gtReport.hasGroundTruthReference)) {
    gtReport = validateAgainstGroundTruth(rawPdbForCheck, referencePdbInput.pdbText, {
      referenceId: referencePdbInput.resolvedId,
      referenceName: referencePdbInput.referenceName ?? referencePdbInput.resolvedId,
      referenceChainId: referencePdbInput.chainId,
      referenceSource: referencePdbInput.source,
      iterationId: metrics?.iterationId ?? 'export',
      parameterVersion: metrics?.parameterVersion ?? 'v1.0.0'
    });
  }

  const effectiveRmsd =
    rmsdComparison ??
    (gtReport && referencePdbInput?.pdbText
      ? buildRmsdComparisonFromGroundTruth(gtReport, referencePdbInput.pdbText)
      : null);

  return exportStandardAllAtomPdb(headerTitle, slicedCodes, allAtom, secondaryTypes, {
    pipelineMode: metrics?.pipelineMode,
    meanMicroFacetErrorDeg: metrics?.meanMicroFacetErrorDeg,
    dMaxActual: metrics?.dMaxActual,
    rgActual: metrics?.rgActual,
    cohesiveEnergy: metrics?.totalCohesiveEnergy,
    rmsdComparison: effectiveRmsd,
    cartesianAlignment: metrics?.cartesianAlignment,
    groundTruthValidation: gtReport
  });
}

// -----------------------------------------------------------------------------
// 7. FUNKCJA POBIERANIA PLIKU .PDB DO PRZEGLĄDARKI
// -----------------------------------------------------------------------------
export function triggerPdbDownload(pdbContent: string, filename: string = 'snaplock_model.pdb') {
  const blob = new Blob([pdbContent], { type: 'chemical/x-pdb;charset=utf-8' });
  const url = URL.createObjectURL(blob);
  const anchor = document.createElement('a');
  anchor.href = url;
  anchor.download = filename.endsWith('.pdb') ? filename : `${filename}.pdb`;
  document.body.appendChild(anchor);
  anchor.click();
  document.body.removeChild(anchor);
  URL.revokeObjectURL(url);
}
