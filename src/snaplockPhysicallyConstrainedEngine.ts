// SNAPLOCK Phase 2 — Physically Constrained Hierarchical Folding Engine
// ============================================================================
// Implements:
// 1. Explicit residue-level internal coordinates (φ_i, ψ_i, ω_i) for i = 1..N
//    with exact Engh & Huber covalent NeRF backbone geometry (N, CA, C, O, H, CB).
// 2. Sequence-driven secondary structure segmentation (α-helices, β-strands,
//    canonical Type I / I' / II / II' β-turns, and flexible loops) using the
//    22-parameter SNAPLOCK physicochemical database + conformational propensities.
// 3. Explicit long-range contact constraints (β-strand pairing, amphipathic
//    hydrophobic core contacts, salt bridges, covalent S-S clamps, and optional
//    explicitly labeled external restraints) strictly separated from Ground-Truth
//    evaluation (zero reference coordinate leakage in de-novo & Held-Out Test).
// 4. Directional backbone hydrogen-bond restraints (O_i ... H_j - N_j) and
//    stereochemical disulfide loop closure (d(Sγ,Sγ) ≈ 2.05 Å, d(Cα,Cα) ≈ 5.50 Å).
// 5. 5-Stage Hierarchical Optimization (Stages A -> B -> C -> D -> E) across
//    multiple topological candidates with topology-preserving atomic relaxation,
//    steric clash removal, and independent Ground-Truth validation.

import * as THREE from 'three';
import { SNAPLOCK_21_BLOCKS } from './snaplock21Blocks';
import { E8_PATCHES } from './snaplockEngine';
import {
  AllAtomStructureResult,
  AlphaFoldRmsdComparisonResult,
  ExplicitResidueBackboneCoords,
  buildAllAtomStructure,
  buildRmsdComparisonFromGroundTruth,
  exportStandardAllAtomPdb,
  placeAtomNeRF
} from './snaplockAllAtomBuilder';
import {
  TARGET_CA_CA_DIST,
  TARGET_SG_SG_DIST,
  getPolarMicroLattice,
  inferUniversalDisulfideTopology,
  reconstructSidechainCbAndSg,
  snapVectorToDiscreteMicroState
} from './snaplockLoopClosureSolver';
import {
  GroundTruthPassCriteria,
  GroundTruthValidationReport,
  ParsedPdbStructure,
  parsePdbForGroundTruth,
  validateAgainstGroundTruth
} from './ground_truth_validation';

const DEG2RAD = Math.PI / 180.0;
const RAD2DEG = 180.0 / Math.PI;

// -----------------------------------------------------------------------------
// 1. ARCHITECTURAL DIAGNOSIS OF LEGACY GENERATOR LIMITATIONS
// -----------------------------------------------------------------------------

export interface GeneratorLimitationDiagnostic {
  id: string;
  component: string;
  limitationTitle: string;
  rootCauseInLegacyCode: string;
  structuralImpact: string;
  phase2Resolution: string;
}

export function diagnoseLegacyGeneratorLimitations(): GeneratorLimitationDiagnostic[] {
  return [
    {
      id: 'LIM-01',
      component: 'snaplockFastaAndPdb.ts (linie 443–518)',
      limitationTitle: 'Brak niezależnych kątów torsyjnych (φ_i, ψ_i) dla poszczególnych reszt',
      rootCauseInLegacyCode:
        'Wszystkie reszty współdzieliły 2 globalne pary skalarów (helixPhiDeg/helixPsiDeg, sheetPhiDeg/sheetPsiDeg) z globalnym przesunięciem (phiOffsetDeg, psiOffsetDeg), przekształcane przez ad-hoc formułę rawBendDeg = |(180 + φ + ψ) mod 120|.',
      structuralImpact:
        'Przestrzeń konformacyjna była zredukowana z 2N stopni swobody do ~6 globalnych skalarów; ponadto φ + ψ występowało sumarycznie, czyniąc phiOffsetDeg i psiOffsetDeg niemal współliniowymi.',
      phase2Resolution:
        'Wprowadzono jawną tablicę ResidueTorsionState[] z niezależnymi (φ_i, ψ_i, ω_i) dla każdej reszty i = 1..N oraz dokładną kinematykę kowalencyjną NeRF dla atomów (N, CA, C, O, H, CB).'
    },
    {
      id: 'LIM-02',
      component: 'snaplockFastaAndPdb.ts (linie 467–475)',
      limitationTitle: 'Naiwny podział struktury drugorzędowej według ułamka pozycji i / N',
      rootCauseInLegacyCode:
        'Warunek normalizedPos < betaRatio wymuszał jeden ciągły odcinek SHEET na N-końcu i jeden ciągły odcinek HELIX w środku łańcucha, ignorując rzeczywisty rozkład motywów α/β w sekwencji.',
      structuralImpact:
        'Uniemożliwiało to odtworzenie przeplatanych topologii α+β (np. ubikwityny β1-β2-α1-β3-β4-β5 czy lizozymu α+β).',
      phase2Resolution:
        'Wdrożono segmentację sekwencyjną SecondaryStructureSegment[] opartą na skłonnościach konformacyjnych 21 aminokwasów, N-capach helis oraz 4-resztowych zwrotach β (typy I, I\', II, II\').'
    },
    {
      id: 'LIM-03',
      component: 'snaplockFastaAndPdb.ts (linie 527–562)',
      limitationTitle: 'Brak jawnych kontaktów dalekiego zasięgu i brak parowania nici arkusza β',
      rootCauseInLegacyCode:
        'Jedyny człon nielokalny przyciągał reszty hydrofobowe do wspólnego izotropowego centroidu z małym krokiem 0.045, bez grafu par kontaktowych (i, j) i bez drabin wiązań wodorowych między nićmi β.',
      structuralImpact:
        'Nici β uciekały od siebie w otwartą wstęgę zamiast tworzyć spinki β-hairpin i arkusze β (RMSD ~13.9–15.9 Å).',
      phase2Resolution:
        'Wprowadzono jawny graf kontaktów dalekiego zasięgu LongRangeContactConstraint[] (parowanie nici β, kontakty amfipatycznego rdzenia hydrofobowego, mostki solne, klamry S-S oraz opcjonalne jawne więzy zewnętrzne) oraz kierunkowy potencjał wiązań wodorowych O...H-N.'
    },
    {
      id: 'LIM-04',
      component: 'snaplockAllAtomBuilder.ts (linie 786–887)',
      limitationTitle: 'Rekonstrukcja szkieletu (N, C, O) post-hoc z heurystycznych normalnych Cα',
      rootCauseInLegacyCode:
        'Atomy N, C, O były dobudowywane dopiero po ustaleniu śladu Cα na podstawie iloczynu wektorowego sąsiednich wektorów Cα, bez powiązania z kątami Ramachandrana (φ_i, ψ_i).',
      structuralImpact:
        'Geometria grup peptydowych nie wynikała z rzeczywistych rotacji torsyjnych wokół wiązań N-Cα i Cα-C, co generowało kolizje steryczne grup karbonylowych i łańcuchów bocznych.',
      phase2Resolution:
        'Szkielet (N, CA, C, O, H, CB) jest budowany bezpośrednio przez kinematykę NeRF z (φ_i, ψ_i, ω_i) ze stałymi długościami wiązań Engha-Hubera i przekazywany wprost do buildAllAtomStructure.'
    },
    {
      id: 'LIM-05',
      component: 'snaplockAutocalibrationLoop.ts (linie 66–76, 578–599)',
      limitationTitle: 'Parametry nieaktywne (ΔX = 0 Å) i odrzucanie prób przez nierelaksowane kolizje ciężkich atomów',
      rootCauseInLegacyCode:
        'Optymalizator próbkował m.in. ssRestraintWeight na białkach bez cystein (1UBQ: ΔX = 0 Å) oraz dMaxRelaxation poniżej progu obwiedni, a próby zagęszczające strukturę były natychmiast odrzucane przez brak etapu relaksacji sterycznej (Stage D).',
      structuralImpact:
        'Ponad 85% prób optymalizatora kończyło się natychmiastowym wycofaniem (22/25 odrzuconych prób dla ubikwityny).',
      phase2Resolution:
        'Wdrożono detekcję parametrów nieaktywnych i redundantnych przed optymalizacją, wieloetapową relaksację steryczną Etapu D usuwającą kolizje bez niszczenia topologii oraz pełne logowanie ΔX, energii fizycznej i kolizji dla każdej próby.'
    }
  ];
}

// -----------------------------------------------------------------------------
// 2. DATA STRUCTURES FOR EXPANDED STRUCTURAL REPRESENTATION
// -----------------------------------------------------------------------------

export type FoldingGeneratorMode =
  | 'LEGACY_E8_L2'
  | 'INDEPENDENT_TORSIONS_ONLY'
  | 'SECONDARY_STRUCTURE_CONSTRAINED'
  | 'LONG_RANGE_CONTACTS_ONLY'
  | 'FULL_HIERARCHICAL_ENGINE';

export type RestraintProvenanceMode =
  | 'DE_NOVO_SEQUENCE_ONLY'
  | 'EVOLUTIONARY_MSA_DCA'
  | 'PREDICTED_CONTACT_MAP'
  | 'EXPERIMENTAL_NOE_OR_DISULFIDE'
  | 'CALIBRATION_TRAIN_SUPERVISED';

export type CandidateSelectionPolicy =
  | 'INTERNAL_PHYSICS_ONLY'
  | 'CALIBRATION_TRAIN_SUPERVISED';

export type BetaTurnSubtype =
  | 'TYPE_I'
  | 'TYPE_I_PRIME'
  | 'TYPE_II'
  | 'TYPE_II_PRIME'
  | 'FLEX_HINGE';

export type RamachandranBasin =
  | 'ALPHA_R'
  | 'BETA'
  | 'PII'
  | 'ALPHA_L'
  | 'TURN_FAVORED'
  | 'PROLINE_RESTRICTED'
  | 'GLYCINE_PERMISSIVE'
  | 'OUTLIER';

export interface ResidueTorsionState {
  resIndex: number;          // 0-based
  resSeq: number;            // 1-based
  resName: string;           // 3-letter AA
  code1: string;             // 1-letter AA
  phiDeg: number;            // [-180, 180]
  psiDeg: number;            // [-180, 180]
  omegaDeg: number;          // Typically 180° (trans) or 0° (cis-PRO)
  secondaryType: 'HELIX' | 'SHEET' | 'TURN' | 'COIL';
  turnSubtype?: BetaTurnSubtype;
  ramachandranBasin: RamachandranBasin;
  ramachandranFavored: boolean;
  ramachandranEnergyKcal: number;
  kinematicFlexibility: number;
  e8NearestPatchId: number;
  l2MicroFacetErrorDeg: number;
}

export interface SecondaryStructureSegment {
  segmentId: string;
  type: 'HELIX' | 'SHEET' | 'TURN' | 'COIL';
  startIdx: number;          // 0-based inclusive
  endIdx: number;            // 0-based inclusive
  length: number;
  meanPropensityScore: number;
  turnSubtype?: BetaTurnSubtype;
  hydrophobicMomentMag: number;
}

export interface BetaSheetPairing {
  strandAId: string;
  strandBId: string;
  strandARange: [number, number];
  strandBRange: [number, number];
  orientation: 'ANTIPARALLEL' | 'PARALLEL';
  residuePairs: [number, number][]; // 0-based paired residue indices
  isHairpin: boolean;
}

export type ContactConstraintCategory =
  | 'DISULFIDE_CLAMP'
  | 'BETA_STRAND_PAIR'
  | 'HYDROPHOBIC_CORE'
  | 'SALT_BRIDGE'
  | 'EXTERNAL_RESTRAINT';

export interface LongRangeContactConstraint {
  contactId: string;
  resIdxA: number;
  resIdxB: number;
  resSeqA: number;
  resSeqB: number;
  sequenceSeparation: number;
  category: ContactConstraintCategory;
  targetCaDistAngstrom: number;
  targetCbDistAngstrom: number;
  toleranceAngstrom: number;
  weight: number;
  provenance: RestraintProvenanceMode;
  actualCaDistAngstrom?: number;
  actualCbDistAngstrom?: number;
  satisfied?: boolean;
}

export interface HydrogenBondConstraint {
  hBondId: string;
  donorResIdx: number;       // Residue providing N-H
  acceptorResIdx: number;    // Residue providing C=O
  hBondType: 'ALPHA_I_PLUS_4' | 'TURN_I_PLUS_3' | 'BETA_ANTIPARALLEL' | 'BETA_PARALLEL';
  targetONDistAngstrom: number;
  actualONDistAngstrom?: number;
  actualNHOAngleDeg?: number;
  energyKcal?: number;
  satisfied?: boolean;
}

export interface ExplicitExternalContactInput {
  resSeqA: number;           // 1-based residue index A
  resSeqB: number;           // 1-based residue index B
  targetDistAngstrom?: number;
  toleranceAngstrom?: number;
  confidenceWeight?: number;
}

export interface ExternalRestraintSource {
  sourceType: Exclude<RestraintProvenanceMode, 'DE_NOVO_SEQUENCE_ONLY'>;
  sourceLabel: string;
  contacts: ExplicitExternalContactInput[];
  secondaryStructure1Letter?: string; // Optional 'H'/'E'/'T'/'C' string of length N
}

export interface PhysicalEnergyWeights {
  ramaWeight: number;
  hBondWeight: number;
  longRangeContactWeight: number;
  betaPairingWeight: number;
  hydrophobicPackingWeight: number;
  compactnessRgWeight: number;
  stericRepulsionWeight: number;
  disulfideClampWeight: number;
  electrostaticWeight: number;
  turnNucleationBias: number;
}

export const DEFAULT_PHYSICAL_ENERGY_WEIGHTS: PhysicalEnergyWeights = {
  ramaWeight: 1.00,
  hBondWeight: 1.25,
  longRangeContactWeight: 1.15,
  betaPairingWeight: 1.35,
  hydrophobicPackingWeight: 1.05,
  compactnessRgWeight: 0.95,
  stericRepulsionWeight: 1.50,
  disulfideClampWeight: 3.80,
  electrostaticWeight: 0.65,
  turnNucleationBias: 1.00
};

export interface PhysicalEnergyBreakdown {
  ramachandranEnergyKcal: number;
  backboneHBondEnergyKcal: number;
  longRangeContactEnergyKcal: number;
  betaPairingEnergyKcal: number;
  hydrophobicCoreCompactnessKcal: number;
  disulfideClampEnergyKcal: number;
  electrostaticEnergyKcal: number;
  stericClashPenaltyKcal: number;
  totalInternalPhysicalEnergyKcal: number; // Explicitly internal physics score — NEVER Ground-Truth accuracy!
}

export interface NeRFBackboneAtomSet {
  N: THREE.Vector3;
  CA: THREE.Vector3;
  C: THREE.Vector3;
  O: THREE.Vector3;
  H: THREE.Vector3 | null;
  CB: THREE.Vector3 | null;
  OXT?: THREE.Vector3;
}

export interface HierarchicalStageTrace {
  stageCode: 'A' | 'B' | 'C' | 'D' | 'E';
  stageName: string;
  evaluationsExecuted: number;
  acceptedMoves: number;
  internalPhysicalEnergyKcal: number;
  rgAngstrom: number;
  satisfiedContactsCount: number;
  satisfiedHBondsCount: number;
  backboneStericOverlapsCount: number;
  heavyAtomClashesUnder155A: number;
  topologyDriftFromPrevStageRmsd: number;
  notes: string;
}

export interface CandidateEvaluationMetrics {
  // 1. Ground-Truth Similarity Metrics (External Reference Only — null when unavailable)
  groundTruth: {
    hasReference: boolean;
    caRmsdAngstrom: number | null;
    independentCaRmsdAngstrom: number | null;
    backboneRmsdAngstrom: number | null;
    allHeavyAtomRmsdAngstrom: number | null;
    tmScore: number | null;
    gdtTsPct: number | null;
    pairedResidues: number;
    coveragePct: number;
  };
  // 2. Local Covalent & Torsional Geometry
  localGeometry: {
    geometryValid: boolean;
    covalentBondViolationsCount: number;
    meanCaCaStepAngstrom: number;
    minCaCaStepAngstrom: number;
    maxCaCaStepAngstrom: number;
    ramachandranFavoredPct: number;
    ramachandranOutliersCount: number;
    meanL2MicroFacetErrorDeg: number;
  };
  // 3. Steric Clash Plausibility
  sterics: {
    nonBondedCaClashesCount: number;
    heavyAtomClashesUnder155A: number;
    softHeavyAtomContactsUnder215A: number;
    collapsedChainDetected: boolean;
  };
  // 4. Contact & Hydrogen-Bond Network Quality
  contacts: {
    restraintProvenance: RestraintProvenanceMode;
    totalTargetLongRangeContacts: number;
    satisfiedLongRangeContacts: number;
    longRangeSatisfactionPct: number;
    totalTargetHBonds: number;
    satisfiedHBonds: number;
    disulfideBridgesCount: number;
    maxDisulfideSgErrorAngstrom: number;
    // Post-hoc native 8.0 Å Cα contact map accuracy (computed ONLY during evaluation if reference exists)
    postHocNativeContactPrecisionPct: number | null;
    postHocNativeContactRecallPct: number | null;
    postHocNativeContactF1Pct: number | null;
  };
  // 5. Internal Physical Model Score (Explicitly separated from Ground-Truth accuracy!)
  internalModelPhysics: PhysicalEnergyBreakdown;
}

export interface HierarchicalCandidateRecord {
  candidateIndex: number;
  candidateId: string;
  topologySeedLabel: string;
  parameterVersion: string;
  selectedAsBest: boolean;
  selectionReason: string;
  residueTorsions: ResidueTorsionState[];
  caPositions: [number, number, number][];
  explicitBackbone: ExplicitResidueBackboneCoords[];
  cbPositions: ([number, number, number] | null)[];
  sgPositions: ([number, number, number] | null)[];
  allAtomStructure: AllAtomStructureResult;
  pdbText: string;
  modelPdbHash: string;
  stageTraces: HierarchicalStageTrace[];
  evaluation: CandidateEvaluationMetrics;
  groundTruthValidation: GroundTruthValidationReport;
  rmsdComparison: AlphaFoldRmsdComparisonResult | null;
}

export interface PhysicallyConstrainedFoldingOptions {
  generatorMode?: FoldingGeneratorMode;
  parameterVersion?: string;
  iterationId?: string;
  weights?: Partial<PhysicalEnergyWeights>;
  cdAlphaTargetPct?: number;
  cdBetaTargetPct?: number;
  dMaxTarget?: number;
  rgTarget?: number;
  explicitDisulfidePairs?: [number, number][];
  rotamerOverrides?: Record<number, number>;
  externalRestraints?: ExternalRestraintSource | null;
  // Ground-Truth reference strictly for Stage E post-generation independent validation!
  groundTruthReference?: {
    pdbText: string;
    resolvedId: string;
    referenceName?: string;
    chainId?: string;
    allowSubsequenceAlignment?: boolean;
    source: 'RCSB_PDB' | 'ALPHAFOLD_DB' | 'CUSTOM_PDB';
  } | null;
  selectionPolicy?: CandidateSelectionPolicy;
  candidateCount?: number;
  stageCMovesBudget?: number;
  criteria?: Partial<GroundTruthPassCriteria>;
  artifactSink?: (candidateId: string, pdbText: string) => string | void;
}

export interface PhysicallyConstrainedFoldingResult {
  generatorMode: FoldingGeneratorMode;
  restraintProvenance: RestraintProvenanceMode;
  usedGroundTruthReferenceDuringGeneration: false;
  usedGroundTruthReferenceForCandidateSelection: boolean;
  deNovoLimitationNotice: string | null;
  parameterVersion: string;
  iterationId: string;
  computeTimeMs: number;
  totalFunctionEvaluations: number;

  // Topology & Restraints Representation
  secondarySegments: SecondaryStructureSegment[];
  betaSheetPairings: BetaSheetPairing[];
  longRangeContacts: LongRangeContactConstraint[];
  hydrogenBondConstraints: HydrogenBondConstraint[];
  disulfidePairs: [number, number][];

  // Winning Verified Candidate
  selectedCandidate: HierarchicalCandidateRecord;
  allCandidates: HierarchicalCandidateRecord[];
  architecturalLimitationsDiagnosed: GeneratorLimitationDiagnostic[];
}

// -----------------------------------------------------------------------------
// 3. DETERMINISTIC PRNG & GEOMETRIC UTILITIES
// -----------------------------------------------------------------------------

function wrapAngleDeg(deg: number): number {
  let a = ((deg + 180.0) % 360.0 + 360.0) % 360.0 - 180.0;
  if (a === -180.0) a = 180.0;
  return a;
}

function angleDiffDeg(a: number, b: number): number {
  return wrapAngleDeg(a - b);
}

class DeterministicRng {
  private state: number;
  constructor(seed: number) {
    this.state = (seed >>> 0) || 0x9e3779b9;
  }
  nextFloat(): number {
    // Mulberry32 deterministic PRNG
    let t = (this.state += 0x6d2b79f5);
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  }
  nextRange(min: number, max: number): number {
    return min + (max - min) * this.nextFloat();
  }
  nextInt(maxExclusive: number): number {
    return Math.floor(this.nextFloat() * Math.max(1, maxExclusive));
  }
}

function hashSequenceSeed(sequenceCodes: string[], salt: number = 0): number {
  let h = (2166136261 ^ salt) >>> 0;
  for (let i = 0; i < sequenceCodes.length; i++) {
    const s = sequenceCodes[i] || 'ALA';
    for (let c = 0; c < s.length; c++) {
      h ^= s.charCodeAt(c);
      h = Math.imul(h, 16777619) >>> 0;
    }
  }
  return h || 1337;
}

const THREE_TO_ONE: Record<string, string> = {
  ALA: 'A', ARG: 'R', ASN: 'N', ASP: 'D', CYS: 'C',
  GLU: 'E', GLN: 'Q', GLY: 'G', HIS: 'H', ILE: 'I',
  LEU: 'L', LYS: 'K', MET: 'M', PHE: 'F', PRO: 'P',
  SER: 'S', THR: 'T', TRP: 'W', TYR: 'Y', VAL: 'V',
  SEC: 'U'
};

// -----------------------------------------------------------------------------
// 4. RAMACHANDRAN BASIN & KINEMATIC FLEXIBILITY MODEL
// -----------------------------------------------------------------------------

export function evaluateRamachandranState(
  code3: string,
  phiDeg: number,
  psiDeg: number
): {
  basin: RamachandranBasin;
  favored: boolean;
  energyKcal: number;
} {
  const phi = wrapAngleDeg(phiDeg);
  const psi = wrapAngleDeg(psiDeg);

  if (code3 === 'GLY') {
    // Achiral Glycine: symmetric basins, only severe C=O...C=O syn-clash around (0, ±180) / (0, 0) penalized
    const nearZeroClash = Math.abs(phi) < 28 && Math.abs(psi) < 35;
    if (nearZeroClash) {
      return { basin: 'OUTLIER', favored: false, energyKcal: 3.8 };
    }
    return { basin: 'GLYCINE_PERMISSIVE', favored: true, energyKcal: -0.35 };
  }

  if (code3 === 'PRO') {
    // Pyrrolidine ring locks φ around -63° ± 18°
    const dPhi = Math.abs(angleDiffDeg(phi, -63.0));
    const dPsiAlpha = Math.abs(angleDiffDeg(psi, -38.0));
    const dPsiPii = Math.abs(angleDiffDeg(psi, 145.0));
    const bestPsiErr = Math.min(dPsiAlpha, dPsiPii);
    if (dPhi <= 22.0 && bestPsiErr <= 42.0) {
      return {
        basin: 'PROLINE_RESTRICTED',
        favored: true,
        energyKcal: -0.45 + 0.002 * (dPhi * dPhi + 0.4 * bestPsiErr * bestPsiErr)
      };
    }
    return {
      basin: 'OUTLIER',
      favored: false,
      energyKcal: 2.5 + 0.025 * dPhi * dPhi + 0.004 * bestPsiErr * bestPsiErr
    };
  }

  // Standard L-amino acid Ramachandran basins:
  // 1) Right-handed α-helix / 3_10 basin: (-62°, -42°)
  const dAlpha = Math.hypot(angleDiffDeg(phi, -62.0) / 28.0, angleDiffDeg(psi, -42.0) / 30.0);
  // 2) Extended β-sheet basin: (-122°, +132°)
  const dBeta = Math.hypot(angleDiffDeg(phi, -122.0) / 36.0, angleDiffDeg(psi, 132.0) / 36.0);
  // 3) Polyproline II (PII) basin: (-75°, +145°)
  const dPii = Math.hypot(angleDiffDeg(phi, -75.0) / 26.0, angleDiffDeg(psi, 145.0) / 28.0);
  // 4) Bridge / Type I-II turn region: (-85°, 0°)
  const dTurn = Math.hypot(angleDiffDeg(phi, -85.0) / 32.0, angleDiffDeg(psi, -5.0) / 32.0);
  // 5) Left-handed α_L turn basin (accessible to ASN, ASP, SER, LYS, HIS in tight turns): (+58°, +42°)
  const alphaLPermissive = code3 === 'ASN' || code3 === 'ASP' || code3 === 'SER' || code3 === 'LYS' || code3 === 'HIS';
  const dAlphaL = Math.hypot(angleDiffDeg(phi, 58.0) / 22.0, angleDiffDeg(psi, 42.0) / 25.0);

  let bestNorm = dAlpha;
  let basin: RamachandranBasin = 'ALPHA_R';
  if (dBeta < bestNorm) {
    bestNorm = dBeta;
    basin = 'BETA';
  }
  if (dPii < bestNorm) {
    bestNorm = dPii;
    basin = 'PII';
  }
  if (dTurn < bestNorm) {
    bestNorm = dTurn;
    basin = 'TURN_FAVORED';
  }
  if (dAlphaL < bestNorm) {
    bestNorm = dAlphaL;
    basin = 'ALPHA_L';
  }

  // Amino-acid intrinsic conformational preference from SNAPLOCK 22-parameter table
  const block = SNAPLOCK_21_BLOCKS[code3] || SNAPLOCK_21_BLOCKS['ALA'];
  let prefBonus = 0;
  if (basin === 'BETA' && block.physico.ramachandranPref.includes('Beta')) prefBonus = -0.30;
  if (basin === 'ALPHA_R' && block.physico.ramachandranPref.includes('Helisa')) prefBonus = -0.30;
  if (basin === 'ALPHA_L' && !alphaLPermissive) prefBonus += 1.35; // β-branched VAL/ILE/THR strongly disfavor α_L

  if (bestNorm <= 1.35) {
    return {
      basin,
      favored: true,
      energyKcal: parseFloat((-0.55 + 0.45 * bestNorm * bestNorm + prefBonus).toFixed(4))
    };
  }
  if (bestNorm <= 1.95) {
    return {
      basin,
      favored: true,
      energyKcal: parseFloat((0.35 + 0.75 * (bestNorm - 1.35) ** 2 + prefBonus).toFixed(4))
    };
  }
  return {
    basin: 'OUTLIER',
    favored: false,
    energyKcal: parseFloat((2.60 + 1.85 * (bestNorm - 1.95) ** 2 + Math.max(0, prefBonus)).toFixed(4))
  };
}

// -----------------------------------------------------------------------------
// 5. EXACT ENGH & HUBER NeRF BACKBONE FORWARD KINEMATICS FROM (φ_i, ψ_i, ω_i)
// -----------------------------------------------------------------------------

export function buildNeRFBackboneFromTorsions(
  sequenceCodes: string[],
  phiDegArr: number[],
  psiDegArr: number[],
  omegaDegArr?: number[]
): NeRFBackboneAtomSet[] {
  const n = sequenceCodes.length;
  if (n === 0) return [];

  const atomsByRes: NeRFBackboneAtomSet[] = new Array(n);

  // Residue 0 canonical frame:
  // |N0 - CA0| = 1.458 Å, |CA0 - C0| = 1.525 Å, ∠(N0 - CA0 - C0) = 111.2°
  const angle0Rad = 111.2 * DEG2RAD;
  const n0 = new THREE.Vector3(1.458 * Math.cos(angle0Rad), 1.458 * Math.sin(angle0Rad), 0.0);
  const ca0 = new THREE.Vector3(0.0, 0.0, 0.0);
  const c0 = new THREE.Vector3(1.525, 0.0, 0.0);
  const h0 = new THREE.Vector3(
    n0.x - 1.010 * Math.cos(35.0 * DEG2RAD),
    n0.y + 1.010 * Math.sin(35.0 * DEG2RAD),
    0.0
  );

  atomsByRes[0] = {
    N: n0,
    CA: ca0,
    C: c0,
    O: new THREE.Vector3(),
    H: sequenceCodes[0] === 'PRO' ? null : h0,
    CB: null
  };

  for (let i = 0; i < n; i++) {
    const cur = atomsByRes[i];
    const codeCur = sequenceCodes[i] || 'ALA';
    const psiI = wrapAngleDeg(psiDegArr[i] ?? -45.0);

    // L-amino acid tetrahedral Cβ placement via CORN rule
    if (codeCur !== 'GLY') {
      const cbAngle = codeCur === 'PRO' ? 103.5 : 110.5;
      const cbImproper = codeCur === 'PRO' ? 118.0 : 122.6;
      cur.CB = placeAtomNeRF(cur.N, cur.C, cur.CA, 1.530, cbAngle, cbImproper);
    }

    if (i < n - 1) {
      const nextCode = sequenceCodes[i + 1] || 'ALA';
      const omegaI = wrapAngleDeg(omegaDegArr?.[i] ?? 180.0);
      const phiNext =
        nextCode === 'PRO'
          ? Math.max(-82.0, Math.min(-46.0, wrapAngleDeg(phiDegArr[i + 1] ?? -63.0)))
          : wrapAngleDeg(phiDegArr[i + 1] ?? -60.0);

      // Place N(i+1) using ψ(i) = τ(N_i, CA_i, C_i, N_{i+1})
      const nNext = placeAtomNeRF(cur.N, cur.CA, cur.C, 1.329, 116.2, psiI);
      // Place carbonyl O(i) planar-trans to N(i+1) around CA_i - C_i
      cur.O = placeAtomNeRF(cur.N, cur.CA, cur.C, 1.231, 120.8, wrapAngleDeg(psiI + 180.0));

      // Place CA(i+1) using ω(i) = τ(CA_i, C_i, N_{i+1}, CA_{i+1})
      const caNext = placeAtomNeRF(cur.CA, cur.C, nNext, 1.458, 121.7, omegaI);
      // Place amide H(i+1) planar-trans to CA(i+1) around C_i - N_{i+1}
      const hNext =
        nextCode === 'PRO'
          ? null
          : placeAtomNeRF(cur.CA, cur.C, nNext, 1.010, 119.2, wrapAngleDeg(omegaI - 180.0));

      // Place C(i+1) using φ(i+1) = τ(C_i, N_{i+1}, CA_{i+1}, C_{i+1})
      const tauNCaC = nextCode === 'PRO' ? 110.2 : 111.2;
      const cNext = placeAtomNeRF(cur.C, nNext, caNext, 1.525, tauNCaC, phiNext);

      atomsByRes[i + 1] = {
        N: nNext,
        CA: caNext,
        C: cNext,
        O: new THREE.Vector3(),
        H: hNext,
        CB: null
      };
    } else {
      // C-terminal carboxylate O and OXT
      cur.O = placeAtomNeRF(cur.N, cur.CA, cur.C, 1.249, 120.2, wrapAngleDeg(psiI + 180.0));
      cur.OXT = placeAtomNeRF(cur.N, cur.CA, cur.C, 1.249, 119.8, psiI);
    }
  }

  // Center entire chain at Cartesian origin (0, 0, 0) via pure rigid translation
  const com = new THREE.Vector3(0, 0, 0);
  for (let i = 0; i < n; i++) {
    com.add(atomsByRes[i].CA);
  }
  com.divideScalar(n);

  for (let i = 0; i < n; i++) {
    const r = atomsByRes[i];
    r.N.sub(com);
    r.CA.sub(com);
    r.C.sub(com);
    r.O.sub(com);
    if (r.H) r.H.sub(com);
    if (r.CB) r.CB.sub(com);
    if (r.OXT) r.OXT.sub(com);
  }

  return atomsByRes;
}

// -----------------------------------------------------------------------------
// 6. SEQUENCE-DRIVEN SECONDARY STRUCTURE & TURN SEGMENTATION
// -----------------------------------------------------------------------------

// Chou-Fasman / Levitt empirical conformational propensities (P_alpha, P_beta, P_turn)
const CONFORMATIONAL_PROPENSITY: Record<string, { pAlpha: number; pBeta: number; pTurn: number }> = {
  ALA: { pAlpha: 1.42, pBeta: 0.83, pTurn: 0.66 },
  ARG: { pAlpha: 0.98, pBeta: 0.93, pTurn: 0.95 },
  ASN: { pAlpha: 0.67, pBeta: 0.89, pTurn: 1.56 },
  ASP: { pAlpha: 1.01, pBeta: 0.54, pTurn: 1.46 },
  CYS: { pAlpha: 0.70, pBeta: 1.19, pTurn: 1.19 },
  GLN: { pAlpha: 1.11, pBeta: 1.10, pTurn: 0.98 },
  GLU: { pAlpha: 1.51, pBeta: 0.37, pTurn: 0.74 },
  GLY: { pAlpha: 0.57, pBeta: 0.75, pTurn: 1.56 },
  HIS: { pAlpha: 1.00, pBeta: 0.87, pTurn: 0.95 },
  ILE: { pAlpha: 1.08, pBeta: 1.60, pTurn: 0.47 },
  LEU: { pAlpha: 1.21, pBeta: 1.30, pTurn: 0.59 },
  LYS: { pAlpha: 1.16, pBeta: 0.74, pTurn: 1.01 },
  MET: { pAlpha: 1.45, pBeta: 1.05, pTurn: 0.60 },
  PHE: { pAlpha: 1.13, pBeta: 1.38, pTurn: 0.60 },
  PRO: { pAlpha: 0.57, pBeta: 0.55, pTurn: 1.52 },
  SER: { pAlpha: 0.77, pBeta: 0.75, pTurn: 1.43 },
  THR: { pAlpha: 0.83, pBeta: 1.19, pTurn: 0.96 },
  TRP: { pAlpha: 1.08, pBeta: 1.37, pTurn: 0.96 },
  TYR: { pAlpha: 0.69, pBeta: 1.47, pTurn: 1.14 },
  VAL: { pAlpha: 1.06, pBeta: 1.70, pTurn: 0.50 },
  SEC: { pAlpha: 0.72, pBeta: 1.18, pTurn: 1.15 }
};

export function predictSecondaryStructureTopology(
  sequenceCodes: string[],
  options?: {
    cdAlphaTargetPct?: number;
    cdBetaTargetPct?: number;
    turnBias?: number;
    externalSS1Letter?: string;
  }
): {
  perResidueTypes: ('HELIX' | 'SHEET' | 'TURN' | 'COIL')[];
  perResidueTurnSubtypes: (BetaTurnSubtype | undefined)[];
  segments: SecondaryStructureSegment[];
  betaPairings: BetaSheetPairing[];
} {
  const n = sequenceCodes.length;
  const perResidueTypes: ('HELIX' | 'SHEET' | 'TURN' | 'COIL')[] = new Array(n).fill('COIL');
  const perResidueTurnSubtypes: (BetaTurnSubtype | undefined)[] = new Array(n).fill(undefined);

  if (n === 0) {
    return { perResidueTypes, perResidueTurnSubtypes, segments: [], betaPairings: [] };
  }

  // Optional explicit external secondary structure restraint if provided
  if (options?.externalSS1Letter && options.externalSS1Letter.length === n) {
    for (let i = 0; i < n; i++) {
      const ch = options.externalSS1Letter[i].toUpperCase();
      if (ch === 'H' || ch === 'G' || ch === 'I') perResidueTypes[i] = 'HELIX';
      else if (ch === 'E' || ch === 'B') perResidueTypes[i] = 'SHEET';
      else if (ch === 'T' || ch === 'S') {
        perResidueTypes[i] = 'TURN';
        perResidueTurnSubtypes[i] = 'TYPE_I';
      } else perResidueTypes[i] = 'COIL';
    }
  } else {
    // Sequence-driven propensity smoothing with mild Far-UV CD composition prior
    const alphaPrior =
      options?.cdAlphaTargetPct !== undefined
        ? Math.max(-0.14, Math.min(0.16, (options.cdAlphaTargetPct - 30.0) * 0.003))
        : 0.0;
    const betaPrior =
      options?.cdBetaTargetPct !== undefined
        ? Math.max(-0.14, Math.min(0.18, (options.cdBetaTargetPct - 25.0) * 0.0035))
        : 0.0;
    const turnBias = options?.turnBias ?? 1.0;

    const alphaScore = new Array<number>(n).fill(0);
    const betaScore = new Array<number>(n).fill(0);
    const turnScore = new Array<number>(n).fill(0);

    for (let i = 0; i < n; i++) {
      // 5-residue window for α-helix, 5-residue window for β-strand
      let sumA = 0, cntA = 0;
      for (let k = Math.max(0, i - 2); k <= Math.min(n - 1, i + 2); k++) {
        const p = CONFORMATIONAL_PROPENSITY[sequenceCodes[k]] || CONFORMATIONAL_PROPENSITY.ALA;
        sumA += p.pAlpha;
        cntA++;
      }
      alphaScore[i] = sumA / Math.max(1, cntA) + alphaPrior;
      if (sequenceCodes[i] === 'PRO' && i > 0) alphaScore[i] -= 0.42; // Internal Proline breaks H-bond i -> i-4
      if (sequenceCodes[i] === 'GLY') alphaScore[i] -= 0.25;

      let sumB = 0, cntB = 0;
      for (let k = Math.max(0, i - 2); k <= Math.min(n - 1, i + 2); k++) {
        const p = CONFORMATIONAL_PROPENSITY[sequenceCodes[k]] || CONFORMATIONAL_PROPENSITY.ALA;
        sumB += p.pBeta;
        cntB++;
      }
      betaScore[i] = sumB / Math.max(1, cntB) + betaPrior;
      if (sequenceCodes[i] === 'PRO') betaScore[i] -= 0.35;
    }

    // 4-residue β-turn nucleation score at i .. i+3
    for (let i = 0; i <= n - 4; i++) {
      const c1 = sequenceCodes[i + 1];
      const c2 = sequenceCodes[i + 2];
      const p0 = CONFORMATIONAL_PROPENSITY[sequenceCodes[i]] || CONFORMATIONAL_PROPENSITY.ALA;
      const p1 = CONFORMATIONAL_PROPENSITY[c1] || CONFORMATIONAL_PROPENSITY.ALA;
      const p2 = CONFORMATIONAL_PROPENSITY[c2] || CONFORMATIONAL_PROPENSITY.ALA;
      const p3 = CONFORMATIONAL_PROPENSITY[sequenceCodes[i + 3]] || CONFORMATIONAL_PROPENSITY.ALA;

      let tVal = (0.20 * p0.pTurn + 0.35 * p1.pTurn + 0.35 * p2.pTurn + 0.10 * p3.pTurn) * turnBias;
      if (c1 === 'PRO') tVal += 0.32;
      if (c2 === 'GLY' || c2 === 'ASN' || c2 === 'ASP') tVal += 0.28;
      if (c1 === 'ASP' || c1 === 'ASN' || c1 === 'SER') tVal += 0.15;
      turnScore[i + 1] = Math.max(turnScore[i + 1], tVal);
      turnScore[i + 2] = Math.max(turnScore[i + 2], tVal);
    }

    // Initial residue classification
    for (let i = 1; i < n - 1; i++) {
      if (turnScore[i] >= 1.24 && turnScore[i] >= alphaScore[i] + 0.08 && turnScore[i] >= betaScore[i] + 0.05) {
        perResidueTypes[i] = 'TURN';
      } else if (alphaScore[i] >= 1.03 && alphaScore[i] >= betaScore[i] - 0.02) {
        perResidueTypes[i] = 'HELIX';
      } else if (betaScore[i] >= 1.04 && betaScore[i] > alphaScore[i]) {
        perResidueTypes[i] = 'SHEET';
      } else {
        perResidueTypes[i] = 'COIL';
      }
    }

    // Enforce minimum physical segment lengths (HELIX >= 4 residues, SHEET >= 3 residues)
    // and split overly long continuous β-strands (> 8 residues) at highest turn propensity site
    let idx = 0;
    while (idx < n) {
      const t = perResidueTypes[idx];
      let end = idx;
      while (end + 1 < n && perResidueTypes[end + 1] === t) end++;
      const len = end - idx + 1;

      if (t === 'HELIX' && len < 4) {
        for (let k = idx; k <= end; k++) perResidueTypes[k] = 'COIL';
      } else if (t === 'SHEET' && len < 3) {
        for (let k = idx; k <= end; k++) perResidueTypes[k] = 'COIL';
      } else if (t === 'SHEET' && len >= 9) {
        // Long β-segments in globular proteins almost always bend via a central β-hairpin turn!
        let bestSplit = Math.floor((idx + end) / 2);
        let bestT = -Infinity;
        for (let k = idx + 3; k <= end - 3; k++) {
          if (turnScore[k] > bestT) {
            bestT = turnScore[k];
            bestSplit = k;
          }
        }
        perResidueTypes[bestSplit] = 'TURN';
        if (bestSplit + 1 <= end - 2) perResidueTypes[bestSplit + 1] = 'TURN';
      }
      idx = end + 1;
    }

    // Assign canonical β-turn subtypes (Type I, I', II) for TURN residues
    for (let i = 1; i < n - 1; i++) {
      if (perResidueTypes[i] === 'TURN') {
        const nextIsTurn = i + 1 < n && perResidueTypes[i + 1] === 'TURN';
        if (nextIsTurn) {
          const c2 = sequenceCodes[i + 1];
          const prevIsSheet = i > 0 && perResidueTypes[i - 1] === 'SHEET';
          const afterIsSheet = i + 2 < n && perResidueTypes[i + 2] === 'SHEET';
          if (c2 === 'GLY') {
            perResidueTurnSubtypes[i] = 'TYPE_II';
            perResidueTurnSubtypes[i + 1] = 'TYPE_II';
          } else if (prevIsSheet && afterIsSheet && (c2 === 'ASN' || c2 === 'ASP' || c2 === 'GLY')) {
            perResidueTurnSubtypes[i] = 'TYPE_I_PRIME';
            perResidueTurnSubtypes[i + 1] = 'TYPE_I_PRIME';
          } else {
            perResidueTurnSubtypes[i] = 'TYPE_I';
            perResidueTurnSubtypes[i + 1] = 'TYPE_I';
          }
          i++; // Handled 2-residue turn core
        } else if (!perResidueTurnSubtypes[i]) {
          perResidueTurnSubtypes[i] = 'FLEX_HINGE';
        }
      }
    }
  }

  // Build contiguous SecondaryStructureSegment[] with Eisenberg hydrophobic moment
  const segments: SecondaryStructureSegment[] = [];
  let sIdx = 0;
  let segNum = 1;
  while (sIdx < n) {
    const sType = perResidueTypes[sIdx];
    let eIdx = sIdx;
    while (eIdx + 1 < n && perResidueTypes[eIdx + 1] === sType) eIdx++;
    const len = eIdx - sIdx + 1;

    // Eisenberg hydrophobic moment |μ_H| at δ = 100° (HELIX) or δ = 160° (SHEET)
    const deltaRad = (sType === 'HELIX' ? 100.0 : sType === 'SHEET' ? 160.0 : 120.0) * DEG2RAD;
    let mx = 0, my = 0, pSum = 0;
    for (let k = sIdx; k <= eIdx; k++) {
      const code = sequenceCodes[k] || 'ALA';
      const h = (SNAPLOCK_21_BLOCKS[code] || SNAPLOCK_21_BLOCKS.ALA).physico.hydropathyKyte;
      mx += h * Math.cos(deltaRad * (k - sIdx));
      my += h * Math.sin(deltaRad * (k - sIdx));
      const cp = CONFORMATIONAL_PROPENSITY[code] || CONFORMATIONAL_PROPENSITY.ALA;
      pSum += sType === 'HELIX' ? cp.pAlpha : sType === 'SHEET' ? cp.pBeta : cp.pTurn;
    }

    segments.push({
      segmentId: `SSE-${segNum++}-${sType}`,
      type: sType,
      startIdx: sIdx,
      endIdx: eIdx,
      length: len,
      meanPropensityScore: parseFloat((pSum / len).toFixed(3)),
      turnSubtype: sType === 'TURN' ? perResidueTurnSubtypes[sIdx] : undefined,
      hydrophobicMomentMag: parseFloat((Math.hypot(mx, my) / len).toFixed(3))
    });
    sIdx = eIdx + 1;
  }

  // Build β-sheet strand pairings (antiparallel β-hairpins for short turn gaps, parallel/antiparallel sheet cradle for distant strands)
  const strandSegments = segments.filter(s => s.type === 'SHEET' && s.length >= 3);
  const betaPairings: BetaSheetPairing[] = [];

  for (let s = 0; s < strandSegments.length - 1; s++) {
    const stA = strandSegments[s];
    const stB = strandSegments[s + 1];
    const gap = stB.startIdx - stA.endIdx - 1;
    const pairLen = Math.min(stA.length, stB.length);
    if (pairLen < 2) continue;

    if (gap <= 7) {
      // Canonical antiparallel β-hairpin pairing outward from the turn!
      const residuePairs: [number, number][] = [];
      for (let k = 0; k < pairLen; k++) {
        const rA = stA.endIdx - k;
        const rB = stB.startIdx + k;
        if (rB - rA >= 4) {
          residuePairs.push([rA, rB]);
        }
      }
      if (residuePairs.length > 0) {
        betaPairings.push({
          strandAId: stA.segmentId,
          strandBId: stB.segmentId,
          strandARange: [stA.startIdx, stA.endIdx],
          strandBRange: [stB.startIdx, stB.endIdx],
          orientation: 'ANTIPARALLEL',
          residuePairs,
          isHairpin: true
        });
      }
    } else if (gap <= 42) {
      // Non-local tertiary β-sheet assembly (e.g., Ubiquitin β1..β5 or β2..β3/β4)
      // Long intervening α-helix crossover favors PARALLEL (β-α-β Rossmann/Ubiquitin motif) when a HELIX sits in between!
      const hasInterveningHelix = segments.some(
        seg => seg.type === 'HELIX' && seg.startIdx > stA.endIdx && seg.endIdx < stB.startIdx && seg.length >= 6
      );
      const orientation: 'ANTIPARALLEL' | 'PARALLEL' = hasInterveningHelix ? 'PARALLEL' : 'ANTIPARALLEL';
      const residuePairs: [number, number][] = [];
      for (let k = 0; k < pairLen; k++) {
        const rA = orientation === 'PARALLEL' ? stA.startIdx + k : stA.endIdx - k;
        const rB = stB.startIdx + k;
        if (rB - rA >= 6) {
          residuePairs.push([rA, rB]);
        }
      }
      if (residuePairs.length > 0) {
        betaPairings.push({
          strandAId: stA.segmentId,
          strandBId: stB.segmentId,
          strandARange: [stA.startIdx, stA.endIdx],
          strandBRange: [stB.startIdx, stB.endIdx],
          orientation,
          residuePairs,
          isHairpin: false
        });
      }
    }
  }

  // Also check N-to-C terminal β-sheet closure (e.g., β1 with C-terminal β_last in β-grasp / ubiquitin-like folds!)
  if (strandSegments.length >= 3) {
    const firstStrand = strandSegments[0];
    const lastStrand = strandSegments[strandSegments.length - 1];
    const alreadyPaired = betaPairings.some(
      bp => bp.strandAId === firstStrand.segmentId && bp.strandBId === lastStrand.segmentId
    );
    if (!alreadyPaired && lastStrand.startIdx - firstStrand.endIdx >= 18) {
      const pairLen = Math.min(firstStrand.length, lastStrand.length, 4);
      const residuePairs: [number, number][] = [];
      for (let k = 0; k < pairLen; k++) {
        residuePairs.push([firstStrand.startIdx + k, lastStrand.startIdx + k]);
      }
      betaPairings.push({
        strandAId: firstStrand.segmentId,
        strandBId: lastStrand.segmentId,
        strandARange: [firstStrand.startIdx, firstStrand.endIdx],
        strandBRange: [lastStrand.startIdx, lastStrand.endIdx],
        orientation: 'PARALLEL',
        residuePairs,
        isHairpin: false
      });
    }
  }

  return {
    perResidueTypes,
    perResidueTurnSubtypes,
    segments,
    betaPairings
  };
}

// -----------------------------------------------------------------------------
// 7. EXPLICIT LONG-RANGE CONTACTS & HYDROGEN-BOND RESTRAINT GRAPH
//    (STRICT SEPARATION OF DE-NOVO / EXTERNAL RESTRAINTS FROM GROUND TRUTH!)
// -----------------------------------------------------------------------------

export function buildPhysicalRestraintGraph(
  sequenceCodes: string[],
  generatorMode: FoldingGeneratorMode,
  segments: SecondaryStructureSegment[],
  betaPairings: BetaSheetPairing[],
  disulfidePairs: [number, number][],
  externalRestraints?: ExternalRestraintSource | null
): {
  restraintProvenance: RestraintProvenanceMode;
  deNovoLimitationNotice: string | null;
  longRangeContacts: LongRangeContactConstraint[];
  hydrogenBonds: HydrogenBondConstraint[];
} {
  const n = sequenceCodes.length;
  const longRangeContacts: LongRangeContactConstraint[] = [];
  const hydrogenBonds: HydrogenBondConstraint[] = [];

  const enableSSE =
    generatorMode === 'SECONDARY_STRUCTURE_CONSTRAINED' ||
    generatorMode === 'FULL_HIERARCHICAL_ENGINE';
  const enableLongRange =
    generatorMode === 'LONG_RANGE_CONTACTS_ONLY' ||
    generatorMode === 'FULL_HIERARCHICAL_ENGINE';

  const restraintProvenance: RestraintProvenanceMode = externalRestraints
    ? externalRestraints.sourceType
    : 'DE_NOVO_SEQUENCE_ONLY';

  const deNovoLimitationNotice = externalRestraints
    ? null
    : 'Tryb DE-NOVO (wyłącznie z sekwencji aminokwasowej): brak zewnętrznych profili ewolucyjnych MSA/DCA oraz brak przewidywanych map kontaktów. Kontakty dalekiego zasięgu wnioskowane są wyłącznie z fizykochemii 21 aminokwasów, parowania arkuszy β, momentów hydrofobowych helis oraz topologii mostków S-S.';

  // 1. Secondary Structure Backbone Hydrogen Bonds (i+4 -> i α-helix, i+3 -> i turns, β-sheet ladders)
  if (enableSSE) {
    let hbIdx = 1;
    for (const seg of segments) {
      if (seg.type === 'HELIX' && seg.length >= 4) {
        for (let i = seg.startIdx; i + 4 <= seg.endIdx; i++) {
          if (sequenceCodes[i + 4] === 'PRO') continue; // Proline lacks amide H
          hydrogenBonds.push({
            hBondId: `HB-${hbIdx++}`,
            donorResIdx: i + 4,
            acceptorResIdx: i,
            hBondType: 'ALPHA_I_PLUS_4',
            targetONDistAngstrom: 2.95
          });
        }
      } else if (seg.type === 'TURN') {
        const turnStart = Math.max(0, seg.startIdx - 1);
        const turnEnd = Math.min(n - 1, turnStart + 3);
        if (turnEnd - turnStart === 3 && sequenceCodes[turnEnd] !== 'PRO') {
          hydrogenBonds.push({
            hBondId: `HB-${hbIdx++}`,
            donorResIdx: turnEnd,
            acceptorResIdx: turnStart,
            hBondType: 'TURN_I_PLUS_3',
            targetONDistAngstrom: 3.05
          });
        }
      }
    }

    // Inter-strand β-sheet hydrogen bond ladders (only when long-range / full hierarchical or local hairpin)
    for (const bp of betaPairings) {
      if (!enableLongRange && !bp.isHairpin) continue;
      for (let p = 0; p < bp.residuePairs.length; p++) {
        const [rA, rB] = bp.residuePairs[p];
        if (bp.orientation === 'ANTIPARALLEL') {
          // Alternating H-bonded pairs in antiparallel β-sheets (every 2nd pair forms 2 symmetric H-bonds)
          if (p % 2 === 0) {
            if (sequenceCodes[rB] !== 'PRO') {
              hydrogenBonds.push({
                hBondId: `HB-${hbIdx++}`,
                donorResIdx: rB,
                acceptorResIdx: rA,
                hBondType: 'BETA_ANTIPARALLEL',
                targetONDistAngstrom: 2.92
              });
            }
            if (sequenceCodes[rA] !== 'PRO') {
              hydrogenBonds.push({
                hBondId: `HB-${hbIdx++}`,
                donorResIdx: rA,
                acceptorResIdx: rB,
                hBondType: 'BETA_ANTIPARALLEL',
                targetONDistAngstrom: 2.92
              });
            }
          }
        } else {
          // Parallel β-sheet H-bond bridge: O(rA) ... H-N(rB) and N(rA+2) ... O(rB)
          if (p % 2 === 0 && sequenceCodes[rB] !== 'PRO') {
            hydrogenBonds.push({
              hBondId: `HB-${hbIdx++}`,
              donorResIdx: rB,
              acceptorResIdx: rA,
              hBondType: 'BETA_PARALLEL',
              targetONDistAngstrom: 2.95
            });
          }
        }
      }
    }
  }

  if (!enableLongRange) {
    return { restraintProvenance, deNovoLimitationNotice, longRangeContacts, hydrogenBonds };
  }

  let cIdx = 1;
  const pairSet = new Set<string>();
  const addContact = (
    a: number,
    b: number,
    category: ContactConstraintCategory,
    targetCa: number,
    targetCb: number,
    tol: number,
    weight: number,
    prov: RestraintProvenanceMode
  ) => {
    const i1 = Math.min(a, b);
    const i2 = Math.max(a, b);
    if (i1 < 0 || i2 >= n || i2 - i1 < 4) return;
    const key = `${i1}:${i2}`;
    if (pairSet.has(key)) return;
    pairSet.add(key);
    longRangeContacts.push({
      contactId: `LRC-${cIdx++}`,
      resIdxA: i1,
      resIdxB: i2,
      resSeqA: i1 + 1,
      resSeqB: i2 + 1,
      sequenceSeparation: i2 - i1,
      category,
      targetCaDistAngstrom: targetCa,
      targetCbDistAngstrom: targetCb,
      toleranceAngstrom: tol,
      weight,
      provenance: prov
    });
  };

  // A. Track 5 Highest-Priority Covalent Disulfide Bridges (S-S)
  for (const [i1, i2] of disulfidePairs) {
    addContact(i1, i2, 'DISULFIDE_CLAMP', TARGET_CA_CA_DIST, 3.85, 0.55, 4.2, restraintProvenance);
  }

  // B. Optional Explicit External Contact Restraints (when user/protocol supplies MSA/DCA/NOE restraints)
  if (externalRestraints && externalRestraints.contacts.length > 0) {
    for (const ext of externalRestraints.contacts) {
      const a = ext.resSeqA - 1;
      const b = ext.resSeqB - 1;
      const tgt = ext.targetDistAngstrom ?? 6.5;
      addContact(
        a,
        b,
        'EXTERNAL_RESTRAINT',
        tgt,
        Math.max(4.5, tgt - 0.6),
        ext.toleranceAngstrom ?? 1.2,
        ext.confidenceWeight ?? 2.2,
        externalRestraints.sourceType
      );
    }
    return { restraintProvenance, deNovoLimitationNotice, longRangeContacts, hydrogenBonds };
  }

  // C. De-Novo β-Sheet Strand Pairing Contacts (only when SSE pairing is active or in full hierarchical mode)
  if (generatorMode === 'FULL_HIERARCHICAL_ENGINE') {
    for (const bp of betaPairings) {
      for (const [rA, rB] of bp.residuePairs) {
        addContact(
          rA,
          rB,
          'BETA_STRAND_PAIR',
          bp.orientation === 'ANTIPARALLEL' ? 4.75 : 4.95,
          bp.orientation === 'ANTIPARALLEL' ? 4.65 : 4.85,
          0.75,
          bp.isHairpin ? 2.1 : 1.55,
          'DE_NOVO_SEQUENCE_ONLY'
        );
      }
    }
  }

  // D. De-Novo Amphipathic Hydrophobic Core Contacts across distinct SSEs / sequence blocks
  const coreAnchors: { idx: number; hydro: number; segIdx: number }[] = [];
  segments.forEach((seg, sIndex) => {
    for (let i = seg.startIdx; i <= seg.endIdx; i++) {
      const h = (SNAPLOCK_21_BLOCKS[sequenceCodes[i]] || SNAPLOCK_21_BLOCKS.ALA).physico.hydropathyKyte;
      if (h >= 1.8) {
        coreAnchors.push({ idx: i, hydro: h, segIdx: sIndex });
      }
    }
  });

  // Cap per-residue contact degree <= 2 so the de-novo contact graph stays 3D-realizable without frustration
  const degreeByRes = new Map<number, number>();
  const candidateHydroPairs: { a: number; b: number; score: number }[] = [];
  for (let u = 0; u < coreAnchors.length; u++) {
    for (let v = u + 1; v < coreAnchors.length; v++) {
      const ca = coreAnchors[u];
      const cb = coreAnchors[v];
      if (ca.segIdx === cb.segIdx) continue;
      const sep = Math.abs(cb.idx - ca.idx);
      if (sep < 7) continue;
      // Favor packing between core residues across neighboring/tertiary structural elements
      const sepBonus = sep >= 8 && sep <= 48 ? 1.25 : 0.85;
      const score = (ca.hydro + cb.hydro) * sepBonus;
      candidateHydroPairs.push({ a: ca.idx, b: cb.idx, score });
    }
  }
  candidateHydroPairs.sort((x, y) => y.score - x.score);

  const maxCoreContacts = Math.max(3, Math.min(24, Math.floor(n * 0.22)));
  let addedCore = 0;
  for (const hp of candidateHydroPairs) {
    if (addedCore >= maxCoreContacts) break;
    const degA = degreeByRes.get(hp.a) ?? 0;
    const degB = degreeByRes.get(hp.b) ?? 0;
    if (degA >= 2 || degB >= 2) continue;
    addContact(hp.a, hp.b, 'HYDROPHOBIC_CORE', 6.80, 6.10, 1.45, 0.95, 'DE_NOVO_SEQUENCE_ONLY');
    degreeByRes.set(hp.a, degA + 1);
    degreeByRes.set(hp.b, degB + 1);
    addedCore++;
  }

  // E. De-Novo Complementary Salt-Bridge Restraints (Track 3 Electrostatics)
  let addedSalt = 0;
  const maxSalt = Math.max(1, Math.min(6, Math.floor(n / 16)));
  for (let i = 0; i < n && addedSalt < maxSalt; i++) {
    const qA = (SNAPLOCK_21_BLOCKS[sequenceCodes[i]] || SNAPLOCK_21_BLOCKS.ALA).physico.chargePH74;
    if (Math.abs(qA) < 0.8) continue;
    for (let j = i + 6; j < Math.min(n, i + 26) && addedSalt < maxSalt; j++) {
      const qB = (SNAPLOCK_21_BLOCKS[sequenceCodes[j]] || SNAPLOCK_21_BLOCKS.ALA).physico.chargePH74;
      if (qA * qB < -0.5 && (degreeByRes.get(i) ?? 0) < 2 && (degreeByRes.get(j) ?? 0) < 2) {
        addContact(i, j, 'SALT_BRIDGE', 7.40, 6.60, 1.60, 0.55, 'DE_NOVO_SEQUENCE_ONLY');
        degreeByRes.set(i, (degreeByRes.get(i) ?? 0) + 1);
        degreeByRes.set(j, (degreeByRes.get(j) ?? 0) + 1);
        addedSalt++;
        break;
      }
    }
  }

  return {
    restraintProvenance,
    deNovoLimitationNotice,
    longRangeContacts,
    hydrogenBonds
  };
}

// -----------------------------------------------------------------------------
// 8. COMPOSITE INTERNAL PHYSICAL ENERGY EVALUATOR (ZERO REFERENCE ACCESS!)
// -----------------------------------------------------------------------------

export function evaluateInternalPhysicalEnergy(
  sequenceCodes: string[],
  phiDeg: number[],
  psiDeg: number[],
  backbone: NeRFBackboneAtomSet[],
  longRangeContacts: LongRangeContactConstraint[],
  hydrogenBonds: HydrogenBondConstraint[],
  disulfidePairs: [number, number][],
  weights: PhysicalEnergyWeights,
  targetRgOverride?: number
): {
  breakdown: PhysicalEnergyBreakdown;
  rgActual: number;
  satisfiedContactsCount: number;
  satisfiedHBondsCount: number;
  backboneStericOverlapsCount: number;
} {
  const n = sequenceCodes.length;
  if (n === 0) {
    return {
      breakdown: {
        ramachandranEnergyKcal: 0,
        backboneHBondEnergyKcal: 0,
        longRangeContactEnergyKcal: 0,
        betaPairingEnergyKcal: 0,
        hydrophobicCoreCompactnessKcal: 0,
        disulfideClampEnergyKcal: 0,
        electrostaticEnergyKcal: 0,
        stericClashPenaltyKcal: 0,
        totalInternalPhysicalEnergyKcal: 0
      },
      rgActual: 0,
      satisfiedContactsCount: 0,
      satisfiedHBondsCount: 0,
      backboneStericOverlapsCount: 0
    };
  }

  // 1. Ramachandran Torsion Energy
  let eRama = 0;
  for (let i = 1; i < n - 1; i++) {
    const r = evaluateRamachandranState(sequenceCodes[i], phiDeg[i], psiDeg[i]);
    eRama += r.energyKcal;
  }
  eRama *= weights.ramaWeight;

  // 2. Directional Backbone Hydrogen Bonds (O_acc ... H_don - N_don)
  let eHBond = 0;
  let satisfiedHBondsCount = 0;
  const hBondPairKeys = new Set<string>();

  for (const hb of hydrogenBonds) {
    const don = backbone[hb.donorResIdx];
    const acc = backbone[hb.acceptorResIdx];
    if (!don || !acc) continue;
    const iMin = Math.min(hb.donorResIdx, hb.acceptorResIdx);
    const iMax = Math.max(hb.donorResIdx, hb.acceptorResIdx);
    hBondPairKeys.add(`${iMin}:${iMax}`);

    const dON = acc.O.distanceTo(don.N);
    let angleFactor = 1.0;
    let angleDeg = 160.0;
    if (don.H) {
      const vNH = don.H.clone().sub(don.N).normalize();
      const vHO = acc.O.clone().sub(don.H).normalize();
      const cosAng = Math.max(-1, Math.min(1, vNH.dot(vHO)));
      // vNH points N->H and vHO points H->O: collinear linear H-bond has cosAng ≈ +1 (∠N-H...O ≈ 180°)
      angleDeg = Math.acos(-cosAng) * RAD2DEG;
      angleFactor = Math.max(0.0, (cosAng + 0.25) / 1.25);
    }

    hb.actualONDistAngstrom = parseFloat(dON.toFixed(3));
    hb.actualNHOAngleDeg = parseFloat(angleDeg.toFixed(1));

    const dev = dON - hb.targetONDistAngstrom;
    if (dON >= 2.55 && dON <= 3.45 && angleFactor > 0.35) {
      const well = -2.25 * angleFactor * Math.exp(-(dev * dev) / 0.22);
      eHBond += well;
      hb.energyKcal = parseFloat((well * weights.hBondWeight).toFixed(3));
      hb.satisfied = true;
      satisfiedHBondsCount++;
    } else {
      // Gentle guiding attraction toward H-bond capture radius without blowing up
      const pull = 0.65 * Math.min(8.0, (Math.max(0, Math.abs(dev) - 0.35)) ** 2 / (1.0 + 0.25 * Math.abs(dev)));
      eHBond += pull;
      hb.energyKcal = parseFloat((pull * weights.hBondWeight).toFixed(3));
      hb.satisfied = false;
    }
  }
  eHBond *= weights.hBondWeight;

  // 3. Explicit Long-Range Contacts, β-Pairing & Disulfide Clamps
  let eContact = 0;
  let eBeta = 0;
  let eDisulfide = 0;
  let eElectro = 0;
  let satisfiedContactsCount = 0;

  for (const c of longRangeContacts) {
    const resA = backbone[c.resIdxA];
    const resB = backbone[c.resIdxB];
    if (!resA || !resB) continue;

    const dCa = resA.CA.distanceTo(resB.CA);
    const ptA = resA.CB ?? resA.CA;
    const ptB = resB.CB ?? resB.CA;
    const dCb = ptA.distanceTo(ptB);

    c.actualCaDistAngstrom = parseFloat(dCa.toFixed(3));
    c.actualCbDistAngstrom = parseFloat(dCb.toFixed(3));

    const errCa = Math.max(0, Math.abs(dCa - c.targetCaDistAngstrom) - c.toleranceAngstrom);
    const errCb = Math.max(0, Math.abs(dCb - c.targetCbDistAngstrom) - c.toleranceAngstrom);
    const withinTol = errCa <= 0.05 && errCb <= 0.15;
    c.satisfied = withinTol;
    if (withinTol) satisfiedContactsCount++;

    if (c.category === 'DISULFIDE_CLAMP') {
      // Strong harmonic + Lorentzian funnel for covalent disulfide bridge reach
      const ssErr = (dCa - TARGET_CA_CA_DIST) ** 2 + 0.85 * (dCb - 3.85) ** 2;
      eDisulfide += c.weight * (withinTol ? -3.8 : 0.95 * ssErr);
    } else if (c.category === 'BETA_STRAND_PAIR') {
      const robustWell = withinTol
        ? -2.10 * c.weight
        : c.weight * ((errCa * errCa + 0.6 * errCb * errCb) / (1.0 + 0.08 * errCa * errCa));
      eBeta += robustWell;
    } else if (c.category === 'SALT_BRIDGE') {
      const electroTerm = withinTol
        ? -1.45 * c.weight
        : c.weight * ((errCb * errCb) / (1.0 + 0.18 * errCb * errCb));
      eElectro += electroTerm;
    } else {
      const lorentz = withinTol
        ? -1.80 * c.weight
        : c.weight * ((errCb * errCb) / (1.0 + 0.14 * errCb * errCb));
      eContact += lorentz;
    }
  }

  // Additional stereochemical Sγ-Sγ disulfide reconstruction term
  if (disulfidePairs.length > 0) {
    const caTuple: [number, number, number][] = backbone.map(b => [b.CA.x, b.CA.y, b.CA.z]);
    const { sgPositions } = reconstructSidechainCbAndSg(sequenceCodes, caTuple, disulfidePairs);
    for (const [i1, i2] of disulfidePairs) {
      const sg1 = sgPositions[i1];
      const sg2 = sgPositions[i2];
      if (sg1 && sg2) {
        const dSg = Math.hypot(sg1[0] - sg2[0], sg1[1] - sg2[1], sg1[2] - sg2[2]);
        const sgErr = Math.abs(dSg - TARGET_SG_SG_DIST);
        eDisulfide += sgErr < 0.35 ? -2.5 : 1.8 * sgErr * sgErr;
      }
    }
  }

  eBeta *= weights.betaPairingWeight;
  eContact *= weights.longRangeContactWeight;
  eDisulfide *= weights.disulfideClampWeight / 3.8;
  eElectro *= weights.electrostaticWeight;

  // 4. Globular Radius of Gyration & Amphipathic Hydrophobic Burial
  //    (Because backbone is centered at (0,0,0), sum(CA^2)/N is exact Rg^2!)
  let sumCaSq = 0;
  let hydroR2Sum = 0, hydroCnt = 0;
  let polarR2Sum = 0, polarCnt = 0;
  for (let i = 0; i < n; i++) {
    const ca = backbone[i].CA;
    const r2 = ca.lengthSq();
    sumCaSq += r2;
    const cb = backbone[i].CB ?? ca;
    const cbR = cb.length();
    const phys = (SNAPLOCK_21_BLOCKS[sequenceCodes[i]] || SNAPLOCK_21_BLOCKS.ALA).physico;
    if (phys.hydropathyKyte >= 1.5) {
      hydroR2Sum += cbR;
      hydroCnt++;
    } else if (phys.hydropathyKyte <= -1.2) {
      polarR2Sum += cbR;
      polarCnt++;
    }
  }
  const rgActual = Math.sqrt(sumCaSq / n);
  // Empirical native globular protein Flory law: Rg ≈ 2.28 * N^0.38 Å
  const rgTarget = targetRgOverride && targetRgOverride > 4.0 ? targetRgOverride : 2.28 * Math.pow(n, 0.38);
  const rgDiff = rgActual - rgTarget;
  let eHydroCore =
    weights.compactnessRgWeight * (rgDiff > 0 ? 1.35 * rgDiff * rgDiff : 2.20 * rgDiff * rgDiff);

  if (hydroCnt > 0 && polarCnt > 0) {
    const meanHydroR = hydroR2Sum / hydroCnt;
    const meanPolarR = polarR2Sum / polarCnt;
    // Reward hydrophobic sidechains sitting deeper inside core than charged/polar sidechains
    eHydroCore += weights.hydrophobicPackingWeight * 1.65 * (meanHydroR - 0.88 * meanPolarR);
  }

  // 5. Hard Van der Waals Steric Exclusion across all non-bonded (N, CA, C, O, CB) atoms
  let eSteric = 0;
  let backboneStericOverlapsCount = 0;

  for (let i = 0; i < n; i++) {
    const bA = backbone[i];
    for (let j = i + 2; j < n; j++) {
      const bB = backbone[j];
      const caDx = bB.CA.x - bA.CA.x;
      const caDy = bB.CA.y - bA.CA.y;
      const caDz = bB.CA.z - bA.CA.z;
      const caD2 = caDx * caDx + caDy * caDy + caDz * caDz;

      // Bounding sphere skip: max distance from CA to O/CB is ~2.45 Å -> 2*2.45 + 3.0 = 7.9 Å
      if (caD2 > 68.0) continue;

      const caDist = Math.sqrt(caD2);
      const minCaAllowed = j - i === 2 ? 4.85 : 3.75;
      if (caDist < minCaAllowed) {
        const pen = minCaAllowed - caDist;
        eSteric += 38.0 * pen * pen + (caDist < 3.10 ? 180.0 * (3.10 - caDist) ** 2 : 0);
        if (caDist < 3.20) backboneStericOverlapsCount++;
      }

      const isHBondPair = hBondPairKeys.has(`${i}:${j}`);
      const isDisulfidePair = disulfidePairs.some(([p1, p2]) => p1 === i && p2 === j);

      const atomsA: { name: string; pt: THREE.Vector3 }[] = [
        { name: 'N', pt: bA.N },
        { name: 'CA', pt: bA.CA },
        { name: 'C', pt: bA.C },
        { name: 'O', pt: bA.O }
      ];
      if (bA.CB) atomsA.push({ name: 'CB', pt: bA.CB });

      const atomsB: { name: string; pt: THREE.Vector3 }[] = [
        { name: 'N', pt: bB.N },
        { name: 'CA', pt: bB.CA },
        { name: 'C', pt: bB.C },
        { name: 'O', pt: bB.O }
      ];
      if (bB.CB) atomsB.push({ name: 'CB', pt: bB.CB });

      for (let u = 0; u < atomsA.length; u++) {
        const a1 = atomsA[u];
        for (let v = 0; v < atomsB.length; v++) {
          const a2 = atomsB[v];
          if (a1.name === 'CA' && a2.name === 'CA') continue;

          const d = a1.pt.distanceTo(a2.pt);
          const isON = (a1.name === 'O' && a2.name === 'N') || (a1.name === 'N' && a2.name === 'O');
          const isCBCB = a1.name === 'CB' && a2.name === 'CB';

          let minD = 2.80;
          if (isON && isHBondPair) minD = 2.58;
          else if (isON) minD = 2.68;
          else if (isCBCB && isDisulfidePair) minD = 3.35;
          else if (isCBCB) minD = 3.30;

          if (d < minD) {
            const ov = minD - d;
            eSteric += 32.0 * ov * ov + (d < 2.15 ? 260.0 * (2.15 - d) ** 2 : 0);
            if (d < 2.20) backboneStericOverlapsCount++;
          }
        }
      }
    }
  }
  eSteric *= weights.stericRepulsionWeight;

  const totalInternalPhysicalEnergyKcal =
    eRama + eHBond + eContact + eBeta + eHydroCore + eDisulfide + eElectro + eSteric;

  return {
    breakdown: {
      ramachandranEnergyKcal: parseFloat(eRama.toFixed(3)),
      backboneHBondEnergyKcal: parseFloat(eHBond.toFixed(3)),
      longRangeContactEnergyKcal: parseFloat(eContact.toFixed(3)),
      betaPairingEnergyKcal: parseFloat(eBeta.toFixed(3)),
      hydrophobicCoreCompactnessKcal: parseFloat(eHydroCore.toFixed(3)),
      disulfideClampEnergyKcal: parseFloat(eDisulfide.toFixed(3)),
      electrostaticEnergyKcal: parseFloat(eElectro.toFixed(3)),
      stericClashPenaltyKcal: parseFloat(eSteric.toFixed(3)),
      totalInternalPhysicalEnergyKcal: parseFloat(totalInternalPhysicalEnergyKcal.toFixed(3))
    },
    rgActual: parseFloat(rgActual.toFixed(3)),
    satisfiedContactsCount,
    satisfiedHBondsCount,
    backboneStericOverlapsCount
  };
}

// -----------------------------------------------------------------------------
// 9. POST-HOC NATIVE CONTACT PRECISION / RECALL (EVALUATION ONLY!)
// -----------------------------------------------------------------------------

function computePostHocNativeContactOverlap(
  modelCaPositions: [number, number, number][],
  parsedReference: ParsedPdbStructure | null
): {
  precisionPct: number | null;
  recallPct: number | null;
  f1Pct: number | null;
} {
  if (!parsedReference || parsedReference.residuesWithCa.length !== modelCaPositions.length) {
    return { precisionPct: null, recallPct: null, f1Pct: null };
  }
  const n = modelCaPositions.length;
  const CONTACT_CUTOFF = 8.0; // Standard CASP 8.0 Å contact threshold for |i - j| >= 6
  let nativeContacts = 0;
  let predictedContacts = 0;
  let truePositives = 0;

  for (let i = 0; i < n; i++) {
    const mA = modelCaPositions[i];
    const rA = parsedReference.residuesWithCa[i].caCoord!;
    for (let j = i + 6; j < n; j++) {
      const mB = modelCaPositions[j];
      const rB = parsedReference.residuesWithCa[j].caCoord!;
      const dMod = Math.hypot(mB[0] - mA[0], mB[1] - mA[1], mB[2] - mA[2]);
      const dRef = Math.hypot(rB[0] - rA[0], rB[1] - rA[1], rB[2] - rA[2]);
      const isMod = dMod <= CONTACT_CUTOFF;
      const isRef = dRef <= CONTACT_CUTOFF;
      if (isMod) predictedContacts++;
      if (isRef) nativeContacts++;
      if (isMod && isRef) truePositives++;
    }
  }

  if (nativeContacts === 0 && predictedContacts === 0) {
    return { precisionPct: 100.0, recallPct: 100.0, f1Pct: 100.0 };
  }
  const prec = predictedContacts > 0 ? (truePositives / predictedContacts) * 100.0 : 0.0;
  const rec = nativeContacts > 0 ? (truePositives / nativeContacts) * 100.0 : 0.0;
  const f1 = prec + rec > 1e-6 ? (2.0 * prec * rec) / (prec + rec) : 0.0;

  return {
    precisionPct: parseFloat(prec.toFixed(2)),
    recallPct: parseFloat(rec.toFixed(2)),
    f1Pct: parseFloat(f1.toFixed(2))
  };
}

function computeFastCaRmsdCentered(a: NeRFBackboneAtomSet[], b: THREE.Vector3[]): number {
  const n = Math.min(a.length, b.length);
  if (n === 0) return 0;
  // Upper-bound / fast Kabsch-free coordinate drift check in common internal frame
  // plus centroid alignment:
  let sq = 0;
  for (let i = 0; i < n; i++) {
    sq += a[i].CA.distanceToSquared(b[i]);
  }
  return Math.sqrt(sq / n);
}

// -----------------------------------------------------------------------------
// 10. FIVE-STAGE HIERARCHICAL FOLDING ENGINE (STAGES A -> B -> C -> D -> E)
// -----------------------------------------------------------------------------

const SEED_TOPOLOGY_HYPOTHESES = [
  { label: 'Compact-Globular-Cradle-A', turnGain: 1.00, strandTwist: 0.0, loopHandedness: 1 },
  { label: 'Antiparallel-Hairpin-Register-B', turnGain: 1.12, strandTwist: -6.0, loopHandedness: -1 },
  { label: 'Amphipathic-Helix-Sheet-Pack-C', turnGain: 0.92, strandTwist: 5.0, loopHandedness: 1 },
  { label: 'Alternative-Turn-Pivot-D', turnGain: 1.06, strandTwist: -3.5, loopHandedness: -1 },
  { label: 'Tight-Core-Fold-E', turnGain: 1.18, strandTwist: 3.0, loopHandedness: 1 },
  { label: 'Open-Clamp-Loop-F', turnGain: 0.88, strandTwist: 0.0, loopHandedness: -1 }
];

export function runPhysicallyConstrainedFolding(
  sequenceCodes: string[],
  options?: PhysicallyConstrainedFoldingOptions
): PhysicallyConstrainedFoldingResult {
  const tStart = performance.now();
  const generatorMode: FoldingGeneratorMode = options?.generatorMode ?? 'FULL_HIERARCHICAL_ENGINE';
  const paramVersion = options?.parameterVersion ?? `v2.0-${generatorMode.toLowerCase()}`;
  const iterationId = options?.iterationId ?? `hier-${Date.now()}`;
  const selectionPolicy: CandidateSelectionPolicy = options?.selectionPolicy ?? 'INTERNAL_PHYSICS_ONLY';

  const weights: PhysicalEnergyWeights = {
    ...DEFAULT_PHYSICAL_ENERGY_WEIGHTS,
    ...(options?.weights ?? {})
  };

  const n = sequenceCodes.length;
  const disulfidePairs = inferUniversalDisulfideTopology(sequenceCodes, options?.explicitDisulfidePairs);

  // Predict sequence-driven secondary structure topology
  const sseTopology = predictSecondaryStructureTopology(sequenceCodes, {
    cdAlphaTargetPct: options?.cdAlphaTargetPct,
    cdBetaTargetPct: options?.cdBetaTargetPct,
    turnBias: weights.turnNucleationBias,
    externalSS1Letter: options?.externalRestraints?.secondaryStructure1Letter
  });

  // Build explicit contact & H-bond constraint graph (ZERO access to groundTruthReference!)
  const restraintGraph = buildPhysicalRestraintGraph(
    sequenceCodes,
    generatorMode,
    sseTopology.segments,
    sseTopology.betaPairings,
    disulfidePairs,
    options?.externalRestraints
  );

  const numCandidates =
    generatorMode === 'FULL_HIERARCHICAL_ENGINE'
      ? Math.max(1, Math.min(6, options?.candidateCount ?? 4))
      : Math.max(1, Math.min(4, options?.candidateCount ?? 2));

  const stageCBudgetPerCandidate =
    options?.stageCMovesBudget ??
    (generatorMode === 'INDEPENDENT_TORSIONS_ONLY'
      ? Math.min(65, Math.max(28, n))
      : generatorMode === 'SECONDARY_STRUCTURE_CONSTRAINED'
      ? Math.min(75, Math.max(36, n))
      : Math.min(135, Math.max(60, Math.round(n * 1.35))));

  const l2Lattice = getPolarMicroLattice(2);
  let totalFunctionEvaluations = 0;

  // Parse reference ONCE strictly for Stage E post-hoc evaluation (NEVER passed to Stages A-D!)
  const parsedRefForPostHoc = options?.groundTruthReference?.pdbText
    ? parsePdbForGroundTruth(options.groundTruthReference.pdbText, {
        chainId: options.groundTruthReference.chainId
      })
    : null;

  const allCandidates: HierarchicalCandidateRecord[] = [];

  for (let candIdx = 0; candIdx < numCandidates; candIdx++) {
    const hyp = SEED_TOPOLOGY_HYPOTHESES[candIdx % SEED_TOPOLOGY_HYPOTHESES.length];
    const rng = new DeterministicRng(hashSequenceSeed(sequenceCodes, 1009 * (candIdx + 1)));
    const candidateId = `${iterationId}-cand-${candIdx + 1}`;
    const stageTraces: HierarchicalStageTrace[] = [];

    // =========================================================================
    // STAGE A: LOCAL CONFORMATIONS CONSISTENT WITH COVALENT CHAIN GEOMETRY
    // =========================================================================
    const phiDeg = new Array<number>(n).fill(-62.0);
    const psiDeg = new Array<number>(n).fill(-42.0);
    const omegaDeg = new Array<number>(n).fill(180.0);
    const flexByRes = new Array<number>(n).fill(0.5);

    for (let i = 0; i < n; i++) {
      const code = sequenceCodes[i] || 'ALA';
      const sec = sseTopology.perResidueTypes[i];
      const turnSub = sseTopology.perResidueTurnSubtypes[i];

      if (generatorMode === 'INDEPENDENT_TORSIONS_ONLY' || generatorMode === 'LONG_RANGE_CONTACTS_ONLY') {
        // Without SSE constraints: initialize from amino-acid intrinsic Ramachandran basin + seed variation
        const aaPref = (SNAPLOCK_21_BLOCKS[code] || SNAPLOCK_21_BLOCKS.ALA).physico.ramachandranPref;
        if (code === 'PRO') {
          phiDeg[i] = -63.0;
          psiDeg[i] = candIdx % 2 === 0 ? 142.0 : -36.0;
          flexByRes[i] = 0.15;
        } else if (code === 'GLY') {
          phiDeg[i] = hyp.loopHandedness * 72.0 + rng.nextRange(-18, 18);
          psiDeg[i] = hyp.loopHandedness * 15.0 + rng.nextRange(-22, 22);
          flexByRes[i] = 1.0;
        } else if (aaPref.includes('Beta')) {
          phiDeg[i] = -122.0 + rng.nextRange(-14, 14);
          psiDeg[i] = 132.0 + rng.nextRange(-14, 14);
          flexByRes[i] = 0.75;
        } else if (aaPref.includes('Zwroty') || aaPref.includes('Pętle')) {
          phiDeg[i] = -72.0 + rng.nextRange(-18, 18);
          psiDeg[i] = (candIdx % 2 === 0 ? 125.0 : -20.0) + rng.nextRange(-18, 18);
          flexByRes[i] = 0.92;
        } else {
          phiDeg[i] = -64.0 + rng.nextRange(-14, 14);
          psiDeg[i] = -41.0 + rng.nextRange(-14, 14);
          flexByRes[i] = 0.75;
        }
      } else {
        // Stage A + B foundation: SSE-aware local Ramachandran conformations
        if (code === 'PRO') {
          phiDeg[i] = -63.0;
          psiDeg[i] = sec === 'HELIX' ? -38.0 : 142.0;
          flexByRes[i] = 0.12;
        } else if (sec === 'HELIX') {
          phiDeg[i] = -58.5 + rng.nextRange(-2.5, 2.5);
          psiDeg[i] = -46.0 + rng.nextRange(-2.5, 2.5);
          flexByRes[i] = 0.22;
        } else if (sec === 'SHEET') {
          phiDeg[i] = -124.0 + hyp.strandTwist + rng.nextRange(-4.0, 4.0);
          psiDeg[i] = 131.0 - hyp.strandTwist + rng.nextRange(-4.0, 4.0);
          flexByRes[i] = 0.30;
        } else if (sec === 'TURN') {
          const isSecondTurnRes = i > 0 && sseTopology.perResidueTypes[i - 1] === 'TURN';
          if (turnSub === 'TYPE_II') {
            phiDeg[i] = isSecondTurnRes ? 80.0 : -60.0;
            psiDeg[i] = isSecondTurnRes ? 0.0 : 120.0;
          } else if (turnSub === 'TYPE_I_PRIME') {
            phiDeg[i] = isSecondTurnRes ? 90.0 : 60.0;
            psiDeg[i] = isSecondTurnRes ? 0.0 : 30.0;
          } else if (turnSub === 'TYPE_I') {
            phiDeg[i] = isSecondTurnRes ? -90.0 * hyp.turnGain : -60.0;
            psiDeg[i] = isSecondTurnRes ? 0.0 : -30.0 * hyp.turnGain;
          } else {
            phiDeg[i] = code === 'GLY' ? hyp.loopHandedness * 78.0 : -75.0 + rng.nextRange(-20, 20);
            psiDeg[i] = hyp.loopHandedness * 25.0 + rng.nextRange(-25, 25);
          }
          flexByRes[i] = 0.95;
        } else {
          // Flexible COIL / linker hinge
          phiDeg[i] = code === 'GLY' ? hyp.loopHandedness * 75.0 : -78.0 + rng.nextRange(-28, 28);
          psiDeg[i] = (i % 2 === 0 ? 132.0 : -28.0) + rng.nextRange(-24, 24);
          flexByRes[i] = 1.00;
        }
      }
    }

    let curBackbone = buildNeRFBackboneFromTorsions(sequenceCodes, phiDeg, psiDeg, omegaDeg);
    let curEval = evaluateInternalPhysicalEnergy(
      sequenceCodes,
      phiDeg,
      psiDeg,
      curBackbone,
      restraintGraph.longRangeContacts,
      restraintGraph.hydrogenBonds,
      disulfidePairs,
      weights,
      options?.rgTarget
    );
    totalFunctionEvaluations++;

    stageTraces.push({
      stageCode: 'A',
      stageName: 'Generowanie lokalnych konformacji NeRF (φ, ψ, ω)',
      evaluationsExecuted: 1,
      acceptedMoves: 1,
      internalPhysicalEnergyKcal: curEval.breakdown.totalInternalPhysicalEnergyKcal,
      rgAngstrom: curEval.rgActual,
      satisfiedContactsCount: curEval.satisfiedContactsCount,
      satisfiedHBondsCount: curEval.satisfiedHBondsCount,
      backboneStericOverlapsCount: curEval.backboneStericOverlapsCount,
      heavyAtomClashesUnder155A: 0,
      topologyDriftFromPrevStageRmsd: 0,
      notes: `Zbudowano kowalencyjny szkielet NeRF dla ${n} reszt (${hyp.label}).`
    });

    // =========================================================================
    // STAGE B: SECONDARY STRUCTURE ASSEMBLY & RELATIVE ORIENTATION
    // =========================================================================
    const prevStageACa = curBackbone.map(b => b.CA.clone());
    let stageBEvals = 0;
    let stageBAccepted = 0;

    if (
      generatorMode === 'SECONDARY_STRUCTURE_CONSTRAINED' ||
      generatorMode === 'FULL_HIERARCHICAL_ENGINE'
    ) {
      // Refine turn & SSE junction pivots to nucleate β-hairpins and α-helix H-bonds cleanly
      const junctionIndices: number[] = [];
      for (const seg of sseTopology.segments) {
        if (seg.type === 'TURN' || seg.type === 'COIL') {
          for (let k = seg.startIdx; k <= seg.endIdx; k++) {
            if (k > 0 && k < n - 1) junctionIndices.push(k);
          }
        } else {
          if (seg.startIdx > 0) junctionIndices.push(seg.startIdx);
          if (seg.endIdx < n - 1) junctionIndices.push(seg.endIdx);
        }
      }

      const stageBWeights: PhysicalEnergyWeights = {
        ...weights,
        hBondWeight: weights.hBondWeight * 1.45,
        betaPairingWeight: weights.betaPairingWeight * 1.35,
        longRangeContactWeight: generatorMode === 'SECONDARY_STRUCTURE_CONSTRAINED' ? 0.0 : weights.longRangeContactWeight * 0.65
      };

      curEval = evaluateInternalPhysicalEnergy(
        sequenceCodes,
        phiDeg,
        psiDeg,
        curBackbone,
        restraintGraph.longRangeContacts,
        restraintGraph.hydrogenBonds,
        disulfidePairs,
        stageBWeights,
        options?.rgTarget
      );
      stageBEvals++;

      const maxBSteps = Math.min(42, Math.max(16, junctionIndices.length * 2));
      for (let step = 0; step < maxBSteps && junctionIndices.length > 0; step++) {
        const pivotIdx = junctionIndices[step % junctionIndices.length];
        const code = sequenceCodes[pivotIdx];
        const savedPhi = phiDeg[pivotIdx];
        const savedPsi = psiDeg[pivotIdx];

        const deltaScale = step < maxBSteps / 2 ? 22.0 : 9.0;
        if (code !== 'PRO') {
          phiDeg[pivotIdx] = wrapAngleDeg(savedPhi + rng.nextRange(-deltaScale, deltaScale));
        }
        psiDeg[pivotIdx] = wrapAngleDeg(savedPsi + rng.nextRange(-deltaScale, deltaScale));

        const trialBb = buildNeRFBackboneFromTorsions(sequenceCodes, phiDeg, psiDeg, omegaDeg);
        const trialEv = evaluateInternalPhysicalEnergy(
          sequenceCodes,
          phiDeg,
          psiDeg,
          trialBb,
          restraintGraph.longRangeContacts,
          restraintGraph.hydrogenBonds,
          disulfidePairs,
          stageBWeights,
          options?.rgTarget
        );
        stageBEvals++;

        if (trialEv.breakdown.totalInternalPhysicalEnergyKcal < curEval.breakdown.totalInternalPhysicalEnergyKcal - 1e-4) {
          curBackbone = trialBb;
          curEval = trialEv;
          stageBAccepted++;
        } else {
          phiDeg[pivotIdx] = savedPhi;
          psiDeg[pivotIdx] = savedPsi;
        }
      }
      totalFunctionEvaluations += stageBEvals;
    }

    // Re-evaluate with full canonical weights at end of Stage B
    curEval = evaluateInternalPhysicalEnergy(
      sequenceCodes,
      phiDeg,
      psiDeg,
      curBackbone,
      restraintGraph.longRangeContacts,
      restraintGraph.hydrogenBonds,
      disulfidePairs,
      weights,
      options?.rgTarget
    );
    totalFunctionEvaluations++;

    stageTraces.push({
      stageCode: 'B',
      stageName: 'Budowa struktur drugorzędowych (SSE) i wzajemnego ułożenia',
      evaluationsExecuted: stageBEvals,
      acceptedMoves: stageBAccepted,
      internalPhysicalEnergyKcal: curEval.breakdown.totalInternalPhysicalEnergyKcal,
      rgAngstrom: curEval.rgActual,
      satisfiedContactsCount: curEval.satisfiedContactsCount,
      satisfiedHBondsCount: curEval.satisfiedHBondsCount,
      backboneStericOverlapsCount: curEval.backboneStericOverlapsCount,
      heavyAtomClashesUnder155A: 0,
      topologyDriftFromPrevStageRmsd: parseFloat(computeFastCaRmsdCentered(curBackbone, prevStageACa).toFixed(3)),
      notes: `Zsegmentowano ${sseTopology.segments.length} bloków SSE i ${sseTopology.betaPairings.length} parowań nici β.`
    });

    // =========================================================================
    // STAGE C: GLOBAL TOPOLOGY SEARCH WITH LONG-RANGE CONTACTS (SA + HINGE MOVES)
    // =========================================================================
    const prevStageBCa = curBackbone.map(b => b.CA.clone());
    let stageCEvals = 0;
    let stageCAccepted = 0;

    let bestStageCPhi = [...phiDeg];
    let bestStageCPsi = [...psiDeg];
    let bestStageCBackbone = curBackbone;
    let bestStageCEval = curEval;

    // Build weighted sampling pool favoring flexible hinges/turns/loops while allowing gentle SSE breathing
    const movableResidues: number[] = [];
    for (let i = 1; i < n - 1; i++) {
      const copies = flexByRes[i] >= 0.8 ? 4 : flexByRes[i] >= 0.5 ? 2 : 1;
      for (let c = 0; c < copies; c++) movableResidues.push(i);
    }

    const tInitial = 2.4;
    const tFinal = 0.08;

    for (let m = 0; m < stageCBudgetPerCandidate && movableResidues.length > 0; m++) {
      const progress = m / Math.max(1, stageCBudgetPerCandidate - 1);
      const temp = tInitial * Math.pow(tFinal / tInitial, progress);

      const idx1 = movableResidues[rng.nextInt(movableResidues.length)];
      // Compensatory crank/pair move on (idx1, idx2) 40% of the time so local loops can fold without whipping entire tail
      const useCompensatoryPair = rng.nextFloat() < 0.42 && idx1 + 2 < n - 1;
      const idx2 = useCompensatoryPair ? Math.min(n - 2, idx1 + 1 + rng.nextInt(3)) : -1;

      const sPhi1 = phiDeg[idx1], sPsi1 = psiDeg[idx1];
      const sPhi2 = idx2 >= 0 ? phiDeg[idx2] : 0;
      const sPsi2 = idx2 >= 0 ? psiDeg[idx2] : 0;

      const amp = (1.0 - 0.65 * progress) * 34.0 * flexByRes[idx1];
      const d1 = rng.nextRange(-amp, amp);
      const d2 = rng.nextRange(-amp, amp);

      if (sequenceCodes[idx1] !== 'PRO') {
        phiDeg[idx1] = wrapAngleDeg(sPhi1 + d1);
      }
      psiDeg[idx1] = wrapAngleDeg(sPsi1 + d2);

      if (idx2 >= 0) {
        const compScale = 0.65 * flexByRes[idx2];
        if (sequenceCodes[idx2] !== 'PRO') {
          phiDeg[idx2] = wrapAngleDeg(sPhi2 - d2 * compScale);
        }
        psiDeg[idx2] = wrapAngleDeg(sPsi2 - d1 * compScale);
      }

      const trialBb = buildNeRFBackboneFromTorsions(sequenceCodes, phiDeg, psiDeg, omegaDeg);
      const trialEv = evaluateInternalPhysicalEnergy(
        sequenceCodes,
        phiDeg,
        psiDeg,
        trialBb,
        restraintGraph.longRangeContacts,
        restraintGraph.hydrogenBonds,
        disulfidePairs,
        weights,
        options?.rgTarget
      );
      stageCEvals++;

      const dE = trialEv.breakdown.totalInternalPhysicalEnergyKcal - curEval.breakdown.totalInternalPhysicalEnergyKcal;
      const acceptMetropolis = dE <= 0 || rng.nextFloat() < Math.exp(-dE / Math.max(0.02, temp));

      if (acceptMetropolis) {
        curBackbone = trialBb;
        curEval = trialEv;
        stageCAccepted++;

        if (
          trialEv.breakdown.totalInternalPhysicalEnergyKcal <
          bestStageCEval.breakdown.totalInternalPhysicalEnergyKcal - 1e-4
        ) {
          bestStageCPhi = [...phiDeg];
          bestStageCPsi = [...psiDeg];
          bestStageCBackbone = trialBb;
          bestStageCEval = trialEv;
        }
      } else {
        phiDeg[idx1] = sPhi1;
        psiDeg[idx1] = sPsi1;
        if (idx2 >= 0) {
          phiDeg[idx2] = sPhi2;
          psiDeg[idx2] = sPsi2;
        }
      }
    }

    // Restore best global topology discovered in Stage C
    for (let i = 0; i < n; i++) {
      phiDeg[i] = bestStageCPhi[i];
      psiDeg[i] = bestStageCPsi[i];
    }
    curBackbone = bestStageCBackbone;
    curEval = bestStageCEval;
    totalFunctionEvaluations += stageCEvals;

    stageTraces.push({
      stageCode: 'C',
      stageName: 'Wyszukiwanie globalnej topologii z kontaktami dalekiego zasięgu',
      evaluationsExecuted: stageCEvals,
      acceptedMoves: stageCAccepted,
      internalPhysicalEnergyKcal: curEval.breakdown.totalInternalPhysicalEnergyKcal,
      rgAngstrom: curEval.rgActual,
      satisfiedContactsCount: curEval.satisfiedContactsCount,
      satisfiedHBondsCount: curEval.satisfiedHBondsCount,
      backboneStericOverlapsCount: curEval.backboneStericOverlapsCount,
      heavyAtomClashesUnder155A: 0,
      topologyDriftFromPrevStageRmsd: parseFloat(computeFastCaRmsdCentered(curBackbone, prevStageBCa).toFixed(3)),
      notes: `Przeszukano globalną topologię torsyjną: spełniono ${curEval.satisfiedContactsCount}/${restraintGraph.longRangeContacts.length} kontaktów dalekiego zasięgu.`
    });

    // =========================================================================
    // STAGE D: TOPOLOGY-PRESERVING ATOMIC RELAXATION & STERIC CLASH REMOVAL
    // =========================================================================
    const stageCTopologyAnchorCa = curBackbone.map(b => b.CA.clone());
    let stageDEvals = 0;
    let stageDAccepted = 0;

    const stageDWeights: PhysicalEnergyWeights = {
      ...weights,
      stericRepulsionWeight: weights.stericRepulsionWeight * 3.5,
      disulfideClampWeight: weights.disulfideClampWeight * 1.8
    };

    curEval = evaluateInternalPhysicalEnergy(
      sequenceCodes,
      phiDeg,
      psiDeg,
      curBackbone,
      restraintGraph.longRangeContacts,
      restraintGraph.hydrogenBonds,
      disulfidePairs,
      stageDWeights,
      options?.rgTarget
    );
    stageDEvals++;

    // Identify Stage D hotspot residues (close steric contacts < 2.95 Å, disulfide hinges, turns/coils, Ramachandran outliers)
    const collectStageDHotspots = (bb: NeRFBackboneAtomSet[]): number[] => {
      const hot = new Set<number>();
      for (let i = 1; i < n - 1; i++) {
        if (sseTopology.perResidueTypes[i] === 'TURN' || sseTopology.perResidueTypes[i] === 'COIL') {
          hot.add(i);
        } else if (!evaluateRamachandranState(sequenceCodes[i], phiDeg[i], psiDeg[i]).favored) {
          hot.add(i);
        }
      }
      for (const [s1, s2] of disulfidePairs) {
        for (let k = Math.max(1, s1 - 2); k <= Math.min(n - 2, s1 + 3); k++) hot.add(k);
        for (let k = Math.max(1, s2 - 3); k <= Math.min(n - 2, s2 + 2); k++) hot.add(k);
      }
      for (let i = 0; i < n; i++) {
        for (let j = i + 2; j < n; j++) {
          if (bb[i].CA.distanceToSquared(bb[j].CA) > 49.0) continue;
          const dMin = Math.min(
            bb[i].CA.distanceTo(bb[j].CA) - 0.9,
            bb[i].O.distanceTo(bb[j].O),
            bb[i].O.distanceTo(bb[j].C),
            bb[i].CB && bb[j].CB ? bb[i].CB!.distanceTo(bb[j].CB!) - 0.4 : 9.9
          );
          if (dMin < 2.90) {
            for (let k = Math.max(1, i - 1); k <= Math.min(n - 2, i + 1); k++) hot.add(k);
            for (let k = Math.max(1, j - 1); k <= Math.min(n - 2, j + 1); k++) hot.add(k);
            const mid = Math.max(1, Math.min(n - 2, Math.floor((i + j) / 2)));
            hot.add(mid);
          }
        }
      }
      return Array.from(hot).sort((a, b) => a - b);
    };

    // Local fine torsion relaxation with strict topology tether (drift <= 1.45 Å from Stage C scaffold)
    const maxDSweeps = curEval.backboneStericOverlapsCount > 0 || disulfidePairs.length > 0 ? 2 : 1;
    for (let sweep = 0; sweep < maxDSweeps; sweep++) {
      const stepDeg = sweep === 0 ? 6.0 : 2.5;
      const hotspots = collectStageDHotspots(curBackbone);
      for (const i of hotspots) {
        for (const isPhi of [true, false]) {
          if (isPhi && sequenceCodes[i] === 'PRO') continue;
          const orig = isPhi ? phiDeg[i] : psiDeg[i];

          for (const sign of [1, -1] as const) {
            if (isPhi) phiDeg[i] = wrapAngleDeg(orig + sign * stepDeg);
            else psiDeg[i] = wrapAngleDeg(orig + sign * stepDeg);

            const trialBb = buildNeRFBackboneFromTorsions(sequenceCodes, phiDeg, psiDeg, omegaDeg);
            const drift = computeFastCaRmsdCentered(trialBb, stageCTopologyAnchorCa);
            if (drift > 1.45) {
              if (isPhi) phiDeg[i] = orig;
              else psiDeg[i] = orig;
              continue;
            }

            const trialEv = evaluateInternalPhysicalEnergy(
              sequenceCodes,
              phiDeg,
              psiDeg,
              trialBb,
              restraintGraph.longRangeContacts,
              restraintGraph.hydrogenBonds,
              disulfidePairs,
              stageDWeights,
              options?.rgTarget
            );
            stageDEvals++;

            const currObj = curEval.breakdown.totalInternalPhysicalEnergyKcal;
            const trialObj = trialEv.breakdown.totalInternalPhysicalEnergyKcal + 0.35 * drift * drift;

            if (trialObj < currObj - 1e-3) {
              curBackbone = trialBb;
              curEval = trialEv;
              stageDAccepted++;
              break;
            } else {
              if (isPhi) phiDeg[i] = orig;
              else psiDeg[i] = orig;
            }
          }
        }
      }
      if (curEval.backboneStericOverlapsCount === 0 && disulfidePairs.length === 0) break;
    }

    // Build all-atom sidechains with continuous off-rotamer χ steric clash relief
    const buildAllAtomFromCurrentBackbone = (bb: NeRFBackboneAtomSet[]) => {
      const caPos: [number, number, number][] = bb.map(r => [
        parseFloat(r.CA.x.toFixed(3)),
        parseFloat(r.CA.y.toFixed(3)),
        parseFloat(r.CA.z.toFixed(3))
      ]);
      const explicitBb: ExplicitResidueBackboneCoords[] = bb.map(r => ({
        N: [parseFloat(r.N.x.toFixed(3)), parseFloat(r.N.y.toFixed(3)), parseFloat(r.N.z.toFixed(3))],
        CA: [parseFloat(r.CA.x.toFixed(3)), parseFloat(r.CA.y.toFixed(3)), parseFloat(r.CA.z.toFixed(3))],
        C: [parseFloat(r.C.x.toFixed(3)), parseFloat(r.C.y.toFixed(3)), parseFloat(r.C.z.toFixed(3))],
        O: [parseFloat(r.O.x.toFixed(3)), parseFloat(r.O.y.toFixed(3)), parseFloat(r.O.z.toFixed(3))],
        OXT: r.OXT
          ? [parseFloat(r.OXT.x.toFixed(3)), parseFloat(r.OXT.y.toFixed(3)), parseFloat(r.OXT.z.toFixed(3))]
          : undefined
      }));
      const cbPos: ([number, number, number] | null)[] = bb.map(r =>
        r.CB ? [parseFloat(r.CB.x.toFixed(3)), parseFloat(r.CB.y.toFixed(3)), parseFloat(r.CB.z.toFixed(3))] : null
      );
      const { sgPositions } = reconstructSidechainCbAndSg(sequenceCodes, caPos, disulfidePairs);
      const allAtom = buildAllAtomStructure({
        sequenceCodes,
        caPositions: caPos,
        cbOverrides: cbPos,
        sgOverrides: sgPositions,
        ssBondPairs: disulfidePairs,
        rotamerOverrides: options?.rotamerOverrides,
        secondaryTypes: sseTopology.perResidueTypes,
        explicitBackboneByRes: explicitBb,
        enableContinuousChiClashRelief: true
      });
      return { caPos, explicitBb, cbPos, sgPositions, allAtom };
    };

    let assembled = buildAllAtomFromCurrentBackbone(curBackbone);

    // Check if any severe heavy-atom pair (|resA - resB| >= 2, non-disulfide) has d < 1.62 Å; if so, nudge hinge between them!
    const findSevereClashingResiduePairs = (allAtom: AllAtomStructureResult): [number, number][] => {
      const clashingPairs: [number, number][] = [];
      const seen = new Set<string>();
      const atoms = allAtom.atoms;
      for (let a = 0; a < atoms.length; a++) {
        const at1 = atoms[a];
        for (let b = a + 1; b < atoms.length; b++) {
          const at2 = atoms[b];
          if (Math.abs(at1.resSeq - at2.resSeq) < 2) continue;
          const isSs =
            (at1.atomName === 'SG' || at1.atomName === 'SE') &&
            (at2.atomName === 'SG' || at2.atomName === 'SE');
          if (isSs) continue;
          const d2 = (at1.x - at2.x) ** 2 + (at1.y - at2.y) ** 2 + (at1.z - at2.z) ** 2;
          if (d2 < 1.62 * 1.62) {
            const r1 = Math.min(at1.resSeq - 1, at2.resSeq - 1);
            const r2 = Math.max(at1.resSeq - 1, at2.resSeq - 1);
            const k = `${r1}:${r2}`;
            if (!seen.has(k)) {
              seen.add(k);
              clashingPairs.push([r1, r2]);
            }
          }
        }
      }
      return clashingPairs;
    };

    let severePairs = findSevereClashingResiduePairs(assembled.allAtom);
    for (let reliefIter = 0; reliefIter < 6 && severePairs.length > 0; reliefIter++) {
      const [r1, r2] = severePairs[0];
      const midHinge = Math.max(1, Math.min(n - 2, Math.floor((r1 + r2) / 2)));
      let resolved = false;
      for (const delta of [5.0, -5.0, 11.0, -11.0, 18.0, -18.0]) {
        const savedPsi = psiDeg[midHinge];
        psiDeg[midHinge] = wrapAngleDeg(savedPsi + delta);
        const trialBb = buildNeRFBackboneFromTorsions(sequenceCodes, phiDeg, psiDeg, omegaDeg);
        const trialAssembled = buildAllAtomFromCurrentBackbone(trialBb);
        stageDEvals++;
        const nextSevere = findSevereClashingResiduePairs(trialAssembled.allAtom);
        if (nextSevere.length < severePairs.length) {
          curBackbone = trialBb;
          assembled = trialAssembled;
          severePairs = nextSevere;
          stageDAccepted++;
          resolved = true;
          break;
        } else {
          psiDeg[midHinge] = savedPsi;
        }
      }
      if (!resolved) break;
    }

    // Recompute final canonical internal physical energy (with standard unboosted weights)
    const finalEnergyEval = evaluateInternalPhysicalEnergy(
      sequenceCodes,
      phiDeg,
      psiDeg,
      curBackbone,
      restraintGraph.longRangeContacts,
      restraintGraph.hydrogenBonds,
      disulfidePairs,
      weights,
      options?.rgTarget
    );
    totalFunctionEvaluations += stageDEvals + 1;

    const stageDDrift = parseFloat(computeFastCaRmsdCentered(curBackbone, stageCTopologyAnchorCa).toFixed(3));
    stageTraces.push({
      stageCode: 'D',
      stageName: 'Relaksacja atomowa i usuwanie kolizji z ochroną topologii',
      evaluationsExecuted: stageDEvals,
      acceptedMoves: stageDAccepted,
      internalPhysicalEnergyKcal: finalEnergyEval.breakdown.totalInternalPhysicalEnergyKcal,
      rgAngstrom: finalEnergyEval.rgActual,
      satisfiedContactsCount: finalEnergyEval.satisfiedContactsCount,
      satisfiedHBondsCount: finalEnergyEval.satisfiedHBondsCount,
      backboneStericOverlapsCount: finalEnergyEval.backboneStericOverlapsCount,
      heavyAtomClashesUnder155A: severePairs.length,
      topologyDriftFromPrevStageRmsd: stageDDrift,
      notes: `Dryf Cα względem Etapu C = ${stageDDrift} Å; ciężkie kolizje (<1.62 Å) po relaksacji rotamerów χ: ${severePairs.length}.`
    });

    // Build detailed per-residue torsion & E8/L2 diagnostics
    const residueTorsions: ResidueTorsionState[] = [];
    let ramaFavoredCount = 0;
    let ramaOutliersCount = 0;
    let l2ErrSum = 0;

    for (let i = 0; i < n; i++) {
      const code = sequenceCodes[i] || 'ALA';
      const rama = evaluateRamachandranState(code, phiDeg[i], psiDeg[i]);
      if (rama.favored) ramaFavoredCount++;
      else ramaOutliersCount++;

      // Map outgoing virtual bond CA(i) -> CA(i+1) onto discrete E8 (248) and L2 (2232) spherical lattices
      let e8PatchId = 248;
      let l2ErrDeg = 0;
      if (i < n - 1) {
        const bondVec = curBackbone[i + 1].CA.clone().sub(curBackbone[i].CA);
        const uBond = bondVec.clone().normalize();
        let bestDot = -Infinity;
        for (const ep of E8_PATCHES) {
          const dot = uBond.x * ep.vector[0] + uBond.y * ep.vector[1] + uBond.z * ep.vector[2];
          if (dot > bestDot) {
            bestDot = dot;
            e8PatchId = ep.id;
          }
        }
        const inAx = i > 0 ? curBackbone[i].CA.clone().sub(curBackbone[i - 1].CA) : new THREE.Vector3(0, 1, 0);
        const snapL2 = snapVectorToDiscreteMicroState(bondVec, inAx, l2Lattice);
        l2ErrDeg = snapL2.angleErrorDeg;
        l2ErrSum += l2ErrDeg;
      }

      residueTorsions.push({
        resIndex: i,
        resSeq: i + 1,
        resName: code,
        code1: THREE_TO_ONE[code] || 'X',
        phiDeg: parseFloat(phiDeg[i].toFixed(2)),
        psiDeg: parseFloat(psiDeg[i].toFixed(2)),
        omegaDeg: parseFloat(omegaDeg[i].toFixed(2)),
        secondaryType: sseTopology.perResidueTypes[i],
        turnSubtype: sseTopology.perResidueTurnSubtypes[i],
        ramachandranBasin: rama.basin,
        ramachandranFavored: rama.favored,
        ramachandranEnergyKcal: rama.energyKcal,
        kinematicFlexibility: parseFloat(flexByRes[i].toFixed(2)),
        e8NearestPatchId: e8PatchId,
        l2MicroFacetErrorDeg: parseFloat(l2ErrDeg.toFixed(3))
      });
    }

    const meanL2Err = n > 1 ? parseFloat((l2ErrSum / (n - 1)).toFixed(3)) : 0;

    // =========================================================================
    // STAGE E: INDEPENDENT GROUND-TRUTH VALIDATION OF CANDIDATE
    // =========================================================================
    const rawCandidatePdb = exportStandardAllAtomPdb(
      `SNAPLOCK PHASE2 ${candidateId} (${paramVersion})`,
      sequenceCodes,
      assembled.allAtom,
      sseTopology.perResidueTypes,
      {
        pipelineMode: generatorMode,
        meanMicroFacetErrorDeg: meanL2Err,
        rgActual: finalEnergyEval.rgActual,
        cohesiveEnergy: finalEnergyEval.breakdown.totalInternalPhysicalEnergyKcal
      }
    );

    let recordedPath: string | undefined;
    if (options?.artifactSink) {
      const p = options.artifactSink(candidateId, rawCandidatePdb);
      if (typeof p === 'string' && p.length > 0) recordedPath = p;
    }

    const gtRef = options?.groundTruthReference;
    const gtValidation = validateAgainstGroundTruth(
      rawCandidatePdb,
      gtRef?.pdbText ?? null,
      {
        referenceId: gtRef?.resolvedId,
        referenceName: gtRef?.referenceName ?? gtRef?.resolvedId,
        referenceChainId: gtRef?.chainId,
        referenceSource: gtRef?.source,
        allowSubsequenceAlignment: gtRef?.allowSubsequenceAlignment,
        iterationId: candidateId,
        parameterVersion: paramVersion,
        modelArtifactPath: recordedPath,
        criteria: options?.criteria
      }
    );

    const rmsdComparison = gtRef?.pdbText
      ? buildRmsdComparisonFromGroundTruth(gtValidation, gtRef.pdbText)
      : null;

    const finalPdbWithVerifiedHeader = exportStandardAllAtomPdb(
      `SNAPLOCK PHASE2 ${candidateId} (${paramVersion})`,
      sequenceCodes,
      assembled.allAtom,
      sseTopology.perResidueTypes,
      {
        pipelineMode: generatorMode,
        meanMicroFacetErrorDeg: meanL2Err,
        rgActual: finalEnergyEval.rgActual,
        cohesiveEnergy: finalEnergyEval.breakdown.totalInternalPhysicalEnergyKcal,
        groundTruthValidation: gtValidation,
        rmsdComparison
      }
    );

    // Compute maximum Sγ-Sγ disulfide error
    let maxSgErr = 0;
    for (const [i1, i2] of disulfidePairs) {
      const sg1 = assembled.sgPositions[i1];
      const sg2 = assembled.sgPositions[i2];
      if (sg1 && sg2) {
        const dSg = Math.hypot(sg1[0] - sg2[0], sg1[1] - sg2[1], sg1[2] - sg2[2]);
        maxSgErr = Math.max(maxSgErr, Math.abs(dSg - TARGET_SG_SG_DIST));
      }
    }

    const nativeOverlap = computePostHocNativeContactOverlap(assembled.caPos, parsedRefForPostHoc);

    stageTraces.push({
      stageCode: 'E',
      stageName: 'Niezależna walidacja kandydata (ground_truth_validation)',
      evaluationsExecuted: 1,
      acceptedMoves: gtValidation.modelGeometry.geometryValid ? 1 : 0,
      internalPhysicalEnergyKcal: finalEnergyEval.breakdown.totalInternalPhysicalEnergyKcal,
      rgAngstrom: finalEnergyEval.rgActual,
      satisfiedContactsCount: finalEnergyEval.satisfiedContactsCount,
      satisfiedHBondsCount: finalEnergyEval.satisfiedHBondsCount,
      backboneStericOverlapsCount: finalEnergyEval.backboneStericOverlapsCount,
      heavyAtomClashesUnder155A: gtValidation.modelGeometry.heavyAtomStericClashesCount,
      topologyDriftFromPrevStageRmsd: 0,
      notes: gtValidation.hasGroundTruthReference
        ? `Cα RMSD = ${gtValidation.caRmsdAngstrom?.toFixed(4) ?? 'ERR'} Å, TM-score = ${gtValidation.tmScore?.toFixed(4) ?? 'N/A'}, geometria fizyczna: ${gtValidation.modelGeometry.geometryValid ? 'OK' : 'NARUSZENIE'}.`
        : `Tryb de-novo bez referencji: geometria fizyczna ${gtValidation.modelGeometry.geometryValid ? 'OK' : 'NARUSZENIE'}.`
    });

    const totalLrc = restraintGraph.longRangeContacts.length;
    const evaluation: CandidateEvaluationMetrics = {
      groundTruth: {
        hasReference: gtValidation.hasGroundTruthReference,
        caRmsdAngstrom: gtValidation.caRmsdAngstrom,
        independentCaRmsdAngstrom: gtValidation.independentCaRmsdAngstrom,
        backboneRmsdAngstrom: gtValidation.backboneRmsdAngstrom,
        allHeavyAtomRmsdAngstrom: gtValidation.allHeavyAtomRmsdAngstrom,
        tmScore: gtValidation.tmScore,
        gdtTsPct: gtValidation.gdtTsPct,
        pairedResidues: gtValidation.pairedResidueCount,
        coveragePct: gtValidation.modelCoveragePct
      },
      localGeometry: {
        geometryValid: gtValidation.modelGeometry.geometryValid,
        covalentBondViolationsCount: gtValidation.modelGeometry.caBondViolations.length,
        meanCaCaStepAngstrom: gtValidation.modelGeometry.meanCaCaDistAngstrom,
        minCaCaStepAngstrom: gtValidation.modelGeometry.minCaCaDistAngstrom,
        maxCaCaStepAngstrom: gtValidation.modelGeometry.maxCaCaDistAngstrom,
        ramachandranFavoredPct: parseFloat(((ramaFavoredCount / Math.max(1, n)) * 100.0).toFixed(1)),
        ramachandranOutliersCount: ramaOutliersCount,
        meanL2MicroFacetErrorDeg: meanL2Err
      },
      sterics: {
        nonBondedCaClashesCount: gtValidation.modelGeometry.nonBondedCaClashes.length,
        heavyAtomClashesUnder155A: gtValidation.modelGeometry.heavyAtomStericClashesCount,
        softHeavyAtomContactsUnder215A: assembled.allAtom.stericClashesCount,
        collapsedChainDetected: gtValidation.modelGeometry.collapsedChainDetected
      },
      contacts: {
        restraintProvenance: restraintGraph.restraintProvenance,
        totalTargetLongRangeContacts: totalLrc,
        satisfiedLongRangeContacts: finalEnergyEval.satisfiedContactsCount,
        longRangeSatisfactionPct:
          totalLrc > 0
            ? parseFloat(((finalEnergyEval.satisfiedContactsCount / totalLrc) * 100.0).toFixed(1))
            : 100.0,
        totalTargetHBonds: restraintGraph.hydrogenBonds.length,
        satisfiedHBonds: finalEnergyEval.satisfiedHBondsCount,
        disulfideBridgesCount: disulfidePairs.length,
        maxDisulfideSgErrorAngstrom: parseFloat(maxSgErr.toFixed(3)),
        postHocNativeContactPrecisionPct: nativeOverlap.precisionPct,
        postHocNativeContactRecallPct: nativeOverlap.recallPct,
        postHocNativeContactF1Pct: nativeOverlap.f1Pct
      },
      internalModelPhysics: finalEnergyEval.breakdown
    };

    allCandidates.push({
      candidateIndex: candIdx + 1,
      candidateId,
      topologySeedLabel: hyp.label,
      parameterVersion: paramVersion,
      selectedAsBest: false,
      selectionReason: '',
      residueTorsions,
      caPositions: assembled.caPos,
      explicitBackbone: assembled.explicitBb,
      cbPositions: assembled.cbPos,
      sgPositions: assembled.sgPositions,
      allAtomStructure: assembled.allAtom,
      pdbText: finalPdbWithVerifiedHeader,
      modelPdbHash: gtValidation.modelPdbHash,
      stageTraces,
      evaluation,
      groundTruthValidation: gtValidation,
      rmsdComparison
    });
  }

  // ===========================================================================
  // CANDIDATE SELECTION (STRICT PHYSICAL GATE + TRANSPARENT POLICY)
  // ===========================================================================
  let bestIdx = 0;
  for (let i = 1; i < allCandidates.length; i++) {
    const curBest = allCandidates[bestIdx];
    const cand = allCandidates[i];

    const bestGeo = curBest.evaluation.localGeometry.geometryValid;
    const candGeo = cand.evaluation.localGeometry.geometryValid;

    // Physical validity ALWAYS dominates candidate ranking!
    if (candGeo && !bestGeo) {
      bestIdx = i;
      continue;
    }
    if (!candGeo && bestGeo) {
      continue;
    }

    if (
      selectionPolicy === 'CALIBRATION_TRAIN_SUPERVISED' &&
      cand.evaluation.groundTruth.caRmsdAngstrom !== null &&
      curBest.evaluation.groundTruth.caRmsdAngstrom !== null
    ) {
      if (cand.evaluation.groundTruth.caRmsdAngstrom < curBest.evaluation.groundTruth.caRmsdAngstrom - 1e-4) {
        bestIdx = i;
      }
    } else {
      // STRICT DE-NOVO / HELD-OUT TEST SELECTION: rank strictly by internal physical energy (ZERO reference leakage!)
      if (
        cand.evaluation.internalModelPhysics.totalInternalPhysicalEnergyKcal <
        curBest.evaluation.internalModelPhysics.totalInternalPhysicalEnergyKcal - 1e-4
      ) {
        bestIdx = i;
      }
    }
  }

  for (let i = 0; i < allCandidates.length; i++) {
    const c = allCandidates[i];
    if (i === bestIdx) {
      c.selectedAsBest = true;
      c.selectionReason =
        selectionPolicy === 'CALIBRATION_TRAIN_SUPERVISED'
          ? `Wybrano kandydata #${c.candidateIndex} (tryb jawnie nadzorowany TRAIN: geometria fizyczna=${c.evaluation.localGeometry.geometryValid ? 'OK' : 'WARN'}, Cα RMSD=${c.evaluation.groundTruth.caRmsdAngstrom?.toFixed(4) ?? 'N/A'} Å).`
          : `Wybrano kandydata #${c.candidateIndex} wyłącznie na podstawie wiarygodności geometrii fizycznej (${c.evaluation.localGeometry.geometryValid ? '0 naruszeń' : 'ostrzeżenie'}) i najniższej wewnętrznej energii fizycznej E_phys = ${c.evaluation.internalModelPhysics.totalInternalPhysicalEnergyKcal.toFixed(2)} kcal/mol (ZERO dostępu do referencji Ground-Truth przy wyborze modelu).`;
    } else {
      c.selectedAsBest = false;
      c.selectionReason = !c.evaluation.localGeometry.geometryValid
        ? `Odrzucono kandydata #${c.candidateIndex}: naruszenie geometrii fizycznej (${c.groundTruthValidation.modelGeometry.reasons.join(' ')}).`
        : `Kandydat #${c.candidateIndex} poprawny geometrycznie, lecz wyższa energia fizyczna E_phys = ${c.evaluation.internalModelPhysics.totalInternalPhysicalEnergyKcal.toFixed(2)} kcal/mol.`;
    }
  }

  const computeTimeMs = parseFloat((performance.now() - tStart).toFixed(2));

  return {
    generatorMode,
    restraintProvenance: restraintGraph.restraintProvenance,
    usedGroundTruthReferenceDuringGeneration: false,
    usedGroundTruthReferenceForCandidateSelection:
      selectionPolicy === 'CALIBRATION_TRAIN_SUPERVISED' && Boolean(options?.groundTruthReference?.pdbText),
    deNovoLimitationNotice: restraintGraph.deNovoLimitationNotice,
    parameterVersion: paramVersion,
    iterationId,
    computeTimeMs,
    totalFunctionEvaluations,
    secondarySegments: sseTopology.segments,
    betaSheetPairings: sseTopology.betaPairings,
    longRangeContacts: restraintGraph.longRangeContacts,
    hydrogenBondConstraints: restraintGraph.hydrogenBonds,
    disulfidePairs,
    selectedCandidate: allCandidates[bestIdx],
    allCandidates,
    architecturalLimitationsDiagnosed: diagnoseLegacyGeneratorLimitations()
  };
}
