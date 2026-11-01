// Closed-Loop Ground-Truth Validation & Derivative-Free Autocalibration Optimizer for SNAPLOCK
// ---------------------------------------------------------------------------------------------
// Architectural Note on Differentiability:
// The SNAPLOCK folding pipeline maps continuous torsion trajectories onto discrete spherical
// E8 (248 vectors) and L2 (2232 micro-facets) Voronoi cells via `argmax_k (v · n_k)`.
// Consequently, the end-to-end mapping θ -> X_model -> RMSD(X_model, X_ref) has zero gradient
// almost everywhere and jump discontinuities across facet boundaries (non-differentiable).
// Therefore, we NEVER claim RMSD is an end-to-end differentiable loss here; instead we use a
// rigorous derivative-free directional pattern search (Compass / Coordinate Pattern Search
// with adaptive step halving) coupled with strict parameter checkpointing and rollback.

import {
  DEFAULT_HIGH_ACCURACY_CRITERIA,
  GroundTruthIterationStatus,
  GroundTruthPassCriteria,
  GroundTruthValidationReport,
  kabschSuperposeSvd,
  parsePdbForGroundTruth
} from './ground_truth_validation';
import {
  BoundaryConditionsConfig,
  DEFAULT_SNAPLOCK_ENGINE_PARAMETERS,
  ExternalReferencePdbInput,
  FoldingMetrics,
  SnaplockEngineParameters,
  foldWithFastaAndBoundaryConstraints
} from './snaplockFastaAndPdb';
import {
  AlphaFoldRmsdComparisonResult
} from './snaplockAllAtomBuilder';
import {
  EXPERIMENTAL_GROUND_TRUTH_PDBS,
  resolveExperimentalReferenceById
} from './snaplockExperimentalReferences';

export type CalibrationDatasetSplit = 'TRAIN' | 'VALIDATION' | 'TEST' | 'SINGLE_TARGET';

export type CalibrationTerminalStatus =
  | 'SKALIBROWANY'
  | 'OSIAGNIETO_PROG_POSREDNI'
  | 'BRAK_ZBIEZNOSCI'
  | 'BLAD_WALIDACJI'
  | 'BRAK_REFERENCJI';

export interface BenchmarkProteinTarget {
  targetId: string;
  proteinName: string;
  pdbId: string;
  chainId: string;
  split: CalibrationDatasetSplit;
  sequenceCodes: string[];
  boundaryConfig: BoundaryConditionsConfig;
  referencePdbText: string;
  referenceSource: 'RCSB_PDB' | 'ALPHAFOLD_DB' | 'CUSTOM_PDB';
  criteria: GroundTruthPassCriteria;
}

export interface TunableParameterSpec {
  key: Exclude<keyof SnaplockEngineParameters, 'parameterVersion'>;
  label: string;
  initialStep: number;
  minStep: number;
  minBound: number;
  maxBound: number;
}

export const SNAPLOCK_TUNABLE_PARAMETER_SPECS: TunableParameterSpec[] = [
  { key: 'phiOffsetDeg', label: 'Kompensacja kąta φ [°]', initialStep: 8.0, minStep: 0.5, minBound: -45.0, maxBound: 45.0 },
  { key: 'psiOffsetDeg', label: 'Kompensacja kąta ψ [°]', initialStep: 8.0, minStep: 0.5, minBound: -45.0, maxBound: 45.0 },
  { key: 'helicalPitchTwistDeg', label: 'Skręt skoku helisy α [°]', initialStep: 6.0, minStep: 0.5, minBound: 72.0, maxBound: 132.0 },
  { key: 'betaPleatTwistDeg', label: 'Skręt harmonijki β [°]', initialStep: 6.0, minStep: 0.5, minBound: 130.0, maxBound: 195.0 },
  { key: 'turnBendGain', label: 'Wzmocnienie zwrotów pętli', initialStep: 0.18, minStep: 0.02, minBound: 0.45, maxBound: 2.10 },
  { key: 'hydrophobicCoreWeight', label: 'Waga rdzenia hydrofobowego', initialStep: 0.25, minStep: 0.03, minBound: 0.20, maxBound: 2.80 },
  { key: 'hBondWeight', label: 'Waga mostków wodorowych', initialStep: 0.25, minStep: 0.03, minBound: 0.20, maxBound: 2.80 },
  { key: 'dMaxRelaxation', label: 'Relaksacja obwiedni SAXS Dmax', initialStep: 0.10, minStep: 0.015, minBound: 0.65, maxBound: 1.45 },
  { key: 'ssRestraintWeight', label: 'Waga klamer disiarczkowych S-S', initialStep: 1.5, minStep: 0.25, minBound: 0.5, maxBound: 18.0 }
];

export interface AutocalibrationIterationRecord {
  trialNumber: number;
  iterationId: string;
  parameterVersion: string;
  tunedParameterKey: string;
  stepDeltaApplied: number;
  parametersBeforeTrial: SnaplockEngineParameters;
  candidateParameters: SnaplockEngineParameters;
  accepted: boolean;
  rolledBack: boolean;
  decisionReason: string;
  modelPdbHash: string;
  referencePdbHash: string | null;
  modelArtifactPath: string | null;
  groundTruthCaRmsdAngstrom: number | null;
  independentCaRmsdAngstrom: number | null;
  groundTruthBackboneRmsdAngstrom: number | null;
  groundTruthAllAtomRmsdAngstrom: number | null;
  tmScore: number | null;
  coordinateDeltaRmsdAngstrom: number;
  geometryViolationsCount: number;
  stericClashCount: number;
  internalSnaplockScoreKcal: number;
  internalL2ResidualAngstrom: number;
  pairedResidueCount: number;
  modelCoveragePct: number;
  referenceCoveragePct: number;
  sequenceIdentityPct: number;
  geometryValid: boolean;
  iterationStatus: GroundTruthIterationStatus;
  bestVerifiedCaRmsdSoFar: number | null;
}

export interface SingleTargetAutocalibrationConfig {
  maxTrials?: number;
  criteria?: Partial<GroundTruthPassCriteria>;
  initialParameters?: Partial<SnaplockEngineParameters>;
  artifactSink?: (iterationId: string, pdbText: string) => string | void;
}

export interface SingleTargetAutocalibrationResult {
  targetId: string;
  proteinName: string;
  referenceId: string;
  referenceChainId: string;
  split: CalibrationDatasetSplit;
  optimizationMethod: 'DERIVATIVE_FREE_DIRECTIONAL_PATTERN_SEARCH';
  isPipelineDifferentiable: false;
  calibrationStatus: CalibrationTerminalStatus;
  passedTargetCriterion: boolean;
  passedIntermediateCriterion: boolean;
  convergenceAchieved: boolean;
  nonConvergenceReason: string | null;
  totalTrialsExecuted: number;
  acceptedImprovementsCount: number;
  rolledBackTrialsCount: number;

  baselineCaRmsdAngstrom: number | null;
  bestVerifiedCaRmsdAngstrom: number | null;
  bestVerifiedBackboneRmsdAngstrom: number | null;
  bestVerifiedAllAtomRmsdAngstrom: number | null;
  rmsdImprovementAngstrom: number | null;

  internalBaselineScoreKcal: number;
  internalBestScoreKcal: number;
  internalBestL2ResidualAngstrom: number;

  initialParameters: SnaplockEngineParameters;
  bestVerifiedParameters: SnaplockEngineParameters;
  bestModelPdbText: string;
  bestValidationReport: GroundTruthValidationReport;
  bestRmsdComparison: AlphaFoldRmsdComparisonResult | null;
  bestFoldingMetrics: FoldingMetrics;
  criteriaApplied: GroundTruthPassCriteria;
  iterationHistory: AutocalibrationIterationRecord[];
}

export interface ProteinSplitEvaluationEntry {
  targetId: string;
  proteinName: string;
  pdbId: string;
  chainId: string;
  split: 'TRAIN' | 'VALIDATION' | 'TEST';
  residueCount: number;
  baselineCaRmsdAngstrom: number | null;
  tunedCaRmsdAngstrom: number | null;
  tunedBackboneRmsdAngstrom: number | null;
  tunedAllAtomRmsdAngstrom: number | null;
  deltaCaRmsdAngstrom: number | null;
  internalSnaplockEnergyKcal: number;
  internalL2ResidualAngstrom: number;
  pairedResidueCount: number;
  modelCoveragePct: number;
  referenceCoveragePct: number;
  sequenceIdentityPct: number;
  geometryValid: boolean;
  iterationStatus: GroundTruthIterationStatus;
  proteinCalibrationStatus: CalibrationTerminalStatus;
  failureReason: string | null;
  validationReport: GroundTruthValidationReport;
  rmsdComparison: AlphaFoldRmsdComparisonResult | null;
}

export interface SplitCheckpointRecord {
  epoch: number;
  parameterVersion: string;
  parameters: SnaplockEngineParameters;
  trainMeanCaRmsdAngstrom: number;
  validationCaRmsdAngstrom: number;
  allGeometriesValid: boolean;
  selectedAsBestValidationCheckpoint: boolean;
}

export interface AntiOverfittingSuiteResult {
  suiteStatus: 'SKALIBROWANY' | 'NIESKALIBROWANY_BRAK_ZBIEZNOSCI';
  suiteFailureSummary: string | null;
  optimizationMethod: 'DERIVATIVE_FREE_DIRECTIONAL_PATTERN_SEARCH';
  antiOverfittingProtocol: string;
  totalTuningEvaluations: number;
  selectedValidationCheckpointEpoch: number;
  initialParameters: SnaplockEngineParameters;
  selectedParametersFromValidation: SnaplockEngineParameters;
  checkpoints: SplitCheckpointRecord[];

  trainMeanBaselineCaRmsd: number | null;
  trainMeanTunedCaRmsd: number | null;
  validationBaselineCaRmsd: number | null;
  validationTunedCaRmsd: number | null;
  testMeanBaselineCaRmsd: number | null;
  testMeanTunedCaRmsd: number | null;

  // Separate explicit per-protein reports (Ubiquitin, Lysozyme, Insulin A/B, Crambin)
  trainResults: ProteinSplitEvaluationEntry[];
  validationResults: ProteinSplitEvaluationEntry[];
  testResults: ProteinSplitEvaluationEntry[];
  allProteinResults: ProteinSplitEvaluationEntry[];
}

// -----------------------------------------------------------------------------
// 1. CANONICAL EXPERIMENTAL BENCHMARK TARGETS (TRAIN / VALIDATION / HELD-OUT TEST)
// -----------------------------------------------------------------------------

function extractBenchmarkTargetFromRcsb(
  targetId: string,
  proteinName: string,
  pdbKey: string,
  chainId: string,
  split: CalibrationDatasetSplit,
  cdAlphaTargetPct: number,
  cdBetaTargetPct: number,
  explicitDisulfidePairs?: [number, number][],
  criteriaOverrides?: Partial<GroundTruthPassCriteria>
): BenchmarkProteinTarget {
  const refEntry = EXPERIMENTAL_GROUND_TRUTH_PDBS[pdbKey] ?? resolveExperimentalReferenceById(pdbKey);
  if (!refEntry) {
    throw new Error(`Brak wbudowanej struktury eksperymentalnej RCSB dla klucza ${pdbKey}.`);
  }
  const parsed = parsePdbForGroundTruth(refEntry.pdbText, { chainId });
  const sequenceCodes = parsed.residuesWithCa.map(r => r.canonicalResName);
  const n = sequenceCodes.length;
  const rgTheo = parseFloat((2.8 * Math.pow(Math.max(1, n), 0.33)).toFixed(2));
  const dMaxTheo = parseFloat((rgTheo * 2.95).toFixed(2));

  return {
    targetId,
    proteinName,
    pdbId: refEntry.id,
    chainId,
    split,
    sequenceCodes,
    boundaryConfig: {
      dMaxTarget: dMaxTheo,
      rgTarget: rgTheo,
      enforceBoundaries: true,
      cdAlphaTargetPct,
      cdBetaTargetPct,
      ssBondTolerance: 0.45,
      ssRestraintWeight: 3.5,
      explicitDisulfidePairs
    },
    referencePdbText: refEntry.pdbText,
    referenceSource: refEntry.source,
    criteria: {
      ...DEFAULT_HIGH_ACCURACY_CRITERIA,
      testId: targetId,
      testName: `${proteinName} (${refEntry.id}:${chainId})`,
      ...(criteriaOverrides ?? {})
    }
  };
}

export function getCanonicalBenchmarkTargets(
  criteriaOverrides?: Partial<GroundTruthPassCriteria>
): BenchmarkProteinTarget[] {
  return [
    // --- ZBIÓR TRENINGOWY (do strojenia parametrów) ---
    extractBenchmarkTargetFromRcsb(
      'crambin_1crn_A',
      'Krambina roślinna (1CRN, 3× S-S)',
      '1CRN',
      'A',
      'TRAIN',
      38,
      18,
      [[2, 39], [3, 31], [15, 25]],
      criteriaOverrides
    ),
    extractBenchmarkTargetFromRcsb(
      'insulin_3i40_A',
      'Insulina ludzka – Łańcuch A (3I40:A)',
      '3I40',
      'A',
      'TRAIN',
      52,
      10,
      [[5, 10], [6, 19]],
      criteriaOverrides
    ),

    // --- ZBIÓR WALIDACYJNY (do wyboru najlepszego zestawu parametrów) ---
    extractBenchmarkTargetFromRcsb(
      'ubiquitin_1ubq_A',
      'Ubikwityna ludzka (1UBQ:A, 76 aa)',
      '1UBQ',
      'A',
      'VALIDATION',
      24,
      36,
      [],
      criteriaOverrides
    ),

    // --- NIEZALEŻNY ZBIÓR TESTOWY (zamrożona ocena końcowa generalizacji) ---
    extractBenchmarkTargetFromRcsb(
      'lysozyme_1lyz_A',
      'Lizozym białka jaja kurzego HEWL (1LYZ:A, 129 aa)',
      '1LYZ',
      'A',
      'TEST',
      40,
      12,
      [[5, 126], [29, 114], [63, 79], [75, 93]],
      criteriaOverrides
    ),
    extractBenchmarkTargetFromRcsb(
      'insulin_3i40_B',
      'Insulina ludzka – Łańcuch B (3I40:B, 30 aa)',
      '3I40',
      'B',
      'TEST',
      45,
      15,
      [],
      criteriaOverrides
    )
  ];
}

// -----------------------------------------------------------------------------
// 2. CLOSED-LOOP SINGLE-TARGET AUTOCALIBRATION (WITH CHECKPOINTING & ROLLBACK)
// -----------------------------------------------------------------------------

function cloneParameters(params: SnaplockEngineParameters, nextVersion?: string): SnaplockEngineParameters {
  return {
    ...params,
    parameterVersion: nextVersion ?? params.parameterVersion
  };
}

function deriveTerminalCalibrationStatus(report: GroundTruthValidationReport): CalibrationTerminalStatus {
  if (!report.hasGroundTruthReference || report.status === 'NO_REFERENCE') {
    return 'BRAK_REFERENCJI';
  }
  if (!report.valid || report.status === 'VALIDATION_ERROR') {
    return 'BLAD_WALIDACJI';
  }
  if (report.passedTargetCriterion) {
    return 'SKALIBROWANY';
  }
  if (report.passedIntermediateCriterion) {
    return 'OSIAGNIETO_PROG_POSREDNI';
  }
  return 'BRAK_ZBIEZNOSCI';
}

export function runSingleTargetAutocalibration(
  target: BenchmarkProteinTarget,
  config?: SingleTargetAutocalibrationConfig
): SingleTargetAutocalibrationResult {
  const maxTrials = Math.max(1, Math.min(200, config?.maxTrials ?? 28));
  const criteria: GroundTruthPassCriteria = {
    ...target.criteria,
    ...(config?.criteria ?? {})
  };

  const initialParameters: SnaplockEngineParameters = {
    ...DEFAULT_SNAPLOCK_ENGINE_PARAMETERS,
    ssRestraintWeight: target.boundaryConfig.ssRestraintWeight,
    ...(config?.initialParameters ?? {}),
    parameterVersion: config?.initialParameters?.parameterVersion ?? 'v1.0.0-init'
  };

  const externalRef: ExternalReferencePdbInput | null = target.referencePdbText
    ? {
        pdbText: target.referencePdbText,
        resolvedId: target.pdbId,
        referenceName: target.proteinName,
        chainId: target.chainId,
        allowSubsequenceAlignment: false,
        source: target.referenceSource
      }
    : null;

  const iterationHistory: AutocalibrationIterationRecord[] = [];

  // --- Trial 0: Baseline Evaluation ---
  const baselineIterId = `${target.targetId}-iter-000`;
  const baselineFold = foldWithFastaAndBoundaryConstraints(
    target.sequenceCodes,
    {
      ...target.boundaryConfig,
      ssRestraintWeight: initialParameters.ssRestraintWeight
    },
    {},
    externalRef,
    initialParameters,
    {
      iterationId: baselineIterId,
      parameterVersion: initialParameters.parameterVersion,
      criteria,
      artifactSink: config?.artifactSink
    }
  );

  const baselineGt = baselineFold.groundTruthValidation;
  const baselineBackbonePhysical =
    baselineGt.modelGeometry.caBondViolations.length === 0 &&
    baselineGt.modelGeometry.nonBondedCaClashes.length === 0 &&
    !baselineGt.modelGeometry.collapsedChainDetected;
  const baselineValidForOptimization =
    baselineGt.hasGroundTruthReference &&
    baselineGt.valid &&
    baselineGt.caRmsdAngstrom !== null &&
    baselineGt.independentVerificationPassed &&
    baselineBackbonePhysical;

  let bestParams = cloneParameters(initialParameters);
  let bestFold = baselineFold;
  let bestCaRmsd: number | null = baselineValidForOptimization ? baselineGt.caRmsdAngstrom : null;

  const isHeldOutTestSplit = target.split === 'TEST';
  const baselineClashCount =
    baselineGt.modelGeometry.nonBondedCaClashes.length +
    baselineGt.modelGeometry.heavyAtomStericClashesCount;
  const baselineGeomViolations =
    baselineGt.modelGeometry.caBondViolations.length +
    baselineGt.modelGeometry.nonBondedCaClashes.length +
    (baselineGt.modelGeometry.collapsedChainDetected ? 1 : 0);

  const baselineAccepted = isHeldOutTestSplit
    ? baselineBackbonePhysical && baselineGt.modelGeometry.geometryValid
    : baselineValidForOptimization;

  iterationHistory.push({
    trialNumber: 0,
    iterationId: baselineIterId,
    parameterVersion: initialParameters.parameterVersion,
    tunedParameterKey: 'BASELINE',
    stepDeltaApplied: 0,
    parametersBeforeTrial: cloneParameters(initialParameters),
    candidateParameters: cloneParameters(initialParameters),
    accepted: baselineAccepted,
    rolledBack: !baselineAccepted,
    decisionReason: isHeldOutTestSplit
      ? `Punkt bazowy Held-Out TEST (optymalizacja wyłącznie wg energii fizycznej bez użycia referencji): E_int = ${baselineFold.metrics.totalCohesiveEnergy.toFixed(2)} kcal/mol.`
      : baselineValidForOptimization
      ? `Pomiar bazowy zweryfikowany: Ground-Truth Cα RMSD = ${baselineGt.caRmsdAngstrom?.toFixed(4)} Å.`
      : `Błąd walidacji lub geometrii w punkcie bazowym: ${baselineGt.failureReason ?? 'brak referencji'}.`,
    modelPdbHash: baselineGt.modelPdbHash,
    referencePdbHash: baselineGt.referencePdbHash,
    modelArtifactPath: baselineGt.modelArtifactPath,
    groundTruthCaRmsdAngstrom: baselineGt.caRmsdAngstrom,
    independentCaRmsdAngstrom: baselineGt.independentCaRmsdAngstrom,
    groundTruthBackboneRmsdAngstrom: baselineGt.backboneRmsdAngstrom,
    groundTruthAllAtomRmsdAngstrom: baselineGt.allHeavyAtomRmsdAngstrom,
    tmScore: baselineGt.tmScore,
    coordinateDeltaRmsdAngstrom: 0,
    geometryViolationsCount: baselineGeomViolations,
    stericClashCount: baselineClashCount,
    internalSnaplockScoreKcal: baselineFold.metrics.totalCohesiveEnergy,
    internalL2ResidualAngstrom: baselineFold.metrics.internalL2BackboneResidualAngstrom,
    pairedResidueCount: baselineGt.pairedResidueCount,
    modelCoveragePct: baselineGt.modelCoveragePct,
    referenceCoveragePct: baselineGt.referenceCoveragePct,
    sequenceIdentityPct: baselineGt.sequenceIdentityPct,
    geometryValid: baselineGt.modelGeometry.geometryValid,
    iterationStatus: baselineGt.status,
    bestVerifiedCaRmsdSoFar: bestCaRmsd
  });

  // If no valid reference on TRAIN/VAL/SINGLE_TARGET (or broken geometry on TEST), halt immediately
  if (
    (!isHeldOutTestSplit && (!baselineGt.hasGroundTruthReference || !baselineGt.valid || bestCaRmsd === null)) ||
    (isHeldOutTestSplit && !baselineBackbonePhysical)
  ) {
    const status = deriveTerminalCalibrationStatus(baselineGt);
    return {
      targetId: target.targetId,
      proteinName: target.proteinName,
      referenceId: target.pdbId,
      referenceChainId: target.chainId,
      split: target.split,
      optimizationMethod: 'DERIVATIVE_FREE_DIRECTIONAL_PATTERN_SEARCH',
      isPipelineDifferentiable: false,
      calibrationStatus: status,
      passedTargetCriterion: false,
      passedIntermediateCriterion: false,
      convergenceAchieved: false,
      nonConvergenceReason:
        baselineGt.failureReason ??
        'Brak poprawnej referencji Ground Truth lub naruszenie geometrii w iteracji bazowej.',
      totalTrialsExecuted: 1,
      acceptedImprovementsCount: 0,
      rolledBackTrialsCount: 1,
      baselineCaRmsdAngstrom: baselineGt.caRmsdAngstrom,
      bestVerifiedCaRmsdAngstrom: null,
      bestVerifiedBackboneRmsdAngstrom: null,
      bestVerifiedAllAtomRmsdAngstrom: null,
      rmsdImprovementAngstrom: null,
      internalBaselineScoreKcal: baselineFold.metrics.totalCohesiveEnergy,
      internalBestScoreKcal: baselineFold.metrics.totalCohesiveEnergy,
      internalBestL2ResidualAngstrom: baselineFold.metrics.internalL2BackboneResidualAngstrom,
      initialParameters,
      bestVerifiedParameters: bestParams,
      bestModelPdbText: baselineFold.modelPdbText,
      bestValidationReport: baselineGt,
      bestRmsdComparison: baselineFold.rmsdComparison ?? null,
      bestFoldingMetrics: baselineFold.metrics,
      criteriaApplied: criteria,
      iterationHistory
    };
  }

  // Stop immediately if baseline already meets target criterion
  let acceptedImprovementsCount = 0;
  let rolledBackTrialsCount = 0;
  let trialCounter = 0;

  const stepSizes = new Map<TunableParameterSpec['key'], number>();
  for (const spec of SNAPLOCK_TUNABLE_PARAMETER_SPECS) {
    stepSizes.set(spec.key, spec.initialStep);
  }

  let paramCursor = 0;
  let directionSign: 1 | -1 = 1;
  let consecutiveRejectsInPass = 0;

  let bestInternalPhysicalScore = baselineFold.metrics.totalCohesiveEnergy;

  while (
    trialCounter < maxTrials &&
    (isHeldOutTestSplit || !bestFold.groundTruthValidation.passedTargetCriterion)
  ) {
    const spec = SNAPLOCK_TUNABLE_PARAMETER_SPECS[paramCursor % SNAPLOCK_TUNABLE_PARAMETER_SPECS.length];
    const currentStep = stepSizes.get(spec.key) ?? spec.initialStep;

    // Check if all step sizes are exhausted below minimum resolution
    const anyStepActive = Array.from(stepSizes.values()).some((s, idx) => {
      const sp = SNAPLOCK_TUNABLE_PARAMETER_SPECS[idx];
      return s >= sp.minStep - 1e-9;
    });
    if (!anyStepActive) {
      break;
    }

    if (currentStep < spec.minStep - 1e-9) {
      paramCursor++;
      directionSign = 1;
      continue;
    }

    // 1. Save exact parameter snapshot BEFORE the trial
    const savedBeforeTrial = cloneParameters(bestParams);
    const caBeforeTrial = bestFold.positions;

    const baseVal = Number(savedBeforeTrial[spec.key]);
    const rawNextVal = baseVal + directionSign * currentStep;
    const clampedNextVal = parseFloat(
      Math.max(spec.minBound, Math.min(spec.maxBound, rawNextVal)).toFixed(4)
    );
    const actualDelta = parseFloat((clampedNextVal - baseVal).toFixed(4));

    if (Math.abs(actualDelta) < 1e-7) {
      if (directionSign === 1) {
        directionSign = -1;
      } else {
        stepSizes.set(spec.key, currentStep * 0.5);
        directionSign = 1;
        paramCursor++;
      }
      continue;
    }

    trialCounter++;
    const iterId = `${target.targetId}-iter-${String(trialCounter).padStart(3, '0')}`;
    const candidateVersion = `v1.1.${trialCounter}-${spec.key}`;
    const candidateParams: SnaplockEngineParameters = {
      ...savedBeforeTrial,
      [spec.key]: clampedNextVal,
      parameterVersion: candidateVersion
    };

    // 2. Generate brand-new 3D structure from modified engine parameters & run independent Ground-Truth validator
    const trialFold = foldWithFastaAndBoundaryConstraints(
      target.sequenceCodes,
      {
        ...target.boundaryConfig,
        ssRestraintWeight: candidateParams.ssRestraintWeight
      },
      {},
      externalRef,
      candidateParams,
      {
        iterationId: iterId,
        parameterVersion: candidateVersion,
        criteria,
        artifactSink: config?.artifactSink
      }
    );

    const gt = trialFold.groundTruthValidation;

    // Compute exact Kabsch-superposed conformational coordinate change ΔX [Å]
    const coordDeltaRmsd =
      caBeforeTrial.length >= 3 && trialFold.positions.length === caBeforeTrial.length
        ? parseFloat(kabschSuperposeSvd(trialFold.positions, caBeforeTrial).caRmsdAngstrom.toFixed(5))
        : 0;

    const geomViolationsCount =
      gt.modelGeometry.caBondViolations.length +
      gt.modelGeometry.nonBondedCaClashes.length +
      (gt.modelGeometry.collapsedChainDetected ? 1 : 0);
    const stericClashCount =
      gt.modelGeometry.nonBondedCaClashes.length +
      gt.modelGeometry.heavyAtomStericClashesCount;

    // 3. Strict validation & physical geometry gate
    const candidateBackbonePhysical =
      gt.modelGeometry.caBondViolations.length === 0 &&
      gt.modelGeometry.nonBondedCaClashes.length === 0 &&
      !gt.modelGeometry.collapsedChainDetected &&
      gt.modelGeometry.heavyAtomStericClashesCount <=
        bestFold.groundTruthValidation.modelGeometry.heavyAtomStericClashesCount;

    const nonZeroCoordinateChange = coordDeltaRmsd >= 1e-5;

    // CRITICAL: On Held-Out TEST split, NEVER use Ground-Truth RMSD for optimization acceptance!
    let improved = false;
    const candidateCaRmsd = gt.caRmsdAngstrom;
    const candidatePhysicalScore = trialFold.metrics.totalCohesiveEnergy;

    if (isHeldOutTestSplit) {
      improved =
        nonZeroCoordinateChange &&
        candidateBackbonePhysical &&
        gt.modelGeometry.geometryValid &&
        candidatePhysicalScore < bestInternalPhysicalScore - 1e-4;
    } else {
      const trialIsPhysicallyAndNumericallyValid =
        nonZeroCoordinateChange &&
        gt.hasGroundTruthReference &&
        gt.valid &&
        gt.independentVerificationPassed &&
        candidateBackbonePhysical &&
        candidateCaRmsd !== null &&
        gt.modelCoveragePct >= criteria.minModelCoveragePct &&
        (!criteria.requireExactSequenceIdentity || gt.mismatchedResidueCount === 0);

      improved =
        trialIsPhysicallyAndNumericallyValid &&
        candidateCaRmsd !== null &&
        bestCaRmsd !== null &&
        candidateCaRmsd < bestCaRmsd - 1e-4;
    }

    let decisionReason = '';
    if (improved) {
      const prevBestRmsd = bestCaRmsd;
      const prevPhys = bestInternalPhysicalScore;
      bestParams = cloneParameters(candidateParams);
      bestFold = trialFold;
      bestInternalPhysicalScore = candidatePhysicalScore;
      if (gt.caRmsdAngstrom !== null) {
        bestCaRmsd = gt.caRmsdAngstrom;
      }
      acceptedImprovementsCount++;
      consecutiveRejectsInPass = 0;
      decisionReason = isHeldOutTestSplit
        ? `ZAAKCEPTOWANO (TEST — bez referencji): spadek energii fizycznej ${prevPhys.toFixed(2)} -> ${candidatePhysicalScore.toFixed(2)} kcal/mol (ΔX=${coordDeltaRmsd.toFixed(4)} Å, kolizje=${stericClashCount}).`
        : `ZAAKCEPTOWANO: Ground-Truth Cα RMSD spadł z ${prevBestRmsd?.toFixed(4)} Å do ${candidateCaRmsd?.toFixed(4)} Å (ΔX=${coordDeltaRmsd.toFixed(4)} Å, ${spec.key} ${actualDelta >= 0 ? '+' : ''}${actualDelta}).`;
      directionSign = 1;
      paramCursor++;
    } else {
      // 4. Rollback to savedBeforeTrial immediately!
      bestParams = savedBeforeTrial;
      rolledBackTrialsCount++;
      consecutiveRejectsInPass++;

      if (!nonZeroCoordinateChange) {
        decisionReason = `COFNIĘTO (PARAMETR NIEAKTYWNY / ΔX = 0.0000 Å): zmiana ${spec.key} (${actualDelta >= 0 ? '+' : ''}${actualDelta}) nie wpłynęła na współrzędne struktury.`;
        // Disable structurally inactive parameter for remaining trials instead of wasting budget!
        if (spec.key === 'ssRestraintWeight' && (target.boundaryConfig.explicitDisulfidePairs?.length ?? 0) === 0) {
          stepSizes.set(spec.key, 0);
        }
      } else if (!isHeldOutTestSplit && (!gt.valid || !gt.independentVerificationPassed)) {
        decisionReason = `COFNIĘTO (BŁĄD WALIDACJI): ${gt.failureReason ?? 'Niezgodność weryfikacji niezależnej SVD/Kwaternion'}. Przywrócono ${savedBeforeTrial.parameterVersion}.`;
      } else if (!candidateBackbonePhysical) {
        decisionReason = `COFNIĘTO (NARUSZENIE GEOMETRII / KOLIZJE=${stericClashCount}, ΔX=${coordDeltaRmsd.toFixed(3)} Å): ${gt.modelGeometry.reasons.join(' ')} Przywrócono ${savedBeforeTrial.parameterVersion}.`;
      } else if (isHeldOutTestSplit) {
        decisionReason = `COFNIĘTO (TEST — brak poprawy energii fizycznej): E_int = ${candidatePhysicalScore.toFixed(2)} >= ${bestInternalPhysicalScore.toFixed(2)} kcal/mol (ΔX=${coordDeltaRmsd.toFixed(3)} Å).`;
      } else {
        decisionReason = `COFNIĘTO (BRAK POPRAWY / EFEKT DŹWIGNI ΔX=${coordDeltaRmsd.toFixed(3)} Å): Cα RMSD = ${candidateCaRmsd?.toFixed(4) ?? 'N/A'} Å >= najlepszy ${bestCaRmsd?.toFixed(4)} Å.`;
      }

      if (directionSign === 1 && nonZeroCoordinateChange) {
        directionSign = -1;
      } else {
        directionSign = 1;
        paramCursor++;
        if (consecutiveRejectsInPass >= SNAPLOCK_TUNABLE_PARAMETER_SPECS.length * 2) {
          for (const sp of SNAPLOCK_TUNABLE_PARAMETER_SPECS) {
            const cur = stepSizes.get(sp.key) ?? sp.initialStep;
            stepSizes.set(sp.key, cur * 0.5);
          }
          consecutiveRejectsInPass = 0;
        }
      }
    }

    iterationHistory.push({
      trialNumber: trialCounter,
      iterationId: iterId,
      parameterVersion: candidateVersion,
      tunedParameterKey: spec.key,
      stepDeltaApplied: actualDelta,
      parametersBeforeTrial: savedBeforeTrial,
      candidateParameters: candidateParams,
      accepted: improved,
      rolledBack: !improved,
      decisionReason,
      modelPdbHash: gt.modelPdbHash,
      referencePdbHash: gt.referencePdbHash,
      modelArtifactPath: gt.modelArtifactPath,
      groundTruthCaRmsdAngstrom: gt.caRmsdAngstrom,
      independentCaRmsdAngstrom: gt.independentCaRmsdAngstrom,
      groundTruthBackboneRmsdAngstrom: gt.backboneRmsdAngstrom,
      groundTruthAllAtomRmsdAngstrom: gt.allHeavyAtomRmsdAngstrom,
      tmScore: gt.tmScore,
      coordinateDeltaRmsdAngstrom: coordDeltaRmsd,
      geometryViolationsCount: geomViolationsCount,
      stericClashCount,
      internalSnaplockScoreKcal: trialFold.metrics.totalCohesiveEnergy,
      internalL2ResidualAngstrom: trialFold.metrics.internalL2BackboneResidualAngstrom,
      pairedResidueCount: gt.pairedResidueCount,
      modelCoveragePct: gt.modelCoveragePct,
      referenceCoveragePct: gt.referenceCoveragePct,
      sequenceIdentityPct: gt.sequenceIdentityPct,
      geometryValid: gt.modelGeometry.geometryValid,
      iterationStatus: gt.status,
      bestVerifiedCaRmsdSoFar: bestCaRmsd
    });
  }

  const finalGt = bestFold.groundTruthValidation;
  const terminalStatus = deriveTerminalCalibrationStatus(finalGt);
  const convergenceAchieved = finalGt.passedTargetCriterion;

  const nonConvergenceReason = convergenceAchieved
    ? null
    : `BRAK ZBIEŻNOŚCI DO CELU WYSOKIEJ DOKŁADNOŚCI (< ${criteria.targetCaRmsdAngstrom.toFixed(2)} Å): po ${trialCounter} próbach optymalizacji bezgradientowej najlepszy rzeczywiście zweryfikowany Ground-Truth Cα RMSD wyniósł ${bestCaRmsd?.toFixed(4)} Å (start: ${baselineGt.caRmsdAngstrom?.toFixed(4)} Å). Dyskretny generator de-novo (248 wektorów E8 / 2232 mikrofasetki L2 bez więzów kontaktowych dalekiego zasięgu MSA/komplementarności 3D) nie osiąga natywnego upakowania trzeciorzędowego poniżej ${criteria.targetCaRmsdAngstrom.toFixed(2)} Å.`;

  return {
    targetId: target.targetId,
    proteinName: target.proteinName,
    referenceId: target.pdbId,
    referenceChainId: target.chainId,
    split: target.split,
    optimizationMethod: 'DERIVATIVE_FREE_DIRECTIONAL_PATTERN_SEARCH',
    isPipelineDifferentiable: false,
    calibrationStatus: terminalStatus,
    passedTargetCriterion: finalGt.passedTargetCriterion,
    passedIntermediateCriterion: finalGt.passedIntermediateCriterion,
    convergenceAchieved,
    nonConvergenceReason,
    totalTrialsExecuted: trialCounter + 1,
    acceptedImprovementsCount,
    rolledBackTrialsCount,
    baselineCaRmsdAngstrom: baselineGt.caRmsdAngstrom,
    bestVerifiedCaRmsdAngstrom: bestCaRmsd,
    bestVerifiedBackboneRmsdAngstrom: finalGt.backboneRmsdAngstrom,
    bestVerifiedAllAtomRmsdAngstrom: finalGt.allHeavyAtomRmsdAngstrom,
    rmsdImprovementAngstrom:
      baselineGt.caRmsdAngstrom !== null && bestCaRmsd !== null
        ? parseFloat((baselineGt.caRmsdAngstrom - bestCaRmsd).toFixed(4))
        : null,
    internalBaselineScoreKcal: baselineFold.metrics.totalCohesiveEnergy,
    internalBestScoreKcal: bestFold.metrics.totalCohesiveEnergy,
    internalBestL2ResidualAngstrom: bestFold.metrics.internalL2BackboneResidualAngstrom,
    initialParameters,
    bestVerifiedParameters: bestParams,
    bestModelPdbText: bestFold.modelPdbText,
    bestValidationReport: finalGt,
    bestRmsdComparison: bestFold.rmsdComparison ?? null,
    bestFoldingMetrics: bestFold.metrics,
    criteriaApplied: criteria,
    iterationHistory
  };
}

// -----------------------------------------------------------------------------
// 3. ANTI-OVERFITTING TRAIN / VALIDATION / HELD-OUT TEST BENCHMARK SUITE
// -----------------------------------------------------------------------------

function evaluateTargetWithParameters(
  target: BenchmarkProteinTarget,
  params: SnaplockEngineParameters,
  iterationTag: string
) {
  const externalRef: ExternalReferencePdbInput = {
    pdbText: target.referencePdbText,
    resolvedId: target.pdbId,
    referenceName: target.proteinName,
    chainId: target.chainId,
    allowSubsequenceAlignment: false,
    source: target.referenceSource
  };

  return foldWithFastaAndBoundaryConstraints(
    target.sequenceCodes,
    {
      ...target.boundaryConfig,
      ssRestraintWeight: params.ssRestraintWeight
    },
    {},
    externalRef,
    params,
    {
      iterationId: `${target.targetId}-${iterationTag}`,
      parameterVersion: params.parameterVersion,
      criteria: target.criteria
    }
  );
}

export function runAntiOverfittingSplitSuite(options?: {
  epochs?: number;
  criteria?: Partial<GroundTruthPassCriteria>;
}): AntiOverfittingSuiteResult {
  const epochs = Math.max(1, Math.min(8, options?.epochs ?? 3));
  const targets = getCanonicalBenchmarkTargets(options?.criteria);

  const trainTargets = targets.filter(t => t.split === 'TRAIN');
  const valTargets = targets.filter(t => t.split === 'VALIDATION');
  const testTargets = targets.filter(t => t.split === 'TEST');

  const initialParams: SnaplockEngineParameters = {
    ...DEFAULT_SNAPLOCK_ENGINE_PARAMETERS,
    parameterVersion: 'v1.0.0-baseline'
  };

  // Baseline evaluation across all splits (for honest before/after comparison; TEST is NOT read by optimizer!)
  const baselineByTarget = new Map<string, ReturnType<typeof evaluateTargetWithParameters>>();
  for (const t of targets) {
    baselineByTarget.set(t.targetId, evaluateTargetWithParameters(t, initialParams, 'baseline'));
  }

  let totalEvaluations = targets.length;

  const meanCaRmsdOf = (subset: BenchmarkProteinTarget[], p: SnaplockEngineParameters, tag: string) => {
    let sum = 0;
    let validCount = 0;
    let allGeoOk = true;
    for (const t of subset) {
      const res = evaluateTargetWithParameters(t, p, tag);
      totalEvaluations++;
      const gt = res.groundTruthValidation;
      const backboneOk =
        gt.modelGeometry.caBondViolations.length === 0 &&
        gt.modelGeometry.nonBondedCaClashes.length === 0 &&
        !gt.modelGeometry.collapsedChainDetected;
      if (!gt.valid || !backboneOk || gt.caRmsdAngstrom === null) {
        allGeoOk = false;
        sum += 999.0;
      } else {
        sum += gt.caRmsdAngstrom;
      }
      validCount++;
    }
    return {
      meanRmsd: validCount > 0 ? parseFloat((sum / validCount).toFixed(4)) : 999.0,
      allGeoOk
    };
  };

  // Epoch 0 checkpoint (baseline on TRAIN + VALIDATION)
  const baseTrainMean =
    trainTargets.reduce((acc, t) => acc + (baselineByTarget.get(t.targetId)?.groundTruthValidation.caRmsdAngstrom ?? 0), 0) /
    Math.max(1, trainTargets.length);
  const baseValMean =
    valTargets.reduce((acc, t) => acc + (baselineByTarget.get(t.targetId)?.groundTruthValidation.caRmsdAngstrom ?? 0), 0) /
    Math.max(1, valTargets.length);

  const checkpoints: SplitCheckpointRecord[] = [
    {
      epoch: 0,
      parameterVersion: initialParams.parameterVersion,
      parameters: cloneParameters(initialParams),
      trainMeanCaRmsdAngstrom: parseFloat(baseTrainMean.toFixed(4)),
      validationCaRmsdAngstrom: parseFloat(baseValMean.toFixed(4)),
      allGeometriesValid: true,
      selectedAsBestValidationCheckpoint: false
    }
  ];

  let currentTrainParams = cloneParameters(initialParams);
  let currentTrainMean = baseTrainMean;

  // Tune strictly on TRAIN targets (Crambin 1CRN + Insulin Chain A 3I40:A)
  const activeSpecs = SNAPLOCK_TUNABLE_PARAMETER_SPECS.slice(0, 6);
  for (let epoch = 1; epoch <= epochs; epoch++) {
    const stepScale = Math.pow(0.6, epoch - 1);

    for (const spec of activeSpecs) {
      const step = spec.initialStep * stepScale;
      const savedBefore = cloneParameters(currentTrainParams);

      let improvedInCoord = false;
      for (const sign of [1, -1] as const) {
        const baseNum = Number(savedBefore[spec.key]);
        const nextVal = parseFloat(
          Math.max(spec.minBound, Math.min(spec.maxBound, baseNum + sign * step)).toFixed(4)
        );
        if (Math.abs(nextVal - baseNum) < 1e-7) continue;

        const cand: SnaplockEngineParameters = {
          ...savedBefore,
          [spec.key]: nextVal,
          parameterVersion: `v1.2.e${epoch}-${spec.key}`
        };

        const evalTrain = meanCaRmsdOf(trainTargets, cand, `e${epoch}-${spec.key}-${sign > 0 ? 'p' : 'm'}`);
        if (evalTrain.allGeoOk && evalTrain.meanRmsd < currentTrainMean - 1e-4) {
          currentTrainParams = cand;
          currentTrainMean = evalTrain.meanRmsd;
          improvedInCoord = true;
          break;
        }
      }

      if (!improvedInCoord) {
        currentTrainParams = savedBefore;
      }
    }

    // Evaluate epoch checkpoint on VALIDATION set (Ubiquitin 1UBQ:A) for early-stopping / model selection
    const epochVersion = `v1.2.epoch-${epoch}`;
    const checkpointParams = cloneParameters(currentTrainParams, epochVersion);
    const valEval = meanCaRmsdOf(valTargets, checkpointParams, `val-e${epoch}`);

    checkpoints.push({
      epoch,
      parameterVersion: epochVersion,
      parameters: checkpointParams,
      trainMeanCaRmsdAngstrom: parseFloat(currentTrainMean.toFixed(4)),
      validationCaRmsdAngstrom: valEval.meanRmsd,
      allGeometriesValid: valEval.allGeoOk,
      selectedAsBestValidationCheckpoint: false
    });
  }

  // Model selection strictly via VALIDATION split (prevents overfitting to TRAIN and zero leakage from TEST)
  let bestCkptIdx = 0;
  for (let i = 1; i < checkpoints.length; i++) {
    if (
      checkpoints[i].allGeometriesValid &&
      checkpoints[i].validationCaRmsdAngstrom < checkpoints[bestCkptIdx].validationCaRmsdAngstrom - 1e-4
    ) {
      bestCkptIdx = i;
    }
  }
  checkpoints[bestCkptIdx].selectedAsBestValidationCheckpoint = true;
  const selectedParameters = cloneParameters(
    checkpoints[bestCkptIdx].parameters,
    `${checkpoints[bestCkptIdx].parameterVersion}-selected`
  );

  // Final frozen evaluation across TRAIN, VALIDATION, and HELD-OUT TEST (Lysozyme 1LYZ:A + Insulin B 3I40:B)
  const buildEntry = (t: BenchmarkProteinTarget): ProteinSplitEvaluationEntry => {
    const baseEval = baselineByTarget.get(t.targetId)!;
    const tunedEval = evaluateTargetWithParameters(t, selectedParameters, 'final-frozen');
    totalEvaluations++;

    const baseRmsd = baseEval.groundTruthValidation.caRmsdAngstrom;
    const gt = tunedEval.groundTruthValidation;
    const tunedRmsd = gt.caRmsdAngstrom;

    return {
      targetId: t.targetId,
      proteinName: t.proteinName,
      pdbId: t.pdbId,
      chainId: t.chainId,
      split: t.split as 'TRAIN' | 'VALIDATION' | 'TEST',
      residueCount: t.sequenceCodes.length,
      baselineCaRmsdAngstrom: baseRmsd,
      tunedCaRmsdAngstrom: tunedRmsd,
      tunedBackboneRmsdAngstrom: gt.backboneRmsdAngstrom,
      tunedAllAtomRmsdAngstrom: gt.allHeavyAtomRmsdAngstrom,
      deltaCaRmsdAngstrom:
        baseRmsd !== null && tunedRmsd !== null
          ? parseFloat((tunedRmsd - baseRmsd).toFixed(4))
          : null,
      internalSnaplockEnergyKcal: tunedEval.metrics.totalCohesiveEnergy,
      internalL2ResidualAngstrom: tunedEval.metrics.internalL2BackboneResidualAngstrom,
      pairedResidueCount: gt.pairedResidueCount,
      modelCoveragePct: gt.modelCoveragePct,
      referenceCoveragePct: gt.referenceCoveragePct,
      sequenceIdentityPct: gt.sequenceIdentityPct,
      geometryValid: gt.modelGeometry.geometryValid,
      iterationStatus: gt.status,
      proteinCalibrationStatus: deriveTerminalCalibrationStatus(gt),
      failureReason: gt.failureReason,
      validationReport: gt,
      rmsdComparison: tunedEval.rmsdComparison ?? null
    };
  };

  const trainResults = trainTargets.map(buildEntry);
  const validationResults = valTargets.map(buildEntry);
  const testResults = testTargets.map(buildEntry);
  const allProteinResults = [...trainResults, ...validationResults, ...testResults];

  const allPassedTarget = allProteinResults.every(
    r => r.proteinCalibrationStatus === 'SKALIBROWANY' && r.validationReport.passedTargetCriterion
  );

  const avgOf = (arr: (number | null)[]) => {
    const nums = arr.filter((v): v is number => v !== null && Number.isFinite(v));
    if (nums.length === 0) return null;
    return parseFloat((nums.reduce((a, b) => a + b, 0) / nums.length).toFixed(4));
  };

  const failedProteins = allProteinResults.filter(r => r.proteinCalibrationStatus !== 'SKALIBROWANY');

  return {
    suiteStatus: allPassedTarget ? 'SKALIBROWANY' : 'NIESKALIBROWANY_BRAK_ZBIEZNOSCI',
    suiteFailureSummary: allPassedTarget
      ? null
      : `Status globalny NIESKALIBROWANY: ${failedProteins.length}/${allProteinResults.length} struktur nie spełniło kryterium Ground-Truth Cα RMSD < ${allProteinResults[0]?.validationReport.criteriaApplied.targetCaRmsdAngstrom.toFixed(2) ?? '1.00'} Å (${failedProteins.map(f => `${f.pdbId}:${f.chainId} = ${f.tunedCaRmsdAngstrom?.toFixed(2) ?? 'ERR'} Å`).join(', ')}). Status „SKALIBROWANY” wymaga rzeczywistego sukcesu wszystkich testów zbioru.`,
    optimizationMethod: 'DERIVATIVE_FREE_DIRECTIONAL_PATTERN_SEARCH',
    antiOverfittingProtocol:
      'Strojenie parametrów wyłącznie na zbiorze TRENINGOWYM (1CRN:A, 3I40:A) -> wybór punktu kontrolnego na zbiorze WALIDACYJNYM (1UBQ:A) -> pojedyncza zamrożona ewaluacja końcowa na niezależnym zbiorze TESTOWYM (1LYZ:A, 3I40:B).',
    totalTuningEvaluations: totalEvaluations,
    selectedValidationCheckpointEpoch: checkpoints[bestCkptIdx].epoch,
    initialParameters: initialParams,
    selectedParametersFromValidation: selectedParameters,
    checkpoints,
    trainMeanBaselineCaRmsd: avgOf(trainResults.map(r => r.baselineCaRmsdAngstrom)),
    trainMeanTunedCaRmsd: avgOf(trainResults.map(r => r.tunedCaRmsdAngstrom)),
    validationBaselineCaRmsd: avgOf(validationResults.map(r => r.baselineCaRmsdAngstrom)),
    validationTunedCaRmsd: avgOf(validationResults.map(r => r.tunedCaRmsdAngstrom)),
    testMeanBaselineCaRmsd: avgOf(testResults.map(r => r.baselineCaRmsdAngstrom)),
    testMeanTunedCaRmsd: avgOf(testResults.map(r => r.tunedCaRmsdAngstrom)),
    trainResults,
    validationResults,
    testResults,
    allProteinResults
  };
}
