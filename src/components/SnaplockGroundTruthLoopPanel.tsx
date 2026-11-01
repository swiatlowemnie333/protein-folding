import React, { useMemo, useState } from 'react';
import {
  Activity,
  AlertTriangle,
  BarChart3,
  CheckCircle2,
  Cpu,
  Download,
  GitCompare,
  Layers,
  Play,
  RefreshCw,
  ShieldAlert,
  ShieldCheck,
  SlidersHorizontal
} from 'lucide-react';
import {
  Bar,
  BarChart,
  CartesianGrid,
  Legend,
  Line,
  LineChart,
  ReferenceLine,
  ResponsiveContainer,
  Tooltip,
  XAxis,
  YAxis
} from 'recharts';
import {
  AntiOverfittingSuiteResult,
  BenchmarkProteinTarget,
  SingleTargetAutocalibrationResult,
  getCanonicalBenchmarkTargets,
  runAntiOverfittingSplitSuite,
  runSingleTargetAutocalibration
} from '../snaplockAutocalibrationLoop';
import {
  AblationBenchmarkSuiteReport,
  HierarchicalFoldingResult,
  ParameterSensitivityReport,
  analyzeParameterSensitivity,
  runAblationBenchmarkSuite,
  runHierarchicalFoldingEngine
} from '../snaplockPhysicalFoldingEngine';
import { triggerPdbDownload } from '../snaplockFastaAndPdb';
import { SnaplockRmsdComparisonChart } from './SnaplockRmsdComparisonChart';
import {
  GLOBAL_PROTEIN_TARGETS,
  GlobalProteinTargetId
} from '../snaplockGlobalTargets';

interface SnaplockGroundTruthLoopPanelProps {
  activeGlobalTarget?: GlobalProteinTargetId;
}

export const SnaplockGroundTruthLoopPanel: React.FC<SnaplockGroundTruthLoopPanelProps> = ({
  activeGlobalTarget = 'UBIQUITIN_76'
}) => {
  const canonicalTargets = useMemo(() => getCanonicalBenchmarkTargets(), []);

  const [selectedTargetId, setSelectedTargetId] = useState<string>(() => {
    if (activeGlobalTarget && activeGlobalTarget !== 'CUSTOM' && GLOBAL_PROTEIN_TARGETS[activeGlobalTarget]) {
      return GLOBAL_PROTEIN_TARGETS[activeGlobalTarget].autocalibrationTargetId;
    }
    return 'ubiquitin_1ubq_A';
  });
  const [simulateMissingReference, setSimulateMissingReference] = useState<boolean>(false);
  const [targetThreshold, setTargetThreshold] = useState<number>(1.0);
  const [enableIntermediateThreshold, setEnableIntermediateThreshold] = useState<boolean>(true);
  const [intermediateThreshold, setIntermediateThreshold] = useState<number>(1.5);
  const [maxTrials, setMaxTrials] = useState<number>(24);

  const [isRunningSingle, setIsRunningSingle] = useState<boolean>(false);
  const [isRunningHierarchical, setIsRunningHierarchical] = useState<boolean>(false);
  const [isRunningAblation, setIsRunningAblation] = useState<boolean>(false);
  const [isRunningSuite, setIsRunningSuite] = useState<boolean>(false);

  const [singleResult, setSingleResult] = useState<SingleTargetAutocalibrationResult | null>(() => {
    const ubq = getCanonicalBenchmarkTargets().find(t => t.targetId === 'ubiquitin_1ubq_A')!;
    return runSingleTargetAutocalibration(ubq, {
      maxTrials: 18,
      criteria: {
        targetCaRmsdAngstrom: 1.0,
        intermediateCaRmsdAngstrom: 1.5
      }
    });
  });

  const [hierarchicalResult, setHierarchicalResult] = useState<HierarchicalFoldingResult | null>(() => {
    const ubq = getCanonicalBenchmarkTargets().find(t => t.targetId === 'ubiquitin_1ubq_A')!;
    return runHierarchicalFoldingEngine(
      ubq.sequenceCodes,
      ubq.boundaryConfig,
      {
        pdbText: ubq.referencePdbText,
        resolvedId: ubq.pdbId,
        referenceName: ubq.proteinName,
        chainId: ubq.chainId,
        source: ubq.referenceSource,
        allowSubsequenceAlignment: false
      },
      {
        mode: 'M5_FULL_HIERARCHICAL_ENGINE',
        splitProtocol: ubq.split,
        candidateCount: 4,
        stageCTrialsPerCandidate: 18,
        stageDRelaxSteps: 6,
        quantizeToL2Grid: true,
        criteria: {
          targetCaRmsdAngstrom: 1.0,
          intermediateCaRmsdAngstrom: 1.5
        },
        modelIdPrefix: `hier-${ubq.targetId}`
      }
    );
  });

  const [ablationReport, setAblationReport] = useState<AblationBenchmarkSuiteReport | null>(() =>
    runAblationBenchmarkSuite({
      stageCTrialsBudget: 16,
      candidateCountM5: 3
    })
  );

  const [suiteResult, setSuiteResult] = useState<AntiOverfittingSuiteResult | null>(() =>
    runAntiOverfittingSplitSuite({
      epochs: 2,
      criteria: {
        targetCaRmsdAngstrom: 1.0,
        intermediateCaRmsdAngstrom: 1.5
      }
    })
  );

  // Sync with Single Global Protein Selector whenever activeGlobalTarget changes!
  React.useEffect(() => {
    if (activeGlobalTarget && activeGlobalTarget !== 'CUSTOM' && GLOBAL_PROTEIN_TARGETS[activeGlobalTarget]) {
      const mappedId = GLOBAL_PROTEIN_TARGETS[activeGlobalTarget].autocalibrationTargetId;
      const found = canonicalTargets.find(t => t.targetId === mappedId);
      if (found) {
        setSelectedTargetId(found.targetId);
        setSingleResult(
          runSingleTargetAutocalibration(found, {
            maxTrials: 18,
            criteria: { targetCaRmsdAngstrom: 1.0, intermediateCaRmsdAngstrom: 1.5 }
          })
        );
        setHierarchicalResult(
          runHierarchicalFoldingEngine(
            found.sequenceCodes,
            found.boundaryConfig,
            {
              pdbText: found.referencePdbText,
              resolvedId: found.pdbId,
              referenceName: found.proteinName,
              chainId: found.chainId,
              source: found.referenceSource,
              allowSubsequenceAlignment: false
            },
            {
              mode: 'M5_FULL_HIERARCHICAL_ENGINE',
              splitProtocol: found.split,
              candidateCount: 4,
              stageCTrialsPerCandidate: 18,
              stageDRelaxSteps: 6,
              quantizeToL2Grid: true,
              criteria: { targetCaRmsdAngstrom: 1.0, intermediateCaRmsdAngstrom: 1.5 },
              modelIdPrefix: `hier-${found.targetId}`
            }
          )
        );
      }
    }
  }, [activeGlobalTarget, canonicalTargets]);

  const activeTarget: BenchmarkProteinTarget = useMemo(() => {
    return canonicalTargets.find(t => t.targetId === selectedTargetId) || canonicalTargets[0];
  }, [canonicalTargets, selectedTargetId]);

  const activeSensitivity: ParameterSensitivityReport = useMemo(() => {
    return analyzeParameterSensitivity(activeTarget);
  }, [activeTarget]);

  const handleRunSingleAutocalibration = () => {
    setIsRunningSingle(true);
    setTimeout(() => {
      try {
        const targetToRun: BenchmarkProteinTarget = simulateMissingReference
          ? {
              ...activeTarget,
              pdbId: 'BRAK-REF',
              referencePdbText: ''
            }
          : activeTarget;

        const res = runSingleTargetAutocalibration(targetToRun, {
          maxTrials,
          criteria: {
            targetCaRmsdAngstrom: targetThreshold,
            intermediateCaRmsdAngstrom: enableIntermediateThreshold ? intermediateThreshold : undefined,
            minModelCoveragePct: 95.0,
            minReferenceCoveragePct: 95.0,
            requireExactSequenceIdentity: true,
            requireValidPhysicalGeometry: true
          }
        });
        setSingleResult(res);
      } finally {
        setIsRunningSingle(false);
      }
    }, 20);
  };

  const handleRunHierarchicalEngine = () => {
    setIsRunningHierarchical(true);
    setTimeout(() => {
      try {
        const extRef = simulateMissingReference
          ? null
          : {
              pdbText: activeTarget.referencePdbText,
              resolvedId: activeTarget.pdbId,
              referenceName: activeTarget.proteinName,
              chainId: activeTarget.chainId,
              source: activeTarget.referenceSource,
              allowSubsequenceAlignment: false
            };

        const res = runHierarchicalFoldingEngine(
          activeTarget.sequenceCodes,
          activeTarget.boundaryConfig,
          extRef,
          {
            mode: 'M5_FULL_HIERARCHICAL_ENGINE',
            splitProtocol: activeTarget.split,
            candidateCount: 4,
            stageCTrialsPerCandidate: Math.max(8, Math.round(maxTrials * 0.85)),
            stageDRelaxSteps: 8,
            quantizeToL2Grid: true,
            criteria: {
              targetCaRmsdAngstrom: targetThreshold,
              intermediateCaRmsdAngstrom: enableIntermediateThreshold ? intermediateThreshold : undefined
            },
            modelIdPrefix: `hier-${activeTarget.targetId}`
          }
        );
        setHierarchicalResult(res);
      } finally {
        setIsRunningHierarchical(false);
      }
    }, 20);
  };

  const handleRunAblationBenchmark = () => {
    setIsRunningAblation(true);
    setTimeout(() => {
      try {
        const rep = runAblationBenchmarkSuite({
          stageCTrialsBudget: Math.max(12, Math.min(32, maxTrials)),
          candidateCountM5: 4
        });
        setAblationReport(rep);
      } finally {
        setIsRunningAblation(false);
      }
    }, 20);
  };

  const handleRunFullSplitSuite = () => {
    setIsRunningSuite(true);
    setTimeout(() => {
      try {
        const res = runAntiOverfittingSplitSuite({
          epochs: 3,
          criteria: {
            targetCaRmsdAngstrom: targetThreshold,
            intermediateCaRmsdAngstrom: enableIntermediateThreshold ? intermediateThreshold : undefined
          }
        });
        setSuiteResult(res);
      } finally {
        setIsRunningSuite(false);
      }
    }, 20);
  };

  const handleDownloadHierarchicalPdb = () => {
    if (!hierarchicalResult) return;
    const sel = hierarchicalResult.selectedCandidate;
    triggerPdbDownload(sel.pdbText, `snaplock_hierarchical_${sel.modelId}.pdb`);
  };

  const handleDownloadVerifiedPdb = () => {
    if (!singleResult) return;
    triggerPdbDownload(
      singleResult.bestModelPdbText,
      `snaplock_autocalib_${singleResult.targetId}_${singleResult.bestVerifiedParameters.parameterVersion}.pdb`
    );
  };

  const handleDownloadAuditJson = () => {
    const payload = {
      exportedAt: new Date().toISOString(),
      hierarchicalEngineM5: hierarchicalResult,
      parameterSensitivity: activeSensitivity,
      ablationStudy: ablationReport,
      singleTargetAutocalibrationM1: singleResult,
      antiOverfittingSplitSuite: suiteResult
    };
    const blob = new Blob([JSON.stringify(payload, null, 2)], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `snaplock_next_phase_audit_${selectedTargetId}.json`;
    a.click();
    URL.revokeObjectURL(url);
  };

  const chartHistoryData = useMemo(() => {
    if (!singleResult) return [];
    return singleResult.iterationHistory.map(item => ({
      trial: `#${item.trialNumber}`,
      trialNum: item.trialNumber,
      trialCaRmsd: item.groundTruthCaRmsdAngstrom,
      bestCaRmsd: item.bestVerifiedCaRmsdSoFar,
      coordDelta: item.coordinateDeltaRmsdAngstrom,
      accepted: item.accepted ? 'TAK' : 'COFNIĘTO'
    }));
  }, [singleResult]);

  const ablationChartData = useMemo(() => {
    if (!ablationReport) return [];
    const targetIds = ['crambin_1crn_A', 'insulin_3i40_A', 'ubiquitin_1ubq_A', 'lysozyme_1lyz_A', 'insulin_3i40_B'];
    const labels: Record<string, string> = {
      crambin_1crn_A: '1CRN:A (Train)',
      insulin_3i40_A: '3I40:A (Train)',
      ubiquitin_1ubq_A: '1UBQ:A (Val)',
      lysozyme_1lyz_A: '1LYZ:A (Test)',
      insulin_3i40_B: '3I40:B (Test)'
    };

    return targetIds.map(tid => {
      const findVal = (m: string) =>
        ablationReport.deNovoRows.find(r => r.targetId === tid && r.methodId === m)?.caRmsdAngstrom ?? null;
      return {
        protein: labels[tid] || tid,
        M1_Baseline: findVal('M1_BASELINE_E8_L2'),
        M2_LocalTorsions: findVal('M2_INDEPENDENT_LOCAL_TORSIONS'),
        M3_SecondaryStruct: findVal('M3_SECONDARY_STRUCTURE_CONSTRAINTS'),
        M4_LongRange: findVal('M4_LONG_RANGE_CONTACTS'),
        M5_Hierarchical: findVal('M5_FULL_HIERARCHICAL_ENGINE')
      };
    });
  }, [ablationReport]);

  const renderStatusBadge = (status: string) => {
    if (status === 'SKALIBROWANY' || status === 'PASSED_TARGET') {
      return (
        <span className="px-2.5 py-0.5 rounded bg-emerald-950/90 border border-emerald-500 text-emerald-300 font-extrabold text-[11px] inline-flex items-center gap-1">
          <CheckCircle2 className="w-3.5 h-3.5 text-emerald-400" />
          <span>SKALIBROWANY</span>
        </span>
      );
    }
    if (status === 'OSIAGNIETO_PROG_POSREDNI' || status === 'PASSED_INTERMEDIATE') {
      return (
        <span className="px-2.5 py-0.5 rounded bg-amber-950/90 border border-amber-500 text-amber-200 font-bold text-[11px] inline-flex items-center gap-1">
          <AlertTriangle className="w-3.5 h-3.5 text-amber-400" />
          <span>PRÓG POŚREDNI (&lt;{intermediateThreshold} Å)</span>
        </span>
      );
    }
    if (status === 'BRAK_REFERENCJI' || status === 'NO_REFERENCE') {
      return (
        <span className="px-2.5 py-0.5 rounded bg-zinc-900 border border-zinc-700 text-zinc-300 font-bold text-[11px] inline-flex items-center gap-1">
          <ShieldAlert className="w-3.5 h-3.5 text-zinc-400" />
          <span>BRAK REFERENCJI GROUND-TRUTH</span>
        </span>
      );
    }
    if (status === 'BLAD_WALIDACJI' || status === 'VALIDATION_ERROR' || status === 'GEOMETRY_VIOLATION') {
      return (
        <span className="px-2.5 py-0.5 rounded bg-rose-950 border border-rose-500 text-rose-200 font-extrabold text-[11px] inline-flex items-center gap-1">
          <AlertTriangle className="w-3.5 h-3.5 text-rose-400" />
          <span>BŁĄD WALIDACJI / GEOMETRII</span>
        </span>
      );
    }
    return (
      <span className="px-2.5 py-0.5 rounded bg-rose-950/85 border border-rose-600/80 text-rose-200 font-extrabold text-[11px] inline-flex items-center gap-1">
        <AlertTriangle className="w-3.5 h-3.5 text-rose-400" />
        <span>NIESKALIBROWANY • BRAK ZBIEŻNOŚCI</span>
      </span>
    );
  };

  const gt = singleResult?.bestValidationReport ?? null;
  const hierSelected = hierarchicalResult?.selectedCandidate ?? null;
  const hierGt = hierSelected?.validationReport ?? null;

  return (
    <div className="flex-1 w-full h-full min-h-0 overflow-y-auto p-4 sm:p-5 space-y-5 font-mono text-xs bg-[#020205]">
      {/* 1. CONTROL & CONFIGURABLE PASS CRITERIA HEADER */}
      <div className="bg-[#05060e] border border-zinc-800 rounded-xl p-4 space-y-3.5">
        <div className="flex flex-wrap items-center justify-between gap-3 border-b border-zinc-900 pb-3">
          <div className="flex items-center gap-3">
            <div className="w-10 h-10 rounded-lg bg-purple-950/90 border border-purple-500/60 flex items-center justify-center text-purple-300">
              <GitCompare className="w-5 h-5" />
            </div>
            <div>
              <div className="flex items-center gap-2 flex-wrap">
                <h2 className="text-sm sm:text-base font-bold text-white">
                  SNAPLOCK Physically Constrained Engine (NeRF φ/ψ + Kontakty) &amp; Walidacja Ground-Truth
                </h2>
                <span className="text-[10px] px-2 py-0.5 rounded bg-cyan-950 text-cyan-300 border border-cyan-700/60 font-bold">
                  Niezależny Kabsch SVD + TM-score + Ablation M1–M5
                </span>
              </div>
              <p className="text-[11px] text-zinc-400 mt-0.5">
                Ścisłe rozdzielenie generowania de-novo od oceny <code className="text-purple-300">ground_truth_validation</code>. Brak użycia referencji podczas zwijania na zbiorze Held-Out Test.
              </p>
            </div>
          </div>

          <div className="flex items-center gap-2 flex-wrap">
            <button
              onClick={handleDownloadHierarchicalPdb}
              disabled={!hierarchicalResult}
              className="px-3 py-1.5 rounded bg-emerald-700 hover:bg-emerald-600 disabled:opacity-40 text-white font-bold flex items-center gap-1.5 cursor-pointer transition"
            >
              <Download className="w-3.5 h-3.5" />
              <span>Eksportuj M5 .pdb</span>
            </button>
            <button
              onClick={handleDownloadVerifiedPdb}
              disabled={!singleResult}
              className="px-3 py-1.5 rounded bg-sky-700 hover:bg-sky-600 disabled:opacity-40 text-white font-bold flex items-center gap-1.5 cursor-pointer transition"
            >
              <Download className="w-3.5 h-3.5" />
              <span>Eksportuj M1 .pdb</span>
            </button>
            <button
              onClick={handleDownloadAuditJson}
              className="px-3 py-1.5 rounded bg-zinc-900 hover:bg-zinc-800 border border-zinc-700 text-cyan-300 font-bold flex items-center gap-1.5 cursor-pointer transition"
            >
              <Download className="w-3.5 h-3.5" />
              <span>Pełny Audyt JSON</span>
            </button>
          </div>
        </div>

        {/* Controls Grid */}
        <div className="grid grid-cols-1 md:grid-cols-5 gap-3 items-end">
          <div className="md:col-span-2 space-y-1">
            <label className="text-[10px] uppercase font-bold text-zinc-400 block">
              Struktura Testowa &amp; Eksperymentalna Referencja RCSB PDB:
            </label>
            <select
              value={selectedTargetId}
              onChange={(e) => setSelectedTargetId(e.target.value)}
              className="w-full bg-black border border-zinc-800 rounded px-2.5 py-1.5 text-white font-bold focus:outline-none focus:border-purple-500"
            >
              {canonicalTargets.map(t => (
                <option key={t.targetId} value={t.targetId}>
                  [{t.split}] {t.proteinName} — Łańcuch {t.chainId} ({t.sequenceCodes.length} reszt)
                </option>
              ))}
            </select>
          </div>

          <div className="space-y-1">
            <label className="text-[10px] uppercase font-bold text-zinc-400 block">
              Cel Cα RMSD [Å] (Wys. dokładność):
            </label>
            <input
              type="number"
              step="0.1"
              min="0.2"
              max="5.0"
              value={targetThreshold}
              onChange={(e) => setTargetThreshold(parseFloat(e.target.value) || 1.0)}
              className="w-full bg-black border border-zinc-800 rounded px-2.5 py-1.5 text-emerald-300 font-bold"
            />
          </div>

          <div className="space-y-1">
            <label className="text-[10px] uppercase font-bold text-zinc-400 flex items-center justify-between">
              <span>Próg pośredni [Å]:</span>
              <input
                type="checkbox"
                checked={enableIntermediateThreshold}
                onChange={(e) => setEnableIntermediateThreshold(e.target.checked)}
                className="accent-amber-500"
              />
            </label>
            <input
              type="number"
              step="0.1"
              min="0.5"
              max="10.0"
              disabled={!enableIntermediateThreshold}
              value={intermediateThreshold}
              onChange={(e) => setIntermediateThreshold(parseFloat(e.target.value) || 1.5)}
              className="w-full bg-black border border-zinc-800 disabled:opacity-40 rounded px-2.5 py-1.5 text-amber-300 font-bold"
            />
          </div>

          <div className="space-y-1">
            <label className="text-[10px] uppercase font-bold text-zinc-400 block">
              Budżet prób optymalizatora:
            </label>
            <input
              type="number"
              step="4"
              min="8"
              max="80"
              value={maxTrials}
              onChange={(e) => setMaxTrials(Math.max(8, Math.min(80, parseInt(e.target.value, 10) || 24)))}
              className="w-full bg-black border border-zinc-800 rounded px-2.5 py-1.5 text-cyan-300 font-bold"
            />
          </div>
        </div>

        <div className="flex flex-wrap items-center justify-between gap-3 pt-1">
          <label className="flex items-center gap-2 text-[11px] text-zinc-300 cursor-pointer">
            <input
              type="checkbox"
              checked={simulateMissingReference}
              onChange={(e) => setSimulateMissingReference(e.target.checked)}
              className="accent-purple-500"
            />
            <span>
              Testuj tryb <strong>bez dostępnej referencji</strong> (czysta kontrola geometrii fizycznej — bez fałszywego 0.0 Å)
            </span>
          </label>

          <div className="flex items-center gap-2 flex-wrap">
            <button
              onClick={handleRunHierarchicalEngine}
              disabled={isRunningHierarchical}
              className="px-3.5 py-1.5 rounded bg-emerald-700 hover:bg-emerald-600 disabled:opacity-50 text-white font-bold flex items-center gap-1.5 cursor-pointer"
            >
              {isRunningHierarchical ? (
                <>
                  <RefreshCw className="w-3.5 h-3.5 animate-spin" />
                  <span>Zwijanie 5-etapowe M5...</span>
                </>
              ) : (
                <>
                  <Cpu className="w-3.5 h-3.5" />
                  <span>Zwiń Hierarchicznie M5 ({activeTarget.pdbId}:{activeTarget.chainId})</span>
                </>
              )}
            </button>

            <button
              onClick={handleRunSingleAutocalibration}
              disabled={isRunningSingle}
              className="px-3.5 py-1.5 rounded bg-purple-700 hover:bg-purple-600 disabled:opacity-50 text-white font-bold flex items-center gap-1.5 cursor-pointer"
            >
              {isRunningSingle ? (
                <>
                  <RefreshCw className="w-3.5 h-3.5 animate-spin" />
                  <span>Strojenie M1...</span>
                </>
              ) : (
                <>
                  <Play className="w-3.5 h-3.5" />
                  <span>Autokalibracja M1 ({activeTarget.pdbId}:{activeTarget.chainId})</span>
                </>
              )}
            </button>

            <button
              onClick={handleRunAblationBenchmark}
              disabled={isRunningAblation}
              className="px-3.5 py-1.5 rounded bg-amber-700 hover:bg-amber-600 disabled:opacity-50 text-white font-bold flex items-center gap-1.5 cursor-pointer"
            >
              {isRunningAblation ? (
                <>
                  <RefreshCw className="w-3.5 h-3.5 animate-spin" />
                  <span>Liczenie Ablation M1–M5...</span>
                </>
              ) : (
                <>
                  <BarChart3 className="w-3.5 h-3.5" />
                  <span>Uruchom Testy Ablation (M1–M5)</span>
                </>
              )}
            </button>

            <button
              onClick={handleRunFullSplitSuite}
              disabled={isRunningSuite}
              className="px-3.5 py-1.5 rounded bg-cyan-800 hover:bg-cyan-700 disabled:opacity-50 text-white font-bold flex items-center gap-1.5 cursor-pointer"
            >
              {isRunningSuite ? (
                <>
                  <RefreshCw className="w-3.5 h-3.5 animate-spin" />
                  <span>Train/Val/Test...</span>
                </>
              ) : (
                <>
                  <Layers className="w-3.5 h-3.5" />
                  <span>Zbiór Train / Val / Test</span>
                </>
              )}
            </button>
          </div>
        </div>
      </div>

      {/* 2. HIERARCHICAL PHYSICALLY CONSTRAINED ENGINE (M5: STAGES A -> B -> C -> D -> E) */}
      {hierarchicalResult && hierSelected && hierGt && (
        <div className="bg-[#05060e] border border-emerald-900/70 rounded-xl p-4 space-y-4">
          <div className="flex flex-wrap items-center justify-between gap-2 border-b border-zinc-900 pb-3">
            <div className="space-y-1">
              <div className="flex items-center gap-2 flex-wrap">
                <Cpu className="w-4 h-4 text-emerald-400" />
                <h3 className="text-sm font-bold text-white">
                  5-Etapowy Hierarchiczny Silnik Fizyczny M5 (Niezależne Kąty φ/ψ NeRF + Wiązania H + Kontakty Dalekiego Zasięgu)
                </h3>
                {renderStatusBadge(hierGt.status)}
                <span className="text-[10px] px-2 py-0.5 rounded bg-emerald-950/80 border border-emerald-700 text-emerald-300 font-bold">
                  Źródło więzów: {hierarchicalResult.restraintProvenance}
                </span>
              </div>
              <p className="text-[11px] text-zinc-400">
                Wybrany kandydat de-novo (wg najniższej energii fizycznej i poprawnej geometrii):{' '}
                <strong className="text-emerald-300">{hierSelected.modelId}</strong> ({hierSelected.topologyLabel}) • Hash PDB:{' '}
                <code className="text-cyan-300">{hierSelected.pdbHash}</code>
              </p>
            </div>

            <div className="text-right text-[11px] text-zinc-400 tabular-nums">
              <div>
                Czas obliczeń: <strong className="text-white">{hierarchicalResult.runtimeMs} ms</strong> • Kandydatów:{' '}
                <strong className="text-emerald-300">{hierarchicalResult.candidates.length}</strong>
              </div>
              <div>
                Prób łącznie: <strong className="text-white">{hierarchicalResult.totalTrialsExecuted}</strong> (Przyjęto:{' '}
                <strong className="text-emerald-400">{hierarchicalResult.acceptedTrialsCount}</strong>, Odrzucono:{' '}
                <strong className="text-amber-400">{hierarchicalResult.rejectedTrialsCount}</strong>, ΔX=0:{' '}
                <strong>{hierarchicalResult.zeroDisplacementTrialsCount}</strong>)
              </div>
            </div>
          </div>

          {/* Provenance & De-Novo Limitations Notice */}
          <div className="p-2.5 rounded bg-zinc-900/80 border border-zinc-800 text-[11px] text-zinc-300">
            <strong className="text-cyan-300">Rozdział generowania od oceny:</strong>{' '}
            {hierarchicalResult.deNovoLimitationsNotice}
          </div>

          {/* 6 Separate Candidate Evaluation Metrics (Section 6 Requirement) */}
          <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-6 gap-2.5 tabular-nums">
            <div className="p-3 rounded-lg bg-black/80 border border-purple-700/80">
              <div className="text-[9.5px] uppercase font-bold text-purple-300">
                1. Ground-Truth Cα RMSD
              </div>
              <div className="text-lg font-extrabold text-white mt-1">
                {hierGt.hasGroundTruthReference && hierGt.caRmsdAngstrom !== null
                  ? `${hierGt.caRmsdAngstrom.toFixed(4)} Å`
                  : 'Niedostępny'}
              </div>
              <div className="text-[9.5px] text-zinc-400 mt-0.5">
                Backbone / All-Atom: {hierGt.backboneRmsdAngstrom?.toFixed(2) ?? '—'} / {hierGt.allHeavyAtomRmsdAngstrom?.toFixed(2) ?? '—'} Å
              </div>
            </div>

            <div className="p-3 rounded-lg bg-black/80 border border-indigo-800/80">
              <div className="text-[9.5px] uppercase font-bold text-indigo-300">
                2. Globalne Podobieństwo (TM-score)
              </div>
              <div className="text-lg font-extrabold text-indigo-200 mt-1">
                {hierGt.tmScore !== null ? hierGt.tmScore.toFixed(4) : 'Niedostępny'}
              </div>
              <div className="text-[9.5px] text-zinc-400 mt-0.5">
                GDT-TS: {hierGt.gdtTsPct !== null ? `${hierGt.gdtTsPct.toFixed(1)}%` : '—'} • Pary: {hierGt.pairedResidueCount}
              </div>
            </div>

            <div className="p-3 rounded-lg bg-black/80 border border-teal-800/80">
              <div className="text-[9.5px] uppercase font-bold text-teal-300">
                3. Lokalna Geometria Łańcucha
              </div>
              <div className={`text-lg font-extrabold mt-1 ${hierGt.modelGeometry.geometryValid ? 'text-emerald-400' : 'text-rose-400'}`}>
                {hierGt.modelGeometry.geometryValid ? 'POPRAWNA' : 'NARUSZONA'}
              </div>
              <div className="text-[9.5px] text-zinc-400 mt-0.5">
                Śr. Cα-Cα: {hierGt.modelGeometry.meanCaCaDistAngstrom} Å • E_rama: {hierSelected.physicalEnergy.ramachandranEnergyKcal}
              </div>
            </div>

            <div className="p-3 rounded-lg bg-black/80 border border-rose-900/70">
              <div className="text-[9.5px] uppercase font-bold text-rose-300">
                4. Kolizje Steryczne (VdW)
              </div>
              <div className="text-lg font-extrabold text-white mt-1">
                {hierSelected.physicalEnergy.stericClashCount} kolizji
              </div>
              <div className="text-[9.5px] text-zinc-400 mt-0.5">
                Ostre kolizje Cα: {hierSelected.physicalEnergy.caClashCount} • E_vdw: {hierSelected.physicalEnergy.stericClashEnergyKcal} kcal
              </div>
            </div>

            <div className="p-3 rounded-lg bg-black/80 border border-cyan-800/80">
              <div className="text-[9.5px] uppercase font-bold text-cyan-300">
                5. Jakość Kontaktów &amp; Wiązań H
              </div>
              <div className="text-lg font-extrabold text-cyan-200 mt-1">
                {hierSelected.physicalEnergy.contactSatisfactionPct}%
              </div>
              <div className="text-[9.5px] text-zinc-400 mt-0.5">
                Kontakty: {hierSelected.physicalEnergy.satisfiedLongRangeContacts}/{hierSelected.physicalEnergy.totalLongRangeContacts} • Mostki H: {hierSelected.physicalEnergy.satisfiedHBondCount}/{hierSelected.physicalEnergy.totalTargetHBonds}
              </div>
            </div>

            <div className="p-3 rounded-lg bg-black/80 border border-amber-800/70">
              <div className="text-[9.5px] uppercase font-bold text-amber-300">
                6. Wewnętrzna Energia Fizyczna
              </div>
              <div className="text-lg font-extrabold text-amber-300 mt-1">
                {hierSelected.physicalEnergy.totalPhysicalEnergyKcal.toFixed(2)} kcal
              </div>
              <div className="text-[9.5px] text-zinc-500 mt-0.5">
                Wewn. score modelu (Rg={hierSelected.physicalEnergy.rgActualAngstrom} Å) — NIE Ground-Truth!
              </div>
            </div>
          </div>

          {/* Multi-Candidate Ensemble Comparison Table */}
          <div className="border border-zinc-800 rounded-lg overflow-hidden bg-black/60">
            <div className="px-3 py-2 bg-zinc-900/80 border-b border-zinc-800 flex items-center justify-between text-[10.5px]">
              <span className="font-bold text-zinc-200 uppercase">
                Porównanie Kandydatów Topologicznych Etapu E (Wybór wg Geometrii + Wewnętrznej Energii Fizycznej — bez podglądania RMSD!)
              </span>
              <span className="text-zinc-400">
                Segmentów SSE: {hierarchicalResult.segments.length} • Aktywnych więzów fizycznych: {hierarchicalResult.activeConstraints.length}
              </span>
            </div>
            <div className="overflow-x-auto">
              <table className="w-full text-left text-[10.5px] tabular-nums">
                <thead className="bg-zinc-950 text-zinc-400 uppercase text-[9.5px]">
                  <tr>
                    <th className="p-2">Kandydat / ID</th>
                    <th className="p-2">Archetyp Topologii (Etap B)</th>
                    <th className="p-2">Wewn. E_fiz [kcal/mol]</th>
                    <th className="p-2">Kontakty %</th>
                    <th className="p-2">Mostki H</th>
                    <th className="p-2">Kolizje</th>
                    <th className="p-2">Geometria</th>
                    <th className="p-2">GT Cα RMSD</th>
                    <th className="p-2">TM-score</th>
                    <th className="p-2">Wybór De-Novo</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-zinc-900">
                  {hierarchicalResult.candidates.map(cand => (
                    <tr
                      key={cand.modelId}
                      className={cand.selectedAsBestDeNovo ? 'bg-emerald-950/20' : 'hover:bg-zinc-900/40'}
                    >
                      <td className="p-2 font-bold text-cyan-300">{cand.modelId}</td>
                      <td className="p-2 text-zinc-300">{cand.topologyLabel}</td>
                      <td className="p-2 font-bold text-amber-300">
                        {cand.physicalEnergy.totalPhysicalEnergyKcal.toFixed(2)}
                      </td>
                      <td className="p-2 text-cyan-200">{cand.physicalEnergy.contactSatisfactionPct}%</td>
                      <td className="p-2 text-zinc-300">
                        {cand.physicalEnergy.satisfiedHBondCount}/{cand.physicalEnergy.totalTargetHBonds}
                      </td>
                      <td className="p-2 text-zinc-300">{cand.physicalEnergy.stericClashCount}</td>
                      <td className="p-2">
                        <span className={cand.geometryValid ? 'text-emerald-400 font-bold' : 'text-rose-400 font-bold'}>
                          {cand.geometryValid ? 'OK' : 'BŁĄD'}
                        </span>
                      </td>
                      <td className="p-2 font-bold text-purple-300">
                        {cand.validationReport.caRmsdAngstrom !== null
                          ? `${cand.validationReport.caRmsdAngstrom.toFixed(4)} Å`
                          : 'Niedostępny'}
                      </td>
                      <td className="p-2 text-indigo-300">
                        {cand.validationReport.tmScore !== null ? cand.validationReport.tmScore.toFixed(4) : '—'}
                      </td>
                      <td className="p-2">
                        {cand.selectedAsBestDeNovo ? (
                          <span className="px-2 py-0.5 rounded bg-emerald-950 border border-emerald-600 text-emerald-300 font-bold text-[9.5px]">
                            WYBRANY DE-NOVO
                          </span>
                        ) : (
                          <span className="text-zinc-500">Kandydat alternatywny</span>
                        )}
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>

          {/* Hierarchical Trial Audit Log (Showing ΔX Coordinate Displacement!) */}
          <div className="border border-zinc-800 rounded-lg overflow-hidden bg-black/60">
            <div className="px-3 py-2 bg-zinc-900/80 border-b border-zinc-800 flex items-center justify-between text-[10.5px]">
              <span className="font-bold text-zinc-200 uppercase">
                Rejestr Prób Etapów C i D (Parametry Przed/Po, Rzeczywista Zmiana Współrzędnych ΔX, Kolizje i Energia Fizyczna)
              </span>
              <span className="text-zinc-400">
                Pokazano pierwsze {Math.min(40, hierarchicalResult.trialAuditLog.length)} z {hierarchicalResult.trialAuditLog.length} zapisanych prób
              </span>
            </div>
            <div className="max-h-56 overflow-y-auto">
              <table className="w-full text-left text-[10.5px] tabular-nums">
                <thead className="bg-zinc-950 text-zinc-500 uppercase text-[9.5px] sticky top-0">
                  <tr>
                    <th className="p-2">Próba</th>
                    <th className="p-2">Etap / Kandydat</th>
                    <th className="p-2">Przed zmianą</th>
                    <th className="p-2">Po zmianie</th>
                    <th className="p-2">Zmiana ΔX [Å]</th>
                    <th className="p-2">Kolizje</th>
                    <th className="p-2">E_fiz (Przed→Po)</th>
                    <th className="p-2">Decyzja i przyczyna</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-zinc-900">
                  {hierarchicalResult.trialAuditLog.slice(0, 40).map(tr => (
                    <tr
                      key={`${tr.candidateId}-${tr.trialNumber}`}
                      className={tr.accepted ? 'bg-emerald-950/15' : 'hover:bg-zinc-900/40'}
                    >
                      <td className="p-2 font-bold text-zinc-300">#{tr.trialNumber}</td>
                      <td className="p-2 text-cyan-300">
                        {tr.stage === 'STAGE_C_GLOBAL_TOPOLOGY' ? 'Etap C' : 'Etap D'} ({tr.candidateId})
                      </td>
                      <td className="p-2 text-zinc-400">{tr.parameterBeforeSummary}</td>
                      <td className="p-2 text-zinc-200">{tr.parameterAfterSummary}</td>
                      <td className="p-2 font-bold text-sky-300">{tr.coordinateDeltaRmsdAngstrom.toFixed(4)} Å</td>
                      <td className="p-2 text-zinc-300">{tr.stericClashCount}</td>
                      <td className="p-2 text-amber-300">
                        {tr.physicalEnergyBeforeKcal.toFixed(1)} → {tr.physicalEnergyAfterKcal.toFixed(1)}
                      </td>
                      <td className="p-2 text-zinc-300">{tr.decisionReason}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>
        </div>
      )}

      {/* 3. PARAMETER SENSITIVITY, INACTIVE & REDUNDANT PARAMETER DIAGNOSTICS */}
      <div className="bg-[#05060e] border border-zinc-800 rounded-xl p-4 space-y-3.5">
        <div className="flex flex-wrap items-center justify-between gap-2 border-b border-zinc-900 pb-2.5">
          <div className="flex items-center gap-2">
            <SlidersHorizontal className="w-4 h-4 text-amber-400" />
            <h3 className="text-sm font-bold text-white">
              Diagnostyka Aktywności, Redundancji i Efektu Dźwigni Parametrów ({activeSensitivity.proteinName})
            </h3>
          </div>
          <div className="flex items-center gap-2 text-[10.5px] tabular-nums">
            <span className="px-2 py-0.5 rounded bg-rose-950/80 border border-rose-700 text-rose-200 font-bold">
              Nieaktywne (ΔX=0): {activeSensitivity.inactiveParametersCount}
            </span>
            <span className="px-2 py-0.5 rounded bg-amber-950/80 border border-amber-700 text-amber-200 font-bold">
              Redundantne: {activeSensitivity.redundantParametersCount}
            </span>
            <span className="px-2 py-0.5 rounded bg-purple-950/80 border border-purple-700 text-purple-200 font-bold">
              Nadmiarowa dźwignia globalna: {activeSensitivity.overcoupledGlobalParametersCount}
            </span>
          </div>
        </div>

        <div className="border border-zinc-800 rounded-lg overflow-hidden bg-black/60">
          <table className="w-full text-left text-[10.5px] tabular-nums">
            <thead className="bg-zinc-950 text-zinc-400 uppercase text-[9.5px]">
              <tr>
                <th className="p-2">Rodzina Silnika</th>
                <th className="p-2">Parametr</th>
                <th className="p-2">Perturbacja</th>
                <th className="p-2">Zmiana Konformacji ΔX [Å]</th>
                <th className="p-2">Max Przesunięcie Cα [Å]</th>
                <th className="p-2">Klasyfikacja</th>
                <th className="p-2">Wyjaśnienie Fizyczne / Geometryczne</th>
              </tr>
            </thead>
            <tbody className="divide-y divide-zinc-900">
              {activeSensitivity.entries.map(entry => (
                <tr key={entry.parameterKey} className="hover:bg-zinc-900/40">
                  <td className="p-2 text-zinc-400">
                    {entry.engineFamily === 'HIERARCHICAL_NERF_TORSIONS' ? 'NeRF M5 (Lokalny)' : 'Bazowy M1 (Globalny)'}
                  </td>
                  <td className="p-2 font-bold text-white">{entry.label}</td>
                  <td className="p-2 text-cyan-300">{entry.perturbationApplied}</td>
                  <td className="p-2 font-bold text-sky-300">{entry.coordinateDeltaRmsdAngstrom.toFixed(5)} Å</td>
                  <td className="p-2 text-zinc-300">{entry.maxResidueShiftAngstrom.toFixed(3)} Å</td>
                  <td className="p-2">
                    <span
                      className={`px-2 py-0.5 rounded text-[9.5px] font-bold border ${
                        entry.classification === 'ACTIVE'
                          ? 'bg-emerald-950/80 border-emerald-700 text-emerald-300'
                          : entry.classification === 'INACTIVE_ZERO_DISPLACEMENT'
                          ? 'bg-rose-950/80 border-rose-700 text-rose-200'
                          : entry.classification === 'REDUNDANT_COLLINEAR'
                          ? 'bg-amber-950/80 border-amber-700 text-amber-200'
                          : 'bg-purple-950/80 border-purple-700 text-purple-200'
                      }`}
                    >
                      {entry.classification}
                    </span>
                  </td>
                  <td className="p-2 text-zinc-300">{entry.diagnosticExplanation}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>

      {/* 4. REPRODUCIBLE 5-WAY ABLATION STUDY (SECTION 7) */}
      {ablationReport && (
        <div className="bg-[#05060e] border border-amber-900/60 rounded-xl p-4 space-y-4">
          <div className="flex flex-wrap items-center justify-between gap-2 border-b border-zinc-900 pb-3">
            <div>
              <div className="flex items-center gap-2">
                <BarChart3 className="w-4 h-4 text-amber-400" />
                <h3 className="text-sm font-bold text-white">
                  Testy Ablation (M1–M5): Wpływ Lokalnych Kątów Torsyjnych, Struktur Drugorzędowych i Kontaktów Dalekiego Zasięgu
                </h3>
              </div>
              <p className="text-[11px] text-zinc-400 mt-1">{ablationReport.scientificConclusionPl}</p>
            </div>
          </div>

          {/* Ablation BarChart across 5 proteins */}
          <div className="bg-black/70 border border-zinc-800 rounded-lg p-3 space-y-2">
            <div className="text-[10.5px] font-bold text-zinc-200 uppercase">
              Porównanie Ground-Truth Cα RMSD [Å] Metod De-Novo M1–M5 na Identycznych Strukturach (niższy słupek = lepsze dopasowanie)
            </div>
            <div className="w-full h-60">
              <ResponsiveContainer width="100%" height="100%">
                <BarChart data={ablationChartData} margin={{ top: 8, right: 20, left: 4, bottom: 4 }}>
                  <CartesianGrid strokeDasharray="3 3" stroke="#18181b" />
                  <XAxis dataKey="protein" tick={{ fill: '#d4d4d8', fontSize: 10 }} />
                  <YAxis tick={{ fill: '#a1a1aa', fontSize: 10 }} unit=" Å" />
                  <Tooltip
                    contentStyle={{ backgroundColor: '#09090b', borderColor: '#3f3f46', fontSize: '11px' }}
                  />
                  <Legend wrapperStyle={{ fontSize: '10.5px' }} />
                  <Bar dataKey="M1_Baseline" name="M1: Bazowy E8/L2" fill="#71717a" />
                  <Bar dataKey="M2_LocalTorsions" name="M2: Lokalne (φ, ψ) NeRF" fill="#38bdf8" />
                  <Bar dataKey="M3_SecondaryStruct" name="M3: + Struktury α/β i H-bond" fill="#a855f7" />
                  <Bar dataKey="M4_LongRange" name="M4: + Kontakty Dalekie" fill="#f59e0b" />
                  <Bar dataKey="M5_Hierarchical" name="M5: Pełny Hierarchiczny (A→E)" fill="#10b981" />
                </BarChart>
              </ResponsiveContainer>
            </div>
          </div>

          {/* Aggregate Method Summary Table */}
          <div className="border border-zinc-800 rounded-lg overflow-hidden bg-black/60">
            <div className="px-3 py-2 bg-zinc-900/80 border-b border-zinc-800 text-[10.5px] font-bold text-zinc-200 uppercase">
              Podsumowanie Zbiorcze Metod Ablation i Kosztu Obliczeniowego
            </div>
            <div className="overflow-x-auto">
              <table className="w-full text-left text-[10.5px] tabular-nums">
                <thead className="bg-zinc-950 text-zinc-400 uppercase text-[9.5px]">
                  <tr>
                    <th className="p-2">Wariant Silnika</th>
                    <th className="p-2">Tryb</th>
                    <th className="p-2">Śr. Train Cα RMSD</th>
                    <th className="p-2">Val (1UBQ:A) Cα RMSD</th>
                    <th className="p-2">Held-Out Test Cα RMSD</th>
                    <th className="p-2">Śr. TM-score</th>
                    <th className="p-2">Kontakty %</th>
                    <th className="p-2">Kolizje (Suma)</th>
                    <th className="p-2">Próby / Czas [ms]</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-zinc-900">
                  {ablationReport.methodSummaries.map(m => (
                    <tr
                      key={m.methodId}
                      className={
                        !m.isDeNovo
                          ? 'bg-amber-950/15'
                          : m.methodId === 'M5_FULL_HIERARCHICAL_ENGINE'
                          ? 'bg-emerald-950/15'
                          : 'hover:bg-zinc-900/40'
                      }
                    >
                      <td className="p-2 font-bold text-white">{m.methodLabel}</td>
                      <td className="p-2">
                        {m.isDeNovo ? (
                          <span className="px-1.5 py-0.5 rounded bg-emerald-950 border border-emerald-700 text-emerald-300 text-[9.5px] font-bold">
                            DE-NOVO
                          </span>
                        ) : (
                          <span className="px-1.5 py-0.5 rounded bg-amber-950 border border-amber-600 text-amber-200 text-[9.5px] font-bold">
                            DIAGNOSTYKA KONTROLNA
                          </span>
                        )}
                      </td>
                      <td className="p-2 text-cyan-300">{m.meanTrainCaRmsd?.toFixed(3) ?? '—'} Å</td>
                      <td className="p-2 font-bold text-purple-300">{m.validationUbqCaRmsd?.toFixed(3) ?? '—'} Å</td>
                      <td className="p-2 font-bold text-emerald-300">
                        {m.meanHeldOutTestCaRmsd !== null ? `${m.meanHeldOutTestCaRmsd.toFixed(3)} Å` : 'Zablokowano (TEST)'}
                      </td>
                      <td className="p-2 text-indigo-300">{m.meanTmScoreAll?.toFixed(4) ?? '—'}</td>
                      <td className="p-2 text-zinc-300">{m.meanContactSatisfactionPct}%</td>
                      <td className="p-2 text-zinc-300">{m.totalStericClashesAcrossSuite}</td>
                      <td className="p-2 text-zinc-400">
                        {m.totalTrialsExecuted} prób / {m.totalRuntimeMs.toFixed(1)} ms
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>
        </div>
      )}

      {/* 5. SINGLE-TARGET M1 AUTOCALIBRATION HISTORY (WITH ΔX & HELD-OUT TEST PROTECTION) */}
      {singleResult && gt && (
        <div className="bg-[#05060e] border border-zinc-800 rounded-xl p-4 space-y-4">
          <div className="flex flex-wrap items-center justify-between gap-2 border-b border-zinc-900 pb-3">
            <div className="space-y-0.5">
              <div className="flex items-center gap-2 flex-wrap">
                <Activity className="w-4 h-4 text-purple-400" />
                <span className="text-sm font-bold text-white">
                  Pętla Autokalibracji M1 ({singleResult.proteinName})
                </span>
                {renderStatusBadge(singleResult.calibrationStatus)}
                <span className="text-[10.5px] px-2 py-0.5 rounded bg-zinc-900 border border-zinc-800 text-zinc-300">
                  Iteracja wybrana: <strong>{gt.iterationId}</strong> ({gt.parameterVersion})
                </span>
              </div>
              <div className="text-[11px] text-zinc-400">
                Referencja: <strong className="text-purple-300">{gt.referenceName ?? 'Brak referencji'}</strong>{' '}
                {gt.referenceId ? `(${gt.referenceId}:${gt.selectedReferenceChain ?? '-'})` : ''} • Hash modelu:{' '}
                <code className="text-cyan-300">{gt.modelPdbHash}</code>
              </div>
            </div>

            <div className="text-right text-[11px] text-zinc-400 tabular-nums">
              <div>Prób wykonanych: <strong className="text-white">{singleResult.totalTrialsExecuted}</strong></div>
              <div>
                Zaakceptowano: <strong className="text-emerald-400">{singleResult.acceptedImprovementsCount}</strong> • Cofnięto (Rollback):{' '}
                <strong className="text-amber-400">{singleResult.rolledBackTrialsCount}</strong>
              </div>
            </div>
          </div>

          {chartHistoryData.length > 1 && gt.hasGroundTruthReference && (
            <div className="bg-black/70 border border-zinc-800 rounded-lg p-3 space-y-2">
              <div className="flex items-center justify-between">
                <span className="font-bold text-zinc-200 uppercase text-[10.5px]">
                  Historia Iteracji M1 (Ground-Truth Cα RMSD vs Próg Celu)
                </span>
                <span className="text-[10px] text-zinc-400">
                  Na zbiorze Held-Out TEST akceptacja używa wyłącznie wewnętrznej energii fizycznej
                </span>
              </div>
              <div className="w-full h-48">
                <ResponsiveContainer width="100%" height="100%">
                  <LineChart data={chartHistoryData} margin={{ top: 8, right: 20, left: 4, bottom: 4 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#18181b" />
                    <XAxis dataKey="trial" tick={{ fill: '#a1a1aa', fontSize: 10 }} />
                    <YAxis tick={{ fill: '#a1a1aa', fontSize: 10 }} unit=" Å" />
                    <Tooltip
                      contentStyle={{ backgroundColor: '#09090b', borderColor: '#3f3f46', fontSize: '11px' }}
                    />
                    <Legend wrapperStyle={{ fontSize: '11px' }} />
                    <ReferenceLine
                      y={targetThreshold}
                      stroke="#10b981"
                      strokeDasharray="4 4"
                      label={{ value: `Cel ${targetThreshold} Å`, fill: '#10b981', fontSize: 10 }}
                    />
                    <Line
                      type="monotone"
                      dataKey="trialCaRmsd"
                      name="Kandydat Cα RMSD [Å]"
                      stroke="#a855f7"
                      strokeWidth={1.6}
                      dot={{ r: 2.5 }}
                    />
                    <Line
                      type="stepAfter"
                      dataKey="bestCaRmsd"
                      name="Wybrany Cα RMSD [Å]"
                      stroke="#22d3ee"
                      strokeWidth={2.4}
                      dot={false}
                    />
                  </LineChart>
                </ResponsiveContainer>
              </div>
            </div>
          )}

          {singleResult.bestRmsdComparison && gt.hasGroundTruthReference && (
            <SnaplockRmsdComparisonChart
              rmsdComparison={singleResult.bestRmsdComparison}
              compact={true}
            />
          )}

          <div className="border border-zinc-800 rounded-lg overflow-hidden bg-black/60">
            <div className="px-3 py-2 bg-zinc-900/80 border-b border-zinc-800 flex items-center justify-between text-[10.5px]">
              <span className="font-bold text-zinc-200 uppercase">
                Rejestr Iteracji M1 z Rzeczywistą Zmianą Współrzędnych ΔX ({singleResult.iterationHistory.length} rekordów)
              </span>
              <span className="text-zinc-400">
                Każda próba rejestruje ΔX [Å], kolizje, TM-score i podpis FNV-64
              </span>
            </div>
            <div className="max-h-60 overflow-y-auto">
              <table className="w-full text-left text-[10.5px] tabular-nums">
                <thead className="bg-zinc-950 text-zinc-500 uppercase text-[9.5px] sticky top-0">
                  <tr>
                    <th className="p-2">Iteracja</th>
                    <th className="p-2">Parametr (Δ)</th>
                    <th className="p-2">Zmiana ΔX [Å]</th>
                    <th className="p-2">Kolizje</th>
                    <th className="p-2">GT Cα RMSD</th>
                    <th className="p-2">TM-score</th>
                    <th className="p-2">Wewn. ΔG</th>
                    <th className="p-2">Decyzja Optymalizatora</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-zinc-900">
                  {singleResult.iterationHistory.map(row => (
                    <tr
                      key={row.iterationId}
                      className={row.accepted ? 'bg-emerald-950/15' : 'hover:bg-zinc-900/40'}
                    >
                      <td className="p-2 font-bold text-zinc-300">#{row.trialNumber}</td>
                      <td className="p-2 text-zinc-300">
                        {row.tunedParameterKey}{' '}
                        {row.trialNumber > 0 ? `(${row.stepDeltaApplied >= 0 ? '+' : ''}${row.stepDeltaApplied})` : ''}
                      </td>
                      <td className="p-2 font-bold text-sky-300">{row.coordinateDeltaRmsdAngstrom.toFixed(4)} Å</td>
                      <td className="p-2 text-zinc-300">{row.stericClashCount}</td>
                      <td className="p-2 font-bold text-purple-300">
                        {row.groundTruthCaRmsdAngstrom !== null
                          ? `${row.groundTruthCaRmsdAngstrom.toFixed(4)} Å`
                          : 'Niedostępny'}
                      </td>
                      <td className="p-2 text-indigo-300">
                        {row.tmScore !== null ? row.tmScore.toFixed(4) : '—'}
                      </td>
                      <td className="p-2 text-amber-300">{row.internalSnaplockScoreKcal} kcal</td>
                      <td className="p-2">
                        <span
                          className={`px-1.5 py-0.5 rounded text-[9.5px] font-bold mr-1.5 ${
                            row.accepted
                              ? 'bg-emerald-950 text-emerald-300 border border-emerald-700'
                              : 'bg-zinc-900 text-amber-300 border border-zinc-700'
                          }`}
                        >
                          {row.accepted ? 'ZAAKCEPTOWANO' : 'ROLLBACK'}
                        </span>
                        <span className="text-zinc-400">{row.decisionReason}</span>
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>
        </div>
      )}

      {/* 6. ANTI-OVERFITTING TRAIN / VALIDATION / HELD-OUT TEST BENCHMARK SUITE */}
      {suiteResult && (
        <div className="bg-[#05060e] border border-zinc-800 rounded-xl p-4 space-y-4">
          <div className="flex flex-wrap items-center justify-between gap-3 border-b border-zinc-900 pb-3">
            <div>
              <div className="flex items-center gap-2.5 flex-wrap">
                <ShieldCheck className="w-4 h-4 text-cyan-400" />
                <h3 className="text-sm font-bold text-white">
                  Protokół Anty-Przeuczeniowy: Rozdział Zbiorów Treningowego, Walidacyjnego i Niezależnego Testowego
                </h3>
                {suiteResult.suiteStatus === 'SKALIBROWANY' ? (
                  <span className="px-2.5 py-0.5 rounded bg-emerald-950 border border-emerald-500 text-emerald-300 font-extrabold text-[11px]">
                    ZBIÓR SKALIBROWANY
                  </span>
                ) : (
                  <span className="px-2.5 py-0.5 rounded bg-rose-950 border border-rose-600 text-rose-200 font-extrabold text-[11px]">
                    STATUS ZBIORU: NIESKALIBROWANY (BRAK ZBIEŻNOŚCI &lt; {targetThreshold.toFixed(1)} Å)
                  </span>
                )}
              </div>
              <p className="text-[11px] text-zinc-400 mt-1">{suiteResult.antiOverfittingProtocol}</p>
            </div>

            <div className="text-right text-[11px] text-zinc-400 tabular-nums">
              <div>Wybrany punkt kontrolny z Walidacji: <strong className="text-cyan-300">Epoka #{suiteResult.selectedValidationCheckpointEpoch}</strong></div>
              <div>Łącznie ewaluacji PDB: <strong className="text-white">{suiteResult.totalTuningEvaluations}</strong></div>
            </div>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-3 gap-3 tabular-nums">
            <div className="p-3 rounded-lg bg-black/75 border border-cyan-900/70">
              <div className="text-[10px] uppercase font-bold text-cyan-400">
                1. Zbiór Treningowy (Strojenie: 1CRN:A, Insulina 3I40:A)
              </div>
              <div className="text-base font-extrabold text-white mt-1">
                Śr. Cα RMSD: {suiteResult.trainMeanBaselineCaRmsd?.toFixed(3)} Å ➔{' '}
                <span className="text-cyan-300">{suiteResult.trainMeanTunedCaRmsd?.toFixed(3)} Å</span>
              </div>
            </div>

            <div className="p-3 rounded-lg bg-black/75 border border-purple-900/70">
              <div className="text-[10px] uppercase font-bold text-purple-400">
                2. Zbiór Walidacyjny (Wybór Parametrów: Ubikwityna 1UBQ:A)
              </div>
              <div className="text-base font-extrabold text-white mt-1">
                Cα RMSD: {suiteResult.validationBaselineCaRmsd?.toFixed(3)} Å ➔{' '}
                <span className="text-purple-300">{suiteResult.validationTunedCaRmsd?.toFixed(3)} Å</span>
              </div>
            </div>

            <div className="p-3 rounded-lg bg-black/75 border border-emerald-900/70">
              <div className="text-[10px] uppercase font-bold text-emerald-400">
                3. Niezależny Zbiór Testowy (Lizozym 1LYZ:A, Insulina 3I40:B)
              </div>
              <div className="text-base font-extrabold text-white mt-1">
                Śr. Cα RMSD: {suiteResult.testMeanBaselineCaRmsd?.toFixed(3)} Å ➔{' '}
                <span className="text-emerald-300">{suiteResult.testMeanTunedCaRmsd?.toFixed(3)} Å</span>
              </div>
            </div>
          </div>

          <div className="border border-zinc-800 rounded-lg overflow-hidden bg-black/60">
            <table className="w-full text-left text-[11px] tabular-nums">
              <thead className="bg-zinc-950 text-zinc-400 uppercase text-[9.5px]">
                <tr>
                  <th className="p-2.5">Zbiór</th>
                  <th className="p-2.5">Białko &amp; Referencja</th>
                  <th className="p-2.5">Łańcuch / Pokrycie</th>
                  <th className="p-2.5">Przed (Bazowy Cα)</th>
                  <th className="p-2.5">Po Strojeniu (GT Cα RMSD)</th>
                  <th className="p-2.5">TM-score</th>
                  <th className="p-2.5">Wewn. Score SNAPLOCK</th>
                  <th className="p-2.5">Geometria</th>
                  <th className="p-2.5">Status Kalibracji Białka</th>
                </tr>
              </thead>
              <tbody className="divide-y divide-zinc-900">
                {suiteResult.allProteinResults.map(entry => (
                  <tr key={entry.targetId} className="hover:bg-zinc-900/40">
                    <td className="p-2.5">
                      <span
                        className={`px-2 py-0.5 rounded text-[9.5px] font-bold border ${
                          entry.split === 'TRAIN'
                            ? 'bg-cyan-950/70 border-cyan-700 text-cyan-200'
                            : entry.split === 'VALIDATION'
                            ? 'bg-purple-950/70 border-purple-700 text-purple-200'
                            : 'bg-emerald-950/70 border-emerald-700 text-emerald-200'
                        }`}
                      >
                        {entry.split}
                      </span>
                    </td>
                    <td className="p-2.5 font-bold text-white">{entry.proteinName}</td>
                    <td className="p-2.5 text-zinc-300">
                      {entry.pdbId}:{entry.chainId} ({entry.pairedResidueCount}/{entry.residueCount} aa • {entry.modelCoveragePct}%)
                    </td>
                    <td className="p-2.5 text-zinc-400">
                      {entry.baselineCaRmsdAngstrom !== null ? `${entry.baselineCaRmsdAngstrom.toFixed(4)} Å` : '-'}
                    </td>
                    <td className="p-2.5 font-extrabold text-purple-300">
                      {entry.tunedCaRmsdAngstrom !== null ? `${entry.tunedCaRmsdAngstrom.toFixed(4)} Å` : 'Niedostępny'}
                      {entry.deltaCaRmsdAngstrom !== null && (
                        <span className={`ml-1.5 text-[10px] ${entry.deltaCaRmsdAngstrom <= 0 ? 'text-emerald-400' : 'text-rose-400'}`}>
                          ({entry.deltaCaRmsdAngstrom >= 0 ? '+' : ''}{entry.deltaCaRmsdAngstrom.toFixed(3)} Å)
                        </span>
                      )}
                    </td>
                    <td className="p-2.5 text-indigo-300">
                      {entry.validationReport.tmScore?.toFixed(4) ?? '—'}
                    </td>
                    <td className="p-2.5 text-amber-300">
                      {entry.internalSnaplockEnergyKcal.toFixed(1)} kcal (L2: {entry.internalL2ResidualAngstrom.toFixed(2)} Å)
                    </td>
                    <td className="p-2.5">
                      <span className={entry.geometryValid ? 'text-emerald-400 font-bold' : 'text-rose-400 font-bold'}>
                        {entry.geometryValid ? 'OK' : 'BŁĄD'}
                      </span>
                    </td>
                    <td className="p-2.5">{renderStatusBadge(entry.proteinCalibrationStatus)}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      )}
    </div>
  );
};
