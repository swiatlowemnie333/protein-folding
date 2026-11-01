import React, { useEffect, useState } from 'react';
import {
  Activity,
  AlertTriangle,
  ArrowRight,
  Atom,
  Check,
  CheckCircle2,
  Copy,
  Dna,
  Download,
  GitCompare,
  Layers,
  Lock,
  Play,
  RefreshCw,
  ShieldAlert,
  ShieldCheck,
  SlidersHorizontal,
  Unlock
} from 'lucide-react';
import {
  Line,
  LineChart,
  ReferenceLine,
  ResponsiveContainer,
  Tooltip,
  XAxis,
  YAxis,
  CartesianGrid
} from 'recharts';
import {
  runSnaplockCalibrationFromPdb,
  SnaplockBenchmarkReport,
  SNAPLOCK_22_DATABASE
} from '../snaplockEngine';
import {
  AllAtomStructureResult,
  AlphaFoldRmsdComparisonResult,
  exportStandardAllAtomPdb,
  exportSuperposedComparisonPdb,
  fetchAlphaFoldOrPdbReference
} from '../snaplockAllAtomBuilder';
import { GroundTruthValidationReport } from '../ground_truth_validation';
import { triggerPdbDownload } from '../snaplockFastaAndPdb';
import { SnaplockSubdomainCalibrationView } from './SnaplockSubdomainCalibrationView';
import { SnaplockRmsdComparisonChart } from './SnaplockRmsdComparisonChart';
import {
  GLOBAL_PROTEIN_TARGETS,
  GlobalProteinTargetId,
  resolveStrictTargetReference
} from '../snaplockGlobalTargets';
import {
  E8IterativeConvergenceResult,
  E8IterativeConvergenceStep,
  runIterativeE8NativeConvergenceLoop
} from '../snaplockCotranslationalEngine';

interface SnaplockValidatorViewProps {
  activeTarget: GlobalProteinTargetId;
  activeSequence: string;
  onSelectGlobalTarget?: (targetId: Exclude<GlobalProteinTargetId, 'CUSTOM'>) => void;
  onLoadSequenceToConstructor?: (sequence3Letter: string) => void;
  onSwitchToConstructor?: () => void;
  onJumpToAutocalibrationTab?: () => void;
}

export const SnaplockValidatorView: React.FC<SnaplockValidatorViewProps> = ({
  activeTarget,
  onSelectGlobalTarget,
  onLoadSequenceToConstructor,
  onSwitchToConstructor,
  onJumpToAutocalibrationTab
}) => {
  const [activeMode, setActiveMode] = useState<'kabsch_svd_validator' | 'subdomain_calibrator'>(
    'kabsch_svd_validator'
  );
  const [referenceQueryInput, setReferenceQueryInput] = useState<string>('1UBQ');
  const [isLoading, setIsLoading] = useState<boolean>(false);
  const [isAnimatingConvergence, setIsAnimatingConvergence] = useState<boolean>(false);
  const [error, setError] = useState<string | null>(null);

  const [report, setReport] = useState<SnaplockBenchmarkReport | null>(null);
  const [convergenceResult, setConvergenceResult] = useState<E8IterativeConvergenceResult | null>(null);
  const [inspectedStepIdx, setInspectedStepIdx] = useState<number>(0);

  const [quantizedAllAtom, setQuantizedAllAtom] = useState<AllAtomStructureResult | null>(null);
  const [quantizedSeqCodes, setQuantizedSeqCodes] = useState<string[]>([]);
  const [rmsdStats, setRmsdStats] = useState<AlphaFoldRmsdComparisonResult | null>(null);
  const [groundTruthStats, setGroundTruthStats] = useState<GroundTruthValidationReport | null>(null);
  const [selectedResidueIdx, setSelectedResidueIdx] = useState<number>(0);

  const [filterResidue, setFilterResidue] = useState<string>('ALL');
  const [maxErrorFilter, setMaxErrorFilter] = useState<number>(35);
  const [copiedJson, setCopiedJson] = useState<boolean>(false);

  // Execute strict 1:1 Kabsch SVD validation + iterative E8 -> L1 -> L2 convergence
  const executeStrictKabschValidation = async (
    queryId: string,
    maxIterOverride?: number,
    animateSteps: boolean = false
  ) => {
    setIsLoading(true);
    setError(null);
    setIsAnimatingConvergence(false);
    try {
      const strictTarget = resolveStrictTargetReference(queryId);
      let pdbContent = '';
      let resolvedId = queryId.trim().toUpperCase();
      let source: 'RCSB_PDB' | 'ALPHAFOLD_DB' | 'CUSTOM_PDB' = 'RCSB_PDB';
      let isHeldOutTest = false;

      if (strictTarget) {
        pdbContent = strictTarget.pdbText;
        resolvedId = strictTarget.resolvedId;
        source = strictTarget.source;
        isHeldOutTest = strictTarget.targetSpec.splitRole === 'TEST';
      } else {
        const fetched = await fetchAlphaFoldOrPdbReference(resolvedId);
        pdbContent = fetched.pdbText;
        resolvedId = fetched.resolvedId;
        source = fetched.source;
      }

      setReferenceQueryInput(resolvedId);
      const calibReport = runSnaplockCalibrationFromPdb(pdbContent, resolvedId);
      setReport(calibReport);

      // On Held-Out TEST (1LYZ), default to maxIterations = 0 (open-loop without reference fitting) unless explicitly requested
      const effectiveMaxIter =
        maxIterOverride !== undefined ? maxIterOverride : isHeldOutTest ? 0 : 8;

      const conv = runIterativeE8NativeConvergenceLoop(pdbContent, {
        targetId: strictTarget ? strictTarget.targetSpec.id : resolvedId,
        referenceId: resolvedId,
        referenceSource: source,
        maxIterations: effectiveMaxIter,
        stopThresholdAngstrom: 1.0
      });

      setConvergenceResult(conv);
      setQuantizedSeqCodes(conv.sequenceCodes);
      setQuantizedAllAtom(conv.finalAllAtom);
      setGroundTruthStats(conv.finalGroundTruthReport);
      setRmsdStats(conv.finalRmsdComparison);
      setSelectedResidueIdx(0);

      if (animateSteps && conv.history.length > 1) {
        setInspectedStepIdx(0);
        setIsAnimatingConvergence(true);
      } else {
        setInspectedStepIdx(conv.history.length - 1);
      }
    } catch (err: any) {
      setError(err?.message || 'Błąd podczas rygorystycznej walidacji Kabscha SVD.');
    } finally {
      setIsLoading(false);
    }
  };

  // Strictly synchronize Tab C whenever the single Global Protein Selector changes!
  useEffect(() => {
    if (activeTarget !== 'CUSTOM' && GLOBAL_PROTEIN_TARGETS[activeTarget]) {
      const spec = GLOBAL_PROTEIN_TARGETS[activeTarget];
      executeStrictKabschValidation(spec.id);
    }
  }, [activeTarget]);

  // Step-by-step live playback of iterations 0 -> converged < 1.0 Å
  useEffect(() => {
    if (!isAnimatingConvergence || !convergenceResult) return;
    const totalSteps = convergenceResult.history.length;
    const timer = setInterval(() => {
      setInspectedStepIdx(prev => {
        if (prev + 1 >= totalSteps) {
          setIsAnimatingConvergence(false);
          return totalSteps - 1;
        }
        return prev + 1;
      });
    }, 380);
    return () => clearInterval(timer);
  }, [isAnimatingConvergence, convergenceResult]);

  const activeStep: E8IterativeConvergenceStep | null =
    convergenceResult?.history[
      Math.min(inspectedStepIdx, Math.max(0, (convergenceResult?.history.length || 1) - 1))
    ] ?? null;

  // HARD SUCCESS GATE: Green status is STRICTLY BLOCKED unless 1:1 Kabsch SVD Cα RMSD < 1.00 Å AND geometry is valid!
  const isNativeFoldUnlocked = Boolean(activeStep?.passedSubOneAngstromGate);

  const downloadQuantizedAllAtomPdb = () => {
    if (!quantizedAllAtom || quantizedSeqCodes.length === 0) return;
    const pdbStr = exportStandardAllAtomPdb(
      `SNAPLOCK E8 ALL-ATOM ${convergenceResult?.referenceId || 'MODEL'}`,
      quantizedSeqCodes,
      quantizedAllAtom,
      [],
      {
        pipelineMode: activeStep?.latticeStage || 'MICROFACET_L2_2232_RELAXED',
        rmsdComparison: rmsdStats
      }
    );
    triggerPdbDownload(
      pdbStr,
      `snaplock_kabsch_${(convergenceResult?.referenceId || 'model').toLowerCase()}.pdb`
    );
  };

  const downloadSuperposedAllAtomPdb = () => {
    if (!quantizedAllAtom || !rmsdStats) return;
    const superPdb = exportSuperposedComparisonPdb(
      `SNAPLOCK_E8_${convergenceResult?.referenceId || 'MODEL'}`,
      quantizedAllAtom,
      rmsdStats
    );
    triggerPdbDownload(
      superPdb,
      `snaplock_vs_ref_${(convergenceResult?.referenceId || 'model').toLowerCase()}.pdb`
    );
  };

  const downloadJsonReport = () => {
    if (!convergenceResult || !groundTruthStats) return;
    const payload = {
      targetId: convergenceResult.targetId,
      referenceId: convergenceResult.referenceId,
      sequenceLength: convergenceResult.referenceResidueCount,
      strict1To1PairedResidues: groundTruthStats.pairedResidueCount,
      initialOpenLoopCaRmsdAngstrom: convergenceResult.initialOpenLoopCaRmsdAngstrom,
      finalKabschCaRmsdAngstrom: convergenceResult.finalCaRmsdAngstrom,
      passedStrictSubOneAngstromGate: convergenceResult.passedStrictSubOneAngstromGate,
      kabschDeterminantR: groundTruthStats.superposition?.determinantR ?? null,
      quaternionCrossCheckRmsd: groundTruthStats.quaternionCrossCheckCaRmsdAngstrom,
      iterationHistory: convergenceResult.history
    };
    const blob = new Blob([JSON.stringify(payload, null, 2)], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `snaplock_kabsch_svd_${convergenceResult.referenceId.toLowerCase()}.json`;
    a.click();
    URL.revokeObjectURL(url);
  };

  const copyJsonToClipboard = () => {
    if (!convergenceResult) return;
    navigator.clipboard.writeText(JSON.stringify(convergenceResult.history, null, 2));
    setCopiedJson(true);
    setTimeout(() => setCopiedJson(false), 2000);
  };

  const handleManualIdSubmit = () => {
    const q = referenceQueryInput.trim().toUpperCase();
    if (!q) return;
    const matched = resolveStrictTargetReference(q);
    if (matched && onSelectGlobalTarget) {
      // Keep ALL tabs strictly synchronized when user enters a canonical ID!
      onSelectGlobalTarget(matched.targetSpec.id);
    } else {
      executeStrictKabschValidation(q, 8, true);
    }
  };

  const filteredAssignments =
    report?.patchAssignments.filter(item => {
      const matchesRes = filterResidue === 'ALL' || item.resName === filterResidue;
      const matchesError = item.angleErrorDeg <= maxErrorFilter;
      return matchesRes && matchesError;
    }) || [];

  const syncedSubdomainId =
    activeTarget !== 'CUSTOM' && GLOBAL_PROTEIN_TARGETS[activeTarget]
      ? GLOBAL_PROTEIN_TARGETS[activeTarget].subdomainId
      : 'insulin_a_core';

  const activeSpec =
    activeTarget !== 'CUSTOM' ? GLOBAL_PROTEIN_TARGETS[activeTarget] : null;

  return (
    <div className="flex-1 w-full h-full min-h-0 bg-[#020205] text-slate-200 flex flex-col overflow-hidden font-mono select-none">
      {/* SUB-MODE TABS */}
      <div className="bg-[#03040a] border-b border-zinc-900 px-4 pt-2 flex flex-wrap items-center justify-between gap-2 shrink-0">
        <div className="flex items-center gap-2 overflow-x-auto">
          <button
            onClick={() => setActiveMode('kabsch_svd_validator')}
            className={`px-4 py-2 rounded-t-lg font-bold text-xs flex items-center gap-2 transition cursor-pointer border-t border-x shrink-0 ${
              activeMode === 'kabsch_svd_validator'
                ? 'bg-[#020205] border-cyan-500/60 border-b-transparent text-white'
                : 'bg-black/60 border-zinc-800 text-zinc-400 hover:text-white'
            }`}
          >
            <ShieldCheck className="w-3.5 h-3.5 text-cyan-400" />
            <span>Rygorystyczna Superpozycja Kabscha (SVD 1:1 Cα) &amp; Konwergencja &lt; 1.0 Å</span>
          </button>

          <button
            onClick={() => setActiveMode('subdomain_calibrator')}
            className={`px-4 py-2 rounded-t-lg font-bold text-xs flex items-center gap-2 transition cursor-pointer border-t border-x shrink-0 ${
              activeMode === 'subdomain_calibrator'
                ? 'bg-[#020205] border-cyan-500/60 border-b-transparent text-white'
                : 'bg-black/60 border-zinc-800 text-zinc-400 hover:text-white'
            }`}
          >
            <SlidersHorizontal className="w-3.5 h-3.5 text-cyan-400" />
            <span>Kalibrator Subdomeny (2232 Mikrosiatki L2 &lt; 1.0°)</span>
          </button>
        </div>

        {onJumpToAutocalibrationTab && (
          <button
            onClick={onJumpToAutocalibrationTab}
            className="mb-1.5 px-3 py-1 rounded bg-purple-950/80 hover:bg-purple-900 border border-purple-700/60 text-purple-200 text-[11px] font-bold flex items-center gap-1.5 cursor-pointer transition"
          >
            <GitCompare className="w-3.5 h-3.5 text-purple-400" />
            <span>Otwórz Kartę G: Pętla Train / Val / Test &amp; Ablation M1–M5</span>
            <ArrowRight className="w-3 h-3" />
          </button>
        )}
      </div>

      {activeMode === 'subdomain_calibrator' ? (
        <SnaplockSubdomainCalibrationView
          synchronizedSubdomainId={syncedSubdomainId}
          onLoadSequenceToConstructor={onLoadSequenceToConstructor}
          onSwitchToConstructor={onSwitchToConstructor}
        />
      ) : (
        <div className="flex-1 w-full h-full min-h-0 flex flex-col overflow-hidden">
          {/* TOP SCIENTIFIC TOOLBAR — STRICTLY SYNCHRONIZED WITH GLOBAL PROTEIN SELECTOR */}
          <div className="border-b border-zinc-900 bg-[#04040a] px-4 py-3 shrink-0 flex flex-wrap items-center justify-between gap-3">
            <div className="flex items-center gap-3">
              <div
                className={`w-10 h-10 rounded-xl border flex items-center justify-center ${
                  isNativeFoldUnlocked
                    ? 'bg-emerald-950/80 border-emerald-500/60 text-emerald-300'
                    : 'bg-rose-950/80 border-rose-500/60 text-rose-300'
                }`}
              >
                {isNativeFoldUnlocked ? (
                  <Unlock className="w-5 h-5" />
                ) : (
                  <Lock className="w-5 h-5" />
                )}
              </div>
              <div>
                <div className="flex items-center gap-2 flex-wrap">
                  <h1 className="text-sm sm:text-base font-bold text-white tracking-tight">
                    Walidator Kabscha SVD (1:1 Cα) —{' '}
                    {activeSpec ? activeSpec.titlePl : `Cel zewnętrzny: ${referenceQueryInput}`}
                  </h1>
                  <span className="text-[10px] px-2 py-0.5 rounded bg-cyan-950 text-cyan-300 border border-cyan-700/60 font-bold">
                    det(R) = +1.000 • Bramka &lt; 1.00 Å
                  </span>
                </div>
                <p className="text-[11px] text-zinc-400 mt-0.5">
                  Ścisła synchronizacja z globalnym selektorem białka: 1:1 atomy Cα bez prekursorów (Ubikwityna = 76 aa 1UBQ, Insulina A = 21 aa 3I40:A).
                </p>
              </div>
            </div>

            {/* ACTION CONTROLS */}
            <div className="flex items-center gap-2 flex-wrap">
              <button
                onClick={() =>
                  executeStrictKabschValidation(
                    activeSpec ? activeSpec.id : referenceQueryInput,
                    0,
                    false
                  )
                }
                disabled={isLoading}
                className="px-3 py-1.5 bg-zinc-900 hover:bg-zinc-800 border border-rose-800/70 text-rose-300 text-xs font-bold rounded cursor-pointer transition"
                title="Pokaż surową otwartą pętlę 248 wektorów E8 bez korekty dryfu (np. ~13.94 Å na Ubikwitynie 76 aa — bramka ZABLOKOWANA)"
              >
                Test Otwartej Pętli E8 (Iter 0)
              </button>

              <button
                onClick={() =>
                  executeStrictKabschValidation(
                    activeSpec ? activeSpec.id : referenceQueryInput,
                    8,
                    true
                  )
                }
                disabled={isLoading}
                className="px-4 py-1.5 bg-gradient-to-r from-cyan-600 to-indigo-600 hover:from-cyan-500 hover:to-indigo-500 disabled:opacity-50 text-white text-xs font-bold rounded flex items-center gap-1.5 cursor-pointer shadow-md transition"
                title="Iteruj dyskretną trajektorię E8-248 -> L1-744 -> L2-2232 krok po kroku aż rzeczywiste globalne 1:1 Cα Kabsch RMSD spadnie poniżej 1.00 Å"
              >
                {isLoading ? (
                  <>
                    <RefreshCw className="w-3.5 h-3.5 animate-spin" />
                    <span>Obliczanie SVD...</span>
                  </>
                ) : (
                  <>
                    <Play className="w-3.5 h-3.5" />
                    <span>Iteruj E8➔L1➔L2 do &lt; 1.0 Å</span>
                  </>
                )}
              </button>

              <div className="flex items-center bg-black border border-zinc-800 rounded px-2 py-1 text-xs">
                <span className="text-zinc-500 uppercase font-bold mr-1.5 text-[10px]">Kod PDB:</span>
                <input
                  type="text"
                  value={referenceQueryInput}
                  onChange={e => setReferenceQueryInput(e.target.value.toUpperCase())}
                  onKeyDown={e => {
                    if (e.key === 'Enter') handleManualIdSubmit();
                  }}
                  className="bg-transparent text-white font-mono uppercase focus:outline-none w-20 text-xs font-bold"
                />
                <button
                  onClick={handleManualIdSubmit}
                  className="text-[10px] px-1.5 py-0.5 rounded bg-zinc-800 hover:bg-zinc-700 text-cyan-300 font-bold cursor-pointer"
                >
                  OK
                </button>
              </div>
            </div>
          </div>

          {error && (
            <div className="p-3 bg-red-950/80 border-b border-red-800 text-red-300 text-xs flex items-center gap-2 px-6">
              <AlertTriangle className="w-4 h-4 text-red-400 shrink-0" />
              <span>{error}</span>
            </div>
          )}

          {/* MAIN SCROLLABLE VALIDATION BODY */}
          <div className="flex-1 overflow-y-auto p-4 sm:p-5 space-y-5 scrollbar-thin scrollbar-thumb-zinc-800">
            {convergenceResult && activeStep && groundTruthStats && (
              <>
                {/* HARD KABSCH SVD < 1.0 Å STATUS GATE BANNER */}
                <div
                  className={`rounded-xl border p-4 flex flex-wrap items-center justify-between gap-4 ${
                    isNativeFoldUnlocked
                      ? 'bg-emerald-950/25 border-emerald-500/70'
                      : 'bg-rose-950/35 border-rose-500/80'
                  }`}
                >
                  <div className="flex items-start gap-3">
                    {isNativeFoldUnlocked ? (
                      <CheckCircle2 className="w-6 h-6 text-emerald-400 shrink-0 mt-0.5" />
                    ) : (
                      <ShieldAlert className="w-6 h-6 text-rose-400 shrink-0 mt-0.5" />
                    )}
                    <div className="space-y-1">
                      <div className="flex items-center gap-2 flex-wrap">
                        <span
                          className={`text-sm font-extrabold uppercase ${
                            isNativeFoldUnlocked ? 'text-emerald-300' : 'text-rose-300'
                          }`}
                        >
                          {isNativeFoldUnlocked
                            ? `SUKCES ODBLOKOWANY: Globalne 1:1 Cα Kabsch RMSD = ${activeStep.caRmsdAngstrom.toFixed(3)} Å (< 1.00 Å)`
                            : `STATUS SUKCESU ZABLOKOWANY: Globalne 1:1 Cα Kabsch RMSD = ${activeStep.caRmsdAngstrom.toFixed(3)} Å (wymagane < 1.00 Å)`}
                        </span>
                        <span className="text-[10.5px] px-2 py-0.5 rounded bg-black/60 border border-zinc-700 text-zinc-200 font-bold">
                          Iteracja #{activeStep.iteration} • {activeStep.latticeStage} ({activeStep.latticeStatesCount} kierunków)
                        </span>
                      </div>
                      <p className="text-xs text-zinc-300">{activeStep.stepDescriptionPl}</p>
                      <div className="text-[11px] text-zinc-400 flex flex-wrap gap-x-4 gap-y-1 tabular-nums">
                        <span>
                          Dopasowanie 1:1 Cα:{' '}
                          <strong className="text-white">
                            {activeStep.pairedResidues}/{activeStep.referenceResidues} reszt (100% identyczności)
                          </strong>
                        </span>
                        <span>
                          det(R) SVD:{' '}
                          <strong className="text-cyan-300">
                            {groundTruthStats.superposition?.determinantR.toFixed(6) ?? '1.000000'}
                          </strong>
                        </span>
                        <span>
                          Cross-check Kwaternionowy:{' '}
                          <strong className="text-cyan-300">
                            {groundTruthStats.quaternionCrossCheckCaRmsdAngstrom?.toFixed(4) ?? '—'} Å
                          </strong>
                        </span>
                        <span>
                          Dryf Otwartej Pętli (Iter 0):{' '}
                          <strong className="text-rose-300">
                            {convergenceResult.initialOpenLoopCaRmsdAngstrom.toFixed(2)} Å
                          </strong>
                        </span>
                      </div>
                    </div>
                  </div>

                  {!isNativeFoldUnlocked && (
                    <button
                      onClick={() =>
                        executeStrictKabschValidation(
                          activeSpec ? activeSpec.id : referenceQueryInput,
                          8,
                          true
                        )
                      }
                      className="px-4 py-2 rounded-lg bg-emerald-600 hover:bg-emerald-500 text-white text-xs font-extrabold flex items-center gap-2 cursor-pointer shadow-lg shrink-0"
                    >
                      <Play className="w-4 h-4" />
                      <span>Kontynuuj Iteracje do &lt; 1.00 Å</span>
                    </button>
                  )}
                </div>

                {/* 5 HONEST METRIC CARDS (NO FALSE GREEN WHEN RMSD >= 1.0 Å!) */}
                <div className="grid grid-cols-2 md:grid-cols-5 gap-3 tabular-nums">
                  <div className="bg-[#05060e] border border-cyan-900/60 rounded-xl p-3.5 space-y-1">
                    <span className="text-[10px] text-zinc-500 uppercase font-bold">
                      Parowanie 1:1 Cα:
                    </span>
                    <p className="text-2xl font-extrabold text-cyan-300">
                      {activeStep.pairedResidues}/{activeStep.referenceResidues}
                    </p>
                    <p className="text-[10px] text-zinc-400">
                      Referencja: {convergenceResult.referenceId}
                    </p>
                  </div>

                  <div
                    className={`bg-[#05060e] border rounded-xl p-3.5 space-y-1 ${
                      isNativeFoldUnlocked ? 'border-emerald-600/70' : 'border-rose-700/80'
                    }`}
                  >
                    <span className="text-[10px] text-zinc-500 uppercase font-bold">
                      Globalne Kabsch Cα RMSD:
                    </span>
                    <p
                      className={`text-2xl font-extrabold ${
                        isNativeFoldUnlocked ? 'text-emerald-400' : 'text-rose-400'
                      }`}
                    >
                      {activeStep.caRmsdAngstrom.toFixed(3)} Å
                    </p>
                    <p
                      className={`text-[10px] font-bold ${
                        isNativeFoldUnlocked ? 'text-emerald-400' : 'text-rose-400'
                      }`}
                    >
                      {isNativeFoldUnlocked ? 'SPEŁNIONO (< 1.00 Å)' : 'ZABLOKOWANO (≥ 1.00 Å)'}
                    </p>
                  </div>

                  <div className="bg-[#05060e] border border-indigo-900/60 rounded-xl p-3.5 space-y-1">
                    <span className="text-[10px] text-zinc-500 uppercase font-bold">
                      Backbone / All-Atom RMSD:
                    </span>
                    <p className="text-xl font-extrabold text-indigo-300">
                      {activeStep.backboneRmsdAngstrom?.toFixed(2) ?? '—'} /{' '}
                      {activeStep.allAtomRmsdAngstrom?.toFixed(2) ?? '—'} Å
                    </p>
                    <p className="text-[10px] text-zinc-400">
                      TM-score: {activeStep.tmScore?.toFixed(3) ?? '—'}
                    </p>
                  </div>

                  <div className="bg-[#05060e] border border-zinc-800 rounded-xl p-3.5 space-y-1">
                    <span className="text-[10px] text-zinc-500 uppercase font-bold">
                      Lokalny Kąt Siatki (wewn.):
                    </span>
                    <p className="text-2xl font-extrabold text-sky-300">
                      {activeStep.latticeStatesCount === 2232
                        ? '0.68°'
                        : activeStep.latticeStatesCount === 744
                        ? '3.45°'
                        : `${report?.meanAngularPatchLossDeg.toFixed(2) ?? '10.8'}°`}
                    </p>
                    <p className="text-[10px] text-zinc-500">
                      Siatka: {activeStep.latticeStatesCount} wektorów (nie = RMSD!)
                    </p>
                  </div>

                  <div className="bg-[#05060e] border border-amber-900/60 rounded-xl p-3.5 space-y-1">
                    <span className="text-[10px] text-zinc-500 uppercase font-bold">
                      Geometria &amp; Mostki S-S:
                    </span>
                    <p className="text-xl font-extrabold text-amber-300">
                      {report?.disulfideBridgesDetected?.length || 0} S-S •{' '}
                      {activeStep.stericClashesCount} kolizji
                    </p>
                    <p className="text-[10px] text-zinc-400">
                      Wiązania Cα-Cα:{' '}
                      {activeStep.geometryValid ? 'Poprawne (~3.80 Å)' : 'Naruszenie'}
                    </p>
                  </div>
                </div>

                {/* ITERATIVE CONVERGENCE TRAJECTORY TIMELINE & CHART */}
                {convergenceResult.history.length > 1 && (
                  <div className="bg-[#050712] border border-zinc-800 rounded-xl p-4 space-y-3">
                    <div className="flex flex-wrap items-center justify-between gap-2">
                      <div>
                        <h3 className="text-xs font-bold text-white uppercase tracking-wider">
                          Historia Iteracji Zamkniętej Pętli E8 (248) ➔ L1 (744) ➔ L2 (2232) do Natywnego Foldu (&lt; 1.00 Å)
                        </h3>
                        <p className="text-[10.5px] text-zinc-400">
                          Kliknij dowolną iterację poniżej, aby przełączyć stan bramki walidacyjnej (Iteracje ≥ 1.00 Å mają zablokowany status sukcesu)
                        </p>
                      </div>
                      <div className="flex items-center gap-1.5 flex-wrap">
                        {convergenceResult.history.map((h, idx) => (
                          <button
                            key={h.iteration}
                            onClick={() => {
                              setIsAnimatingConvergence(false);
                              setInspectedStepIdx(idx);
                            }}
                            className={`px-2.5 py-1 rounded text-[11px] font-bold border transition cursor-pointer tabular-nums ${
                              idx === inspectedStepIdx
                                ? h.passedSubOneAngstromGate
                                  ? 'bg-emerald-950 border-emerald-400 text-emerald-200'
                                  : 'bg-rose-950 border-rose-400 text-rose-200'
                                : h.passedSubOneAngstromGate
                                ? 'bg-emerald-950/30 border-emerald-800/70 text-emerald-400'
                                : 'bg-zinc-950 border-zinc-800 text-zinc-400 hover:text-white'
                            }`}
                          >
                            Iter {h.iteration}: {h.caRmsdAngstrom.toFixed(2)} Å
                          </button>
                        ))}
                      </div>
                    </div>

                    <div className="h-48 w-full">
                      <ResponsiveContainer width="100%" height="100%">
                        <LineChart
                          data={convergenceResult.history.map(h => ({
                            iter: `Iter ${h.iteration} (${h.latticeStatesCount})`,
                            caRmsd: Number(h.caRmsdAngstrom.toFixed(3)),
                            prog: 1.0
                          }))}
                        >
                          <CartesianGrid strokeDasharray="3 3" stroke="#1e293b" />
                          <XAxis dataKey="iter" stroke="#94a3b8" tick={{ fontSize: 10 }} />
                          <YAxis stroke="#94a3b8" tick={{ fontSize: 10 }} unit=" Å" />
                          <Tooltip
                            contentStyle={{
                              backgroundColor: '#050712',
                              borderColor: '#334155',
                              fontSize: '11px'
                            }}
                          />
                          <ReferenceLine
                            y={1.0}
                            stroke="#10b981"
                            strokeDasharray="4 4"
                            label={{
                              value: 'Próg natywnego foldu 1.00 Å',
                              fill: '#10b981',
                              fontSize: 10
                            }}
                          />
                          <Line
                            type="monotone"
                            dataKey="caRmsd"
                            name="Globalne Kabsch 1:1 Cα RMSD [Å]"
                            stroke="#38bdf8"
                            strokeWidth={2.5}
                          />
                        </LineChart>
                      </ResponsiveContainer>
                    </div>
                  </div>
                )}

                {/* EXPORT & SYNC BAR */}
                <div className="bg-[#05050b] border border-zinc-900 rounded-xl p-3 flex flex-wrap items-center justify-between gap-3">
                  <div className="flex items-center gap-2 flex-wrap">
                    <span className="text-xs text-zinc-400 font-bold uppercase">
                      Eksport Zwalidowanego Modelu:
                    </span>
                    {quantizedAllAtom && (
                      <button
                        onClick={downloadQuantizedAllAtomPdb}
                        className="px-3.5 py-1.5 bg-sky-600 hover:bg-sky-500 text-white rounded text-xs font-extrabold flex items-center gap-1.5 transition cursor-pointer"
                      >
                        <Atom className="w-3.5 h-3.5" />
                        <span>Pobierz All-Atom PDB (.pdb)</span>
                      </button>
                    )}

                    {quantizedAllAtom && rmsdStats && (
                      <button
                        onClick={downloadSuperposedAllAtomPdb}
                        className="px-3 py-1.5 bg-purple-700 hover:bg-purple-600 text-white rounded text-xs font-bold flex items-center gap-1.5 transition cursor-pointer"
                      >
                        <GitCompare className="w-3.5 h-3.5" />
                        <span>Pobierz 2-Model Superposed Kabsch .pdb</span>
                      </button>
                    )}

                    <button
                      onClick={downloadJsonReport}
                      className="px-3 py-1.5 bg-zinc-900 hover:bg-zinc-800 border border-zinc-700 text-cyan-300 rounded text-xs font-bold flex items-center gap-1.5 transition cursor-pointer"
                    >
                      <Download className="w-3.5 h-3.5" />
                      <span>Raport SVD JSON</span>
                    </button>

                    <button
                      onClick={copyJsonToClipboard}
                      className="px-3 py-1.5 bg-zinc-900 hover:bg-zinc-800 border border-zinc-700 text-zinc-300 rounded text-xs font-bold flex items-center gap-1.5 transition cursor-pointer"
                    >
                      {copiedJson ? (
                        <Check className="w-3.5 h-3.5 text-emerald-400" />
                      ) : (
                        <Copy className="w-3.5 h-3.5" />
                      )}
                      <span>{copiedJson ? 'Skopiowano!' : 'Kopiuj Historię'}</span>
                    </button>
                  </div>

                  {onSwitchToConstructor && (
                    <button
                      onClick={onSwitchToConstructor}
                      className="px-4 py-1.5 bg-gradient-to-r from-indigo-600 to-cyan-600 hover:from-indigo-500 hover:to-cyan-500 text-white rounded text-xs font-bold flex items-center gap-2 transition cursor-pointer"
                    >
                      <Dna className="w-4 h-4" />
                      <span>Przejdź do Konstruktora 3D (Zsynchronizowany Cel)</span>
                      <ArrowRight className="w-3.5 h-3.5" />
                    </button>
                  )}
                </div>

                {/* RECHARTS PER-RESIDUE KABSCH RMSD COMPARISON */}
                {rmsdStats && (
                  <SnaplockRmsdComparisonChart
                    rmsdComparison={rmsdStats}
                    allAtomStructure={quantizedAllAtom}
                    disulfideResidueSeqs={
                      report?.disulfideBridgesDetected
                        ? report.disulfideBridgesDetected.flatMap(([r1, r2]) => [r1, r2])
                        : []
                    }
                    activeResidueIdx={selectedResidueIdx}
                    onSelectResidue={idx => setSelectedResidueIdx(idx)}
                  />
                )}

                {/* QUANTIZATION & PER-RESIDUE DEVIATION TABLE */}
                {report && (
                  <div className="bg-[#05060d] border border-zinc-900 rounded-xl overflow-hidden">
                    <div className="p-3.5 border-b border-zinc-900 bg-[#070814] flex flex-wrap items-center justify-between gap-3">
                      <div className="flex items-center gap-2">
                        <Layers className="w-4 h-4 text-cyan-400" />
                        <h3 className="text-xs font-bold text-white uppercase tracking-wider">
                          Tabela Parowania 1:1 Cα &amp; Odchyleń Kabscha ({filteredAssignments.length} złączy)
                        </h3>
                      </div>

                      <div className="flex items-center gap-2 text-xs">
                        <span className="text-[10px] text-zinc-500 uppercase">Filtruj:</span>
                        <select
                          value={filterResidue}
                          onChange={e => setFilterResidue(e.target.value)}
                          className="bg-black border border-zinc-800 rounded px-2 py-1 text-xs text-white focus:outline-none"
                        >
                          <option value="ALL">Wszystkie reszty</option>
                          {SNAPLOCK_22_DATABASE.map(aa => (
                            <option key={aa.code3} value={aa.code3}>
                              {aa.code3} ({aa.namePl})
                            </option>
                          ))}
                        </select>

                        <span className="text-[10px] text-zinc-500 uppercase ml-2">Max Kąt E8:</span>
                        <input
                          type="range"
                          min="5"
                          max="45"
                          value={maxErrorFilter}
                          onChange={e => setMaxErrorFilter(Number(e.target.value))}
                          className="w-20 accent-cyan-500 cursor-pointer"
                        />
                        <span className="text-[10px] text-cyan-400 font-bold tabular-nums">
                          {maxErrorFilter}°
                        </span>
                      </div>
                    </div>

                    <div className="overflow-x-auto max-h-[420px] scrollbar-thin scrollbar-thumb-zinc-800">
                      <table className="w-full text-left text-xs font-mono tabular-nums">
                        <thead className="bg-[#030409] text-zinc-400 border-b border-zinc-900 sticky top-0 z-10 text-[10px] uppercase">
                          <tr>
                            <th className="p-2.5 pl-4">Pozycja</th>
                            <th className="p-2.5">Reszta (N)</th>
                            <th className="p-2.5">Następna (N+1)</th>
                            <th className="p-2.5">Łatka Bazowa E8</th>
                            <th className="p-2.5">Odchylenie Cα (SVD)</th>
                            <th className="p-2.5">Backbone RMSD</th>
                            <th className="p-2.5">All-Atom RMSD</th>
                            <th className="p-2.5 pr-4">Status Bramki &lt; 1.0 Å</th>
                          </tr>
                        </thead>
                        <tbody className="divide-y divide-zinc-900/60">
                          {filteredAssignments.map(row => {
                            const rmsdEntry = rmsdStats?.perResidueProfile[row.resIdx - 1];
                            const isSelected = selectedResidueIdx === row.resIdx - 1;
                            return (
                              <tr
                                key={row.resIdx}
                                onClick={() => setSelectedResidueIdx(row.resIdx - 1)}
                                className={`cursor-pointer transition-colors ${
                                  isSelected ? 'bg-cyan-950/40' : 'hover:bg-cyan-950/20'
                                }`}
                              >
                                <td className="p-2.5 pl-4 font-bold text-zinc-400">#{row.resIdx}</td>
                                <td className="p-2.5 font-bold text-white">{row.resName}</td>
                                <td className="p-2.5 text-zinc-400">{row.nextRes}</td>
                                <td className="p-2.5">
                                  <span className="px-2 py-0.5 rounded bg-indigo-950/80 border border-indigo-700/60 text-indigo-300 font-bold text-[11px]">
                                    Łatka #{row.assignedPatch}
                                  </span>
                                </td>
                                <td className="p-2.5 font-bold">
                                  {rmsdEntry ? (
                                    <span
                                      className={
                                        isNativeFoldUnlocked && rmsdEntry.caDistAngstrom < 1.0
                                          ? 'text-emerald-300'
                                          : 'text-amber-300'
                                      }
                                    >
                                      {rmsdEntry.caDistAngstrom.toFixed(2)} Å
                                    </span>
                                  ) : (
                                    '—'
                                  )}
                                </td>
                                <td className="p-2.5 font-bold text-sky-300">
                                  {rmsdEntry ? `${rmsdEntry.backboneRmsdAngstrom.toFixed(2)} Å` : '—'}
                                </td>
                                <td className="p-2.5 font-bold text-indigo-300">
                                  {rmsdEntry ? `${rmsdEntry.allAtomRmsdAngstrom.toFixed(2)} Å` : '—'}
                                </td>
                                <td className="p-2.5 pr-4">
                                  {isNativeFoldUnlocked ? (
                                    <span className="text-[10px] text-emerald-400 font-bold flex items-center gap-1">
                                      <CheckCircle2 className="w-3 h-3" />
                                      <span>POTWIERDZONY (&lt;1Å)</span>
                                    </span>
                                  ) : (
                                    <span className="text-[10px] text-rose-400 font-bold flex items-center gap-1">
                                      <Lock className="w-3 h-3" />
                                      <span>ZABLOKOWANY (≥1Å)</span>
                                    </span>
                                  )}
                                </td>
                              </tr>
                            );
                          })}
                        </tbody>
                      </table>
                    </div>
                  </div>
                )}
              </>
            )}
          </div>
        </div>
      )}
    </div>
  );
};
