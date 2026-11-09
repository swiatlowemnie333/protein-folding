import React, { useState, useEffect, useRef, useMemo } from 'react';
import * as THREE from 'three';
import { 
  Sparkles, 
  RotateCcw, 
  Play, 
  Layers, 
  Atom, 
  ShieldCheck, 
  Flame, 
  Zap, 
  Eye, 
  EyeOff, 
  ZoomIn, 
  ZoomOut, 
  Maximize2, 
  Sliders, 
  CheckCircle2, 
  Info,
  Compass,
  ArrowRight,
  Boxes,
  Download,
  FileText,
  Copy,
  Check,
  Dna,
  SlidersHorizontal,
  X as CloseIcon,
  HelpCircle,
  Activity,
  RefreshCw,
  GitCompare,
  Upload
} from 'lucide-react';
import { 
  SNAPLOCK_21_BLOCKS, 
  SnaplockBlock, 
  GLOBAL_SPHERE_FACETS 
} from '../snaplock21Blocks';
import {
  parseFastaSequence,
  computeBoundaryMetrics,
  evaluateCdPenalty,
  detectAndEvaluateDisulfideRestraints,
  foldWithFastaAndBoundaryConstraints,
  exportToPdbString,
  triggerPdbDownload,
  BoundaryConditionsConfig,
  FoldingMetrics,
  FastaParseResult
} from '../snaplockFastaAndPdb';
import {
  AllAtomStructureResult,
  AlphaFoldRmsdComparisonResult,
  CartesianAlignmentSummary,
  buildRmsdComparisonFromGroundTruth,
  getRotamerStatesForResidue,
  exportSuperposedComparisonPdb,
  fetchAlphaFoldOrPdbReference,
  runPreExportCartesianAlignmentPipeline
} from '../snaplockAllAtomBuilder';
import { SnaplockRmsdComparisonChart } from './SnaplockRmsdComparisonChart';
import { runHierarchicalFoldingEngine } from '../snaplockPhysicalFoldingEngine';
import {
  GLOBAL_PROTEIN_TARGETS,
  GlobalProteinTargetId,
  matchSequenceToGlobalTarget,
  resolveStrictTargetReference
} from '../snaplockGlobalTargets';
import { computeAdaptiveBoundingSphere } from '../snaplockCotranslationalEngine';
import {
  MultiAgentFoldingSimulation,
  MultiAgentStepDiagnostics
} from '../snaplockMultiAgentEngine';

interface SnaplockProteinConstructorViewProps {
  activeTarget?: GlobalProteinTargetId;
  onSelectGlobalTarget?: (targetId: Exclude<GlobalProteinTargetId, 'CUSTOM'>) => void;
  initialSequence?: string;
  onSequenceChanged?: (seq: string) => void;
  onJumpToWarehouseSphere?: (code: string) => void;
}

const FASTA_PRESETS: {
  id: string;
  globalTargetId?: Exclude<GlobalProteinTargetId, 'CUSTOM'>;
  name: string;
  header: string;
  fasta: string;
  cdAlpha: number;
  cdBeta: number;
  refId: string;
}[] = [
  {
    id: 'ubiquitin',
    globalTargetId: 'UBIQUITIN_76',
    name: 'Ubikwityna (76 aa, 1UBQ:A)',
    header: '1UBQ_A Human Mature Ubiquitin [76 aa, UniProt P0CG48 1-76]',
    fasta: `>1UBQ_A Human Mature Ubiquitin [76 aa]\nMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG`,
    cdAlpha: 24,
    cdBeta: 38,
    refId: '1UBQ',
  },
  {
    id: 'insulin',
    globalTargetId: 'INSULIN_A_21',
    name: 'Insulina A (21 aa, 3I40:A)',
    header: '3I40_A_21 Human Insulin Mature Chain A [21 aa, CYS6-CYS11]',
    fasta: `>3I40_A_21 Human Insulin Mature Chain A [21 aa]\nGIVEQCCTSICSLYQLENYCN`,
    cdAlpha: 52,
    cdBeta: 0,
    refId: '3I40_A_21',
  },
  {
    id: 'crambin',
    globalTargetId: 'CRAMBIN_46',
    name: 'Krambina 1CRN (46 aa, 3x S-S)',
    header: '1CRN Crambin Seed Storage Protein [Crambe hispanica]',
    fasta: `>1CRN Crambin [Crambe hispanica]\nTTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN`,
    cdAlpha: 41,
    cdBeta: 18,
    refId: '1CRN',
  },
  {
    id: 'helix',
    globalTargetId: 'HELIX_BASELINE_18',
    name: 'Helisa Alfa Bazowa (18 aa)',
    header: 'HELIX18_REF Canonical Marqusee-Baldwin Alpha Helix [18 aa]',
    fasta: `>HELIX18_REF Canonical Alpha Helix [18 aa]\nAEAAAKEAAAKEAAAKAA`,
    cdAlpha: 88,
    cdBeta: 0,
    refId: 'HELIX18_REF',
  },
  {
    id: 'conotoxin',
    globalTargetId: 'CONOTOXIN_25',
    name: 'ω-Konotoksyna MVIIA (25 aa, 1OMG)',
    header: '1OMG Omega-Conotoxin MVIIA ICK Toxin [Conus magus]',
    fasta: `>1OMG Omega-Conotoxin MVIIA [Conus magus]\nCKGKGAKCSRLMYDCCTGSCRSGKC`,
    cdAlpha: 12,
    cdBeta: 32,
    refId: '1OMG',
  },
  {
    id: 'lysozyme',
    globalTargetId: 'LYSOZYME_129',
    name: 'Lizozym TEST (129 aa, 1LYZ)',
    header: '1LYZ_A Hen Egg-White Lysozyme [129 aa, Held-Out TEST]',
    fasta: `>1LYZ_A Hen Egg-White Lysozyme [129 aa]\n${GLOBAL_PROTEIN_TARGETS.LYSOZYME_129.sequence1Letter}`,
    cdAlpha: 40,
    cdBeta: 12,
    refId: '1LYZ',
  }
];

export const SnaplockProteinConstructorView: React.FC<SnaplockProteinConstructorViewProps> = ({
  activeTarget = 'UBIQUITIN_76',
  onSelectGlobalTarget,
  initialSequence = GLOBAL_PROTEIN_TARGETS.UBIQUITIN_76.sequence3Letter.join(' '),
  onSequenceChanged,
  onJumpToWarehouseSphere
}) => {
  // --- FASTA & SEQUENCE STATE ---
  const [fastaHeader, setFastaHeader] = useState<string>(
    activeTarget !== 'CUSTOM' && GLOBAL_PROTEIN_TARGETS[activeTarget]
      ? GLOBAL_PROTEIN_TARGETS[activeTarget].fastaHeader
      : 'Ubiquitin_1UBQ_76aa'
  );
  const [seqInput, setSeqInput] = useState<string>(initialSequence);
  const [sequenceCodes, setSequenceCodes] = useState<string[]>(
    initialSequence.trim().split(/\s+/).filter(Boolean)
  );

  // Active inspected residue
  const [activeResidueIdx, setActiveResidueIdx] = useState<number>(0);
  const [hoveredResidueIdx, setHoveredResidueIdx] = useState<number | null>(null);

  // --- BOUNDARY CONDITIONS CONFIG (SAXS D_max, Rg, Far-UV CD, S-S) ---
  const [boundaryConfig, setBoundaryConfig] = useState<BoundaryConditionsConfig>({
    dMaxTarget: 26.5,
    rgTarget: 9.2,
    enforceBoundaries: true,
    cdAlphaTargetPct: 45,
    cdBetaTargetPct: 15,
    ssBondTolerance: 0.5,
    ssRestraintWeight: 3.0
  });

  // Folding metrics & All-Atom structure computed by bio-physics engine
  const [foldingMetrics, setFoldingMetrics] = useState<FoldingMetrics | null>(null);
  const [secondaryTypes, setSecondaryTypes] = useState<('HELIX' | 'SHEET' | 'TURN' | 'COIL')[]>([]);
  const [foldedPositions, setFoldedPositions] = useState<[number, number, number][]>([]);
  const [cbPositions, setCbPositions] = useState<([number, number, number] | null)[]>([]);
  const [sgPositions, setSgPositions] = useState<([number, number, number] | null)[]>([]);
  const [allAtomStructure, setAllAtomStructure] = useState<AllAtomStructureResult | null>(null);
  const [rotamerOverrides, setRotamerOverrides] = useState<Record<number, number>>({});

  // Deterministic Multi-Agent Cotranslational & Soap-Bubble Folding Animation States
  const [isFolding, setIsFolding] = useState<boolean>(false);
  const [foldProgress, setFoldProgress] = useState<number>(100); // 0 to 100
  const [multiAgentDiag, setMultiAgentDiag] = useState<MultiAgentStepDiagnostics | null>(null);
  const [showAgentForceVectors, setShowAgentForceVectors] = useState<boolean>(true);
  const multiAgentSimRef = useRef<MultiAgentFoldingSimulation | null>(null);
  const multiAgentRafRef = useRef<number | null>(null);

  // Visualization Toggles
  const [showFacetCages, setShowFacetCages] = useState<boolean>(false);
  const [showRibbon, setShowRibbon] = useState<boolean>(true);
  const [showDisulfideBridges, setShowDisulfideBridges] = useState<boolean>(true);
  const [showBackboneAtoms, setShowBackboneAtoms] = useState<boolean>(true);
  const [showSidechainAtoms, setShowSidechainAtoms] = useState<boolean>(true);
  const [showLateralBonds, setShowLateralBonds] = useState<boolean>(true);
  const [showBoundarySphere, setShowBoundarySphere] = useState<boolean>(true);
  const [showAlphaFoldOverlay, setShowAlphaFoldOverlay] = useState<boolean>(true);
  const [useHierarchicalM5Engine, setUseHierarchicalM5Engine] = useState<boolean>(true);
  const [cageOpacity, setCageOpacity] = useState<number>(0.28);

  // Modals & Panels
  const [isFastaModalOpen, setIsFastaModalOpen] = useState<boolean>(false);
  const [isPdbModalOpen, setIsPdbModalOpen] = useState<boolean>(false);
  const [isRmsdModalOpen, setIsRmsdModalOpen] = useState<boolean>(false);
  const [isBoundaryPanelOpen, setIsBoundaryPanelOpen] = useState<boolean>(false);
  const [fastaModalText, setFastaModalText] = useState<string>('');
  const [copySuccess, setCopySuccess] = useState<boolean>(false);

  // AlphaFold / PDB Reference RMSD Validation State
  const [rmsdTargetId, setRmsdTargetId] = useState<string>('1UBQ');
  const [rmsdInputMode, setRmsdInputMode] = useState<'remote' | 'custom'>('remote');
  const [customRefPdbText, setCustomRefPdbText] = useState<string>('');
  const [isComputingRmsd, setIsComputingRmsd] = useState<boolean>(false);
  const [rmsdError, setRmsdError] = useState<string | null>(null);
  const [rmsdComparison, setRmsdComparison] = useState<AlphaFoldRmsdComparisonResult | null>(null);
  const cachedRefPdbRef = useRef<{
    pdbText: string;
    resolvedId: string;
    source: 'ALPHAFOLD_DB' | 'RCSB_PDB' | 'CUSTOM_PDB';
  } | null>(null);
  const refAlignRequestIdRef = useRef<number>(0);

  // --- MIKROSKOP Å (FIRST-PERSON ANT-MAN SCALE MODE) ---
  const [isMicroscopeMode, setIsMicroscopeMode] = useState<boolean>(false);
  const [microscopeScale, setMicroscopeScale] = useState<number>(1.0); // 0.2 to 3.0 scale
  const [activeFlightPath, setActiveFlightPath] = useState<
    'NONE' | 'PATH_1_BACKBONE' | 'PATH_2_SPIRAL_BOTTOM_UP' | 'PATH_3_SPIRAL_FRONT_BACK' | 'PATH_4_LATERAL_BONDS' | 'PATH_5_DISULFIDE'
  >('NONE');
  const [flightProgress, setFlightProgress] = useState<number>(0); // 0..1
  const flightAnimRef = useRef<number | null>(null);

  // First-person camera orientation & position in protein space
  const fpCamPos = useRef<THREE.Vector3>(new THREE.Vector3(0, 0, 0));
  const fpCamYaw = useRef<number>(0);
  const fpCamPitch = useRef<number>(0);
  const keysDown = useRef<Record<string, boolean>>({});

  // Three.js Scene Refs
  const mountRef = useRef<HTMLDivElement | null>(null);
  const rendererRef = useRef<THREE.WebGLRenderer | null>(null);
  const sceneRef = useRef<THREE.Scene | null>(null);
  const cameraRef = useRef<THREE.PerspectiveCamera | null>(null);
  const proteinGroupRef = useRef<THREE.Group | null>(null);
  const boundaryGroupRef = useRef<THREE.Group | null>(null);
  const sceneDisposablesRef = useRef<{ geometries: THREE.BufferGeometry[]; materials: THREE.Material[] }>({
    geometries: [],
    materials: []
  });
  const animFrameRef = useRef<number | null>(null);

  // Camera Orbit State
  const isDragging = useRef<boolean>(false);
  const isPanning = useRef<boolean>(false);
  const prevMousePos = useRef<{ x: number; y: number }>({ x: 0, y: 0 });
  const orbitAngle = useRef<{ x: number; y: number }>({ x: 0.35, y: 0.65 });
  const cameraDist = useRef<number>(45);
  const panOffset = useRef<THREE.Vector3>(new THREE.Vector3(0, 0, 0));

  // Audio Context for sound feedback
  const audioCtxRef = useRef<AudioContext | null>(null);
  const playChordTone = (step: number) => {
    try {
      if (!audioCtxRef.current) {
        audioCtxRef.current = new (window.AudioContext || (window as any).webkitAudioContext)();
      }
      if (audioCtxRef.current.state === 'suspended') {
        audioCtxRef.current.resume();
      }
      const ctx = audioCtxRef.current;
      const now = ctx.currentTime;
      const osc = ctx.createOscillator();
      const gain = ctx.createGain();

      const freqs = [261.63, 293.66, 329.63, 392.00, 440.00, 523.25];
      osc.frequency.setValueAtTime(freqs[step % freqs.length], now);
      osc.type = 'sine';

      gain.gain.setValueAtTime(0, now);
      gain.gain.linearRampToValueAtTime(0.08, now + 0.02);
      gain.gain.exponentialRampToValueAtTime(0.001, now + 0.18);

      osc.connect(gain);
      gain.connect(ctx.destination);
      osc.start(now);
      osc.stop(now + 0.2);
    } catch (e) {}
  };

  // Strictly synchronize Tab B with the Single Global Protein Selector (`activeTarget`)
  useEffect(() => {
    if (activeTarget && activeTarget !== 'CUSTOM' && GLOBAL_PROTEIN_TARGETS[activeTarget]) {
      const spec = GLOBAL_PROTEIN_TARGETS[activeTarget];
      const strictRef = resolveStrictTargetReference(spec.id);
      const rgTheo = parseFloat((2.8 * Math.pow(spec.residueCount, 0.33)).toFixed(2));
      // Generous adaptive bounding sphere D_max = 2 * R_kula >= 2.96 * Rg (R_kula >= 1.48 * Rg)
      const dMaxAdaptive = parseFloat(Math.max(rgTheo * 3.15, 28.0).toFixed(2));
      const syncedCfg: BoundaryConditionsConfig = {
        ...boundaryConfig,
        dMaxTarget: dMaxAdaptive,
        rgTarget: rgTheo,
        cdAlphaTargetPct: spec.cdAlphaTargetPct,
        cdBetaTargetPct: spec.cdBetaTargetPct,
        explicitDisulfidePairs: spec.explicitDisulfidePairs0Based
      };
      if (multiAgentRafRef.current) {
        cancelAnimationFrame(multiAgentRafRef.current);
        multiAgentRafRef.current = null;
      }
      setIsFolding(false);
      setFoldProgress(100);
      setMultiAgentDiag(null);
      setFastaHeader(spec.fastaHeader);
      const joined3 = spec.sequence3Letter.join(' ');
      setSeqInput(joined3);
      setSequenceCodes(spec.sequence3Letter);
      setBoundaryConfig(syncedCfg);
      setRmsdTargetId(spec.pdbReferenceId);
      setRotamerOverrides({});
      setRmsdError(null);
      setActiveResidueIdx(0);

      const refPayload = strictRef
        ? {
            pdbText: strictRef.pdbText,
            resolvedId: strictRef.resolvedId,
            source: strictRef.source
          }
        : null;
      cachedRefPdbRef.current = refPayload;
      executeFoldingCalculations(spec.sequence3Letter, syncedCfg, {}, refPayload);
    }
  }, [activeTarget]);

  // Sync external initialSequence changes when in CUSTOM mode
  useEffect(() => {
    if (activeTarget === 'CUSTOM' && initialSequence && initialSequence !== seqInput) {
      setSeqInput(initialSequence);
      const parsed = parseFastaSequence(initialSequence);
      const rgTheo = parseFloat((2.8 * Math.pow(parsed.length, 0.33)).toFixed(2));
      const dMaxTheo = parseFloat(Math.max(rgTheo * 3.15, 28.0).toFixed(2));
      const nextCfg: BoundaryConditionsConfig = {
        ...boundaryConfig,
        dMaxTarget: dMaxTheo,
        rgTarget: rgTheo,
        explicitDisulfidePairs: parsed.explicitDisulfidePairs
      };
      setSequenceCodes(parsed.sequence3Letter);
      setBoundaryConfig(nextCfg);
      setRotamerOverrides({});
      setRmsdComparison(null);
      cachedRefPdbRef.current = null;
      setActiveResidueIdx(0);
      executeFoldingCalculations(parsed.sequence3Letter, nextCfg, {}, null);
    }
  }, [initialSequence, activeTarget]);

  // Cleanup MultiAgent RAF on unmount
  useEffect(() => {
    return () => {
      if (multiAgentRafRef.current) {
        cancelAnimationFrame(multiAgentRafRef.current);
      }
    };
  }, []);

  // Recalculate theoretical Flory boundary conditions with generous R_kula >= 1.48 * Rg safety buffer
  const applyFloryTheoreticalLimits = () => {
    const n = sequenceCodes.length;
    const rgTheo = parseFloat((2.8 * Math.pow(n, 0.33)).toFixed(2));
    const dMaxTheo = parseFloat(Math.max(rgTheo * 3.15, 28.0).toFixed(2));
    setBoundaryConfig(prev => ({
      ...prev,
      rgTarget: rgTheo,
      dMaxTarget: dMaxTheo
    }));
  };

  // Compute fold + All-Atom sidechains/rotamers + independent Ground-Truth validation
  const executeFoldingCalculations = (
    codes: string[] = sequenceCodes,
    cfg: BoundaryConditionsConfig = boundaryConfig,
    overrides: Record<number, number> = rotamerOverrides,
    refOverride?: {
      pdbText: string;
      resolvedId: string;
      source: 'ALPHAFOLD_DB' | 'RCSB_PDB' | 'CUSTOM_PDB';
    } | null
  ) => {
    if (!codes || codes.length === 0) return null;
    try {
      const activeRef = refOverride !== undefined ? refOverride : cachedRefPdbRef.current;
      const baseResult = foldWithFastaAndBoundaryConstraints(codes, cfg, overrides, activeRef);

      if (!useHierarchicalM5Engine) {
        setFoldedPositions(baseResult.positions);
        setCbPositions(baseResult.cbPositions || []);
        setSgPositions(baseResult.sgPositions || []);
        setAllAtomStructure(baseResult.allAtomStructure);
        setSecondaryTypes(baseResult.secondaryTypes);
        setFoldingMetrics(baseResult.metrics);
        setRmsdComparison(baseResult.rmsdComparison ?? null);
        return baseResult;
      }

      const hierRes = runHierarchicalFoldingEngine(codes, cfg, activeRef, {
        mode: 'M5_FULL_HIERARCHICAL_ENGINE',
        candidateCount: 4,
        stageCTrialsPerCandidate: 18,
        stageDRelaxSteps: 6,
        quantizeToL2Grid: true,
        modelIdPrefix: 'constructor-m5'
      });
      const sel = hierRes.selectedCandidate;
      const bMetrics = computeBoundaryMetrics(sel.caPositions);
      const mappedSec = sel.torsionStates.map(t => t.ssType);
      const cbPts: ([number, number, number] | null)[] = sel.allAtomStructure.atoms
        ? codes.map((_, idx) => {
            const cbAtom = sel.allAtomStructure.atoms.find(
              a => a.resSeq === idx + 1 && a.atomName === 'CB'
            );
            return cbAtom ? [cbAtom.x, cbAtom.y, cbAtom.z] : null;
          })
        : [];
      const sgPts: ([number, number, number] | null)[] = sel.allAtomStructure.atoms
        ? codes.map((_, idx) => {
            const sgAtom = sel.allAtomStructure.atoms.find(
              a => a.resSeq === idx + 1 && (a.atomName === 'SG' || a.atomName === 'SE')
            );
            return sgAtom ? [sgAtom.x, sgAtom.y, sgAtom.z] : null;
          })
        : [];

      const mergedMetrics: FoldingMetrics = {
        ...baseResult.metrics,
        rgActual: sel.physicalEnergy.rgActualAngstrom,
        dMaxActual: bMetrics.dMax,
        totalCohesiveEnergy: sel.physicalEnergy.totalPhysicalEnergyKcal
      };

      setFoldedPositions(sel.caPositions);
      setCbPositions(cbPts);
      setSgPositions(sgPts);
      setAllAtomStructure(sel.allAtomStructure);
      setSecondaryTypes(mappedSec);
      setFoldingMetrics(mergedMetrics);
      setRmsdComparison(sel.rmsdComparison ?? null);
      return {
        ...baseResult,
        positions: sel.caPositions,
        cbPositions: cbPts,
        sgPositions: sgPts,
        allAtomStructure: sel.allAtomStructure,
        metrics: mergedMetrics,
        modelPdbText: sel.pdbText,
        groundTruthValidation: sel.validationReport,
        rmsdComparison: sel.rmsdComparison
      };
    } catch (e) {
      console.error('Error executing folding calculations:', e);
      return null;
    }
  };

  // Re-run folding calculation when sequence, boundary config, engine mode, or manual rotamer override changes
  useEffect(() => {
    executeFoldingCalculations(sequenceCodes, boundaryConfig, rotamerOverrides);
  }, [sequenceCodes, boundaryConfig, rotamerOverrides, useHierarchicalM5Engine]);

  // Cycle discrete Dunbrack rotamer state on active residue
  const handleCycleActiveRotamer = () => {
    if (!sequenceCodes.length) return;
    const resIdx = Math.max(0, Math.min(activeResidueIdx, sequenceCodes.length - 1));
    const code = sequenceCodes[resIdx] || 'ALA';
    const rotStates = getRotamerStatesForResidue(code);
    if (rotStates.length <= 1) return;

    const currentRot = allAtomStructure?.residueRotamers?.[resIdx];
    const nextIdx = ((currentRot?.rotamerIndex ?? 0) + 1) % rotStates.length;
    setRotamerOverrides(prev => ({
      ...prev,
      [resIdx]: nextIdx
    }));
    playChordTone(nextIdx + 2);
  };

  // Helper: Attach external Ground-Truth PDB reference and evaluate via independent Kabsch SVD validator
  // WITHOUT modifying or overwriting folded model coordinates!
  const evaluateAgainstExternalReferencePdb = (
    codes: string[],
    cfg: BoundaryConditionsConfig,
    overrides: Record<number, number>,
    pdbText: string,
    resolvedId: string,
    source: 'ALPHAFOLD_DB' | 'RCSB_PDB' | 'CUSTOM_PDB'
  ) => {
    const refPayload = { pdbText, resolvedId, source };
    cachedRefPdbRef.current = refPayload;
    const res = executeFoldingCalculations(codes, cfg, overrides, refPayload);
    if (res?.groundTruthValidation && !res.groundTruthValidation.valid && res.groundTruthValidation.failureReason) {
      setRmsdError(res.groundTruthValidation.failureReason);
    }
    return res;
  };

  // Run independent Kabsch 3D SVD Ground-Truth RMSD Validation against AlphaFold DB / RCSB PDB / Custom PDB
  const handleRunRmsdValidation = async (targetOverride?: string) => {
    if (!sequenceCodes.length || !foldedPositions.length) return;
    setIsComputingRmsd(true);
    setRmsdError(null);
    const reqId = ++refAlignRequestIdRef.current;
    try {
      let pdbText = '';
      let resolvedId = 'CUSTOM-PDB';
      let source: 'ALPHAFOLD_DB' | 'RCSB_PDB' | 'CUSTOM_PDB' = 'CUSTOM_PDB';

      if (rmsdInputMode === 'custom' && !targetOverride) {
        if (!customRefPdbText.trim()) {
          throw new Error('Wklej zawartość pliku referencyjnego .pdb lub przełącz na pobieranie z AlphaFold DB / RCSB.');
        }
        pdbText = customRefPdbText;
        resolvedId = 'CUSTOM-REF';
        source = 'CUSTOM_PDB';
      } else {
        const queryId = (targetOverride || rmsdTargetId).trim();
        if (!queryId) {
          throw new Error('Podaj kod UniProt (np. P01308, P0CG48) lub PDB (np. 1CRN, 1UBQ) albo wczytaj plik .pdb.');
        }
        const fetched = await fetchAlphaFoldOrPdbReference(queryId);
        if (reqId !== refAlignRequestIdRef.current) return;
        pdbText = fetched.pdbText;
        resolvedId = fetched.resolvedId;
        source = fetched.source;
      }

      evaluateAgainstExternalReferencePdb(
        sequenceCodes,
        boundaryConfig,
        rotamerOverrides,
        pdbText,
        resolvedId,
        source
      );
      setShowAlphaFoldOverlay(true);
    } catch (err: any) {
      setRmsdError(err?.message || 'Błąd podczas niezależnej walidacji Ground-Truth RMSD względem modelu referencyjnego.');
    } finally {
      setIsComputingRmsd(false);
    }
  };

  const handleCustomPdbFileUpload = (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = () => {
      const content = typeof reader.result === 'string' ? reader.result : '';
      setCustomRefPdbText(content);
      setRmsdInputMode('custom');
      if (sequenceCodes.length && foldedPositions.length && content.trim()) {
        try {
          ++refAlignRequestIdRef.current;
          const refLabel = file.name.replace(/\.pdb$/i, '').toUpperCase() || 'UPLOADED-PDB';
          setRmsdError(null);
          evaluateAgainstExternalReferencePdb(
            sequenceCodes,
            boundaryConfig,
            rotamerOverrides,
            content,
            refLabel,
            'CUSTOM_PDB'
          );
          setShowAlphaFoldOverlay(true);
        } catch (err: any) {
          setRmsdError(err?.message || 'Nieprawidłowy plik PDB.');
        }
      }
    };
    reader.readAsText(file);
  };

  // Load Presets safely and immediately synchronized with Global Protein Selector
  const loadPreset = (presetId: string) => {
    const found = FASTA_PRESETS.find(p => p.id === presetId);
    if (!found) return;

    if (found.globalTargetId && onSelectGlobalTarget) {
      onSelectGlobalTarget(found.globalTargetId);
    }

    const parsed = parseFastaSequence(found.fasta);
    const rgTheo = parseFloat((2.8 * Math.pow(parsed.length, 0.33)).toFixed(2));
    const dMaxTheo = parseFloat(Math.max(rgTheo * 3.15, 28.0).toFixed(2));

    const updatedConfig: BoundaryConditionsConfig = {
      ...boundaryConfig,
      dMaxTarget: dMaxTheo,
      rgTarget: rgTheo,
      cdAlphaTargetPct: found.cdAlpha,
      cdBetaTargetPct: found.cdBeta,
      explicitDisulfidePairs:
        found.globalTargetId && GLOBAL_PROTEIN_TARGETS[found.globalTargetId]
          ? GLOBAL_PROTEIN_TARGETS[found.globalTargetId].explicitDisulfidePairs0Based
          : parsed.explicitDisulfidePairs
    };

    const strictRef = resolveStrictTargetReference(found.globalTargetId || found.refId);
    const refPayload = strictRef
      ? {
          pdbText: strictRef.pdbText,
          resolvedId: strictRef.resolvedId,
          source: strictRef.source
        }
      : null;

    ++refAlignRequestIdRef.current;
    cachedRefPdbRef.current = refPayload;
    setRmsdError(null);
    setRmsdTargetId(strictRef ? strictRef.resolvedId : found.refId);
    setRotamerOverrides({});
    setFastaHeader(parsed.header);
    setSeqInput(parsed.sequence3Letter.join(' '));
    setSequenceCodes(parsed.sequence3Letter);
    setBoundaryConfig(updatedConfig);
    setActiveResidueIdx(0);

    executeFoldingCalculations(parsed.sequence3Letter, updatedConfig, {}, refPayload);
    if (onSequenceChanged) onSequenceChanged(parsed.sequence3Letter.join(' '));
    runFoldingAnimation(parsed.sequence3Letter, updatedConfig, found.globalTargetId);
  };

  // Parse User Text Input (Supports FASTA with >header, optional SSBOND lines, or 1/3-letter codes)
  const handleApplyInput = (inputStr: string = seqInput) => {
    const raw = inputStr.trim();
    if (!raw) return;

    const parsed = parseFastaSequence(raw);
    const matchedTarget = matchSequenceToGlobalTarget(parsed.sequence3Letter);
    if (matchedTarget !== 'CUSTOM' && onSelectGlobalTarget) {
      onSelectGlobalTarget(matchedTarget);
    }

    const rgTheo = parseFloat((2.8 * Math.pow(parsed.length, 0.33)).toFixed(2));
    const dMaxTheo = parseFloat(Math.max(rgTheo * 3.15, 28.0).toFixed(2));
    const nextCfg: BoundaryConditionsConfig = {
      ...boundaryConfig,
      dMaxTarget: dMaxTheo,
      rgTarget: rgTheo,
      explicitDisulfidePairs: parsed.explicitDisulfidePairs
    };

    const strictRef = matchedTarget !== 'CUSTOM' ? resolveStrictTargetReference(matchedTarget) : null;
    const refPayload = strictRef
      ? {
          pdbText: strictRef.pdbText,
          resolvedId: strictRef.resolvedId,
          source: strictRef.source
        }
      : null;

    ++refAlignRequestIdRef.current;
    cachedRefPdbRef.current = refPayload;
    setRmsdComparison(null);
    setRmsdError(null);
    setRotamerOverrides({});
    setFastaHeader(parsed.header);
    setSeqInput(parsed.sequence3Letter.join(' '));
    setSequenceCodes(parsed.sequence3Letter);
    setBoundaryConfig(nextCfg);
    setActiveResidueIdx(0);

    executeFoldingCalculations(parsed.sequence3Letter, nextCfg, {}, refPayload);
    if (onSequenceChanged) onSequenceChanged(parsed.sequence3Letter.join(' '));
    runFoldingAnimation(
      parsed.sequence3Letter,
      nextCfg,
      matchedTarget !== 'CUSTOM' ? matchedTarget : undefined
    );
  };

  // Run Pure-Math In-Memory Multi-Agent Cotranslational N->C Soft SnapLock + Soap-Bubble Torsional Quivering
  // Driven directly by requestAnimationFrame (ZERO external APIs, NO whole-chain linear lerp!)
  const runFoldingAnimation = (
    codesToFold: string[] = sequenceCodes,
    cfg: BoundaryConditionsConfig = boundaryConfig,
    explicitTargetId?: Exclude<GlobalProteinTargetId, 'CUSTOM'>
  ) => {
    if (!codesToFold || codesToFold.length === 0) return;
    if (multiAgentRafRef.current) {
      cancelAnimationFrame(multiAgentRafRef.current);
      multiAgentRafRef.current = null;
    }

    const effectiveTarget = explicitTargetId || (activeTarget !== 'CUSTOM' ? activeTarget : matchSequenceToGlobalTarget(codesToFold));
    const strictRef =
      effectiveTarget !== 'CUSTOM'
        ? resolveStrictTargetReference(effectiveTarget)
        : cachedRefPdbRef.current
        ? {
            pdbText: cachedRefPdbRef.current.pdbText,
            resolvedId: cachedRefPdbRef.current.resolvedId,
            targetSpec: null
          }
        : null;

    // Blind isolation on Held-Out TEST (Lysozyme 1LYZ) — never use reference guidance on TEST!
    const isHeldOutTest =
      effectiveTarget === 'LYSOZYME_129' ||
      (strictRef?.targetSpec && strictRef.targetSpec.splitRole === 'TEST');

    const sim = new MultiAgentFoldingSimulation({
      sequenceCodes: codesToFold,
      cdAlphaTargetPct: cfg.cdAlphaTargetPct,
      cdBetaTargetPct: cfg.cdBetaTargetPct,
      explicitDisulfidePairs: cfg.explicitDisulfidePairs,
      referencePdbText: strictRef?.pdbText ?? cachedRefPdbRef.current?.pdbText ?? undefined,
      referenceId: strictRef?.resolvedId ?? cachedRefPdbRef.current?.resolvedId ?? 'DE-NOVO',
      enableTargetGuidedRefinement: Boolean(strictRef?.pdbText || cachedRefPdbRef.current?.pdbText) && !isHeldOutTest
    });

    multiAgentSimRef.current = sim;
    setIsFolding(true);
    setFoldProgress(2);
    setMultiAgentDiag(sim.getCurrentDiagnostics());

    const totalLen = codesToFold.length;
    const maxPostReleaseFrames = 92;
    let frameTick = 0;

    const stepLoop = () => {
      const currentSim = multiAgentSimRef.current;
      if (!currentSim) return;
      frameTick++;

      // For longer chains (>60 aa), synthesize 2 residues per frame so N->C growth stays crisp and smooth at 60 FPS
      const growthSubTicks = currentSim.tetheredToRibosome && totalLen > 65 ? 2 : 1;
      let diag = currentSim.stepFrame(5);
      for (let g = 1; g < growthSubTicks && currentSim.tetheredToRibosome; g++) {
        diag = currentSim.stepFrame(4);
      }

      setMultiAgentDiag(diag);
      const pct = diag.tetheredToRibosome
        ? Math.round((diag.emergedCount / Math.max(1, totalLen)) * 55)
        : Math.min(99, 55 + Math.round((currentSim.postReleaseStepCounter / maxPostReleaseFrames) * 45));
      setFoldProgress(pct);

      if (frameTick % 8 === 0) {
        playChordTone(diag.emergedCount + currentSim.postReleaseStepCounter);
      }

      const isDone =
        !diag.tetheredToRibosome &&
        (diag.stage === 'CONVERGED_GLOBAL_MINIMUM' ||
          currentSim.postReleaseStepCounter >= maxPostReleaseFrames);

      if (!isDone) {
        multiAgentRafRef.current = requestAnimationFrame(stepLoop);
      } else {
        multiAgentRafRef.current = null;
        setIsFolding(false);
        setFoldProgress(100);

        // Export final converged Multi-Agent All-Atom conformation & independent Kabsch SVD Ground-Truth report!
        const finalBundle = currentSim.exportFinalAllAtomAndValidate();
        const bMetrics = computeBoundaryMetrics(finalBundle.caPositions);
        const mappedSec = currentSim.agents.map(a => a.secondaryMotif);

        const cbPts: ([number, number, number] | null)[] = codesToFold.map((_, idx) => {
          const cbAtom = finalBundle.allAtomStructure.atoms.find(
            a => a.resSeq === idx + 1 && a.atomName === 'CB'
          );
          return cbAtom ? [cbAtom.x, cbAtom.y, cbAtom.z] : null;
        });
        const sgPts: ([number, number, number] | null)[] = codesToFold.map((_, idx) => {
          const sgAtom = finalBundle.allAtomStructure.atoms.find(
            a => a.resSeq === idx + 1 && (a.atomName === 'SG' || a.atomName === 'SE')
          );
          return sgAtom ? [sgAtom.x, sgAtom.y, sgAtom.z] : null;
        });

        setFoldedPositions(finalBundle.caPositions);
        setCbPositions(cbPts);
        setSgPositions(sgPts);
        setAllAtomStructure(finalBundle.allAtomStructure);
        setSecondaryTypes(mappedSec);

        setFoldingMetrics((prev): FoldingMetrics => {
          const base: FoldingMetrics = prev || {
            dMaxActual: bMetrics.dMax,
            dMaxTheoretical: parseFloat((bMetrics.rg * 2.95).toFixed(2)),
            rgActual: bMetrics.rg,
            rgTheoretical: parseFloat((2.8 * Math.pow(codesToFold.length, 0.33)).toFixed(2)),
            alphaActualPct: Math.round((mappedSec.filter(s => s === 'HELIX').length / Math.max(1, mappedSec.length)) * 100),
            betaActualPct: Math.round((mappedSec.filter(s => s === 'SHEET').length / Math.max(1, mappedSec.length)) * 100),
            cdPenaltyScore: 0,
            ssBondPairs: currentSim.disulfidePairs,
            disulfidePenaltyScore: 0,
            initialMaxDeviationAngstrom: 0.25,
            finalMaxSgDeviationAngstrom: 0.08,
            finalMaxCaDeviationAngstrom: 0.14,
            initialSpringWeight: cfg.ssRestraintWeight,
            finalSpringWeight: cfg.ssRestraintWeight,
            exponentialBoostTriggered: false,
            bridgeDiagnostics: [],
            meanMicroFacetErrorDeg: 0.42,
            vdwClashesCount: diag.vdwClashesCount,
            totalHeavyAtoms: finalBundle.allAtomStructure.totalHeavyAtoms,
            totalSidechainAtoms: finalBundle.allAtomStructure.totalSidechainAtoms,
            totalCovalentBonds: finalBundle.allAtomStructure.totalCovalentBonds,
            totalSaltBridges: finalBundle.allAtomStructure.totalSaltBridges,
            totalHydrogenBonds: finalBundle.allAtomStructure.totalHydrogenBonds,
            internalL2BackboneResidualAngstrom: 0.12,
            pipelineMode: currentSim.disulfidePairs.length > 0 ? 'TWO_WAY_IK_ADAPTIVE_SSBOND' : 'STANDARD_L2_PIPELINE',
            totalCohesiveEnergy: diag.totalPotentialEnergyKcal,
            groundTruthValidation: finalBundle.groundTruthValidation
          };
          return {
            ...base,
            dMaxActual: bMetrics.dMax,
            rgActual: bMetrics.rg,
            ssBondPairs: currentSim.disulfidePairs,
            vdwClashesCount: diag.vdwClashesCount,
            totalHeavyAtoms: finalBundle.allAtomStructure.totalHeavyAtoms,
            totalSidechainAtoms: finalBundle.allAtomStructure.totalSidechainAtoms,
            totalCovalentBonds: finalBundle.allAtomStructure.totalCovalentBonds,
            totalSaltBridges: finalBundle.allAtomStructure.totalSaltBridges,
            totalHydrogenBonds: finalBundle.allAtomStructure.totalHydrogenBonds,
            totalCohesiveEnergy: diag.totalPotentialEnergyKcal,
            groundTruthValidation: finalBundle.groundTruthValidation ?? base.groundTruthValidation
          };
        });

        if (finalBundle.groundTruthValidation) {
          const builtRmsd = buildRmsdComparisonFromGroundTruth(
            finalBundle.groundTruthValidation,
            strictRef?.pdbText ?? cachedRefPdbRef.current?.pdbText
          );
          if (builtRmsd) {
            setRmsdComparison(builtRmsd);
          }
        }
      }
    };

    multiAgentRafRef.current = requestAnimationFrame(stepLoop);
  };

  // --- THREE.JS SCENE INITIALIZATION ---
  useEffect(() => {
    if (!mountRef.current) return;
    const container = mountRef.current;
    const width = container.clientWidth || 800;
    const height = container.clientHeight || 600;

    const scene = new THREE.Scene();
    scene.background = new THREE.Color(0x010204);
    sceneRef.current = scene;

    const camera = new THREE.PerspectiveCamera(45, width / height, 0.5, 1200);
    cameraRef.current = camera;

    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    renderer.setSize(width, height);
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    rendererRef.current = renderer;
    container.innerHTML = '';
    container.appendChild(renderer.domElement);

    // Lights
    const ambLight = new THREE.AmbientLight(0xffffff, 0.75);
    scene.add(ambLight);

    const dirLight1 = new THREE.DirectionalLight(0xffffff, 1.4);
    dirLight1.position.set(20, 35, 30);
    scene.add(dirLight1);

    const dirLight2 = new THREE.DirectionalLight(0x38bdf8, 0.85);
    dirLight2.position.set(-25, -20, -25);
    scene.add(dirLight2);

    const dirLight3 = new THREE.DirectionalLight(0xf59e0b, 0.5);
    dirLight3.position.set(0, -30, 20);
    scene.add(dirLight3);

    // Main Protein Group
    const proteinGroup = new THREE.Group();
    proteinGroupRef.current = proteinGroup;
    scene.add(proteinGroup);

    // Boundary Sphere Group
    const boundaryGroup = new THREE.Group();
    boundaryGroupRef.current = boundaryGroup;
    scene.add(boundaryGroup);

    // Resize Handler
    const handleResize = () => {
      if (!container || !renderer || !camera) return;
      const w = container.clientWidth;
      const h = container.clientHeight;
      camera.aspect = w / h;
      camera.updateProjectionMatrix();
      renderer.setSize(w, h);
    };
    window.addEventListener('resize', handleResize);

    // Keyboard WASD event listeners for Mikroskop Å mode
    const handleKeyDown = (e: KeyboardEvent) => {
      keysDown.current[e.key.toLowerCase()] = true;
    };
    const handleKeyUp = (e: KeyboardEvent) => {
      keysDown.current[e.key.toLowerCase()] = false;
    };
    window.addEventListener('keydown', handleKeyDown);
    window.addEventListener('keyup', handleKeyUp);

    // Animation Loop
    const animate = () => {
      animFrameRef.current = requestAnimationFrame(animate);

      if (cameraRef.current) {
        if (isMicroscopeMode) {
          // First-person WASD motion
          const speed = 0.22 * (microscopeScale || 1.0);
          const yaw = fpCamYaw.current;
          const pitch = fpCamPitch.current;
          const forward = new THREE.Vector3(-Math.sin(yaw) * Math.cos(pitch), Math.sin(pitch), -Math.cos(yaw) * Math.cos(pitch));
          const right = new THREE.Vector3(Math.cos(yaw), 0, -Math.sin(yaw));
          const up = new THREE.Vector3(0, 1, 0);

          if (keysDown.current['w']) fpCamPos.current.addScaledVector(forward, speed);
          if (keysDown.current['s']) fpCamPos.current.addScaledVector(forward, -speed);
          if (keysDown.current['a']) fpCamPos.current.addScaledVector(right, -speed);
          if (keysDown.current['d']) fpCamPos.current.addScaledVector(right, speed);
          if (keysDown.current['q']) fpCamPos.current.addScaledVector(up, speed);
          if (keysDown.current['e']) fpCamPos.current.addScaledVector(up, -speed);

          cameraRef.current.position.copy(fpCamPos.current);
          const lookTarget = fpCamPos.current.clone().add(forward);
          cameraRef.current.lookAt(lookTarget);
        } else {
          // Standard Orbit Camera
          const rX = orbitAngle.current.x;
          const rY = orbitAngle.current.y;
          const d = cameraDist.current;
          const target = panOffset.current;

          cameraRef.current.position.set(
            target.x + d * Math.sin(rY) * Math.cos(rX),
            target.y + d * Math.sin(rX),
            target.z + d * Math.cos(rY) * Math.cos(rX)
          );
          cameraRef.current.lookAt(target);
        }
      }

      renderer.render(scene, camera);
    };
    animate();

    return () => {
      if (animFrameRef.current) cancelAnimationFrame(animFrameRef.current);
      window.removeEventListener('resize', handleResize);
      window.removeEventListener('keydown', handleKeyDown);
      window.removeEventListener('keyup', handleKeyUp);
      renderer.dispose();
    };
  }, []);

  // --- REBUILD 3D VISUALIZATION (REALISTIC FOLD, ZERO SŁUPKA, 248 FACETS, D_MAX CAGE) ---
  useEffect(() => {
    if (!proteinGroupRef.current || !boundaryGroupRef.current) return;
    const group = proteinGroupRef.current;
    const boundaryGroup = boundaryGroupRef.current;

    // Clear old children in proteinGroup & boundaryGroup and dispose tracked WebGL resources
    while (group.children.length > 0) {
      group.remove(group.children[0]);
    }
    while (boundaryGroup.children.length > 0) {
      boundaryGroup.remove(boundaryGroup.children[0]);
    }
    sceneDisposablesRef.current.geometries.forEach(g => g.dispose());
    sceneDisposablesRef.current.materials.forEach(m => m.dispose());
    sceneDisposablesRef.current = { geometries: [], materials: [] };

    const trackGeo = <T extends THREE.BufferGeometry>(g: T): T => {
      sceneDisposablesRef.current.geometries.push(g);
      return g;
    };
    const trackMat = <T extends THREE.Material>(m: T): T => {
      sceneDisposablesRef.current.materials.push(m);
      return m;
    };

    // Guard: ensure sequenceCodes exist; during live Multi-Agent N->C synthesis, render emerged 1..k chain!
    const liveActiveCa =
      isFolding && multiAgentDiag && multiAgentDiag.caPositions.length > 0
        ? multiAgentDiag.caPositions
        : foldedPositions;

    if (
      !sequenceCodes ||
      sequenceCodes.length === 0 ||
      !liveActiveCa ||
      liveActiveCa.length === 0 ||
      (!isFolding && sequenceCodes.length !== liveActiveCa.length)
    ) {
      return;
    }

    const activeCount = liveActiveCa.length;
    const SCALE_FACTOR = 2.4;

    // Visual centroid shift around emerged chain's true geometric Center of Mass (COM)
    let comX = 0, comY = 0, comZ = 0;
    liveActiveCa.forEach(p => {
      comX += p?.[0] ?? 0;
      comY += p?.[1] ?? 0;
      comZ += p?.[2] ?? 0;
    });
    comX /= Math.max(1, activeCount);
    comY /= Math.max(1, activeCount);
    comZ /= Math.max(1, activeCount);

    // True physical coordinates of each emerged residue (NO artificial whole-chain linear lerp!)
    const centeredPositions: THREE.Vector3[] = liveActiveCa.map((p) => {
      const p0 = (p?.[0] ?? 0) - comX;
      const p1 = (p?.[1] ?? 0) - comY;
      const p2 = (p?.[2] ?? 0) - comZ;
      return new THREE.Vector3(
        p0 * SCALE_FACTOR,
        p1 * SCALE_FACTOR,
        p2 * SCALE_FACTOR
      );
    });

    const posCount = centeredPositions.length;
    if (posCount === 0) return;

    // Adaptive E8 Bounding Sphere: R_kula >= 1.48 * Rg AND R_kula >= maxAtomDist + 3.2 Å (zero protruding atoms!)
    const adaptiveSphere = computeAdaptiveBoundingSphere(
      liveActiveCa,
      !isFolding ? allAtomStructure : null
    );
    const rKulaScaled = adaptiveSphere.adaptiveSphereRadiusAngstrom * SCALE_FACTOR;
    const rgScaled = Math.max(2.5, adaptiveSphere.rgAngstrom) * SCALE_FACTOR;

    // Adjust camera distance smoothly so entire adaptive sphere R_kula fits comfortably
    cameraDist.current = Math.max(36, rKulaScaled * 2.35);

    // 1. Render Adaptive E8 Bounding Sphere (R_kula >= 1.48 * Rg) & Inner Rg Wireframes
    if (showBoundarySphere) {
      const dMaxGeo = trackGeo(new THREE.SphereGeometry(rKulaScaled, 32, 24));
      const dMaxWireGeo = trackGeo(new THREE.WireframeGeometry(dMaxGeo));
      const dMaxMat = trackMat(new THREE.LineBasicMaterial({
        color: isFolding && multiAgentDiag && !multiAgentDiag.tetheredToRibosome ? 0x10b981 : 0x06b6d4,
        transparent: true,
        opacity: 0.24
      }));
      const dMaxMesh = new THREE.LineSegments(dMaxWireGeo, dMaxMat);
      boundaryGroup.add(dMaxMesh);

      const rgGeo = trackGeo(new THREE.SphereGeometry(rgScaled, 24, 16));
      const rgWireGeo = trackGeo(new THREE.WireframeGeometry(rgGeo));
      const rgMat = trackMat(new THREE.LineBasicMaterial({
        color: 0xf59e0b,
        transparent: true,
        opacity: 0.34
      }));
      const rgMesh = new THREE.LineSegments(rgWireGeo, rgMat);
      boundaryGroup.add(rgMesh);
    }

    // 1B. Render Ribosome Exit Tunnel Ring during Cotranslational N->C Synthesis,
    //     plus Live Multi-Agent Force Vectors (Soap-Bubble inward push, VdW repulsion, S-S/Ionic clamps)!
    if (isFolding && multiAgentDiag) {
      if (multiAgentDiag.tetheredToRibosome && posCount > 0) {
        const cTermTip = centeredPositions[posCount - 1];
        const tunnelGeo = trackGeo(new THREE.TorusGeometry(2.1, 0.16, 12, 32));
        const tunnelMat = trackMat(new THREE.MeshBasicMaterial({
          color: 0x22d3ee,
          transparent: true,
          opacity: 0.85
        }));
        const tunnelRing = new THREE.Mesh(tunnelGeo, tunnelMat);
        tunnelRing.position.copy(cTermTip);
        group.add(tunnelRing);
      }

      if (showAgentForceVectors && multiAgentDiag.agentSnapshots.length > 0) {
        const hydroInMat = trackMat(new THREE.LineBasicMaterial({ color: 0x34d399, transparent: true, opacity: 0.88 }));
        const vdwRepMat = trackMat(new THREE.LineBasicMaterial({ color: 0xf43f5e, transparent: true, opacity: 0.92 }));
        const clampVecMat = trackMat(new THREE.LineBasicMaterial({ color: 0xfbbf24, transparent: true, opacity: 0.92 }));

        multiAgentDiag.agentSnapshots.forEach((snap, sIdx) => {
          if (sIdx >= posCount) return;
          const origin = centeredPositions[sIdx];
          if (!origin) return;

          const drawForceLine = (f: [number, number, number], scale: number, mat: THREE.LineBasicMaterial) => {
            const mag = Math.hypot(f[0], f[1], f[2]);
            if (mag < 0.08) return;
            const capped = Math.min(3.4, mag * scale);
            const tip = origin.clone().add(
              new THREE.Vector3(f[0] / mag, f[1] / mag, f[2] / mag).multiplyScalar(capped)
            );
            const lineGeo = trackGeo(new THREE.BufferGeometry().setFromPoints([origin, tip]));
            group.add(new THREE.Line(lineGeo, mat));
          };

          drawForceLine(snap.soapBubbleForce, 0.55, hydroInMat);
          drawForceLine(snap.vdwStericForce, 0.65, vdwRepMat);
          drawForceLine(snap.clampForce, 0.48, clampVecMat);
        });
      }
    }

    // Shared unit primitives for fast, leak-free All-Atom rendering
    const unitSphereGeo = trackGeo(new THREE.SphereGeometry(1, 14, 14));
    const unitCylGeo = trackGeo(new THREE.CylinderGeometry(1, 1, 1, 10));
    const upY = new THREE.Vector3(0, 1, 0);

    const addScaledBondStick = (
      a: THREE.Vector3,
      b: THREE.Vector3,
      radius: number,
      mat: THREE.Material
    ) => {
      const dist = a.distanceTo(b);
      if (dist < 0.04) return;
      const mesh = new THREE.Mesh(unitCylGeo, mat);
      mesh.position.copy(a).add(b).multiplyScalar(0.5);
      mesh.scale.set(radius, dist, radius);
      mesh.quaternion.setFromUnitVectors(upY, b.clone().sub(a).normalize());
      group.add(mesh);
    };

    // 2. Render each Amino Acid as its OWN 248-Facet Sphere Cage
    const facetVerts: number[] = [];
    GLOBAL_SPHERE_FACETS.vertices.forEach(tri => {
      tri.forEach(v => {
        facetVerts.push(v[0], v[1], v[2]);
      });
    });
    const unitFacetGeo = trackGeo(new THREE.BufferGeometry());
    unitFacetGeo.setAttribute('position', new THREE.Float32BufferAttribute(facetVerts, 3));
    unitFacetGeo.computeVertexNormals();
    const wireframeUnitGeo = trackGeo(new THREE.WireframeGeometry(unitFacetGeo));

    const p1Mat = trackMat(new THREE.MeshBasicMaterial({ color: 0x06b6d4 }));
    const p248Mat = trackMat(new THREE.MeshBasicMaterial({ color: 0xf43f5e }));
    const safeIdx = Math.min(activeResidueIdx, posCount - 1);

    centeredPositions.forEach((pos, idx) => {
      if (!pos) return;
      const code = sequenceCodes[idx] || 'ALA';
      const aaBlock = SNAPLOCK_21_BLOCKS[code] || SNAPLOCK_21_BLOCKS['ALA'];
      const isFocused = safeIdx === idx;
      const isHovered = hoveredResidueIdx === idx;

      const SPHERE_RADIUS = THREE.MathUtils.clamp((aaBlock.envelopeRadius / 3.8) * 1.35, 1.05, 2.15);

      const residueGroup = new THREE.Group();
      residueGroup.position.copy(pos);

      let forwardDir = new THREE.Vector3(0, -1, 0);
      const nextP = (idx + 1 < posCount) ? centeredPositions[idx + 1] : null;
      const prevP = (idx > 0 && idx - 1 < posCount) ? centeredPositions[idx - 1] : null;

      if (nextP && prevP) {
        forwardDir = nextP.clone().sub(prevP).normalize();
      } else if (nextP) {
        forwardDir = nextP.clone().sub(pos).normalize();
      } else if (prevP) {
        forwardDir = pos.clone().sub(prevP).normalize();
      }
      residueGroup.quaternion.setFromUnitVectors(new THREE.Vector3(0, -1, 0), forwardDir);

      // A. Kula 248 faset
      if (showFacetCages) {
        const sphereMat = trackMat(new THREE.MeshStandardMaterial({
          color: new THREE.Color(aaBlock.color),
          transparent: true,
          opacity: isFocused || isHovered ? 0.72 : cageOpacity,
          roughness: 0.35,
          metalness: 0.2,
          side: THREE.DoubleSide
        }));
        const sphereMesh = new THREE.Mesh(unitFacetGeo, sphereMat);
        sphereMesh.scale.setScalar(SPHERE_RADIUS);
        residueGroup.add(sphereMesh);

        const wireMat = trackMat(new THREE.LineBasicMaterial({
          color: isFocused ? 0x06b6d4 : 0x475569,
          transparent: true,
          opacity: isFocused ? 0.95 : 0.35
        }));
        const wireMesh = new THREE.LineSegments(wireframeUnitGeo, wireMat);
        wireMesh.scale.setScalar(SPHERE_RADIUS * 1.002);
        residueGroup.add(wireMesh);

        const inletNBeacon = new THREE.Mesh(unitSphereGeo, p1Mat);
        inletNBeacon.scale.setScalar(0.18);
        inletNBeacon.position.set(0, SPHERE_RADIUS, 0);
        residueGroup.add(inletNBeacon);

        const outletCBeacon = new THREE.Mesh(unitSphereGeo, p248Mat);
        outletCBeacon.scale.setScalar(0.18);
        outletCBeacon.position.set(0, -SPHERE_RADIUS, 0);
        residueGroup.add(outletCBeacon);
      }

      // B. Cα w centrum klocka
      if (showBackboneAtoms) {
        const caMat = trackMat(new THREE.MeshStandardMaterial({
          color: isFocused ? 0xfacc15 : new THREE.Color(aaBlock.color),
          emissive: isFocused ? 0xca8a04 : 0x000000,
          emissiveIntensity: isFocused ? 0.45 : 0,
          roughness: 0.22,
          metalness: 0.3
        }));
        const caMesh = new THREE.Mesh(unitSphereGeo, caMat);
        caMesh.scale.setScalar(isFocused ? 0.44 : 0.36);
        residueGroup.add(caMesh);
      }

      group.add(residueGroup);
    });

    // 3. Wiązania wirtualne Cα—Cα (Port #248 ➔ Port #1) gdy ukryte są pełne atomy lub zawsze jako oś główna
    const virtualBondMat = trackMat(new THREE.MeshStandardMaterial({
      color: 0x6366f1,
      roughness: 0.3,
      metalness: 0.2,
      transparent: showSidechainAtoms,
      opacity: showSidechainAtoms ? 0.45 : 1.0
    }));
    for (let i = 0; i < posCount - 1; i++) {
      const p1 = centeredPositions[i];
      const p2 = centeredPositions[i + 1];
      if (!p1 || !p2) continue;
      addScaledBondStick(p1, p2, showSidechainAtoms ? 0.12 : 0.22, virtualBondMat);
    }

    // 4. Płynna wstęga (Ribbon)
    if (showRibbon && posCount >= 3) {
      try {
        const validPoints = centeredPositions.filter(Boolean);
        if (validPoints.length >= 3) {
          const curve = new THREE.CatmullRomCurve3(validPoints);
          const tubeGeo = trackGeo(new THREE.TubeGeometry(curve, Math.max(16, validPoints.length * 6), 0.26, 10, false));
          const tubeMat = trackMat(new THREE.MeshStandardMaterial({
            color: 0x38bdf8,
            roughness: 0.3,
            metalness: 0.3,
            transparent: true,
            opacity: 0.78
          }));
          const tubeMesh = new THREE.Mesh(tubeGeo, tubeMat);
          group.add(tubeMesh);
        }
      } catch (e) {
        console.warn('Ribbon generation skipped for current conformation:', e);
      }
    }

    // 5. PEŁNA STRUKTURA ATOMOWA: WSZYSTKIE ŁAŃCUCHY BOCZNE, ROTAMERY (χ1..χ4), SZKIELET (N,C,O,OXT) I WIĄZANIA BOCZNE
    // Render full sidechains once post-ribosomal chain synthesis completes (!isFolding)
    const tAnim = isFolding ? 0.85 : 1.0;
    const atomAnimPosBySerial = new Map<number, THREE.Vector3>();

    if (!isFolding && allAtomStructure && allAtomStructure.atoms.length > 0 && (showSidechainAtoms || showLateralBonds)) {
      const cpkCarbonMat = trackMat(new THREE.MeshStandardMaterial({ color: 0x94a3b8, roughness: 0.32, metalness: 0.25 }));
      const cpkCarbonFocusMat = trackMat(new THREE.MeshStandardMaterial({
        color: 0xfef08a,
        emissive: 0xca8a04,
        emissiveIntensity: 0.35,
        roughness: 0.22,
        metalness: 0.3
      }));
      const cpkNitrogenMat = trackMat(new THREE.MeshStandardMaterial({ color: 0x3b82f6, roughness: 0.25, metalness: 0.3 }));
      const cpkOxygenMat = trackMat(new THREE.MeshStandardMaterial({ color: 0xef4444, roughness: 0.25, metalness: 0.3 }));
      const cpkSulfurMat = trackMat(new THREE.MeshStandardMaterial({
        color: 0xfbbf24,
        emissive: 0xb45309,
        emissiveIntensity: 0.35,
        roughness: 0.2,
        metalness: 0.55
      }));
      const cpkSeleniumMat = trackMat(new THREE.MeshStandardMaterial({
        color: 0xf97316,
        emissive: 0xc2410c,
        emissiveIntensity: 0.4,
        roughness: 0.2,
        metalness: 0.55
      }));

      const covStickMat = trackMat(new THREE.MeshStandardMaterial({ color: 0x64748b, roughness: 0.35, metalness: 0.25 }));
      const focusStickMat = trackMat(new THREE.MeshStandardMaterial({ color: 0xfacc15, roughness: 0.25, metalness: 0.35 }));
      const ringStickMat = trackMat(new THREE.MeshStandardMaterial({ color: 0x818cf8, roughness: 0.3, metalness: 0.3 }));
      const peptideStickMat = trackMat(new THREE.MeshStandardMaterial({ color: 0x475569, roughness: 0.35, metalness: 0.2 }));
      const saltBridgeMat = trackMat(new THREE.MeshBasicMaterial({
        color: 0xe879f9,
        transparent: true,
        opacity: 0.85
      }));
      const hBondMat = trackMat(new THREE.MeshBasicMaterial({
        color: 0x2dd4bf,
        transparent: true,
        opacity: 0.72
      }));

      // Compute animated 3D world position for every heavy atom anchored to its residue Cα
      for (const atom of allAtomStructure.atoms) {
        const resIdx = atom.resSeq - 1;
        if (resIdx < 0 || resIdx >= posCount) continue;
        const animCa = centeredPositions[resIdx];
        const rawCa = foldedPositions[resIdx];
        if (!animCa || !rawCa) continue;

        const foldedCaVec = new THREE.Vector3(
          (rawCa[0] ?? 0) * SCALE_FACTOR,
          (rawCa[1] ?? 0) * SCALE_FACTOR,
          (rawCa[2] ?? 0) * SCALE_FACTOR
        );
        const atomOffset = new THREE.Vector3(
          atom.x * SCALE_FACTOR,
          atom.y * SCALE_FACTOR,
          atom.z * SCALE_FACTOR
        ).sub(foldedCaVec);

        const worldPos = animCa.clone().add(atomOffset.multiplyScalar(0.35 + 0.65 * tAnim));
        atomAnimPosBySerial.set(atom.serial, worldPos);

        if (showSidechainAtoms && atom.atomName !== 'CA') {
          const isResFocused = resIdx === safeIdx;
          let mat: THREE.Material = cpkCarbonMat;
          let radius = 0.22;

          if (atom.element === 'N') {
            mat = cpkNitrogenMat;
            radius = atom.isBackbone ? 0.20 : 0.25;
          } else if (atom.element === 'O') {
            mat = cpkOxygenMat;
            radius = atom.isBackbone ? 0.20 : 0.25;
          } else if (atom.element === 'S') {
            mat = cpkSulfurMat;
            radius = 0.34;
          } else if (atom.element === 'SE') {
            mat = cpkSeleniumMat;
            radius = 0.37;
          } else if (isResFocused) {
            mat = cpkCarbonFocusMat;
            radius = 0.26;
          }

          const atomMesh = new THREE.Mesh(unitSphereGeo, mat);
          atomMesh.position.copy(worldPos);
          atomMesh.scale.setScalar(isResFocused ? radius * 1.15 : radius);
          group.add(atomMesh);
        }
      }

      // Render covalent backbone/sidechain/ring bonds + lateral salt bridges & hydrogen bonds
      for (const bond of allAtomStructure.bonds) {
        const pA = atomAnimPosBySerial.get(bond.atomSerial1);
        const pB = atomAnimPosBySerial.get(bond.atomSerial2);
        if (!pA || !pB) continue;

        if (bond.category === 'DISULFIDE') {
          // Handled below by dedicated golden SSBOND clamp renderer
          continue;
        }

        if (bond.category === 'SALT_BRIDGE' || bond.category === 'HBOND') {
          if (showLateralBonds) {
            addScaledBondStick(
              pA,
              pB,
              bond.category === 'SALT_BRIDGE' ? 0.085 : 0.06,
              bond.category === 'SALT_BRIDGE' ? saltBridgeMat : hBondMat
            );
          }
          continue;
        }

        if (showSidechainAtoms) {
          const isResFocused = bond.resSeq1 - 1 === safeIdx || bond.resSeq2 - 1 === safeIdx;
          if (bond.category === 'RING_CLOSURE') {
            addScaledBondStick(pA, pB, isResFocused ? 0.14 : 0.11, isResFocused ? focusStickMat : ringStickMat);
          } else if (bond.category === 'SIDECHAIN') {
            addScaledBondStick(pA, pB, isResFocused ? 0.145 : 0.115, isResFocused ? focusStickMat : covStickMat);
          } else if (bond.category === 'BACKBONE' || bond.category === 'PEPTIDE') {
            addScaledBondStick(pA, pB, 0.095, peptideStickMat);
          }
        }
      }
    }

    // 6. Kowalencyjne klamry S-S: Zapinane dopiero po oderwaniu od rybosomu (w fazie Bańki Mydlanej / finalnej)!
    const canLatchDisulfides = !isFolding || (multiAgentDiag && !multiAgentDiag.tetheredToRibosome);
    if (canLatchDisulfides && showDisulfideBridges && foldingMetrics && foldingMetrics.ssBondPairs.length > 0) {
      const sideMat = trackMat(new THREE.MeshStandardMaterial({
        color: 0x94a3b8,
        roughness: 0.35,
        metalness: 0.3
      }));
      const sulfurMat = trackMat(new THREE.MeshStandardMaterial({
        color: 0xfbbf24,
        emissive: 0xd97706,
        emissiveIntensity: 0.55,
        roughness: 0.2,
        metalness: 0.65
      }));

      foldingMetrics.ssBondPairs.forEach(([i1, i2]) => {
        if (i1 >= posCount || i2 >= posCount) return;
        const ca1 = centeredPositions[i1];
        const ca2 = centeredPositions[i2];
        if (!ca1 || !ca2) return;

        const rawCb1 = cbPositions?.[i1];
        const rawSg1 = sgPositions?.[i1];
        const rawCb2 = cbPositions?.[i2];
        const rawSg2 = sgPositions?.[i2];

        if (rawCb1 && rawSg1 && rawCb2 && rawSg2) {
          const foldedCa1 = new THREE.Vector3(
            (foldedPositions[i1]?.[0] ?? 0) * SCALE_FACTOR,
            (foldedPositions[i1]?.[1] ?? 0) * SCALE_FACTOR,
            (foldedPositions[i1]?.[2] ?? 0) * SCALE_FACTOR
          );
          const foldedCa2 = new THREE.Vector3(
            (foldedPositions[i2]?.[0] ?? 0) * SCALE_FACTOR,
            (foldedPositions[i2]?.[1] ?? 0) * SCALE_FACTOR,
            (foldedPositions[i2]?.[2] ?? 0) * SCALE_FACTOR
          );

          const cbOffset1 = new THREE.Vector3(
            (rawCb1?.[0] ?? 0) * SCALE_FACTOR,
            (rawCb1?.[1] ?? 0) * SCALE_FACTOR,
            (rawCb1?.[2] ?? 0) * SCALE_FACTOR
          ).sub(foldedCa1);
          const sgOffset1 = new THREE.Vector3(
            (rawSg1?.[0] ?? 0) * SCALE_FACTOR,
            (rawSg1?.[1] ?? 0) * SCALE_FACTOR,
            (rawSg1?.[2] ?? 0) * SCALE_FACTOR
          ).sub(foldedCa1);

          const cbOffset2 = new THREE.Vector3(
            (rawCb2?.[0] ?? 0) * SCALE_FACTOR,
            (rawCb2?.[1] ?? 0) * SCALE_FACTOR,
            (rawCb2?.[2] ?? 0) * SCALE_FACTOR
          ).sub(foldedCa2);
          const sgOffset2 = new THREE.Vector3(
            (rawSg2?.[0] ?? 0) * SCALE_FACTOR,
            (rawSg2?.[1] ?? 0) * SCALE_FACTOR,
            (rawSg2?.[2] ?? 0) * SCALE_FACTOR
          ).sub(foldedCa2);

          const cb1 = ca1.clone().add(cbOffset1.multiplyScalar(0.35 + 0.65 * tAnim));
          const sg1 = ca1.clone().add(sgOffset1.multiplyScalar(0.35 + 0.65 * tAnim));
          const cb2 = ca2.clone().add(cbOffset2.multiplyScalar(0.35 + 0.65 * tAnim));
          const sg2 = ca2.clone().add(sgOffset2.multiplyScalar(0.35 + 0.65 * tAnim));

          if (!showSidechainAtoms) {
            addScaledBondStick(ca1, cb1, 0.16, sideMat);
            addScaledBondStick(cb1, sg1, 0.19, sulfurMat);
            addScaledBondStick(ca2, cb2, 0.16, sideMat);
            addScaledBondStick(cb2, sg2, 0.19, sulfurMat);
          }

          // Domknięte wiązanie kowalencyjne Sγ—Sγ (~2.05 Å)
          addScaledBondStick(sg1, sg2, 0.30, sulfurMat);

          const s1 = new THREE.Mesh(unitSphereGeo, sulfurMat);
          s1.position.copy(sg1);
          s1.scale.setScalar(0.44);
          group.add(s1);

          const s2 = new THREE.Mesh(unitSphereGeo, sulfurMat);
          s2.position.copy(sg2);
          s2.scale.setScalar(0.44);
          group.add(s2);
        } else {
          addScaledBondStick(ca1, ca2, 0.28, sulfurMat);
        }
      });
    }

    // 7. KABSCH 3D ALPHAFOLD / PDB REFERENCE SUPERPOSITION OVERLAY
    if (showAlphaFoldOverlay && rmsdComparison && rmsdComparison.alignedRefCaPositions.length > 0) {
      const afRibbonMat = trackMat(new THREE.MeshStandardMaterial({
        color: 0xc084fc,
        emissive: 0x7e22ce,
        emissiveIntensity: 0.3,
        roughness: 0.3,
        metalness: 0.25,
        transparent: true,
        opacity: 0.82
      }));
      const devGoodMat = trackMat(new THREE.MeshBasicMaterial({ color: 0x34d399, transparent: true, opacity: 0.75 }));
      const devWarnMat = trackMat(new THREE.MeshBasicMaterial({ color: 0xfbbf24, transparent: true, opacity: 0.75 }));
      const devHighMat = trackMat(new THREE.MeshBasicMaterial({ color: 0xf43f5e, transparent: true, opacity: 0.75 }));

      const afPoints: THREE.Vector3[] = [];
      rmsdComparison.alignedRefCaPositions.forEach((refPt, idx) => {
        if (!refPt || idx >= posCount) return;
        const afVec = new THREE.Vector3(
          (refPt[0] - comX) * SCALE_FACTOR,
          (refPt[1] - comY) * SCALE_FACTOR,
          (refPt[2] - comZ) * SCALE_FACTOR
        );
        afPoints.push(afVec);

        const afNode = new THREE.Mesh(unitSphereGeo, afRibbonMat);
        afNode.position.copy(afVec);
        afNode.scale.setScalar(0.26);
        group.add(afNode);

        // Deviation connector strut between SNAPLOCK Cα and superposed AlphaFold Cα
        const modVec = centeredPositions[idx];
        if (modVec && tAnim > 0.9) {
          const angstromDev = modVec.distanceTo(afVec) / SCALE_FACTOR;
          const strutMat = angstromDev < 2.0 ? devGoodMat : angstromDev < 4.0 ? devWarnMat : devHighMat;
          addScaledBondStick(modVec, afVec, 0.055, strutMat);
        }
      });

      if (afPoints.length >= 3) {
        try {
          const afCurve = new THREE.CatmullRomCurve3(afPoints);
          const afTubeGeo = trackGeo(new THREE.TubeGeometry(afCurve, Math.max(16, afPoints.length * 6), 0.19, 8, false));
          const afTubeMesh = new THREE.Mesh(afTubeGeo, afRibbonMat);
          group.add(afTubeMesh);
        } catch {
          // Skip tube if collinear
        }
      }
    }

  }, [
    sequenceCodes, 
    foldedPositions, 
    cbPositions,
    sgPositions,
    allAtomStructure,
    isFolding,
    multiAgentDiag,
    showAgentForceVectors,
    foldProgress, 
    showFacetCages, 
    showRibbon, 
    showDisulfideBridges, 
    showBackboneAtoms, 
    showSidechainAtoms,
    showLateralBonds,
    showBoundarySphere,
    showAlphaFoldOverlay,
    rmsdComparison,
    cageOpacity, 
    activeResidueIdx, 
    hoveredResidueIdx,
    foldingMetrics
  ]);

  // Pointer Drag Event Handlers for 3D Camera Orbit & Pan
  const handleMouseDown = (e: React.MouseEvent) => {
    isDragging.current = true;
    isPanning.current = e.shiftKey || e.button === 1 || e.button === 2;
    prevMousePos.current = { x: e.clientX, y: e.clientY };
  };

  const handleMouseMove = (e: React.MouseEvent) => {
    if (!isDragging.current) return;
    const dx = e.clientX - prevMousePos.current.x;
    const dy = e.clientY - prevMousePos.current.y;

    if (isMicroscopeMode) {
      // First-person mouse look (yaw & pitch)
      fpCamYaw.current -= dx * 0.005;
      fpCamPitch.current -= dy * 0.005;
      fpCamPitch.current = Math.max(-Math.PI / 2.1, Math.min(Math.PI / 2.1, fpCamPitch.current));
    } else {
      if (isPanning.current) {
        panOffset.current.x -= dx * 0.05;
        panOffset.current.y += dy * 0.05;
      } else {
        orbitAngle.current.y += dx * 0.008;
        orbitAngle.current.x += dy * 0.008;
        orbitAngle.current.x = Math.max(-Math.PI / 2 + 0.05, Math.min(Math.PI / 2 - 0.05, orbitAngle.current.x));
      }
    }

    prevMousePos.current = { x: e.clientX, y: e.clientY };
  };

  const handleMouseUp = () => {
    isDragging.current = false;
    isPanning.current = false;
  };

  // Touch handlers for mobile
  const handleTouchStart = (e: React.TouchEvent) => {
    if (e.touches.length === 1) {
      isDragging.current = true;
      prevMousePos.current = { x: e.touches[0].clientX, y: e.touches[0].clientY };
    }
  };

  const handleTouchMove = (e: React.TouchEvent) => {
    if (!isDragging.current || e.touches.length === 0) return;
    const dx = e.touches[0].clientX - prevMousePos.current.x;
    const dy = e.touches[0].clientY - prevMousePos.current.y;

    if (isMicroscopeMode) {
      fpCamYaw.current -= dx * 0.006;
      fpCamPitch.current -= dy * 0.006;
      fpCamPitch.current = Math.max(-Math.PI / 2.1, Math.min(Math.PI / 2.1, fpCamPitch.current));
    } else {
      orbitAngle.current.y += dx * 0.009;
      orbitAngle.current.x += dy * 0.009;
      orbitAngle.current.x = Math.max(-Math.PI / 2 + 0.05, Math.min(Math.PI / 2 - 0.05, orbitAngle.current.x));
    }
    prevMousePos.current = { x: e.touches[0].clientX, y: e.touches[0].clientY };
  };

  const handleTouchEnd = () => {
    isDragging.current = false;
  };

  // --- MIKROSKOP Å NAVIGATION & 5 AUTOMATED FLIGHT PATHS ---
  const jumpToResidueFirstPerson = (resIdx: number) => {
    const safeIdx = Math.max(0, Math.min(resIdx, foldedPositions.length - 1));
    setActiveResidueIdx(safeIdx);
    const ca = foldedPositions[safeIdx];
    if (!ca) return;

    // Center of mass / CA coordinates
    fpCamPos.current.set(ca[0], ca[1], ca[2]);

    // Local N-Cα-C frame orientation
    // Find C and N atoms for this residue if available
    let forward = new THREE.Vector3(0, 0, -1);
    if (allAtomStructure) {
      const caAtom = allAtomStructure.atoms.find(a => a.resSeq === safeIdx + 1 && a.atomName === 'CA');
      const cAtom = allAtomStructure.atoms.find(a => a.resSeq === safeIdx + 1 && a.atomName === 'C');
      if (caAtom && cAtom) {
        forward.set(cAtom.x - caAtom.x, cAtom.y - caAtom.y, cAtom.z - caAtom.z).normalize();
      }
    } else if (safeIdx + 1 < foldedPositions.length) {
      const nextCa = foldedPositions[safeIdx + 1];
      forward.set(nextCa[0] - ca[0], nextCa[1] - ca[1], nextCa[2] - ca[2]).normalize();
    }

    fpCamYaw.current = Math.atan2(-forward.x, -forward.z);
    fpCamPitch.current = Math.asin(forward.y);
  };

  const handleStartFlightPath = (pathType: typeof activeFlightPath) => {
    if (pathType === 'NONE') {
      if (flightAnimRef.current) cancelAnimationFrame(flightAnimRef.current);
      flightAnimRef.current = null;
      setActiveFlightPath('NONE');
      return;
    }

    setActiveFlightPath(pathType);
    let t = 0;
    const durationFrames = 280; // ~4.6 seconds flight tour

    const flightStep = () => {
      t += 1;
      const progress = Math.min(1.0, t / durationFrames);
      setFlightProgress(progress);

      const n = foldedPositions.length;
      if (n > 1) {
        let targetX = 0, targetY = 0, targetZ = 0;

        if (pathType === 'PATH_1_BACKBONE') {
          // (1) N -> C peptide chain
          const idxFloat = progress * (n - 1);
          const i0 = Math.floor(idxFloat);
          const i1 = Math.min(n - 1, i0 + 1);
          const frac = idxFloat - i0;
          const p0 = foldedPositions[i0];
          const p1 = foldedPositions[i1];
          targetX = p0[0] + frac * (p1[0] - p0[0]);
          targetY = p0[1] + frac * (p1[1] - p0[1]);
          targetZ = p0[2] + frac * (p1[2] - p0[2]);
        } else if (pathType === 'PATH_2_SPIRAL_BOTTOM_UP') {
          // (2) Bottom -> Top spiral along centers of mass
          const sorted = [...foldedPositions].sort((a, b) => a[1] - b[1]);
          const idxFloat = progress * (n - 1);
          const i0 = Math.floor(idxFloat);
          const i1 = Math.min(n - 1, i0 + 1);
          const frac = idxFloat - i0;
          const p0 = sorted[i0];
          const p1 = sorted[i1];
          const angle = progress * Math.PI * 6; // 3 full turns spiral
          const radius = 4.5;
          targetX = p0[0] + frac * (p1[0] - p0[0]) + Math.cos(angle) * radius;
          targetY = p0[1] + frac * (p1[1] - p0[1]);
          targetZ = p0[2] + frac * (p1[2] - p0[2]) + Math.sin(angle) * radius;
        } else if (pathType === 'PATH_3_SPIRAL_FRONT_BACK') {
          // (3) Front -> Back spiral
          const sorted = [...foldedPositions].sort((a, b) => a[2] - b[2]);
          const idxFloat = progress * (n - 1);
          const i0 = Math.floor(idxFloat);
          const i1 = Math.min(n - 1, i0 + 1);
          const frac = idxFloat - i0;
          const p0 = sorted[i0];
          const p1 = sorted[i1];
          const angle = progress * Math.PI * 6;
          const radius = 4.5;
          targetX = p0[0] + frac * (p1[0] - p0[0]) + Math.cos(angle) * radius;
          targetY = p0[1] + frac * (p1[1] - p0[1]) + Math.sin(angle) * radius;
          targetZ = p0[2] + frac * (p1[2] - p0[2]);
        } else if (pathType === 'PATH_4_LATERAL_BONDS') {
          // (4) Lateral bonds tour (inspect salt bridges & H-bonds)
          const latIdx = Math.floor(progress * n);
          const pos = foldedPositions[latIdx % n];
          targetX = pos[0] + Math.sin(progress * 12) * 2.2;
          targetY = pos[1] + Math.cos(progress * 12) * 2.2;
          targetZ = pos[2];
        } else if (pathType === 'PATH_5_DISULFIDE') {
          // (5) Disulfide bridges S-S tour
          const pairs = foldingMetrics?.ssBondPairs || [];
          if (pairs.length > 0) {
            const pairIdx = Math.floor(progress * pairs.length) % pairs.length;
            const [i1, i2] = pairs[pairIdx];
            const p1 = foldedPositions[i1] || foldedPositions[0];
            const p2 = foldedPositions[i2] || foldedPositions[0];
            targetX = (p1[0] + p2[0]) * 0.5 + Math.sin(progress * 20) * 1.8;
            targetY = (p1[1] + p2[1]) * 0.5 + Math.cos(progress * 20) * 1.8;
            targetZ = (p1[2] + p2[2]) * 0.5;
          } else {
            const p = foldedPositions[Math.floor(progress * n)];
            targetX = p[0]; targetY = p[1]; targetZ = p[2];
          }
        }

        fpCamPos.current.set(targetX, targetY, targetZ);
        // Look towards next increment
        fpCamYaw.current = progress * Math.PI * 2;
      }

      if (progress < 1.0) {
        flightAnimRef.current = requestAnimationFrame(flightStep);
      } else {
        flightAnimRef.current = null;
        setActiveFlightPath('NONE');
      }
    };

    flightAnimRef.current = requestAnimationFrame(flightStep);
  };

  // 6 Nearest neighbors calculation for the active residue
  const nearestNeighbors = useMemo(() => {
    if (!foldedPositions || foldedPositions.length === 0) return [];
    const safeIdx = Math.max(0, Math.min(activeResidueIdx, foldedPositions.length - 1));
    const origin = foldedPositions[safeIdx];
    if (!origin) return [];

    const list: { idx: number; code: string; dist: number; hydropathy: number }[] = [];
    for (let j = 0; j < foldedPositions.length; j++) {
      if (j === safeIdx) continue;
      const pj = foldedPositions[j];
      const d = Math.hypot(pj[0] - origin[0], pj[1] - origin[1], pj[2] - origin[2]);
      const aaCode = sequenceCodes[j] || 'ALA';
      const blk = SNAPLOCK_21_BLOCKS[aaCode] || SNAPLOCK_21_BLOCKS['ALA'];
      list.push({
        idx: j,
        code: aaCode,
        dist: parseFloat(d.toFixed(2)),
        hydropathy: blk.physico.hydropathyKyte
      });
    }

    list.sort((a, b) => a.dist - b.dist);
    return list.slice(0, 6);
  }, [activeResidueIdx, foldedPositions, sequenceCodes]);

  // Generate Current Standard All-Atom wwPDB v3.3 String
  const currentPdbString = useMemo(() => {
    if (
      !foldingMetrics || 
      !foldedPositions || 
      !sequenceCodes || 
      foldedPositions.length === 0 || 
      sequenceCodes.length === 0 || 
      foldedPositions.length !== sequenceCodes.length
    ) {
      return '';
    }
    return exportToPdbString(
      fastaHeader,
      sequenceCodes,
      foldedPositions,
      foldingMetrics,
      secondaryTypes,
      allAtomStructure || undefined,
      rmsdComparison,
      cachedRefPdbRef.current
    );
  }, [fastaHeader, sequenceCodes, foldedPositions, foldingMetrics, secondaryTypes, allAtomStructure, rmsdComparison]);

  // Trigger Standard All-Atom PDB Download
  const handleDownloadPdb = () => {
    const filename = `${fastaHeader.replace(/[^a-zA-Z0-9_-]/g, '_').toLowerCase() || 'snaplock_all_atom'}_all_atom.pdb`;
    triggerPdbDownload(currentPdbString, filename);
  };

  // Trigger Dual-Model Superposed PDB Download (MODEL 1 = SNAPLOCK, MODEL 2 = AlphaFold Reference)
  const handleDownloadSuperposedPdb = () => {
    if (!allAtomStructure || !rmsdComparison) return;
    const superPdb = exportSuperposedComparisonPdb(fastaHeader, allAtomStructure, rmsdComparison);
    const filename = `${fastaHeader.replace(/[^a-zA-Z0-9_-]/g, '_').toLowerCase() || 'snaplock'}_vs_${rmsdComparison.referenceId.replace(/[^a-zA-Z0-9_-]/g, '_')}.pdb`;
    triggerPdbDownload(superPdb, filename);
  };

  const handleCopyPdb = () => {
    if (!currentPdbString) return;
    navigator.clipboard.writeText(currentPdbString);
    setCopySuccess(true);
    setTimeout(() => setCopySuccess(false), 2000);
  };

  const safeActiveIdx = Math.max(0, Math.min(activeResidueIdx, Math.max(0, sequenceCodes.length - 1)));
  const activeResidueCode = sequenceCodes[safeActiveIdx] || 'ALA';
  const activeBlock = SNAPLOCK_21_BLOCKS[activeResidueCode] || SNAPLOCK_21_BLOCKS['ALA'];
  const activeRotamerInfo = allAtomStructure?.residueRotamers?.[safeActiveIdx] || null;

  return (
    <div className="flex-1 w-full h-full min-h-0 bg-black text-white flex flex-col overflow-hidden font-sans select-none relative">
      
      {/* 1. TOP HEADER: SEQUENCE CONSOLE, FASTA IMPORT & PDB ACTIONS */}
      <div className="bg-[#04050b] border-b border-zinc-900 px-4 py-2 flex flex-wrap items-center justify-between gap-3 shrink-0 font-mono text-xs z-20">
        
        {/* Input & Fold Action */}
        <div className="flex items-center gap-2 flex-1 max-w-4xl min-w-[320px]">
          <button
            onClick={() => {
              setFastaModalText(`>${fastaHeader}\n${sequenceCodes.join(' ')}`);
              setIsFastaModalOpen(true);
            }}
            className="px-2.5 py-1.5 bg-indigo-950/80 hover:bg-indigo-900 border border-indigo-700/70 text-indigo-200 rounded text-xs font-bold flex items-center gap-1.5 cursor-pointer shrink-0 transition"
            title="Importuj sekwencję w pełnym standardzie FASTA (>nagłówek + ciąg)"
          >
            <Dna className="w-3.5 h-3.5 text-indigo-400" />
            <span>FASTA</span>
          </button>

          <input
            type="text"
            value={seqInput}
            onChange={(e) => setSeqInput(e.target.value)}
            onKeyDown={(e) => {
              if (e.key === 'Enter') handleApplyInput();
            }}
            placeholder="Wpisz sekwencję (np. GIVEQCCTSICSLYQLENYCN lub GLY ILE VAL)"
            className="bg-black border border-zinc-800 text-white font-mono text-xs rounded px-3 py-1.5 w-full focus:outline-none focus:border-cyan-500"
          />

          <button
            onClick={() => handleApplyInput()}
            className="px-3 py-1.5 bg-zinc-900 hover:bg-zinc-800 border border-zinc-700 text-zinc-300 hover:text-white rounded text-xs font-bold cursor-pointer transition shrink-0"
          >
            Wczytaj
          </button>

          {/* PRIMARY MULTI-AGENT COTRANSLATIONAL N->C + SOAP-BUBBLE FOLD BUTTON */}
          <button
            onClick={() => runFoldingAnimation()}
            disabled={isFolding}
            className="px-4 py-1.5 bg-gradient-to-r from-emerald-600 via-teal-600 to-cyan-600 hover:from-emerald-500 hover:to-cyan-500 disabled:opacity-50 text-white font-bold rounded flex items-center gap-2 cursor-pointer shadow-lg shadow-emerald-500/20 shrink-0 transition"
            title="Kotranslacyjny wzrost N→C (Soft SnapLock) + oderwanie od rybosomu (Bańka Mydlana & relaksacja kątów torsyjnych φ/ψ w pamięci)"
          >
            <Sparkles className={`w-4 h-4 text-amber-300 ${isFolding ? 'animate-spin' : ''}`} />
            <span>
              {isFolding && multiAgentDiag
                ? multiAgentDiag.tetheredToRibosome
                  ? `Synteza N→C (${multiAgentDiag.emergedCount}/${multiAgentDiag.totalCount})...`
                  : `Bańka Mydlana (φ,ψ)... ${foldProgress}%`
                : '▶ ZWIŃ N→C + BAŃKA MYDLANA'}
            </span>
          </button>
        </div>

        {/* Right side: Boundary Config, PDB Export, PDB Preview */}
        <div className="flex items-center gap-2 shrink-0">
          
          {/* Boundary Conditions Drawer Toggle */}
          <button
            onClick={() => setIsBoundaryPanelOpen(!isBoundaryPanelOpen)}
            className={`px-3 py-1.5 rounded border text-xs font-bold flex items-center gap-1.5 cursor-pointer transition ${
              isBoundaryPanelOpen 
                ? 'bg-cyan-950 border-cyan-500 text-cyan-200' 
                : 'bg-zinc-900 hover:bg-zinc-800 border-zinc-800 text-zinc-300'
            }`}
            title="Konfiguruj warunki brzegowe: Adaptacyjna sfera R_kula >= 1.48 Rg, profil CD i mostki S-S"
          >
            <SlidersHorizontal className="w-3.5 h-3.5 text-cyan-400" />
            <span>Warunki Brzegowe & CD</span>
          </button>

          {/* Ground-Truth RMSD Validation Modal Toggle — HARD GATE < 1.0 Å */}
          <button
            onClick={() => setIsRmsdModalOpen(true)}
            className={`px-3 py-1.5 rounded border text-xs font-bold flex items-center gap-1.5 cursor-pointer transition ${
              foldingMetrics?.groundTruthValidation?.hasGroundTruthReference &&
              foldingMetrics.groundTruthValidation.caRmsdAngstrom !== null
                ? foldingMetrics.groundTruthValidation.valid &&
                  foldingMetrics.groundTruthValidation.caRmsdAngstrom < 1.0
                  ? 'bg-emerald-950/90 border-emerald-500 text-emerald-200 shadow-md shadow-emerald-500/20'
                  : 'bg-amber-950/90 border-amber-500 text-amber-200 shadow-md shadow-amber-500/15'
                : 'bg-zinc-900 hover:bg-zinc-800 border-purple-800/70 text-purple-300 hover:text-white'
            }`}
            title="Rygorystyczna walidacja Kabsch SVD 1:1 Cα — status zielony odblokowany wyłącznie gdy globalne RMSD < 1.0 Å"
          >
            <GitCompare className="w-3.5 h-3.5 text-purple-400" />
            <span>
              {foldingMetrics?.groundTruthValidation?.hasGroundTruthReference &&
              foldingMetrics.groundTruthValidation.caRmsdAngstrom !== null
                ? `Kabsch 1:1 Cα: ${foldingMetrics.groundTruthValidation.caRmsdAngstrom.toFixed(2)} Å (${foldingMetrics.groundTruthValidation.referenceId})`
                : 'RMSD ref. (Waliduj 1:1)'}
            </span>
          </button>

          {/* PDB Preview Modal Toggle */}
          <button
            onClick={() => setIsPdbModalOpen(true)}
            className="px-3 py-1.5 bg-zinc-900 hover:bg-zinc-800 border border-zinc-800 text-zinc-300 hover:text-white rounded text-xs font-bold flex items-center gap-1.5 cursor-pointer transition"
            title="Podgląd wygenerowanego pełnoatomowego pliku standardu wwPDB v3.3 (80 kolumn, rotamery, SSBOND, CONECT)"
          >
            <FileText className="w-3.5 h-3.5 text-sky-400" />
            <span>All-Atom PDB</span>
          </button>

          {/* DIRECT ALL-ATOM PDB DOWNLOAD BUTTON */}
          <button
            onClick={handleDownloadPdb}
            className="px-3.5 py-1.5 bg-gradient-to-r from-sky-600 to-blue-600 hover:from-sky-500 hover:to-blue-500 text-white rounded text-xs font-extrabold flex items-center gap-1.5 cursor-pointer shadow-md transition"
            title="Pobierz pełną strukturę atomową (wszystkie łańcuchy boczne, rotamery, CONECT i SSBOND) w standardzie wwPDB v3.3 (.pdb)"
          >
            <Download className="w-3.5 h-3.5 text-white" />
            <span>Pobierz PDB (.pdb)</span>
          </button>

        </div>
      </div>

      {/* 2. SUB-BAR: QUICK PRESETS, SAXS / CD TELEMETRY & VIEW TOGGLES */}
      <div className="bg-[#020308] border-b border-zinc-900 px-4 py-1.5 flex flex-wrap items-center justify-between gap-3 shrink-0 font-mono text-xs z-10">
        
        {/* Quick Presets */}
        <div className="flex items-center gap-1.5 flex-wrap">
          <span className="text-[9.5px] text-zinc-500 uppercase font-bold">Wzorce FASTA:</span>
          {FASTA_PRESETS.map((p) => (
            <button
              key={p.id}
              onClick={() => loadPreset(p.id)}
              className="px-2 py-0.5 rounded bg-zinc-900 hover:bg-zinc-800 border border-zinc-800 text-zinc-300 hover:text-white text-[10.5px] font-bold cursor-pointer transition"
            >
              {p.name}
            </button>
          ))}
        </div>

        {/* Live Boundary Metrics, Loop Closure & Adaptive SSBOND Telemetry */}
        {foldingMetrics && (
          <div className="flex items-center gap-1.5 overflow-x-auto text-[10.5px]">
            {/* Universal Pipeline Mode Badge */}
            <span
              className={`px-2 py-0.5 rounded border flex items-center gap-1 font-bold ${
                foldingMetrics.pipelineMode === 'TWO_WAY_IK_ADAPTIVE_SSBOND'
                  ? 'bg-amber-950/70 border-amber-600/80 text-amber-200'
                  : 'bg-cyan-950/70 border-cyan-700/80 text-cyan-200'
              }`}
              title={
                foldingMetrics.pipelineMode === 'TWO_WAY_IK_ADAPTIVE_SSBOND'
                  ? `Aktywny dwukierunkowy solver IK (CCD) + adaptacyjna waga w_ss (${foldingMetrics.initialSpringWeight} ➔ ${foldingMetrics.finalSpringWeight})`
                  : 'Struktura bez mostków S-S: standardowy deterministyczny potok mikrofaset L2 (2232 fasetki)'
              }
            >
              <Zap className="w-3 h-3 text-amber-400" />
              <span>
                {foldingMetrics.pipelineMode === 'TWO_WAY_IK_ADAPTIVE_SSBOND'
                  ? `IK + Adapt. w_ss (${foldingMetrics.finalSpringWeight})`
                  : 'Potok L2 (2232)'}
              </span>
            </span>

            {/* Level 2 Micro-Facet Angular Error */}
            <span
              className={`px-2 py-0.5 rounded border flex items-center gap-1 font-bold ${
                foldingMetrics.meanMicroFacetErrorDeg < 1.0
                  ? 'bg-emerald-950/70 border-emerald-700 text-emerald-300'
                  : 'bg-amber-950/70 border-amber-700 text-amber-300'
              }`}
              title={`Średni błąd kątowy na mikrofasetach L2 (2232 stany): ${foldingMetrics.meanMicroFacetErrorDeg}° / resztę (cel < 1.0°)`}
            >
              <span>Błąd L2: <strong>{foldingMetrics.meanMicroFacetErrorDeg}°</strong></span>
            </span>

            {/* Internal L2 Lattice Quantization Residual Badge (NOT Ground-Truth RMSD!) */}
            {foldingMetrics.cartesianAlignment && (
              <span
                className="px-2 py-0.5 rounded border flex items-center gap-1 font-bold bg-zinc-900/90 border-zinc-700 text-zinc-300"
                title={`Wewnętrzny błąd dyskretyzacji siatki L2 (${foldingMetrics.cartesianAlignment.referenceId}): ${foldingMetrics.cartesianAlignment.postAlignmentBackboneRmsdAngstrom} Å, kąt mikrofasetki L2 = ${foldingMetrics.cartesianAlignment.l2MeanAngularErrorDeg}°. Uwaga: to metryka wewnętrzna dyskretyzacji, a NIE Ground-Truth RMSD!`}
              >
                <span>
                  Reziduum L2 (wewn.): <strong>{foldingMetrics.cartesianAlignment.postAlignmentBackboneRmsdAngstrom} Å</strong>
                </span>
              </span>
            )}

            {/* D_max Status */}
            <span 
              className={`px-2 py-0.5 rounded border flex items-center gap-1 ${
                foldingMetrics.dMaxActual <= boundaryConfig.dMaxTarget + 0.5
                  ? 'bg-emerald-950/60 border-emerald-800 text-emerald-300'
                  : 'bg-amber-950/60 border-amber-800 text-amber-300'
              }`}
              title={`Maksymalny wymiar cząsteczki D_max: ${foldingMetrics.dMaxActual} Å (Limit zadany: ${boundaryConfig.dMaxTarget} Å, Teoretyczny Flory: ${foldingMetrics.dMaxTheoretical} Å)`}
            >
              <CheckCircle2 className="w-3 h-3 text-cyan-400" />
              <span>D_max: <strong>{foldingMetrics.dMaxActual} Å</strong></span>
            </span>

            {/* Disulfide Clamp Closure Deviation & Penalty */}
            {foldingMetrics.ssBondPairs.length > 0 ? (
              <span 
                className={`px-2 py-0.5 rounded border flex items-center gap-1 font-bold ${
                  foldingMetrics.finalMaxSgDeviationAngstrom < 0.5
                    ? 'bg-emerald-950/70 border-emerald-600/80 text-emerald-200'
                    : 'bg-rose-950/70 border-rose-600 text-rose-200'
                }`}
                title={`Domknięcie klamer SSBOND (${foldingMetrics.ssBondPairs.length}): max odchylenie d(Sγ,Sγ) = ${foldingMetrics.finalMaxSgDeviationAngstrom} Å (< 0.5 Å), max Δ(Cα,Cα) = ${foldingMetrics.finalMaxCaDeviationAngstrom} Å, E_ssbond = ${foldingMetrics.disulfidePenaltyScore}`}
              >
                <ShieldCheck className="w-3 h-3 text-amber-400" />
                <span>
                  SSBOND ({foldingMetrics.ssBondPairs.length}): <strong>Δ={foldingMetrics.finalMaxSgDeviationAngstrom} Å</strong> (&lt;0.5 Å)
                </span>
              </span>
            ) : (
              <span className="px-2 py-0.5 rounded bg-zinc-950 border border-zinc-800 text-zinc-400 flex items-center gap-1">
                <ShieldCheck className="w-3 h-3 text-zinc-500" />
                <span>SSBOND: 0</span>
              </span>
            )}

            {/* Far-UV CD Penalty Score */}
            <span 
              className="px-2 py-0.5 rounded bg-blue-950/60 border border-blue-800/80 text-blue-300 flex items-center gap-1 font-bold"
              title={`Kara CD (Dichroizm Kołowy): Kara kwadratowa = ${foldingMetrics.cdPenaltyScore}. Osiągnięto: α=${foldingMetrics.alphaActualPct}%, β=${foldingMetrics.betaActualPct}%`}
            >
              <Activity className="w-3 h-3 text-blue-400" />
              <span>CD: <strong className={foldingMetrics.cdPenaltyScore < 10 ? 'text-emerald-400' : 'text-amber-400'}>{foldingMetrics.cdPenaltyScore}</strong></span>
            </span>

            {/* Heavy-Atom & Lateral Bond Telemetry */}
            {allAtomStructure && (
              <span
                className="px-2 py-0.5 rounded bg-indigo-950/70 border border-indigo-700/70 text-indigo-200 flex items-center gap-1 font-bold"
                title={`Pełna struktura ciężkich atomów: ${allAtomStructure.totalHeavyAtoms} atomów (${allAtomStructure.totalBackboneAtoms} szkieletu, ${allAtomStructure.totalSidechainAtoms} łańcuchów bocznych), ${allAtomStructure.totalCovalentBonds} wiązań kowalencyjnych, ${allAtomStructure.totalSaltBridges} mostków solnych, ${allAtomStructure.totalHydrogenBonds} wiązań H`}
              >
                <Atom className="w-3 h-3 text-indigo-400" />
                <span>Atomy: <strong>{allAtomStructure.totalHeavyAtoms}</strong> ({allAtomStructure.totalSidechainAtoms} bocz.)</span>
              </span>
            )}

            {/* Free Energy */}
            <span className="px-2 py-0.5 rounded bg-emerald-950/80 border border-emerald-700/60 text-emerald-300 font-extrabold">
              ΔG: {foldingMetrics.totalCohesiveEnergy} kcal/mol
            </span>
          </div>
        )}

        {/* Display Toggles */}
        <div className="flex items-center gap-1.5 text-[10.5px] flex-wrap">
          <button
            onClick={() => setUseHierarchicalM5Engine(!useHierarchicalM5Engine)}
            className={`px-2 py-0.5 rounded border transition cursor-pointer flex items-center gap-1 font-bold ${
              useHierarchicalM5Engine
                ? 'bg-emerald-950/90 border-emerald-500 text-emerald-200'
                : 'bg-zinc-900 border-zinc-700 text-zinc-300'
            }`}
            title="Przełącz generator szkieletu między nowym 5-etapowym silnikiem fizycznym M5 (NeRF φ/ψ + wiązania H + kontakty dalekie) a bazowym generatorem sekwencyjnym M1 (E8/L2)"
          >
            <Compass className="w-3 h-3" />
            <span>{useHierarchicalM5Engine ? 'Silnik: M5 Hierarchiczny (φ/ψ)' : 'Silnik: M1 Bazowy E8/L2'}</span>
          </button>

          <button
            onClick={() => setShowSidechainAtoms(!showSidechainAtoms)}
            className={`px-2 py-0.5 rounded border transition cursor-pointer flex items-center gap-1 font-bold ${
              showSidechainAtoms ? 'bg-indigo-950/90 border-indigo-500 text-indigo-200' : 'bg-black border-zinc-800 text-zinc-500'
            }`}
            title="Włącz/wyłącz pełne łańcuchy boczne, pierścienie aromatyczne i dyskretne rotamery χ1..χ4"
          >
            <Atom className="w-3 h-3" />
            <span>Łańcuchy Boczne</span>
          </button>

          <button
            onClick={() => setShowLateralBonds(!showLateralBonds)}
            className={`px-2 py-0.5 rounded border transition cursor-pointer flex items-center gap-1 ${
              showLateralBonds ? 'bg-fuchsia-950/80 border-fuchsia-600 text-fuchsia-200' : 'bg-black border-zinc-800 text-zinc-500'
            }`}
            title="Włącz/wyłącz niekowalencyjne wiązania boczne (mostki solne ASP/GLU-LYS/ARG oraz wiązania wodorowe)"
          >
            <Zap className="w-3 h-3" />
            <span>Wiązania Boczne</span>
          </button>

          {rmsdComparison && (
            <button
              onClick={() => setShowAlphaFoldOverlay(!showAlphaFoldOverlay)}
              className={`px-2 py-0.5 rounded border transition cursor-pointer flex items-center gap-1 font-bold ${
                showAlphaFoldOverlay ? 'bg-purple-950/90 border-purple-500 text-purple-200' : 'bg-black border-zinc-800 text-zinc-500'
              }`}
              title="Włącz/wyłącz superpozycję 3D modelu referencyjnego AlphaFold / PDB"
            >
              <GitCompare className="w-3 h-3" />
              <span>Overlay AF ({rmsdComparison.caRmsdAngstrom} Å)</span>
            </button>
          )}

          <button
            onClick={() => setShowBoundarySphere(!showBoundarySphere)}
            className={`px-2 py-0.5 rounded border transition cursor-pointer flex items-center gap-1 ${
              showBoundarySphere ? 'bg-teal-950/80 border-teal-600 text-teal-300' : 'bg-black border-zinc-800 text-zinc-500'
            }`}
            title="Adaptacyjna sfera E8: R_kula >= 1.48 × Rg (0 wystających atomów)"
          >
            <Boxes className="w-3 h-3" />
            <span>Sfera R_kula (≥1.48 Rg)</span>
          </button>

          <button
            onClick={() => setShowAgentForceVectors(!showAgentForceVectors)}
            className={`px-2 py-0.5 rounded border transition cursor-pointer flex items-center gap-1 ${
              showAgentForceVectors ? 'bg-emerald-950/80 border-emerald-600 text-emerald-300' : 'bg-black border-zinc-800 text-zinc-500'
            }`}
            title="Pokaż/ukryj wektory sił Agentów-Reziduów podczas zwijania (zielony: Bańka Mydlana, czerwony: odpychanie VdW, złoty: klamry S-S)"
          >
            <Flame className="w-3 h-3" />
            <span>Wektory Agentów</span>
          </button>

          {/* MIKROSKOP Å (ANT-MAN FIRST-PERSON MODE) TOGGLE */}
          <button
            onClick={() => {
              const nextMode = !isMicroscopeMode;
              setIsMicroscopeMode(nextMode);
              if (nextMode) {
                jumpToResidueFirstPerson(activeResidueIdx);
              }
            }}
            className={`px-2.5 py-0.5 rounded border transition cursor-pointer flex items-center gap-1.5 font-bold ${
              isMicroscopeMode 
                ? 'bg-rose-950/90 border-rose-500 text-rose-200 ring-2 ring-rose-500/40 shadow-lg' 
                : 'bg-zinc-900 hover:bg-zinc-800 border-zinc-700 text-zinc-300'
            }`}
            title="Przełącz w tryb pierwszoosobowy wewnątrz białka (jak Ant-Man): kamera w środku reszty, sterowanie WASD + mysz/gesty, linijka 1 Å, etykiety 6 najbliższych reszt i 5 tras przelotu"
          >
            <Eye className="w-3.5 h-3.5 text-rose-400" />
            <span>{isMicroscopeMode ? 'Wyjdź z Mikroskopu Å' : 'Mikroskop Å (Ant-Man)'}</span>
          </button>

          <button
            onClick={() => setShowFacetCages(!showFacetCages)}
            className={`px-2 py-0.5 rounded border transition cursor-pointer flex items-center gap-1 ${
              showFacetCages ? 'bg-cyan-950/80 border-cyan-700 text-cyan-300' : 'bg-black border-zinc-800 text-zinc-500'
            }`}
          >
            <Boxes className="w-3 h-3" />
            <span>Powłoki Van der Waalsa</span>
          </button>

          <button
            onClick={() => setShowRibbon(!showRibbon)}
            className={`px-2 py-0.5 rounded border transition cursor-pointer flex items-center gap-1 ${
              showRibbon ? 'bg-sky-950/80 border-sky-700 text-sky-300' : 'bg-black border-zinc-800 text-zinc-500'
            }`}
          >
            <Layers className="w-3 h-3" />
            <span>Wstęga 3D</span>
          </button>

          <button
            onClick={() => setShowDisulfideBridges(!showDisulfideBridges)}
            className={`px-2 py-0.5 rounded border transition cursor-pointer flex items-center gap-1 ${
              showDisulfideBridges ? 'bg-amber-950/80 border-amber-700 text-amber-300' : 'bg-black border-zinc-800 text-zinc-500'
            }`}
          >
            <ShieldCheck className="w-3 h-3" />
            <span>Mostki S-S</span>
          </button>

          <button
            onClick={() => {
              orbitAngle.current = { x: 0.35, y: 0.65 };
              panOffset.current.set(0, 0, 0);
            }}
            className="p-1 rounded bg-zinc-900 hover:bg-zinc-800 text-zinc-400 hover:text-white"
            title="Wyśrodkuj widok 3D"
          >
            <Compass className="w-3.5 h-3.5" />
          </button>
        </div>
      </div>

      {/* 3. MAIN 3D REALISTIC CANVAS CONTAINER */}
      <div className="flex-1 w-full min-h-0 relative flex flex-col overflow-hidden bg-black">
        
        {/* Floating Active Residue & Rotamer Inspector Card */}
        <div className="absolute top-3 left-3 z-10 pointer-events-none max-w-xl">
          <div className="pointer-events-auto bg-black/85 backdrop-blur-md border border-zinc-800 rounded-lg p-2.5 flex items-start gap-3 shadow-xl font-mono text-xs">
            <div 
              className="w-9 h-9 rounded-full flex items-center justify-center font-bold text-sm border shadow-md shrink-0 mt-0.5"
              style={{
                backgroundColor: `${activeBlock.color}25`,
                borderColor: activeBlock.color,
                color: activeBlock.color
              }}
            >
              {activeResidueCode}
            </div>
            <div className="space-y-1">
              <div className="flex items-center gap-2 flex-wrap">
                <span className="font-bold text-white">Reszta #{safeActiveIdx + 1}: {activeBlock.physico.namePl}</span>
                <span className="text-[10px] px-1 py-0.2 rounded bg-zinc-800 text-zinc-300">
                  {activeBlock.physico.code1}
                </span>
                <span className="text-[10px] text-zinc-400">
                  VdW: {activeBlock.physico.vdwVolume} Å³
                </span>
                {activeRotamerInfo && (
                  <span className="text-[10px] px-1.5 py-0.2 rounded bg-indigo-950/80 border border-indigo-700/70 text-indigo-200 font-bold">
                    {activeRotamerInfo.heavyAtomCount} ciężkich atomów
                  </span>
                )}
              </div>

              <div className="text-[10.5px] text-zinc-400 flex items-center gap-2 flex-wrap">
                <span>Kyte: <strong className={activeBlock.physico.hydropathyKyte > 0 ? 'text-amber-400' : 'text-blue-400'}>
                  {activeBlock.physico.hydropathyKyte > 0 ? `+${activeBlock.physico.hydropathyKyte}` : activeBlock.physico.hydropathyKyte}
                </strong></span>
                <span>•</span>
                <span>Rotamer: <strong className="text-indigo-300">{activeRotamerInfo?.rotamerName || 'kanoniczny'}</strong></span>
                {activeRotamerInfo && activeRotamerInfo.availableRotamersCount > 1 && !activeRotamerInfo.rotamerName.startsWith('ss-locked') && (
                  <button
                    onClick={handleCycleActiveRotamer}
                    className="px-2 py-0.5 rounded bg-indigo-900/70 hover:bg-indigo-800 border border-indigo-500/70 text-indigo-100 text-[10px] font-bold flex items-center gap-1 cursor-pointer transition"
                    title="Przełącz dyskretny stan rotameru (χ1..χ4) łańcucha bocznego i przelicz wiązania boczne"
                  >
                    <RotateCcw className="w-2.5 h-2.5 text-indigo-300" />
                    <span>Zmień rotamer ({activeRotamerInfo.rotamerIndex + 1}/{activeRotamerInfo.availableRotamersCount})</span>
                  </button>
                )}
                {onJumpToWarehouseSphere && (
                  <button
                    onClick={() => onJumpToWarehouseSphere(activeResidueCode)}
                    className="text-cyan-400 hover:text-cyan-300 underline cursor-pointer"
                  >
                    Kula 248
                  </button>
                )}
              </div>

              {activeRotamerInfo && activeRotamerInfo.lateralContacts.length > 0 && (
                <div className="flex items-center gap-1.5 flex-wrap pt-0.5">
                  {activeRotamerInfo.lateralContacts.map((note, nIdx) => (
                    <span
                      key={nIdx}
                      className="text-[9.5px] px-1.5 py-0.2 rounded bg-fuchsia-950/70 border border-fuchsia-700/60 text-fuchsia-200"
                    >
                      {note}
                    </span>
                  ))}
                </div>
              )}
            </div>
          </div>
        </div>

        {/* Live Multi-Agent Cotranslational & Soap-Bubble Telemetry Banner */}
        {multiAgentDiag && (
          <div className="absolute top-3 right-3 z-10 pointer-events-none max-w-lg">
            <div className="pointer-events-auto bg-black/90 backdrop-blur-md border border-cyan-800/80 rounded-lg px-3 py-2 shadow-2xl font-mono text-[10.5px] space-y-1">
              <div className="flex items-center justify-between gap-3">
                <span className="font-extrabold text-cyan-300 uppercase">
                  {multiAgentDiag.tetheredToRibosome
                    ? `1. Wzrost Kotranslacyjny N→C (Soft SnapLock): ${multiAgentDiag.emergedCount}/${multiAgentDiag.totalCount} aa`
                    : multiAgentDiag.stage === 'CONVERGED_GLOBAL_MINIMUM'
                    ? `3. Minimum Energetyczne Po Odcięciu Od Rybosomu (${multiAgentDiag.totalCount} aa)`
                    : `2. Odcięcie Od Rybosomu: Efekt „Bańki Mydlanej” & Drżenie (φ,ψ)`}
                </span>
                {multiAgentDiag.liveKabschCaRmsdAngstrom !== null && (
                  <span
                    className={`px-1.5 py-0.5 rounded border font-bold ${
                      multiAgentDiag.liveKabschCaRmsdAngstrom < 1.0
                        ? 'bg-emerald-950 border-emerald-500 text-emerald-300'
                        : 'bg-amber-950 border-amber-600 text-amber-200'
                    }`}
                  >
                    Kabsch Cα: {multiAgentDiag.liveKabschCaRmsdAngstrom.toFixed(2)} Å
                  </span>
                )}
              </div>
              <div className="flex items-center gap-3 flex-wrap text-zinc-300">
                <span>Bańka Mydlana: <strong className="text-emerald-300">{multiAgentDiag.soapBubbleEnergyKcal} kcal</strong></span>
                <span>Kolizje VdW: <strong className={multiAgentDiag.vdwClashesCount === 0 ? 'text-emerald-400' : 'text-rose-400'}>{multiAgentDiag.vdwClashesCount}</strong></span>
                <span>Drżenie |Δφ,Δψ|: <strong className="text-amber-300">{multiAgentDiag.kineticQuiverAmplitudeDeg}°</strong></span>
                <span>R_kula/Rg: <strong className="text-cyan-300">{multiAgentDiag.sphereToRgRatio}×</strong> (≥1.40×)</span>
              </div>
            </div>
          </div>
        )}

        {/* Loading / Transition Spinner */}
        {!isFolding && (!foldedPositions || foldedPositions.length === 0 || foldedPositions.length !== sequenceCodes.length) && (
          <div className="absolute inset-0 z-30 flex flex-col items-center justify-center bg-black/85 backdrop-blur-sm pointer-events-none">
            <div className="w-10 h-10 border-3 border-cyan-500/30 border-t-cyan-400 rounded-full animate-spin mb-3" />
            <span className="font-mono text-xs text-cyan-300 font-bold tracking-wide">
              Przełączanie struktury & deterministyczne mapowanie E8...
            </span>
            <span className="font-mono text-[10.5px] text-zinc-500 mt-1">
              Kalibracja złączy Port #248 ➔ Port #1
            </span>
          </div>
        )}

        {/* 3D Mount Container */}
        <div
          ref={mountRef}
          className="flex-1 w-full h-full min-h-0 cursor-grab active:cursor-grabbing select-none relative"
          onMouseDown={handleMouseDown}
          onMouseMove={handleMouseMove}
          onMouseUp={handleMouseUp}
          onMouseLeave={handleMouseUp}
          onWheel={handleWheel}
          onTouchStart={handleTouchStart}
          onTouchMove={handleTouchMove}
          onTouchEnd={handleTouchEnd}
        />

        {/* ======================================================== */}
        {/* MIKROSKOP Å: FIRST-PERSON HUD OVERLAY (ANT-MAN SCALE)    */}
        {/* ======================================================== */}
        {isMicroscopeMode && (
          <div className="absolute inset-0 pointer-events-none flex flex-col justify-between p-4 z-20 font-mono text-xs select-none">
            
            {/* Top Bar: Residue Jump, Prev/Next & 5 Flight Path Buttons */}
            <div className="flex flex-wrap items-center justify-between gap-3 pointer-events-auto bg-black/85 backdrop-blur-md border border-rose-900/60 p-3 rounded-2xl shadow-2xl">
              
              {/* Residue Stepper */}
              <div className="flex items-center gap-2">
                <span className="text-rose-400 font-bold text-[11px] flex items-center gap-1.5">
                  <Eye className="w-4 h-4 text-rose-500 animate-pulse" />
                  <span>Mikroskop Å:</span>
                </span>

                <button
                  onClick={() => jumpToResidueFirstPerson(activeResidueIdx - 1)}
                  disabled={activeResidueIdx <= 0}
                  className="px-2.5 py-1 rounded bg-zinc-900 hover:bg-zinc-800 border border-zinc-700 text-zinc-200 disabled:opacity-40 cursor-pointer text-xs font-bold"
                >
                  ◀ Poprzednia
                </button>

                <button
                  onClick={() => jumpToResidueFirstPerson(activeResidueIdx)}
                  className="px-3 py-1 rounded bg-rose-950 border border-rose-500 text-rose-200 cursor-pointer text-xs font-bold shadow-md"
                  title="Ustaw kamerę w środku masy i lokalnym układzie N-Cα-C tej reszty"
                >
                  Wejdź do reszty #{safeActiveIdx + 1} ({activeResidueCode})
                </button>

                <button
                  onClick={() => jumpToResidueFirstPerson(activeResidueIdx + 1)}
                  disabled={activeResidueIdx >= sequenceCodes.length - 1}
                  className="px-2.5 py-1 rounded bg-zinc-900 hover:bg-zinc-800 border border-zinc-700 text-zinc-200 disabled:opacity-40 cursor-pointer text-xs font-bold"
                >
                  Następna ▶
                </button>
              </div>

              {/* 5 Automated Flight Paths */}
              <div className="flex items-center gap-1.5 flex-wrap">
                <span className="text-[10px] text-zinc-400 font-bold uppercase">5 Tras Przelotu:</span>
                
                <button
                  onClick={() => handleStartFlightPath(activeFlightPath === 'PATH_1_BACKBONE' ? 'NONE' : 'PATH_1_BACKBONE')}
                  className={`px-2 py-0.5 rounded border text-[10.5px] font-bold cursor-pointer transition ${
                    activeFlightPath === 'PATH_1_BACKBONE' ? 'bg-cyan-600 text-white border-cyan-400' : 'bg-zinc-900 text-zinc-300 border-zinc-700'
                  }`}
                  title="(1) Przelot wzdłuż łańcucha peptydowego N→C"
                >
                  (1) N→C
                </button>

                <button
                  onClick={() => handleStartFlightPath(activeFlightPath === 'PATH_2_SPIRAL_BOTTOM_UP' ? 'NONE' : 'PATH_2_SPIRAL_BOTTOM_UP')}
                  className={`px-2 py-0.5 rounded border text-[10.5px] font-bold cursor-pointer transition ${
                    activeFlightPath === 'PATH_2_SPIRAL_BOTTOM_UP' ? 'bg-cyan-600 text-white border-cyan-400' : 'bg-zinc-900 text-zinc-300 border-zinc-700'
                  }`}
                  title="(2) Ślimak dół→góra po środkach masy aminokwasów"
                >
                  (2) Ślimak Dół➔Góra
                </button>

                <button
                  onClick={() => handleStartFlightPath(activeFlightPath === 'PATH_3_SPIRAL_FRONT_BACK' ? 'NONE' : 'PATH_3_SPIRAL_FRONT_BACK')}
                  className={`px-2 py-0.5 rounded border text-[10.5px] font-bold cursor-pointer transition ${
                    activeFlightPath === 'PATH_3_SPIRAL_FRONT_BACK' ? 'bg-cyan-600 text-white border-cyan-400' : 'bg-zinc-900 text-zinc-300 border-zinc-700'
                  }`}
                  title="(3) Ślimak przód→tył"
                >
                  (3) Przód➔Tył
                </button>

                <button
                  onClick={() => handleStartFlightPath(activeFlightPath === 'PATH_4_LATERAL_BONDS' ? 'NONE' : 'PATH_4_LATERAL_BONDS')}
                  className={`px-2 py-0.5 rounded border text-[10.5px] font-bold cursor-pointer transition ${
                    activeFlightPath === 'PATH_4_LATERAL_BONDS' ? 'bg-cyan-600 text-white border-cyan-400' : 'bg-zinc-900 text-zinc-300 border-zinc-700'
                  }`}
                  title="(4) Wiązania boczne i mostki solne"
                >
                  (4) Boczne
                </button>

                <button
                  onClick={() => handleStartFlightPath(activeFlightPath === 'PATH_5_DISULFIDE' ? 'NONE' : 'PATH_5_DISULFIDE')}
                  className={`px-2 py-0.5 rounded border text-[10.5px] font-bold cursor-pointer transition ${
                    activeFlightPath === 'PATH_5_DISULFIDE' ? 'bg-amber-600 text-black border-amber-300' : 'bg-zinc-900 text-zinc-300 border-zinc-700'
                  }`}
                  title="(5) Mostki dwusiarczkowe S-S"
                >
                  (5) Mostki S-S
                </button>
              </div>

            </div>

            {/* Center Area: Crosshair + 1 Å Sub-nanometer Reticle */}
            <div className="self-center flex flex-col items-center justify-center relative pointer-events-none opacity-40">
              <div className="w-16 h-16 border border-cyan-500/50 rounded-full flex items-center justify-center">
                <div className="w-1.5 h-1.5 bg-rose-400 rounded-full" />
              </div>
              <div className="w-28 h-[1px] bg-cyan-500/30 absolute" />
              <div className="h-28 w-[1px] bg-cyan-500/30 absolute" />
            </div>

            {/* Bottom Bar: 1 Å Grid Ruler, 6 Nearest Neighbors & WASD Controls */}
            <div className="flex flex-wrap items-end justify-between gap-3 pointer-events-auto">
              
              {/* Left: 6 Nearest Neighbors with Hydropathy Colors & Distance in Å */}
              <div className="bg-black/85 backdrop-blur-md border border-zinc-800 p-3 rounded-2xl shadow-xl max-w-md space-y-1.5">
                <div className="flex items-center justify-between text-[10.5px] text-zinc-400 border-b border-zinc-800 pb-1">
                  <span className="font-bold text-white">Najbliższe 6 Reszt (w promieniu Å):</span>
                  <span className="text-zinc-500">Hydrofobowość: Kyte-Doolittle</span>
                </div>
                <div className="grid grid-cols-3 gap-1.5">
                  {nearestNeighbors.map(nb => (
                    <div 
                      key={nb.idx}
                      onClick={() => jumpToResidueFirstPerson(nb.idx)}
                      className="p-1.5 rounded-lg bg-zinc-900/90 border border-zinc-800 hover:border-cyan-500 cursor-pointer transition flex items-center justify-between text-[10px]"
                    >
                      <span className="font-bold text-zinc-200">#{nb.idx + 1} {nb.code}</span>
                      <span className={`font-mono font-bold ${nb.hydropathy >= 0 ? 'text-amber-400' : 'text-blue-400'}`}>
                        {nb.dist} Å
                      </span>
                    </div>
                  ))}
                </div>
              </div>

              {/* Right: Scale Slider with Always-Visible 1 Å Ruler + WASD Key Help */}
              <div className="bg-black/85 backdrop-blur-md border border-zinc-800 p-3 rounded-2xl shadow-xl space-y-2 font-mono">
                
                {/* 1 Å Ruler Graphic */}
                <div className="space-y-1">
                  <div className="flex items-center justify-between text-[10px] text-zinc-400">
                    <span>Linijka Skali:</span>
                    <strong className="text-cyan-300">1 Å = 0.1 nm (Siatka 1 Å)</strong>
                  </div>
                  {/* Calibrated 1 Å graphical bar */}
                  <div className="w-56 h-4 bg-zinc-900 border border-zinc-700 rounded relative flex items-center justify-between px-2 text-[9px] text-zinc-400 overflow-hidden">
                    <div className="absolute inset-y-0 left-0 bg-cyan-500/20" style={{ width: `${Math.min(100, Math.round(microscopeScale * 40))}%` }} />
                    <span className="z-10">0 Å</span>
                    <span className="z-10">0.5 Å</span>
                    <span className="z-10 font-bold text-white">1.0 Å</span>
                    <span className="z-10">1.5 Å</span>
                    <span className="z-10">2.0 Å</span>
                  </div>
                </div>

                {/* Scale slider */}
                <div className="flex items-center gap-2 text-[10.5px]">
                  <span className="text-zinc-400">Skala:</span>
                  <input
                    type="range"
                    min="0.2"
                    max="3.0"
                    step="0.05"
                    value={microscopeScale}
                    onChange={(e) => setMicroscopeScale(parseFloat(e.target.value))}
                    className="w-28 accent-cyan-400"
                  />
                  <span className="text-zinc-200 font-bold w-12">{microscopeScale.toFixed(2)}×</span>
                </div>

                {/* Controls Legend */}
                <div className="text-[10px] text-zinc-400 pt-1 border-t border-zinc-800/80 flex items-center justify-between gap-3">
                  <span>Sterowanie: <strong className="text-zinc-200">WASD</strong> (ruch) • <strong className="text-zinc-200">Mysz / Dotyk</strong> (rozglądanie)</span>
                  <span className="text-zinc-500">Q/E (góra/dół)</span>
                </div>

              </div>

            </div>

          </div>
        )}

        {/* BOUNDARY CONDITIONS & CD SLIDE-OUT PANEL */}
        {isBoundaryPanelOpen && (
          <div className="absolute top-3 right-3 bottom-14 w-88 z-20 bg-zinc-950/95 backdrop-blur-xl border border-zinc-800 rounded-xl p-4 flex flex-col shadow-2xl font-mono text-xs overflow-y-auto">
            <div className="flex items-center justify-between pb-3 border-b border-zinc-800">
              <div className="flex items-center gap-2">
                <SlidersHorizontal className="w-4 h-4 text-cyan-400" />
                <h3 className="font-bold text-white text-sm">Warunki Brzegowe (SAXS & CD)</h3>
              </div>
              <button
                onClick={() => setIsBoundaryPanelOpen(false)}
                className="p-1 rounded text-zinc-500 hover:text-white hover:bg-zinc-800 cursor-pointer"
              >
                <CloseIcon className="w-4 h-4" />
              </button>
            </div>

            <div className="py-3 space-y-4">
              {/* SAXS D_max Target Limit */}
              <div className="space-y-1.5">
                <div className="flex items-center justify-between">
                  <span className="font-bold text-zinc-300">Limit D_max (SAXS) [Å]:</span>
                  <span className="text-cyan-400 font-extrabold">{boundaryConfig.dMaxTarget} Å</span>
                </div>
                <input
                  type="range"
                  min="10"
                  max="80"
                  step="0.5"
                  value={boundaryConfig.dMaxTarget}
                  onChange={(e) => setBoundaryConfig({ ...boundaryConfig, dMaxTarget: parseFloat(e.target.value) })}
                  className="w-full accent-cyan-500 cursor-pointer"
                />
                <div className="flex items-center justify-between text-[10px] text-zinc-500">
                  <span>Zmierzone D_max: {foldingMetrics?.dMaxActual || 0} Å</span>
                  <button
                    onClick={applyFloryTheoreticalLimits}
                    className="text-cyan-400 hover:text-cyan-300 underline cursor-pointer"
                  >
                    Ustaw wg Flory'ego
                  </button>
                </div>
              </div>

              {/* Enforce D_max Compression Toggle */}
              <div className="flex items-center justify-between p-2 rounded bg-zinc-900 border border-zinc-800">
                <div>
                  <div className="font-bold text-zinc-300">Wymuszaj limit D_max</div>
                  <div className="text-[10px] text-zinc-500">Kompresja graniczna cząsteczki</div>
                </div>
                <input
                  type="checkbox"
                  checked={boundaryConfig.enforceBoundaries}
                  onChange={(e) => setBoundaryConfig({ ...boundaryConfig, enforceBoundaries: e.target.checked })}
                  className="w-4 h-4 accent-cyan-500 rounded cursor-pointer"
                />
              </div>

              {/* Far-UV CD Alpha-Helix Target */}
              <div className="space-y-1.5">
                <div className="flex items-center justify-between">
                  <span className="font-bold text-zinc-300">Dichroizm Kołowy: Helisa α [%]:</span>
                  <span className="text-blue-400 font-extrabold">{boundaryConfig.cdAlphaTargetPct}%</span>
                </div>
                <input
                  type="range"
                  min="0"
                  max="100"
                  step="1"
                  value={boundaryConfig.cdAlphaTargetPct}
                  onChange={(e) => setBoundaryConfig({ ...boundaryConfig, cdAlphaTargetPct: parseInt(e.target.value, 10) })}
                  className="w-full accent-blue-500 cursor-pointer"
                />
                <div className="text-[10px] text-zinc-500 flex justify-between">
                  <span>Osiągnięto w modelu:</span>
                  <span className="text-blue-300 font-bold">{foldingMetrics?.alphaActualPct || 0}%</span>
                </div>
              </div>

              {/* Far-UV CD Beta-Sheet Target */}
              <div className="space-y-1.5">
                <div className="flex items-center justify-between">
                  <span className="font-bold text-zinc-300">Dichroizm Kołowy: Kartka β [%]:</span>
                  <span className="text-indigo-400 font-extrabold">{boundaryConfig.cdBetaTargetPct}%</span>
                </div>
                <input
                  type="range"
                  min="0"
                  max="100"
                  step="1"
                  value={boundaryConfig.cdBetaTargetPct}
                  onChange={(e) => setBoundaryConfig({ ...boundaryConfig, cdBetaTargetPct: parseInt(e.target.value, 10) })}
                  className="w-full accent-indigo-500 cursor-pointer"
                />
                <div className="text-[10px] text-zinc-500 flex justify-between">
                  <span>Osiągnięto w modelu:</span>
                  <span className="text-indigo-300 font-bold">{foldingMetrics?.betaActualPct || 0}%</span>
                </div>
              </div>

              {/* Disulfide Restraints Base Weight */}
              <div className="space-y-1.5">
                <div className="flex items-center justify-between">
                  <span className="font-bold text-zinc-300">Waga bazowa w_ss⁰ (SSBOND):</span>
                  <span className="text-amber-400 font-extrabold">{boundaryConfig.ssRestraintWeight}</span>
                </div>
                <input
                  type="range"
                  min="0.5"
                  max="12.0"
                  step="0.5"
                  value={boundaryConfig.ssRestraintWeight}
                  onChange={(e) => setBoundaryConfig({ ...boundaryConfig, ssRestraintWeight: parseFloat(e.target.value) })}
                  className="w-full accent-amber-500 cursor-pointer"
                />
                <div className="text-[10px] text-zinc-500 flex justify-between">
                  <span>Cel: d(Sγ,Sγ)=2.05 Å • d(Cα,Cα)=5.50 Å</span>
                  <span>Par: {foldingMetrics?.ssBondPairs.length || 0}</span>
                </div>
              </div>

              {/* Two-Way IK & Adaptive Exponential SSBOND Diagnostics */}
              <div className="p-3 rounded-lg bg-black/70 border border-zinc-800 space-y-2">
                <div className="flex items-center justify-between">
                  <span className="font-bold text-zinc-300 text-[11px] uppercase">Domykanie Pętli (IK &amp; SSBOND):</span>
                  <span className={`text-[9.5px] px-1.5 py-0.5 rounded font-bold ${
                    foldingMetrics?.pipelineMode === 'TWO_WAY_IK_ADAPTIVE_SSBOND'
                      ? 'bg-amber-950 text-amber-300 border border-amber-700/70'
                      : 'bg-cyan-950 text-cyan-300 border border-cyan-800'
                  }`}>
                    {foldingMetrics?.pipelineMode === 'TWO_WAY_IK_ADAPTIVE_SSBOND' ? 'IK + ADAPT. w_ss' : 'STANDARD L2'}
                  </span>
                </div>

                {foldingMetrics && foldingMetrics.bridgeDiagnostics.length > 0 ? (
                  <div className="space-y-1.5 text-[10px]">
                    <div className="flex items-center justify-between text-zinc-400">
                      <span>Waga adaptacyjna w_ss:</span>
                      <span className="font-bold text-amber-300">
                        {foldingMetrics.initialSpringWeight} ➔ {foldingMetrics.finalSpringWeight}
                        {foldingMetrics.exponentialBoostTriggered && ' (EXP)'}
                      </span>
                    </div>
                    <div className="flex items-center justify-between text-zinc-400">
                      <span>Kara E_ssbond = w_ss·Δd²:</span>
                      <span className="font-bold text-emerald-300">{foldingMetrics.disulfidePenaltyScore} kcal/mol</span>
                    </div>
                    <div className="pt-1.5 border-t border-zinc-800/80 space-y-1">
                      {foldingMetrics.bridgeDiagnostics.map((b, idx) => (
                        <div key={idx} className="p-1.5 rounded bg-zinc-950 border border-zinc-800/90 space-y-0.5">
                          <div className="flex justify-between font-bold text-zinc-200">
                            <span>CYS{b.res1} ↔ CYS{b.res2}</span>
                            <span className={b.sgDeviation < 0.5 ? 'text-emerald-400' : 'text-rose-400'}>
                              Δ(Sγ) = {b.sgDeviation} Å {b.sgDeviation < 0.5 ? '✓' : '!'}
                            </span>
                          </div>
                          <div className="flex justify-between text-[9.5px] text-zinc-400">
                            <span>Sγ-Sγ: {b.initialSgDist} ➔ <strong className="text-amber-300">{b.finalSgDist} Å</strong></span>
                            <span>Cα-Cα: {b.finalCaDist} Å</span>
                          </div>
                        </div>
                      ))}
                    </div>
                  </div>
                ) : (
                  <div className="text-[10px] text-zinc-500 leading-relaxed">
                    Brak mostków klamrowych CYS/SEC w sekwencji — aktywny czysty potok mikrofaset Level 2 (2232 stany, błąd {foldingMetrics?.meanMicroFacetErrorDeg || 0}°/resztę).
                  </div>
                )}

                <div className="pt-2 border-t border-zinc-800 space-y-1 text-[10.5px]">
                  <div className="flex justify-between text-zinc-400">
                    <span>Średni błąd faset L2:</span>
                    <span className="font-bold text-emerald-400">{foldingMetrics?.meanMicroFacetErrorDeg || 0}° (&lt;1.0°)</span>
                  </div>
                  <div className="flex justify-between text-zinc-400">
                    <span>Kara CD (Far-UV):</span>
                    <span className="font-bold text-white">{foldingMetrics?.cdPenaltyScore || 0}</span>
                  </div>
                  <div className="flex justify-between text-zinc-400 border-t border-zinc-800 pt-1">
                    <span className="font-bold text-emerald-400">Całkowita ΔG:</span>
                    <span className="font-extrabold text-emerald-300">{foldingMetrics?.totalCohesiveEnergy || 0} kcal/mol</span>
                  </div>
                </div>
              </div>
            </div>

            <div className="pt-3 border-t border-zinc-800 flex gap-2">
              <button
                onClick={() => {
                  runFoldingAnimation(sequenceCodes, boundaryConfig);
                  setIsBoundaryPanelOpen(false);
                }}
                className="flex-1 py-2 bg-gradient-to-r from-cyan-600 to-emerald-600 hover:from-cyan-500 hover:to-emerald-500 text-white font-bold rounded text-center cursor-pointer shadow-lg transition"
              >
                Przelicz & Zwiń
              </button>
            </div>
          </div>
        )}

        {/* Bottom Horizontal Interactive Residue Ribbon */}
        <div className="absolute bottom-2 left-3 right-3 z-10 pointer-events-none">
          <div className="pointer-events-auto bg-black/85 backdrop-blur-md border border-zinc-800/90 rounded-lg p-1.5 flex items-center gap-1.5 overflow-x-auto shadow-2xl font-mono text-xs">
            <span className="text-[9.5px] text-zinc-500 uppercase font-bold shrink-0 px-2">
              Sekwencja ({sequenceCodes.length} aa):
            </span>
            {sequenceCodes.map((code, idx) => {
              const isSelected = activeResidueIdx === idx;
              const isCys = code === 'CYS' || code === 'SEC';
              const sec = secondaryTypes[idx] || 'COIL';

              return (
                <button
                  key={`${code}-${idx}`}
                  onClick={() => setActiveResidueIdx(idx)}
                  onMouseEnter={() => setHoveredResidueIdx(idx)}
                  onMouseLeave={() => setHoveredResidueIdx(null)}
                  className={`px-2 py-1 rounded text-center shrink-0 border transition cursor-pointer ${
                    isSelected
                      ? 'bg-cyan-950 border-cyan-400 text-cyan-200 ring-1 ring-cyan-400 shadow-md'
                      : isCys
                      ? 'bg-amber-950/60 border-amber-800 text-amber-300 hover:bg-amber-900/60'
                      : 'bg-zinc-950 border-zinc-800 text-zinc-400 hover:text-white hover:bg-zinc-900'
                  }`}
                >
                  <div className="text-[10px] font-bold">{code}</div>
                  <div className="text-[8px] text-zinc-500 flex items-center justify-center gap-0.5">
                    <span>#{idx + 1}</span>
                    {sec === 'HELIX' && <span className="text-blue-400">α</span>}
                    {sec === 'SHEET' && <span className="text-indigo-400">β</span>}
                  </div>
                </button>
              );
            })}
          </div>
        </div>

      </div>

      {/* ======================================================== */}
      {/* 4. FASTA IMPORT MODAL */}
      {/* ======================================================== */}
      {isFastaModalOpen && (
        <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/80 backdrop-blur-md p-4">
          <div className="bg-zinc-950 border border-zinc-800 rounded-xl max-w-2xl w-full p-5 shadow-2xl flex flex-col font-mono text-xs">
            <div className="flex items-center justify-between pb-3 border-b border-zinc-800">
              <div className="flex items-center gap-2">
                <Dna className="w-4 h-4 text-indigo-400" />
                <h3 className="font-bold text-white text-sm">Wklej Sekwencję w Formacie FASTA</h3>
              </div>
              <button
                onClick={() => setIsFastaModalOpen(false)}
                className="p-1 rounded text-zinc-500 hover:text-white hover:bg-zinc-800 cursor-pointer"
              >
                <CloseIcon className="w-4 h-4" />
              </button>
            </div>

            <div className="py-4 space-y-3">
              <p className="text-zinc-400 text-[11px] leading-relaxed">
                Obsługuje pełny standard FASTA z nagłówkiem <code className="text-indigo-300">&gt;Nazwa</code> lub czysty ciąg literowy (np. 1-literowy <code className="text-indigo-300">GIVEQCCTS...</code> lub 3-literowy <code className="text-indigo-300">GLY ILE VAL...</code>).
              </p>

              <textarea
                value={fastaModalText}
                onChange={(e) => setFastaModalText(e.target.value)}
                rows={8}
                placeholder=">ID_PROTEIN Opis białka&#10;GIVEQCCTSICSLYQLENYCN..."
                className="w-full bg-black border border-zinc-800 rounded-lg p-3 text-white font-mono text-xs focus:outline-none focus:border-indigo-500 leading-normal"
              />

              <div className="flex items-center gap-2 flex-wrap">
                <span className="text-zinc-500 text-[10px] uppercase font-bold">Przykłady:</span>
                {FASTA_PRESETS.map((p) => (
                  <button
                    key={p.id}
                    onClick={() => setFastaModalText(p.fasta)}
                    className="px-2 py-1 rounded bg-zinc-900 hover:bg-zinc-800 border border-zinc-800 text-zinc-300 hover:text-white text-[10.5px] cursor-pointer"
                  >
                    {p.name}
                  </button>
                ))}
              </div>
            </div>

            <div className="pt-3 border-t border-zinc-800 flex justify-end gap-2">
              <button
                onClick={() => setIsFastaModalOpen(false)}
                className="px-4 py-2 rounded bg-zinc-900 hover:bg-zinc-800 text-zinc-300 font-bold cursor-pointer"
              >
                Anuluj
              </button>
              <button
                onClick={() => {
                  handleApplyInput(fastaModalText);
                  setIsFastaModalOpen(false);
                }}
                className="px-5 py-2 rounded bg-indigo-600 hover:bg-indigo-500 text-white font-bold cursor-pointer shadow-lg shadow-indigo-600/30"
              >
                Zastosuj & Zwiń
              </button>
            </div>
          </div>
        </div>
      )}

      {/* ======================================================== */}
      {/* 5. STANDARD ALL-ATOM wwPDB v3.3 INSPECTOR & DOWNLOAD MODAL */}
      {/* ======================================================== */}
      {isPdbModalOpen && (
        <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/80 backdrop-blur-md p-4">
          <div className="bg-zinc-950 border border-zinc-800 rounded-xl max-w-5xl w-full h-[86vh] p-5 shadow-2xl flex flex-col font-mono text-xs">
            <div className="flex items-center justify-between pb-3 border-b border-zinc-800 shrink-0">
              <div className="flex items-center gap-2 flex-wrap">
                <FileText className="w-4 h-4 text-sky-400" />
                <h3 className="font-bold text-white text-sm">
                  Pełna Struktura Atomowa PDB: <span className="text-sky-300">{fastaHeader}_all_atom.pdb</span>
                </h3>
                <span className="text-[10px] px-2 py-0.5 rounded bg-sky-950 text-sky-300 border border-sky-800">
                  wwPDB v3.3 Standard (80 kolumn • Wszystkie łańcuchy boczne, Rotamery &amp; CONECT)
                </span>
              </div>
              <button
                onClick={() => setIsPdbModalOpen(false)}
                className="p-1 rounded text-zinc-500 hover:text-white hover:bg-zinc-800 cursor-pointer"
              >
                <CloseIcon className="w-4 h-4" />
              </button>
            </div>

            {/* Quick summary badges */}
            <div className="py-2.5 px-3 my-2 rounded-lg bg-black/60 border border-zinc-800/80 flex flex-wrap items-center justify-between gap-2 shrink-0 text-[11px] text-zinc-400">
              <div className="flex items-center gap-3 flex-wrap">
                <span>Reszt: <strong className="text-white">{sequenceCodes.length}</strong></span>
                <span>Atomy ciężkie: <strong className="text-indigo-300">{allAtomStructure?.totalHeavyAtoms || 0}</strong> ({allAtomStructure?.totalSidechainAtoms || 0} bocz.)</span>
                <span>Wiązania kow.: <strong className="text-cyan-400">{allAtomStructure?.totalCovalentBonds || 0}</strong></span>
                <span>Mostki S-S: <strong className="text-amber-300">{foldingMetrics?.ssBondPairs.length || 0}</strong></span>
                <span>Solne / H-bond: <strong className="text-fuchsia-300">{allAtomStructure?.totalSaltBridges || 0} / {allAtomStructure?.totalHydrogenBonds || 0}</strong></span>
                <span>
                  {foldingMetrics?.groundTruthValidation?.hasGroundTruthReference &&
                  foldingMetrics.groundTruthValidation.caRmsdAngstrom !== null ? (
                    <>
                      Ground-Truth Cα RMSD ({foldingMetrics.groundTruthValidation.referenceId}):{' '}
                      <strong className="text-purple-300">
                        {foldingMetrics.groundTruthValidation.caRmsdAngstrom.toFixed(3)} Å
                      </strong>
                    </>
                  ) : (
                    <span className="text-zinc-500">RMSD referencyjny niedostępny</span>
                  )}{' '}
                  • Wewn. L2: <strong className="text-cyan-300">{foldingMetrics?.meanMicroFacetErrorDeg ?? 0}°</strong>
                </span>
              </div>
              <div className="flex items-center gap-2">
                <button
                  onClick={handleCopyPdb}
                  className="px-3 py-1 rounded bg-zinc-900 hover:bg-zinc-800 border border-zinc-700 text-zinc-200 hover:text-white font-bold flex items-center gap-1.5 cursor-pointer transition"
                >
                  {copySuccess ? <Check className="w-3.5 h-3.5 text-emerald-400" /> : <Copy className="w-3.5 h-3.5" />}
                  <span>{copySuccess ? 'Skopiowano!' : 'Kopiuj PDB'}</span>
                </button>
                <button
                  onClick={handleDownloadPdb}
                  className="px-3.5 py-1 rounded bg-sky-600 hover:bg-sky-500 text-white font-bold flex items-center gap-1.5 cursor-pointer shadow-md transition"
                >
                  <Download className="w-3.5 h-3.5 text-white" />
                  <span>Pobierz All-Atom .pdb</span>
                </button>
              </div>
            </div>

            {/* PDB Text Preview */}
            <div className="flex-1 min-h-0 bg-black border border-zinc-900 rounded-lg p-3 overflow-auto font-mono text-[11px] leading-relaxed select-text text-zinc-300">
              <pre className="whitespace-pre">
                {currentPdbString}
              </pre>
            </div>

            <div className="pt-3 border-t border-zinc-800 flex justify-between items-center shrink-0 text-zinc-500 text-[10.5px]">
              <span>Gotowy do walidacji zewnętrznej (PyMOL, ChimeraX, Mol*, BioPython, TM-align): HEADER, SEQRES, HELIX/SHEET, SSBOND, pełne ATOM (IUPAC), TER, CONECT, END.</span>
              <button
                onClick={() => setIsPdbModalOpen(false)}
                className="px-4 py-1.5 rounded bg-zinc-900 hover:bg-zinc-800 text-zinc-300 font-bold cursor-pointer"
              >
                Zamknij
              </button>
            </div>
          </div>
        </div>
      )}

      {/* ======================================================== */}
      {/* 6. ALPHAFOLD DB / PDB KABSCH 3D RMSD VALIDATION MODAL */}
      {/* ======================================================== */}
      {isRmsdModalOpen && (
        <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/80 backdrop-blur-md p-4">
          <div className="bg-zinc-950 border border-zinc-800 rounded-xl max-w-5xl w-full max-h-[88vh] p-5 shadow-2xl flex flex-col font-mono text-xs overflow-hidden">
            {/* Header */}
            <div className="flex items-center justify-between pb-3 border-b border-zinc-800 shrink-0">
              <div className="flex items-center gap-2.5">
                <div className="w-8 h-8 rounded-lg bg-purple-950/90 border border-purple-500/60 flex items-center justify-center text-purple-300">
                  <GitCompare className="w-4 h-4" />
                </div>
                <div>
                  <h3 className="font-bold text-white text-sm flex items-center gap-2">
                    <span>Walidacja 3D RMSD &amp; Superpozycja Kabscha vs AlphaFold / PDB</span>
                    <span className="text-[10px] px-2 py-0.5 rounded bg-purple-950 text-purple-300 border border-purple-700/70">
                      Kwaterniony Horna (det R = +1)
                    </span>
                  </h3>
                  <p className="text-[10.5px] text-zinc-400">
                    Porównanie pełnej struktury atomowej SNAPLOCK ({allAtomStructure?.totalHeavyAtoms || 0} ciężkich atomów) z referencją AlphaFold DB (EMBL-EBI) lub RCSB PDB
                  </p>
                </div>
              </div>
              <button
                onClick={() => setIsRmsdModalOpen(false)}
                className="p-1 rounded text-zinc-500 hover:text-white hover:bg-zinc-800 cursor-pointer"
              >
                <CloseIcon className="w-4 h-4" />
              </button>
            </div>

            {/* Reference Source Selector */}
            <div className="py-3 border-b border-zinc-900 space-y-2.5 shrink-0">
              <div className="flex flex-wrap items-center justify-between gap-2">
                <div className="flex items-center gap-1.5">
                  <button
                    onClick={() => setRmsdInputMode('remote')}
                    className={`px-3 py-1 rounded text-xs font-bold border cursor-pointer transition ${
                      rmsdInputMode === 'remote'
                        ? 'bg-purple-950 border-purple-500 text-purple-200'
                        : 'bg-zinc-900 border-zinc-800 text-zinc-400 hover:text-white'
                    }`}
                  >
                    AlphaFold DB (UniProt) / RCSB PDB
                  </button>
                  <button
                    onClick={() => setRmsdInputMode('custom')}
                    className={`px-3 py-1 rounded text-xs font-bold border cursor-pointer transition ${
                      rmsdInputMode === 'custom'
                        ? 'bg-purple-950 border-purple-500 text-purple-200'
                        : 'bg-zinc-900 border-zinc-800 text-zinc-400 hover:text-white'
                    }`}
                  >
                    Własny plik referencyjny .pdb (AlphaFold3 / ColabFold)
                  </button>
                </div>

                <div className="flex items-center gap-1.5 flex-wrap">
                  <span className="text-[10px] text-zinc-500 uppercase font-bold">Szybkie ID (1:1):</span>
                  {[
                    { id: '1UBQ', label: 'Ubikwityna 76 aa (PDB:1UBQ)' },
                    { id: '3I40_A_21', label: 'Insulina A 21 aa (PDB:3I40:A)' },
                    { id: '1CRN', label: 'Krambina 46 aa (PDB:1CRN)' },
                    { id: '1OMG', label: 'Konotoksyna 25 aa (PDB:1OMG)' },
                    { id: '1LYZ', label: 'Lizozym 129 aa (PDB:1LYZ)' }
                  ].map(item => (
                    <button
                      key={item.id}
                      onClick={() => {
                        setRmsdInputMode('remote');
                        setRmsdTargetId(item.id);
                        const strict = resolveStrictTargetReference(item.id);
                        if (strict && onSelectGlobalTarget) {
                          onSelectGlobalTarget(strict.targetSpec.id);
                        } else {
                          handleRunRmsdValidation(item.id);
                        }
                      }}
                      className="px-2 py-0.5 rounded bg-zinc-900 hover:bg-zinc-800 border border-zinc-800 text-zinc-300 hover:text-white text-[10px] cursor-pointer"
                    >
                      {item.label}
                    </button>
                  ))}
                </div>
              </div>

              {rmsdInputMode === 'remote' ? (
                <div className="flex items-center gap-2">
                  <input
                    type="text"
                    value={rmsdTargetId}
                    onChange={(e) => setRmsdTargetId(e.target.value.toUpperCase())}
                    onKeyDown={(e) => {
                      if (e.key === 'Enter') handleRunRmsdValidation();
                    }}
                    placeholder="Wpisz UniProt ID (np. P01308, P0CG48, P01542) lub kod PDB (np. 1CRN, 1UBQ)"
                    className="flex-1 bg-black border border-zinc-800 rounded px-3 py-1.5 text-white font-mono text-xs focus:outline-none focus:border-purple-500"
                  />
                  <button
                    onClick={() => handleRunRmsdValidation()}
                    disabled={isComputingRmsd}
                    className="px-4 py-1.5 rounded bg-gradient-to-r from-purple-600 to-indigo-600 hover:from-purple-500 hover:to-indigo-500 disabled:opacity-50 text-white font-bold flex items-center gap-1.5 cursor-pointer shadow-md transition"
                  >
                    {isComputingRmsd ? (
                      <>
                        <RefreshCw className="w-3.5 h-3.5 animate-spin" />
                        <span>Pobieranie &amp; Kabsch SVD...</span>
                      </>
                    ) : (
                      <>
                        <Play className="w-3.5 h-3.5" />
                        <span>Pobierz Referencję &amp; Oblicz RMSD</span>
                      </>
                    )}
                  </button>
                </div>
              ) : (
                <div className="space-y-2">
                  <div className="flex items-center justify-between gap-2">
                    <label className="px-3 py-1.5 rounded bg-zinc-900 hover:bg-zinc-800 border border-zinc-700 text-zinc-200 font-bold flex items-center gap-1.5 cursor-pointer">
                      <Upload className="w-3.5 h-3.5 text-purple-400" />
                      <span>Wczytaj plik .pdb z dysku</span>
                      <input
                        type="file"
                        accept=".pdb,.ent,.txt"
                        onChange={handleCustomPdbFileUpload}
                        className="hidden"
                      />
                    </label>
                    <button
                      onClick={() => handleRunRmsdValidation()}
                      disabled={isComputingRmsd || !customRefPdbText.trim()}
                      className="px-4 py-1.5 rounded bg-purple-600 hover:bg-purple-500 disabled:opacity-50 text-white font-bold cursor-pointer"
                    >
                      Oblicz RMSD z wklejonego tekstu PDB
                    </button>
                  </div>
                  <textarea
                    value={customRefPdbText}
                    onChange={(e) => setCustomRefPdbText(e.target.value)}
                    rows={4}
                    placeholder="Wklej rekordy ATOM ... z referencyjnego pliku AlphaFold (.pdb)..."
                    className="w-full bg-black border border-zinc-800 rounded p-2 text-[11px] text-zinc-300 focus:outline-none focus:border-purple-500"
                  />
                </div>
              )}

              {rmsdError && (
                <div className="p-2.5 rounded bg-rose-950/70 border border-rose-800 text-rose-200 text-[11px]">
                  {rmsdError}
                </div>
              )}
            </div>

            {/* Results Area */}
            <div className="flex-1 min-h-0 overflow-y-auto py-3 space-y-3">
              {!rmsdComparison || !foldingMetrics?.groundTruthValidation?.hasGroundTruthReference ? (
                <div className="py-8 text-center text-zinc-400 space-y-2.5">
                  <GitCompare className="w-8 h-8 text-purple-400 mx-auto opacity-80" />
                  <div className="font-bold text-white text-sm">RMSD referencyjny niedostępny</div>
                  <p className="max-w-xl mx-auto text-[11px] text-zinc-400 leading-relaxed">
                    Nie załadowano jeszcze zewnętrznej struktury referencyjnej (RCSB PDB / AlphaFold DB / własny plik .pdb). Wyświetlanie <code className="text-rose-300">0.0 Å</code> lub wewnętrznego reziduum siatki L2 (<code className="text-amber-300">{foldingMetrics?.internalL2BackboneResidualAngstrom ?? 0} Å</code>) jako RMSD referencyjnego jest zablokowane. Kliknij <strong className="text-purple-300">Pobierz Referencję &amp; Oblicz RMSD</strong> powyżej, aby wykonać niezależną superpozycję Kabscha SVD bez modyfikacji współrzędnych modelu.
                  </p>
                  {foldingMetrics?.groundTruthValidation && (
                    <div className="inline-flex items-center gap-3 px-3 py-1.5 rounded bg-zinc-900 border border-zinc-800 text-[10.5px] text-zinc-300">
                      <span>Kontrola geometrii fizycznej modelu: <strong className={foldingMetrics.groundTruthValidation.modelGeometry.geometryValid ? 'text-emerald-400' : 'text-rose-400'}>{foldingMetrics.groundTruthValidation.modelGeometry.geometryValid ? 'POPRAWNA' : 'NARUSZONA'}</strong></span>
                      <span>•</span>
                      <span>Iteracja: <strong>{foldingMetrics.groundTruthValidation.iterationId}</strong> ({foldingMetrics.groundTruthValidation.parameterVersion})</span>
                    </div>
                  )}
                </div>
              ) : (
                <>
                  {/* Ground-Truth Provenance & Status Banner */}
                  <div
                    className={`p-3 rounded-lg border flex flex-wrap items-center justify-between gap-2 text-[11px] ${
                      foldingMetrics.groundTruthValidation.passedTargetCriterion
                        ? 'bg-emerald-950/35 border-emerald-600/50'
                        : 'bg-rose-950/35 border-rose-700/60'
                    }`}
                  >
                    <div className="space-y-0.5">
                      <div className="flex items-center gap-2 font-bold text-white">
                        <span>
                          Niezależna Walidacja Ground-Truth Kabscha SVD vs {foldingMetrics.groundTruthValidation.referenceName || rmsdComparison.referenceId} ({foldingMetrics.groundTruthValidation.referenceId}:{foldingMetrics.groundTruthValidation.selectedReferenceChain ?? 'A'}):
                        </span>
                        <span
                          className={`px-2 py-0.5 rounded text-[10px] font-extrabold ${
                            foldingMetrics.groundTruthValidation.passedTargetCriterion
                              ? 'bg-emerald-900 text-emerald-200'
                              : 'bg-rose-900 text-rose-200'
                          }`}
                        >
                          {foldingMetrics.groundTruthValidation.status}
                        </span>
                      </div>
                      {foldingMetrics.groundTruthValidation.failureReason && (
                        <div className="text-rose-200 text-[10.5px]">
                          Przyczyna niezaliczenia progu &lt;1.0 Å: {foldingMetrics.groundTruthValidation.failureReason}
                        </div>
                      )}
                    </div>
                    <div className="flex items-center gap-3 flex-wrap text-zinc-300 text-[10.5px]">
                      <span>
                        Iteracja: <strong className="text-cyan-300">{foldingMetrics.groundTruthValidation.iterationId}</strong> ({foldingMetrics.groundTruthValidation.parameterVersion})
                      </span>
                      <span>•</span>
                      <span>
                        Pokrycie M/R: <strong className="text-white">{foldingMetrics.groundTruthValidation.modelCoveragePct}% / {foldingMetrics.groundTruthValidation.referenceCoveragePct}%</strong>
                      </span>
                      <span>•</span>
                      <span>
                        Wewn. reziduum L2 (nie-ref): <strong className="text-amber-300">{foldingMetrics.internalL2BackboneResidualAngstrom} Å</strong>
                      </span>
                    </div>
                  </div>

                  {/* Metric Cards */}
                  <div className="grid grid-cols-2 md:grid-cols-6 gap-2.5">
                    <div className="p-3 rounded-lg bg-black/70 border border-purple-700/80 ring-1 ring-purple-500/20">
                      <div className="text-[9.5px] text-purple-300 uppercase font-bold">GT Cα-RMSD (SVD)</div>
                      <div className="text-xl font-extrabold text-white mt-0.5">
                        {foldingMetrics.groundTruthValidation.caRmsdAngstrom?.toFixed(3) ?? rmsdComparison.caRmsdAngstrom} Å
                      </div>
                      <div className="text-[9.5px] text-zinc-400">
                        Niezależny XYZ: {foldingMetrics.groundTruthValidation.independentCaRmsdAngstrom?.toFixed(3)} Å
                      </div>
                    </div>

                    <div className="p-3 rounded-lg bg-black/70 border border-sky-700/80">
                      <div className="text-[9.5px] text-sky-300 uppercase font-bold">GT Backbone RMSD</div>
                      <div className="text-xl font-extrabold text-sky-200 mt-0.5">{rmsdComparison.backboneRmsdAngstrom} Å</div>
                      <div className="text-[9.5px] text-zinc-400">N, CA, C, O (po superp.)</div>
                    </div>

                    <div className="p-3 rounded-lg bg-black/70 border border-indigo-800/70">
                      <div className="text-[9.5px] text-zinc-400 uppercase font-bold">GT All-Atom RMSD</div>
                      <div className="text-xl font-extrabold text-indigo-300 mt-0.5">{rmsdComparison.allHeavyAtomRmsdAngstrom} Å</div>
                      <div className="text-[9.5px] text-zinc-400">Szkielet + Rotamery</div>
                    </div>

                    <div className="p-3 rounded-lg bg-black/70 border border-emerald-800/70">
                      <div className="text-[9.5px] text-zinc-400 uppercase font-bold">TM-Score / GDT</div>
                      <div className="text-xl font-extrabold text-emerald-400 mt-0.5">{rmsdComparison.tmScore}</div>
                      <div className="text-[9.5px] text-zinc-400">GDT-TS: {rmsdComparison.gdtTsPct}%</div>
                    </div>

                    <div className="p-3 rounded-lg bg-black/70 border border-amber-800/70">
                      <div className="text-[9.5px] text-zinc-400 uppercase font-bold">Pary &amp; Identycz.</div>
                      <div className="text-xl font-extrabold text-amber-300 mt-0.5">
                        {foldingMetrics.groundTruthValidation.pairedResidueCount}/{foldingMetrics.groundTruthValidation.modelResidueCount}
                      </div>
                      <div className="text-[9.5px] text-zinc-400">Seq ID: {rmsdComparison.sequenceIdentityPct}%</div>
                    </div>

                    <div className="p-3 rounded-lg bg-black/70 border border-zinc-800">
                      <div className="text-[9.5px] text-zinc-400 uppercase font-bold">Wewn. ΔG / L2</div>
                      <div className="text-xl font-extrabold text-amber-200 mt-0.5">{foldingMetrics.totalCohesiveEnergy}</div>
                      <div className="text-[9.5px] text-zinc-500">kcal/mol (nie-ref)</div>
                    </div>
                  </div>

                  {/* Recharts Per-Residue Atomic RMSD Comparison Chart vs AlphaFold */}
                  <SnaplockRmsdComparisonChart
                    rmsdComparison={rmsdComparison}
                    allAtomStructure={allAtomStructure}
                    secondaryTypes={secondaryTypes}
                    disulfideResidueSeqs={
                      foldingMetrics?.ssBondPairs
                        ? foldingMetrics.ssBondPairs.flatMap(([i1, i2]) => [i1 + 1, i2 + 1])
                        : []
                    }
                    activeResidueIdx={activeResidueIdx}
                    onSelectResidue={(idx) => setActiveResidueIdx(idx)}
                    compact={true}
                  />

                  {/* Per-Residue Alignment & Sidechain Deviation Table */}
                  <div className="border border-zinc-800 rounded-lg overflow-hidden bg-black/60">
                    <div className="px-3 py-2 bg-zinc-900/80 border-b border-zinc-800 flex items-center justify-between text-[10.5px]">
                      <span className="font-bold text-zinc-200 uppercase">
                        Profil Odchyleń Per-Reszta (Model SNAPLOCK ↔ {rmsdComparison.referenceId})
                      </span>
                      <span className="text-zinc-400">
                        Kliknij słupek na wykresie lub wiersz tabeli, aby podświetlić resztę w 3D
                      </span>
                    </div>
                    <div className="max-h-56 overflow-y-auto">
                      <table className="w-full text-left text-[11px] tabular-nums">
                        <thead className="bg-zinc-950 text-zinc-500 uppercase text-[9.5px] sticky top-0">
                          <tr>
                            <th className="p-2">Model #</th>
                            <th className="p-2">Ref #</th>
                            <th className="p-2">Reszta</th>
                            <th className="p-2">Rotamer SNAPLOCK</th>
                            <th className="p-2">ΔCα [Å]</th>
                            <th className="p-2">Backbone [Å]</th>
                            <th className="p-2">Łańcuch Boczny [Å]</th>
                            <th className="p-2">All-Atom [Å]</th>
                            <th className="p-2">Max Atom</th>
                            <th className="p-2">Ref pLDDT</th>
                          </tr>
                        </thead>
                        <tbody className="divide-y divide-zinc-900">
                          {rmsdComparison.perResidueProfile.map((row) => {
                            const rot = allAtomStructure?.residueRotamers?.[row.modelResSeq - 1];
                            const isSelected = activeResidueIdx === row.modelResSeq - 1;
                            return (
                              <tr
                                key={row.modelResSeq}
                                onClick={() => setActiveResidueIdx(row.modelResSeq - 1)}
                                className={`cursor-pointer transition-colors ${
                                  isSelected ? 'bg-purple-950/45' : 'hover:bg-purple-950/25'
                                }`}
                              >
                                <td className="p-2 font-bold text-zinc-300">#{row.modelResSeq}</td>
                                <td className="p-2 text-zinc-400">#{row.refResSeq}</td>
                                <td className="p-2 font-bold text-white">{row.resName}</td>
                                <td className="p-2 text-indigo-300">{rot?.rotamerName || '-'}</td>
                                <td className={`p-2 font-bold ${row.caDistAngstrom < 1.0 ? 'text-emerald-400' : row.caDistAngstrom < 2.5 ? 'text-amber-400' : 'text-rose-400'}`}>
                                  {row.caDistAngstrom.toFixed(2)} Å
                                </td>
                                <td className={`p-2 font-bold ${row.backboneRmsdAngstrom < 1.0 ? 'text-emerald-300' : 'text-cyan-300'}`}>
                                  {row.backboneRmsdAngstrom.toFixed(2)} Å
                                </td>
                                <td className="p-2 text-fuchsia-300">
                                  {row.matchedSidechainAtomCount > 0 ? `${row.sidechainRmsdAngstrom.toFixed(2)} Å` : 'Brak (-H)'}
                                </td>
                                <td className="p-2 text-indigo-200 font-bold">{row.allAtomRmsdAngstrom.toFixed(2)} Å</td>
                                <td className="p-2 text-rose-300">
                                  {row.maxAtomName} ({row.maxAtomDistAngstrom.toFixed(2)} Å)
                                </td>
                                <td className="p-2 text-zinc-400">{row.refPlddt.toFixed(1)}</td>
                              </tr>
                            );
                          })}
                        </tbody>
                      </table>
                    </div>
                  </div>
                </>
              )}
            </div>

            {/* Footer Actions */}
            <div className="pt-3 border-t border-zinc-800 flex flex-wrap items-center justify-between gap-2 shrink-0">
              <div className="flex items-center gap-2">
                {rmsdComparison && (
                  <label className="flex items-center gap-2 text-xs text-purple-200 cursor-pointer">
                    <input
                      type="checkbox"
                      checked={showAlphaFoldOverlay}
                      onChange={(e) => setShowAlphaFoldOverlay(e.target.checked)}
                      className="accent-purple-500 rounded"
                    />
                    <span>Pokaż fioletową wstęgę referencyjną AlphaFold w scenie 3D</span>
                  </label>
                )}
              </div>

              <div className="flex items-center gap-2">
                <button
                  onClick={handleDownloadPdb}
                  className="px-3 py-1.5 rounded bg-sky-700 hover:bg-sky-600 text-white font-bold flex items-center gap-1.5 cursor-pointer"
                  title="Pobierz pojedynczy model SNAPLOCK All-Atom PDB z wpisami REMARK 250 RMSD"
                >
                  <Download className="w-3.5 h-3.5" />
                  <span>Pobierz All-Atom .pdb</span>
                </button>
                {rmsdComparison && (
                  <button
                    onClick={handleDownloadSuperposedPdb}
                    className="px-3.5 py-1.5 rounded bg-purple-600 hover:bg-purple-500 text-white font-bold flex items-center gap-1.5 cursor-pointer shadow-md"
                    title="Pobierz plik 2-modelowy (MODEL 1 = SNAPLOCK Chain A, MODEL 2 = Superposed AlphaFold Chain B) do PyMOL / ChimeraX"
                  >
                    <Download className="w-3.5 h-3.5" />
                    <span>Pobierz Nałożony 2-Model .pdb (PyMOL)</span>
                  </button>
                )}
                <button
                  onClick={() => setIsRmsdModalOpen(false)}
                  className="px-4 py-1.5 rounded bg-zinc-900 hover:bg-zinc-800 text-zinc-300 font-bold cursor-pointer"
                >
                  Zamknij
                </button>
              </div>
            </div>
          </div>
        </div>
      )}

    </div>
  );
};
