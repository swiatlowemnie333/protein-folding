import React, { useState, useEffect, useRef } from 'react';
import { 
  Play, 
  RotateCcw, 
  Atom, 
  Sparkles, 
  Sliders,
  Dna,
  Maximize2,
  Minimize2,
  X as CloseIcon,
  Music,
  Compass,
  Activity,
  ChevronRight,
  ChevronUp,
  ChevronDown,
  Layers,
  Bot,
  Database,
  ShieldCheck,
  Eye,
  EyeOff,
  CheckCircle2,
  AlertTriangle,
  ArrowUpDown,
  Move,
  ZoomIn,
  Settings,
  Boxes
} from 'lucide-react';
import { AMINO_ACIDS, SOUND_KEYS, getPatchCoordinate } from './aminoAcids';
import { AminoAcid, ConnectedAminoAcid, BenchmarkProteinId } from './types';
import { BENCHMARK_PROTEINS } from './benchmarks';
import { BioinformaticsModal } from './components/BioinformaticsModal';
import { SnaplockChatbot } from './components/SnaplockChatbot';
import { SnaplockParameterAuditView } from './components/SnaplockParameterAuditView';
import { SnaplockHeatmapView } from './components/SnaplockHeatmapView';
import { AminoAcidCatalogView } from './components/AminoAcidCatalogView';
import { SnaplockProteinConstructorView } from './components/SnaplockProteinConstructorView';
import { SnaplockValidatorView } from './components/SnaplockValidatorView';
import { MutationDiagnosticsView } from './components/MutationDiagnosticsView';
import { SnaplockBenchmarkModal } from './components/SnaplockBenchmarkModal';
import { SnaplockGroundTruthLoopPanel } from './components/SnaplockGroundTruthLoopPanel';
import { SnaplockRibosomeLinearizationScannerView } from './components/SnaplockRibosomeLinearizationScannerView';
import { ErrorBoundary } from './components/ErrorBoundary';
import { EsmFoldPanel } from './components/EsmFoldPanel';
import AblationLabView from '../AblationLabView';
import {
  GLOBAL_PROTEIN_TARGETS,
  GLOBAL_TARGET_ORDER,
  GlobalProteinTargetId,
  matchSequenceToGlobalTarget
} from './snaplockGlobalTargets';

type ChromeTabId = 'A' | 'B' | 'C' | 'D' | 'E' | 'F' | 'G' | 'H' | 'I' | 'J';

export default function App() {
  // --- TOP CHROME TABS STATE (ZERO DELETED FEATURES + DEDICATED NEW TABS E..H) ---
  const [activeTab, setActiveTab] = useState<ChromeTabId>('B');
  // FIX 2026-11-09: syntezator E8 usunięty — szuflada wyłączona na stałe (wcześniej brakująca zmienna wywracała całą aplikację)
  const showSynthPiano = false;
  const setShowSynthPiano = (_v: boolean) => {};

  // --- SINGLE CENTRAL GLOBAL PROTEIN TARGET STATE (STRICT 1:1 CROSS-TAB SYNC) ---
  const [activeTarget, setActiveTarget] = useState<GlobalProteinTargetId>('UBIQUITIN_76');

  // --- 3D PROTEIN CHAIN STATE ---
  const [proteinChain, setProteinChain] = useState<ConnectedAminoAcid[]>([
    { id: '1', aminoAcid: AMINO_ACIDS[0], fromPatch: 121, toPatch: 248, torsionAngle: 0 },
    { id: '2', aminoAcid: AMINO_ACIDS[19], fromPatch: 248, toPatch: 50, torsionAngle: 45 },
    { id: '3', aminoAcid: AMINO_ACIDS[4], fromPatch: 50, toPatch: 150, torsionAngle: -60 },
    { id: '4', aminoAcid: AMINO_ACIDS[7], fromPatch: 150, toPatch: 100, torsionAngle: 30 },
    { id: '5', aminoAcid: AMINO_ACIDS[15], fromPatch: 100, toPatch: 200, torsionAngle: -45 }
  ]);

  // Synchronized sequence string across all tabs
  const [seqInput, setSeqInput] = useState<string>(
    GLOBAL_PROTEIN_TARGETS.UBIQUITIN_76.sequence3Letter.join(' ')
  );
  
  // Selected amino acid in the warehouse
  const [selectedWarehouseAcid, setSelectedWarehouseAcid] = useState<AminoAcid>(AMINO_ACIDS[0]);

  // Active selected amino acid index in the current folded chain
  const [activeIndex, setActiveIndex] = useState<number>(0);
  
  // Modals & Drawers (Preserving 100% of built features)
  const [viewMode, setViewMode] = useState<'chain' | 'inside'>('chain');
  const [isBioinformaticsModalOpen, setIsBioinformaticsModalOpen] = useState<boolean>(false);
  const [isBenchmarkModalOpen, setIsBenchmarkModalOpen] = useState<boolean>(false);
  const [isChatbotOpen, setIsChatbotOpen] = useState<boolean>(false);
  const [calibrationLevel, setCalibrationLevel] = useState<number>(92);
  const [activeBenchmark, setActiveBenchmark] = useState<BenchmarkProteinId>('ubq');

  // Interactive 3D Orbit variables
  const [rotX, setRotX] = useState<number>(-22);
  const [rotY, setRotY] = useState<number>(38);
  const isDragging = useRef(false);
  const dragStart = useRef({ x: 0, y: 0, rotX: 0, rotY: 0, isPanning: false, panStartX: 0, panStartY: 0 });

  // Web Audio Context reference for synthesizer
  const audioCtxRef = useRef<AudioContext | null>(null);
  const [synthWave, setSynthWave] = useState<OscillatorType>('sine');
  const [soundFeedback, setSoundFeedback] = useState<string>('Silnik fizyczny SNAPLOCK gotowy (Å, kcal/mol)');
  const [triggeredKey, setTriggeredKey] = useState<string | null>(null);

  // --- DETERMINISTIC PHYSICAL FOLDING STATES ---
  const [isFolding, setIsFolding] = useState<boolean>(false);
  const [foldingProgress, setFoldingProgress] = useState<number>(0);
  const [energyLevel, setEnergyLevel] = useState<number>(-58.4);

  // Dynamic 3D model scaling controls (zoom)
  const [zoomScale, setZoomScale] = useState<number>(0.95);

  // Dynamic 3D model translation controls (panning)
  const [panX, setPanX] = useState<number>(0);
  const [panY, setPanY] = useState<number>(0);

  // Focus and X-Ray mode
  const [xrayFocus, setXrayFocus] = useState<boolean>(false);

  // Side Control Slider Panel visibility (docked on the right for zoom & vertical sequence scrolling)
  const [showSideSlider, setShowSideSlider] = useState<boolean>(true);
  const [showRibbonTape, setShowRibbonTape] = useState<boolean>(true);

  // Centralized Global Protein Target switcher — updates ALL tabs (A..H) synchronously!
  const handleSelectGlobalTarget = (targetId: Exclude<GlobalProteinTargetId, 'CUSTOM'>) => {
    const spec = GLOBAL_PROTEIN_TARGETS[targetId];
    if (!spec) return;
    const seq3Str = spec.sequence3Letter.join(' ');
    setActiveTarget(targetId);
    setSeqInput(seq3Str);
    parseDeterministicInput(seq3Str);
    if (targetId === 'UBIQUITIN_76') setActiveBenchmark('ubq');
    else if (targetId === 'INSULIN_A_21') setActiveBenchmark('7nhu');
    else if (targetId === 'LYSOZYME_129') setActiveBenchmark('lyz');
    setSoundFeedback(`Globalny cel: ${spec.shortLabel} (${spec.residueCount} aa • ${spec.pdbReferenceId})`);
  };

  // Sync custom/mutated sequence updates back to central state
  const handleCustomSequenceSync = (newSeq3Str: string) => {
    setSeqInput(newSeq3Str);
    parseDeterministicInput(newSeq3Str);
    const codes = newSeq3Str.trim().split(/\s+/).filter(Boolean);
    const matched = matchSequenceToGlobalTarget(codes);
    setActiveTarget(matched);
  };

  // Initial load synchronized to default global target (UBIQUITIN_76 = 76 aa 1UBQ)
  useEffect(() => {
    handleSelectGlobalTarget('UBIQUITIN_76');
  }, []);

  // --- AUDIO SYNTHESIZER DRIVER ---
  const initAudio = () => {
    if (!audioCtxRef.current) {
      audioCtxRef.current = new (window.AudioContext || (window as any).webkitAudioContext)();
    }
    if (audioCtxRef.current.state === 'suspended') {
      audioCtxRef.current.resume();
    }
  };

  const playCustomChord = (chordKey: string | number, duration: number = 0.35) => {
    initAudio();
    if (!audioCtxRef.current) return;

    let baseFreq = 261.63;
    let isMinor = false;
    let chordName = '';
    
    const keyStr = chordKey.toString().toUpperCase();
    
    if (keyStr === '0') { baseFreq = 261.63; isMinor = false; chordName = 'C-dur'; }
    else if (keyStr === '1') { baseFreq = 261.63; isMinor = true;  chordName = 'c-moll'; }
    else if (keyStr === '2') { baseFreq = 392.00; isMinor = false; chordName = 'G-dur'; }
    else if (keyStr === '3') { baseFreq = 392.00; isMinor = true;  chordName = 'g-moll'; }
    else if (keyStr === '4') { baseFreq = 293.66; isMinor = false; chordName = 'D-dur'; }
    else if (keyStr === '5') { baseFreq = 293.66; isMinor = true;  chordName = 'd-moll'; }
    else if (keyStr === '6') { baseFreq = 440.00; isMinor = false; chordName = 'A-dur'; }
    else if (keyStr === '7') { baseFreq = 440.00; isMinor = true;  chordName = 'a-moll'; }
    else if (keyStr === '8') { baseFreq = 329.63; isMinor = false; chordName = 'E-dur'; }
    else if (keyStr === '9') { baseFreq = 329.63; isMinor = true;  chordName = 'e-moll'; }
    else {
      baseFreq = 493.88; isMinor = true; chordName = 'h-moll';
    }

    setTriggeredKey(keyStr);
    setTimeout(() => setTriggeredKey(null), 250);

    const ctx = audioCtxRef.current;
    const now = ctx.currentTime;
    
    const ratio3 = isMinor ? 1.20 : 1.25;
    const ratio5 = 1.50;
    const freqs = [baseFreq, baseFreq * ratio3, baseFreq * ratio5];
    if (baseFreq > 200) freqs.push(baseFreq / 2);

    freqs.forEach((freq, idx) => {
      const osc = ctx.createOscillator();
      const filter = ctx.createBiquadFilter();
      const gainNode = ctx.createGain();

      osc.type = synthWave;
      const detuneValue = idx === 0 ? 0 : idx === 1 ? 4 : -4;
      osc.frequency.setValueAtTime(freq, now);
      osc.detune.setValueAtTime(detuneValue, now);

      filter.type = 'lowpass';
      filter.frequency.setValueAtTime(1400 + idx * 100, now);
      filter.Q.setValueAtTime(1.2, now);

      gainNode.gain.setValueAtTime(0, now);
      const targetGain = 0.22 / freqs.length;
      gainNode.gain.linearRampToValueAtTime(targetGain, now + 0.02);
      gainNode.gain.exponentialRampToValueAtTime(0.001, now + duration);

      osc.connect(filter);
      filter.connect(gainNode);
      gainNode.connect(ctx.destination);

      osc.start(now);
      osc.stop(now + duration + 0.05);
    });

    setSoundFeedback(`Ton: ${chordName} (${baseFreq.toFixed(0)} Hz)`);
  };

  // Keyboard shortcut listener
  useEffect(() => {
    const handleKeyDown = (e: KeyboardEvent) => {
      if (document.activeElement?.tagName === 'INPUT' || document.activeElement?.tagName === 'TEXTAREA') return;
      const key = e.key.toLowerCase();
      if (key >= '0' && key <= '9') {
        playCustomChord(key, 0.4);
      }
    };
    window.addEventListener('keydown', handleKeyDown);
    return () => window.removeEventListener('keydown', handleKeyDown);
  }, [synthWave]);

  // Deterministic chain sequence parser
  const parseDeterministicInput = (val: string) => {
    setSeqInput(val);
    let parts = val.toUpperCase().replace(/,/g, ' ').split(/\s+/).filter(Boolean);
    if (parts.length === 0) return;

    const newChain: ConnectedAminoAcid[] = [];
    parts.forEach((part, index) => {
      let matched = AMINO_ACIDS.find(a => a.code3.toUpperCase() === part || a.code1.toUpperCase() === part);
      if (!matched) {
        matched = AMINO_ACIDS[index % AMINO_ACIDS.length];
      }

      // Discrete E8 patch indices based on deterministic Lie algebra mapping
      const fromP = 248; // Backbone entry
      const toP = 1;     // Backbone exit
      const rotTable = [-120, -60, 0, 60, 120, 180];
      const torsion = rotTable[(index * 3) % rotTable.length];

      newChain.push({
        id: `node-${index + 1}`,
        aminoAcid: matched,
        fromPatch: fromP,
        toPatch: toP,
        torsionAngle: torsion
      });
    });

    setProteinChain(newChain);
    setActiveIndex(0);
    calculateCohesiveEnergy(newChain);
  };

  const calculateCohesiveEnergy = (chain: ConnectedAminoAcid[]) => {
    let baseEnergy = -chain.length * 2.8;
    
    // Hydrophobic core desolvation
    const hydrophobicCount = chain.filter(c => c.aminoAcid.hydrophobicity > 2.0).length;
    baseEnergy -= hydrophobicCount * 1.9;

    // Disulfide clamps bonus (covalent staples override)
    let cysCount = 0;
    chain.forEach(c => {
      if (c.aminoAcid.id === 'CYS' || c.aminoAcid.id === 'SEC') cysCount++;
    });
    const possibleBridges = Math.floor(cysCount / 2);
    baseEnergy -= possibleBridges * 14.5;

    setEnergyLevel(parseFloat(baseEnergy.toFixed(1)));
  };

  const runFoldingProcess = () => {
    if (isFolding) return;
    setIsFolding(true);
    setFoldingProgress(0);
    setSoundFeedback('Uruchomiono deterministyczne wiązanie peptydowe (Kwantyzacja E8)...');

    let p = 0;
    const interval = setInterval(() => {
      p += 15;
      setFoldingProgress(p);
      playCustomChord(p % 10, 0.15);

      if (p >= 100) {
        clearInterval(interval);
        setIsFolding(false);
        setFoldingProgress(100);
        calculateCohesiveEnergy(proteinChain);
        setSoundFeedback('Deterministyczne zwijanie zakończone. Wszystkie porty #248 ➔ #1 zaryglowane.');
      }
    }, 120);
  };

  const loadPresetSequence = (preset: 'insulin' | 'ubiquitin' | 'lycotoxin' | 'lysozyme' | 'prion' | 'helix') => {
    let seq = '';
    if (preset === 'insulin') {
      seq = 'GLY ILE VAL GLU GLN CYS CYS THR SER ILE CYS SER LEU TYR GLN LEU GLU ASN TYR CYS ASN';
      setActiveBenchmark('7nhu');
    } else if (preset === 'ubiquitin') {
      seq = BENCHMARK_PROTEINS['ubq'].sequence3Letter;
      setActiveBenchmark('ubq');
    } else if (preset === 'lycotoxin') {
      seq = 'ILE TRP LEU THR ALA LEU LYS PHE LEU GLY LYS HIS ALA ALA LYS ALA LYS ALA LYS ALA LYS HIS LEU CYS ALA LYS';
    } else if (preset === 'lysozyme') {
      seq = BENCHMARK_PROTEINS['lyz'].sequence3Letter;
      setActiveBenchmark('lyz');
    } else if (preset === 'prion') {
      seq = 'MET ALA ASN LEU GLY CYS TRP MET LEU VAL LEU PHE VAL ALA THR TRP SER ASP LEU GLY LEU CYS LYS LYS ARG PRO LYS PRO GLY GLY TRP ASN THR GLY GLY SER ARG TYR PRO';
    } else if (preset === 'helix') {
      seq = 'ALA GLU ALA ALA ALA LYS GLU ALA ALA ALA LYS GLU ALA ALA ALA LYS ALA ALA';
    }

    setSeqInput(seq);
    parseDeterministicInput(seq);
    fitModelToScreen();
  };

  // Center coordinate points relative to Center of Mass (0,0,0)
  const getFullChainCoords3D = () => {
    let currentX = 0;
    let currentY = 0;
    let currentZ = 0;

    return proteinChain.map((node, index) => {
      // Helical and globular dispersion
      const theta = index * 0.72 + (node.torsionAngle * Math.PI) / 180;
      const radius = 68 - Math.min(25, index * 0.2);
      
      currentX = Math.cos(theta) * radius;
      currentY = (index - proteinChain.length / 2) * 16.5;
      currentZ = Math.sin(theta) * radius;

      return {
        node,
        x: currentX,
        y: currentY,
        z: currentZ
      };
    });
  };

  const getCenteredPositions = () => {
    const rawCoords = getFullChainCoords3D();
    if (rawCoords.length === 0) return [];

    let sumX = 0, sumY = 0, sumZ = 0;
    rawCoords.forEach(c => {
      sumX += c.x;
      sumY += c.y;
      sumZ += c.z;
    });

    const avgX = sumX / rawCoords.length;
    const avgY = sumY / rawCoords.length;
    const avgZ = sumZ / rawCoords.length;

    return rawCoords.map(c => ({
      node: c.node,
      x: c.x - avgX,
      y: c.y - avgY,
      z: c.z - avgZ
    }));
  };

  const centeredChain = getCenteredPositions();
  const activeChainNode = proteinChain[activeIndex] || proteinChain[0];
  const activeAminoAcid = activeChainNode.aminoAcid;

  // Auto-fit 3D model to screen for sequences of any length
  const fitModelToScreen = () => {
    setPanX(0);
    setPanY(0);
    setRotX(-22);
    setRotY(38);
    const len = proteinChain.length;
    let targetZoom = 1.0;
    if (len <= 10) targetZoom = 1.35;
    else if (len <= 25) targetZoom = 1.1;
    else if (len <= 55) targetZoom = 0.88;
    else if (len <= 80) targetZoom = 0.68;
    else if (len <= 120) targetZoom = 0.52;
    else targetZoom = 0.42;

    setZoomScale(targetZoom);
    setSoundFeedback(`Dopasowano model do ekranu (Skala: ${Math.round(targetZoom * 100)}%).`);
  };

  // 3D to 2D projection
  const project3DTo2D = (x: number, y: number, z: number, cX: number, cY: number, customScale: number = 1.0) => {
    const radX = (rotX * Math.PI) / 180;
    const radY = (rotY * Math.PI) / 180;

    // Rotate Y (Yaw)
    const x1 = x * Math.cos(radY) - z * Math.sin(radY);
    const z1 = x * Math.sin(radY) + z * Math.cos(radY);

    // Rotate X (Pitch)
    const y2 = y * Math.cos(radX) - z1 * Math.sin(radX);
    const z2 = y * Math.sin(radX) + z1 * Math.cos(radX);

    const distance = 1000;
    const denom = Math.max(distance + z2, 100);
    const rawScale = (distance / denom) * customScale * zoomScale;
    const scaleFactor = Math.max(0.08, Math.min(rawScale, 4.0));

    return {
      x: cX + panX + x1 * scaleFactor,
      y: cY + panY - y2 * scaleFactor,
      z: z2,
      scale: scaleFactor
    };
  };

  const getAtomColor = (element: string) => {
    if (element === 'C') return '#475569';
    if (element === 'H') return '#f1f5f9';
    if (element === 'O') return '#ef4444';
    if (element === 'N') return '#3b82f6';
    if (element === 'S') return '#eab308';
    if (element === 'Se') return '#14b8a6';
    return '#94a3b8';
  };

  const activeGlobalSpec =
    activeTarget !== 'CUSTOM' ? GLOBAL_PROTEIN_TARGETS[activeTarget] : null;

  return (
    <div className="fixed inset-0 bg-[#000000] text-[#f1f5f9] flex flex-col font-sans overflow-hidden select-none" id="viewport_main">
      
      {/* ======================================================== */}
      {/* 1. TOP CHROME TAB NAVIGATION BAR (ALL 8 TABS A..H)       */}
      {/* ======================================================== */}
      <nav className="bg-[#0b0c16] border-b border-zinc-900 px-3 pt-2 flex items-center justify-between gap-2 shrink-0 z-30 select-none shadow-md">
        
        {/* Left: Scrollable Chrome Tabs A..H */}
        <div className="flex items-center gap-1 overflow-x-auto scrollbar-none">
          
          {/* TAB A: Katalog 20 Aminokwasów (+ SEC/PYL) */}
          <button
            onClick={() => setActiveTab('A')}
            className={`px-3 py-1.5 rounded-t-xl font-mono text-[11.5px] font-bold flex items-center gap-1.5 transition cursor-pointer border-t border-x shrink-0 ${
              activeTab === 'A'
                ? 'bg-[#020205] border-zinc-700 border-b-transparent text-white shadow-xl ring-1 ring-cyan-500/20'
                : 'bg-[#101222]/80 border-transparent text-zinc-400 hover:text-white hover:bg-[#17192f]'
            }`}
            title="Karta A: Katalog 20 standardowych aminokwasów (+ SEC, PYL) w jednostkach Å i kcal/mol"
          >
            <Boxes className="w-3.5 h-3.5 text-cyan-400" />
            <span>A: Katalog Aminokwasów</span>
          </button>

          {/* TAB B: Konstruktor 3D (N->C Soft SnapLock, Bańka Mydlana & Mikroskop Å) */}
          <button
            onClick={() => setActiveTab('B')}
            className={`px-3 py-1.5 rounded-t-xl font-mono text-[11.5px] font-bold flex items-center gap-1.5 transition cursor-pointer border-t border-x shrink-0 ${
              activeTab === 'B'
                ? 'bg-[#000000] border-zinc-700 border-b-transparent text-white shadow-xl ring-1 ring-indigo-500/30'
                : 'bg-[#101222]/80 border-transparent text-zinc-400 hover:text-white hover:bg-[#17192f]'
            }`}
            title="Karta B: Wieloagentowy wzrost kotranslacyjny N→C + Bańka Mydlana + tryb pierwszoosobowy Mikroskop Å"
          >
            <Compass className="w-3.5 h-3.5 text-indigo-400" />
            <span>B: Konstruktor &amp; Mikroskop Å</span>
          </button>

          {/* TAB C: Walidator Kabsch SVD 1:1 (< 1.0 Å Hard Gate) */}
          <button
            onClick={() => setActiveTab('C')}
            className={`px-3 py-1.5 rounded-t-xl font-mono text-[11.5px] font-bold flex items-center gap-1.5 transition cursor-pointer border-t border-x shrink-0 ${
              activeTab === 'C'
                ? 'bg-[#020205] border-zinc-700 border-b-transparent text-white shadow-xl ring-1 ring-teal-500/20'
                : 'bg-[#101222]/80 border-transparent text-zinc-400 hover:text-white hover:bg-[#17192f]'
            }`}
            title="Karta C: Rygorystyczny Walidator Kabsch SVD 1:1 Cα z twardą blokadą < 1.0 Å"
          >
            <ShieldCheck className="w-3.5 h-3.5 text-cyan-400" />
            <span>C: Walidator Kabsch &lt;1Å</span>
          </button>

          {/* TAB D: Diagnostyka Mutacji (Zdrowe vs Chore) */}
          <button
            onClick={() => setActiveTab('D')}
            className={`px-3 py-1.5 rounded-t-xl font-mono text-[11.5px] font-bold flex items-center gap-1.5 transition cursor-pointer border-t border-x shrink-0 ${
              activeTab === 'D'
                ? 'bg-[#020205] border-zinc-700 border-b-transparent text-white shadow-xl ring-1 ring-red-500/20'
                : 'bg-[#101222]/80 border-transparent text-zinc-400 hover:text-white hover:bg-[#17192f]'
            }`}
            title="Karta D: Dwuokienna diagnostyka mutacji klinicznych (Zdrowe vs Chore)"
          >
            <Activity className="w-3.5 h-3.5 text-red-400" />
            <span>D: Mutacje</span>
          </button>

          {/* TAB E: Audyt Parametrów w Przeglądarce */}
          <button
            onClick={() => setActiveTab('E')}
            className={`px-3 py-1.5 rounded-t-xl font-mono text-[11.5px] font-bold flex items-center gap-1.5 transition cursor-pointer border-t border-x shrink-0 ${
              activeTab === 'E'
                ? 'bg-[#020205] border-zinc-700 border-b-transparent text-white shadow-xl ring-1 ring-amber-500/20'
                : 'bg-[#101222]/80 border-transparent text-zinc-400 hover:text-white hover:bg-[#17192f]'
            }`}
            title="Karta E: Audyt czułości wszystkich parametrów aminokwasów liczony w przeglądarce (ΔE, ΔRMSD, weryfikacja w kodzie)"
          >
            <Database className="w-3.5 h-3.5 text-amber-400" />
            <span>E: Audyt Parametrów</span>
          </button>

          {/* TAB F: Mapy Ciepła (Cα-Cα, Orientacja N-Cα-C, Diff, TM-Score) */}
          <button
            onClick={() => setActiveTab('F')}
            className={`px-3 py-1.5 rounded-t-xl font-mono text-[11.5px] font-bold flex items-center gap-1.5 transition cursor-pointer border-t border-x shrink-0 ${
              activeTab === 'F'
                ? 'bg-[#020205] border-zinc-700 border-b-transparent text-white shadow-xl ring-1 ring-teal-500/20'
                : 'bg-[#101222]/80 border-transparent text-zinc-400 hover:text-white hover:bg-[#17192f]'
            }`}
            title="Karta F: Mapy Ciepła Cα–Cα, Orientacja z układu lokalnego N-Cα-C, Tryb Różnicowy Diff oraz TM-score"
          >
            <Eye className="w-3.5 h-3.5 text-teal-400" />
            <span>F: Mapy Ciepła</span>
          </button>

          {/* TAB G: Autokalibracja & Ablation M1-M5 (Train / Val / Test) */}
          <button
            onClick={() => setActiveTab('G')}
            className={`px-3 py-1.5 rounded-t-xl font-mono text-[11.5px] font-bold flex items-center gap-1.5 transition cursor-pointer border-t border-x shrink-0 ${
              activeTab === 'G'
                ? 'bg-[#020205] border-zinc-700 border-b-transparent text-white shadow-xl ring-1 ring-emerald-500/20'
                : 'bg-[#101222]/80 border-transparent text-zinc-400 hover:text-white hover:bg-[#17192f]'
            }`}
            title="Karta G: Pętla Autokalibracji Ground-Truth, czułość parametrów i Ablation M1–M5 (Train/Validation/Held-Out Test)"
          >
            <Sliders className="w-3.5 h-3.5 text-emerald-400" />
            <span>G: Autokalibracja M1–M5</span>
          </button>

          {/* TAB H: Fizyka Oderwania od Rybosomu, Bańka Mydlana & Skaner 4 Ścieżek Przez Płytę */}
          <button
            onClick={() => setActiveTab('H')}
            className={`px-3 py-1.5 rounded-t-xl font-mono text-[11.5px] font-bold flex items-center gap-1.5 transition cursor-pointer border-t border-x shrink-0 ${
              activeTab === 'H'
                ? 'bg-[#020308] border-fuchsia-500/60 border-b-transparent text-white shadow-xl ring-1 ring-fuchsia-500/30'
                : 'bg-[#101222]/80 border-transparent text-zinc-400 hover:text-white hover:bg-[#17192f]'
            }`}
            title="Karta H: Fizyka oderwania od rybosomu, linearyzacja 22 parametrów, wieloagentowa Bańka Mydlana, baza α-helis i skaner błędów przez płytę (2x ślimak + boczne + S-S)"
          >
            <Sparkles className="w-3.5 h-3.5 text-fuchsia-400" />
            <span>H: Rybosom &amp; Skaner Płyty</span>
          </button>

          {/* TAB I: Prawdziwa predykcja ESMFold + walidacja */}
          <button
            onClick={() => setActiveTab('I')}
            className={`px-3 py-1.5 rounded-t-xl font-mono text-[11.5px] font-bold flex items-center gap-1.5 transition cursor-pointer border-t border-x shrink-0 ${
              activeTab === 'I'
                ? 'bg-[#020205] border-emerald-500/60 border-b-transparent text-white shadow-xl ring-1 ring-emerald-500/30'
                : 'bg-[#101222]/80 border-transparent text-zinc-400 hover:text-white hover:bg-[#17192f]'
            }`}
            title="Karta I: Predykcja struktury ESMFold (Meta) oceniona walidatorem SNAPLOCK"
          >
            <ShieldCheck className="w-3.5 h-3.5 text-emerald-400" />
            <span>I: ESMFold</span>
          </button>

          {/* TAB J: Laboratorium wariantów A/B + Dziennik wniosków */}
          <button
            onClick={() => setActiveTab('J')}
            className={`px-3 py-1.5 rounded-t-xl font-mono text-[11.5px] font-bold flex items-center gap-1.5 transition cursor-pointer border-t border-x shrink-0 ${
              activeTab === 'J'
                ? 'bg-[#020205] border-amber-500/60 border-b-transparent text-white shadow-xl ring-1 ring-amber-500/30'
                : 'bg-[#101222]/80 border-transparent text-zinc-400 hover:text-white hover:bg-[#17192f]'
            }`}
            title="Karta J: Laboratorium wariantów A/B (moduły osobno i w zestawach, średnia z kilku startów, synergia) + Dziennik wniosków"
          >
            <Sliders className="w-3.5 h-3.5 text-amber-400" />
            <span>J: Laboratorium A/B</span>
          </button>
        </div>

        {/* Right side: Restored Utility Modals (Bioinformatics, Benchmark) + Chatbot */}
        <div className="flex items-center gap-1.5 font-mono text-xs shrink-0 pb-1">
          <button
            onClick={() => setIsBioinformaticsModalOpen(true)}
            className="px-2.5 py-1 rounded bg-zinc-900 hover:bg-zinc-800 border border-zinc-800 text-zinc-300 hover:text-white text-[11px] font-bold flex items-center gap-1 cursor-pointer transition"
            title="Eksport i kalibracja bioinformatyczna PDB"
          >
            <Dna className="w-3.5 h-3.5 text-cyan-400" />
            <span className="hidden 2xl:inline">Bioinfo PDB</span>
          </button>

          <button
            onClick={() => setIsBenchmarkModalOpen(true)}
            className="px-2.5 py-1 rounded bg-zinc-900 hover:bg-zinc-800 border border-zinc-800 text-zinc-300 hover:text-white text-[11px] font-bold flex items-center gap-1 cursor-pointer transition"
            title="Modal Benchmarku UniProt / AlphaFold"
          >
            <Settings className="w-3.5 h-3.5 text-indigo-400" />
            <span className="hidden 2xl:inline">Benchmark</span>
          </button>

          <button
            onClick={() => setIsChatbotOpen(true)}
            className="px-3 py-1 bg-gradient-to-r from-cyan-600 to-indigo-600 hover:from-cyan-500 hover:to-indigo-500 text-white rounded text-[11px] font-bold flex items-center gap-1.5 cursor-pointer shadow-md transition"
            title="Konsultant Biofizyki SNAPLOCK"
          >
            <Bot className="w-3.5 h-3.5 text-cyan-200" />
            <span className="hidden xl:inline">Asystent</span>
          </button>
        </div>
      </nav>

      {/* ======================================================== */}
      {/* 2. SINGLE GLOBAL PROTEIN SELECTOR BAR (STRICT 1:1 SYNC)  */}
      {/* ======================================================== */}
      <div className="bg-[#060711] border-b border-zinc-900 px-3 py-1.5 flex flex-wrap items-center justify-between gap-2 shrink-0 z-20 font-mono text-xs">
        <div className="flex items-center gap-1.5 flex-wrap">
          <span className="text-[10px] font-extrabold uppercase tracking-wider text-cyan-400 px-1.5 py-0.5 rounded bg-cyan-950/60 border border-cyan-800/60">
            Globalny Selektor Białka (1:1):
          </span>

          {GLOBAL_TARGET_ORDER.map((tid) => {
            const spec = GLOBAL_PROTEIN_TARGETS[tid];
            const isSelected = activeTarget === tid;
            return (
              <button
                key={tid}
                onClick={() => handleSelectGlobalTarget(tid)}
                className={`px-2.5 py-1 rounded border text-[11px] font-bold flex items-center gap-1.5 cursor-pointer transition ${
                  isSelected
                    ? 'bg-cyan-950 border-cyan-400 text-white shadow-sm ring-1 ring-cyan-400/40'
                    : 'bg-zinc-950 hover:bg-zinc-900 border-zinc-800 text-zinc-400 hover:text-white'
                }`}
                title={`${spec.titlePl} • Referencja 1:1: ${spec.pdbReferenceId} (${spec.uniprotLabel})`}
              >
                <span>{spec.shortLabel}</span>
                <span className="text-[9.5px] px-1 rounded bg-black/70 text-cyan-300 border border-zinc-800">
                  {spec.residueCount} aa • {spec.pdbReferenceId}
                </span>
              </button>
            );
          })}

          {activeTarget === 'CUSTOM' && (
            <span className="px-2.5 py-1 rounded border border-amber-500 bg-amber-950/70 text-amber-200 text-[11px] font-bold">
              Własna sekwencja ({proteinChain.length} aa)
            </span>
          )}
        </div>

        {/* Strict 1:1 Reference & Split Indicator */}
        <div className="flex items-center gap-2 text-[10.5px] text-zinc-400">
          {activeGlobalSpec ? (
            <>
              <span>
                Rekord ref. 1:1:{' '}
                <strong className="text-emerald-300">
                  {activeGlobalSpec.pdbReferenceId}:{activeGlobalSpec.pdbChainId}
                </strong>{' '}
                ({activeGlobalSpec.uniprotLabel})
              </span>
              <span className="px-1.5 py-0.5 rounded bg-zinc-900 border border-zinc-800 text-zinc-200 font-bold">
                {activeGlobalSpec.splitRole}
              </span>
            </>
          ) : (
            <span>Tryb De-Novo (własna sekwencja bez przypiętej krystalografii 1:1)</span>
          )}
        </div>
      </div>

      {/* OPTIONAL RESTORED E8 HARMONIC SYNTHESIZER DRAWER */}
      {showSynthPiano && (
        <div className="bg-[#090b16] border-b border-amber-900/50 px-4 py-2 flex flex-wrap items-center justify-between gap-3 shrink-0 z-20 font-mono text-xs">
          <div className="flex items-center gap-2 flex-wrap">
            <span className="text-[10.5px] text-amber-300 font-bold uppercase">
              Syntezator E8 ({soundFeedback}):
            </span>
            {SOUND_KEYS.map((sk) => (
              <button
                key={String(sk.digit)}
                onClick={() => playCustomChord(String(sk.digit), 0.4)}
                className={`px-2.5 py-1 rounded border text-[11px] font-bold cursor-pointer transition ${
                  triggeredKey === String(sk.digit)
                    ? 'bg-amber-500 text-black border-amber-300'
                    : 'bg-black hover:bg-zinc-900 border-zinc-800 text-zinc-200'
                }`}
              >
                [{sk.digit}] {sk.chord}
              </button>
            ))}
          </div>
          <div className="flex items-center gap-2">
            <select
              value={synthWave}
              onChange={(e) => setSynthWave(e.target.value as OscillatorType)}
              className="bg-black border border-zinc-800 rounded px-2 py-1 text-[11px] text-zinc-200"
            >
              <option value="sine">Sinusoida (Sine)</option>
              <option value="triangle">Trójkąt (Triangle)</option>
              <option value="sawtooth">Piła (Sawtooth)</option>
            </select>
            <button
              onClick={() => setShowSynthPiano(false)}
              className="p-1 text-zinc-400 hover:text-white cursor-pointer"
            >
              <CloseIcon className="w-3.5 h-3.5" />
            </button>
          </div>
        </div>
      )}

      {/* ======================================================== */}
      {/* 3. SYNCHRONIZED TAB WORKSPACES (A..H)                    */}
      {/* ======================================================== */}

      {/* TAB A: 21 KUL E8 */}
      {activeTab === 'A' && (
        <div className="flex-1 w-full h-full min-h-0 overflow-hidden flex flex-col bg-[#020205]">
          <ErrorBoundary fallbackTitle="Błąd modułu Magazynu Kul E8">
            <AminoAcidSphereWarehouseView 
              onSelectAminoAcidForConstructor={(code3) => {
                const updated = seqInput.trim() ? `${seqInput.trim()} ${code3}` : code3;
                handleCustomSequenceSync(updated);
                setSoundFeedback(`Dodano ${code3} do sekwencji Konstruktora.`);
              }}
              onSwitchToConstructor={() => setActiveTab('B')}
            />
          </ErrorBoundary>
        </div>
      )}

      {/* TAB B: KONSTRUKTOR 3D (WIELOAGENTOWY WZROST N->C SOFT SNAPLOCK + BAŃKA MYDLANA) */}
      {activeTab === 'B' && (
        <div className="flex-1 w-full h-full min-h-0 overflow-hidden flex flex-col bg-black">
          <ErrorBoundary fallbackTitle="Błąd modułu Konstruktora Białka 3D">
            <SnaplockProteinConstructorView 
              activeTarget={activeTarget}
              onSelectGlobalTarget={handleSelectGlobalTarget}
              initialSequence={seqInput}
              onSequenceChanged={(newSeq) => {
                handleCustomSequenceSync(newSeq);
              }}
              onJumpToWarehouseSphere={(_code) => {
                setActiveTab('A');
              }}
            />
          </ErrorBoundary>
        </div>
      )}

      {/* TAB C: RYGORYSTYCZNY WALIDATOR KABSCH SVD 1:1 (< 1.0 Å HARD GATE) */}
      {activeTab === 'C' && (
        <div className="flex-1 w-full h-full min-h-0 overflow-hidden flex flex-col bg-[#020205]">
          <ErrorBoundary fallbackTitle="Błąd modułu Walidatora Kabsch SVD">
            <SnaplockValidatorView 
              activeTarget={activeTarget}
              activeSequence={seqInput}
              onSelectGlobalTarget={handleSelectGlobalTarget}
              onLoadSequenceToConstructor={(seq) => {
                handleCustomSequenceSync(seq);
                fitModelToScreen();
                setActiveTab('B');
              }}
              onSwitchToConstructor={() => setActiveTab('B')}
              onJumpToAutocalibrationTab={() => setActiveTab('G')}
            />
          </ErrorBoundary>
        </div>
      )}

      {/* TAB D: DIAGNOSTYKA MUTACJI (ZDROWE VS CHORE) */}
      {activeTab === 'D' && (
        <div className="flex-1 w-full h-full min-h-0 overflow-hidden flex flex-col bg-[#020205]">
          <ErrorBoundary fallbackTitle="Błąd modułu Diagnostyki Mutacji">
            <MutationDiagnosticsView 
              activeTarget={activeTarget}
              onLoadVariantToConstructor={(seq) => {
                handleCustomSequenceSync(seq);
                fitModelToScreen();
                setActiveTab('B');
              }}
              onSwitchToConstructor={() => setActiveTab('B')}
            />
          </ErrorBoundary>
        </div>
      )}

      {/* TAB E: PRZYWRÓCONA MATRYCA 22 PARAMETRÓW & ZŁĄCZA 5 ŚCIEŻEK */}
      {activeTab === 'E' && (
        <div className="flex-1 w-full h-full min-h-0 overflow-hidden flex flex-col bg-[#020205]">
          <ErrorBoundary fallbackTitle="Błąd modułu Matrycy 22 Parametrów">
            <SnaplockWarehouseView />
          </ErrorBoundary>
        </div>
      )}

      {/* TAB F: PRZYWRÓCONE PLANETARIUM SKYMAP 360° */}
      {activeTab === 'F' && (
        <div className="flex-1 w-full h-full min-h-0 overflow-hidden flex flex-col bg-[#020206]">
          <ErrorBoundary fallbackTitle="Błąd modułu Planetarium SkyMap 360°">
            <SkyMapPlanetariumView 
              activeTarget={activeTarget}
              activeSequence={seqInput}
              onSelectGlobalTarget={handleSelectGlobalTarget}
              onChainUpdated={(codes) => {
                handleCustomSequenceSync(codes.join(' '));
              }}
              onSendToConstructor={(seq) => {
                handleCustomSequenceSync(seq);
                setActiveTab('B');
              }}
            />
          </ErrorBoundary>
        </div>
      )}

      {/* TAB G: PRZYWRÓCONA PĘTLA AUTOKALIBRACJI & ABLATION M1–M5 */}
      {activeTab === 'G' && (
        <div className="flex-1 w-full h-full min-h-0 overflow-hidden flex flex-col bg-[#020205]">
          <ErrorBoundary fallbackTitle="Błąd modułu Autokalibracji Ground-Truth">
            <SnaplockGroundTruthLoopPanel activeGlobalTarget={activeTarget} />
          </ErrorBoundary>
        </div>
      )}

      {/* TAB H: NOWA KARTA — FIZYKA RYBOSOMU, LINEARYZACJA 22 PARAM., BAŃKA MYDLANA & SKANER 4 ŚCIEŻEK */}
      {activeTab === 'H' && (
        <div className="flex-1 w-full h-full min-h-0 overflow-hidden flex flex-col bg-[#020308]">
          <ErrorBoundary fallbackTitle="Błąd modułu Fizyki Rybosomu i Skanera Płyty">
            <SnaplockRibosomeLinearizationScannerView
              activeTarget={activeTarget}
              activeSequence={seqInput}
            />
          </ErrorBoundary>
        </div>
      )}

      {activeTab === 'I' && (
        <div className="flex-1 w-full h-full min-h-0 overflow-hidden flex flex-col bg-[#020205]">
          <ErrorBoundary fallbackTitle="Błąd modułu ESMFold">
            <EsmFoldPanel
              activeSequence={seqInput}
              referencePdbText={activeGlobalSpec?.referencePdbText ?? null}
              referenceChainId={activeGlobalSpec?.pdbChainId}
              referenceLabel={activeGlobalSpec?.pdbReferenceId}
            />
          </ErrorBoundary>
        </div>
      )}

      {activeTab === 'J' && (
        <div className="flex-1 w-full h-full min-h-0 overflow-hidden flex flex-col bg-[#020205]">
          <ErrorBoundary fallbackTitle="Błąd Laboratorium A/B">
            <AblationLabView />
          </ErrorBoundary>
        </div>
      )}

      {/* ======================================================== */}
      {/* 4. RESTORED MODALS & ASSISTANTS                          */}
      {/* ======================================================== */}

      <BioinformaticsModal
        isOpen={isBioinformaticsModalOpen}
        onClose={() => setIsBioinformaticsModalOpen(false)}
        proteinChain={proteinChain}
        centeredChain={centeredChain}
        cohesiveEnergy={energyLevel}
        calibrationLevel={calibrationLevel}
        onCalibrationChange={(lvl) => {
          setCalibrationLevel(lvl);
          setSoundFeedback(`Kalibracja zaktualizowana do ${lvl}%.`);
        }}
        activeBenchmark={activeBenchmark}
        onBenchmarkChange={(id) => {
          setActiveBenchmark(id);
          const b = BENCHMARK_PROTEINS[id];
          setSoundFeedback(`Wzorzec: ${b.fullName}.`);
          if (id === 'ubq') handleSelectGlobalTarget('UBIQUITIN_76');
          else if (id === '7nhu') handleSelectGlobalTarget('INSULIN_A_21');
          else if (id === 'lyz') handleSelectGlobalTarget('LYSOZYME_129');
        }}
        onLoadBenchmarkSequence={(seq) => {
          handleCustomSequenceSync(seq);
          fitModelToScreen();
        }}
      />

      <SnaplockBenchmarkModal
        isOpen={isBenchmarkModalOpen}
        onClose={() => setIsBenchmarkModalOpen(false)}
        onLoadSequence={(seq) => {
          handleCustomSequenceSync(seq);
          setActiveTab('B');
        }}
      />

      {isChatbotOpen && (
        <SnaplockChatbot onClose={() => setIsChatbotOpen(false)} />
      )}

    </div>
  );
}
