import React, { useState, useMemo } from 'react';
import { 
  Grid, 
  Layers, 
  Eye, 
  ArrowRight, 
  Activity, 
  Sliders, 
  Sparkles,
  Info,
  Maximize2,
  RefreshCw,
  GitCompare,
  TrendingDown
} from 'lucide-react';
import { 
  GLOBAL_PROTEIN_TARGETS, 
  GlobalProteinTargetId 
} from '../snaplockGlobalTargets';
import { 
  runHierarchicalFoldingEngine 
} from '../snaplockPhysicalFoldingEngine';
import { 
  parsePdbForGroundTruth, 
  validateAgainstGroundTruth,
  Vec3Tuple,
  THREE_TO_ONE_MAP
} from '../ground_truth_validation';

interface SnaplockHeatmapViewProps {
  activeTarget: GlobalProteinTargetId;
  activeSequence: string;
  onSelectGlobalTarget?: (targetId: Exclude<GlobalProteinTargetId, 'CUSTOM'>) => void;
  onLoadSequenceToConstructor?: (seq: string) => void;
  onSwitchToConstructor?: () => void;
}

export type HeatmapMode = 'DISTANCE_CA' | 'ORIENTATION_FRAME' | 'DIFF_COMPARISON';

interface LocalFrame3D {
  ux: [number, number, number]; // along CA -> C
  uy: [number, number, number]; // in plane
  uz: [number, number, number]; // normal
  ca: [number, number, number];
}

function normalize(v: [number, number, number]): [number, number, number] {
  const d = Math.hypot(v[0], v[1], v[2]);
  if (d < 1e-9) return [1, 0, 0];
  return [v[0] / d, v[1] / d, v[2] / d];
}

function cross(a: [number, number, number], b: [number, number, number]): [number, number, number] {
  return [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0]
  ];
}

function dot(a: [number, number, number], b: [number, number, number]): number {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

export const SnaplockHeatmapView: React.FC<SnaplockHeatmapViewProps> = ({
  activeTarget,
  activeSequence,
  onSelectGlobalTarget,
  onLoadSequenceToConstructor,
  onSwitchToConstructor
}) => {
  const [heatmapMode, setHeatmapMode] = useState<HeatmapMode>('DISTANCE_CA');
  const [compareTarget, setCompareTarget] = useState<'PREDICTED_VS_PDB' | 'HEALTHY_VS_MUTANT'>('PREDICTED_VS_PDB');
  const [mutantSite, setMutantSite] = useState<number>(48); // K48R in Ubiquitin
  const [mutantNewAA, setMutantNewAA] = useState<string>('ARG');
  const [distanceCutoff, setDistanceCutoff] = useState<number>(20.0); // max Å scale
  const [hoveredCell, setHoveredCell] = useState<{ i: number; j: number; val1: number; val2?: number; diff?: number } | null>(null);

  const activeSpec = activeTarget !== 'CUSTOM' ? GLOBAL_PROTEIN_TARGETS[activeTarget] : null;
  const codes = useMemo(() => activeSequence.trim().split(/\s+/).filter(Boolean), [activeSequence]);
  const N = codes.length;

  // 1. Fold predicted structure in memory
  const predictedModel = useMemo(() => {
    if (N === 0) return null;
    const res = runHierarchicalFoldingEngine(codes, 45, 15);
    return res;
  }, [codes, N]);

  // 2. Reference PDB parse (Ground Truth)
  const refPdbData = useMemo(() => {
    if (!activeSpec || !activeSpec.pdbReferenceContent) return null;
    try {
      return parsePdbForGroundTruth(activeSpec.pdbReferenceContent, {
        targetChainId: activeSpec.targetChainId,
        enforceSingleChain: true,
        stripHydrogens: true
      });
    } catch {
      return null;
    }
  }, [activeSpec]);

  // 3. Ground-truth validation (Kabsch superposition + TM-Score)
  const validationReport = useMemo(() => {
    if (!predictedModel || !refPdbData) return null;
    try {
      const predCa = predictedModel.caPositions as Vec3Tuple[];
      return validateAgainstGroundTruth(predCa, refPdbData, {
        allowPartialAlignment: true,
        targetCaRmsdAngstrom: 1.0,
        iterationId: `heatmap_${activeTarget}`
      });
    } catch {
      return null;
    }
  }, [predictedModel, refPdbData, activeTarget]);

  // 4. Mutated variant folding (Healthy vs Mutant)
  const mutantCodes = useMemo(() => {
    const copy = [...codes];
    const safeSite = Math.max(0, Math.min(mutantSite - 1, copy.length - 1));
    if (copy[safeSite]) {
      copy[safeSite] = mutantNewAA;
    }
    return copy;
  }, [codes, mutantSite, mutantNewAA]);

  const mutantModel = useMemo(() => {
    if (compareTarget !== 'HEALTHY_VS_MUTANT' || mutantCodes.length === 0) return null;
    return runHierarchicalFoldingEngine(mutantCodes, 45, 15);
  }, [compareTarget, mutantCodes]);

  // 5. Compute Local Frames (N, CA, C) for orientation heatmap
  // ux = CA -> C, plane defined by N - CA - C, uz = ux x (N - CA), uy = uz x ux
  const localFrames = useMemo(() => {
    if (!predictedModel || !predictedModel.allAtomStructure) return [];
    const atoms = predictedModel.allAtomStructure.atoms;
    const frames: (LocalFrame3D | null)[] = [];

    for (let i = 1; i <= N; i++) {
      const caAtom = atoms.find(a => a.resSeq === i && a.atomName === 'CA');
      const cAtom = atoms.find(a => a.resSeq === i && a.atomName === 'C');
      const nAtom = atoms.find(a => a.resSeq === i && a.atomName === 'N');

      if (caAtom && cAtom && nAtom) {
        const ca: [number, number, number] = [caAtom.x, caAtom.y, caAtom.z];
        const c: [number, number, number] = [cAtom.x, cAtom.y, cAtom.z];
        const n: [number, number, number] = [nAtom.x, nAtom.y, nAtom.z];

        const ux = normalize([c[0] - ca[0], c[1] - ca[1], c[2] - ca[2]]);
        const vN = [n[0] - ca[0], n[1] - ca[1], n[2] - ca[2]] as [number, number, number];
        const uz = normalize(cross(ux, vN));
        const uy = normalize(cross(uz, ux));

        frames.push({ ux, uy, uz, ca });
      } else if (predictedModel.caPositions[i - 1]) {
        // Fallback using CA sequence derivative
        const ca = predictedModel.caPositions[i - 1];
        const nextCa = predictedModel.caPositions[i] || ca;
        const ux = normalize([nextCa[0] - ca[0], nextCa[1] - ca[1], nextCa[2] - ca[2]]);
        const uz: [number, number, number] = [0, 0, 1];
        const uy = normalize(cross(uz, ux));
        frames.push({ ux, uy, uz, ca });
      } else {
        frames.push(null);
      }
    }
    return frames;
  }, [predictedModel, N]);

  // Distance Matrices
  const predDistMatrix = useMemo(() => {
    if (!predictedModel) return [];
    const ca = predictedModel.caPositions;
    const mat: number[][] = [];
    for (let i = 0; i < N; i++) {
      mat[i] = [];
      for (let j = 0; j < N; j++) {
        const d = Math.hypot(ca[i][0] - ca[j][0], ca[i][1] - ca[j][1], ca[i][2] - ca[j][2]);
        mat[i][j] = d;
      }
    }
    return mat;
  }, [predictedModel, N]);

  const refDistMatrix = useMemo(() => {
    if (!refPdbData) return [];
    const ca = refPdbData.residues.map(r => r.caCoord);
    const len = Math.min(N, ca.length);
    const mat: number[][] = [];
    for (let i = 0; i < len; i++) {
      mat[i] = [];
      for (let j = 0; j < len; j++) {
        const d = Math.hypot(ca[i][0] - ca[j][0], ca[i][1] - ca[j][1], ca[i][2] - ca[j][2]);
        mat[i][j] = d;
      }
    }
    return mat;
  }, [refPdbData, N]);

  const mutantDistMatrix = useMemo(() => {
    if (!mutantModel) return [];
    const ca = mutantModel.caPositions;
    const mat: number[][] = [];
    for (let i = 0; i < N; i++) {
      mat[i] = [];
      for (let j = 0; j < N; j++) {
        const d = Math.hypot(ca[i][0] - ca[j][0], ca[i][1] - ca[j][1], ca[i][2] - ca[j][2]);
        mat[i][j] = d;
      }
    }
    return mat;
  }, [mutantModel, N]);

  // Orientation Matrix (angle in degrees of residue j relative to residue i's local frame)
  const orientationMatrix = useMemo(() => {
    if (localFrames.length === 0) return [];
    const mat: number[][] = [];
    for (let i = 0; i < N; i++) {
      mat[i] = [];
      const fi = localFrames[i];
      for (let j = 0; j < N; j++) {
        if (i === j || !fi || !localFrames[j]) {
          mat[i][j] = 0;
          continue;
        }
        const fj = localFrames[j]!;
        const dij: [number, number, number] = [
          fj.ca[0] - fi.ca[0],
          fj.ca[1] - fi.ca[1],
          fj.ca[2] - fi.ca[2]
        ];
        const dist = Math.hypot(dij[0], dij[1], dij[2]);
        if (dist < 1e-4) {
          mat[i][j] = 0;
          continue;
        }
        // Project dij into frame i
        const xProj = dot(dij, fi.ux);
        const cosTheta = Math.max(-1.0, Math.min(1.0, xProj / dist));
        const thetaDeg = (Math.acos(cosTheta) * 180) / Math.PI;
        mat[i][j] = thetaDeg;
      }
    }
    return mat;
  }, [localFrames, N]);

  // Color mapper for distance (Ångström scale)
  const getDistColor = (val: number, maxD: number = distanceCutoff): string => {
    const clamped = Math.max(0, Math.min(val, maxD));
    const t = clamped / maxD;
    // Blue (close <5Å) -> Teal -> Yellow -> Red -> Dark purple (>20Å)
    if (t < 0.25) {
      // 0 to 5 Å: deep blue to cyan
      const sub = t / 0.25;
      const r = Math.round(15 + sub * 20);
      const g = Math.round(60 + sub * 150);
      const b = Math.round(200 + sub * 50);
      return `rgb(${r},${g},${b})`;
    } else if (t < 0.5) {
      // 5 to 10 Å: cyan to yellow-green
      const sub = (t - 0.25) / 0.25;
      const r = Math.round(35 + sub * 180);
      const g = Math.round(210 + sub * 30);
      const b = Math.round(250 - sub * 200);
      return `rgb(${r},${g},${b})`;
    } else if (t < 0.75) {
      // 10 to 15 Å: yellow to red
      const sub = (t - 0.5) / 0.25;
      const r = Math.round(215 + sub * 40);
      const g = Math.round(240 - sub * 180);
      const b = Math.round(50 - sub * 30);
      return `rgb(${r},${g},${b})`;
    } else {
      // 15 to 20+ Å: red to dark purple/slate
      const sub = (t - 0.75) / 0.25;
      const r = Math.round(255 - sub * 190);
      const g = Math.round(60 - sub * 40);
      const b = Math.round(20 + sub * 70);
      return `rgb(${r},${g},${b})`;
    }
  };

  // Color mapper for orientation angle (0° to 180°)
  const getAngleColor = (angleDeg: number): string => {
    const t = Math.max(0, Math.min(angleDeg, 180)) / 180;
    // Violet (0°) -> Cyan (60°) -> Lime (120°) -> Amber (180°)
    const hue = Math.round((1 - t) * 280);
    return `hsl(${hue}, 85%, 50%)`;
  };

  // Color mapper for difference matrix (|A - B| in Å)
  const getDiffColor = (diffVal: number): string => {
    const clamped = Math.min(12, Math.abs(diffVal));
    const t = clamped / 12.0;
    // 0 = Dark charcoal -> Yellow -> Intense Neon Red (>8 Å)
    if (t < 0.1) return '#0d1117';
    if (t < 0.4) {
      const sub = (t - 0.1) / 0.3;
      return `rgb(${Math.round(40 + sub * 160)}, ${Math.round(40 + sub * 160)}, 40)`;
    }
    const sub = (t - 0.4) / 0.6;
    return `rgb(${Math.round(200 + sub * 55)}, ${Math.round(200 - sub * 170)}, ${Math.round(30 - sub * 20)})`;
  };

  // Metric displays
  const rmsdVal = validationReport?.caRmsdAngstrom ?? null;
  const tmScoreVal = validationReport?.tmScore ?? null;

  return (
    <div className="flex-1 w-full h-full min-h-0 bg-[#04050b] text-zinc-100 flex flex-col font-sans overflow-hidden">
      
      {/* Top Header & Metrics Bar */}
      <div className="bg-[#090b14] border-b border-zinc-800 px-5 py-3 flex flex-wrap items-center justify-between gap-4 shrink-0">
        <div className="flex items-center gap-3">
          <div className="w-8 h-8 rounded-lg bg-teal-950/70 border border-teal-500/40 flex items-center justify-center">
            <Grid className="w-5 h-5 text-teal-400" />
          </div>
          <div>
            <h1 className="text-sm font-bold text-white tracking-wide flex items-center gap-2">
              Karta Map Ciepła: Cα–Cα, Orientacja N-Cα-C &amp; Różnica (Diff)
              <span className="text-[10px] px-2 py-0.5 rounded bg-zinc-800 text-teal-300 font-mono border border-zinc-700">
                Macierz {N}×{N}
              </span>
            </h1>
            <p className="text-[11px] text-zinc-400">
              Topologia odległościowa, lokalne układy współrzędnych i dokładne porównanie modelu z krystalografią PDB / mutacją punktową.
            </p>
          </div>
        </div>

        {/* Global Protein & TM-Score / RMSD Display */}
        <div className="flex items-center gap-3 font-mono text-xs">
          <div className="px-3 py-1.5 rounded-lg bg-black/60 border border-zinc-800 flex items-center gap-3">
            <div>
              <span className="text-zinc-500 text-[10px] block">GLOBALNE Cα RMSD:</span>
              <span className={`text-sm font-extrabold ${rmsdVal && rmsdVal < 1.0 ? 'text-emerald-400' : 'text-amber-400'}`}>
                {rmsdVal !== null ? `${rmsdVal.toFixed(3)} Å` : 'Obliczanie...'}
              </span>
            </div>
            <div className="h-7 w-[1px] bg-zinc-800" />
            <div>
              <span className="text-zinc-500 text-[10px] block">TM-SCORE (0..1):</span>
              <span className={`text-sm font-extrabold ${tmScoreVal && tmScoreVal > 0.5 ? 'text-emerald-400' : 'text-cyan-400'}`}>
                {tmScoreVal !== null ? tmScoreVal.toFixed(4) : '—'}
              </span>
            </div>
          </div>

          {onSwitchToConstructor && (
            <button
              onClick={onSwitchToConstructor}
              className="px-3 py-1.5 bg-indigo-600 hover:bg-indigo-500 text-white rounded text-xs font-bold flex items-center gap-1.5 cursor-pointer transition"
            >
              <span>Konstruktor 3D</span>
              <ArrowRight className="w-3.5 h-3.5" />
            </button>
          )}
        </div>
      </div>

      {/* Subnav & Mode Selectors */}
      <div className="bg-[#0b0c16] border-b border-zinc-800/80 px-5 py-2.5 flex flex-wrap items-center justify-between gap-3 shrink-0 text-xs font-mono">
        <div className="flex items-center gap-1.5">
          <button
            onClick={() => setHeatmapMode('DISTANCE_CA')}
            className={`px-3 py-1 rounded-md font-bold transition cursor-pointer flex items-center gap-1.5 ${
              heatmapMode === 'DISTANCE_CA'
                ? 'bg-teal-500 text-black shadow'
                : 'bg-zinc-900 hover:bg-zinc-800 text-zinc-300'
            }`}
          >
            <span>(a) Odległości Cα–Cα</span>
          </button>

          <button
            onClick={() => setHeatmapMode('ORIENTATION_FRAME')}
            className={`px-3 py-1 rounded-md font-bold transition cursor-pointer flex items-center gap-1.5 ${
              heatmapMode === 'ORIENTATION_FRAME'
                ? 'bg-purple-500 text-white shadow'
                : 'bg-zinc-900 hover:bg-zinc-800 text-zinc-300'
            }`}
          >
            <span>(b) Orientacja Lokalna N-Cα-C</span>
          </button>

          <button
            onClick={() => setHeatmapMode('DIFF_COMPARISON')}
            className={`px-3 py-1 rounded-md font-bold transition cursor-pointer flex items-center gap-1.5 ${
              heatmapMode === 'DIFF_COMPARISON'
                ? 'bg-amber-500 text-black shadow'
                : 'bg-zinc-900 hover:bg-zinc-800 text-zinc-300'
            }`}
          >
            <GitCompare className="w-3.5 h-3.5" />
            <span>(c) Tryb Porównania &amp; Diff</span>
          </button>
        </div>

        {/* Dynamic Controls based on selected mode */}
        <div className="flex items-center gap-3">
          {heatmapMode === 'DIFF_COMPARISON' && (
            <div className="flex items-center gap-2 bg-black/60 px-2.5 py-1 rounded border border-zinc-800">
              <span className="text-zinc-500 text-[11px]">Cel:</span>
              <button
                onClick={() => setCompareTarget('PREDICTED_VS_PDB')}
                className={`px-2 py-0.5 rounded text-[11px] font-bold ${
                  compareTarget === 'PREDICTED_VS_PDB' ? 'bg-cyan-600 text-white' : 'text-zinc-400 hover:text-white'
                }`}
              >
                Model vs Referencja PDB
              </button>
              <button
                onClick={() => setCompareTarget('HEALTHY_VS_MUTANT')}
                className={`px-2 py-0.5 rounded text-[11px] font-bold ${
                  compareTarget === 'HEALTHY_VS_MUTANT' ? 'bg-rose-600 text-white' : 'text-zinc-400 hover:text-white'
                }`}
              >
                Zdrowe vs Mutacja
              </button>
            </div>
          )}

          {heatmapMode === 'DIFF_COMPARISON' && compareTarget === 'HEALTHY_VS_MUTANT' && (
            <div className="flex items-center gap-1.5 bg-black/60 px-2 py-1 rounded border border-zinc-800">
              <span className="text-[10px] text-zinc-500">Mutacja w:</span>
              <input
                type="number"
                min={1}
                max={N}
                value={mutantSite}
                onChange={(e) => setMutantSite(parseInt(e.target.value) || 1)}
                className="w-12 bg-zinc-900 border border-zinc-700 px-1 py-0.5 text-xs text-center rounded text-white"
              />
              <span className="text-zinc-400 text-xs">➔</span>
              <select
                value={mutantNewAA}
                onChange={(e) => setMutantNewAA(e.target.value)}
                className="bg-zinc-900 border border-zinc-700 text-xs px-1.5 py-0.5 rounded text-white"
              >
                {['ARG', 'ALA', 'CYS', 'GLU', 'LYS', 'PRO', 'SER', 'TRP', 'VAL'].map(aa => (
                  <option key={aa} value={aa}>{aa}</option>
                ))}
              </select>
            </div>
          )}

          <div className="flex items-center gap-2">
            <span className="text-zinc-500 text-[10px]">Skala odległości:</span>
            <input
              type="range"
              min={8}
              max={32}
              value={distanceCutoff}
              onChange={(e) => setDistanceCutoff(parseFloat(e.target.value))}
              className="w-20 accent-teal-400"
            />
            <span className="text-zinc-300 text-[11px] w-10">{distanceCutoff.toFixed(0)} Å</span>
          </div>
        </div>
      </div>

      {/* Main Heatmap Canvas Area */}
      <div className="flex-1 min-h-0 flex flex-col lg:flex-row overflow-hidden p-4 gap-4 bg-[#05060b]">
        
        {/* Heatmap Matrix Display */}
        <div className="flex-1 flex flex-col bg-[#070912] border border-zinc-800 rounded-2xl p-4 overflow-hidden shadow-2xl relative min-h-0">
          
          {/* Header of the current heatmap view */}
          <div className="flex items-center justify-between mb-3 shrink-0 font-mono text-xs">
            <div className="flex items-center gap-2 text-zinc-300">
              {heatmapMode === 'DISTANCE_CA' && (
                <span>Mapa Odległości Cα(i)—Cα(j) [Å] • Cała Struktura 3D</span>
              )}
              {heatmapMode === 'ORIENTATION_FRAME' && (
                <span>Kąt Orientacji Reszty j w Układzie Lokalnym N-Cα-C Reszty i [stopnie 0°..180°]</span>
              )}
              {heatmapMode === 'DIFF_COMPARISON' && (
                <span>
                  {compareTarget === 'PREDICTED_VS_PDB'
                    ? `Tryb Różnicowy: |D_model - D_pdb| względem ${refPdbData?.pdbId || 'PDB'}`
                    : `Tryb Różnicowy: |D_wild - D_mutant(${mutantCodes[mutantSite-1]}${mutantSite}${mutantNewAA})|`}
                </span>
              )}
            </div>

            {/* Hovered cell tooltip */}
            {hoveredCell ? (
              <div className="px-2.5 py-1 rounded bg-black border border-cyan-500/40 text-[11px] text-cyan-300 font-bold flex items-center gap-2">
                <span>Reszty #{hoveredCell.i + 1} ({codes[hoveredCell.i]}) ↔ #{hoveredCell.j + 1} ({codes[hoveredCell.j]})</span>
                <span className="text-white font-extrabold">
                  {heatmapMode === 'ORIENTATION_FRAME'
                    ? `${hoveredCell.val1.toFixed(1)}°`
                    : `${hoveredCell.val1.toFixed(2)} Å`}
                </span>
                {hoveredCell.diff !== undefined && (
                  <span className="text-amber-400">Δ = {hoveredCell.diff.toFixed(2)} Å</span>
                )}
              </div>
            ) : (
              <span className="text-[11px] text-zinc-500">Najedź kursorem na komórkę, aby odczytać wartość</span>
            )}
          </div>

          {/* Grid View */}
          <div className="flex-1 min-h-0 flex items-center justify-center overflow-auto p-2">
            {heatmapMode === 'DIFF_COMPARISON' ? (
              // Three-matrix comparison layout (A, B, Diff)
              <div className="flex items-center justify-center gap-4 flex-wrap max-w-full max-h-full">
                
                {/* Matrix 1: Model / Healthy */}
                <div className="flex flex-col items-center">
                  <span className="text-[10px] font-mono text-zinc-400 mb-1">
                    {compareTarget === 'PREDICTED_VS_PDB' ? '1. Przewidywana (Model)' : '1. Typ Dziki (Zdrowe)'}
                  </span>
                  <div 
                    className="grid border border-zinc-700 bg-black"
                    style={{
                      gridTemplateColumns: `repeat(${N}, minmax(0, 1fr))`,
                      width: N <= 46 ? '240px' : '280px',
                      height: N <= 46 ? '240px' : '280px'
                    }}
                  >
                    {predDistMatrix.map((row, i) =>
                      row.map((val, j) => (
                        <div
                          key={`p-${i}-${j}`}
                          onMouseEnter={() => setHoveredCell({ i, j, val1: val })}
                          onMouseLeave={() => setHoveredCell(null)}
                          style={{ backgroundColor: getDistColor(val) }}
                          className="w-full h-full"
                        />
                      ))
                    )}
                  </div>
                </div>

                {/* Matrix 2: PDB Reference / Mutant */}
                <div className="flex flex-col items-center">
                  <span className="text-[10px] font-mono text-zinc-400 mb-1">
                    {compareTarget === 'PREDICTED_VS_PDB' ? '2. Referencja PDB (Krystalografia)' : `2. Wariant Mutacyjny (${mutantNewAA})`}
                  </span>
                  <div 
                    className="grid border border-zinc-700 bg-black"
                    style={{
                      gridTemplateColumns: `repeat(${N}, minmax(0, 1fr))`,
                      width: N <= 46 ? '240px' : '280px',
                      height: N <= 46 ? '240px' : '280px'
                    }}
                  >
                    {(compareTarget === 'PREDICTED_VS_PDB' ? refDistMatrix : mutantDistMatrix).map((row, i) =>
                      row.map((val, j) => (
                        <div
                          key={`r-${i}-${j}`}
                          onMouseEnter={() => setHoveredCell({ i, j, val1: val })}
                          onMouseLeave={() => setHoveredCell(null)}
                          style={{ backgroundColor: getDistColor(val) }}
                          className="w-full h-full"
                        />
                      ))
                    )}
                  </div>
                </div>

                {/* Matrix 3: Difference (Δ Diff) */}
                <div className="flex flex-col items-center">
                  <span className="text-[10px] font-mono text-amber-300 font-bold mb-1 flex items-center gap-1">
                    <TrendingDown className="w-3 h-3 text-amber-400" />
                    <span>3. Różnica Bezwzględna |Δ|</span>
                  </span>
                  <div 
                    className="grid border border-amber-600/70 bg-black shadow-lg shadow-amber-500/10"
                    style={{
                      gridTemplateColumns: `repeat(${N}, minmax(0, 1fr))`,
                      width: N <= 46 ? '240px' : '280px',
                      height: N <= 46 ? '240px' : '280px'
                    }}
                  >
                    {predDistMatrix.map((row, i) =>
                      row.map((val1, j) => {
                        const targetMat = compareTarget === 'PREDICTED_VS_PDB' ? refDistMatrix : mutantDistMatrix;
                        const val2 = targetMat[i]?.[j] ?? val1;
                        const diff = Math.abs(val1 - val2);
                        return (
                          <div
                            key={`d-${i}-${j}`}
                            onMouseEnter={() => setHoveredCell({ i, j, val1, val2, diff })}
                            onMouseLeave={() => setHoveredCell(null)}
                            style={{ backgroundColor: getDiffColor(diff) }}
                            className="w-full h-full cursor-pointer hover:ring-1 hover:ring-white"
                          />
                        );
                      })
                    )}
                  </div>
                </div>

              </div>
            ) : (
              // Single full matrix
              <div 
                className="grid border border-zinc-700 bg-black shadow-2xl relative"
                style={{
                  gridTemplateColumns: `repeat(${N}, minmax(0, 1fr))`,
                  width: 'min(90vw, 520px)',
                  height: 'min(90vw, 520px)'
                }}
              >
                {(heatmapMode === 'DISTANCE_CA' ? predDistMatrix : orientationMatrix).map((row, i) =>
                  row.map((val, j) => {
                    const bgColor = heatmapMode === 'DISTANCE_CA' ? getDistColor(val) : getAngleColor(val);
                    return (
                      <div
                        key={`cell-${i}-${j}`}
                        onMouseEnter={() => setHoveredCell({ i, j, val1: val })}
                        onMouseLeave={() => setHoveredCell(null)}
                        style={{ backgroundColor: bgColor }}
                        className="w-full h-full hover:ring-2 hover:ring-white hover:z-10 cursor-pointer transition-none"
                      />
                    );
                  })
                )}
              </div>
            )}
          </div>

          {/* Interactive Colorbar / Legend */}
          <div className="mt-2 pt-2 border-t border-zinc-800/80 flex flex-wrap items-center justify-between gap-3 font-mono text-[11px] shrink-0">
            {heatmapMode === 'DISTANCE_CA' && (
              <div className="flex items-center gap-2">
                <span className="text-zinc-500">Legenda odległości:</span>
                <div className="flex items-center gap-1">
                  <div className="w-4 h-3 rounded" style={{ backgroundColor: getDistColor(2) }} />
                  <span className="text-zinc-400">Blisko (&lt;4 Å)</span>
                  <div className="w-4 h-3 rounded ml-2" style={{ backgroundColor: getDistColor(8) }} />
                  <span className="text-zinc-400">Kontakt (8 Å)</span>
                  <div className="w-4 h-3 rounded ml-2" style={{ backgroundColor: getDistColor(14) }} />
                  <span className="text-zinc-400">Średnie (14 Å)</span>
                  <div className="w-4 h-3 rounded ml-2" style={{ backgroundColor: getDistColor(22) }} />
                  <span className="text-zinc-400">Dalekie (&gt;20 Å)</span>
                </div>
              </div>
            )}

            {heatmapMode === 'ORIENTATION_FRAME' && (
              <div className="flex items-center gap-2">
                <span className="text-zinc-500">Kąt w układzie N-Cα-C:</span>
                <div className="flex items-center gap-1">
                  <div className="w-4 h-3 rounded" style={{ backgroundColor: getAngleColor(0) }} />
                  <span className="text-zinc-400">0° (przód Cα→C)</span>
                  <div className="w-4 h-3 rounded ml-2" style={{ backgroundColor: getAngleColor(90) }} />
                  <span className="text-zinc-400">90° (prostopadły)</span>
                  <div className="w-4 h-3 rounded ml-2" style={{ backgroundColor: getAngleColor(180) }} />
                  <span className="text-zinc-400">180° (tył)</span>
                </div>
              </div>
            )}

            {heatmapMode === 'DIFF_COMPARISON' && (
              <div className="flex items-center gap-2">
                <span className="text-zinc-500">Skala różnicy Δ:</span>
                <div className="flex items-center gap-1">
                  <div className="w-4 h-3 rounded bg-[#0d1117] border border-zinc-700" />
                  <span className="text-zinc-400">0 Å (idealna zgodność)</span>
                  <div className="w-4 h-3 rounded ml-2 bg-yellow-400" />
                  <span className="text-zinc-400">~4 Å</span>
                  <div className="w-4 h-3 rounded ml-2 bg-rose-600" />
                  <span className="text-zinc-400">&gt;8 Å (rozbieżność)</span>
                </div>
              </div>
            )}

            <span className="text-zinc-500 text-[10px]">
              Oś X: Reszta j (1..{N}) • Oś Y: Reszta i (1..{N})
            </span>
          </div>

        </div>

        {/* Right Panel: Structural Topology Diagnostics & TM-Score Summary */}
        <div className="w-full lg:w-80 flex flex-col gap-3 font-mono text-xs shrink-0 overflow-y-auto">
          
          <div className="p-4 rounded-xl bg-[#090b14] border border-zinc-800 space-y-3">
            <h3 className="text-xs font-bold text-zinc-300 uppercase tracking-wider flex items-center gap-2">
              <Sparkles className="w-4 h-4 text-teal-400" />
              <span>Metryki Topologii Globalnej</span>
            </h3>

            <div className="p-3 rounded-lg bg-black/60 border border-zinc-800 space-y-2">
              <div className="flex items-center justify-between">
                <span className="text-zinc-500">TM-SCORE:</span>
                <span className={`font-extrabold ${tmScoreVal && tmScoreVal > 0.5 ? 'text-emerald-400' : 'text-zinc-300'}`}>
                  {tmScoreVal !== null ? tmScoreVal.toFixed(4) : '—'}
                </span>
              </div>
              <div className="text-[10px] text-zinc-500 leading-relaxed">
                TM-score &gt; 0.5 oznacza identyczny fold topologiczny. Niezależny od wielkości białka.
              </div>
              <div className="flex items-center justify-between pt-1 border-t border-zinc-800/80">
                <span className="text-zinc-500">Cα RMSD (KABSCH):</span>
                <span className={`font-extrabold ${rmsdVal && rmsdVal < 1.0 ? 'text-emerald-400' : 'text-amber-400'}`}>
                  {rmsdVal !== null ? `${rmsdVal.toFixed(3)} Å` : '—'}
                </span>
              </div>
              {validationReport && (
                <div className="flex items-center justify-between">
                  <span className="text-zinc-500">GDT-TS:</span>
                  <span className="text-zinc-300 font-bold">
                    {validationReport.gdtTsPct.toFixed(1)}%
                  </span>
                </div>
              )}
            </div>

            <div className="p-3 rounded-lg bg-black/60 border border-zinc-800 space-y-1.5 text-[11px]">
              <div className="text-zinc-400 font-bold">Lokalny Układ N-Cα-C:</div>
              <p className="text-zinc-500 text-[10px] leading-relaxed">
                Każda reszta i tworzy układ kartezjański z wektorem <strong>u_x</strong> wzdłuż wiązania Cα→C oraz płaszczyzną zdefiniowaną przez atom N. Mapa orientacji ujawnia relacje kątowe helis i harmonijek.
              </p>
            </div>
          </div>

          <div className="p-4 rounded-xl bg-[#090b14] border border-zinc-800 space-y-2 text-[11px]">
            <span className="text-zinc-400 font-bold block">Kontakty dalekiego zasięgu (|i - j| &gt; 5, d &lt; 8 Å):</span>
            <div className="p-2 rounded bg-black/60 border border-zinc-800 text-[10px] text-zinc-300">
              Liczba par w kontakcie: {
                predDistMatrix.reduce((acc, row, i) => 
                  acc + row.filter((d, j) => j > i + 5 && d < 8.0).length, 0
                )
              }
            </div>
          </div>

        </div>

      </div>

    </div>
  );
};
