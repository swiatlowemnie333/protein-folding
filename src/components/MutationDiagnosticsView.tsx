import React, { useEffect, useMemo, useState } from 'react';
import { 
  ShieldCheck, 
  AlertTriangle, 
  Activity, 
  RotateCcw, 
  Dna, 
  Flame, 
  Zap, 
  Sliders, 
  ArrowRight, 
  CheckCircle2, 
  XCircle,
  Eye,
  Layers,
  Sparkles,
  RefreshCw
} from 'lucide-react';
import { SNAPLOCK_22_DATABASE, E8_PATCHES } from '../snaplockEngine';
import {
  GLOBAL_PROTEIN_TARGETS,
  GlobalProteinTargetId,
  TargetMutationVariant
} from '../snaplockGlobalTargets';

interface MutationDiagnosticsViewProps {
  activeTarget?: GlobalProteinTargetId;
  onLoadVariantToConstructor?: (sequence: string) => void;
  onSwitchToConstructor?: () => void;
}

interface ClinicalMutationPreset extends TargetMutationVariant {
  healthyEnergy: number;
  mutatedEnergy: number;
}

function estimateChainCohesiveEnergy(seq: string[], ssPairs: [number, number][], mutatedPos1Based?: number): number {
  let e = -seq.length * 1.45;
  for (const code of seq) {
    const row = SNAPLOCK_22_DATABASE.find(a => a.code3 === code);
    if (row && row.hydropathyKyte > 1.2) {
      e -= row.hydropathyKyte * 0.85;
    }
  }
  for (const [i1, i2] of ssPairs) {
    const intact =
      (!mutatedPos1Based || (mutatedPos1Based !== i1 + 1 && mutatedPos1Based !== i2 + 1)) &&
      (seq[i1] === 'CYS' || seq[i1] === 'SEC') &&
      (seq[i2] === 'CYS' || seq[i2] === 'SEC');
    if (intact) {
      e -= 14.2;
    } else {
      e += 8.6;
    }
  }
  if (mutatedPos1Based) {
    e += 12.4;
  }
  return parseFloat(e.toFixed(1));
}

export const MutationDiagnosticsView: React.FC<MutationDiagnosticsViewProps> = ({
  activeTarget = 'UBIQUITIN_76',
  onLoadVariantToConstructor,
  onSwitchToConstructor
}) => {
  const targetSpec = useMemo(() => {
    const key = activeTarget && activeTarget !== 'CUSTOM' ? activeTarget : 'UBIQUITIN_76';
    return GLOBAL_PROTEIN_TARGETS[key];
  }, [activeTarget]);

  const mutationPresets: ClinicalMutationPreset[] = useMemo(() => {
    const baseHealthyE = estimateChainCohesiveEnergy(
      targetSpec.sequence3Letter,
      targetSpec.explicitDisulfidePairs0Based
    );
    return targetSpec.mutations.map(m => {
      const mutSeq = targetSpec.sequence3Letter.map((r, idx) =>
        idx + 1 === m.position ? m.mutatedResidue : r
      );
      const mutE = estimateChainCohesiveEnergy(
        mutSeq,
        targetSpec.explicitDisulfidePairs0Based,
        m.position
      );
      return {
        ...m,
        healthyEnergy: baseHealthyE,
        mutatedEnergy: mutE
      };
    });
  }, [targetSpec]);

  const [selectedMutation, setSelectedMutation] = useState<ClinicalMutationPreset>(
    () => mutationPresets[0]
  );

  // Immediately sync selected mutation whenever global activeTarget switches!
  useEffect(() => {
    if (mutationPresets.length > 0) {
      setSelectedMutation(mutationPresets[0]);
    }
  }, [mutationPresets]);

  const [syncCamera, setSyncCamera] = useState<boolean>(true);
  
  // Independent / Synchronized Orbit Rotations
  const [rotX1, setRotX1] = useState<number>(-20);
  const [rotY1, setRotY1] = useState<number>(35);
  const [rotX2, setRotX2] = useState<number>(-20);
  const [rotY2, setRotY2] = useState<number>(35);

  const [zoomScale, setZoomScale] = useState<number>(1.0);
  const [hoveredResidue, setHoveredResidue] = useState<number | null>(null);

  // Generate sequences strictly from active global target
  const healthySequence = targetSpec.sequence3Letter;
  const mutatedSequence = healthySequence.map((res, idx) => {
    if (idx + 1 === selectedMutation.position) {
      return selectedMutation.mutatedResidue;
    }
    return res;
  });

  const deltaEnergy = selectedMutation.mutatedEnergy - selectedMutation.healthyEnergy;

  // Simple 3D projection for dual viewport
  const project3D = (x: number, y: number, z: number, rX: number, rY: number, cx: number, cy: number, scale: number) => {
    const radX = (rX * Math.PI) / 180;
    const radY = (rY * Math.PI) / 180;

    // Rotate Y
    const x1 = x * Math.cos(radY) + z * Math.sin(radY);
    const z1 = -x * Math.sin(radY) + z * Math.cos(radY);

    // Rotate X
    const y2 = y * Math.cos(radX) - z1 * Math.sin(radX);
    const z2 = y * Math.sin(radX) + z1 * Math.cos(radX);

    const fov = 450;
    const distance = 420;
    const projScale = (fov / (fov + z2 + distance)) * scale;

    return {
      x: cx + x1 * projScale,
      y: cy + y2 * projScale,
      z: z2,
      scale: projScale
    };
  };

  // Generate 3D coordinates scaled dynamically to any protein length (18..129 aa)
  const generateChainCoordinates = (seq: string[], isMutated: boolean) => {
    const len = Math.max(1, seq.length);
    const yStep = Math.min(14, 240 / len);
    return seq.map((res, i) => {
      const turns = Math.max(2.4, len / 9.5);
      const theta = (i / len) * Math.PI * 2 * turns;
      const radius = 84 - (i % 4) * 10;
      let x = Math.cos(theta) * radius;
      let y = (i - len / 2) * yStep;
      let z = Math.sin(theta) * radius;

      if (isMutated && Math.abs(i + 1 - selectedMutation.position) <= 4) {
        const repulsion = (5 - Math.abs(i + 1 - selectedMutation.position)) * 18;
        x += Math.cos(theta) * repulsion;
        z += Math.sin(theta) * repulsion;
        y += (i % 2 === 0 ? 1 : -1) * repulsion * 0.7;
      }

      return {
        index: i + 1,
        code: res,
        x,
        y,
        z,
        isMutatedNode: isMutated && i + 1 === selectedMutation.position,
        isDisulfideCys: res === 'CYS'
      };
    });
  };

  const healthyNodes = generateChainCoordinates(healthySequence, false);
  const mutatedNodes = generateChainCoordinates(mutatedSequence, true);

  const handleRotate = (dx: number, dy: number, isRightPanel: boolean) => {
    if (syncCamera) {
      setRotY1(prev => (prev + dx * 0.5) % 360);
      setRotX1(prev => Math.max(-85, Math.min(85, prev - dy * 0.5)));
      setRotY2(prev => (prev + dx * 0.5) % 360);
      setRotX2(prev => Math.max(-85, Math.min(85, prev - dy * 0.5)));
    } else if (isRightPanel) {
      setRotY2(prev => (prev + dx * 0.5) % 360);
      setRotX2(prev => Math.max(-85, Math.min(85, prev - dy * 0.5)));
    } else {
      setRotY1(prev => (prev + dx * 0.5) % 360);
      setRotX1(prev => Math.max(-85, Math.min(85, prev - dy * 0.5)));
    }
  };

  const render3DViewport = (
    nodes: ReturnType<typeof generateChainCoordinates>, 
    isPatient: boolean, 
    rotX: number, 
    rotY: number
  ) => {
    const width = 500;
    const height = 360;
    const cx = width / 2;
    const cy = height / 2;

    const projected = nodes.map(n => {
      const p = project3D(n.x, n.y, n.z, rotX, rotY, cx, cy, zoomScale * 1.5);
      return { ...n, px: p.x, py: p.y, pz: p.z, pScale: p.scale };
    });

    // Sort depth
    const sorted = [...projected].sort((a, b) => a.pz - b.pz);

    return (
      <svg 
        viewBox={`0 0 ${width} ${height}`} 
        className="w-full h-full cursor-grab active:cursor-grabbing select-none"
        onMouseDown={(e) => {
          const startX = e.clientX;
          const startY = e.clientY;
          const onMove = (me: MouseEvent) => {
            const dx = me.clientX - startX;
            const dy = me.clientY - startY;
            handleRotate(dx, dy, isPatient);
          };
          const onUp = () => {
            window.removeEventListener('mousemove', onMove);
            window.removeEventListener('mouseup', onUp);
          };
          window.addEventListener('mousemove', onMove);
          window.addEventListener('mouseup', onUp);
        }}
      >
        {/* Subtle grid sphere */}
        <circle cx={cx} cy={cy} r={140} fill="none" stroke="#121324" strokeWidth="1" strokeDasharray="3,4" />
        <circle cx={cx} cy={cy} r={70} fill="none" stroke="#0b0c16" strokeWidth="1" strokeDasharray="2,4" />

        {/* Backbone connections (peptides) */}
        {projected.map((n, i) => {
          if (i === projected.length - 1) return null;
          const next = projected[i + 1];
          const isNearMutation = isPatient && (n.isMutatedNode || next.isMutatedNode);

          return (
            <line
              key={`backbone-${i}`}
              x1={n.px}
              y1={n.py}
              x2={next.px}
              y2={next.py}
              stroke={isNearMutation ? '#ef4444' : isPatient ? '#64748b' : '#38bdf8'}
              strokeWidth={isNearMutation ? '3' : '2'}
              strokeDasharray={isNearMutation ? '4,3' : 'none'}
              opacity={isNearMutation ? 0.9 : 0.65}
            />
          );
        })}

        {/* Dynamic Disulfide Bridges from active global target */}
        {targetSpec.explicitDisulfidePairs0Based.map(([idxA, idxB], bIdx) => {
          const pA = projected[idxA];
          const pB = projected[idxB];
          if (!pA || !pB) return null;
          const isBroken =
            isPatient &&
            (selectedMutation.position === idxA + 1 || selectedMutation.position === idxB + 1);
          return (
            <g key={`ss-${bIdx}`}>
              <line
                x1={pA.px}
                y1={pA.py}
                x2={pB.px}
                y2={pB.py}
                stroke={isBroken ? '#ef4444' : '#10b981'}
                strokeWidth={isBroken ? '2.5' : '3.6'}
                strokeDasharray={isBroken ? '4,4' : 'none'}
                strokeLinecap="round"
              />
              <circle
                cx={(pA.px + pB.px) / 2}
                cy={(pA.py + pB.py) / 2}
                r={isBroken ? 6 : 8}
                fill={isBroken ? '#7f1d1d' : '#065f46'}
                stroke={isBroken ? '#ef4444' : '#34d399'}
                strokeWidth="1.5"
              />
              <text
                x={(pA.px + pB.px) / 2}
                y={(pA.py + pB.py) / 2 + 3}
                textAnchor="middle"
                className="fill-white font-mono text-[7px] font-bold"
              >
                {isBroken ? '✕ S-S' : 'S-S'}
              </text>
            </g>
          );
        })}

        {/* Amino Acid Nodes */}
        {sorted.map((n) => {
          const rad = Math.max(7, Math.min(18, 11 * n.pScale));
          const isMut = n.isMutatedNode;
          const isCys = n.isDisulfideCys;

          let fill = '#3b82f6';
          if (isMut) fill = '#ef4444';
          else if (isCys) fill = '#10b981';
          else if (['VAL', 'LEU', 'ILE', 'PHE', 'MET'].includes(n.code)) fill = '#f59e0b';
          else if (['ASP', 'GLU'].includes(n.code)) fill = '#ec4899';
          else if (['LYS', 'ARG', 'HIS'].includes(n.code)) fill = '#8b5cf6';

          return (
            <g 
              key={`node-${n.index}`}
              onMouseEnter={() => setHoveredResidue(n.index)}
              onMouseLeave={() => setHoveredResidue(null)}
              className="cursor-pointer"
            >
              {/* Outer glow if mutated or active bridge */}
              {(isMut || isCys) && (
                <circle cx={n.px} cy={n.py} r={rad + 5} fill={isMut ? '#ef4444' : '#10b981'} opacity="0.25" />
              )}
              
              <circle
                cx={n.px}
                cy={n.py}
                r={rad}
                fill={fill}
                stroke={isMut ? '#fee2e2' : '#000'}
                strokeWidth={isMut ? '2' : '1.2'}
              />
              
              <text
                x={n.px}
                y={n.py + 3}
                textAnchor="middle"
                className="fill-black font-mono text-[8.5px] font-extrabold select-none pointer-events-none"
              >
                {n.code}
              </text>

              <text
                x={n.px}
                y={n.py - rad - 3}
                textAnchor="middle"
                className={`font-mono text-[8px] font-bold select-none pointer-events-none ${
                  isMut ? 'fill-red-400' : 'fill-zinc-400'
                }`}
              >
                #{n.index}
              </text>
            </g>
          );
        })}
      </svg>
    );
  };

  return (
    <div className="flex-1 w-full h-full min-h-0 bg-[#020205] text-slate-200 flex flex-col overflow-hidden font-mono select-none">
      
      {/* HEADER BAR */}
      <div className="border-b border-zinc-900 bg-[#04040a] p-4 shrink-0 flex flex-wrap items-center justify-between gap-3">
        <div className="flex items-center gap-3">
          <div className="w-10 h-10 rounded-xl bg-red-950/80 border border-red-500/50 flex items-center justify-center text-red-400 shadow-inner">
            <Activity className="w-5 h-5 text-red-300" />
          </div>
          <div>
            <div className="flex items-center gap-2">
              <h1 className="text-base sm:text-lg font-bold text-white tracking-tight flex items-center gap-2">
                Diagnostyka Mutacji: Zdrowe vs Chore (Widok Dwuokienny)
              </h1>
              <span className="text-[10px] px-2 py-0.5 rounded bg-red-950 text-red-300 border border-red-700/60 font-bold uppercase">
                Klamry S-S & Destabilizacja Rdzenia
              </span>
            </div>
            <p className="text-xs text-zinc-400 mt-0.5">
              Zsynchronizowany cel: <strong className="text-cyan-300">{targetSpec.titlePl}</strong> ({targetSpec.residueCount} aa • {targetSpec.pdbReferenceId})
            </p>
          </div>
        </div>

        {/* MUTATION PRESET DROPDOWN */}
        <div className="flex items-center gap-2 flex-wrap">
          <span className="text-xs text-zinc-400 font-bold uppercase">Wybierz Wariant Mutacji:</span>
          <select
            value={selectedMutation.id}
            onChange={(e) => {
              const found = mutationPresets.find(m => m.id === e.target.value);
              if (found) setSelectedMutation(found);
            }}
            className="bg-black border border-red-900/80 text-white font-mono text-xs rounded px-3 py-1.5 focus:outline-none focus:border-red-500"
          >
            {mutationPresets.map(m => (
              <option key={m.id} value={m.id}>{m.name}</option>
            ))}
          </select>

          <button
            onClick={() => setSyncCamera(prev => !prev)}
            className={`px-3 py-1.5 rounded text-xs font-bold border transition cursor-pointer flex items-center gap-1.5 ${
              syncCamera 
                ? 'bg-indigo-950/80 border-indigo-700 text-indigo-300 shadow-sm' 
                : 'bg-zinc-900 border-zinc-800 text-zinc-400 hover:text-white'
            }`}
            title="Synchronizuj obrót kamery 3D pomiędzy oknem zdrowym i chorym"
          >
            <RotateCcw className="w-3.5 h-3.5" />
            <span>{syncCamera ? 'Kamery: Zsynchronizowane' : 'Kamery: Niezależne'}</span>
          </button>

          {onLoadVariantToConstructor && (
            <button
              onClick={() => {
                onLoadVariantToConstructor(mutatedSequence.join(' '));
                if (onSwitchToConstructor) onSwitchToConstructor();
              }}
              className="px-3 py-1.5 bg-gradient-to-r from-red-600 to-amber-600 hover:from-red-500 hover:to-amber-500 text-white rounded text-xs font-bold flex items-center gap-1.5 cursor-pointer shadow-md transition"
              title="Przenieś ten zmutowany łańcuch do pełnego Konstruktora 3D"
            >
              <Dna className="w-3.5 h-3.5" />
              <span>Wczytaj do Konstruktora 3D</span>
            </button>
          )}
        </div>
      </div>

      {/* CLINICAL INFO BANNER */}
      <div className="bg-[#05040a] border-b border-zinc-900/80 px-4 py-2.5 flex items-center justify-between text-xs font-mono">
        <div className="flex items-center gap-2">
          <AlertTriangle className="w-4 h-4 text-amber-400 shrink-0" />
          <span className="text-zinc-400 font-bold uppercase">{selectedMutation.gene}:</span>
          <span className="text-slate-300">{selectedMutation.clinicalImpact}</span>
        </div>
        <div className="text-zinc-500 text-[11px] shrink-0 hidden md:block">
          Pozycja: <strong className="text-white">#{selectedMutation.position}</strong> ({selectedMutation.healthyResidue} ➔ {selectedMutation.mutatedResidue})
        </div>
      </div>

      {/* DUAL WINDOW COMPARISON AREA */}
      <div className="flex-1 grid grid-cols-1 md:grid-cols-2 divide-y md:divide-y-0 md:divide-x divide-zinc-900 min-h-0 overflow-hidden">
        
        {/* LEFT WINDOW: HEALTHY / WILD TYPE */}
        <div className="flex flex-col h-full min-h-0 bg-[#020308] relative overflow-hidden">
          <div className="p-3 bg-[#030510] border-b border-zinc-900 flex items-center justify-between shrink-0">
            <div className="flex items-center gap-2">
              <span className="w-2.5 h-2.5 rounded-full bg-emerald-400" />
              <h2 className="text-xs font-bold text-emerald-300 uppercase tracking-wider">
                1. Wariant Natywny (Zdrowy / Wild-Type • {healthySequence.length} aa)
              </h2>
            </div>
            <span className="text-[10px] px-2 py-0.5 rounded bg-emerald-950/80 border border-emerald-700/60 text-emerald-300 font-bold">
              Konformacja Stabilna
            </span>
          </div>

          {/* 3D CANVAS */}
          <div className="flex-1 w-full h-full min-h-0 relative flex items-center justify-center overflow-hidden">
            {render3DViewport(healthyNodes, false, rotX1, rotY1)}

            {/* Overlaid Badges */}
            <div className="absolute top-3 left-3 bg-black/80 backdrop-blur-sm border border-emerald-900/60 rounded-lg p-2 text-[11px] space-y-1 pointer-events-none">
              <div className="flex justify-between gap-4">
                <span className="text-zinc-400">Energia Kohezyjna:</span>
                <strong className="text-emerald-400 font-extrabold">{selectedMutation.healthyEnergy} kcal/mol</strong>
              </div>
              <div className="flex justify-between gap-4">
                <span className="text-zinc-400">Klamry S-S (Mostki):</span>
                <strong className="text-emerald-300">{targetSpec.explicitDisulfidePairs0Based.length} nienaruszone</strong>
              </div>
              <div className="flex justify-between gap-4">
                <span className="text-zinc-400">Rdzeń Hydrofobowy:</span>
                <strong className="text-emerald-300">Zwarty ({targetSpec.pdbReferenceId})</strong>
              </div>
            </div>

            <div className="absolute bottom-2 left-3 text-[9.5px] text-zinc-500 font-mono pointer-events-none">
              Przeciągnij myszką, aby obrócić model 3D
            </div>
          </div>

          {/* Track verification mini footer */}
          <div className="p-2.5 bg-[#03040a] border-t border-zinc-900 text-[10.5px] flex items-center justify-between text-zinc-400">
            <span className="flex items-center gap-1 text-emerald-400 font-bold">
              <CheckCircle2 className="w-3.5 h-3.5" />
              Track 5: {targetSpec.explicitDisulfidePairs0Based.length > 0 ? 'Klamry kowalencyjne S-S aktywne' : 'Rdzeń hydrofobowy i wiązania H aktywne'}
            </span>
            <span className="text-zinc-500">
              Reszta #{selectedMutation.position}: <strong className="text-emerald-300">{selectedMutation.healthyResidue}</strong> (Natywna)
            </span>
          </div>
        </div>

        {/* RIGHT WINDOW: MUTATED / SICK PATIENT */}
        <div className="flex flex-col h-full min-h-0 bg-[#080204] relative overflow-hidden">
          <div className="p-3 bg-[#130307] border-b border-zinc-900 flex items-center justify-between shrink-0">
            <div className="flex items-center gap-2">
              <span className="w-2.5 h-2.5 rounded-full bg-red-500" />
              <h2 className="text-xs font-bold text-red-400 uppercase tracking-wider">
                2. Wariant Zmutowany (Chory / Pacjent)
              </h2>
            </div>
            <span className="text-[10px] px-2 py-0.5 rounded bg-red-950/80 border border-red-700/60 text-red-300 font-bold">
              ALERT: Destabilizacja
            </span>
          </div>

          {/* 3D CANVAS */}
          <div className="flex-1 w-full h-full min-h-0 relative flex items-center justify-center overflow-hidden">
            {render3DViewport(mutatedNodes, true, rotX2, rotY2)}

            {/* Overlaid Badges */}
            <div className="absolute top-3 left-3 bg-black/80 backdrop-blur-sm border border-red-900/60 rounded-lg p-2 text-[11px] space-y-1 pointer-events-none">
              <div className="flex justify-between gap-4">
                <span className="text-zinc-400">Energia Kohezyjna:</span>
                <strong className="text-red-400 font-extrabold">{selectedMutation.mutatedEnergy} kcal/mol</strong>
              </div>
              <div className="flex justify-between gap-4">
                <span className="text-zinc-400">Destabilizacja ΔΔG:</span>
                <strong className="text-red-400 font-extrabold">+{deltaEnergy.toFixed(1)} kcal/mol</strong>
              </div>
              <div className="flex justify-between gap-4">
                <span className="text-zinc-400">Stan Mostków S-S:</span>
                <strong className="text-red-300">{selectedMutation.brokenBridges[0]}</strong>
              </div>
            </div>

            <div className="absolute bottom-2 left-3 text-[9.5px] text-red-400/80 font-mono pointer-events-none">
              Czerwone linie i węzły wskazują strefę uszkodzenia mechanicznego szkieletu
            </div>
          </div>

          {/* Track verification mini footer */}
          <div className="p-2.5 bg-[#0b0306] border-t border-zinc-900 text-[10.5px] flex items-center justify-between text-zinc-400">
            <span className="flex items-center gap-1 text-red-400 font-bold">
              <XCircle className="w-3.5 h-3.5" />
              Track 5: Zerwanie wiązania kowalencyjnego
            </span>
            <span className="text-red-300 font-bold">
              Reszta #{selectedMutation.position}: {selectedMutation.mutatedResidue} (Zmutowana)
            </span>
          </div>
        </div>
      </div>

      {/* BOTTOM COMPARATIVE METRICS & 5 TRACKS TELEMETRY */}
      <div className="border-t border-zinc-900 bg-[#030308] p-3 shrink-0 flex flex-wrap items-center justify-between gap-3 text-xs font-mono">
        <div className="flex items-center gap-4 flex-wrap">
          <div className="flex items-center gap-1.5">
            <span className="text-zinc-500 uppercase font-bold text-[10px]">Weryfikacja 5 Ścieżek SNAPLOCK:</span>
          </div>

          <div className="flex items-center gap-1.5 px-2.5 py-1 rounded bg-black border border-zinc-800 text-[11px]">
            <span className="text-zinc-500">1. Szkielet (#248➔#1):</span>
            <span className="text-emerald-400 font-bold">Koaksjalny</span>
          </div>

          <div className="flex items-center gap-1.5 px-2.5 py-1 rounded bg-black border border-zinc-800 text-[11px]">
            <span className="text-zinc-500">2. VdW Exclusion:</span>
            <span className={selectedMutation.coreStatus === 'COMPACT' ? 'text-emerald-400 font-bold' : 'text-amber-400 font-bold'}>
              {selectedMutation.coreStatus === 'COMPACT' ? 'Zgodny (Brak kolizji)' : 'Deformacja steryczna'}
            </span>
          </div>

          <div className="flex items-center gap-1.5 px-2.5 py-1 rounded bg-black border border-zinc-800 text-[11px]">
            <span className="text-zinc-500">3. Elektrostatyka:</span>
            <span className="text-emerald-400 font-bold">Zbilansowana (pH 7.4)</span>
          </div>

          <div className="flex items-center gap-1.5 px-2.5 py-1 rounded bg-black border border-zinc-800 text-[11px]">
            <span className="text-zinc-500">4. Rdzeń Hydrofobowy:</span>
            <span className={selectedMutation.coreStatus === 'COMPACT' ? 'text-emerald-400 font-bold' : 'text-red-400 font-bold'}>
              {selectedMutation.coreStatus === 'COMPACT' ? 'Zwinięty' : 'Odsłonięty / Zdestabilizowany'}
            </span>
          </div>

          <div className="flex items-center gap-1.5 px-2.5 py-1 rounded bg-black border border-zinc-800 text-[11px]">
            <span className="text-zinc-500">5. Klamry S-S:</span>
            <span className={selectedMutation.track5Status !== 'BROKEN' ? 'text-emerald-400 font-bold' : 'text-red-400 font-extrabold'}>
              {selectedMutation.track5Status === 'LOCKED'
                ? `Locked (${targetSpec.explicitDisulfidePairs0Based.length} aktywne)`
                : selectedMutation.track5Status === 'NONE'
                ? 'Brak S-S (Rdzeń β-Grasp/α)'
                : 'BROKEN (Utrata mostka!)'}
            </span>
          </div>
        </div>

        {/* ZOOM & CONTROLS */}
        <div className="flex items-center gap-2">
          <span className="text-[10px] text-zinc-500 uppercase">Powiększenie 3D:</span>
          <input
            type="range"
            min="0.5"
            max="2.0"
            step="0.05"
            value={zoomScale}
            onChange={(e) => setZoomScale(parseFloat(e.target.value))}
            className="w-24 accent-red-500 cursor-pointer"
          />
          <span className="text-zinc-400 font-bold text-[10px]">{(zoomScale * 100).toFixed(0)}%</span>
          
          <button
            onClick={() => {
              setRotX1(-20);
              setRotY1(35);
              setRotX2(-20);
              setRotY2(35);
              setZoomScale(1.0);
            }}
            className="px-2 py-0.5 rounded bg-zinc-900 hover:bg-zinc-800 border border-zinc-700 text-zinc-300 hover:text-white text-[10px] cursor-pointer"
          >
            Reset Widoku
          </button>
        </div>
      </div>
    </div>
  );
};
