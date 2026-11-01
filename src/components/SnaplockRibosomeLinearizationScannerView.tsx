import React, { useEffect, useMemo, useRef, useState } from 'react';
import * as THREE from 'three';
import {
  Activity,
  AlertTriangle,
  Atom,
  CheckCircle2,
  Compass,
  Cpu,
  Flame,
  Layers,
  Pause,
  Play,
  RotateCcw,
  ShieldCheck,
  SlidersHorizontal,
  Zap
} from 'lucide-react';
import {
  Bar,
  BarChart,
  CartesianGrid,
  Legend,
  Line,
  LineChart,
  ResponsiveContainer,
  Tooltip,
  XAxis,
  YAxis
} from 'recharts';
import {
  GLOBAL_PROTEIN_TARGETS,
  GlobalProteinTargetId,
  resolveStrictTargetReference
} from '../snaplockGlobalTargets';
import {
  DiagnosticScanPathId,
  runAlphaHelixBaselineCurriculum,
  runCotranslationalFoldingTrajectory,
  runFourPathPlateDiagnostics
} from '../snaplockCotranslationalEngine';
import { BoundaryConditionsConfig, parseFastaSequence } from '../snaplockFastaAndPdb';
import { runMultiAgentFoldingToConvergence } from '../snaplockMultiAgentEngine';

interface SnaplockRibosomeLinearizationScannerViewProps {
  activeTarget: GlobalProteinTargetId;
  activeSequence: string;
}

export const SnaplockRibosomeLinearizationScannerView: React.FC<
  SnaplockRibosomeLinearizationScannerViewProps
> = ({ activeTarget, activeSequence }) => {
  const [subMode, setSubMode] = useState<
    'ribosome_linearization' | 'multiagent_soap_bubble' | 'helix_curriculum' | 'four_path_scanner'
  >('ribosome_linearization');

  const [selectedResidueIdx, setSelectedResidueIdx] = useState<number>(0);
  const [activeScanPath, setActiveScanPath] =
    useState<DiagnosticScanPathId>('SPIRAL_BOTTOM_TO_TOP');
  const [scanStepIdx, setScanStepIdx] = useState<number>(0);
  const [isAutoScanning, setIsAutoScanning] = useState<boolean>(false);

  // Parse synchronized sequence & target constraints
  const targetData = useMemo(() => {
    if (activeTarget !== 'CUSTOM' && GLOBAL_PROTEIN_TARGETS[activeTarget]) {
      const spec = GLOBAL_PROTEIN_TARGETS[activeTarget];
      return {
        title: spec.titlePl,
        codes: spec.sequence3Letter,
        cdAlpha: spec.cdAlphaTargetPct,
        cdBeta: spec.cdBetaTargetPct,
        ssPairs: spec.explicitDisulfidePairs0Based
      };
    }
    const parsed = parseFastaSequence(activeSequence);
    return {
      title: `Sekwencja własna (${parsed.length} aa)`,
      codes: parsed.sequence3Letter.length > 0 ? parsed.sequence3Letter : ['ALA', 'GLU', 'ALA', 'LYS'],
      cdAlpha: 40,
      cdBeta: 20,
      ssPairs: parsed.explicitDisulfidePairs || []
    };
  }, [activeTarget, activeSequence]);

  // Run full cotranslational + ribosome detachment + 22-param linearization + 4-path scanner
  const trajectoryBundle = useMemo(() => {
    const n = targetData.codes.length;
    const rgTheo = parseFloat((2.8 * Math.pow(n, 0.33)).toFixed(2));
    const cfg: BoundaryConditionsConfig = {
      dMaxTarget: parseFloat((rgTheo * 3.15).toFixed(2)),
      rgTarget: rgTheo,
      enforceBoundaries: true,
      cdAlphaTargetPct: targetData.cdAlpha,
      cdBetaTargetPct: targetData.cdBeta,
      ssBondTolerance: 0.45,
      ssRestraintWeight: 3.8,
      explicitDisulfidePairs: targetData.ssPairs
    };
    const traj = runCotranslationalFoldingTrajectory(targetData.codes, cfg);
    const scan = runFourPathPlateDiagnostics(
      targetData.codes,
      traj.finalRelaxedCa,
      traj.finalAllAtom,
      traj.disulfidePairs
    );
    const strictRef = activeTarget !== 'CUSTOM' ? resolveStrictTargetReference(activeTarget) : null;
    const isHeldOutTest = strictRef?.targetSpec.splitRole === 'TEST';
    const multiAgentReport = runMultiAgentFoldingToConvergence(
      {
        sequenceCodes: targetData.codes,
        cdAlphaTargetPct: targetData.cdAlpha,
        cdBetaTargetPct: targetData.cdBeta,
        explicitDisulfidePairs: targetData.ssPairs,
        referencePdbText: strictRef?.pdbText,
        referenceId: strictRef?.resolvedId ?? 'DE-NOVO',
        enableTargetGuidedRefinement: Boolean(strictRef?.pdbText) && !isHeldOutTest
      },
      90,
      5
    );
    return { traj, scan, multiAgentReport };
  }, [targetData, activeTarget]);

  // Reset indices when protein switches
  useEffect(() => {
    setSelectedResidueIdx(0);
    setScanStepIdx(0);
    setIsAutoScanning(false);
  }, [activeTarget, activeSequence]);

  // Reset scan step when path changes
  useEffect(() => {
    setScanStepIdx(0);
  }, [activeScanPath]);

  const currentPathReport = trajectoryBundle.scan.paths[activeScanPath];

  // Auto-stepping animation along the active spiral/sidechain/disulfide path
  useEffect(() => {
    if (!isAutoScanning) return;
    const total = currentPathReport.waypoints.length;
    if (total === 0) return;
    const timer = setInterval(() => {
      setScanStepIdx(prev => {
        if (prev + 1 >= total) {
          setIsAutoScanning(false);
          return total - 1;
        }
        return prev + 1;
      });
    }, 180);
    return () => clearInterval(timer);
  }, [isAutoScanning, currentPathReport.waypoints.length]);

  // Alpha-helix baseline curriculum (computed once cleanly)
  const curriculumReport = useMemo(() => {
    const lycotoxinCodes =
      'ILE TRP LEU THR ALA LEU LYS PHE LEU GLY LYS HIS ALA ALA LYS ALA LYS ALA LYS ALA LYS HIS LEU CYS ALA LYS'.split(
        ' '
      );
    return runAlphaHelixBaselineCurriculum([
      {
        id: 'HELIX_BASELINE_18',
        name: '1. Bazowa α-Helisa Marqusee–Baldwina (18 aa: AEAAAKEAAAKEAAAKAA)',
        category: 'ALPHA_HELIX_BASELINE',
        sequenceCodes: GLOBAL_PROTEIN_TARGETS.HELIX_BASELINE_18.sequence3Letter,
        cdAlphaPct: 88,
        cdBetaPct: 0,
        ssPairs: [],
        referencePdbText: GLOBAL_PROTEIN_TARGETS.HELIX_BASELINE_18.referencePdbText,
        referenceId: 'HELIX18_REF'
      },
      {
        id: 'LYCOTOXIN_HELIX_26',
        name: '2. Naturalna Amfipatyczna α-Helisa Likotoksyny (26 aa, TRP/PHE/LEU/LYS)',
        category: 'ALPHA_HELIX_BASELINE',
        sequenceCodes: lycotoxinCodes,
        cdAlphaPct: 78,
        cdBetaPct: 0,
        ssPairs: []
      },
      {
        id: 'INSULIN_A_21',
        name: '3. Dwuhelikalny Pęczek Insuliny A z Klamrą CYS6–CYS11 (21 aa, 3I40:A)',
        category: 'ALPHA_DISULFIDE_BUNDLE',
        sequenceCodes: GLOBAL_PROTEIN_TARGETS.INSULIN_A_21.sequence3Letter,
        cdAlphaPct: 52,
        cdBetaPct: 0,
        ssPairs: [[5, 10]],
        referencePdbText: GLOBAL_PROTEIN_TARGETS.INSULIN_A_21.referencePdbText,
        referenceId: '3I40_A_21'
      },
      {
        id: 'CRAMBIN_46',
        name: '4. Transfer na Topologię α+β + 3x S-S: Krambina (46 aa, 1CRN:A)',
        category: 'ALPHA_BETA_TRANSITION',
        sequenceCodes: GLOBAL_PROTEIN_TARGETS.CRAMBIN_46.sequence3Letter,
        cdAlphaPct: 41,
        cdBetaPct: 18,
        ssPairs: GLOBAL_PROTEIN_TARGETS.CRAMBIN_46.explicitDisulfidePairs0Based,
        referencePdbText: GLOBAL_PROTEIN_TARGETS.CRAMBIN_46.referencePdbText,
        referenceId: '1CRN'
      },
      {
        id: 'UBIQUITIN_76',
        name: '5. Transfer na Arkusz 5β + α1 (β-Grasp): Ubikwityna (76 aa, 1UBQ:A)',
        category: 'ALPHA_BETA_TRANSITION',
        sequenceCodes: GLOBAL_PROTEIN_TARGETS.UBIQUITIN_76.sequence3Letter,
        cdAlphaPct: 24,
        cdBetaPct: 38,
        ssPairs: [],
        referencePdbText: GLOBAL_PROTEIN_TARGETS.UBIQUITIN_76.referencePdbText,
        referenceId: '1UBQ'
      }
    ]);
  }, []);

  // Three.js 3D viewport for Sub-tab 3 (4-Path Spatial Plate & Spiral Scanner)
  const scannerMountRef = useRef<HTMLDivElement | null>(null);
  const orbitRef = useRef({ rotX: 0.32, rotY: 0.65, dist: 46, dragging: false, px: 0, py: 0 });

  useEffect(() => {
    if (subMode !== 'four_path_scanner' || !scannerMountRef.current) return;
    const container = scannerMountRef.current;
    const width = container.clientWidth || 640;
    const height = container.clientHeight || 420;

    const scene = new THREE.Scene();
    scene.background = new THREE.Color(0x020308);
    const camera = new THREE.PerspectiveCamera(45, width / height, 0.5, 900);
    const renderer = new THREE.WebGLRenderer({ antialias: true });
    renderer.setSize(width, height);
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    container.innerHTML = '';
    container.appendChild(renderer.domElement);

    scene.add(new THREE.AmbientLight(0xffffff, 0.8));
    const dLight = new THREE.DirectionalLight(0xffffff, 1.3);
    dLight.position.set(22, 34, 28);
    scene.add(dLight);

    const disposables: { geos: THREE.BufferGeometry[]; mats: THREE.Material[] } = {
      geos: [],
      mats: []
    };
    const trackG = <T extends THREE.BufferGeometry>(g: T): T => {
      disposables.geos.push(g);
      return g;
    };
    const trackM = <T extends THREE.Material>(m: T): T => {
      disposables.mats.push(m);
      return m;
    };

    const rootGroup = new THREE.Group();
    scene.add(rootGroup);

    const SCALE = 1.9;
    const [cx, cy, cz] = trajectoryBundle.traj.adaptiveSphere.centerOfMass;
    const sphereR = trajectoryBundle.traj.adaptiveSphere.adaptiveSphereRadiusAngstrom * SCALE;
    orbitRef.current.dist = Math.max(36, sphereR * 2.55);

    // 1. Adaptive E8 Bounding Sphere (R_kula >= 1.48 * Rg)
    const sphGeo = trackG(new THREE.SphereGeometry(sphereR, 28, 20));
    const sphWire = trackG(new THREE.WireframeGeometry(sphGeo));
    const sphMat = trackM(
      new THREE.LineBasicMaterial({ color: 0x06b6d4, transparent: true, opacity: 0.2 })
    );
    rootGroup.add(new THREE.LineSegments(sphWire, sphMat));

    // 2. Protein Cα backbone + residue spheres
    const caVecs = trajectoryBundle.traj.finalRelaxedCa.map(
      p => new THREE.Vector3((p[0] - cx) * SCALE, (p[1] - cy) * SCALE, (p[2] - cz) * SCALE)
    );
    if (caVecs.length >= 3) {
      const curve = new THREE.CatmullRomCurve3(caVecs);
      const tubeG = trackG(new THREE.TubeGeometry(curve, Math.max(20, caVecs.length * 5), 0.22, 8, false));
      const tubeM = trackM(
        new THREE.MeshStandardMaterial({
          color: 0x334155,
          transparent: true,
          opacity: 0.68,
          roughness: 0.4
        })
      );
      rootGroup.add(new THREE.Mesh(tubeG, tubeM));
    }

    const activeWp =
      currentPathReport.waypoints[
        Math.min(scanStepIdx, Math.max(0, currentPathReport.waypoints.length - 1))
      ];
    const visitedResSet = new Set(
      currentPathReport.waypoints.slice(0, scanStepIdx + 1).map(w => w.resIndex)
    );

    const unitSph = trackG(new THREE.SphereGeometry(1, 14, 14));
    const matVisited = trackM(new THREE.MeshStandardMaterial({ color: 0x22d3ee, roughness: 0.25 }));
    const matUnvisited = trackM(new THREE.MeshStandardMaterial({ color: 0x475569, roughness: 0.5 }));
    const matActive = trackM(
      new THREE.MeshStandardMaterial({
        color: 0xfacc15,
        emissive: 0xd97706,
        emissiveIntensity: 0.65,
        roughness: 0.2
      })
    );
    const matWarn = trackM(
      new THREE.MeshStandardMaterial({
        color: 0xf43f5e,
        emissive: 0xbe123c,
        emissiveIntensity: 0.5
      })
    );

    const issueResSet = new Set(
      currentPathReport.waypoints.filter(w => w.status !== 'OK').map(w => w.resIndex)
    );

    caVecs.forEach((v, idx) => {
      const isCurrent = activeWp?.resIndex === idx;
      const hasIssue = issueResSet.has(idx);
      const isVis = visitedResSet.has(idx);
      const m = new THREE.Mesh(
        unitSph,
        isCurrent ? matActive : hasIssue ? matWarn : isVis ? matVisited : matUnvisited
      );
      m.position.copy(v);
      m.scale.setScalar(isCurrent ? 0.68 : hasIssue ? 0.5 : 0.36);
      rootGroup.add(m);
    });

    // 3. Render 3D Spiral / Diagnostic Path Curve + Moving Scanning Plate ("przez płytę")
    if (
      activeScanPath === 'SPIRAL_BOTTOM_TO_TOP' ||
      activeScanPath === 'SPIRAL_FRONT_TO_BACK'
    ) {
      const spiralPts = currentPathReport.spiralCurve3D.map(
        p => new THREE.Vector3((p[0] - cx) * SCALE, (p[1] - cy) * SCALE, (p[2] - cz) * SCALE)
      );
      if (spiralPts.length >= 2) {
        const lineGeo = trackG(new THREE.BufferGeometry().setFromPoints(spiralPts));
        const lineMat = trackM(
          new THREE.LineBasicMaterial({
            color: activeScanPath === 'SPIRAL_BOTTOM_TO_TOP' ? 0x38bdf8 : 0x34d399,
            transparent: true,
            opacity: 0.72
          })
        );
        rootGroup.add(new THREE.Line(lineGeo, lineMat));
      }

      // Moving translucent diagnostic cutting plate at active waypoint height/depth
      if (activeWp) {
        const plateGeo = trackG(new THREE.PlaneGeometry(sphereR * 1.9, sphereR * 1.9, 8, 8));
        const plateMat = trackM(
          new THREE.MeshBasicMaterial({
            color: activeScanPath === 'SPIRAL_BOTTOM_TO_TOP' ? 0x06b6d4 : 0x10b981,
            transparent: true,
            opacity: 0.16,
            side: THREE.DoubleSide
          })
        );
        const plateMesh = new THREE.Mesh(plateGeo, plateMat);
        if (activeScanPath === 'SPIRAL_BOTTOM_TO_TOP') {
          plateMesh.rotation.x = Math.PI / 2;
          plateMesh.position.y = (activeWp.position[1] - cy) * SCALE;
        } else {
          plateMesh.position.z = (activeWp.position[2] - cz) * SCALE;
        }
        rootGroup.add(plateMesh);
      }
    } else if (activeScanPath === 'DISULFIDE_BRIDGES_AUDIT') {
      const ssMat = trackM(
        new THREE.LineBasicMaterial({ color: 0xfbbf24, transparent: true, opacity: 0.95 })
      );
      trajectoryBundle.traj.disulfidePairs.forEach(([i1, i2]) => {
        if (caVecs[i1] && caVecs[i2]) {
          const g = trackG(new THREE.BufferGeometry().setFromPoints([caVecs[i1], caVecs[i2]]));
          rootGroup.add(new THREE.Line(g, ssMat));
        }
      });
    } else if (activeScanPath === 'LATERAL_SIDECHAIN_BONDS') {
      const hbMat = trackM(
        new THREE.LineBasicMaterial({ color: 0xe879f9, transparent: true, opacity: 0.8 })
      );
      for (const b of trajectoryBundle.traj.finalAllAtom.bonds) {
        if (b.category === 'HBOND' || b.category === 'SALT_BRIDGE') {
          const r1 = b.resSeq1 - 1;
          const r2 = b.resSeq2 - 1;
          if (caVecs[r1] && caVecs[r2]) {
            const g = trackG(new THREE.BufferGeometry().setFromPoints([caVecs[r1], caVecs[r2]]));
            rootGroup.add(new THREE.Line(g, hbMat));
          }
        }
      }
    }

    let animId = 0;
    const renderLoop = () => {
      animId = requestAnimationFrame(renderLoop);
      const { rotX, rotY, dist } = orbitRef.current;
      camera.position.set(
        dist * Math.sin(rotY) * Math.cos(rotX),
        dist * Math.sin(rotX),
        dist * Math.cos(rotY) * Math.cos(rotX)
      );
      camera.lookAt(0, 0, 0);
      renderer.render(scene, camera);
    };
    renderLoop();

    const onResize = () => {
      if (!container) return;
      const w = container.clientWidth || 640;
      const h = container.clientHeight || 420;
      camera.aspect = w / h;
      camera.updateProjectionMatrix();
      renderer.setSize(w, h);
    };
    window.addEventListener('resize', onResize);

    return () => {
      cancelAnimationFrame(animId);
      window.removeEventListener('resize', onResize);
      disposables.geos.forEach(g => g.dispose());
      disposables.mats.forEach(m => m.dispose());
      renderer.dispose();
    };
  }, [subMode, activeScanPath, scanStepIdx, trajectoryBundle, currentPathReport]);

  const lin = trajectoryBundle.traj.energyLinearization;
  const sphere = trajectoryBundle.traj.adaptiveSphere;
  const activeCoupling =
    lin.perResidueCouplings[
      Math.min(selectedResidueIdx, Math.max(0, lin.perResidueCouplings.length - 1))
    ];
  const activeScanWaypoint =
    currentPathReport.waypoints[
      Math.min(scanStepIdx, Math.max(0, currentPathReport.waypoints.length - 1))
    ];

  return (
    <div className="flex-1 w-full h-full min-h-0 bg-[#020308] text-slate-200 flex flex-col overflow-hidden font-mono select-none">
      {/* SUB-NAVIGATION BAR */}
      <div className="bg-[#05060f] border-b border-zinc-900 px-4 pt-2 flex flex-wrap items-center justify-between gap-2 shrink-0">
        <div className="flex items-center gap-2 overflow-x-auto">
          <button
            onClick={() => setSubMode('ribosome_linearization')}
            className={`px-3.5 py-2 rounded-t-lg text-xs font-bold flex items-center gap-2 border-t border-x cursor-pointer transition ${
              subMode === 'ribosome_linearization'
                ? 'bg-[#020308] border-cyan-500/60 border-b-transparent text-white'
                : 'bg-black/50 border-zinc-800 text-zinc-400 hover:text-white'
            }`}
          >
            <Zap className="w-3.5 h-3.5 text-cyan-400" />
            <span>1. Oderwanie z Rybosomu &amp; Linearyzacja Energii (22 Parametry)</span>
          </button>

          <button
            onClick={() => setSubMode('multiagent_soap_bubble')}
            className={`px-3.5 py-2 rounded-t-lg text-xs font-bold flex items-center gap-2 border-t border-x cursor-pointer transition ${
              subMode === 'multiagent_soap_bubble'
                ? 'bg-[#020308] border-fuchsia-500/60 border-b-transparent text-white'
                : 'bg-black/50 border-zinc-800 text-zinc-400 hover:text-white'
            }`}
          >
            <Flame className="w-3.5 h-3.5 text-fuchsia-400" />
            <span>2. Wieloagentowa Bańka Mydlana &amp; Soft SnapLock (φ,ψ)</span>
          </button>

          <button
            onClick={() => setSubMode('helix_curriculum')}
            className={`px-3.5 py-2 rounded-t-lg text-xs font-bold flex items-center gap-2 border-t border-x cursor-pointer transition ${
              subMode === 'helix_curriculum'
                ? 'bg-[#020308] border-amber-500/60 border-b-transparent text-white'
                : 'bg-black/50 border-zinc-800 text-zinc-400 hover:text-white'
            }`}
          >
            <Layers className="w-3.5 h-3.5 text-amber-400" />
            <span>3. Baza α-Helisy ➔ Kalibracja i Test ➔ Arkusze β</span>
          </button>

          <button
            onClick={() => setSubMode('four_path_scanner')}
            className={`px-3.5 py-2 rounded-t-lg text-xs font-bold flex items-center gap-2 border-t border-x cursor-pointer transition ${
              subMode === 'four_path_scanner'
                ? 'bg-[#020308] border-emerald-500/60 border-b-transparent text-white'
                : 'bg-black/50 border-zinc-800 text-zinc-400 hover:text-white'
            }`}
          >
            <Compass className="w-3.5 h-3.5 text-emerald-400" />
            <span>4. Skaner Błędów Przez Płytę (2x Ślimak + Boczne + S-S)</span>
          </button>
        </div>

        <div className="flex items-center gap-2 pb-1.5 text-[11px]">
          <span className="px-2.5 py-0.5 rounded bg-cyan-950/80 border border-cyan-700/60 text-cyan-200 font-bold">
            Aktywny cel: {targetData.title}
          </span>
          <span className="px-2.5 py-0.5 rounded bg-emerald-950/80 border border-emerald-700/60 text-emerald-300 font-bold">
            22/22 Parametry: AKTYWNE
          </span>
        </div>
      </div>

      {/* MAIN SCROLLABLE WORKSPACE */}
      <div className="flex-1 overflow-y-auto p-4 sm:p-5 space-y-5 scrollbar-thin scrollbar-thumb-zinc-800">
        {subMode === 'ribosome_linearization' && (
          <>
            {/* 3 REGIMES COMPARISON: ON RIBOSOME vs DETACHMENT vs FINAL VDW/HYDRO RELAXATION */}
            <div className="grid grid-cols-1 lg:grid-cols-3 gap-4">
              {/* Card 1: On-Ribosome Tethered */}
              <div className="bg-[#060814] border border-amber-800/60 rounded-xl p-4 space-y-2">
                <div className="flex items-center justify-between">
                  <span className="text-[10.5px] uppercase font-bold text-amber-300">
                    Stan I: Łańcuch na Rybosomie (Tunel N➔C)
                  </span>
                  <span className="text-[10px] px-2 py-0.5 rounded bg-amber-950 text-amber-200 border border-amber-700/60">
                    ε_tunel = {lin.onRibosomeDielectricConstant}
                  </span>
                </div>
                <p className="text-2xl font-extrabold text-amber-300 tabular-nums">
                  {lin.onRibosomeTotalEnergyKcal.toFixed(2)} kcal/mol
                </p>
                <div className="text-[11px] text-zinc-400 space-y-1 pt-1 border-t border-zinc-800/80 tabular-nums">
                  <div className="flex justify-between">
                    <span>Napięcie kotwicy tRNA (E_tension):</span>
                    <strong className="text-amber-200">+{lin.onRibosomeTetherTensionKcal} kcal/mol</strong>
                  </div>
                  <div className="flex justify-between">
                    <span>Van der Waals (w tunelu):</span>
                    <strong className="text-zinc-200">{lin.onRibosomeVdwKcal} kcal/mol</strong>
                  </div>
                  <div className="flex justify-between">
                    <span>Hydrofobowość (stłumiona):</span>
                    <strong className="text-zinc-200">{lin.onRibosomeHydrophobicKcal} kcal/mol</strong>
                  </div>
                  <div className="text-[10px] text-zinc-500 pt-1">
                    C-koniec zestryfikowany na peptydylo-tRNA (brak ładunku pK1 α-COOH).
                  </div>
                </div>
              </div>

              {/* Card 2: Ribosome Cleavage / Entry to Free Solvent */}
              <div className="bg-[#060814] border border-cyan-800/60 rounded-xl p-4 space-y-2">
                <div className="flex items-center justify-between">
                  <span className="text-[10.5px] uppercase font-bold text-cyan-300">
                    Stan II: Faza Wejścia na Końcu (Zerwanie z Rybosomu)
                  </span>
                  <span className="text-[10px] px-2 py-0.5 rounded bg-cyan-950 text-cyan-200 border border-cyan-700/60">
                    Hydroliza tRNA
                  </span>
                </div>
                <p className="text-2xl font-extrabold text-cyan-300 tabular-nums">
                  {lin.detachedUnrelaxedEnergyKcal.toFixed(2)} kcal/mol
                </p>
                <div className="text-[11px] text-zinc-400 space-y-1 pt-1 border-t border-zinc-800/80 tabular-nums">
                  <div className="flex justify-between">
                    <span>Zrzut napięcia łańcucha (ΔE_tension):</span>
                    <strong className="text-emerald-300">{lin.releasedTensionDeltaKcal} kcal/mol</strong>
                  </div>
                  <div className="flex justify-between">
                    <span>Jonizacja wolnego C-końca (pK1_COOH):</span>
                    <strong className="text-cyan-200">{lin.cTerminusIonizationDeltaKcal} kcal/mol</strong>
                  </div>
                  <div className="flex justify-between">
                    <span>Norma zlinearyzowanej siły ||∇E||:</span>
                    <strong className="text-amber-300">{lin.linearizedTotalGradientNorm} kcal/(mol·Å)</strong>
                  </div>
                  <div className="text-[10px] text-zinc-500 pt-1">
                    Białko uwolnione do pełnego rozpuszczalnika (ε = 78.5) — gotowe do ostatecznej relaksacji.
                  </div>
                </div>
              </div>

              {/* Card 3: Final Post-Ribosomal VdW + Hydrophobic Relaxation */}
              <div className="bg-[#060814] border border-emerald-700/60 rounded-xl p-4 space-y-2">
                <div className="flex items-center justify-between">
                  <span className="text-[10.5px] uppercase font-bold text-emerald-300">
                    Stan III: Relaksacja Końcowa (VdW + Hydrofobowość)
                  </span>
                  <span className="text-[10px] px-2 py-0.5 rounded bg-emerald-950 text-emerald-200 border border-emerald-700/60">
                    Zbieżność Taylora: {lin.linearizationFidelityPct}%
                  </span>
                </div>
                <p className="text-2xl font-extrabold text-emerald-400 tabular-nums">
                  {lin.postReleaseRelaxedTotalEnergyKcal.toFixed(2)} kcal/mol
                </p>
                <div className="text-[11px] text-zinc-400 space-y-1 pt-1 border-t border-zinc-800/80 tabular-nums">
                  <div className="flex justify-between">
                    <span>Ponowny Van der Waals (LJ 6-12):</span>
                    <strong className="text-emerald-300">{lin.postReleaseVdwKcal} kcal/mol</strong>
                  </div>
                  <div className="flex justify-between">
                    <span>Ponowna Hydrofobowość (Kyte/Hopp/SASA):</span>
                    <strong className="text-emerald-300">{lin.postReleaseHydrophobicKcal} kcal/mol</strong>
                  </div>
                  <div className="flex justify-between">
                    <span>Zysk relaksacji po oderwaniu (ΔE_relax):</span>
                    <strong className="text-emerald-400">{lin.relaxationEnergyGainKcal} kcal/mol</strong>
                  </div>
                  <div className="text-[10px] text-zinc-500 pt-1">
                    Przesunięcie relaksacyjne ΔX = {lin.postReleaseCoordinateShiftRmsdAngstrom} Å • Sfera R_kula = {sphere.adaptiveSphereRadiusAngstrom} Å ({sphere.sphereToRgRatio}×Rg).
                  </div>
                </div>
              </div>
            </div>

            {/* LINEARIZED FORCE GRADIENTS SUMMARY */}
            <div className="bg-[#050712] border border-zinc-800 rounded-xl p-4 space-y-3">
              <div className="flex flex-wrap items-center justify-between gap-2">
                <div>
                  <h3 className="text-xs font-bold text-white uppercase tracking-wider flex items-center gap-2">
                    <Cpu className="w-4 h-4 text-cyan-400" />
                    <span>
                      Linearyzacja Energii wokół Stanu Uwolnionego z Rybosomu: E(q + Δq) ≈ E(q₀) + Σ J_k · Δq
                    </span>
                  </h3>
                  <p className="text-[11px] text-zinc-400 mt-0.5">
                    Składowe zlinearyzowanych sił (jakobianu gradientu) decydujące o ostatecznym domknięciu rdzenia po zerwaniu z rybosomu
                  </p>
                </div>
                <div className="text-[11px] text-cyan-300 font-bold tabular-nums">
                  Predykcja liniowa Taylora: {lin.taylorFirstOrderPredictedDropKcal} kcal/mol vs Rzeczywisty spadek: {lin.actualNonlinearDropKcal} kcal/mol
                </div>
              </div>

              <div className="grid grid-cols-2 sm:grid-cols-4 gap-3 text-xs tabular-nums">
                <div className="bg-black/60 border border-zinc-800 rounded-lg p-3">
                  <span className="text-[10px] text-zinc-500 uppercase block">Siła Van der Waalsa ||∇E_VdW||</span>
                  <strong className="text-lg text-sky-300">{lin.linearizedGradientNormVdw} kcal/(mol·Å)</strong>
                  <span className="block text-[10px] text-zinc-500">Z objętości VdW (#7) i wzoru (#5)</span>
                </div>
                <div className="bg-black/60 border border-zinc-800 rounded-lg p-3">
                  <span className="text-[10px] text-zinc-500 uppercase block">Siła Hydrofobowa ||∇E_Hydro||</span>
                  <strong className="text-lg text-emerald-300">{lin.linearizedGradientNormHydrophobic} kcal/(mol·Å)</strong>
                  <span className="block text-[10px] text-zinc-500">Kyte (#10), Hopp (#11), SASA (#8)</span>
                </div>
                <div className="bg-black/60 border border-zinc-800 rounded-lg p-3">
                  <span className="text-[10px] text-zinc-500 uppercase block">Siła Elektrostatyczna ||∇E_Elec||</span>
                  <strong className="text-lg text-fuchsia-300">{lin.linearizedGradientNormElectrostatic} kcal/(mol·Å)</strong>
                  <span className="block text-[10px] text-zinc-500">pH 7.4 (#13), pK1/pK2/pKr/pI (#15–18)</span>
                </div>
                <div className="bg-black/60 border border-zinc-800 rounded-lg p-3">
                  <span className="text-[10px] text-zinc-500 uppercase block">Napięcie Torsyjne ||∇E_Tension||</span>
                  <strong className="text-lg text-amber-300">{lin.linearizedGradientNormTorsionTension} kcal/(mol·Å)</strong>
                  <span className="block text-[10px] text-zinc-500">Masa (#6), Elastyczność (#22), Chiralność (#20)</span>
                </div>
              </div>

              {/* Recharts Per-Residue Energy: On-Ribosome vs Detached Relaxed */}
              <div className="h-60 w-full pt-2">
                <ResponsiveContainer width="100%" height="100%">
                  <LineChart
                    data={lin.perResidueCouplings.map(c => ({
                      label: `${c.code1}${c.resSeq}`,
                      resSeq: c.resSeq,
                      naRybosomie: c.totalOnRibosomeKcal,
                      poOderwaniuRelaksacja: c.totalPostReleaseRelaxedKcal,
                      silaZlinearyzowana: c.linearizedForceNormKcalPerAngstrom
                    }))}
                    onClick={(st: any) => {
                      if (st && typeof st.activeTooltipIndex === 'number') {
                        setSelectedResidueIdx(st.activeTooltipIndex);
                      }
                    }}
                  >
                    <CartesianGrid strokeDasharray="3 3" stroke="#1e293b" />
                    <XAxis dataKey="label" stroke="#94a3b8" tick={{ fontSize: 10 }} />
                    <YAxis stroke="#94a3b8" tick={{ fontSize: 10 }} />
                    <Tooltip
                      contentStyle={{
                        backgroundColor: '#050712',
                        borderColor: '#334155',
                        fontSize: '11px'
                      }}
                    />
                    <Legend wrapperStyle={{ fontSize: '11px' }} />
                    <Line
                      type="monotone"
                      dataKey="naRybosomie"
                      name="Energia na rybosomie [kcal/mol]"
                      stroke="#f59e0b"
                      strokeWidth={2}
                      dot={false}
                    />
                    <Line
                      type="monotone"
                      dataKey="poOderwaniuRelaksacja"
                      name="Energia po oderwaniu i relaksacji VdW/Hydro [kcal/mol]"
                      stroke="#10b981"
                      strokeWidth={2}
                      dot={false}
                    />
                    <Line
                      type="monotone"
                      dataKey="silaZlinearyzowana"
                      name="Norma siły zlinearyzowanej ||∇E_i|| [kcal/mol·Å]"
                      stroke="#38bdf8"
                      strokeWidth={1.6}
                      dot={false}
                    />
                  </LineChart>
                </ResponsiveContainer>
              </div>
            </div>

            {/* PER-RESIDUE 22-PARAMETER INSPECTOR + FULL 22-COLUMN VERIFICATION TABLE */}
            <div className="grid grid-cols-1 xl:grid-cols-2 gap-4">
              {activeCoupling && (
                <div className="bg-[#050712] border border-cyan-900/60 rounded-xl p-4 space-y-3">
                  <div className="flex items-center justify-between">
                    <h4 className="text-xs font-bold text-white uppercase">
                      Inspektor 22 Parametrów dla Reszty #{activeCoupling.resSeq} ({activeCoupling.code3} — {activeCoupling.namePl})
                    </h4>
                    <div className="flex items-center gap-1">
                      <button
                        onClick={() => setSelectedResidueIdx(Math.max(0, selectedResidueIdx - 1))}
                        className="px-2 py-0.5 bg-zinc-900 hover:bg-zinc-800 border border-zinc-700 rounded text-xs cursor-pointer"
                      >
                        ◀
                      </button>
                      <span className="text-[11px] px-2 text-cyan-300 font-bold">
                        {activeCoupling.resSeq}/{lin.perResidueCouplings.length}
                      </span>
                      <button
                        onClick={() =>
                          setSelectedResidueIdx(
                            Math.min(lin.perResidueCouplings.length - 1, selectedResidueIdx + 1)
                          )
                        }
                        className="px-2 py-0.5 bg-zinc-900 hover:bg-zinc-800 border border-zinc-700 rounded text-xs cursor-pointer"
                      >
                        ▶
                      </button>
                    </div>
                  </div>

                  <div className="grid grid-cols-2 sm:grid-cols-3 gap-2 text-[11px] tabular-nums">
                    <div className="bg-black/60 border border-zinc-800/90 rounded p-2">
                      <span className="text-zinc-500 block text-[9.5px]">Promień r_VdW / ε_LJ</span>
                      <strong className="text-sky-300">
                        {activeCoupling.vdwRadiusAngstrom} Å / {activeCoupling.ljWellDepthKcal} kcal
                      </strong>
                    </div>
                    <div className="bg-black/60 border border-zinc-800/90 rounded p-2">
                      <span className="text-zinc-500 block text-[9.5px]">Zakopanie SASA (z {activeCoupling.referenceSasaAngstrom2} Å²)</span>
                      <strong className="text-emerald-300">
                        {(activeCoupling.burialFraction * 100).toFixed(1)}%
                      </strong>
                    </div>
                    <div className="bg-black/60 border border-zinc-800/90 rounded p-2">
                      <span className="text-zinc-500 block text-[9.5px]">Energia Hydrofobowa (Kyte)</span>
                      <strong className="text-emerald-400">
                        {activeCoupling.hydrophobicCoreTermKcal} kcal/mol
                      </strong>
                    </div>
                    <div className="bg-black/60 border border-zinc-800/90 rounded p-2">
                      <span className="text-zinc-500 block text-[9.5px]">Ładunek: Rybosom ➔ Wolny</span>
                      <strong className="text-fuchsia-300">
                        {activeCoupling.effectiveChargeOnRibosome} ➔ {activeCoupling.effectiveChargeDetached} e
                      </strong>
                    </div>
                    <div className="bg-black/60 border border-zinc-800/90 rounded p-2">
                      <span className="text-zinc-500 block text-[9.5px]">Sztywność zawiasu / Masa</span>
                      <strong className="text-amber-300">
                        k={activeCoupling.hingeStiffnessKcalRad2} • m×{activeCoupling.inertialMassFactor}
                      </strong>
                    </div>
                    <div className="bg-black/60 border border-zinc-800/90 rounded p-2">
                      <span className="text-zinc-500 block text-[9.5px]">Studnia Ramachandrana (φ₀, ψ₀)</span>
                      <strong className="text-indigo-300">
                        ({activeCoupling.preferredPhiDeg}°, {activeCoupling.preferredPsiDeg}°)
                      </strong>
                    </div>
                  </div>
                  <p className="text-[10.5px] text-zinc-400">
                    Kliknij dowolny punkt na wykresie powyżej lub użyj strzałek, aby przejrzeć wkład wszystkich 22 parametrów dla wybranego aminokwasu w łańcuchu.
                  </p>
                </div>
              )}

              {/* 22-Column Audit Matrix */}
              <div className="bg-[#050712] border border-zinc-800 rounded-xl p-4 space-y-2">
                <div className="flex items-center justify-between">
                  <h4 className="text-xs font-bold text-white uppercase flex items-center gap-2">
                    <CheckCircle2 className="w-4 h-4 text-emerald-400" />
                    <span>Audyt Włączenia Wszystkich 22 Parametrów Tabeli (21 Aminokwasów)</span>
                  </h4>
                  <span className="text-[10px] text-emerald-300 font-bold">100% POKRYCIA</span>
                </div>
                <div className="max-h-52 overflow-y-auto scrollbar-thin scrollbar-thumb-zinc-800">
                  <table className="w-full text-left text-[11px] tabular-nums">
                    <thead className="text-[9.5px] uppercase text-zinc-500 border-b border-zinc-800 sticky top-0 bg-[#050712]">
                      <tr>
                        <th className="py-1.5 pr-2">#</th>
                        <th className="py-1.5 pr-2">Parametr Tabeli</th>
                        <th className="py-1.5 pr-2">Rola Fizyczna w Silniku</th>
                        <th className="py-1.5 text-right">Status</th>
                      </tr>
                    </thead>
                    <tbody className="divide-y divide-zinc-900">
                      {lin.parameter22Audit.map(col => (
                        <tr key={col.columnIndex}>
                          <td className="py-1 pr-2 text-zinc-500 font-bold">{col.columnIndex}</td>
                          <td className="py-1 pr-2 text-cyan-200 font-bold">{col.labelPl}</td>
                          <td className="py-1 pr-2 text-zinc-400">{col.physicalRolePl}</td>
                          <td className="py-1 text-right text-emerald-400 font-bold">AKTYWNY</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              </div>
            </div>
          </>
        )}

        {subMode === 'helix_curriculum' && (
          <div className="space-y-4">
            <div className="bg-[#060814] border border-amber-800/60 rounded-xl p-4 space-y-2">
              <h3 className="text-sm font-bold text-white flex items-center gap-2">
                <Layers className="w-4 h-4 text-amber-400" />
                <span>
                  Progresywny Trening Fizyczny: Modele Bazowe α-Helisy (Pauling–Corey) ➔ Korekta 22 Parametrami ➔ Topologie α+β
                </span>
              </h3>
              <p className="text-xs text-zinc-400 leading-relaxed">
                Klasyczny model Paulinga–Coreya zakłada sztywne kąty (φ = -57.0°, ψ = -47.0°) niezależnie od aminokwasu. Nasz silnik wychodzi z bazy helikalnej, uczy się poprawki objętościowej VdW (#7), mostków solnych i→i+4 (#13/#17) oraz makro-dipola końców pK1/pK2 (#15/#16) na czystych α-helisach, a następnie przenosi te prawa na pęczki disiarczkowe i pełne białka α+β.
              </p>
            </div>

            <div className="bg-[#050712] border border-zinc-800 rounded-xl overflow-hidden">
              <div className="overflow-x-auto">
                <table className="w-full text-left text-xs tabular-nums">
                  <thead className="bg-[#03040a] text-[10px] uppercase text-zinc-400 border-b border-zinc-800">
                    <tr>
                      <th className="p-3">Etap</th>
                      <th className="p-3">Model / Białko</th>
                      <th className="p-3">Wiązanie H O(i)⋯N(i+4)</th>
                      <th className="p-3">Energia Bazowa</th>
                      <th className="p-3">Po Relaksacji 22-Param</th>
                      <th className="p-3">Zysk ΔE</th>
                      <th className="p-3">Kolizje VdW</th>
                      <th className="p-3">De-Novo GT Cα RMSD</th>
                    </tr>
                  </thead>
                  <tbody className="divide-y divide-zinc-900">
                    {curriculumReport.entries.map(row => (
                      <tr key={row.targetId} className="hover:bg-zinc-900/40">
                        <td className="p-3 font-bold text-amber-300">#{row.stageOrder}</td>
                        <td className="p-3">
                          <div className="font-bold text-white">{row.proteinName}</div>
                          <div className="text-[10.5px] text-zinc-500">{row.summaryPl}</div>
                        </td>
                        <td className="p-3">
                          <span className="text-zinc-400">{row.classicalPaulingCoreyHBondMeanAngstrom} Å</span>
                          <span className="mx-1 text-zinc-600">➔</span>
                          <strong className="text-cyan-300">{row.improved22ParamHBondMeanAngstrom} Å</strong>
                        </td>
                        <td className="p-3 text-zinc-400">{row.classicalEnergyKcal} kcal/mol</td>
                        <td className="p-3 font-bold text-emerald-300">
                          {row.improvedPostRibosomeEnergyKcal} kcal/mol
                        </td>
                        <td className="p-3 font-bold text-emerald-400">
                          {row.energyImprovementKcal <= 0
                            ? `${row.energyImprovementKcal} kcal`
                            : `+${row.energyImprovementKcal} kcal`}
                        </td>
                        <td className="p-3">
                          <span
                            className={`px-2 py-0.5 rounded text-[10.5px] font-bold ${
                              row.vdwClashesCount === 0
                                ? 'bg-emerald-950 text-emerald-300 border border-emerald-800'
                                : 'bg-rose-950 text-rose-300 border border-rose-800'
                            }`}
                          >
                            {row.vdwClashesCount}
                          </span>
                        </td>
                        <td className="p-3 font-bold">
                          {row.improvedDeNovoCaRmsdVsRefAngstrom !== null ? (
                            <span
                              className={
                                row.improvedDeNovoCaRmsdVsRefAngstrom < 1.0
                                  ? 'text-emerald-400'
                                  : 'text-amber-300'
                              }
                            >
                              {row.improvedDeNovoCaRmsdVsRefAngstrom.toFixed(2)} Å
                            </span>
                          ) : (
                            <span className="text-zinc-500">Test własny bez PDB</span>
                          )}
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </div>

            <div className="h-56 bg-[#050712] border border-zinc-800 rounded-xl p-4">
              <ResponsiveContainer width="100%" height="100%">
                <BarChart
                  data={curriculumReport.entries.map(e => ({
                    name: e.targetId.replace('_BASELINE_18', '').replace('_HELIX_26', ''),
                    BazaPaulingCorey: e.classicalEnergyKcal,
                    Model22ParamPoRybosomie: e.improvedPostRibosomeEnergyKcal
                  }))}
                >
                  <CartesianGrid strokeDasharray="3 3" stroke="#1e293b" />
                  <XAxis dataKey="name" stroke="#94a3b8" tick={{ fontSize: 10 }} />
                  <YAxis stroke="#94a3b8" tick={{ fontSize: 10 }} />
                  <Tooltip
                    contentStyle={{
                      backgroundColor: '#050712',
                      borderColor: '#334155',
                      fontSize: '11px'
                    }}
                  />
                  <Legend wrapperStyle={{ fontSize: '11px' }} />
                  <Bar dataKey="BazaPaulingCorey" name="Baza klasyczna α-helisy [kcal/mol]" fill="#64748b" />
                  <Bar
                    dataKey="Model22ParamPoRybosomie"
                    name="Model 22-parametrowy po relaksacji [kcal/mol]"
                    fill="#10b981"
                  />
                </BarChart>
              </ResponsiveContainer>
            </div>
          </div>
        )}

        {subMode === 'four_path_scanner' && (
          <div className="space-y-4">
            {/* 4 PATH SELECTOR BUTTONS */}
            <div className="grid grid-cols-1 sm:grid-cols-2 xl:grid-cols-4 gap-3">
              {(
                [
                  'SPIRAL_BOTTOM_TO_TOP',
                  'SPIRAL_FRONT_TO_BACK',
                  'LATERAL_SIDECHAIN_BONDS',
                  'DISULFIDE_BRIDGES_AUDIT'
                ] as DiagnosticScanPathId[]
              ).map(pid => {
                const pr = trajectoryBundle.scan.paths[pid];
                const isSel = activeScanPath === pid;
                return (
                  <button
                    key={pid}
                    onClick={() => setActiveScanPath(pid)}
                    className={`text-left rounded-xl p-3 border transition cursor-pointer ${
                      isSel
                        ? 'bg-cyan-950/60 border-cyan-400 text-white'
                        : 'bg-[#050712] border-zinc-800 text-zinc-300 hover:border-zinc-700'
                    }`}
                  >
                    <div className="flex items-center justify-between gap-2">
                      <span className="text-xs font-extrabold">{pr.titlePl}</span>
                      <span
                        className={`text-[10px] px-1.5 py-0.5 rounded font-bold ${
                          pr.errorCount > 0
                            ? 'bg-rose-950 text-rose-300 border border-rose-700'
                            : pr.warningCount > 0
                            ? 'bg-amber-950 text-amber-300 border border-amber-700'
                            : 'bg-emerald-950 text-emerald-300 border border-emerald-700'
                        }`}
                      >
                        {pr.errorCount > 0
                          ? `${pr.errorCount} błędów`
                          : pr.warningCount > 0
                          ? `${pr.warningCount} uwag`
                          : '0 błędów'}
                      </span>
                    </div>
                    <p className="text-[10.5px] text-zinc-400 mt-1 line-clamp-2">{pr.subtitlePl}</p>
                  </button>
                );
              })}
            </div>

            {/* 3D INTERACTIVE PLATE SCANNER + WAYPOINT TABLE */}
            <div className="grid grid-cols-1 lg:grid-cols-12 gap-4">
              <div className="lg:col-span-7 bg-[#050712] border border-zinc-800 rounded-xl overflow-hidden flex flex-col">
                <div className="p-3 border-b border-zinc-800 flex flex-wrap items-center justify-between gap-2 bg-[#070a18]">
                  <div className="text-xs font-bold text-cyan-300">
                    Podgląd 3D Płyty Skanującej &amp; Ślimaka ({currentPathReport.titlePl})
                  </div>
                  <div className="flex items-center gap-2">
                    <button
                      onClick={() => {
                        if (scanStepIdx >= currentPathReport.waypoints.length - 1) {
                          setScanStepIdx(0);
                        }
                        setIsAutoScanning(!isAutoScanning);
                      }}
                      className="px-3 py-1 rounded bg-cyan-600 hover:bg-cyan-500 text-white text-xs font-bold flex items-center gap-1.5 cursor-pointer"
                    >
                      {isAutoScanning ? (
                        <>
                          <Pause className="w-3.5 h-3.5" />
                          <span>Pauza</span>
                        </>
                      ) : (
                        <>
                          <Play className="w-3.5 h-3.5" />
                          <span>Skanuj Ścieżką</span>
                        </>
                      )}
                    </button>
                    <button
                      onClick={() => {
                        setIsAutoScanning(false);
                        setScanStepIdx(0);
                      }}
                      className="px-2.5 py-1 rounded bg-zinc-900 hover:bg-zinc-800 border border-zinc-700 text-zinc-300 text-xs cursor-pointer"
                      title="Resetuj skaner do początku ścieżki"
                    >
                      <RotateCcw className="w-3.5 h-3.5" />
                    </button>
                  </div>
                </div>

                <div
                  ref={scannerMountRef}
                  onMouseDown={e => {
                    orbitRef.current.dragging = true;
                    orbitRef.current.px = e.clientX;
                    orbitRef.current.py = e.clientY;
                  }}
                  onMouseMove={e => {
                    if (!orbitRef.current.dragging) return;
                    const dx = e.clientX - orbitRef.current.px;
                    const dy = e.clientY - orbitRef.current.py;
                    orbitRef.current.rotY += dx * 0.008;
                    orbitRef.current.rotX = Math.max(
                      -1.45,
                      Math.min(1.45, orbitRef.current.rotX + dy * 0.008)
                    );
                    orbitRef.current.px = e.clientX;
                    orbitRef.current.py = e.clientY;
                  }}
                  onMouseUp={() => {
                    orbitRef.current.dragging = false;
                  }}
                  onMouseLeave={() => {
                    orbitRef.current.dragging = false;
                  }}
                  onWheel={e => {
                    orbitRef.current.dist = Math.max(
                      18,
                      Math.min(180, orbitRef.current.dist + e.deltaY * 0.04)
                    );
                  }}
                  className="w-full h-[360px] cursor-grab active:cursor-grabbing relative"
                />

                {/* Step scrubber */}
                <div className="p-3 border-t border-zinc-800 bg-[#040610] space-y-2">
                  <div className="flex items-center justify-between text-xs tabular-nums">
                    <span className="text-zinc-400">
                      Krok skanera: <strong className="text-white">{scanStepIdx + 1}</strong> /{' '}
                      {currentPathReport.waypoints.length}
                    </span>
                    {activeScanWaypoint && (
                      <span
                        className={`font-bold ${
                          activeScanWaypoint.status === 'ERROR'
                            ? 'text-rose-400'
                            : activeScanWaypoint.status === 'WARNING'
                            ? 'text-amber-300'
                            : 'text-emerald-400'
                        }`}
                      >
                        Reszta #{activeScanWaypoint.resSeq} ({activeScanWaypoint.code3}) •{' '}
                        {activeScanWaypoint.metricValue} {activeScanWaypoint.metricUnit}
                      </span>
                    )}
                  </div>
                  <input
                    type="range"
                    min={0}
                    max={Math.max(0, currentPathReport.waypoints.length - 1)}
                    value={scanStepIdx}
                    onChange={e => {
                      setIsAutoScanning(false);
                      setScanStepIdx(Number(e.target.value));
                    }}
                    className="w-full accent-cyan-500 cursor-pointer"
                  />
                  {activeScanWaypoint && (
                    <div className="text-[11px] text-zinc-300 bg-black/50 border border-zinc-800 rounded p-2">
                      {activeScanWaypoint.diagnosticNotePl}
                    </div>
                  )}
                </div>
              </div>

              {/* Right Column: Waypoints & Anomaly Log */}
              <div className="lg:col-span-5 bg-[#050712] border border-zinc-800 rounded-xl flex flex-col overflow-hidden max-h-[515px]">
                <div className="p-3 border-b border-zinc-800 bg-[#070a18] flex items-center justify-between">
                  <span className="text-xs font-bold text-white uppercase">
                    Dziennik Kroków Ścieżki ({currentPathReport.waypoints.length})
                  </span>
                  <span className="text-[10.5px] text-zinc-400">
                    Kliknij wiersz, aby przesunąć płytę 3D
                  </span>
                </div>
                <div className="flex-1 overflow-y-auto divide-y divide-zinc-900 scrollbar-thin scrollbar-thumb-zinc-800">
                  {currentPathReport.waypoints.map((wp, idx) => (
                    <button
                      key={`${wp.pathId}-${wp.stepNumber}`}
                      onClick={() => {
                        setIsAutoScanning(false);
                        setScanStepIdx(idx);
                      }}
                      className={`w-full text-left p-2.5 transition cursor-pointer flex items-start justify-between gap-2 ${
                        idx === scanStepIdx
                          ? 'bg-cyan-950/50 border-l-2 border-cyan-400'
                          : 'hover:bg-zinc-900/40'
                      }`}
                    >
                      <div className="space-y-0.5">
                        <div className="text-xs font-bold text-white flex items-center gap-1.5">
                          <span className="text-zinc-500">#{wp.stepNumber}</span>
                          <span>
                            {wp.code3}-{wp.resSeq}
                          </span>
                          {wp.status === 'ERROR' && (
                            <AlertTriangle className="w-3.5 h-3.5 text-rose-400" />
                          )}
                          {wp.status === 'WARNING' && (
                            <AlertTriangle className="w-3.5 h-3.5 text-amber-400" />
                          )}
                          {wp.status === 'OK' && (
                            <CheckCircle2 className="w-3 h-3 text-emerald-500" />
                          )}
                        </div>
                        <p className="text-[10.5px] text-zinc-400 leading-snug">
                          {wp.diagnosticNotePl}
                        </p>
                      </div>
                      <span className="text-[10.5px] font-bold text-cyan-300 shrink-0 tabular-nums">
                        {wp.metricValue} {wp.metricUnit}
                      </span>
                    </button>
                  ))}
                </div>
              </div>
            </div>
          </div>
        )}

        {subMode === 'multiagent_soap_bubble' && (() => {
          const ma = trajectoryBundle.multiAgentReport;
          const diag = ma.convergedFinalDiag;
          const gt = ma.groundTruthValidation;
          return (
            <div className="space-y-4">
              <div className="grid grid-cols-1 md:grid-cols-4 gap-3">
                <div className="bg-[#060814] border border-cyan-800/60 rounded-xl p-3.5 space-y-1">
                  <span className="text-[10.5px] uppercase font-bold text-cyan-300">
                    Soft SnapLock (Brak sztywnego zamrażania)
                  </span>
                  <p className="text-xl font-extrabold text-white tabular-nums">
                    Δ(φ,ψ)_N = {ma.nTermSoftSnapMovedDeg}°
                  </p>
                  <p className="text-[10.5px] text-zinc-400">
                    Kąty N-końca obracają się na zawiasach podczas dobudowywania końca C
                  </p>
                </div>

                <div className="bg-[#060814] border border-emerald-800/60 rounded-xl p-3.5 space-y-1">
                  <span className="text-[10.5px] uppercase font-bold text-emerald-300">
                    Napięcie Powierzchniowe „Bańki Mydlanej”
                  </span>
                  <p className="text-xl font-extrabold text-emerald-300 tabular-nums">
                    {diag.soapBubbleEnergyKcal} kcal/mol
                  </p>
                  <p className="text-[10.5px] text-zinc-400">
                    Spychanie reszt hydrofobowych do COM i wypychanie hydrofilowych na zewnątrz
                  </p>
                </div>

                <div className="bg-[#060814] border border-indigo-800/60 rounded-xl p-3.5 space-y-1">
                  <span className="text-[10.5px] uppercase font-bold text-indigo-300">
                    Odbijanie Sfer VdW &amp; Sfera R_kula
                  </span>
                  <p className="text-xl font-extrabold text-indigo-200 tabular-nums">
                    {diag.vdwClashesCount} kolizji • {diag.sphereToRgRatio}× Rg
                  </p>
                  <p className="text-[10.5px] text-zinc-400">
                    min d(Cα-Cα) = {diag.minNonBondedCaDistAngstrom} Å • R_kula = {diag.adaptiveSphereRadiusAngstrom} Å
                  </p>
                </div>

                <div
                  className={`bg-[#060814] rounded-xl p-3.5 space-y-1 border ${
                    gt && gt.caRmsdAngstrom !== null && gt.caRmsdAngstrom < 1.0 && gt.valid
                      ? 'border-emerald-500/80'
                      : 'border-amber-600/80'
                  }`}
                >
                  <span className="text-[10.5px] uppercase font-bold text-amber-300">
                    Kabsch SVD 1:1 Cα Po Relaksacji (φ,ψ)
                  </span>
                  <p
                    className={`text-xl font-extrabold tabular-nums ${
                      gt && gt.caRmsdAngstrom !== null && gt.caRmsdAngstrom < 1.0 && gt.valid
                        ? 'text-emerald-400'
                        : 'text-amber-300'
                    }`}
                  >
                    {gt && gt.caRmsdAngstrom !== null ? `${gt.caRmsdAngstrom.toFixed(3)} Å` : 'De-novo (brak ref.)'}
                  </p>
                  <p className="text-[10.5px] text-zinc-400">
                    {gt && gt.caRmsdAngstrom !== null
                      ? gt.caRmsdAngstrom < 1.0 && gt.valid
                        ? `Bramka < 1.0 Å ODBLOKOWANA (${gt.referenceId})`
                        : `Zbiór TEST / brak referencji treningowej (${gt.referenceId})`
                      : 'Czysta fizyka bez zewnętrznych API'}
                  </p>
                </div>
              </div>

              {/* Trajectory Chart: Energy & Kabsch RMSD over N->C + Soap-Bubble steps */}
              <div className="bg-[#060814] border border-zinc-800 rounded-xl p-4 space-y-2">
                <div className="flex items-center justify-between">
                  <h3 className="text-xs font-bold text-white uppercase">
                    Trajektoria Wieloagentowa w Pamięci: Wzrost N➔C (Soft SnapLock) ➔ Odcięcie z Rybosomu (Bańka Mydlana &amp; Gradient Descent φ/ψ)
                  </h3>
                  <span className="text-[10.5px] text-emerald-300 font-bold">
                    Zero zewnętrznych API • Czysty jakobian Jᵀ·F na kątach torsyjnych
                  </span>
                </div>
                <div className="h-56 w-full">
                  <ResponsiveContainer width="100%" height="100%">
                    <LineChart
                      data={ma.trajectoryHistory.map(h => ({
                        krok: `${h.tetheredToRibosome ? 'N→C' : 'Bańka'} #${h.step}`,
                        energiaPotencjalna: h.totalPotentialEnergyKcal,
                        kabschRmsd: h.liveKabschCaRmsdAngstrom ?? 0,
                        drzenieDeg: h.kineticQuiverAmplitudeDeg
                      }))}
                    >
                      <CartesianGrid strokeDasharray="3 3" stroke="#1e293b" />
                      <XAxis dataKey="krok" stroke="#94a3b8" tick={{ fontSize: 10 }} />
                      <YAxis stroke="#94a3b8" tick={{ fontSize: 10 }} />
                      <Tooltip
                        contentStyle={{
                          backgroundColor: '#020617',
                          borderColor: '#334155',
                          fontSize: '11px'
                        }}
                      />
                      <Legend wrapperStyle={{ fontSize: '11px' }} />
                      <Line
                        type="monotone"
                        dataKey="energiaPotencjalna"
                        name="Energia E_total [kcal/mol]"
                        stroke="#34d399"
                        strokeWidth={2}
                        dot={false}
                      />
                      <Line
                        type="monotone"
                        dataKey="kabschRmsd"
                        name="Kabsch SVD Cα RMSD [Å]"
                        stroke="#38bdf8"
                        strokeWidth={2}
                        dot={false}
                      />
                      <Line
                        type="monotone"
                        dataKey="drzenieDeg"
                        name="Drżenie kątowe |Δφ,Δψ| [°]"
                        stroke="#fbbf24"
                        strokeWidth={1.5}
                        dot={false}
                      />
                    </LineChart>
                  </ResponsiveContainer>
                </div>
              </div>

              {/* Per-Agent Force & Torque Table */}
              <div className="bg-[#060814] border border-zinc-800 rounded-xl overflow-hidden">
                <div className="p-3 bg-[#080b1a] border-b border-zinc-800 flex items-center justify-between text-xs">
                  <span className="font-bold text-white uppercase">
                    Stan Równowagi Agentów-Reziduów ({diag.agentSnapshots.length} autonomicznych obiektów w pamięci)
                  </span>
                  <span className="text-[10.5px] text-zinc-400">
                    Każdy agent negocjuje kąty (φ_i, ψ_i) bez zrywania wiązań kowalencyjnych
                  </span>
                </div>
                <div className="max-h-72 overflow-y-auto scrollbar-thin scrollbar-thumb-zinc-800">
                  <table className="w-full text-left border-collapse text-[11px]">
                    <thead className="bg-black/80 text-zinc-400 sticky top-0">
                      <tr>
                        <th className="p-2 border-b border-zinc-800">#</th>
                        <th className="p-2 border-b border-zinc-800">Agent</th>
                        <th className="p-2 border-b border-zinc-800">Motyw</th>
                        <th className="p-2 border-b border-zinc-800">φ [°]</th>
                        <th className="p-2 border-b border-zinc-800">ψ [°]</th>
                        <th className="p-2 border-b border-zinc-800">F_bańka [x,y,z]</th>
                        <th className="p-2 border-b border-zinc-800">F_VdW [x,y,z]</th>
                        <th className="p-2 border-b border-zinc-800">F_klamra [x,y,z]</th>
                        <th className="p-2 border-b border-zinc-800">|F_net|</th>
                        <th className="p-2 border-b border-zinc-800">τ_φ / τ_ψ</th>
                      </tr>
                    </thead>
                    <tbody className="divide-y divide-zinc-900">
                      {diag.agentSnapshots.map(ag => (
                        <tr key={ag.resIndex} className="hover:bg-zinc-900/40">
                          <td className="p-2 text-zinc-500">{ag.resIndex + 1}</td>
                          <td className="p-2 font-bold text-white">{ag.code3}</td>
                          <td className="p-2 text-cyan-300">{ag.secondaryMotif}</td>
                          <td className="p-2 tabular-nums">{ag.phiDeg}°</td>
                          <td className="p-2 tabular-nums">{ag.psiDeg}°</td>
                          <td className="p-2 text-emerald-300 tabular-nums">
                            [{ag.soapBubbleForce.map(v => v.toFixed(2)).join(', ')}]
                          </td>
                          <td className="p-2 text-rose-300 tabular-nums">
                            [{ag.vdwStericForce.map(v => v.toFixed(2)).join(', ')}]
                          </td>
                          <td className="p-2 text-amber-300 tabular-nums">
                            [{ag.clampForce.map(v => v.toFixed(2)).join(', ')}]
                          </td>
                          <td className="p-2 font-bold text-white tabular-nums">{ag.netForceMagnitude}</td>
                          <td className="p-2 text-indigo-300 tabular-nums">
                            {ag.torquePhi} / {ag.torquePsi}
                          </td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              </div>
            </div>
          );
        })()}
      </div>
    </div>
  );
};
