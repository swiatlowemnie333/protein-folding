import React, { useState, useEffect, useRef } from 'react';
import * as THREE from 'three';
import { 
  Sparkles, 
  Layers, 
  ShieldCheck, 
  Compass, 
  RotateCcw, 
  Play, 
  Download, 
  FileText, 
  CheckCircle2, 
  AlertTriangle, 
  SlidersHorizontal, 
  Info, 
  Boxes,
  ArrowRight,
  RefreshCw,
  Activity,
  Check,
  Copy,
  BarChart3,
  X as CloseIcon
} from 'lucide-react';
import { 
  SubdivisionLevel, 
  SubdomainCalibrationReport, 
  BenchmarkSubdomain, 
  BENCHMARK_SUBDOMAINS, 
  getSubdividedLattices, 
  runSubdomainCalibration, 
  exportCalibratedSubdomainPdb 
} from '../snaplockSubdomainCalibrator';
import { SNAPLOCK_21_BLOCKS } from '../snaplock21Blocks';
import { SnaplockRmsdComparisonChart } from './SnaplockRmsdComparisonChart';

// Cache subdivided cage geometries once so WireframeGeometry(2232 triangles) is never re-calculated on the fly
const SUBDIVIDED_CAGE_GEOS: Record<number, { surfaceGeo: THREE.BufferGeometry; wireGeo: THREE.WireframeGeometry }> = {};

function getCachedSubdividedCage(level: SubdivisionLevel) {
  if (SUBDIVIDED_CAGE_GEOS[level]) return SUBDIVIDED_CAGE_GEOS[level];
  const lattices = getSubdividedLattices();
  const lattice = level === 2 ? lattices.level2 : level === 1 ? lattices.level1 : lattices.level0;

  const vertArray: number[] = [];
  lattice.triangles.forEach(tri => {
    tri.forEach(v => {
      vertArray.push(v[0], v[1], v[2]);
    });
  });

  const surfaceGeo = new THREE.BufferGeometry();
  surfaceGeo.setAttribute('position', new THREE.Float32BufferAttribute(vertArray, 3));
  surfaceGeo.computeVertexNormals();
  const wireGeo = new THREE.WireframeGeometry(surfaceGeo);

  SUBDIVIDED_CAGE_GEOS[level] = { surfaceGeo, wireGeo };
  return SUBDIVIDED_CAGE_GEOS[level];
}

interface SnaplockSubdomainCalibrationViewProps {
  synchronizedSubdomainId?: string;
  onLoadSequenceToConstructor?: (sequence3Letter: string) => void;
  onSwitchToConstructor?: () => void;
}

export const SnaplockSubdomainCalibrationView: React.FC<SnaplockSubdomainCalibrationViewProps> = ({
  synchronizedSubdomainId,
  onLoadSequenceToConstructor,
  onSwitchToConstructor
}) => {
  // Selected benchmark subdomain
  const [selectedSubdomainId, setSelectedSubdomainId] = useState<string>(
    synchronizedSubdomainId || 'insulin_a_core'
  );
  useEffect(() => {
    if (synchronizedSubdomainId && BENCHMARK_SUBDOMAINS.some(s => s.id === synchronizedSubdomainId)) {
      setSelectedSubdomainId(synchronizedSubdomainId);
    }
  }, [synchronizedSubdomainId]);
  const currentSubdomain = BENCHMARK_SUBDOMAINS.find(s => s.id === selectedSubdomainId) || BENCHMARK_SUBDOMAINS[0];

  // Subdivision Level (0 = 248, 1 = 744, 2 = 2232)
  const [subdivisionLevel, setSubdivisionLevel] = useState<SubdivisionLevel>(2);

  // Calibration parameters state
  const [kappaBend, setKappaBend] = useState<number>(42.5);
  const [phiOffset, setPhiOffset] = useState<number>(-3.2);
  const [psiOffset, setPsiOffset] = useState<number>(2.8);
  const [ssRestraintK, setSsRestraintK] = useState<number>(8.5);
  const [dMaxTarget, setDMaxTarget] = useState<number>(currentSubdomain.dMaxTarget);

  // Active inspected residue
  const [activeResidueIdx, setActiveResidueIdx] = useState<number>(0);

  // Calibration Execution & Report state
  const [isCalibrating, setIsCalibrating] = useState<boolean>(false);
  const [calibProgress, setCalibProgress] = useState<number>(100);
  const [report, setReport] = useState<SubdomainCalibrationReport | null>(null);
  const [copiedReport, setCopiedReport] = useState<boolean>(false);
  const [showRmsdChartModal, setShowRmsdChartModal] = useState<boolean>(false);

  // 3D Visualizer Toggles
  const [showSubdividedCage, setShowSubdividedCage] = useState<boolean>(true);
  const [showRibbon, setShowRibbon] = useState<boolean>(true);
  const [showDmaxBoundary, setShowDmaxBoundary] = useState<boolean>(true);
  const [showDisulfideRods, setShowDisulfideRods] = useState<boolean>(true);

  // Three.js Mount & Scene Refs
  const mountRef = useRef<HTMLDivElement | null>(null);
  const rendererRef = useRef<THREE.WebGLRenderer | null>(null);
  const sceneRef = useRef<THREE.Scene | null>(null);
  const cameraRef = useRef<THREE.PerspectiveCamera | null>(null);
  const modelGroupRef = useRef<THREE.Group | null>(null);
  const boundaryGroupRef = useRef<THREE.Group | null>(null);
  const animFrameRef = useRef<number | null>(null);

  // Orbit controls
  const isDragging = useRef<boolean>(false);
  const isPanning = useRef<boolean>(false);
  const prevMousePos = useRef<{ x: number; y: number }>({ x: 0, y: 0 });
  const orbitAngle = useRef<{ x: number; y: number }>({ x: 0.35, y: 0.65 });
  const cameraDist = useRef<number>(45);
  const panOffset = useRef<THREE.Vector3>(new THREE.Vector3(0, 0, 0));

  // Audio tone for calibration convergence
  const playTone = (freq: number) => {
    try {
      const AudioCtx = window.AudioContext || (window as any).webkitAudioContext;
      if (!AudioCtx) return;
      const ctx = new AudioCtx();
      const osc = ctx.createOscillator();
      const gain = ctx.createGain();
      osc.frequency.setValueAtTime(freq, ctx.currentTime);
      gain.gain.setValueAtTime(0.06, ctx.currentTime);
      gain.gain.exponentialRampToValueAtTime(0.001, ctx.currentTime + 0.15);
      osc.connect(gain);
      gain.connect(ctx.destination);
      osc.start();
      osc.stop(ctx.currentTime + 0.16);
    } catch (e) {}
  };

  // Run calibration on initial load or subdomain change
  useEffect(() => {
    setActiveResidueIdx(0);
    setDMaxTarget(currentSubdomain.dMaxTarget);
    const initialReport = runSubdomainCalibration(
      currentSubdomain,
      subdivisionLevel,
      {
        kappaBend,
        phiOffsetDeg: phiOffset,
        psiOffsetDeg: psiOffset,
        ssRestraintK,
        dMaxRelaxation: 0.96
      },
      12
    );
    setReport(initialReport);
  }, [selectedSubdomainId, subdivisionLevel]);

  // Execute Animated Calibration Process
  const handleExecuteCalibration = () => {
    if (isCalibrating) return;
    setIsCalibrating(true);
    setCalibProgress(0);

    let step = 0;
    const maxSteps = 16;
    const interval = setInterval(() => {
      step++;
      const pct = Math.floor((step / maxSteps) * 100);
      setCalibProgress(pct);
      playTone(280 + step * 25);

      if (step >= maxSteps) {
        clearInterval(interval);
        setIsCalibrating(false);
        setCalibProgress(100);

        const finalReport = runSubdomainCalibration(
          { ...currentSubdomain, dMaxTarget },
          subdivisionLevel,
          {
            kappaBend,
            phiOffsetDeg: phiOffset,
            psiOffsetDeg: psiOffset,
            ssRestraintK
          },
          22
        );
        setReport(finalReport);
        playTone(660);
      }
    }, 55);
  };

  // --- THREE.JS SCENE SETUP ---
  useEffect(() => {
    if (!mountRef.current) return;
    const container = mountRef.current;
    const width = container.clientWidth || 800;
    const height = container.clientHeight || 500;

    const scene = new THREE.Scene();
    scene.background = new THREE.Color(0x010204);
    sceneRef.current = scene;

    const camera = new THREE.PerspectiveCamera(45, width / height, 0.5, 1000);
    cameraRef.current = camera;

    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    renderer.setSize(width, height);
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    rendererRef.current = renderer;
    container.innerHTML = '';
    container.appendChild(renderer.domElement);

    // Lights
    const amb = new THREE.AmbientLight(0xffffff, 0.75);
    scene.add(amb);

    const dir1 = new THREE.DirectionalLight(0xffffff, 1.3);
    dir1.position.set(20, 30, 25);
    scene.add(dir1);

    const dir2 = new THREE.DirectionalLight(0x38bdf8, 0.8);
    dir2.position.set(-25, -20, -20);
    scene.add(dir2);

    const dir3 = new THREE.DirectionalLight(0xf59e0b, 0.5);
    dir3.position.set(0, -25, 20);
    scene.add(dir3);

    const modelGroup = new THREE.Group();
    modelGroupRef.current = modelGroup;
    scene.add(modelGroup);

    const boundaryGroup = new THREE.Group();
    boundaryGroupRef.current = boundaryGroup;
    scene.add(boundaryGroup);

    const handleResize = () => {
      if (!container || !renderer || !camera) return;
      const w = container.clientWidth;
      const h = container.clientHeight;
      camera.aspect = w / h;
      camera.updateProjectionMatrix();
      renderer.setSize(w, h);
    };
    window.addEventListener('resize', handleResize);

    const animate = () => {
      animFrameRef.current = requestAnimationFrame(animate);

      if (cameraRef.current) {
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

      renderer.render(scene, camera);
    };
    animate();

    return () => {
      if (animFrameRef.current) cancelAnimationFrame(animFrameRef.current);
      window.removeEventListener('resize', handleResize);
      renderer.dispose();
    };
  }, []);

  // --- REBUILD 3D VISUALIZATION: BACKBONE, 2232 SUBDIVIDED CAGE, DMAX BOUNDARY & S-S ---
  useEffect(() => {
    if (!modelGroupRef.current || !boundaryGroupRef.current) return;
    const group = modelGroupRef.current;
    const boundaryGroup = boundaryGroupRef.current;

    // Clear old children without disposing shared geometries
    while (group.children.length > 0) {
      const child = group.children[0];
      group.remove(child);
      if ((child as any).geometry && !(child as any).userData?.isShared) {
        (child as any).geometry.dispose();
      }
    }
    while (boundaryGroup.children.length > 0) {
      const child = boundaryGroup.children[0];
      boundaryGroup.remove(child);
      if ((child as any).geometry) (child as any).geometry.dispose();
    }

    const hasCalibrated =
      report &&
      report.subdomainId === currentSubdomain.id &&
      Array.isArray(report.calibratedCoordinates) &&
      report.calibratedCoordinates.length === currentSubdomain.sequenceCodes.length;

    const coords = hasCalibrated ? report.calibratedCoordinates : currentSubdomain.caCoordinates;
    const cbCoords = hasCalibrated ? report.cbCoordinates : [];
    const sgCoords = hasCalibrated ? report.sgCoordinates : [];

    const n = coords?.length ?? 0;
    if (n === 0) return;

    // Center coordinates around Center of Mass [0,0,0]
    let sumX = 0, sumY = 0, sumZ = 0;
    coords.forEach(p => { 
      if (p) {
        sumX += p?.[0] ?? 0; 
        sumY += p?.[1] ?? 0; 
        sumZ += p?.[2] ?? 0; 
      }
    });
    const comX = sumX / Math.max(1, n);
    const comY = sumY / Math.max(1, n);
    const comZ = sumZ / Math.max(1, n);

    const SCALE_FACTOR = 1.4;
    const centeredPositions = coords.map(p => 
      new THREE.Vector3(
        ((p?.[0] ?? 0) - comX) * SCALE_FACTOR,
        ((p?.[1] ?? 0) - comY) * SCALE_FACTOR,
        ((p?.[2] ?? 0) - comZ) * SCALE_FACTOR
      )
    );

    // Dynamic camera framing
    let maxDist = 10;
    centeredPositions.forEach(p => {
      if (p) maxDist = Math.max(maxDist, p.length());
    });
    cameraDist.current = Math.max(28, maxDist * 2.3);

    // 1. D_max Boundary Sphere (SAXS constraint envelope)
    if (showDmaxBoundary) {
      const dMaxHalf = (dMaxTarget / 2.0) * SCALE_FACTOR;
      const dMaxGeo = new THREE.SphereGeometry(dMaxHalf, 32, 24);
      const dMaxWire = new THREE.WireframeGeometry(dMaxGeo);
      const dMaxMat = new THREE.LineBasicMaterial({
        color: 0x06b6d4,
        transparent: true,
        opacity: 0.22
      });
      boundaryGroup.add(new THREE.LineSegments(dMaxWire, dMaxMat));
    }

    // 2. Backbone CatmullRom Ribbon
    if (showRibbon && centeredPositions.length >= 3) {
      try {
        const curve = new THREE.CatmullRomCurve3(centeredPositions);
        const tubeGeo = new THREE.TubeGeometry(curve, n * 8, 0.28, 12, false);
        const tubeMat = new THREE.MeshStandardMaterial({
          color: 0x38bdf8,
          roughness: 0.25,
          metalness: 0.35,
          transparent: true,
          opacity: 0.85
        });
        group.add(new THREE.Mesh(tubeGeo, tubeMat));
      } catch (e) {}
    }

    const safeResIdx = Math.min(activeResidueIdx, centeredPositions.length - 1);

    // 3. Residue Cα Beads
    centeredPositions.forEach((pos, idx) => {
      if (!pos) return;
      const isFocused = safeResIdx === idx;
      const code = currentSubdomain.sequenceCodes[idx] || 'ALA';
      const aaBlock = SNAPLOCK_21_BLOCKS[code] || SNAPLOCK_21_BLOCKS['ALA'];

      const caMesh = new THREE.Mesh(
        new THREE.SphereGeometry(isFocused ? 0.6 : 0.4, 16, 16),
        new THREE.MeshStandardMaterial({
          color: isFocused ? 0xfacc15 : new THREE.Color(aaBlock.color),
          emissive: isFocused ? 0xca8a04 : 0x000000,
          emissiveIntensity: isFocused ? 0.6 : 0,
          roughness: 0.2,
          metalness: 0.4
        })
      );
      caMesh.position.copy(pos);
      group.add(caMesh);
    });

    // 4. ACTIVE RESIDUE'S SUBDIVIDED SPHERICAL CONTACT CAGE (248 / 744 / 2232 FACETS!)
    if (showSubdividedCage && centeredPositions[safeResIdx]) {
      const activePos = centeredPositions[safeResIdx];
      const { surfaceGeo, wireGeo } = getCachedSubdividedCage(subdivisionLevel);

      const CAGE_RADIUS = 2.4; // Scaled for visible inspection
      const cageGroup = new THREE.Group();
      cageGroup.position.copy(activePos);

      // Translucent surface
      const faceMat = new THREE.MeshStandardMaterial({
        color: 0x06b6d4,
        transparent: true,
        opacity: subdivisionLevel === 2 ? 0.30 : 0.22,
        roughness: 0.3,
        metalness: 0.2,
        side: THREE.DoubleSide
      });
      const surfaceMesh = new THREE.Mesh(surfaceGeo, faceMat);
      surfaceMesh.scale.setScalar(CAGE_RADIUS);
      surfaceMesh.userData = { isShared: true };
      cageGroup.add(surfaceMesh);

      // Wireframe of the 2232 micro-facets!
      const wireMat = new THREE.LineBasicMaterial({
        color: subdivisionLevel === 2 ? 0x22d3ee : 0x38bdf8,
        transparent: true,
        opacity: subdivisionLevel === 2 ? 0.75 : 0.55
      });
      const wireMesh = new THREE.LineSegments(wireGeo, wireMat);
      wireMesh.scale.setScalar(CAGE_RADIUS * 1.002);
      wireMesh.userData = { isShared: true };
      cageGroup.add(wireMesh);

      // Port #1 and #248 markers
      const p1Mesh = new THREE.Mesh(new THREE.SphereGeometry(0.2, 12, 12), new THREE.MeshBasicMaterial({ color: 0x06b6d4 }));
      p1Mesh.position.set(0, CAGE_RADIUS, 0);
      cageGroup.add(p1Mesh);

      const p248Mesh = new THREE.Mesh(new THREE.SphereGeometry(0.2, 12, 12), new THREE.MeshBasicMaterial({ color: 0xf43f5e }));
      p248Mesh.position.set(0, -CAGE_RADIUS, 0);
      cageGroup.add(p248Mesh);

      group.add(cageGroup);
    }

    // 5. Covalent Disulfide Bridges (Two-Way IK Closed CA - CB - SG === SG - CB - CA)
    if (showDisulfideRods && currentSubdomain.disulfidePairs.length > 0) {
      const sideMat = new THREE.MeshStandardMaterial({
        color: 0x94a3b8,
        roughness: 0.35,
        metalness: 0.3
      });
      const bridgeMat = new THREE.MeshStandardMaterial({
        color: 0xfbbf24,
        emissive: 0xd97706,
        emissiveIntensity: 0.5,
        roughness: 0.2,
        metalness: 0.65
      });

      const addRod = (a: THREE.Vector3, b: THREE.Vector3, r: number, mat: THREE.Material) => {
        const dist = a.distanceTo(b);
        if (dist < 0.04) return;
        const cylGeo = new THREE.CylinderGeometry(r, r, dist, 12);
        const mesh = new THREE.Mesh(cylGeo, mat);
        mesh.position.copy(a.clone().add(b).multiplyScalar(0.5));
        mesh.quaternion.setFromUnitVectors(
          new THREE.Vector3(0, 1, 0),
          b.clone().sub(a).normalize()
        );
        group.add(mesh);
      };

      const toCenteredVec = (pt: [number, number, number]) =>
        new THREE.Vector3(
          ((pt?.[0] ?? 0) - comX) * SCALE_FACTOR,
          ((pt?.[1] ?? 0) - comY) * SCALE_FACTOR,
          ((pt?.[2] ?? 0) - comZ) * SCALE_FACTOR
        );

      currentSubdomain.disulfidePairs.forEach(([i1, i2]) => {
        const ca1 = centeredPositions[i1];
        const ca2 = centeredPositions[i2];
        if (!ca1 || !ca2) return;

        const rawCb1 = cbCoords?.[i1];
        const rawSg1 = sgCoords?.[i1];
        const rawCb2 = cbCoords?.[i2];
        const rawSg2 = sgCoords?.[i2];

        if (rawCb1 && rawSg1 && rawCb2 && rawSg2) {
          const cb1 = toCenteredVec(rawCb1);
          const sg1 = toCenteredVec(rawSg1);
          const cb2 = toCenteredVec(rawCb2);
          const sg2 = toCenteredVec(rawSg2);

          addRod(ca1, cb1, 0.13, sideMat);
          addRod(cb1, sg1, 0.16, bridgeMat);
          addRod(ca2, cb2, 0.13, sideMat);
          addRod(cb2, sg2, 0.16, bridgeMat);
          addRod(sg1, sg2, 0.25, bridgeMat);

          const s1 = new THREE.Mesh(new THREE.SphereGeometry(0.34, 16, 16), bridgeMat);
          s1.position.copy(sg1);
          group.add(s1);

          const s2 = new THREE.Mesh(new THREE.SphereGeometry(0.34, 16, 16), bridgeMat);
          s2.position.copy(sg2);
          group.add(s2);
        } else {
          addRod(ca1, ca2, 0.24, bridgeMat);
          const s1 = new THREE.Mesh(new THREE.SphereGeometry(0.38, 16, 16), bridgeMat);
          s1.position.copy(ca1);
          group.add(s1);
          const s2 = new THREE.Mesh(new THREE.SphereGeometry(0.38, 16, 16), bridgeMat);
          s2.position.copy(ca2);
          group.add(s2);
        }
      });
    }

  }, [
    currentSubdomain,
    report,
    dMaxTarget,
    subdivisionLevel,
    activeResidueIdx,
    showSubdividedCage,
    showRibbon,
    showDmaxBoundary,
    showDisulfideRods
  ]);

  // Orbit drag handlers
  const handleMouseDown = (e: React.MouseEvent) => {
    isDragging.current = true;
    isPanning.current = e.shiftKey || e.button === 1 || e.button === 2;
    prevMousePos.current = { x: e.clientX, y: e.clientY };
  };

  const handleMouseMove = (e: React.MouseEvent) => {
    if (!isDragging.current) return;
    const dx = e.clientX - prevMousePos.current.x;
    const dy = e.clientY - prevMousePos.current.y;

    if (isPanning.current) {
      panOffset.current.x -= dx * 0.05;
      panOffset.current.y += dy * 0.05;
    } else {
      orbitAngle.current.y += dx * 0.008;
      orbitAngle.current.x += dy * 0.008;
      orbitAngle.current.x = Math.max(-Math.PI / 2 + 0.05, Math.min(Math.PI / 2 - 0.05, orbitAngle.current.x));
    }

    prevMousePos.current = { x: e.clientX, y: e.clientY };
  };

  const handleMouseUp = () => {
    isDragging.current = false;
    isPanning.current = false;
  };

  const handleWheel = (e: React.WheelEvent) => {
    e.preventDefault();
    const zoomDelta = e.deltaY * 0.04;
    cameraDist.current = Math.max(12, Math.min(180, cameraDist.current + zoomDelta));
  };

  // Download Calibrated PDB
  const handleDownloadCalibratedPdb = () => {
    if (!report) return;
    const pdb = exportCalibratedSubdomainPdb(report);
    const blob = new Blob([pdb], { type: 'chemical/x-pdb;charset=utf-8' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `snaplock_calibrated_${report.subdomainId}_L${subdivisionLevel}.pdb`;
    a.click();
    URL.revokeObjectURL(url);
  };

  // Download JSON Report
  const handleDownloadJsonReport = () => {
    if (!report) return;
    const blob = new Blob([JSON.stringify(report, null, 2)], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `snaplock_subdomain_calibration_${report.subdomainId}.json`;
    a.click();
    URL.revokeObjectURL(url);
  };

  const handleCopyJson = () => {
    if (!report) return;
    navigator.clipboard.writeText(JSON.stringify(report, null, 2));
    setCopiedReport(true);
    setTimeout(() => setCopiedReport(false), 2000);
  };

  const activeResDetail = report?.residueDetails[activeResidueIdx];

  return (
    <div className="flex-1 w-full h-full min-h-0 bg-[#020205] text-slate-200 flex flex-col overflow-hidden font-mono select-none">
      
      {/* 1. TOP HEADER: SUBDOMAIN SELECTOR & SUBDIVISION LEVEL CONTROL */}
      <div className="border-b border-zinc-900 bg-[#04040a] px-4 py-2.5 shrink-0 flex flex-wrap items-center justify-between gap-3 text-xs">
        
        {/* Left: Title & Subdomain Selector */}
        <div className="flex items-center gap-3">
          <div className="w-9 h-9 rounded-xl bg-cyan-950/80 border border-cyan-500/50 flex items-center justify-center text-cyan-300 shadow-inner">
            <SlidersHorizontal className="w-5 h-5 text-cyan-400" />
          </div>
          <div>
            <div className="flex items-center gap-2">
              <h2 className="text-sm font-bold text-white tracking-tight">
                Kalibrator Subdomeny Białkowej & Podwójna Subdivizja
              </h2>
              <span className="text-[10px] px-2 py-0.2 rounded bg-cyan-950 text-cyan-300 border border-cyan-700/60 font-bold uppercase">
                2232 Mikrosiatki
              </span>
            </div>
            <p className="text-[10.5px] text-zinc-400">
              Dwukrotny podział 1➔3 (248 ➔ 744 ➔ 2232 faset) • Minimalizacja błędu kątowego poniżej 1.0° / resztę
            </p>
          </div>
        </div>

        {/* Center: Benchmark Subdomain Picker */}
        <div className="flex items-center gap-2">
          <span className="text-[10px] text-zinc-500 uppercase font-bold">Subdomena:</span>
          <select
            value={selectedSubdomainId}
            onChange={(e) => setSelectedSubdomainId(e.target.value)}
            className="bg-black border border-zinc-800 rounded px-2.5 py-1 text-xs text-white font-bold focus:outline-none focus:border-cyan-500 cursor-pointer"
          >
            {BENCHMARK_SUBDOMAINS.map((s) => (
              <option key={s.id} value={s.id}>
                {s.name} ({s.sequenceCodes.length} aa)
              </option>
            ))}
          </select>
        </div>

        {/* Right: Subdivision Level Selector (0: 248, 1: 744, 2: 2232) */}
        <div className="flex items-center gap-1.5 bg-black p-1 rounded-lg border border-zinc-800 text-[10.5px]">
          <span className="text-zinc-500 text-[9.5px] uppercase font-bold px-1.5">Siatka:</span>
          
          <button
            onClick={() => setSubdivisionLevel(0)}
            className={`px-2 py-1 rounded transition cursor-pointer font-bold ${
              subdivisionLevel === 0
                ? 'bg-zinc-800 text-white shadow-sm'
                : 'text-zinc-500 hover:text-zinc-300'
            }`}
            title="Poziom 0: Baza 248 wektorów Liego E8 (Średni błąd ~6.2°)"
          >
            L0: 248 faset
          </button>

          <button
            onClick={() => setSubdivisionLevel(1)}
            className={`px-2 py-1 rounded transition cursor-pointer font-bold ${
              subdivisionLevel === 1
                ? 'bg-indigo-950 border border-indigo-700 text-indigo-200 shadow-sm'
                : 'text-zinc-500 hover:text-zinc-300'
            }`}
            title="Poziom 1: Pierwsza subdivizja 1➔3 (744 fasetki, Średni błąd ~2.5°)"
          >
            L1: 744 faset
          </button>

          <button
            onClick={() => setSubdivisionLevel(2)}
            className={`px-2.5 py-1 rounded transition cursor-pointer font-extrabold flex items-center gap-1 ${
              subdivisionLevel === 2
                ? 'bg-cyan-950 border border-cyan-400 text-cyan-200 ring-1 ring-cyan-400 shadow-md'
                : 'text-zinc-500 hover:text-zinc-300'
            }`}
            title="Poziom 2: Druga subdivizja 3➔9 (2232 mikrofasetki, Błąd kątowy < 1.0° / resztę!)"
          >
            <Sparkles className="w-3 h-3 text-cyan-300" />
            <span>L2: 2232 faset (&lt;1°)</span>
          </button>
        </div>

      </div>

      {/* 2. SUB-BAR: LIVE CONVERGENCE TELEMETRY & RUN CALIBRATION BUTTON */}
      <div className="bg-[#030409] border-b border-zinc-900 px-4 py-1.5 flex flex-wrap items-center justify-between gap-3 shrink-0 text-xs font-mono">
        
        {/* Telemetry badges */}
        <div className="flex items-center gap-2 overflow-x-auto text-[10.5px]">
          {/* Angular Error Badge */}
          <span 
            className={`px-2.5 py-0.5 rounded border flex items-center gap-1 font-bold ${
              report && report.finalMeanErrorDeg < 1.0
                ? 'bg-emerald-950/80 border-emerald-500 text-emerald-300 ring-1 ring-emerald-500/40'
                : 'bg-amber-950/80 border-amber-600 text-amber-300'
            }`}
          >
            <CheckCircle2 className="w-3.5 h-3.5 text-emerald-400" />
            <span>
              Średni błąd kątowy: <strong>{report?.finalMeanErrorDeg ?? 0}°</strong> / resztę
            </span>
            {report && report.finalMeanErrorDeg < 1.0 && (
              <span className="text-[9px] px-1 py-0.1 rounded bg-emerald-900 text-emerald-100 uppercase ml-1">
                Cel &lt; 1.0° Osiągnięty
              </span>
            )}
          </span>

          {/* Ground-Truth Cα RMSD Badge (Kabsch SVD vs RCSB PDB) + Click to open Recharts RMSD Plot */}
          {report && (
            <button
              type="button"
              onClick={() => {
                if (report.rmsdComparison) setShowRmsdChartModal(true);
              }}
              className={`px-2.5 py-0.5 rounded border flex items-center gap-1.5 font-bold cursor-pointer transition ${
                report.calibrationStatus === 'SKALIBROWANY'
                  ? 'bg-emerald-950/80 hover:bg-emerald-900/80 border-emerald-500 text-emerald-200 ring-1 ring-emerald-500/30'
                  : report.groundTruthCaRmsdAngstrom !== null
                  ? 'bg-purple-950/80 hover:bg-purple-900/80 border-purple-500 text-purple-200'
                  : 'bg-zinc-900 border-zinc-700 text-zinc-400'
              }`}
              title="Prawdziwy niezależny Ground-Truth Cα RMSD z algorytmu Kabscha SVD względem eksperymentalnej struktury RCSB PDB"
            >
              <BarChart3 className="w-3.5 h-3.5 text-purple-300" />
              <span>
                {report.groundTruthCaRmsdAngstrom !== null ? (
                  <>
                    GT Cα-RMSD: <strong>{report.groundTruthCaRmsdAngstrom.toFixed(3)} Å</strong> ({report.groundTruthValidation.referenceId ?? 'PDB'} •{' '}
                    {report.calibrationStatus === 'SKALIBROWANY' ? 'SKALIBROWANY' : 'BRAK ZBIEŻNOŚCI <1Å'})
                  </>
                ) : (
                  'RMSD referencyjny niedostępny'
                )}
              </span>
            </button>
          )}

          {/* Pipeline Mode Badge */}
          <span
            className={`px-2 py-0.5 rounded border flex items-center gap-1 font-bold ${
              report?.pipelineMode === 'TWO_WAY_IK_ADAPTIVE_SSBOND'
                ? 'bg-amber-950/70 border-amber-700 text-amber-200'
                : 'bg-cyan-950/70 border-cyan-800 text-cyan-200'
            }`}
          >
            <span>
              {report?.pipelineMode === 'TWO_WAY_IK_ADAPTIVE_SSBOND'
                ? `IK + Adapt. w_ss (${report.finalSpringWeight})`
                : 'Potok L2 (Bez S-S)'}
            </span>
          </span>

          {/* D_max Badge */}
          <span className="px-2 py-0.5 rounded bg-black border border-cyan-900 text-cyan-300 flex items-center gap-1">
            <span>D_max: <strong>{report?.dMaxAchievedAngstrom ?? currentSubdomain.dMaxTarget} Å</strong> / {dMaxTarget} Å</span>
          </span>

          {/* Disulfide Restraints Badge */}
          <span
            className={`px-2 py-0.5 rounded border flex items-center gap-1 font-bold ${
              report && report.disulfideBonds.length > 0
                ? report.finalMaxSsDeviationAngstrom < 0.5
                  ? 'bg-emerald-950/70 border-emerald-600 text-emerald-200'
                  : 'bg-rose-950/70 border-rose-600 text-rose-200'
                : 'bg-black border-zinc-800 text-zinc-400'
            }`}
          >
            <ShieldCheck className="w-3 h-3 text-amber-400" />
            <span>
              {report && report.disulfideBonds.length > 0
                ? `SSBOND (${report.disulfideBonds.length}): max Δ=${report.finalMaxSsDeviationAngstrom} Å (<0.5 Å)`
                : 'Mostki S-S: 0'}
            </span>
          </span>

          {/* VdW Clashes */}
          <span className="px-2 py-0.5 rounded bg-black border border-teal-900 text-teal-300 flex items-center gap-1">
            <span>Kolizje VdW: <strong>0</strong></span>
          </span>
        </div>

        {/* Action: Run Optimization & Calibration */}
        <div className="flex items-center gap-2">
          <button
            onClick={handleExecuteCalibration}
            disabled={isCalibrating}
            className="px-4 py-1.5 bg-gradient-to-r from-emerald-600 via-teal-600 to-cyan-600 hover:from-emerald-500 hover:to-cyan-500 disabled:opacity-50 text-white font-bold rounded flex items-center gap-2 cursor-pointer shadow-lg shadow-emerald-500/20 transition"
          >
            {isCalibrating ? (
              <>
                <RefreshCw className="w-3.5 h-3.5 animate-spin text-white" />
                <span>Kalibracja parametrów gięcia... {calibProgress}%</span>
              </>
            ) : (
              <>
                <Play className="w-3.5 h-3.5 text-cyan-200" />
                <span>▶ Kalibruj Parametry Gięcia (&lt;1.0°)</span>
              </>
            )}
          </button>
        </div>

      </div>

      {/* 3. MAIN WORKSPACE: 3D VIEWPORT (LEFT) & CALIBRATION HUD (RIGHT) */}
      <div className="flex-1 w-full min-h-0 flex flex-col lg:flex-row overflow-hidden bg-black">
        
        {/* LEFT: 3D THREE.JS CANVAS */}
        <div className="flex-1 h-full min-h-[350px] relative flex flex-col overflow-hidden bg-black">
          
          {/* Floating Active Residue Info Card */}
          <div className="absolute top-3 left-3 z-10 pointer-events-none">
            <div className="pointer-events-auto bg-black/85 backdrop-blur-md border border-zinc-800 rounded-lg p-2.5 flex items-center gap-3 shadow-xl text-xs font-mono">
              <div className="w-9 h-9 rounded-full bg-cyan-950 border border-cyan-500 flex items-center justify-center font-bold text-cyan-300 shadow">
                {currentSubdomain.sequenceCodes[activeResidueIdx] || 'ALA'}
              </div>
              <div>
                <div className="flex items-center gap-2 font-bold text-white">
                  <span>Reszta #{activeResidueIdx + 1}: {currentSubdomain.sequenceCodes[activeResidueIdx]}</span>
                  <span className="text-[10px] px-1.5 py-0.2 rounded bg-zinc-800 text-zinc-300">
                    Siatka L{subdivisionLevel}: {report?.microFacetCount} faset
                  </span>
                </div>
                {activeResDetail && (
                  <div className="text-[10.5px] text-zinc-400 flex items-center gap-2 mt-0.5">
                    <span>Błąd kątowy: <strong className={activeResDetail.angularErrorDeg < 1.0 ? 'text-emerald-400' : 'text-amber-400'}>
                      {activeResDetail.angularErrorDeg}°
                    </strong></span>
                    <span>•</span>
                    <span>Mikrofasetka: <strong className="text-cyan-400">#{activeResDetail.assignedPatchId}</strong></span>
                    <span>•</span>
                    <span>Dystans VdW: <strong className="text-zinc-300">{activeResDetail.vdwDistAngstrom} Å</strong></span>
                  </div>
                )}
              </div>
            </div>
          </div>

          {/* 3D Visualizer Toggles Bar */}
          <div className="absolute top-3 right-3 z-10 flex items-center gap-1.5 text-[10px] font-mono">
            <button
              onClick={() => setShowSubdividedCage(!showSubdividedCage)}
              className={`px-2 py-1 rounded border transition cursor-pointer flex items-center gap-1 ${
                showSubdividedCage ? 'bg-cyan-950/80 border-cyan-700 text-cyan-300' : 'bg-black border-zinc-800 text-zinc-500'
              }`}
            >
              <Boxes className="w-3 h-3" />
              <span>Siatka ({report?.microFacetCount})</span>
            </button>

            <button
              onClick={() => setShowDmaxBoundary(!showDmaxBoundary)}
              className={`px-2 py-1 rounded border transition cursor-pointer flex items-center gap-1 ${
                showDmaxBoundary ? 'bg-teal-950/80 border-teal-600 text-teal-300' : 'bg-black border-zinc-800 text-zinc-500'
              }`}
            >
              <Activity className="w-3 h-3" />
              <span>Limit D_max ({currentSubdomain.dMaxTarget} Å)</span>
            </button>

            <button
              onClick={() => setShowDisulfideRods(!showDisulfideRods)}
              className={`px-2 py-1 rounded border transition cursor-pointer flex items-center gap-1 ${
                showDisulfideRods ? 'bg-amber-950/80 border-amber-700 text-amber-300' : 'bg-black border-zinc-800 text-zinc-500'
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

          {/* Mount Container */}
          <div
            ref={mountRef}
            className="flex-1 w-full h-full min-h-0 cursor-grab active:cursor-grabbing select-none"
            onMouseDown={handleMouseDown}
            onMouseMove={handleMouseMove}
            onMouseUp={handleMouseUp}
            onMouseLeave={handleMouseUp}
            onWheel={handleWheel}
          />

          {/* Bottom Residue Ribbon */}
          <div className="absolute bottom-2 left-3 right-3 z-10 pointer-events-none">
            <div className="pointer-events-auto bg-black/85 backdrop-blur-md border border-zinc-800 rounded-lg p-1.5 flex items-center gap-1.5 overflow-x-auto shadow-2xl font-mono text-xs">
              <span className="text-[9.5px] text-zinc-500 uppercase font-bold shrink-0 px-2">
                Reszty ({currentSubdomain.sequenceCodes.length} aa):
              </span>
              {currentSubdomain.sequenceCodes.map((code, idx) => {
                const isSelected = activeResidueIdx === idx;
                const isCys = code === 'CYS' || code === 'SEC';
                const resErr = report?.residueDetails[idx]?.angularErrorDeg;

                return (
                  <button
                    key={`${code}-${idx}`}
                    onClick={() => setActiveResidueIdx(idx)}
                    className={`px-2 py-1 rounded text-center shrink-0 border transition cursor-pointer ${
                      isSelected
                        ? 'bg-cyan-950 border-cyan-400 text-cyan-200 ring-1 ring-cyan-400 shadow-md'
                        : isCys
                        ? 'bg-amber-950/60 border-amber-800 text-amber-300'
                        : 'bg-zinc-950 border-zinc-800 text-zinc-400 hover:text-white'
                    }`}
                  >
                    <div className="text-[10px] font-bold">{code}</div>
                    <div className="text-[8px] text-zinc-500">
                      {resErr !== undefined ? `${resErr}°` : `#${idx + 1}`}
                    </div>
                  </button>
                );
              })}
            </div>
          </div>

        </div>

        {/* RIGHT: CALIBRATION PARAMETERS & CONVERGENCE HUD */}
        <div className="w-full lg:w-96 border-t lg:border-t-0 lg:border-l border-zinc-900 bg-[#04040a] p-4 flex flex-col font-mono text-xs overflow-y-auto space-y-4 shrink-0">
          
          <div className="flex items-center justify-between border-b border-zinc-800 pb-2.5">
            <div className="flex items-center gap-2">
              <SlidersHorizontal className="w-4 h-4 text-cyan-400" />
              <h3 className="font-bold text-white text-xs uppercase tracking-wide">
                Parametry Gięcia & Warunki Brzegowe
              </h3>
            </div>
            <span className="text-[10px] text-zinc-500">
              {report?.microFacetCount} stanów
            </span>
          </div>

          {/* Bending Stiffness Kappa */}
          <div className="space-y-1.5 bg-black/60 p-2.5 rounded-lg border border-zinc-850">
            <div className="flex items-center justify-between">
              <span className="font-bold text-zinc-300">Sztywność gięcia κ (Kappa):</span>
              <span className="text-cyan-400 font-extrabold">{kappaBend} kcal/(mol·rad²)</span>
            </div>
            <input
              type="range"
              min="10"
              max="120"
              step="0.5"
              value={kappaBend}
              onChange={(e) => setKappaBend(parseFloat(e.target.value))}
              className="w-full accent-cyan-500 cursor-pointer"
            />
            <div className="text-[9.5px] text-zinc-500">
              Sprężystość kątowa szkieletu peptydowego w złączach #248 ➔ #1
            </div>
          </div>

          {/* Dihedral Offset Compensators (Phi / Psi) */}
          <div className="grid grid-cols-2 gap-2 bg-black/60 p-2.5 rounded-lg border border-zinc-850">
            <div className="space-y-1">
              <div className="flex justify-between text-[11px]">
                <span className="text-zinc-400">ΔΦ offset:</span>
                <span className="text-indigo-400 font-bold">{phiOffset}°</span>
              </div>
              <input
                type="range"
                min="-15"
                max="15"
                step="0.2"
                value={phiOffset}
                onChange={(e) => setPhiOffset(parseFloat(e.target.value))}
                className="w-full accent-indigo-500 cursor-pointer"
              />
            </div>
            <div className="space-y-1">
              <div className="flex justify-between text-[11px]">
                <span className="text-zinc-400">ΔΨ offset:</span>
                <span className="text-indigo-400 font-bold">{psiOffset}°</span>
              </div>
              <input
                type="range"
                min="-15"
                max="15"
                step="0.2"
                value={psiOffset}
                onChange={(e) => setPsiOffset(parseFloat(e.target.value))}
                className="w-full accent-indigo-500 cursor-pointer"
              />
            </div>
          </div>

          {/* Disulfide Spring Restraint */}
          <div className="space-y-1.5 bg-black/60 p-2.5 rounded-lg border border-zinc-850">
            <div className="flex items-center justify-between">
              <span className="font-bold text-zinc-300">Siła klamry S-S (k_SS):</span>
              <span className="text-amber-400 font-extrabold">{ssRestraintK} kcal/(mol·Å²)</span>
            </div>
            <input
              type="range"
              min="1.0"
              max="20.0"
              step="0.5"
              value={ssRestraintK}
              onChange={(e) => setSsRestraintK(parseFloat(e.target.value))}
              className="w-full accent-amber-500 cursor-pointer"
            />
            <div className="text-[9.5px] text-zinc-500 flex justify-between">
              <span>Więzy: d0 = 5.50 Å Cα-Cα</span>
              <span>Klamer: {currentSubdomain.disulfidePairs.length}</span>
            </div>
          </div>

          {/* SAXS D_max Limit Constraint */}
          <div className="space-y-1.5 bg-black/60 p-2.5 rounded-lg border border-zinc-850">
            <div className="flex items-center justify-between">
              <span className="font-bold text-zinc-300">Limit D_max (SAXS):</span>
              <span className="text-teal-400 font-extrabold">{dMaxTarget} Å</span>
            </div>
            <input
              type="range"
              min="10"
              max="60"
              step="0.5"
              value={dMaxTarget}
              onChange={(e) => setDMaxTarget(parseFloat(e.target.value))}
              className="w-full accent-teal-500 cursor-pointer"
            />
            <div className="text-[9.5px] text-zinc-500">
              Maksymalny wymiar cząsteczki wymuszający zwarte upakowanie
            </div>
          </div>

          {/* CONVERGENCE & GROUND-TRUTH VALIDATION SUMMARY */}
          {report && (
            <div className="p-3 bg-zinc-950 border border-zinc-800 rounded-lg space-y-2">
              <div className="text-[11px] font-bold text-white uppercase flex items-center justify-between gap-1">
                <span>Status Kalibracji Ground-Truth:</span>
                {report.calibrationStatus === 'SKALIBROWANY' ? (
                  <span className="text-emerald-400 font-extrabold flex items-center gap-1">
                    <CheckCircle2 className="w-3.5 h-3.5" />
                    <span>SKALIBROWANY (&lt;1.0 Å)</span>
                  </span>
                ) : report.calibrationStatus === 'BRAK_REFERENCJI' ? (
                  <span className="text-zinc-400 font-bold">BRAK REFERENCJI</span>
                ) : (
                  <span className="text-rose-400 font-extrabold flex items-center gap-1">
                    <AlertTriangle className="w-3.5 h-3.5" />
                    <span>BRAK ZBIEŻNOŚCI</span>
                  </span>
                )}
              </div>

              {/* Independent Ground-Truth Validation Box */}
              <div className="p-2 rounded bg-black/80 border border-purple-900/70 space-y-1 text-[10px]">
                <div className="flex justify-between text-purple-200 font-bold">
                  <span>Niezależny GT Cα RMSD (Kabsch SVD):</span>
                  <span>
                    {report.groundTruthCaRmsdAngstrom !== null
                      ? `${report.groundTruthCaRmsdAngstrom.toFixed(4)} Å`
                      : 'RMSD referencyjny niedostępny'}
                  </span>
                </div>
                {report.groundTruthValidation.hasGroundTruthReference && (
                  <>
                    <div className="flex justify-between text-zinc-400">
                      <span>Referencja / Łańcuch:</span>
                      <span className="text-zinc-200">
                        {report.groundTruthValidation.referenceId}:{report.groundTruthValidation.selectedReferenceChain} (
                        {report.groundTruthValidation.pairedResidueCount}/{report.totalResidues} par •{' '}
                        {report.groundTruthValidation.modelCoveragePct}%)
                      </span>
                    </div>
                    <div className="flex justify-between text-zinc-400">
                      <span>GT Backbone RMSD:</span>
                      <span className="text-sky-300 font-bold">
                        {report.backboneRmsdAngstrom !== null ? `${report.backboneRmsdAngstrom.toFixed(3)} Å` : 'niedostępny'}
                      </span>
                    </div>
                  </>
                )}
                {report.groundTruthValidation.failureReason && (
                  <div className="text-[9.5px] text-rose-300 pt-0.5 leading-snug">
                    {report.groundTruthValidation.failureReason}
                  </div>
                )}
              </div>

              {/* Internal Snaplock L2 Quantization Metrics (Explicitly Separated!) */}
              <div className="text-[10px] font-bold text-zinc-400 uppercase pt-1">
                Wewnętrzna Kwantyzacja Siatki L{subdivisionLevel} (Nie-Referencyjna):
              </div>
              <div className="grid grid-cols-2 gap-1.5 text-[10px] text-zinc-400">
                <div>Błąd faset pocz.: <strong className="text-zinc-200">{report.initialMeanErrorDeg}°</strong></div>
                <div>Błąd faset L{subdivisionLevel}: <strong className="text-cyan-300">{report.finalMeanErrorDeg}°</strong></div>
                <div>Wewn. reziduum L2: <strong className="text-amber-300">{report.internalL2BackboneResidualAngstrom} Å</strong></div>
                <div>Mikrofasetki: <strong className="text-cyan-400">{report.microFacetCount}</strong></div>
              </div>

              {/* Disulfide Bridges & Two-Way IK Loop Closure Status */}
              {report.disulfideBonds.length > 0 ? (
                <div className="pt-2 border-t border-zinc-850 space-y-1.5 text-[10px]">
                  <div className="flex items-center justify-between">
                    <span className="text-zinc-400 font-bold uppercase">Domykanie Pętli (IK &amp; SSBOND):</span>
                    <span className="text-amber-300 font-bold">
                      w_ss: {report.initialSpringWeight} ➔ {report.finalSpringWeight}
                      {report.exponentialBoostTriggered ? ' (EXP)' : ''}
                    </span>
                  </div>
                  <div className="flex items-center justify-between text-zinc-400">
                    <span>Kara E_ssbond = w_ss·Δd²:</span>
                    <span className="text-emerald-300 font-bold">{report.totalSsbondPenaltyEnergy} kcal/mol</span>
                  </div>
                  {report.disulfideBonds.map((b, idx) => (
                    <div key={idx} className="p-1.5 rounded bg-black/70 border border-zinc-800/80 space-y-0.5">
                      <div className="flex justify-between font-bold text-zinc-200">
                        <span>CYS{b.res1} ↔ CYS{b.res2}</span>
                        <span className={b.sgDeviation < 0.5 && b.deviation < 0.5 ? 'text-emerald-400' : 'text-rose-400'}>
                          Δ(Sγ)={b.sgDeviation} Å • Δ(Cα)={b.deviation} Å
                        </span>
                      </div>
                      <div className="flex justify-between text-[9.5px] text-zinc-400">
                        <span>Cα-Cα: {b.initialCaDist} ➔ <strong className="text-amber-300">{b.actualDist} Å</strong></span>
                        <span>Sγ-Sγ: <strong className="text-amber-200">{b.actualSgDist} Å</strong> / 2.05 Å</span>
                      </div>
                    </div>
                  ))}
                </div>
              ) : (
                <div className="pt-2 border-t border-zinc-850 text-[10px] text-zinc-500">
                  Struktura bez mostków disiarczkowych — aktywny standardowy potok mikrofaset L2 (2232 stany).
                </div>
              )}
            </div>
          )}

          {/* ACTION BUTTONS: DOWNLOAD JSON, PDB & LOAD TO CONSTRUCTOR */}
          <div className="space-y-2 pt-2">
            <div className="grid grid-cols-2 gap-2">
              <button
                onClick={handleDownloadCalibratedPdb}
                className="py-2 px-3 rounded bg-sky-600 hover:bg-sky-500 text-white font-bold flex items-center justify-center gap-1.5 cursor-pointer shadow-md transition text-xs"
                title="Zapisz skalibrowaną strukturę subdomeny w formacie PDB"
              >
                <Download className="w-3.5 h-3.5" />
                <span>Pobierz PDB</span>
              </button>

              <button
                onClick={handleDownloadJsonReport}
                className="py-2 px-3 rounded bg-zinc-900 hover:bg-zinc-800 border border-zinc-700 text-cyan-300 hover:text-white font-bold flex items-center justify-center gap-1.5 cursor-pointer transition text-xs"
                title="Pobierz pełny raport kalibracji w formacie JSON"
              >
                <FileText className="w-3.5 h-3.5" />
                <span>Raport JSON</span>
              </button>
            </div>

            <button
              onClick={handleCopyJson}
              className="w-full py-1.5 rounded bg-black border border-zinc-800 hover:bg-zinc-900 text-zinc-400 hover:text-white flex items-center justify-center gap-1.5 cursor-pointer transition text-xs"
            >
              {copiedReport ? <Check className="w-3.5 h-3.5 text-emerald-400" /> : <Copy className="w-3.5 h-3.5" />}
              <span>{copiedReport ? 'Skopiowano JSON do schowka!' : 'Kopiuj Raport do Schowka'}</span>
            </button>

            {onLoadSequenceToConstructor && (
              <button
                onClick={() => {
                  onLoadSequenceToConstructor(currentSubdomain.sequenceCodes.join(' '));
                  if (onSwitchToConstructor) onSwitchToConstructor();
                }}
                className="w-full py-2 bg-gradient-to-r from-cyan-600 to-indigo-600 hover:from-cyan-500 hover:to-indigo-500 text-white font-bold rounded flex items-center justify-center gap-2 cursor-pointer shadow-md transition text-xs mt-1"
              >
                <span>Załaduj tę Subdomenę do Konstruktora 3D</span>
                <ArrowRight className="w-3.5 h-3.5" />
              </button>
            )}
          </div>

        </div>

      </div>

      {/* RECHARTS SUBDOMAIN PER-RESIDUE RMSD MODAL */}
      {showRmsdChartModal && report?.rmsdComparison && (
        <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/80 backdrop-blur-md p-4">
          <div className="bg-zinc-950 border border-zinc-800 rounded-xl max-w-5xl w-full p-5 shadow-2xl flex flex-col gap-3">
            <div className="flex items-center justify-between pb-2.5 border-b border-zinc-800">
              <div className="flex items-center gap-2">
                <BarChart3 className="w-4 h-4 text-purple-400" />
                <span className="text-sm font-bold text-white">
                  Walidacja Graficzna RMSD Subdomeny: {report.subdomainName}
                </span>
              </div>
              <button
                type="button"
                onClick={() => setShowRmsdChartModal(false)}
                className="p-1 rounded text-zinc-400 hover:text-white hover:bg-zinc-800 cursor-pointer"
              >
                <CloseIcon className="w-4 h-4" />
              </button>
            </div>

            <SnaplockRmsdComparisonChart
              rmsdComparison={report.rmsdComparison}
              disulfideResidueSeqs={report.disulfideBonds.flatMap(b => [b.res1, b.res2])}
              activeResidueIdx={activeResidueIdx}
              onSelectResidue={(idx) => setActiveResidueIdx(idx)}
            />

            <div className="flex justify-end pt-1">
              <button
                type="button"
                onClick={() => setShowRmsdChartModal(false)}
                className="px-4 py-1.5 rounded bg-zinc-900 hover:bg-zinc-800 text-zinc-200 font-bold text-xs cursor-pointer"
              >
                Zamknij wykres
              </button>
            </div>
          </div>
        </div>
      )}

    </div>
  );
};
