import React, { useState, useEffect, useRef } from 'react';
import * as THREE from 'three';
import { 
  Compass, 
  Atom, 
  Sparkles, 
  Layers, 
  RotateCcw, 
  ShieldCheck, 
  Zap, 
  Flame, 
  Info, 
  Play, 
  Square as StopIcon,
  ChevronRight,
  ChevronLeft,
  Maximize2,
  Minimize2,
  Crosshair,
  Eye,
  Sliders,
  Plus,
  RefreshCw,
  Boxes,
  Activity,
  CheckCircle2
} from 'lucide-react';
import { 
  SNAPLOCK_21_BLOCKS, 
  SnaplockBlock, 
  Facet248Detail, 
  computeNaturalFoldedTrajectory 
} from '../snaplock21Blocks';
import {
  GLOBAL_PROTEIN_TARGETS,
  GlobalProteinTargetId
} from '../snaplockGlobalTargets';

interface SkyMapPlanetariumViewProps {
  activeTarget?: GlobalProteinTargetId;
  activeSequence?: string;
  onSelectGlobalTarget?: (targetId: Exclude<GlobalProteinTargetId, 'CUSTOM'>) => void;
  onChainUpdated?: (chainCodes: string[]) => void;
  onSendToConstructor?: (sequence: string) => void;
}

export const SkyMapPlanetariumView: React.FC<SkyMapPlanetariumViewProps> = ({
  activeTarget = 'UBIQUITIN_76',
  activeSequence,
  onSelectGlobalTarget,
  onChainUpdated,
  onSendToConstructor
}) => {
  // Current sequence being assembled (strictly synchronized with global activeTarget!)
  const [assembledChain, setAssembledChain] = useState<string[]>(() => {
    if (activeTarget && activeTarget !== 'CUSTOM' && GLOBAL_PROTEIN_TARGETS[activeTarget]) {
      return GLOBAL_PROTEIN_TARGETS[activeTarget].sequence3Letter;
    }
    return GLOBAL_PROTEIN_TARGETS.UBIQUITIN_76.sequence3Letter;
  });

  useEffect(() => {
    if (activeTarget && activeTarget !== 'CUSTOM' && GLOBAL_PROTEIN_TARGETS[activeTarget]) {
      setAssembledChain(GLOBAL_PROTEIN_TARGETS[activeTarget].sequence3Letter);
      setActiveResidueIndex(0);
    } else if (activeSequence && activeSequence.trim()) {
      const codes = activeSequence.trim().split(/\s+/).filter(Boolean);
      if (codes.length > 0) {
        setAssembledChain(codes);
        setActiveResidueIndex(0);
      }
    }
  }, [activeTarget, activeSequence]);

  // Active inspected block index in chain
  const [activeResidueIndex, setActiveResidueIndex] = useState<number>(0);
  const activeCode = assembledChain[activeResidueIndex] || 'ALA';
  const currentBlock: SnaplockBlock = SNAPLOCK_21_BLOCKS[activeCode] || SNAPLOCK_21_BLOCKS['ALA'];

  // View modes:
  // 'skymap_inside': camera inside [0,0,0] Center of Mass of active amino acid looking 360° out at its 248 facets
  // 'natural_fold_3d': global 3D orbit around the REAL natural folded tertiary protein composed of 248-facet spheres
  const [viewMode, setViewMode] = useState<'skymap_inside' | 'natural_fold_3d'>('skymap_inside');

  // Visualization toggles for 3D fold
  const [showFacetCages, setShowFacetCages] = useState<boolean>(true);
  const [showAtoms, setShowAtoms] = useState<boolean>(true);
  const [showBridges, setShowBridges] = useState<boolean>(true);
  const [showRibbon, setShowRibbon] = useState<boolean>(true);

  // Selected & hovered facet in SkyMap
  const [selectedFacetId, setSelectedFacetId] = useState<number>(1);
  const [hoveredFacetId, setHoveredFacetId] = useState<number | null>(null);

  // Astronomical look-around angles for SkyMap (deg)
  const [yawDeg, setYawDeg] = useState<number>(0);
  const [pitchDeg, setPitchDeg] = useState<number>(0);
  const [fovZoom, setFovZoom] = useState<number>(75);

  // Global orbit camera angles for 3D fold
  const [orbitRotX, setOrbitRotX] = useState<number>(20);
  const [orbitRotY, setOrbitRotY] = useState<number>(45);
  const [orbitDist, setOrbitDist] = useState<number>(38);

  // Canvas & Three.js refs
  const containerRef = useRef<HTMLDivElement | null>(null);
  const rendererRef = useRef<THREE.WebGLRenderer | null>(null);
  const sceneRef = useRef<THREE.Scene | null>(null);
  const cameraRef = useRef<THREE.PerspectiveCamera | null>(null);
  const animationFrameId = useRef<number | null>(null);

  // Scene groups
  const skyMapGroupRef = useRef<THREE.Group | null>(null);
  const naturalFoldGroupRef = useRef<THREE.Group | null>(null);

  // Pointer drag state
  const isDragging = useRef<boolean>(false);
  const prevMousePos = useRef<{ x: number; y: number }>({ x: 0, y: 0 });

  // Web Audio Context for harmonic tones
  const audioCtxRef = useRef<AudioContext | null>(null);

  const playTone = (freqHz: number = 440) => {
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

      osc.type = 'sine';
      osc.frequency.setValueAtTime(freqHz, now);

      gain.gain.setValueAtTime(0, now);
      gain.gain.linearRampToValueAtTime(0.12, now + 0.02);
      gain.gain.exponentialRampToValueAtTime(0.001, now + 0.35);

      osc.connect(gain);
      gain.connect(ctx.destination);
      osc.start(now);
      osc.stop(now + 0.38);
    } catch (e) {}
  };

  // Active inspected facet
  const activeFacet: Facet248Detail = currentBlock.facets.find(f => f.id === (hoveredFacetId || selectedFacetId)) || currentBlock.facets[0];

  // --- THREE.JS INITIALIZATION ---
  useEffect(() => {
    if (!containerRef.current) return;
    const container = containerRef.current;
    const width = container.clientWidth || 800;
    const height = container.clientHeight || 600;

    const scene = new THREE.Scene();
    scene.background = new THREE.Color(0x020206);
    sceneRef.current = scene;

    const camera = new THREE.PerspectiveCamera(fovZoom, width / height, 0.1, 1000);
    cameraRef.current = camera;

    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true, powerPreference: 'high-performance' });
    renderer.setSize(width, height);
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    renderer.outputColorSpace = THREE.SRGBColorSpace;
    container.appendChild(renderer.domElement);
    rendererRef.current = renderer;

    // Lighting
    const ambLight = new THREE.AmbientLight(0xffffff, 0.9);
    scene.add(ambLight);

    const dirLight1 = new THREE.DirectionalLight(0xffffff, 0.8);
    dirLight1.position.set(20, 40, 20);
    scene.add(dirLight1);

    const dirLight2 = new THREE.DirectionalLight(0x38bdf8, 0.5);
    dirLight2.position.set(-20, -30, -20);
    scene.add(dirLight2);

    const skyGroup = new THREE.Group();
    skyMapGroupRef.current = skyGroup;
    scene.add(skyGroup);

    const foldGroup = new THREE.Group();
    naturalFoldGroupRef.current = foldGroup;
    scene.add(foldGroup);

    // Animation Loop
    const animate = () => {
      animationFrameId.current = requestAnimationFrame(animate);

      if (cameraRef.current) {
        if (viewMode === 'skymap_inside') {
          // Camera rigidly in [0, 0, 0] Center of Mass, rotating 360° Pitch/Yaw
          cameraRef.current.position.set(0, 0, 0);
          const phi = THREE.MathUtils.degToRad(90 - pitchDeg);
          const theta = THREE.MathUtils.degToRad(yawDeg);

          const target = new THREE.Vector3(
            Math.sin(phi) * Math.cos(theta),
            Math.cos(phi),
            Math.sin(phi) * Math.sin(theta)
          );
          cameraRef.current.lookAt(target);
        } else {
          // Global 3D Orbit around natural folded protein
          const radX = THREE.MathUtils.degToRad(orbitRotX);
          const radY = THREE.MathUtils.degToRad(orbitRotY);

          const cx = Math.sin(radY) * Math.cos(radX) * orbitDist;
          const cy = Math.sin(radX) * orbitDist;
          const cz = Math.cos(radY) * Math.cos(radX) * orbitDist;

          cameraRef.current.position.set(cx, cy, cz);
          cameraRef.current.lookAt(0, 0, 0);
        }
      }

      renderer.render(scene, cameraRef.current!);
    };
    animate();

    const handleResize = () => {
      if (!containerRef.current || !rendererRef.current || !cameraRef.current) return;
      const w = containerRef.current.clientWidth;
      const h = containerRef.current.clientHeight;
      cameraRef.current.aspect = w / h;
      cameraRef.current.updateProjectionMatrix();
      rendererRef.current.setSize(w, h);
    };
    window.addEventListener('resize', handleResize);

    return () => {
      window.removeEventListener('resize', handleResize);
      if (animationFrameId.current) cancelAnimationFrame(animationFrameId.current);
      if (rendererRef.current?.domElement) {
        rendererRef.current.domElement.remove();
        rendererRef.current.dispose();
      }
    };
  }, [viewMode]);

  // Update FOV on zoom change
  useEffect(() => {
    if (cameraRef.current) {
      cameraRef.current.fov = fovZoom;
      cameraRef.current.updateProjectionMatrix();
    }
  }, [fovZoom]);

  // --- REBUILD SKYMAP 360° PLANETARIUM DOME (INSIDE ACTIVE AMINO ACID) ---
  useEffect(() => {
    if (!skyMapGroupRef.current) return;
    const group = skyMapGroupRef.current;

    while (group.children.length > 0) {
      group.remove(group.children[0]);
    }

    if (viewMode !== 'skymap_inside') return;

    const SPHERE_RADIUS = 13.0; // Radius of internal vector dome around [0,0,0]

    // 1. Celestial wireframe dome
    const domeGeom = new THREE.IcosahedronGeometry(SPHERE_RADIUS, 3);
    const domeMat = new THREE.MeshBasicMaterial({
      color: 0x1e1b4b,
      wireframe: true,
      transparent: true,
      opacity: 0.18,
      side: THREE.BackSide
    });
    group.add(new THREE.Mesh(domeGeom, domeMat));

    // 2. Build 248 TRIANGULAR FACETS on the sphere
    currentBlock.facets.forEach((facet) => {
      const [v1, v2, v3] = facet.vertices;

      const p1 = new THREE.Vector3(...v1).multiplyScalar(SPHERE_RADIUS);
      const p2 = new THREE.Vector3(...v2).multiplyScalar(SPHERE_RADIUS);
      const p3 = new THREE.Vector3(...v3).multiplyScalar(SPHERE_RADIUS);

      const triGeom = new THREE.BufferGeometry();
      const posArray = new Float32Array([
        p1.x, p1.y, p1.z,
        p2.x, p2.y, p2.z,
        p3.x, p3.y, p3.z
      ]);
      triGeom.setAttribute('position', new THREE.BufferAttribute(posArray, 3));
      triGeom.computeVertexNormals();

      // Facet coloring based on Port #1, Port #248, sidechain R or polarity
      let facetColor = 0x1e293b;
      let opacity = 0.35;
      let emissive = 0x000000;

      if (facet.isPort1) {
        // Port #1: Zenit (+Y) N-inlet
        facetColor = 0x10b981; // Bright emerald
        opacity = 0.9;
        emissive = 0x059669;
      } else if (facet.isPort248) {
        // Port #248: Nadir (-Y) C-outlet
        facetColor = 0xf59e0b; // Bright gold
        opacity = 0.9;
        emissive = 0xd97706;
      } else if (facet.isSidechainR) {
        // Sidechain R sector
        facetColor = new THREE.Color(currentBlock.color).getHex();
        opacity = 0.65;
      } else if (facet.isCardinal) {
        facetColor = 0x38bdf8;
        opacity = 0.7;
      } else if (facet.hydrophobicity > 1.5) {
        facetColor = 0xd97706; // Amber
        opacity = 0.45;
      } else if (facet.chargeAtPH74 !== 0) {
        facetColor = facet.chargeAtPH74 > 0 ? 0x8b5cf6 : 0xec4899;
        opacity = 0.55;
      }

      const triMat = new THREE.MeshStandardMaterial({
        color: facetColor,
        emissive,
        emissiveIntensity: facet.isPort1 || facet.isPort248 ? 0.6 : 0.1,
        side: THREE.DoubleSide,
        transparent: true,
        opacity,
        roughness: 0.3
      });

      const facetMesh = new THREE.Mesh(triGeom, triMat);
      facetMesh.userData = { facetId: facet.id };
      group.add(facetMesh);

      // Triangular wireframe border
      const wireGeom = new THREE.BufferGeometry().setFromPoints([p1, p2, p3, p1]);
      const wireMat = new THREE.LineBasicMaterial({
        color: facet.isPort1 ? 0x34d399 : facet.isPort248 ? 0xfbbf24 : 0x334155,
        transparent: true,
        opacity: 0.6
      });
      group.add(new THREE.Line(wireGeom, wireMat));

      // Node point in the center of facet
      const centerPos = new THREE.Vector3(...facet.normal).multiplyScalar(SPHERE_RADIUS);
      const nodeGeom = new THREE.SphereGeometry(facet.isPort1 || facet.isPort248 ? 0.45 : 0.22, 16, 16);
      const nodeMat = new THREE.MeshBasicMaterial({
        color: facet.isPort1 ? 0x10b981 : facet.isPort248 ? 0xf59e0b : facetColor
      });
      const nodeMesh = new THREE.Mesh(nodeGeom, nodeMat);
      nodeMesh.position.copy(centerPos);
      nodeMesh.userData = { facetId: facet.id };
      group.add(nodeMesh);

      // Beacon ring for Port 1 and Port 248
      if (facet.isPort1 || facet.isPort248) {
        const ringGeom = new THREE.TorusGeometry(0.85, 0.06, 12, 32);
        const ringMat = new THREE.MeshBasicMaterial({
          color: facet.isPort1 ? 0x34d399 : 0xfbbf24
        });
        const ringMesh = new THREE.Mesh(ringGeom, ringMat);
        ringMesh.position.copy(centerPos);
        ringMesh.rotation.x = Math.PI / 2;
        group.add(ringMesh);
      }
    });

    // 3. Render atoms of the active amino acid inside the sphere
    const ATOM_SCALE = 2.4;
    const atomPosMap = new Map<string, THREE.Vector3>();

    currentBlock.atoms.forEach((atom) => {
      const pos = new THREE.Vector3(atom.x * ATOM_SCALE, atom.y * ATOM_SCALE, atom.z * ATOM_SCALE);
      atomPosMap.set(atom.id, pos);

      let color = 0x64748b;
      let rad = 0.48;
      if (atom.element === 'N') { color = 0x2563eb; rad = 0.52; }
      else if (atom.element === 'O') { color = 0xdc2626; rad = 0.48; }
      else if (atom.element === 'S') { color = 0xeab308; rad = 0.65; }
      else if (atom.element === 'Se') { color = 0x0d9488; rad = 0.72; }
      else if (atom.element === 'H') { color = 0xf8fafc; rad = 0.28; }
      else if (atom.element === 'C') { color = 0x334155; rad = 0.54; }

      const atomGeom = new THREE.SphereGeometry(rad, 24, 24);
      const atomMat = new THREE.MeshStandardMaterial({
        color,
        roughness: 0.3,
        metalness: 0.2,
        transparent: true,
        opacity: 0.9
      });
      const aMesh = new THREE.Mesh(atomGeom, atomMat);
      aMesh.position.copy(pos);
      group.add(aMesh);
    });

    // 4. Chemical bonds connecting atoms
    currentBlock.bonds.forEach((bond) => {
      const p1 = atomPosMap.get(bond.from);
      const p2 = atomPosMap.get(bond.to);
      if (!p1 || !p2) return;

      const dist = p1.distanceTo(p2);
      const mid = new THREE.Vector3().addVectors(p1, p2).multiplyScalar(0.5);
      const cylGeom = new THREE.CylinderGeometry(0.08, 0.08, dist, 12);
      const cylMat = new THREE.MeshStandardMaterial({ color: 0x64748b, roughness: 0.4 });
      const cyl = new THREE.Mesh(cylGeom, cylMat);
      cyl.position.copy(mid);
      cyl.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), new THREE.Vector3().subVectors(p2, p1).normalize());
      group.add(cyl);
    });

  }, [currentBlock, viewMode]);

  // --- REBUILD GLOBAL NATURAL 3D FOLD (COMPOSED OF 248-FACET SPHERES, NOT A COLUMN!) ---
  useEffect(() => {
    if (!naturalFoldGroupRef.current) return;
    const group = naturalFoldGroupRef.current;

    while (group.children.length > 0) {
      group.remove(group.children[0]);
    }

    if (viewMode !== 'natural_fold_3d') return;

    // Calculate real, natural folded 3D trajectory
    const trajectory = computeNaturalFoldedTrajectory(assembledChain);
    const n = assembledChain.length;

    const blockPositions: THREE.Vector3[] = [];

    for (let i = 0; i < n; i++) {
      const code = assembledChain[i];
      const block = SNAPLOCK_21_BLOCKS[code] || SNAPLOCK_21_BLOCKS['ALA'];
      const pos = new THREE.Vector3(...trajectory.positions[i]);
      const quat = trajectory.orientations[i];
      blockPositions.push(pos);

      const isCurrentActive = i === activeResidueIndex;

      // 1. INDIVIDUAL 248-FACET SPHERE / CAGE FOR THIS AMINO ACID
      if (showFacetCages) {
        const sphereRadius = 1.65; // physical envelope size
        const cageGroup = new THREE.Group();
        cageGroup.position.copy(pos);
        cageGroup.quaternion.copy(quat);

        // Render polyhedron facets
        const geom = new THREE.IcosahedronGeometry(sphereRadius, 2);
        const mat = new THREE.MeshStandardMaterial({
          color: new THREE.Color(block.color),
          wireframe: true,
          transparent: true,
          opacity: isCurrentActive ? 0.95 : 0.45,
          roughness: 0.3
        });
        const cageMesh = new THREE.Mesh(geom, mat);
        cageGroup.add(cageMesh);

        // Core translucent sphere
        const coreMat = new THREE.MeshStandardMaterial({
          color: new THREE.Color(block.color),
          transparent: true,
          opacity: isCurrentActive ? 0.65 : 0.25,
          roughness: 0.4
        });
        const coreMesh = new THREE.Mesh(new THREE.SphereGeometry(sphereRadius * 0.95, 16, 16), coreMat);
        coreMesh.userData = { residueIndex: i, code };
        cageGroup.add(coreMesh);

        // Port #1 (+Y): N-inlet emerald dot
        const p1Mesh = new THREE.Mesh(new THREE.SphereGeometry(0.32, 16, 16), new THREE.MeshBasicMaterial({ color: 0x10b981 }));
        p1Mesh.position.set(0, sphereRadius, 0);
        cageGroup.add(p1Mesh);

        // Port #248 (-Y): C-outlet golden dot
        const p248Mesh = new THREE.Mesh(new THREE.SphereGeometry(0.32, 16, 16), new THREE.MeshBasicMaterial({ color: 0xf59e0b }));
        p248Mesh.position.set(0, -sphereRadius, 0);
        cageGroup.add(p248Mesh);

        group.add(cageGroup);
      }

      // 2. ATOMS INSIDE EACH SPHERE
      if (showAtoms) {
        const atomScale = 0.35;
        block.atoms.forEach((a) => {
          const aPosLocal = new THREE.Vector3(a.x * atomScale, a.y * atomScale, a.z * atomScale);
          aPosLocal.applyQuaternion(quat);
          const aPosWorld = new THREE.Vector3().addVectors(pos, aPosLocal);

          let color = 0x64748b;
          let rad = 0.22;
          if (a.element === 'N') { color = 0x2563eb; rad = 0.25; }
          else if (a.element === 'O') { color = 0xdc2626; rad = 0.24; }
          else if (a.element === 'S') { color = 0xeab308; rad = 0.32; }

          const aMesh = new THREE.Mesh(
            new THREE.SphereGeometry(rad, 12, 12),
            new THREE.MeshStandardMaterial({ color, roughness: 0.3 })
          );
          aMesh.position.copy(aPosWorld);
          group.add(aMesh);
        });
      }

      // 3. PEPTIDE BOND CONNECTION: Port #248 of (i-1) <-> Port #1 of (i)
      if (i > 0) {
        const prevPos = blockPositions[i - 1];
        const lineGeom = new THREE.BufferGeometry().setFromPoints([prevPos, pos]);
        const lineMat = new THREE.LineBasicMaterial({ color: 0x38bdf8, linewidth: 3 });
        group.add(new THREE.Line(lineGeom, lineMat));
      }
    }

    // 4. CONTINUOUS BACKBONE RIBBON
    if (showRibbon && blockPositions.length > 2) {
      const curve = new THREE.CatmullRomCurve3(blockPositions);
      const tubeGeom = new THREE.TubeGeometry(curve, n * 6, 0.4, 8, false);
      const tubeMat = new THREE.MeshStandardMaterial({
        color: 0x4f46e5,
        roughness: 0.3,
        metalness: 0.2,
        transparent: true,
        opacity: 0.85
      });
      group.add(new THREE.Mesh(tubeGeom, tubeMat));
    }

    // 5. DISULFIDE BRIDGES (CYS-CYS COVALENT CLAMPS)
    if (showBridges) {
      const cysIndices: number[] = [];
      assembledChain.forEach((c, idx) => {
        if (c === 'CYS' || c === 'SEC') cysIndices.push(idx);
      });

      // Pair up cysteines (e.g. 6-11 in Insulin A-chain)
      for (let c = 0; c < cysIndices.length; c += 2) {
        if (c + 1 < cysIndices.length) {
          const idx1 = cysIndices[c];
          const idx2 = cysIndices[c + 1];
          const pos1 = blockPositions[idx1];
          const pos2 = blockPositions[idx2];

          const bridgeLine = new THREE.Line(
            new THREE.BufferGeometry().setFromPoints([pos1, pos2]),
            new THREE.LineDashedMaterial({ color: 0x10b981, dashSize: 0.4, gapSize: 0.2 })
          );
          group.add(bridgeLine);

          // Covalent clamp golden badge
          const mid = new THREE.Vector3().addVectors(pos1, pos2).multiplyScalar(0.5);
          const badge = new THREE.Mesh(
            new THREE.SphereGeometry(0.65, 16, 16),
            new THREE.MeshStandardMaterial({ color: 0x059669, emissive: 0x10b981, emissiveIntensity: 0.5 })
          );
          badge.position.copy(mid);
          group.add(badge);
        }
      }
    }

  }, [assembledChain, viewMode, showFacetCages, showAtoms, showBridges, showRibbon, activeResidueIndex]);

  // --- MOUSE LOOK & DRAG HANDLERS ---
  const handleMouseDown = (e: React.MouseEvent) => {
    isDragging.current = true;
    prevMousePos.current = { x: e.clientX, y: e.clientY };
  };

  const handleMouseMove = (e: React.MouseEvent) => {
    if (isDragging.current) {
      const dx = e.clientX - prevMousePos.current.x;
      const dy = e.clientY - prevMousePos.current.y;
      prevMousePos.current = { x: e.clientX, y: e.clientY };

      if (viewMode === 'skymap_inside') {
        const SENSITIVITY = 0.25;
        setYawDeg(prev => (prev - dx * SENSITIVITY + 360) % 360);
        setPitchDeg(prev => Math.max(-88, Math.min(88, prev + dy * SENSITIVITY)));
      } else {
        // Global 3D Orbit
        setOrbitRotY(prev => (prev - dx * 0.5 + 360) % 360);
        setOrbitRotX(prev => Math.max(-85, Math.min(85, prev + dy * 0.5)));
      }
    }

    // Raycast in SkyMap mode to highlight facets
    if (viewMode === 'skymap_inside' && containerRef.current && cameraRef.current && sceneRef.current) {
      const rect = containerRef.current.getBoundingClientRect();
      const mouseX = ((e.clientX - rect.left) / rect.width) * 2 - 1;
      const mouseY = -((e.clientY - rect.top) / rect.height) * 2 + 1;

      const raycaster = new THREE.Raycaster();
      raycaster.setFromCamera(new THREE.Vector2(mouseX, mouseY), cameraRef.current);

      if (skyMapGroupRef.current) {
        const hits = raycaster.intersectObjects(skyMapGroupRef.current.children);
        if (hits.length > 0 && hits[0].object.userData?.facetId) {
          setHoveredFacetId(hits[0].object.userData.facetId);
          return;
        }
      }
      setHoveredFacetId(null);
    }
  };

  const handleMouseUp = () => {
    isDragging.current = false;
  };

  const handleClickCanvas = () => {
    if (viewMode === 'skymap_inside' && hoveredFacetId) {
      setSelectedFacetId(hoveredFacetId);
      playTone(280 + hoveredFacetId * 2.5);
    }
  };

  const handleWheel = (e: React.WheelEvent) => {
    e.preventDefault();
    if (viewMode === 'skymap_inside') {
      setFovZoom(prev => Math.max(30, Math.min(100, prev + e.deltaY * 0.05)));
    } else {
      setOrbitDist(prev => Math.max(12, Math.min(120, prev + e.deltaY * 0.05)));
    }
  };

  // --- CHAIN OPERATIONS ---
  const handleSnapNext = (code: string) => {
    const next = [...assembledChain, code];
    setAssembledChain(next);
    const newIdx = next.length - 1;
    setActiveResidueIndex(newIdx);
    playTone(523.25);
    if (onChainUpdated) onChainUpdated(next);
  };

  const handleStepBack = () => {
    if (assembledChain.length <= 1) return;
    const next = assembledChain.slice(0, -1);
    setAssembledChain(next);
    const newIdx = next.length - 1;
    setActiveResidueIndex(newIdx);
    if (onChainUpdated) onChainUpdated(next);
  };

  const handleLoadPreset = (preset: 'insulin' | 'ubiquitin' | 'crambin' | 'pure_helix') => {
    let targetId: Exclude<GlobalProteinTargetId, 'CUSTOM'> = 'UBIQUITIN_76';
    if (preset === 'insulin') targetId = 'INSULIN_A_21';
    else if (preset === 'ubiquitin') targetId = 'UBIQUITIN_76';
    else if (preset === 'crambin') targetId = 'CRAMBIN_46';
    else if (preset === 'pure_helix') targetId = 'HELIX_BASELINE_18';

    const seq = GLOBAL_PROTEIN_TARGETS[targetId].sequence3Letter;
    setAssembledChain(seq);
    setActiveResidueIndex(0);
    playTone(440);
    if (onSelectGlobalTarget) onSelectGlobalTarget(targetId);
    if (onChainUpdated) onChainUpdated(seq);
  };

  return (
    <div className="flex-1 w-full h-full min-h-0 bg-[#020206] text-slate-100 flex flex-col overflow-hidden font-mono select-none relative">
      
      {/* 1. TOP STATUS & MODE SELECTOR BAR */}
      <div className="bg-[#03040b] border-b border-zinc-900 px-4 py-2 flex flex-wrap items-center justify-between gap-3 shrink-0 z-20">
        
        {/* Title */}
        <div className="flex items-center gap-3">
          <div className="w-9 h-9 rounded-xl bg-indigo-950/80 border border-indigo-500/50 flex items-center justify-center text-cyan-400 shadow-inner">
            <Boxes className="w-5 h-5 text-cyan-300" />
          </div>
          <div>
            <div className="flex items-center gap-2">
              <h1 className="text-sm sm:text-base font-bold text-white tracking-tight flex items-center gap-2">
                SNAPLOCK • 21 Klocków ze Sferami 248 Faset
              </h1>
              <span className="text-[10px] px-2 py-0.5 rounded bg-emerald-950 text-emerald-300 border border-emerald-700/60 font-bold uppercase">
                Naturalny Kształt 3D (Zero Słupka)
              </span>
            </div>
            <p className="text-[11px] text-zinc-400 mt-0.5">
              Każdy aminokwas ma własną kulę 248 trójkątów: Faseta #1 (Wejście N) łączy się sztywno z Fasetą #248 (Wyjście C).
            </p>
          </div>
        </div>

        {/* View Switcher: SkyMap 360° vs Natural 3D Fold */}
        <div className="flex items-center gap-2 flex-wrap text-xs">
          <div className="inline-flex rounded-lg border border-zinc-800 bg-black p-0.5">
            <button
              onClick={() => setViewMode('skymap_inside')}
              className={`px-3.5 py-1.5 rounded-md flex items-center gap-2 font-bold cursor-pointer transition ${
                viewMode === 'skymap_inside'
                  ? 'bg-gradient-to-r from-cyan-600 to-indigo-600 text-white shadow-lg'
                  : 'text-zinc-400 hover:text-white'
              }`}
              title="Kamera w punkcie [0,0,0] Centrum Masy aminokwasu, rozglądasz się 360° wewnątrz sfery 248 faset"
            >
              <Eye className="w-3.5 h-3.5" />
              <span>SkyMap 360° (Wnętrze Kuli)</span>
            </button>

            <button
              onClick={() => setViewMode('natural_fold_3d')}
              className={`px-3.5 py-1.5 rounded-md flex items-center gap-2 font-bold cursor-pointer transition ${
                viewMode === 'natural_fold_3d'
                  ? 'bg-gradient-to-r from-emerald-600 to-teal-600 text-white shadow-lg'
                  : 'text-zinc-400 hover:text-white'
              }`}
              title="Widok 3D całej cząsteczki: naturalny splot globularny złożony z indywidualnych sfer 248 faset"
            >
              <Boxes className="w-3.5 h-3.5" />
              <span>Naturalny Splot 3D (Klocki 248)</span>
            </button>
          </div>

          {/* Quick presets */}
          <div className="hidden xl:flex items-center gap-1 pl-2 border-l border-zinc-800 text-[11px]">
            <span className="text-[10px] text-zinc-500 uppercase mr-1">Wzorce:</span>
            <button
              onClick={() => handleLoadPreset('insulin')}
              className="px-2 py-1 rounded bg-amber-950/70 hover:bg-amber-900 border border-amber-600/60 text-amber-300 font-bold cursor-pointer"
              title="Insulina (21 aa): globularny zwój z 3 mostkami S-S"
            >
              Insulina (21 aa)
            </button>
            <button
              onClick={() => handleLoadPreset('ubiquitin')}
              className="px-2 py-1 rounded bg-indigo-950/70 hover:bg-indigo-900 border border-indigo-600/60 text-indigo-300 font-bold cursor-pointer"
              title="Ubikwityna (76 aa): harmonijka beta + helisa alfa"
            >
              Ubikwityna (76 aa)
            </button>
            <button
              onClick={() => handleLoadPreset('pure_helix')}
              className="px-2 py-1 rounded bg-cyan-950/70 hover:bg-cyan-900 border border-cyan-600/60 text-cyan-300 font-bold cursor-pointer"
              title="Czysta Helisa Alfa (18 aa): skręt 3.6 reszty/obrót"
            >
              Helisa Alfa
            </button>
          </div>

          {onSendToConstructor && (
            <button
              onClick={() => onSendToConstructor(assembledChain.join(' '))}
              className="px-3 py-1.5 bg-gradient-to-r from-emerald-600 to-teal-600 hover:from-emerald-500 text-white rounded-lg font-bold flex items-center gap-1.5 cursor-pointer shadow-md transition ml-1"
            >
              <span>Wyślij do Konstruktora</span>
              <ChevronRight className="w-3.5 h-3.5" />
            </button>
          )}
        </div>
      </div>

      {/* 2. MAIN 3D THREE.JS CANVAS VIEWPORT */}
      <div
        ref={containerRef}
        className="flex-1 w-full h-full min-h-0 relative cursor-grab active:cursor-grabbing overflow-hidden"
        onMouseDown={handleMouseDown}
        onMouseMove={handleMouseMove}
        onMouseUp={handleMouseUp}
        onClick={handleClickCanvas}
        onWheel={handleWheel}
      >
        {/* CROSSHAIR IN SKYMAP INSIDE MODE */}
        {viewMode === 'skymap_inside' && (
          <div className="absolute inset-0 pointer-events-none flex items-center justify-center">
            <div className="relative flex items-center justify-center">
              <div className="w-8 h-8 rounded-full border border-cyan-400/40" />
              <div className="w-2 h-2 rounded-full bg-cyan-400/60" />
              <div className="absolute w-12 h-[1px] bg-cyan-400/30" />
              <div className="absolute h-12 w-[1px] bg-cyan-400/30" />
            </div>
          </div>
        )}

        {/* TOP-LEFT TELEMETRY HUD */}
        <div className="absolute top-4 left-4 z-10 bg-black/80 backdrop-blur-md border border-zinc-800 rounded-xl p-3 text-xs space-y-2 pointer-events-none shadow-xl w-64">
          <div className="flex items-center justify-between border-b border-zinc-800 pb-1.5">
            <span className="text-[10px] text-zinc-400 uppercase font-bold flex items-center gap-1.5">
              <Compass className="w-3.5 h-3.5 text-cyan-400" />
              <span>{viewMode === 'skymap_inside' ? 'SkyMap Telemetria' : 'Splot Naturalny 3D'}</span>
            </span>
            <span className="text-cyan-300 font-extrabold text-xs">
              #{activeResidueIndex + 1}. {currentBlock.code3}
            </span>
          </div>

          {viewMode === 'skymap_inside' ? (
            <div className="grid grid-cols-2 gap-2 text-[11px]">
              <div>
                <span className="text-zinc-500 block text-[9.5px]">Azymut (Yaw):</span>
                <strong className="text-white">{yawDeg.toFixed(1)}°</strong>
              </div>
              <div>
                <span className="text-zinc-500 block text-[9.5px]">Wysokość (Pitch):</span>
                <strong className="text-white">{pitchDeg.toFixed(1)}°</strong>
              </div>
              <div>
                <span className="text-zinc-500 block text-[9.5px]">Środek Masy:</span>
                <strong className="text-emerald-400">[0.0, 0.0, 0.0]</strong>
              </div>
              <div>
                <span className="text-zinc-500 block text-[9.5px]">Pole FOV:</span>
                <strong className="text-cyan-400">{fovZoom}°</strong>
              </div>
            </div>
          ) : (
            <div className="space-y-1 text-[11px]">
              <div className="flex justify-between">
                <span className="text-zinc-500">Długość łańcucha:</span>
                <strong className="text-white">{assembledChain.length} klocków</strong>
              </div>
              <div className="flex justify-between">
                <span className="text-zinc-500">Konformacja:</span>
                <strong className="text-emerald-400">Globularna / Helikalna</strong>
              </div>
              <div className="flex justify-between">
                <span className="text-zinc-500">Klamry kowalencyjne:</span>
                <strong className="text-amber-300">S-S CYS Aktywne</strong>
              </div>
            </div>
          )}

          <div className="pt-1.5 border-t border-zinc-850 text-[10px] text-zinc-400 leading-tight">
            {viewMode === 'skymap_inside'
              ? 'Przeciągnij myszką, by obrócić się o 360° wewnątrz sfery 248 faset.'
              : 'Obracaj model 3D myszką. Kliknij dowolną kulę, by wejść do jej środka masy!'}
          </div>
        </div>

        {/* TOP-RIGHT BIOCHEMICAL DETAILS HUD (SKYMAP MODE) */}
        {viewMode === 'skymap_inside' && (
          <div className="absolute top-4 right-4 z-10 bg-[#040612]/92 backdrop-blur-md border border-cyan-900/80 rounded-2xl p-4 text-xs space-y-3 shadow-2xl w-80 pointer-events-auto">
            
            <div className="flex items-center justify-between border-b border-zinc-800 pb-2">
              <div className="flex items-center gap-2">
                <span className={`w-3 h-3 rounded-full ${
                  activeFacet.isPort1 ? 'bg-emerald-400 animate-ping' : activeFacet.isPort248 ? 'bg-amber-400 animate-ping' : 'bg-cyan-400'
                }`} />
                <h3 className="text-xs font-bold text-white uppercase tracking-wider">
                  Faseta #{activeFacet.id} z 248
                </h3>
              </div>
              <span className={`text-[10px] px-2 py-0.5 rounded font-bold uppercase ${
                activeFacet.isPort1 ? 'bg-emerald-950 text-emerald-300 border border-emerald-700' :
                activeFacet.isPort248 ? 'bg-amber-950 text-amber-300 border border-amber-700' :
                activeFacet.isSidechainR ? 'bg-purple-950 text-purple-300 border border-purple-700' :
                'bg-zinc-900 text-zinc-400 border border-zinc-800'
              }`}>
                {activeFacet.isPort1 ? 'PORT #1 (N)' :
                 activeFacet.isPort248 ? 'PORT #248 (C)' :
                 activeFacet.isSidechainR ? 'ŁAŃCUCH R' : 'FASETA E8'}
              </span>
            </div>

            {/* Normal vector */}
            <div className="bg-black/75 border border-zinc-850 rounded-lg p-2.5 text-[11px] space-y-1">
              <div className="flex justify-between">
                <span className="text-zinc-400">Wektor Normalny:</span>
                <strong className="text-cyan-300 font-mono">
                  [{activeFacet.normal[0].toFixed(2)}, {activeFacet.normal[1].toFixed(2)}, {activeFacet.normal[2].toFixed(2)}]
                </strong>
              </div>
              {activeFacet.cardinalLabel && (
                <p className="text-[10px] text-amber-300 font-bold pt-1 border-t border-zinc-900">
                  {activeFacet.cardinalLabel}
                </p>
              )}
            </div>

            {/* Physicochemical mapping */}
            <div className="space-y-2">
              <span className="text-[10px] text-zinc-400 uppercase font-bold block">
                Cechy Klocka {currentBlock.code3} na tej Fasecie:
              </span>

              <div className="grid grid-cols-2 gap-2 text-[11px]">
                <div className="bg-[#070b1d] border border-zinc-800 p-2 rounded">
                  <span className="text-[9.5px] text-zinc-400 block">Hydropatia (Kyte):</span>
                  <strong className={activeFacet.hydrophobicity >= 0 ? 'text-amber-400' : 'text-cyan-400'}>
                    {activeFacet.hydrophobicity > 0 ? `+${activeFacet.hydrophobicity}` : activeFacet.hydrophobicity}
                  </strong>
                </div>

                <div className="bg-[#070b1d] border border-zinc-800 p-2 rounded">
                  <span className="text-[9.5px] text-zinc-400 block">Ładunek w pH 7.4:</span>
                  <strong className={activeFacet.chargeAtPH74 > 0 ? 'text-purple-400' : activeFacet.chargeAtPH74 < 0 ? 'text-rose-400' : 'text-zinc-300'}>
                    {activeFacet.chargeAtPH74 > 0 ? `+${activeFacet.chargeAtPH74}` : activeFacet.chargeAtPH74}
                  </strong>
                </div>

                <div className="bg-[#070b1d] border border-zinc-800 p-2 rounded">
                  <span className="text-[9.5px] text-zinc-400 block">Wiązania Wodorowe:</span>
                  <strong className="text-emerald-300">
                    {activeFacet.hBondRole}
                  </strong>
                </div>

                <div className="bg-[#070b1d] border border-zinc-800 p-2 rounded">
                  <span className="text-[9.5px] text-zinc-400 block">Promień VdW:</span>
                  <strong className="text-white">
                    {activeFacet.stericRadius} Å
                  </strong>
                </div>
              </div>

              {/* Nearest atom */}
              <div className="bg-black/60 border border-zinc-800 p-2.5 rounded-lg text-[11px] space-y-1">
                <div className="flex justify-between">
                  <span className="text-zinc-400">Najbliższy Atom:</span>
                  <strong className="text-white flex items-center gap-1">
                    <span>{activeFacet.nearestAtom.id}</span>
                    <span className="text-cyan-400 font-extrabold">({activeFacet.nearestAtom.element})</span>
                  </strong>
                </div>
                <div className="flex justify-between">
                  <span className="text-zinc-400">Dystans od środka masy:</span>
                  <span className="text-emerald-400 font-bold">{activeFacet.nearestAtom.distAngstrom} Å</span>
                </div>
              </div>
            </div>

            {/* Quick jump to Port #1 / Port #248 */}
            <div className="pt-2 border-t border-zinc-850 grid grid-cols-2 gap-2">
              <button
                onClick={() => {
                  setPitchDeg(88); setYawDeg(0);
                  setSelectedFacetId(1);
                  playTone(523.25);
                }}
                className="py-1 rounded bg-emerald-950/80 hover:bg-emerald-900 border border-emerald-600 text-emerald-300 font-bold text-center cursor-pointer text-[10.5px]"
              >
                Zenit: Port #1 (N)
              </button>
              <button
                onClick={() => {
                  setPitchDeg(-88); setYawDeg(0);
                  setSelectedFacetId(248);
                  playTone(392);
                }}
                className="py-1 rounded bg-amber-950/80 hover:bg-amber-900 border border-amber-600 text-amber-300 font-bold text-center cursor-pointer text-[10.5px]"
              >
                Nadir: Port #248 (C)
              </button>
            </div>
          </div>
        )}

        {/* TOP-RIGHT VIEWPORT CONTROLS (NATURAL 3D FOLD MODE) */}
        {viewMode === 'natural_fold_3d' && (
          <div className="absolute top-4 right-4 z-10 bg-black/85 backdrop-blur-md border border-zinc-800 rounded-2xl p-3 text-xs space-y-2 shadow-2xl pointer-events-auto w-64">
            <span className="text-[10px] text-zinc-400 uppercase font-bold block border-b border-zinc-800 pb-1">
              Warstwy Wizualne Splotu 3D:
            </span>

            <label className="flex items-center justify-between text-[11px] cursor-pointer hover:text-white">
              <span className="text-zinc-300">Klatki 248 Faset (Sfery):</span>
              <input
                type="checkbox"
                checked={showFacetCages}
                onChange={(e) => setShowFacetCages(e.target.checked)}
                className="accent-cyan-500 cursor-pointer"
              />
            </label>

            <label className="flex items-center justify-between text-[11px] cursor-pointer hover:text-white">
              <span className="text-zinc-300">Atomy wewnątrz kulek:</span>
              <input
                type="checkbox"
                checked={showAtoms}
                onChange={(e) => setShowAtoms(e.target.checked)}
                className="accent-cyan-500 cursor-pointer"
              />
            </label>

            <label className="flex items-center justify-between text-[11px] cursor-pointer hover:text-white">
              <span className="text-zinc-300">Mostki S-S (Klamry CYS):</span>
              <input
                type="checkbox"
                checked={showBridges}
                onChange={(e) => setShowBridges(e.target.checked)}
                className="accent-emerald-500 cursor-pointer"
              />
            </label>

            <label className="flex items-center justify-between text-[11px] cursor-pointer hover:text-white">
              <span className="text-zinc-300">Wstęga szkieletu Cα:</span>
              <input
                type="checkbox"
                checked={showRibbon}
                onChange={(e) => setShowRibbon(e.target.checked)}
                className="accent-indigo-500 cursor-pointer"
              />
            </label>

            <div className="pt-2 border-t border-zinc-850 flex justify-between items-center text-[10.5px]">
              <span className="text-zinc-400">Oddalenie:</span>
              <input
                type="range"
                min="15"
                max="90"
                value={orbitDist}
                onChange={(e) => setOrbitDist(parseInt(e.target.value, 10))}
                className="w-24 accent-teal-500 cursor-pointer"
              />
            </div>
          </div>
        )}

        {/* BOTTOM-LEFT ACTIVE RESIDUE CARD & ASSEMBLY CONTROLS */}
        <div className="absolute bottom-4 left-4 z-10 bg-black/85 backdrop-blur-md border border-zinc-800 rounded-xl p-3 text-xs flex flex-wrap items-center gap-4 shadow-xl">
          <div className="flex items-center gap-3">
            <div
              className="w-10 h-10 rounded-xl flex items-center justify-center font-bold text-white shadow-lg text-sm"
              style={{ backgroundColor: currentBlock.color }}
            >
              {currentBlock.code3}
            </div>
            <div>
              <div className="flex items-center gap-2">
                <span className="text-white font-extrabold text-sm">{currentBlock.namePl}</span>
                <span className="text-[10px] px-1.5 py-0.2 rounded bg-zinc-800 text-zinc-300">
                  {currentBlock.sidechainFormula}
                </span>
              </div>
              <p className="text-[11px] text-zinc-400 mt-0.5">
                Kloczek <strong className="text-cyan-300">#{activeResidueIndex + 1}</strong> z {assembledChain.length} • Wiązanie C-N: <strong className="text-emerald-400">1.33 Å</strong> (Zamrożone)
              </p>
            </div>
          </div>

          <div className="flex items-center gap-2 pl-3 border-l border-zinc-800">
            <button
              onClick={() => handleSnapNext(currentBlock.code3)}
              className="px-3.5 py-1.5 bg-gradient-to-r from-emerald-600 to-teal-600 hover:from-emerald-500 text-white font-bold rounded-lg text-xs flex items-center gap-1.5 cursor-pointer shadow-md transition"
              title="Zatrzaśnij kolejny kloczek: Port #1 następnika łączy się sztywno z Portem #248 poprzednika"
            >
              <Plus className="w-3.5 h-3.5" />
              <span>Zatrzaśnij Kloczek (+1)</span>
            </button>

            <button
              onClick={handleStepBack}
              disabled={assembledChain.length <= 1}
              className="p-1.5 bg-zinc-900 hover:bg-zinc-800 disabled:opacity-30 border border-zinc-700 text-zinc-300 hover:text-white rounded-lg cursor-pointer"
              title="Cofnij ostatni kloczek"
            >
              <RotateCcw className="w-3.5 h-3.5" />
            </button>

            <button
              onClick={() => {
                setAssembledChain(['GLY']);
                setActiveResidueIndex(0);
                if (onChainUpdated) onChainUpdated(['GLY']);
              }}
              className="p-1.5 bg-zinc-900 hover:bg-zinc-800 border border-zinc-700 text-red-400 hover:text-white rounded-lg cursor-pointer"
              title="Resetuj łańcuch do pierwszego klocka"
            >
              <RefreshCw className="w-3.5 h-3.5" />
            </button>
          </div>
        </div>

      </div>

      {/* 3. WAREHOUSE SELECTION RIBBON (21 KLOCKÓW SNAPLOCK) */}
      <div className="bg-[#030308] border-t border-zinc-900 p-2.5 shrink-0 z-20 flex flex-col gap-1.5 font-mono select-none">
        
        <div className="flex justify-between items-center px-1 text-[11px]">
          <div className="flex items-center gap-2">
            <span className="w-2 h-2 rounded-full bg-cyan-400 animate-pulse" />
            <span className="text-white font-extrabold uppercase tracking-wide">
              Magazyn 21 Klocków SNAPLOCK (Kliknij kloczek, by przejść do jego sfery 248):
            </span>
          </div>
          <div className="flex items-center gap-2 text-[10.5px] text-zinc-400">
            <span>Łańcuch białkowy: <strong className="text-emerald-400">{assembledChain.length} reszt</strong></span>
            <span>•</span>
            <span className="text-amber-300">Port #248 ➔ Port #1 (1.33 Å)</span>
          </div>
        </div>

        {/* Carousel of 21 Amino Acids */}
        <div className="flex items-center gap-1.5 overflow-x-auto pb-1 scrollbar-thin scrollbar-thumb-zinc-800">
          {Object.values(SNAPLOCK_21_BLOCKS).map((blk) => {
            const isSelected = activeCode === blk.code3;
            return (
              <button
                key={blk.id}
                onClick={() => {
                  // Switch active inspected block
                  const curChain = [...assembledChain];
                  curChain[activeResidueIndex] = blk.code3;
                  setAssembledChain(curChain);
                  playTone(260 + blk.mass * 0.8);
                  if (onChainUpdated) onChainUpdated(curChain);
                }}
                className={`px-2.5 py-1.5 rounded-lg border text-left flex items-center gap-2 shrink-0 cursor-pointer transition-all ${
                  isSelected
                    ? 'bg-gradient-to-r from-indigo-950 to-cyan-950 border-cyan-400 ring-2 ring-cyan-500/30 scale-[1.03] shadow-lg'
                    : 'bg-black border-zinc-900 hover:bg-zinc-900/80 hover:border-zinc-750 text-zinc-400'
                }`}
                title={`${blk.namePl} (${blk.nameEn}) • Masa: ${blk.mass} g/mol`}
              >
                <span className="w-2.5 h-2.5 rounded-full shrink-0" style={{ backgroundColor: blk.color }} />
                <div className="leading-none">
                  <span className="text-xs font-black text-white">{blk.code3}</span>
                  <span className="text-[8.5px] text-zinc-500 block">{blk.sidechainFormula}</span>
                </div>
              </button>
            );
          })}
        </div>

        {/* Chain Sequence Navigation Strip */}
        <div className="flex items-center gap-1 overflow-x-auto pt-1 border-t border-zinc-900/80 scrollbar-none text-[10.5px]">
          <span className="text-zinc-500 shrink-0 uppercase font-bold mr-1">Sekwencja Łańcucha:</span>
          {assembledChain.map((code, idx) => {
            const isCur = idx === activeResidueIndex;
            const b = SNAPLOCK_21_BLOCKS[code] || SNAPLOCK_21_BLOCKS['ALA'];
            return (
              <button
                key={idx}
                onClick={() => {
                  setActiveResidueIndex(idx);
                  playTone(300 + idx * 10);
                }}
                className={`px-2 py-0.5 rounded border flex items-center gap-1 shrink-0 cursor-pointer transition ${
                  isCur
                    ? 'bg-indigo-950 border-indigo-500 text-white font-extrabold shadow-sm'
                    : 'bg-zinc-950 border-zinc-850 text-zinc-400 hover:text-white'
                }`}
              >
                <span className="text-[9px] text-zinc-500">{idx + 1}.</span>
                <span className="w-1.5 h-1.5 rounded-full" style={{ backgroundColor: b.color }} />
                <span>{b.code3}</span>
              </button>
            );
          })}
        </div>

      </div>

    </div>
  );
};
