import React, { useState, useEffect, useRef } from 'react';
import * as THREE from 'three';
import { 
  Boxes, 
  Atom, 
  RotateCcw, 
  ZoomIn, 
  ZoomOut, 
  Maximize2, 
  Eye, 
  Layers, 
  ShieldCheck, 
  Zap, 
  Flame, 
  Sparkles, 
  Info, 
  Compass, 
  ChevronRight, 
  Search, 
  CheckCircle2,
  ArrowRight,
  Table,
  Check,
  Scale
} from 'lucide-react';
import { 
  SNAPLOCK_21_BLOCKS, 
  SnaplockBlock, 
  Facet248Detail, 
  SNAPLOCK_22_DATABASE 
} from '../snaplock21Blocks';

interface AminoAcidSphereWarehouseViewProps {
  onSelectAminoAcidForConstructor?: (code3: string) => void;
  onSwitchToConstructor?: () => void;
}

export const AminoAcidSphereWarehouseView: React.FC<AminoAcidSphereWarehouseViewProps> = ({
  onSelectAminoAcidForConstructor,
  onSwitchToConstructor
}) => {
  // Currently selected amino acid (default: CYS #12 or ALA #2)
  const [selectedCode, setSelectedCode] = useState<string>('CYS');
  const activeBlock: SnaplockBlock = SNAPLOCK_21_BLOCKS[selectedCode] || SNAPLOCK_21_BLOCKS['ALA'];

  // Inspector mode: 'FACET' (badanie właściwości na trójkąciku kuli) vs 'TABLE_22' (pełna tabela 22 parametrów)
  const [inspectorTab, setInspectorTab] = useState<'FACET' | 'TABLE_22'>('FACET');

  // Category filter for the 21 balls
  const [categoryFilter, setCategoryFilter] = useState<'ALL' | 'HYDROPHOBIC' | 'POLAR' | 'CHARGED' | 'COVALENT'>('ALL');
  const [searchQuery, setSearchQuery] = useState<string>('');

  // Selected & Hovered Facet on the 248-facet sphere
  const [selectedFacetId, setSelectedFacetId] = useState<number>(1);
  const [hoveredFacetId, setHoveredFacetId] = useState<number | null>(null);

  // Visualization toggles
  const [sphereOpacity, setSphereOpacity] = useState<number>(0.38);
  const [showAtoms, setShowAtoms] = useState<boolean>(true);
  const [showPortBeacons, setShowPortBeacons] = useState<boolean>(true);
  const [autoRotate, setAutoRotate] = useState<boolean>(false);
  const [zoomLevel, setZoomLevel] = useState<number>(1.0);

  // Canvas & Three.js refs
  const mountRef = useRef<HTMLDivElement | null>(null);
  const rendererRef = useRef<THREE.WebGLRenderer | null>(null);
  const sceneRef = useRef<THREE.Scene | null>(null);
  const cameraRef = useRef<THREE.PerspectiveCamera | null>(null);
  const sphereGroupRef = useRef<THREE.Group | null>(null);
  const facetMeshRef = useRef<THREE.Mesh | null>(null);
  const animFrameRef = useRef<number | null>(null);

  // Drag interaction
  const isDragging = useRef<boolean>(false);
  const prevMousePos = useRef<{ x: number; y: number }>({ x: 0, y: 0 });
  const rotAngle = useRef<{ x: number; y: number }>({ x: 0.2, y: 0.4 });

  // Web Audio Context for harmonic tones
  const audioCtxRef = useRef<AudioContext | null>(null);
  const playPatchTone = (patchId: number) => {
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

      const freq = 200 + (patchId / 248) * 600;
      osc.type = 'sine';
      osc.frequency.setValueAtTime(freq, now);

      gain.gain.setValueAtTime(0, now);
      gain.gain.linearRampToValueAtTime(0.08, now + 0.02);
      gain.gain.exponentialRampToValueAtTime(0.001, now + 0.22);

      osc.connect(gain);
      gain.connect(ctx.destination);
      osc.start(now);
      osc.stop(now + 0.25);
    } catch (e) {}
  };

  // List of all 21 amino acids in canonical order
  const all21List = Object.values(SNAPLOCK_21_BLOCKS);

  // Filtered list
  const filteredList = all21List.filter(aa => {
    const p = aa.physico;
    const matchesSearch = 
      p.code3.toLowerCase().includes(searchQuery.toLowerCase()) || 
      p.namePl.toLowerCase().includes(searchQuery.toLowerCase()) ||
      p.code1.toLowerCase().includes(searchQuery.toLowerCase());
    if (!matchesSearch) return false;

    if (categoryFilter === 'HYDROPHOBIC') return p.hydropathyKyte > 1.5;
    if (categoryFilter === 'POLAR') return p.polarity.toLowerCase().includes('polar') && p.chargePH74 === 0;
    if (categoryFilter === 'CHARGED') return p.chargePH74 !== 0;
    if (categoryFilter === 'COVALENT') return p.covalentBridges.includes('TAK');
    return true;
  });

  const activeFacet: Facet248Detail = activeBlock.facets.find(f => f.id === (hoveredFacetId || selectedFacetId)) || activeBlock.facets[0];

  // --- THREE.JS INITIALIZATION FOR SELECTED SPHERE ---
  useEffect(() => {
    if (!mountRef.current) return;
    const container = mountRef.current;
    const width = container.clientWidth || 600;
    const height = container.clientHeight || 550;

    const scene = new THREE.Scene();
    scene.background = new THREE.Color(0x020206);
    sceneRef.current = scene;

    const camera = new THREE.PerspectiveCamera(45, width / height, 0.1, 150);
    const initialCamDist = Math.max(9, activeBlock.envelopeRadius * 2.7) / zoomLevel;
    camera.position.set(0, 0, initialCamDist);
    cameraRef.current = camera;

    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    renderer.setSize(width, height);
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    rendererRef.current = renderer;
    container.innerHTML = '';
    container.appendChild(renderer.domElement);

    // Lights
    const ambientLight = new THREE.AmbientLight(0xffffff, 0.8);
    scene.add(ambientLight);

    const dirLight1 = new THREE.DirectionalLight(0xffffff, 1.3);
    dirLight1.position.set(15, 25, 20);
    scene.add(dirLight1);

    const dirLight2 = new THREE.DirectionalLight(0x38bdf8, 0.7);
    dirLight2.position.set(-15, -15, -15);
    scene.add(dirLight2);

    const dirLight3 = new THREE.DirectionalLight(0xf59e0b, 0.4);
    dirLight3.position.set(0, -20, 15);
    scene.add(dirLight3);

    // Main Sphere Group
    const sphereGroup = new THREE.Group();
    sphereGroupRef.current = sphereGroup;
    scene.add(sphereGroup);

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

    // Animation Loop
    const animate = () => {
      animFrameRef.current = requestAnimationFrame(animate);
      if (sphereGroupRef.current) {
        if (autoRotate && !isDragging.current) {
          rotAngle.current.y += 0.005;
        }
        sphereGroupRef.current.rotation.x = rotAngle.current.x;
        sphereGroupRef.current.rotation.y = rotAngle.current.y;
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

  // Update Camera distance dynamically based on activeBlock size and zoomLevel
  useEffect(() => {
    if (cameraRef.current) {
      const dist = Math.max(9, activeBlock.envelopeRadius * 2.7) / zoomLevel;
      cameraRef.current.position.set(0, 0, dist);
    }
  }, [selectedCode, zoomLevel]);

  // --- REBUILD 3D SPHERE AND ATOMS ---
  // PROMIEN SFERY JEST DOKŁADNIE DOPASOWANY DO FIZYCZNEJ POWŁOKI VDW (ENVELOPE RADIUS)
  // ZERO ATOMÓW NIE WYSTAJE POZA KULĘ! 100% W ŚRODKU!
  useEffect(() => {
    if (!sphereGroupRef.current) return;
    const group = sphereGroupRef.current;
    
    // Clear old children
    while (group.children.length > 0) {
      const child = group.children[0];
      group.remove(child);
      if ((child as any).geometry) (child as any).geometry.dispose();
      if ((child as any).material) {
        if (Array.isArray((child as any).material)) {
          (child as any).material.forEach((m: any) => m.dispose());
        } else {
          (child as any).material.dispose();
        }
      }
    }

    // Promień sfery 248 trójkątów: ściśle równy maksymalnemu zasięgowi chmury VdW tego aminokwasu!
    const SPHERE_RADIUS = activeBlock.envelopeRadius;

    // 1. Build the 248 Triangular Facets Mesh
    const positions: number[] = [];
    const colors: number[] = [];

    activeBlock.facets.forEach((f) => {
      const tri = f.vertices;
      tri.forEach(v => {
        positions.push(v[0] * SPHERE_RADIUS, v[1] * SPHERE_RADIUS, v[2] * SPHERE_RADIUS);
      });

      // Kolor faset wg orientacji biochemicznej
      let col = new THREE.Color(0x1a263d);
      if (f.isPort1) {
        col = new THREE.Color(0x06b6d4); // Bright Cyan for Port #1 (N-inlet)
      } else if (f.isPort248) {
        col = new THREE.Color(0xf43f5e); // Bright Rose for Port #248 (C-outlet)
      } else if (f.id === selectedFacetId) {
        col = new THREE.Color(0xfacc15); // Selected yellow highlight
      } else if (f.isSidechainR) {
        if (f.chargeAtPH74 > 0) col = new THREE.Color(0x3b82f6); // Basic (+) Blue
        else if (f.chargeAtPH74 < 0) col = new THREE.Color(0xef4444); // Acidic (-) Red
        else if (f.hydrophobicity > 2.0) col = new THREE.Color(0xf59e0b); // Hydrophobic Amber
        else col = new THREE.Color(0x10b981); // Polar Green
      } else if (f.isCardinal) {
        col = new THREE.Color(0xa855f7); // Purple for Cardinal (#50, #150, #100, #200)
      }

      for (let k = 0; k < 3; k++) {
        colors.push(col.r, col.g, col.b);
      }
    });

    const facetGeometry = new THREE.BufferGeometry();
    facetGeometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
    facetGeometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
    facetGeometry.computeVertexNormals();

    const facetMaterial = new THREE.MeshStandardMaterial({
      vertexColors: true,
      side: THREE.DoubleSide,
      transparent: true,
      opacity: sphereOpacity,
      roughness: 0.35,
      metalness: 0.15
    });

    const facetMesh = new THREE.Mesh(facetGeometry, facetMaterial);
    facetMeshRef.current = facetMesh;
    group.add(facetMesh);

    // 2. Facet Wireframe outline (248 triangular borders)
    const wireframeGeo = new THREE.WireframeGeometry(facetGeometry);
    const wireframeMat = new THREE.LineBasicMaterial({
      color: 0x475569,
      transparent: true,
      opacity: 0.4
    });
    const wireframeMesh = new THREE.LineSegments(wireframeGeo, wireframeMat);
    group.add(wireframeMesh);

    // 3. Highlight on Selected Facet with Vector Line from CoM
    const activeTri = activeFacet.vertices;
    const activeCenter = new THREE.Vector3(
      (activeTri[0][0] + activeTri[1][0] + activeTri[2][0]) / 3 * SPHERE_RADIUS,
      (activeTri[0][1] + activeTri[1][1] + activeTri[2][1]) / 3 * SPHERE_RADIUS,
      (activeTri[0][2] + activeTri[1][2] + activeTri[2][2]) / 3 * SPHERE_RADIUS
    );

    // Vector Ray from [0,0,0] to facet surface
    const rayGeo = new THREE.BufferGeometry().setFromPoints([
      new THREE.Vector3(0, 0, 0),
      activeCenter.clone().multiplyScalar(1.15)
    ]);
    const rayMat = new THREE.LineBasicMaterial({ color: 0xfacc15, linewidth: 2 });
    const rayLine = new THREE.Line(rayGeo, rayMat);
    group.add(rayLine);

    // Selected Facet Ring Beacon
    const ringMesh = new THREE.Mesh(
      new THREE.RingGeometry(0.2, 0.35, 16),
      new THREE.MeshBasicMaterial({ color: 0xfacc15, side: THREE.DoubleSide })
    );
    ringMesh.position.copy(activeCenter.clone().multiplyScalar(1.02));
    ringMesh.quaternion.setFromUnitVectors(new THREE.Vector3(0, 0, 1), new THREE.Vector3(...activeFacet.normal));
    group.add(ringMesh);

    // 4. Highlighted Port #1 (+Y Zenith) and Port #248 (-Y Nadir) Beacon Rings
    if (showPortBeacons) {
      // Port #1 (+Y) North Pole Beacon (N-terminus peptide bond inlet)
      const p1Mesh = new THREE.Mesh(
        new THREE.TorusGeometry(0.48, 0.08, 16, 32),
        new THREE.MeshBasicMaterial({ color: 0x06b6d4 })
      );
      p1Mesh.position.set(0, SPHERE_RADIUS + 0.05, 0);
      p1Mesh.rotation.x = Math.PI / 2;
      group.add(p1Mesh);

      // Port #248 (-Y) South Pole Beacon (C-terminus peptide bond outlet)
      const p248Mesh = new THREE.Mesh(
        new THREE.TorusGeometry(0.48, 0.08, 16, 32),
        new THREE.MeshBasicMaterial({ color: 0xf43f5e })
      );
      p248Mesh.position.set(0, -SPHERE_RADIUS - 0.05, 0);
      p248Mesh.rotation.x = Math.PI / 2;
      group.add(p248Mesh);

      // Cardinal Beacons for #100 (+X), #200 (-X), #50 (+Z), #150 (-Z)
      const cardinalMarkers: { pos: THREE.Vector3; color: number }[] = [
        { pos: new THREE.Vector3(SPHERE_RADIUS, 0, 0), color: 0xef4444 }, // +X #100
        { pos: new THREE.Vector3(-SPHERE_RADIUS, 0, 0), color: 0xef4444 }, // -X #200
        { pos: new THREE.Vector3(0, 0, SPHERE_RADIUS), color: 0x3b82f6 }, // +Z #50
        { pos: new THREE.Vector3(0, 0, -SPHERE_RADIUS), color: 0x3b82f6 }  // -Z #150
      ];

      cardinalMarkers.forEach(cm => {
        const dot = new THREE.Mesh(
          new THREE.SphereGeometry(0.12, 12, 12),
          new THREE.MeshBasicMaterial({ color: cm.color })
        );
        dot.position.copy(cm.pos);
        group.add(dot);
      });
    }

    // 5. Center of Mass [0,0,0] Indicator (Golden Gizmo)
    const comMesh = new THREE.Mesh(
      new THREE.SphereGeometry(0.25, 16, 16),
      new THREE.MeshBasicMaterial({ color: 0xfacc15, wireframe: true })
    );
    group.add(comMesh);

    // 6. Inside: 3D Ball-and-Stick Structural Model of the Amino Acid
    // Wszystkie współrzędne leżą w 100% wewnątrz kuli envelopeRadius
    if (showAtoms && activeBlock.atoms) {
      const atomGroup = new THREE.Group();
      const atomMap: Record<string, THREE.Vector3> = {};

      activeBlock.atoms.forEach((atom) => {
        let atomColor = 0x94a3b8; // C
        let atomRadius = 0.38;

        if (atom.element === 'N') {
          atomColor = 0x3b82f6; // N Blue
          atomRadius = 0.40;
        } else if (atom.element === 'O') {
          atomColor = 0xef4444; // O Red
          atomRadius = 0.38;
        } else if (atom.element === 'S') {
          atomColor = 0xeab308; // S Yellow
          atomRadius = 0.50;
        } else if (atom.element === 'Se') {
          atomColor = 0xf97316; // Se Orange
          atomRadius = 0.55;
        } else if (atom.element === 'H') {
          atomColor = 0xffffff; // H White
          atomRadius = 0.24;
        }

        const atomGeo = new THREE.SphereGeometry(atomRadius, 24, 24);
        const atomMat = new THREE.MeshStandardMaterial({
          color: atomColor,
          roughness: 0.25,
          metalness: 0.2
        });

        const atomMesh = new THREE.Mesh(atomGeo, atomMat);
        const atomPos = new THREE.Vector3(atom.x, atom.y, atom.z);
        atomMesh.position.copy(atomPos);
        atomGroup.add(atomMesh);

        atomMap[atom.id] = atomPos;
      });

      // Bonds as Cylinders
      if (activeBlock.bonds) {
        activeBlock.bonds.forEach((bond) => {
          const p1 = atomMap[bond.from];
          const p2 = atomMap[bond.to];
          if (!p1 || !p2) return;

          const dist = p1.distanceTo(p2);
          const bondGeo = new THREE.CylinderGeometry(0.10, 0.10, dist, 12);
          const bondMat = new THREE.MeshStandardMaterial({
            color: 0x64748b,
            roughness: 0.4,
            metalness: 0.1
          });

          const bondMesh = new THREE.Mesh(bondGeo, bondMat);
          bondMesh.position.copy(p1.clone().add(p2).multiplyScalar(0.5));
          bondMesh.quaternion.setFromUnitVectors(
            new THREE.Vector3(0, 1, 0),
            p2.clone().sub(p1).normalize()
          );

          atomGroup.add(bondMesh);
        });
      }

      group.add(atomGroup);
    }
  }, [selectedCode, sphereOpacity, showAtoms, showPortBeacons, selectedFacetId]);

  // Pointer drag event handlers for rotating the sphere
  const handleMouseDown = (e: React.MouseEvent) => {
    isDragging.current = true;
    prevMousePos.current = { x: e.clientX, y: e.clientY };
  };

  const handleMouseMove = (e: React.MouseEvent) => {
    if (!isDragging.current) return;
    const dx = e.clientX - prevMousePos.current.x;
    const dy = e.clientY - prevMousePos.current.y;

    rotAngle.current.y += dx * 0.008;
    rotAngle.current.x += dy * 0.008;
    rotAngle.current.x = Math.max(-Math.PI / 2 + 0.05, Math.min(Math.PI / 2 - 0.05, rotAngle.current.x));

    prevMousePos.current = { x: e.clientX, y: e.clientY };
  };

  const handleMouseUp = () => {
    isDragging.current = false;
  };

  const handleWheel = (e: React.WheelEvent) => {
    e.preventDefault();
    const zoomDelta = e.deltaY * -0.0015;
    setZoomLevel(prev => Math.max(0.4, Math.min(2.5, prev + zoomDelta)));
  };

  return (
    <div className="flex-1 w-full h-full min-h-0 bg-[#020206] text-white flex flex-col overflow-hidden font-sans select-none">
      
      {/* 1. TOP HEADER & CATEGORY TABS */}
      <div className="bg-[#05060d] border-b border-zinc-800/80 px-4 py-2 flex flex-wrap items-center justify-between gap-3 shrink-0">
        <div className="flex items-center gap-3">
          <div className="flex items-center gap-2">
            <span className="w-2.5 h-2.5 rounded-full bg-cyan-400 animate-pulse" />
            <h1 className="text-sm font-bold tracking-wide uppercase font-mono text-zinc-100 flex items-center gap-1.5">
              <span>Magazyn 21 Kul E8</span>
              <span className="text-xs px-2 py-0.5 rounded bg-cyan-950/80 border border-cyan-800/60 text-cyan-300 font-mono">
                248 Łatek / Kula
              </span>
            </h1>
          </div>
          <span className="text-xs text-zinc-400 hidden lg:inline">
            Środek ciężkości każdego aminokwasu leży ściśle w [0,0,0], a promień kuli R_max zabezpiecza przed wystawaniem atomów
          </span>
        </div>

        {/* Filter categories */}
        <div className="flex items-center gap-1.5 flex-wrap font-mono text-xs">
          <span className="text-[10px] text-zinc-500 uppercase font-bold mr-1">Filtr:</span>
          {(['ALL', 'HYDROPHOBIC', 'POLAR', 'CHARGED', 'COVALENT'] as const).map(cat => (
            <button
              key={cat}
              onClick={() => setCategoryFilter(cat)}
              className={`px-2.5 py-1 rounded text-[11px] font-bold cursor-pointer transition ${
                categoryFilter === cat
                  ? 'bg-cyan-500 text-black shadow-sm'
                  : 'bg-zinc-900 text-zinc-400 hover:text-white hover:bg-zinc-800'
              }`}
            >
              {cat === 'ALL' ? 'Wszystkie (21)' :
               cat === 'HYDROPHOBIC' ? 'Hydrofobowe' :
               cat === 'POLAR' ? 'Polarne' :
               cat === 'CHARGED' ? 'Naładowane' : 'Kowalencyjne S-S/Se'}
            </button>
          ))}
        </div>
      </div>

      {/* 2. MAIN 2-COLUMN WORKSPACE: LEFT = 3D SPHERE STAGE, RIGHT = FACET & PHYSICOCHEMICAL INSPECTOR */}
      <div className="flex-1 w-full min-h-0 flex flex-col lg:flex-row overflow-hidden relative">
        
        {/* LEFT COLUMN: 3D THREE.JS CANVAS OF THE INDIVIDUAL 248-FACET SPHERE */}
        <div className="flex-1 h-full min-h-0 relative flex flex-col bg-black overflow-hidden">
          
          {/* Top Floating Controls Overlay */}
          <div className="absolute top-3 left-3 right-3 z-10 flex items-center justify-between pointer-events-none">
            {/* Active Amino Acid Dimensional Badge */}
            <div className="pointer-events-auto bg-black/85 backdrop-blur-md border border-zinc-800 rounded-lg p-2.5 flex items-center gap-3 shadow-xl">
              <div 
                className="w-11 h-11 rounded-full flex items-center justify-center font-mono font-extrabold text-base border-2 shadow-md shrink-0"
                style={{ 
                  backgroundColor: `${activeBlock.color}25`, 
                  borderColor: activeBlock.color,
                  color: activeBlock.color 
                }}
              >
                {activeBlock.physico.code3}
              </div>
              <div>
                <div className="flex items-center gap-2">
                  <span className="text-sm font-bold text-white font-mono">{activeBlock.physico.namePl}</span>
                  <span className="text-[10px] px-1.5 py-0.2 rounded bg-zinc-800 text-zinc-300 font-mono">
                    {activeBlock.physico.code1}
                  </span>
                  <span className="text-[10px] px-1.5 py-0.2 rounded bg-emerald-950 text-emerald-300 border border-emerald-800/60 font-mono font-bold flex items-center gap-1">
                    <Check className="w-3 h-3 text-emerald-400" />
                    <span>CoM: [0.0, 0.0, 0.0] 100%</span>
                  </span>
                </div>
                <div className="text-[11px] text-zinc-400 font-mono flex items-center gap-2 mt-0.5">
                  <span className="text-amber-300 font-bold">R_max: {activeBlock.envelopeRadius} Å</span>
                  <span>•</span>
                  <span>Średnica: <strong className="text-zinc-200">{activeBlock.collisionDiameter} Å</strong></span>
                  <span>•</span>
                  <span>VdW: <strong className="text-zinc-200">{activeBlock.physico.vdwVolume} Å³</strong></span>
                </div>
              </div>
            </div>

            {/* Quick Actions (Rotate, Opacity, Atoms, Zoom) */}
            <div className="pointer-events-auto flex items-center gap-2 bg-black/80 backdrop-blur-md border border-zinc-800 p-1.5 rounded-lg shadow-xl font-mono text-xs">
              <button
                onClick={() => setAutoRotate(!autoRotate)}
                className={`px-2 py-1 rounded flex items-center gap-1 cursor-pointer transition ${
                  autoRotate ? 'bg-cyan-600 text-white' : 'bg-zinc-900 text-zinc-400 hover:text-white'
                }`}
                title="Automatyczny obrót kuli 360°"
              >
                <RotateCcw className={`w-3.5 h-3.5 ${autoRotate ? 'animate-spin' : ''}`} />
                <span className="text-[10px]">Obrót</span>
              </button>

              <button
                onClick={() => setShowAtoms(!showAtoms)}
                className={`px-2 py-1 rounded flex items-center gap-1 cursor-pointer transition ${
                  showAtoms ? 'bg-indigo-600 text-white' : 'bg-zinc-900 text-zinc-400 hover:text-white'
                }`}
                title="Pokaż/Ukryj model atomowy w środku kuli"
              >
                <Atom className="w-3.5 h-3.5" />
                <span className="text-[10px]">Atomy</span>
              </button>

              <button
                onClick={() => setSphereOpacity(prev => prev > 0.5 ? 0.25 : 0.65)}
                className="px-2 py-1 rounded bg-zinc-900 hover:bg-zinc-800 text-zinc-300 flex items-center gap-1 cursor-pointer"
                title="Zmień przezroczystość kuli 248 trójkątów"
              >
                <Layers className="w-3.5 h-3.5" />
                <span className="text-[10px]">Klatka {Math.round(sphereOpacity * 100)}%</span>
              </button>

              <div className="flex items-center gap-1 border-l border-zinc-800 pl-2">
                <button
                  onClick={() => setZoomLevel(prev => Math.max(0.4, prev - 0.2))}
                  className="p-1 rounded bg-zinc-900 hover:bg-zinc-800 text-zinc-300"
                  title="Oddal"
                >
                  <ZoomOut className="w-3.5 h-3.5" />
                </button>
                <span className="text-[10px] text-zinc-400 min-w-[28px] text-center font-mono">
                  {Math.round(zoomLevel * 100)}%
                </span>
                <button
                  onClick={() => setZoomLevel(prev => Math.min(2.5, prev + 0.2))}
                  className="p-1 rounded bg-zinc-900 hover:bg-zinc-800 text-zinc-300"
                  title="Przybliż"
                >
                  <ZoomIn className="w-3.5 h-3.5" />
                </button>
              </div>

              {/* Reset view */}
              <button
                onClick={() => {
                  rotAngle.current = { x: 0.2, y: 0.4 };
                  setZoomLevel(1.0);
                }}
                className="p-1 rounded bg-zinc-900 hover:bg-zinc-800 text-zinc-400 hover:text-white"
                title="Zresetuj kamerę"
              >
                <Compass className="w-3.5 h-3.5" />
              </button>
            </div>
          </div>

          {/* 3D Three.js Mounting Container */}
          <div 
            ref={mountRef}
            className="flex-1 w-full h-full min-h-0 cursor-grab active:cursor-grabbing select-none"
            onMouseDown={handleMouseDown}
            onMouseMove={handleMouseMove}
            onMouseUp={handleMouseUp}
            onMouseLeave={handleMouseUp}
            onWheel={handleWheel}
          />

          {/* Bottom Port Legend & Cardinal Vectors Indicator */}
          <div className="absolute bottom-3 left-3 right-3 z-10 flex flex-wrap items-center justify-between gap-2 pointer-events-none">
            <div className="pointer-events-auto bg-black/85 backdrop-blur-md border border-zinc-800 rounded-lg px-3 py-1.5 flex items-center gap-3 text-[11px] font-mono shadow-xl">
              <div className="flex items-center gap-1.5">
                <span className="w-2.5 h-2.5 rounded-full bg-cyan-400 shadow-sm shadow-cyan-400/50" />
                <span className="text-zinc-300 font-bold">Łatka #1 (+Y):</span>
                <span className="text-cyan-300">Wejście N wiązania peptydowego</span>
              </div>
              <span className="text-zinc-600">|</span>
              <div className="flex items-center gap-1.5">
                <span className="w-2.5 h-2.5 rounded-full bg-rose-500 shadow-sm shadow-rose-500/50" />
                <span className="text-zinc-300 font-bold">Łatka #248 (-Y):</span>
                <span className="text-rose-300">Wyjście C wiązania peptydowego</span>
              </div>
              <span className="text-zinc-600">|</span>
              <div className="flex items-center gap-1.5 text-zinc-400">
                <span className="w-2 h-2 rounded-full bg-yellow-400" />
                <span className="text-yellow-400 font-bold">Środek: [0,0,0]</span>
              </div>
            </div>

            {/* Action: Send Block to Constructor */}
            <div className="pointer-events-auto">
              <button
                onClick={() => {
                  if (onSelectAminoAcidForConstructor) {
                    onSelectAminoAcidForConstructor(selectedCode);
                  }
                  if (onSwitchToConstructor) {
                    onSwitchToConstructor();
                  }
                }}
                className="px-4 py-2 bg-gradient-to-r from-emerald-600 to-teal-600 hover:from-emerald-500 hover:to-teal-500 text-white font-mono text-xs font-bold rounded-lg shadow-lg shadow-emerald-500/20 flex items-center gap-2 cursor-pointer transition"
                title="Wstaw ten klocek sferyczny do konstruktora zwijania na Karcie B"
              >
                <span>Użyj w Konstruktorze (Karta B)</span>
                <ArrowRight className="w-4 h-4 text-emerald-200" />
              </button>
            </div>
          </div>
        </div>

        {/* RIGHT COLUMN: DEDICATED INSPECTOR (FACET E8 OR COMPLETE 22-PARAMETER TABLE) */}
        <div className="w-full lg:w-[440px] shrink-0 border-t lg:border-t-0 lg:border-l border-zinc-800 bg-[#06070f] flex flex-col overflow-y-auto font-mono text-xs">
          
          {/* Header & Mode Switcher */}
          <div className="p-2.5 border-b border-zinc-800/80 bg-[#090b17] flex items-center justify-between gap-2">
            <div className="inline-flex rounded-lg border border-zinc-800 bg-black p-0.5 text-[11px]">
              <button
                onClick={() => setInspectorTab('FACET')}
                className={`px-3 py-1 rounded-md font-bold transition cursor-pointer flex items-center gap-1.5 ${
                  inspectorTab === 'FACET'
                    ? 'bg-cyan-600 text-white shadow-sm'
                    : 'text-zinc-400 hover:text-white'
                }`}
              >
                <Compass className="w-3.5 h-3.5" />
                <span>Łatka #{activeFacet.id}</span>
              </button>
              <button
                onClick={() => setInspectorTab('TABLE_22')}
                className={`px-3 py-1 rounded-md font-bold transition cursor-pointer flex items-center gap-1.5 ${
                  inspectorTab === 'TABLE_22'
                    ? 'bg-cyan-600 text-white shadow-sm'
                    : 'text-zinc-400 hover:text-white'
                }`}
              >
                <Table className="w-3.5 h-3.5" />
                <span>Tabela 22 Parametrów</span>
              </button>
            </div>

            <span className="text-[10.5px] text-zinc-400">
              {activeBlock.physico.code3} • {activeBlock.envelopeRadius} Å
            </span>
          </div>

          {/* VIEW MODE 1: 248 FACET DETAIL */}
          {inspectorTab === 'FACET' && (
            <div className="p-3.5 space-y-3.5">
              
              {/* Quick Facet Selector (Number buttons / jump) */}
              <div className="p-3 rounded-lg bg-black/60 border border-zinc-800 space-y-1.5">
                <div className="flex items-center justify-between text-[10.5px] text-zinc-400">
                  <span>Wybierz łatkę na kuli (1..248):</span>
                  <span className="text-zinc-500">Bieguny: #1 (N), #248 (C)</span>
                </div>
                <div className="flex items-center gap-1.5 flex-wrap">
                  {[1, 50, 100, 150, 200, 248].map(pid => (
                    <button
                      key={pid}
                      onClick={() => {
                        setSelectedFacetId(pid);
                        playPatchTone(pid);
                      }}
                      className={`px-2 py-0.5 rounded text-[10.5px] font-bold cursor-pointer transition ${
                        selectedFacetId === pid
                          ? 'bg-cyan-500 text-black shadow-md'
                          : pid === 1 ? 'bg-cyan-950 text-cyan-300 border border-cyan-800'
                          : pid === 248 ? 'bg-rose-950 text-rose-300 border border-rose-800'
                          : 'bg-zinc-900 hover:bg-zinc-800 text-zinc-300'
                      }`}
                    >
                      #{pid} {pid === 1 ? '(N)' : pid === 248 ? '(C)' : ''}
                    </button>
                  ))}
                  <input
                    type="number"
                    min="1"
                    max="248"
                    value={selectedFacetId}
                    onChange={(e) => {
                      const val = Math.max(1, Math.min(248, parseInt(e.target.value) || 1));
                      setSelectedFacetId(val);
                    }}
                    className="w-16 px-1.5 py-0.5 bg-black border border-zinc-700 rounded text-center text-zinc-200 text-xs focus:outline-none focus:border-cyan-500"
                  />
                </div>
              </div>

              {/* Vector Coordinates & Distance from CoM */}
              <div className="p-3 rounded-lg bg-black/60 border border-zinc-800 space-y-2">
                <div className="flex items-center justify-between text-[11px]">
                  <span className="text-zinc-400 font-bold uppercase">Wektor Normalny (Sfera E8):</span>
                  <span className="text-[10px] text-amber-300 font-bold">
                    d = {activeBlock.envelopeRadius} Å
                  </span>
                </div>
                <div className="grid grid-cols-3 gap-1.5 text-center text-[11px]">
                  <div className="p-1 rounded bg-zinc-900 border border-zinc-800">
                    <span className="text-red-400 font-bold">X:</span> {activeFacet.normal[0].toFixed(3)}
                  </div>
                  <div className="p-1 rounded bg-zinc-900 border border-zinc-800">
                    <span className="text-emerald-400 font-bold">Y:</span> {activeFacet.normal[1].toFixed(3)}
                  </div>
                  <div className="p-1 rounded bg-zinc-900 border border-zinc-800">
                    <span className="text-blue-400 font-bold">Z:</span> {activeFacet.normal[2].toFixed(3)}
                  </div>
                </div>
                {activeFacet.cardinalLabel && (
                  <div className="text-[10px] text-cyan-300 bg-cyan-950/60 p-1.5 rounded border border-cyan-800/60">
                    {activeFacet.cardinalLabel}
                  </div>
                )}
              </div>

              {/* Biochemical Properties on this patch */}
              <div className="p-3 rounded-lg bg-black/60 border border-zinc-800 space-y-2">
                <span className="text-zinc-400 font-bold uppercase text-[11px] block">
                  Właściwości Biochemiczne w Kierunku Łatki #{activeFacet.id}:
                </span>

                <div className="space-y-1.5 text-[11px]">
                  <div className="flex items-center justify-between p-1.5 rounded bg-zinc-900/80">
                    <span className="text-zinc-400">Sektor cząsteczki:</span>
                    <span className={`font-bold ${activeFacet.isSidechainR ? 'text-amber-300' : 'text-indigo-300'}`}>
                      {activeFacet.isSidechainR ? 'Łańcuch boczny (Grupa R)' : 'Szkielet peptydowy (Backbone)'}
                    </span>
                  </div>

                  <div className="flex items-center justify-between p-1.5 rounded bg-zinc-900/80">
                    <span className="text-zinc-400">Lokalna hydropatia (Kyte):</span>
                    <span className={`font-bold ${activeFacet.hydrophobicity > 0 ? 'text-amber-400' : 'text-blue-400'}`}>
                      {activeFacet.hydrophobicity > 0 ? `+${activeFacet.hydrophobicity}` : activeFacet.hydrophobicity}
                    </span>
                  </div>

                  <div className="flex items-center justify-between p-1.5 rounded bg-zinc-900/80">
                    <span className="text-zinc-400">Lokalny ładunek (pH 7.4):</span>
                    <span className={`font-bold ${activeFacet.chargeAtPH74 > 0 ? 'text-blue-400' : activeFacet.chargeAtPH74 < 0 ? 'text-red-400' : 'text-zinc-400'}`}>
                      {activeFacet.chargeAtPH74 > 0 ? `+${activeFacet.chargeAtPH74}` : activeFacet.chargeAtPH74}
                    </span>
                  </div>

                  <div className="flex items-center justify-between p-1.5 rounded bg-zinc-900/80">
                    <span className="text-zinc-400">Wiązania wodorowe (H-Bonds):</span>
                    <span className="font-bold text-zinc-200">
                      {activeFacet.hBondRole === 'DONOR' ? 'DONOR (Grupa NH)' :
                       activeFacet.hBondRole === 'ACCEPTOR' ? 'AKCEPTOR (Grupa C=O)' :
                       activeFacet.hBondRole === 'BOTH' ? 'DONOR / AKCEPTOR' : 'Obojętny'}
                    </span>
                  </div>

                  <div className="flex items-center justify-between p-1.5 rounded bg-zinc-900/80">
                    <span className="text-zinc-400">Promień zewn. sfery (R_max):</span>
                    <span className="font-bold text-emerald-400">
                      {activeBlock.envelopeRadius} Å
                    </span>
                  </div>

                  {activeFacet.nearestAtom && (
                    <div className="flex items-center justify-between p-1.5 rounded bg-zinc-900/80">
                      <span className="text-zinc-400">Najbliższy rzutowany atom:</span>
                      <span className="font-bold text-cyan-300">
                        {activeFacet.nearestAtom.id} ({activeFacet.nearestAtom.element}) • d={activeFacet.nearestAtom.distAngstrom} Å
                      </span>
                    </div>
                  )}
                </div>
              </div>

            </div>
          )}

          {/* VIEW MODE 2: COMPLETE 22-PARAMETERS TABLE */}
          {inspectorTab === 'TABLE_22' && (
            <div className="p-3.5 space-y-3">
              <div className="flex items-center justify-between pb-2 border-b border-zinc-800">
                <span className="font-bold text-zinc-200 uppercase text-[11px]">
                  Baza Fizykochemiczna SNAPLOCK: {activeBlock.physico.namePl} ({activeBlock.physico.code3})
                </span>
                <span className="text-[10px] px-1.5 py-0.2 rounded bg-zinc-800 text-zinc-400">
                  Nr #{activeBlock.physico.nr} / 21
                </span>
              </div>

              {/* The 22 Physical Parameters Table */}
              <div className="space-y-1 text-[11px]">
                <div className="flex justify-between p-1.5 rounded bg-zinc-900/70 border border-zinc-800/80">
                  <span className="text-zinc-400">1. Nr porządkowy:</span>
                  <span className="font-bold text-white">#{activeBlock.physico.nr}</span>
                </div>

                <div className="flex justify-between p-1.5 rounded bg-zinc-900/70 border border-zinc-800/80">
                  <span className="text-zinc-400">2. Kod (3-lit / 1-lit):</span>
                  <span className="font-bold text-cyan-300">{activeBlock.physico.code3} / {activeBlock.physico.code1}</span>
                </div>

                <div className="flex justify-between p-1.5 rounded bg-zinc-900/70 border border-zinc-800/80">
                  <span className="text-zinc-400">3. Nazwa polska:</span>
                  <span className="font-bold text-white">{activeBlock.physico.namePl}</span>
                </div>

                <div className="flex justify-between p-1.5 rounded bg-zinc-900/70 border border-zinc-800/80">
                  <span className="text-zinc-400">4. Wzór łańcucha bocznego:</span>
                  <span className="font-bold text-amber-300">{activeBlock.physico.sidechainFormula}</span>
                </div>

                <div className="flex justify-between p-1.5 rounded bg-zinc-900/70 border border-zinc-800/80">
                  <span className="text-zinc-400">5. Masa cząsteczkowa:</span>
                  <span className="font-bold text-white">{activeBlock.physico.mass} g/mol</span>
                </div>

                <div className="flex justify-between p-1.5 rounded bg-zinc-900/70 border border-zinc-800/80">
                  <span className="text-zinc-400">6. Objętość Van der Waalsa:</span>
                  <span className="font-bold text-emerald-400">{activeBlock.physico.vdwVolume} Å³</span>
                </div>

                <div className="flex justify-between p-1.5 rounded bg-zinc-900/70 border border-zinc-800/80">
                  <span className="text-zinc-400">7. Dostępność solwatacyjna (SASA):</span>
                  <span className="font-bold text-white">{activeBlock.physico.sasa} Å²</span>
                </div>

                <div className="flex justify-between p-1.5 rounded bg-zinc-900/70 border border-zinc-800/80">
                  <span className="text-zinc-400">8. Polarność:</span>
                  <span className="font-bold text-white">{activeBlock.physico.polarity}</span>
                </div>

                <div className="flex justify-between p-1.5 rounded bg-zinc-900/70 border border-zinc-800/80">
                  <span className="text-zinc-400">9. Hydropatia (Kyte):</span>
                  <span className={`font-bold ${activeBlock.physico.hydropathyKyte > 0 ? 'text-amber-400' : 'text-blue-400'}`}>
                    {activeBlock.physico.hydropathyKyte > 0 ? `+${activeBlock.physico.hydropathyKyte}` : activeBlock.physico.hydropathyKyte}
                  </span>
                </div>

                <div className="flex justify-between p-1.5 rounded bg-zinc-900/70 border border-zinc-800/80">
                  <span className="text-zinc-400">10. Hydrofilowość (Hopp):</span>
                  <span className="font-bold text-white">{activeBlock.physico.hydrophilicityHopp > 0 ? `+${activeBlock.physico.hydrophilicityHopp}` : activeBlock.physico.hydrophilicityHopp}</span>
                </div>

                <div className="flex justify-between p-1.5 rounded bg-zinc-900/70 border border-zinc-800/80">
                  <span className="text-zinc-400">11. Rozpuszczalność:</span>
                  <span className="font-bold text-white">{activeBlock.physico.solubility}</span>
                </div>

                <div className="flex justify-between p-1.5 rounded bg-zinc-900/70 border border-zinc-800/80">
                  <span className="text-zinc-400">12. Ładunek w pH 7.4:</span>
                  <span className={`font-bold ${activeBlock.physico.chargePH74 > 0 ? 'text-blue-400' : activeBlock.physico.chargePH74 < 0 ? 'text-red-400' : 'text-zinc-300'}`}>
                    {activeBlock.physico.chargePH74 > 0 ? `+${activeBlock.physico.chargePH74}` : activeBlock.physico.chargePH74}
                  </span>
                </div>

                <div className="flex justify-between p-1.5 rounded bg-zinc-900/70 border border-zinc-800/80">
                  <span className="text-zinc-400">13. Jonizacja boczna:</span>
                  <span className="font-bold text-white">{activeBlock.physico.sidechainIonization ? 'Tak' : 'Nie'}</span>
                </div>

                <div className="flex justify-between p-1.5 rounded bg-zinc-900/70 border border-zinc-800/80">
                  <span className="text-zinc-400">14. pK1 (alfa-COOH):</span>
                  <span className="font-bold text-white">{activeBlock.physico.pK1_COOH}</span>
                </div>

                <div className="flex justify-between p-1.5 rounded bg-zinc-900/70 border border-zinc-800/80">
                  <span className="text-zinc-400">15. pK2 (alfa-NH3):</span>
                  <span className="font-bold text-white">{activeBlock.physico.pK2_NH3}</span>
                </div>

                <div className="flex justify-between p-1.5 rounded bg-zinc-900/70 border border-zinc-800/80">
                  <span className="text-zinc-400">16. pKr (łańcuch boczny):</span>
                  <span className="font-bold text-white">{activeBlock.physico.pKr_side !== null ? activeBlock.physico.pKr_side : 'B.d.'}</span>
                </div>

                <div className="flex justify-between p-1.5 rounded bg-zinc-900/70 border border-zinc-800/80">
                  <span className="text-zinc-400">17. Punkt izoelektryczny (pI):</span>
                  <span className="font-bold text-white">{activeBlock.physico.pI}</span>
                </div>

                <div className="flex justify-between p-1.5 rounded bg-zinc-900/70 border border-zinc-800/80">
                  <span className="text-zinc-400">18. Mostki kowalencyjne:</span>
                  <span className={`font-bold ${activeBlock.physico.covalentBridges.includes('TAK') ? 'text-amber-400' : 'text-zinc-400'}`}>
                    {activeBlock.physico.covalentBridges}
                  </span>
                </div>

                <div className="flex justify-between p-1.5 rounded bg-zinc-900/70 border border-zinc-800/80">
                  <span className="text-zinc-400">19. Chiralność:</span>
                  <span className="font-bold text-white">{activeBlock.physico.chirality}</span>
                </div>

                <div className="flex justify-between p-1.5 rounded bg-zinc-900/70 border border-zinc-800/80">
                  <span className="text-zinc-400">20. Preferencja Ramachandrana:</span>
                  <span className="font-bold text-white">{activeBlock.physico.ramachandranPref}</span>
                </div>

                <div className="flex justify-between p-1.5 rounded bg-zinc-900/70 border border-zinc-800/80">
                  <span className="text-zinc-400">21. Elastyczność łańcucha:</span>
                  <span className="font-bold text-white">{activeBlock.physico.chainFlexibility}</span>
                </div>

                <div className="flex justify-between p-1.5 rounded bg-zinc-900/70 border border-zinc-800/80">
                  <span className="text-zinc-400">22. Preferencje Rotamerowe:</span>
                  <span className="font-bold text-white">{activeBlock.physico.rotamerPreferences}</span>
                </div>

                {/* Spatial Envelope additions */}
                <div className="flex justify-between p-1.5 rounded bg-cyan-950/40 border border-cyan-800/60 font-bold">
                  <span className="text-cyan-300">23. Promień zewn. R_max (CoM do VdW):</span>
                  <span className="text-amber-300">{activeBlock.envelopeRadius} Å</span>
                </div>

                <div className="flex justify-between p-1.5 rounded bg-cyan-950/40 border border-cyan-800/60 font-bold">
                  <span className="text-cyan-300">24. Średnica ochronna zderzeń:</span>
                  <span className="text-amber-300">{activeBlock.collisionDiameter} Å</span>
                </div>
              </div>
            </div>
          )}

        </div>
      </div>

      {/* 3. BOTTOM RIBBON: THE 21 PIŁECZEK / KULEK MAGAZYNU */}
      <div className="bg-[#05060d] border-t border-zinc-800/90 px-3 py-2 shrink-0 flex flex-col gap-1.5">
        <div className="flex items-center justify-between text-[11px] font-mono px-1">
          <div className="flex items-center gap-2">
            <span className="text-zinc-400 font-bold uppercase tracking-wider">
              Wybierz jedną z 21 Kul E8:
            </span>
            <span className="text-[10px] text-zinc-500 hidden sm:inline">
              (każda kula ma promień dopasowany do swojej powłoki elektronowej VdW, zero wystawania atomów)
            </span>
          </div>
          <div className="flex items-center gap-1.5">
            <Search className="w-3.5 h-3.5 text-zinc-500" />
            <input
              type="text"
              placeholder="Szukaj (np. CYS, ALA)..."
              value={searchQuery}
              onChange={(e) => setSearchQuery(e.target.value)}
              className="bg-black/60 border border-zinc-800 rounded px-2 py-0.5 text-[10.5px] w-40 text-white focus:outline-none focus:border-cyan-500"
            />
          </div>
        </div>

        {/* 21 Interactive Sphere Buttons */}
        <div className="flex items-center gap-2 overflow-x-auto pb-1 scrollbar-thin scrollbar-thumb-zinc-800">
          {filteredList.map((aa) => {
            const isSelected = selectedCode === aa.physico.code3;
            return (
              <button
                key={aa.physico.code3}
                onClick={() => {
                  setSelectedCode(aa.physico.code3);
                  playPatchTone(1);
                }}
                className={`flex items-center gap-2 px-3 py-1.5 rounded-xl border text-left shrink-0 transition-all cursor-pointer ${
                  isSelected
                    ? 'bg-[#101935] border-cyan-500 ring-2 ring-cyan-500/30 shadow-lg shadow-cyan-500/10'
                    : 'bg-[#090b14] border-zinc-800/80 hover:bg-[#121424] hover:border-zinc-700'
                }`}
              >
                {/* 3D Ball Icon */}
                <div 
                  className="w-8 h-8 rounded-full flex items-center justify-center font-mono font-bold text-xs shadow-md shrink-0 border relative overflow-hidden"
                  style={{
                    backgroundColor: `${aa.color}20`,
                    borderColor: isSelected ? '#06b6d4' : `${aa.color}80`,
                    color: isSelected ? '#fff' : aa.color
                  }}
                >
                  <div className="absolute top-1 left-1 w-2.5 h-2.5 rounded-full bg-white/40 blur-[0.5px]" />
                  <span>{aa.physico.code3}</span>
                </div>

                {/* Details */}
                <div className="min-w-[65px]">
                  <div className="flex items-center gap-1">
                    <span className="font-mono text-xs font-bold text-zinc-100">{aa.physico.namePl}</span>
                  </div>
                  <div className="text-[9.5px] text-zinc-400 font-mono flex items-center gap-1">
                    <span>R: <strong className="text-amber-300">{aa.envelopeRadius}Å</strong></span>
                    {aa.physico.covalentBridges.includes('TAK') && (
                      <span className="text-[8px] px-1 py-0 rounded bg-amber-950 text-amber-300 border border-amber-800/60 font-bold">
                        S-S
                      </span>
                    )}
                  </div>
                </div>
              </button>
            );
          })}
        </div>
      </div>

    </div>
  );
};
