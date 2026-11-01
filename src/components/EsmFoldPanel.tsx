// Karta I: prawdziwa predykcja struktury (ESMFold, Meta) + ocena walidatorem SNAPLOCK
import React, { useEffect, useMemo, useRef, useState } from 'react';
import * as THREE from 'three';
import { LineChart, Line, XAxis, YAxis, Tooltip, ResponsiveContainer, ReferenceLine } from 'recharts';
import { validateAgainstGroundTruth } from '../ground_truth_validation';

const T: Record<string, string> = {
  ALA: 'A', ARG: 'R', ASN: 'N', ASP: 'D', CYS: 'C', GLN: 'Q', GLU: 'E', GLY: 'G', HIS: 'H', ILE: 'I',
  LEU: 'L', LYS: 'K', MET: 'M', PHE: 'F', PRO: 'P', SER: 'S', THR: 'T', TRP: 'W', TYR: 'Y', VAL: 'V'
};

function toOneLetter(input: string): string {
  const s = input.trim();
  if (!s) return '';
  const parts = s.split(/[\s,;-]+/).filter(Boolean);
  if (parts.length > 1 && parts.every((p) => p.length === 3)) {
    return parts.map((p) => T[p.toUpperCase()] || 'X').join('');
  }
  return s.replace(/[^A-Za-z]/g, '').toUpperCase();
}

interface CaAtom { resSeq: number; x: number; y: number; z: number; plddt: number; }

function parseCa(pdb: string): CaAtom[] {
  const out: CaAtom[] = [];
  for (const l of pdb.split('\n')) {
    if (l.startsWith('ATOM') && l.substring(12, 16).trim() === 'CA') {
      const b = parseFloat(l.substring(60, 66));
      out.push({
        resSeq: parseInt(l.substring(22, 26), 10),
        x: parseFloat(l.substring(30, 38)),
        y: parseFloat(l.substring(38, 46)),
        z: parseFloat(l.substring(46, 54)),
        plddt: Number.isNaN(b) ? 0 : b <= 1 ? b * 100 : b
      });
    }
  }
  return out;
}

function plddtColor(p: number): number {
  if (p >= 90) return 0x0053d6;
  if (p >= 70) return 0x65cbf3;
  if (p >= 50) return 0xffdb13;
  return 0xff7d45;
}

interface Props {
  activeSequence: string;
  referencePdbText?: string | null;
  referenceChainId?: string;
  referenceLabel?: string;
}

export const EsmFoldPanel: React.FC<Props> = ({ activeSequence, referencePdbText, referenceChainId, referenceLabel }) => {
  const [seq, setSeq] = useState<string>(() => toOneLetter(activeSequence));
  const [pdb, setPdb] = useState<string>('');
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState('');
  const mountRef = useRef<HTMLDivElement>(null);

  useEffect(() => { setSeq(toOneLetter(activeSequence)); }, [activeSequence]);

  const ca = useMemo(() => (pdb ? parseCa(pdb) : []), [pdb]);
  const meanPlddt = ca.length ? ca.reduce((a, c) => a + c.plddt, 0) / ca.length : null;

  const validation = useMemo(() => {
    if (!pdb || !referencePdbText) return null;
    try {
      return validateAgainstGroundTruth(pdb, referencePdbText, {
        modelChainId: 'A',
        referenceChainId: referenceChainId || 'A',
        referenceId: referenceLabel,
        allowSubsequenceAlignment: false
      });
    } catch {
      return null;
    }
  }, [pdb, referencePdbText, referenceChainId, referenceLabel]);

  const sequenceMatchesReference = validation?.sequenceMapping?.sequenceIdentityPct === 100;

  async function predict() {
    setError('');
    const clean = toOneLetter(seq);
    if (clean.length < 10 || clean.length > 400) { setError('Sekwencja musi mieć 10–400 aminokwasów.'); return; }
    if (/[^ACDEFGHIKLMNPQRSTVWY]/.test(clean)) { setError('Tylko 20 standardowych aminokwasów (kody 1-literowe).'); return; }
    setBusy(true);
    try {
      const resp = await fetch('/api/esmfold', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ sequence: clean })
      });
      const data = await resp.json();
      if (!resp.ok || !data.pdb) throw new Error(data.error || `HTTP ${resp.status}`);
      setPdb(data.pdb);
    } catch (e: any) {
      setError(`ESMFold: ${e.message || e}`);
    } finally {
      setBusy(false);
    }
  }

  function downloadPdb() {
    const blob = new Blob([pdb], { type: 'chemical/x-pdb' });
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = `esmfold_${toOneLetter(seq).slice(0, 12)}.pdb`;
    a.click();
  }

  useEffect(() => {
    const mount = mountRef.current;
    if (!mount || ca.length < 2) return;
    const w = mount.clientWidth || 600;
    const h = mount.clientHeight || 420;
    const renderer = new THREE.WebGLRenderer({ antialias: true });
    renderer.setSize(w, h);
    renderer.setClearColor(0x05060c);
    mount.appendChild(renderer.domElement);
    const scene = new THREE.Scene();
    const camera = new THREE.PerspectiveCamera(45, w / h, 0.1, 2000);
    scene.add(new THREE.AmbientLight(0xffffff, 0.6));
    const dir = new THREE.DirectionalLight(0xffffff, 0.9);
    dir.position.set(1, 1, 2);
    scene.add(dir);

    const cx = ca.reduce((a, c) => a + c.x, 0) / ca.length;
    const cy = ca.reduce((a, c) => a + c.y, 0) / ca.length;
    const cz = ca.reduce((a, c) => a + c.z, 0) / ca.length;
    const group = new THREE.Group();
    let radius = 0;
    const pts = ca.map((c) => {
      const v = new THREE.Vector3(c.x - cx, c.y - cy, c.z - cz);
      radius = Math.max(radius, v.length());
      return v;
    });
    pts.forEach((p, i) => {
      const s = new THREE.Mesh(
        new THREE.SphereGeometry(0.7, 12, 12),
        new THREE.MeshPhongMaterial({ color: plddtColor(ca[i].plddt) })
      );
      s.position.copy(p);
      group.add(s);
      if (i > 0) {
        const a = pts[i - 1];
        const len = a.distanceTo(p);
        const cyl = new THREE.Mesh(
          new THREE.CylinderGeometry(0.3, 0.3, len, 8),
          new THREE.MeshPhongMaterial({ color: plddtColor(ca[i].plddt) })
        );
        cyl.position.copy(a.clone().add(p).multiplyScalar(0.5));
        cyl.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), p.clone().sub(a).normalize());
        group.add(cyl);
      }
    });
    scene.add(group);
    camera.position.set(0, 0, radius * 2.6 + 10);

    let dragging = false, lx = 0, ly = 0, frame = 0;
    const down = (e: PointerEvent) => { dragging = true; lx = e.clientX; ly = e.clientY; };
    const up = () => { dragging = false; };
    const move = (e: PointerEvent) => {
      if (!dragging) return;
      group.rotation.y += (e.clientX - lx) * 0.01;
      group.rotation.x += (e.clientY - ly) * 0.01;
      lx = e.clientX; ly = e.clientY;
    };
    const wheel = (e: WheelEvent) => { e.preventDefault(); camera.position.z = Math.max(10, camera.position.z + e.deltaY * 0.05); };
    renderer.domElement.addEventListener('pointerdown', down);
    window.addEventListener('pointerup', up);
    window.addEventListener('pointermove', move);
    renderer.domElement.addEventListener('wheel', wheel, { passive: false });
    const loop = () => {
      if (!dragging) group.rotation.y += 0.003;
      renderer.render(scene, camera);
      frame = requestAnimationFrame(loop);
    };
    loop();
    return () => {
      cancelAnimationFrame(frame);
      window.removeEventListener('pointerup', up);
      window.removeEventListener('pointermove', move);
      renderer.dispose();
      mount.removeChild(renderer.domElement);
    };
  }, [ca]);

  const fmt = (v: number | null | undefined, d = 2) => (v === null || v === undefined ? '—' : v.toFixed(d));

  return (
    <div className="flex-1 min-h-0 overflow-auto p-4 font-mono text-sm text-zinc-200">
      <h2 className="text-lg font-bold text-white mb-1">I: Predykcja ESMFold + walidacja SNAPLOCK</h2>
      <p className="text-zinc-400 text-xs mb-3">
        Prawdziwa predykcja struktury z sekwencji (ESMFold, Meta, publiczne API, do 400 aa). Kolor = pewność pLDDT:
        <span className="text-[#0053d6]"> ≥90</span>, <span className="text-[#65cbf3]">70–90</span>,
        <span className="text-[#ffdb13]"> 50–70</span>, <span className="text-[#ff7d45]">&lt;50</span>.
      </p>

      <textarea
        value={seq}
        onChange={(e) => setSeq(e.target.value)}
        rows={3}
        className="w-full bg-[#0b0d18] border border-zinc-700 rounded p-2 text-xs"
        placeholder="Sekwencja (kody 1-literowe lub 3-literowe oddzielone spacjami)"
      />
      <div className="flex gap-2 mt-2 items-center">
        <button onClick={predict} disabled={busy}
          className="px-4 py-1.5 rounded bg-emerald-600 hover:bg-emerald-500 disabled:opacity-50 text-white font-bold">
          {busy ? 'Liczę… (do ~1 min)' : 'Przewiduj (ESMFold)'}
        </button>
        {pdb && (
          <button onClick={downloadPdb} className="px-3 py-1.5 rounded bg-zinc-700 hover:bg-zinc-600 text-white">
            Pobierz PDB
          </button>
        )}
        <span className="text-xs text-zinc-500">{toOneLetter(seq).length} aa</span>
      </div>
      {error && <div className="mt-2 text-red-400 text-xs">{error}</div>}

      {pdb && (
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-3 mt-4">
          <div ref={mountRef} className="lg:col-span-2 h-[420px] rounded border border-zinc-800 overflow-hidden" />
          <div className="space-y-3">
            <div className="rounded border border-zinc-800 p-3">
              <div className="text-zinc-400 text-xs">Średnia pewność (pLDDT)</div>
              <div className="text-2xl font-bold text-white">{fmt(meanPlddt, 1)}</div>
              <div className="text-xs text-zinc-500">{ca.length} reszt</div>
            </div>
            <div className="rounded border border-zinc-800 p-3">
              <div className="text-zinc-400 text-xs mb-1">Porównanie z eksperymentem {referenceLabel ? `(${referenceLabel})` : ''}</div>
              {!referencePdbText && <div className="text-xs text-zinc-500">Brak struktury referencyjnej dla tej sekwencji.</div>}
              {referencePdbText && !validation && <div className="text-xs text-zinc-500">Nie udało się porównać.</div>}
              {validation && !sequenceMatchesReference && (
                <div className="text-xs text-amber-400">Sekwencja różni się od referencji — porównanie niewiarygodne.</div>
              )}
              {validation && sequenceMatchesReference && (
                <div className="text-xs space-y-1">
                  <div>RMSD Cα: <b className="text-white">{fmt(validation.caRmsdAngstrom)} Å</b></div>
                  <div>TM-score: <b className="text-white">{fmt(validation.tmScore, 3)}</b> <span className="text-zinc-500">(&gt;0,5 = ten sam fałd)</span></div>
                  <div>Sparowane reszty: {validation.pairedResidueCount}</div>
                </div>
              )}
            </div>
          </div>
          <div className="lg:col-span-3 h-48 rounded border border-zinc-800 p-2">
            <ResponsiveContainer width="100%" height="100%">
              <LineChart data={ca.map((c) => ({ res: c.resSeq, plddt: Number(c.plddt.toFixed(1)) }))}>
                <XAxis dataKey="res" stroke="#71717a" fontSize={10} />
                <YAxis domain={[0, 100]} stroke="#71717a" fontSize={10} />
                <Tooltip contentStyle={{ background: '#0b0d18', border: '1px solid #3f3f46' }} />
                <ReferenceLine y={70} stroke="#65cbf3" strokeDasharray="3 3" />
                <Line type="monotone" dataKey="plddt" stroke="#10b981" dot={false} strokeWidth={2} />
              </LineChart>
            </ResponsiveContainer>
          </div>
        </div>
      )}
    </div>
  );
};
