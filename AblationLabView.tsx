import React, { useEffect, useRef, useState } from 'react';
import { runMultiAgentFoldingToConvergence } from './src/snaplockMultiAgentEngine';
import { GLOBAL_PROTEIN_TARGETS } from './src/snaplockGlobalTargets';

// =====================================================================================
// KARTA J: Laboratorium wariantów A/B + Dziennik wniosków
// A = ścieżka bazowa. Każdy moduł (B, C, D, E) dokładany osobno i w zestawach.
// Każdy wariant liczony z kilku startów (średnia ± odchylenie), bo zwijanie jest chaotyczne.
// Wszystko na ślepo: silnik NIE widzi struktury docelowej (enableTargetGuidedRefinement: false).
// =====================================================================================

type ModuleKey = 'B' | 'C' | 'D' | 'E' | 'F';
const MODULES: { key: ModuleKey; name: string; desc: string }[] = [
  { key: 'B', name: 'Kontakty MJ', desc: 'Energie kontaktu Miyazawa-Jernigan 1996 (210 par aminokwasów)' },
  { key: 'C', name: 'Wiązania H', desc: 'Wiązania wodorowe szkieletu O…N (2,9 Å)' },
  { key: 'D', name: 'GOR', desc: 'Struktura drugorzędowa wyuczona z 691 białek PDB (Q3 62,6%)' },
  { key: 'E', name: 'Monte Carlo', desc: 'Wyżarzanie Monte Carlo w przestrzeni kątów (φ, ψ)' },
  { key: 'F', name: 'Kontakty (test kontrolny)', desc: 'UWAGA: NIE na ślepo. Kontakty z prawdziwej struktury = górna granica tego, co da koewolucja/DCA. Włącza też Monte Carlo.' }
];

// Kontakty CA–CA < 8 Å, |i-j| >= 6, z pliku PDB (tylko do testu kontrolnego F)
function nativeContacts(pdb: string, n: number): { i: number; j: number }[] {
  const ca: number[][] = []; const seen = new Set<number>();
  for (const l of (pdb || '').split('\n')) {
    if (!l.startsWith('ATOM') || l.slice(12, 16).trim() !== 'CA' || !(l[16] === ' ' || l[16] === 'A')) continue;
    const r = +l.slice(22, 26); if (seen.has(r)) continue; seen.add(r);
    ca.push([+l.slice(30, 38), +l.slice(38, 46), +l.slice(46, 54)]);
  }
  const m = Math.min(n, ca.length), out: { i: number; j: number }[] = [];
  for (let i = 0; i < m; i++) for (let j = i + 6; j < m; j++) if (Math.hypot(ca[i][0] - ca[j][0], ca[i][1] - ca[j][1], ca[i][2] - ca[j][2]) < 8) out.push({ i, j });
  return out;
}
const HISTORY_KEY = 'snaplock_historia_przebiegow_v1';

const PROTEINS: { key: string; label: string }[] = [
  { key: 'INSULIN_A_21', label: 'Insulina A (21)' },
  { key: 'CRAMBIN_46', label: 'Krambina (46)' },
  { key: 'UBIQUITIN_76', label: 'Ubikwityna (76)' },
  { key: 'LYSOZYME_129', label: 'Lizozym (129)' }
];

interface CellStats { rmsdMean: number; rmsdSd: number; tmMean: number; tmSd: number; best: number; n: number; }
interface GridRow { cfg: string; cells: Record<string, CellStats>; }
interface JournalEntry { date: string; text: string; }

const JOURNAL_KEY = 'snaplock_dziennik_wnioskow_v1';
const DEFAULT_JOURNAL: JournalEntry[] = [
  { date: '2026-11-08', text: 'Tabela 24 właściwości: większość wartości zgodna z literaturą. Poprawione: LYS objętość 135 Å³ (Richards), SASA wg Millera 1987, HIS ładunek 0,04, chiralność GLY. Wpływ na RMSD: w granicach szumu.' },
  { date: '2026-11-08', text: 'BŁĄD KRYTYCZNY naprawiony: budowa łańcucha (NeRF) liczyła ψ względem C(i-2) zamiast N(i-1). Odbudowa prawdziwej krambiny z jej kątów: 7,84 Å -> 1,33 Å.' },
  { date: '2026-11-08', text: 'Test energii: prawdziwa struktura ma NIŻSZĄ energię niż wynik silnika (krambina, ubikwityna, lizozym). Wzory energii idą w dobrą stronę; wąskie gardło = przeszukiwanie.' },
  { date: '2026-11-09', text: 'GOR (struktura drugorzędowa z 691 białek PDB): Q3 62,6% vs 53,5% stara reguła vs 46,1% Chou-Fasman. Sufit metod z jednej sekwencji ~60-65%. Do 3D samo GOR nie wystarcza.' },
  { date: '2026-11-09', text: 'Siatka A/AB/AC/.../ABCD, 5 startów: różnice między wariantami mniejsze niż rozrzut (1-4 Å). Najlepszy trend: lizozym ABCD 14,5±1,3 Å vs A 15,8±1,4 Å; ubikwityna ACD 11,2±1,6 vs A 14,1±2,3. Żaden moduł nie jest jeszcze przełomem. Następny krok: koewolucja (kontakty z rodzin białek).' },
  { date: '2026-11-09', text: 'Moduł F, test kontrolny: nawet PRAWDZIWE kontakty w Monte Carlo dają mało (krambina 8,7 -> 7,2 Å; ubikwityna 10,6 -> 9,4 Å; realistyczne DCA 70%: 7,9 i 10,2 Å). Monte Carlo po jednym kącie nie umie spełnić kontaktów daleko w łańcuchu.' },
  { date: '2026-11-09', text: 'Moduł G (geometria odległości, poza aplikacją): pełna macierz odległości -> 0,3 Å (optymalizator działa). Same kontakty ~7 Å. Kontakty + pary DALEKO + prawdziwa struktura drugorzędowa -> 4,4 Å. WNIOSEK: kolejność = 1) lepsza struktura drugorzędowa (profil z rodziny białek), 2) przewidywanie odległości, 3) dopiero potem DCA.' }
];

function allConfigs(selected: ModuleKey[]): string[] {
  const out: string[] = [];
  const n = selected.length;
  for (let mask = 0; mask < 1 << n; mask++) {
    out.push('A' + selected.filter((_, i) => mask & (1 << i)).join(''));
  }
  return out;
}

const mean = (a: number[]) => a.reduce((x, y) => x + y, 0) / Math.max(1, a.length);
const sd = (a: number[]) => { const m = mean(a); return Math.sqrt(mean(a.map(x => (x - m) ** 2))); };

export default function AblationLabView() {
  const [selectedModules, setSelectedModules] = useState<ModuleKey[]>(['B', 'C', 'D']);
  const [selectedProteins, setSelectedProteins] = useState<string[]>(['CRAMBIN_46', 'UBIQUITIN_76', 'LYSOZYME_129']);
  const [seeds, setSeeds] = useState(5);
  const [mcSteps, setMcSteps] = useState(5000);
  const [frames, setFrames] = useState(150);
  const [rows, setRows] = useState<GridRow[]>([]);
  const [running, setRunning] = useState(false);
  const [progress, setProgress] = useState('');
  const stopRef = useRef(false);

  const [journal, setJournal] = useState<JournalEntry[]>(() => {
    try { const raw = localStorage.getItem(JOURNAL_KEY); if (raw) return JSON.parse(raw); } catch { /* ignore */ }
    return DEFAULT_JOURNAL;
  });
  // dopisz nowe wnioski domyślne, których jeszcze nie ma w zapisanym dzienniku
  useEffect(() => { setJournal(j => { const missing = DEFAULT_JOURNAL.filter(d => !j.some(x => x.text === d.text)); return missing.length ? [...missing.reverse(), ...j] : j; }); }, []);
  const [newNote, setNewNote] = useState('');
  useEffect(() => { try { localStorage.setItem(JOURNAL_KEY, JSON.stringify(journal)); } catch { /* ignore */ } }, [journal]);

  const toggle = <T,>(list: T[], v: T) => (list.includes(v) ? list.filter(x => x !== v) : [...list, v]);

  const runGrid = async () => {
    stopRef.current = false;
    setRunning(true);
    setRows([]);
    const configs = allConfigs(selectedModules);
    const total = configs.length * selectedProteins.length * seeds;
    let done = 0;
    const result: GridRow[] = [];
    for (const cfg of configs) {
      const row: GridRow = { cfg, cells: {} };
      for (const key of selectedProteins) {
        const spec: any = (GLOBAL_PROTEIN_TARGETS as any)[key];
        if (!spec) continue;
        const r: number[] = [], t: number[] = [];
        for (let seed = 1; seed <= seeds; seed++) {
          if (stopRef.current) break;
          setProgress(`${cfg} · ${key} · start ${seed}/${seeds}  (${done}/${total})`);
          await new Promise(res => setTimeout(res, 0)); // oddaj sterowanie przeglądarce
          const run = runMultiAgentFoldingToConvergence(
            {
              sequenceCodes: spec.sequence3Letter,
              cdAlphaTargetPct: spec.cdAlphaTargetPct,
              cdBetaTargetPct: spec.cdBetaTargetPct,
              explicitDisulfidePairs: spec.explicitDisulfidePairs0Based ?? [],
              referencePdbText: spec.referencePdbText,
              referenceId: key,
              enableTargetGuidedRefinement: false,
              randomSeed: seed,
              initialNoiseDeg: 15,
              mjContactScale: cfg.includes('B') ? 1.0 : 0,
              hBondScale: cfg.includes('C') ? 1.5 : 0,
              secondaryStructurePredictor: cfg.includes('D') ? 'GOR' : 'LEGACY',
              monteCarloSteps: (cfg.includes('E') || cfg.includes('F')) ? mcSteps : 0,
              contactRestraints: cfg.includes('F') ? nativeContacts(spec.referencePdbText, spec.sequence3Letter.length) : undefined,
              contactWeight: 2.0
            },
            frames,
            3
          );
          const g: any = run.groundTruthValidation;
          if (g && g.caRmsdAngstrom != null) { r.push(g.caRmsdAngstrom); t.push(g.tmScore ?? 0); }
          done++;
        }
        row.cells[key] = { rmsdMean: mean(r), rmsdSd: sd(r), tmMean: mean(t), tmSd: sd(t), best: r.length ? Math.min(...r) : NaN, n: r.length };
      }
      result.push(row);
      setRows([...result]);
      if (stopRef.current) break;
    }
    setProgress(stopRef.current ? 'Zatrzymano.' : `Gotowe: ${result.length} wariantów × ${selectedProteins.length} białek × ${seeds} startów.`);
    // UCZENIE: każdy pełny przebieg zostaje w historii i dopisuje automatyczny wniosek do dziennika
    if (!stopRef.current && result.length) {
      const date = new Date().toISOString().slice(0, 10);
      try { const h = JSON.parse(localStorage.getItem(HISTORY_KEY) || '[]'); h.push({ date, seeds, mcSteps, frames, rows: result }); localStorage.setItem(HISTORY_KEY, JSON.stringify(h.slice(-200))); } catch { /* ignore */ }
      const a = result.find(r => r.cfg === 'A');
      const parts = selectedProteins.map(key => {
        const best = [...result].filter(r => r.cells[key]?.n).sort((x, y) => x.cells[key].rmsdMean - y.cells[key].rmsdMean)[0];
        if (!best) return '';
        const b = best.cells[key], base = a?.cells[key];
        let tag = '';
        if (base && best.cfg !== 'A') { const noise = Math.sqrt((base.rmsdSd ** 2 + b.rmsdSd ** 2) / Math.max(1, Math.min(base.n, b.n))); const diff = base.rmsdMean - b.rmsdMean; tag = diff > 2 * noise ? ' (LEPIEJ niż A)' : ' (w granicach szumu)'; }
        return `${key}: najlepszy ${best.cfg} ${b.rmsdMean.toFixed(2)}±${b.rmsdSd.toFixed(2)} Å` + (base ? ` vs A ${base.rmsdMean.toFixed(2)}` : '') + tag;
      }).filter(Boolean);
      setJournal(j => [{ date, text: `[auto] Siatka ${result.map(r => r.cfg).join('/')}, ${seeds} startów, MC ${mcSteps}: ` + parts.join('; ') }, ...j]);
    }
    setRunning(false);
  };

  const base = rows.find(r => r.cfg === 'A');
  const verdict = (row: GridRow, key: string) => {
    const a = base?.cells[key], c = row.cells[key];
    if (!a || !c || row.cfg === 'A') return '';
    const diff = a.rmsdMean - c.rmsdMean; // dodatnie = lepiej niż A
    const noise = Math.sqrt((a.rmsdSd ** 2 + c.rmsdSd ** 2) / Math.max(1, Math.min(a.n, c.n)));
    if (diff > 2 * noise) return 'lepiej';
    if (diff < -2 * noise) return 'gorzej';
    return 'szum';
  };

  // Synergia: efekt(XY) - [efekt(X) + efekt(Y)], gdzie efekt = poprawa RMSD względem A (Å). Dodatnia = moduły się wzmacniają.
  const synergy = (row: GridRow, key: string): number | null => {
    const mods = row.cfg.slice(1).split('');
    if (mods.length < 2 || !base?.cells[key] || !row.cells[key]) return null;
    const a = base.cells[key].rmsdMean;
    let sumSingles = 0;
    for (const m of mods) { const single = rows.find(r => r.cfg === 'A' + m)?.cells[key]; if (!single) return null; sumSingles += a - single.rmsdMean; }
    return (a - row.cells[key].rmsdMean) - sumSingles;
  };

  const exportCsv = () => {
    const head = ['wariant', ...selectedProteins.flatMap(k => [`${k}_rmsd`, `${k}_sd`, `${k}_tm`])].join(',');
    const lines = rows.map(r => [r.cfg, ...selectedProteins.flatMap(k => { const c = r.cells[k]; return c ? [c.rmsdMean.toFixed(3), c.rmsdSd.toFixed(3), c.tmMean.toFixed(4)] : ['', '', '']; })].join(','));
    const blob = new Blob([[head, ...lines].join('\n')], { type: 'text/csv' });
    const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = 'siatka_wariantow.csv'; a.click();
  };

  return (
    <div className="flex-1 w-full h-full overflow-auto p-4 text-zinc-200 font-mono text-[12px] space-y-4">
      <div>
        <h2 className="text-lg font-bold text-white">J: Laboratorium wariantów A/B</h2>
        <p className="text-zinc-400">A = ścieżka bazowa. Każdy moduł dokładany osobno i w zestawach. Każdy wariant z kilku startów — „lepiej/gorzej” tylko gdy różnica &gt; 2× szum. „syn” = synergia: czy moduły razem dają więcej (+) czy mniej (−) niż suma osobno. Wszystko na ślepo.</p>
      </div>

      <div className="grid md:grid-cols-3 gap-3">
        <div className="border border-zinc-800 rounded p-3 space-y-1">
          <div className="font-bold text-white mb-1">Moduły</div>
          {MODULES.map(m => (
            <label key={m.key} className="flex gap-2 items-start cursor-pointer">
              <input type="checkbox" checked={selectedModules.includes(m.key)} onChange={() => setSelectedModules(toggle(selectedModules, m.key))} disabled={running} />
              <span><b>{m.key}</b>: {m.name} <span className="text-zinc-500">— {m.desc}</span></span>
            </label>
          ))}
        </div>
        <div className="border border-zinc-800 rounded p-3 space-y-1">
          <div className="font-bold text-white mb-1">Białka</div>
          {PROTEINS.map(p => (
            <label key={p.key} className="flex gap-2 items-center cursor-pointer">
              <input type="checkbox" checked={selectedProteins.includes(p.key)} onChange={() => setSelectedProteins(toggle(selectedProteins, p.key))} disabled={running} />
              <span>{p.label}</span>
            </label>
          ))}
        </div>
        <div className="border border-zinc-800 rounded p-3 space-y-2">
          <div className="font-bold text-white">Ustawienia</div>
          <label className="flex justify-between gap-2">Starty na wariant <input className="w-16 bg-zinc-900 px-1" type="number" min={1} max={20} value={seeds} onChange={e => setSeeds(+e.target.value)} disabled={running} /></label>
          <label className="flex justify-between gap-2">Klatki symulacji <input className="w-16 bg-zinc-900 px-1" type="number" min={30} max={2000} value={frames} onChange={e => setFrames(+e.target.value)} disabled={running} /></label>
          <label className="flex justify-between gap-2">Kroki Monte Carlo (E) <input className="w-20 bg-zinc-900 px-1" type="number" min={0} max={200000} step={1000} value={mcSteps} onChange={e => setMcSteps(+e.target.value)} disabled={running} /></label>
          <div className="flex gap-2 pt-1">
            {!running ? (
              <button className="px-3 py-1 rounded bg-emerald-700 hover:bg-emerald-600 text-white font-bold" onClick={runGrid}>Uruchom siatkę ({allConfigs(selectedModules).length} wariantów)</button>
            ) : (
              <button className="px-3 py-1 rounded bg-rose-700 hover:bg-rose-600 text-white font-bold" onClick={() => (stopRef.current = true)}>Stop</button>
            )}
            {rows.length > 0 && !running && <button className="px-3 py-1 rounded bg-zinc-800 hover:bg-zinc-700" onClick={exportCsv}>CSV</button>}
          </div>
          <div className="text-zinc-400 min-h-[1.2em]">{progress}</div>
        </div>
      </div>

      {rows.length > 0 && (
        <div className="overflow-auto border border-zinc-800 rounded">
          <table className="w-full text-left">
            <thead className="bg-zinc-900 text-zinc-300">
              <tr>
                <th className="px-2 py-1">Wariant</th>
                {selectedProteins.map(k => <th key={k} className="px-2 py-1">{PROTEINS.find(p => p.key === k)?.label ?? k}<br /><span className="text-zinc-500 font-normal">RMSD Å (śr ± odch.) · TM · najlepszy</span></th>)}
              </tr>
            </thead>
            <tbody>
              {rows.map(r => (
                <tr key={r.cfg} className="border-t border-zinc-800">
                  <td className="px-2 py-1 font-bold text-white">{r.cfg}</td>
                  {selectedProteins.map(k => {
                    const c = r.cells[k]; const v = verdict(r, k);
                    return (
                      <td key={k} className="px-2 py-1">
                        {c ? <>
                          {c.rmsdMean.toFixed(2)} ± {c.rmsdSd.toFixed(2)} · {c.tmMean.toFixed(3)} · {c.best.toFixed(2)}
                          {v && <span className={`ml-2 px-1 rounded ${v === 'lepiej' ? 'bg-emerald-800 text-emerald-200' : v === 'gorzej' ? 'bg-rose-800 text-rose-200' : 'bg-zinc-800 text-zinc-400'}`}>{v}</span>}
                          {(() => { const syn = synergy(r, k); return syn === null ? null : <span className={`ml-2 ${syn > 0 ? 'text-emerald-300' : 'text-rose-300'}`} title="Synergia: o ile Å kombinacja jest lepsza (+) lub gorsza (−) niż suma pojedynczych modułów">syn {syn > 0 ? '+' : ''}{syn.toFixed(2)}</span>; })()}
                        </> : '—'}
                      </td>
                    );
                  })}
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}

      <div className="border border-zinc-800 rounded p-3 space-y-2">
        <div className="font-bold text-white">Dziennik wniosków (zapisujemy, nie wdrażamy bez dowodu)</div>
        <div className="flex gap-2">
          <input className="flex-1 bg-zinc-900 px-2 py-1 rounded" placeholder="Nowy wniosek z badania…" value={newNote} onChange={e => setNewNote(e.target.value)} />
          <button className="px-3 py-1 rounded bg-indigo-700 hover:bg-indigo-600 text-white" onClick={() => { if (newNote.trim()) { setJournal([{ date: new Date().toISOString().slice(0, 10), text: newNote.trim() }, ...journal]); setNewNote(''); } }}>Dodaj</button>
        </div>
        <ul className="space-y-1">
          {journal.map((j, i) => (
            <li key={i} className="border-t border-zinc-800 pt-1"><span className="text-zinc-500">{j.date}</span> — {j.text}</li>
          ))}
        </ul>
      </div>
    </div>
  );
}
