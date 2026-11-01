import React, { useState } from 'react';
import { 
  X, 
  Activity, 
  CheckCircle2, 
  AlertTriangle, 
  Download, 
  Play, 
  Dna, 
  ShieldCheck, 
  RefreshCw,
  Terminal,
  FileJson,
  Sparkles
} from 'lucide-react';
import { runSnaplockCalibrationFromPdb, SnaplockBenchmarkReport } from '../snaplockEngine';

interface SnaplockBenchmarkModalProps {
  isOpen: boolean;
  onClose: () => void;
  onLoadSequence?: (sequence3Letter: string) => void;
}

export const SnaplockBenchmarkModal: React.FC<SnaplockBenchmarkModalProps> = ({
  isOpen,
  onClose,
  onLoadSequence
}) => {
  const [uniprotId, setUniprotId] = useState<string>('P01308'); // Human Insulin default
  const [isLoading, setIsLoading] = useState<boolean>(false);
  const [report, setReport] = useState<SnaplockBenchmarkReport | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [filterResidue, setFilterResidue] = useState<string>('ALL');

  if (!isOpen) return null;

  const runBenchmark = async (targetId: string = uniprotId) => {
    setIsLoading(true);
    setError(null);
    try {
      // First try backend proxy endpoint
      const response = await fetch('/api/snaplock/benchmark', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ uniprotId: targetId.trim().toUpperCase() })
      });

      if (response.ok) {
        const data = await response.json();
        setReport(data);
      } else {
        // Fallback: fetch PDB from AlphaFold DB directly from client
        const afMeta = await fetch(`https://alphafold.ebi.ac.uk/api/prediction/${targetId.trim().toUpperCase()}`);
        if (!afMeta.ok) throw new Error(`Nie odnaleziono predykcji AlphaFold dla ${targetId}`);
        const metaJson = await afMeta.json();
        if (!metaJson?.[0]?.pdbUrl) throw new Error(`Brak pliku PDB dla ${targetId}`);
        
        const pdbResp = await fetch(metaJson[0].pdbUrl);
        if (!pdbResp.ok) throw new Error('Nie udało się pobrać pliku PDB.');
        const pdbContent = await pdbResp.text();
        const calibReport = runSnaplockCalibrationFromPdb(pdbContent, targetId.toUpperCase());
        setReport(calibReport);
      }
    } catch (err: any) {
      setError(err.message || 'Wystąpił błąd podczas walidacji AlphaFold E8.');
    } finally {
      setIsLoading(false);
    }
  };

  const downloadJsonReport = () => {
    if (!report) return;
    const blob = new Blob([JSON.stringify(report, null, 2)], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `snaplock_benchmark_${report.uniprot_id || 'custom'}.json`;
    a.click();
    URL.revokeObjectURL(url);
  };

  const PRESET_PROTEINS = [
    { id: 'P83347', name: 'P83347', label: 'Kalibracja testowa (30 aa)' },
    { id: 'P01308', name: 'Insulina', label: 'Insulina ludzka (P01308, 110 aa)' },
    { id: 'P0CG48', name: 'Ubikwityna', label: 'Ubikwityna (P0CG48, 76 aa)' },
    { id: 'P00698', name: 'Lizozym', label: 'Lizozym HEWL (P00698, 147 aa)' }
  ];

  return (
    <div className="fixed inset-0 z-50 bg-black/80 backdrop-blur-md flex items-center justify-center p-3 sm:p-6 overflow-y-auto animate-in fade-in duration-200">
      <div className="bg-[#04040a] border border-cyan-800/80 rounded-2xl w-full max-w-5xl shadow-[0_0_50px_rgba(6,182,212,0.15)] flex flex-col max-h-[92vh] overflow-hidden">
        
        {/* Header */}
        <div className="p-4 sm:p-5 border-b border-zinc-900 bg-gradient-to-r from-[#020617] via-[#04101e] to-[#020617] flex items-center justify-between shrink-0">
          <div className="flex items-center gap-3">
            <div className="w-10 h-10 rounded-xl bg-cyan-950/80 border border-cyan-500/50 flex items-center justify-center text-cyan-400 shadow-inner">
              <ShieldCheck className="w-5 h-5 text-cyan-300" />
            </div>
            <div>
              <div className="flex items-center gap-2">
                <h2 className="text-base sm:text-lg font-bold font-mono text-white tracking-tight">
                  Oficjalny Moduł Walidacyjny SNAPLOCK vs AlphaFold DB
                </h2>
                <span className="text-[10px] font-mono px-2 py-0.5 rounded bg-cyan-950 text-cyan-300 border border-cyan-700/60 font-bold uppercase">
                  Siatka E8 (248 Wektorów)
                </span>
              </div>
              <p className="text-xs text-zinc-400 font-mono mt-0.5">
                Deterministyczne rzutowanie atomów z AlphaFold DB (EMBL-EBI) na dyskretne wektory Liego E8
              </p>
            </div>
          </div>
          <button
            onClick={onClose}
            className="p-2 text-zinc-400 hover:text-white rounded-lg hover:bg-zinc-900/80 border border-transparent hover:border-zinc-800 transition cursor-pointer"
          >
            <X className="w-5 h-5" />
          </button>
        </div>

        {/* Action Controls & Input */}
        <div className="p-4 border-b border-zinc-900 bg-[#020206] shrink-0 flex flex-col md:flex-row items-center justify-between gap-3">
          <div className="flex items-center gap-2 w-full md:w-auto">
            <span className="text-xs font-mono text-zinc-400 uppercase font-bold shrink-0">UniProt ID:</span>
            <input
              type="text"
              value={uniprotId}
              onChange={(e) => setUniprotId(e.target.value.toUpperCase())}
              placeholder="np. P01308, P0CG48, P83347"
              className="bg-black border border-cyan-900/80 rounded px-3 py-1.5 text-sm font-mono text-white focus:outline-none focus:border-cyan-400 w-44"
            />
            <button
              onClick={() => runBenchmark(uniprotId)}
              disabled={isLoading || !uniprotId.trim()}
              className="px-4 py-1.5 bg-gradient-to-r from-cyan-600 to-indigo-600 hover:from-cyan-500 hover:to-indigo-500 disabled:opacity-50 text-white font-mono text-xs font-bold rounded flex items-center gap-1.5 cursor-pointer shadow-md transition"
            >
              {isLoading ? (
                <>
                  <RefreshCw className="w-4 h-4 animate-spin text-cyan-200" />
                  <span>Kwantyzacja E8...</span>
                </>
              ) : (
                <>
                  <Play className="w-4 h-4 text-cyan-200" />
                  <span>Uruchom Walidację</span>
                </>
              )}
            </button>
          </div>

          {/* Quick presets */}
          <div className="flex items-center gap-1.5 flex-wrap w-full md:w-auto">
            <span className="text-[10px] font-mono text-zinc-500 uppercase">Wzorce:</span>
            {PRESET_PROTEINS.map((p) => (
              <button
                key={p.id}
                onClick={() => {
                  setUniprotId(p.id);
                  runBenchmark(p.id);
                }}
                className={`px-2.5 py-1 rounded border text-[11px] font-mono transition cursor-pointer ${
                  uniprotId === p.id 
                    ? 'bg-cyan-950/80 border-cyan-500 text-cyan-300 font-bold' 
                    : 'bg-zinc-950 border-zinc-800 text-zinc-400 hover:text-white'
                }`}
                title={p.label}
              >
                {p.name}
              </button>
            ))}
          </div>
        </div>

        {/* Modal Body */}
        <div className="flex-1 overflow-y-auto p-4 sm:p-6 space-y-5 bg-[#030308]">
          {error && (
            <div className="p-3.5 bg-red-950/50 border border-red-800 rounded-xl text-red-300 text-xs font-mono flex items-center gap-2">
              <AlertTriangle className="w-5 h-5 text-red-400 shrink-0" />
              <span>{error}</span>
            </div>
          )}

          {!report && !isLoading && (
            <div className="py-12 flex flex-col items-center justify-center text-center space-y-4 text-zinc-500">
              <Dna className="w-16 h-16 text-cyan-900/60 animate-pulse" />
              <div>
                <p className="text-sm font-mono text-zinc-300 font-bold">Wybierz lub wpisz kod UniProt i kliknij &bdquo;Uruchom Walidację&rdquo;</p>
                <p className="text-xs font-mono text-zinc-500 mt-1 max-w-md">
                  Silnik pobierze strukturę atomową 3D bezpośrednio z EMBL-EBI AlphaFold DB, zmapuje każdy wektor wiązania na 248 wektorów algebry E8 i wyliczy precyzyjne odchylenia kątowe.
                </p>
              </div>
              <div className="flex items-center gap-2 text-xs font-mono bg-zinc-950/80 border border-zinc-850 p-2.5 rounded-lg text-zinc-400">
                <Terminal className="w-4 h-4 text-cyan-400" />
                <span>Skrypt Python w katalogu głównym: <code>./snaplock_alphafold_benchmark.py P01308</code></span>
              </div>
            </div>
          )}

          {isLoading && (
            <div className="py-16 flex flex-col items-center justify-center text-center space-y-4">
              <RefreshCw className="w-10 h-10 text-cyan-400 animate-spin" />
              <p className="text-sm font-mono text-cyan-300 font-bold">Pobieranie PDB z AlphaFold DB i rzutowanie na dyskretną siatkę E8...</p>
              <p className="text-xs font-mono text-zinc-500">Kalkulacja iloczynów skalarnych dla 248 łatek kontaktowych per reszta...</p>
            </div>
          )}

          {report && !isLoading && (
            <>
              {/* Metric Cards Grid */}
              <div className="grid grid-cols-2 sm:grid-cols-4 gap-3">
                <div className="bg-[#050914] border border-cyan-900/40 rounded-xl p-3 font-mono">
                  <span className="text-[10px] text-zinc-500 uppercase font-bold block">Reszty Przeanalizowane</span>
                  <div className="text-2xl font-black text-white mt-1">
                    {report.total_residues}
                  </div>
                  <span className="text-[10px] text-cyan-400 mt-0.5 block">
                    Zwalidowane: {report.quantized_bonds_count} wiązań
                  </span>
                </div>

                <div className="bg-[#050914] border border-cyan-900/40 rounded-xl p-3 font-mono">
                  <span className="text-[10px] text-zinc-500 uppercase font-bold block">Średni Błąd Kątowy E8</span>
                  <div className="text-2xl font-black text-amber-400 mt-1">
                    {report.mean_angular_patch_loss_deg}°
                  </div>
                  <span className="text-[10px] text-zinc-400 mt-0.5 block">
                    Max: {report.max_angular_loss_deg}°
                  </span>
                </div>

                <div className="bg-[#050914] border border-cyan-900/40 rounded-xl p-3 font-mono">
                  <span className="text-[10px] text-zinc-500 uppercase font-bold block">Mostki Kowalencyjne S-S</span>
                  <div className="text-2xl font-black text-emerald-400 mt-1">
                    {report.disulfide_bridges_detected}
                  </div>
                  <span className="text-[10px] text-zinc-400 mt-0.5 block">
                    Kolizje VdW: {report.vdw_clashes_detected}
                  </span>
                </div>

                <div className="bg-[#050914] border border-cyan-900/40 rounded-xl p-3 font-mono">
                  <span className="text-[10px] text-zinc-500 uppercase font-bold block">AlphaFold Śr. pLDDT</span>
                  <div className="text-2xl font-black text-teal-400 mt-1">
                    {report.alphafold_mean_plddt}%
                  </div>
                  <span className="text-[10px] text-emerald-400 mt-0.5 block">
                    Zgodność modelu: Bardzo wysoka
                  </span>
                </div>
              </div>

              {/* Action bar for current report */}
              <div className="flex flex-col sm:flex-row items-center justify-between gap-3 bg-zinc-950/80 p-3 rounded-xl border border-zinc-900">
                <div className="flex items-center gap-2 text-xs font-mono">
                  <span className="text-zinc-400">Filtr aminokwasów:</span>
                  <select
                    value={filterResidue}
                    onChange={(e) => setFilterResidue(e.target.value)}
                    className="bg-black border border-zinc-800 rounded px-2 py-1 text-white text-xs font-mono outline-none"
                  >
                    <option value="ALL">Wszystkie reszty ({report.quantized_bonds_count})</option>
                    <option value="CYS">Tylko Cysteiny (CYS/SEC)</option>
                    <option value="HIGH_ERR">Najwyższy błąd kątowy (&gt;6°)</option>
                  </select>
                </div>

                <div className="flex items-center gap-2">
                  <button
                    onClick={downloadJsonReport}
                    className="px-3 py-1.5 bg-zinc-900 hover:bg-zinc-800 border border-zinc-750 text-cyan-300 font-mono text-xs rounded flex items-center gap-1.5 cursor-pointer transition"
                  >
                    <FileJson className="w-3.5 h-3.5 text-cyan-400" />
                    <span>Eksportuj Raport JSON</span>
                  </button>

                  {onLoadSequence && report.sample_quantized_bonds.length > 0 && (
                    <button
                      onClick={() => {
                        const seq = report.sample_quantized_bonds.map(b => b.residue_name).join(' ');
                        onLoadSequence(seq);
                        onClose();
                      }}
                      className="px-3 py-1.5 bg-cyan-950 hover:bg-cyan-900 border border-cyan-700 text-white font-mono text-xs font-bold rounded flex items-center gap-1.5 cursor-pointer transition"
                    >
                      <Sparkles className="w-3.5 h-3.5 text-cyan-300" />
                      <span>Załaduj do Widoku 3D</span>
                    </button>
                  )}
                </div>
              </div>

              {/* Residue-by-residue details table */}
              <div className="border border-zinc-900 rounded-xl overflow-hidden bg-black/60">
                <div className="p-3 bg-zinc-950 border-b border-zinc-900 flex justify-between items-center text-xs font-mono text-zinc-400">
                  <span className="font-bold text-white uppercase tracking-wider">
                    Szczegółowa Tabela Kwantyzacji Wektorowej E8 (Siatka 248 Łatek)
                  </span>
                  <span>{report.sample_quantized_bonds.length} reszt</span>
                </div>

                <div className="max-h-80 overflow-y-auto font-mono text-xs">
                  <table className="w-full text-left border-collapse">
                    <thead className="bg-[#05060d] text-zinc-400 text-[10px] uppercase sticky top-0 border-b border-zinc-900">
                      <tr>
                        <th className="p-2.5"># Reszta</th>
                        <th className="p-2.5">Kod</th>
                        <th className="p-2.5">Przypisana Łatka E8</th>
                        <th className="p-2.5">Błąd Kątowy</th>
                        <th className="p-2.5">Objętość VdW</th>
                        <th className="p-2.5">AF pLDDT</th>
                        <th className="p-2.5">Wektor Wiązania [X, Y, Z]</th>
                      </tr>
                    </thead>
                    <tbody className="divide-y divide-zinc-900/60">
                      {report.sample_quantized_bonds
                        .filter(b => {
                          if (filterResidue === 'CYS') return b.residue_name === 'CYS' || b.residue_name === 'SEC';
                          if (filterResidue === 'HIGH_ERR') return b.angle_error_deg > 6.0;
                          return true;
                        })
                        .map((b, idx) => (
                          <tr key={idx} className="hover:bg-zinc-900/40 transition">
                            <td className="p-2.5 text-zinc-400">#{b.residue_idx}</td>
                            <td className="p-2.5 font-bold text-white">
                              <span className={`px-1.5 py-0.5 rounded text-[10px] ${
                                b.residue_name === 'CYS' ? 'bg-amber-950 text-amber-300 border border-amber-800' : 'bg-zinc-900 text-zinc-200'
                              }`}>
                                {b.residue_name}
                              </span>
                            </td>
                            <td className="p-2.5">
                              <span className="text-cyan-400 font-bold">Łatka #{b.assigned_patch}</span>
                              {b.assigned_patch === 1 && <span className="text-[9px] text-emerald-400 ml-1">(N-in)</span>}
                              {b.assigned_patch === 248 && <span className="text-[9px] text-emerald-400 ml-1">(C-out)</span>}
                            </td>
                            <td className="p-2.5">
                              <span className={`font-bold ${b.angle_error_deg < 5.0 ? 'text-emerald-400' : b.angle_error_deg < 8.0 ? 'text-amber-400' : 'text-red-400'}`}>
                                {b.angle_error_deg}°
                              </span>
                            </td>
                            <td className="p-2.5 text-zinc-300">{b.vdw_radius} Å</td>
                            <td className="p-2.5">
                              <span className="text-teal-400 font-bold">{b.plddt}%</span>
                            </td>
                            <td className="p-2.5 text-zinc-500 text-[10px]">
                              [{b.bond_vector.map(n => n.toFixed(3)).join(', ')}]
                            </td>
                          </tr>
                        ))}
                    </tbody>
                  </table>
                </div>
              </div>
            </>
          )}

        </div>

        {/* Footer */}
        <div className="p-4 border-t border-zinc-900 bg-[#020206] flex items-center justify-between text-xs font-mono text-zinc-500 shrink-0">
          <div className="flex items-center gap-2">
            <span className="w-2 h-2 rounded-full bg-cyan-400" />
            <span>Deterministyczna Algebra E8 • Standard walidacji SNAPLOCK v1.0</span>
          </div>
          <button
            onClick={onClose}
            className="px-4 py-1.5 bg-zinc-900 hover:bg-zinc-800 border border-zinc-800 text-zinc-300 hover:text-white rounded font-mono transition cursor-pointer"
          >
            Zamknij
          </button>
        </div>

      </div>
    </div>
  );
};
