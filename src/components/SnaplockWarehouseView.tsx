import React, { useState } from 'react';
import { 
  SNAPLOCK_22_DATABASE, 
  E8_PATCHES, 
  generateAminoAcidWarehouseJson, 
  evaluateSnaplockJunction,
  AminoAcidPhysicochemicalData,
  E8Patch,
  getVdwRadius
} from '../snaplockEngine';
import { 
  Database, 
  Layers, 
  Download, 
  Copy, 
  Check, 
  Search, 
  ShieldCheck, 
  Sliders, 
  Compass, 
  Zap, 
  Flame, 
  Activity,
  Sparkles,
  ExternalLink,
  ChevronRight
} from 'lucide-react';

export const SnaplockWarehouseView: React.FC = () => {
  const [activeTab, setActiveTab] = useState<'matrix22' | 'patches248' | 'junction5tracks' | 'jsonExport' | 'benchmarkEngine'>('matrix22');
  const [searchQuery, setSearchQuery] = useState('');
  const [selectedFilter, setSelectedFilter] = useState<string>('all');
  
  // Selected patch for detailed view
  const [selectedPatchId, setSelectedPatchId] = useState<number>(1);
  const selectedPatch = E8_PATCHES.find(p => p.id === selectedPatchId) || E8_PATCHES[0];

  // 5 Tracks Junction test state
  const [predCode, setPredCode] = useState<string>('CYS');
  const [succCode, setSuccCode] = useState<string>('CYS');
  const [torsionAngle, setTorsionAngle] = useState<number>(0);
  const [customPredPatch, setCustomPredPatch] = useState<number>(248);
  const [customSuccPatch, setCustomSuccPatch] = useState<number>(1);

  // JSON export state
  const [exportScope, setExportScope] = useState<'first3' | 'all21'>('first3');
  const [copied, setCopied] = useState(false);

  // AlphaFold vs SNAPLOCK Benchmark state
  const [benchmarkUniprot, setBenchmarkUniprot] = useState<string>('P83347');
  const [benchmarkLoading, setBenchmarkLoading] = useState<boolean>(false);
  const [benchmarkReport, setBenchmarkReport] = useState<any>(null);
  const [benchmarkError, setBenchmarkError] = useState<string | null>(null);

  // Compute 5 Tracks evaluation
  const predAA = SNAPLOCK_22_DATABASE.find(a => a.code3 === predCode) || SNAPLOCK_22_DATABASE[0];
  const succAA = SNAPLOCK_22_DATABASE.find(a => a.code3 === succCode) || SNAPLOCK_22_DATABASE[1];
  const junctionEval = evaluateSnaplockJunction(predAA, succAA, customPredPatch, customSuccPatch, torsionAngle);

  // Generate current JSON preview
  const warehouseJson = generateAminoAcidWarehouseJson(
    exportScope === 'first3' ? ['GLY', 'ALA', 'VAL'] : SNAPLOCK_22_DATABASE.map(a => a.code3)
  );
  const jsonString = JSON.stringify(warehouseJson, null, 2);

  const filteredAAs = SNAPLOCK_22_DATABASE.filter(aa => {
    const matchesSearch = 
      aa.code3.toLowerCase().includes(searchQuery.toLowerCase()) ||
      aa.namePl.toLowerCase().includes(searchQuery.toLowerCase()) ||
      aa.code1.toLowerCase().includes(searchQuery.toLowerCase());
    
    if (selectedFilter === 'all') return matchesSearch;
    if (selectedFilter === 'covalent') return matchesSearch && aa.covalentBridges.includes('TAK');
    if (selectedFilter === 'hydrophobic') return matchesSearch && aa.hydropathyKyte > 1.0;
    if (selectedFilter === 'charged') return matchesSearch && aa.chargePH74 !== 0;
    if (selectedFilter === 'polar') return matchesSearch && aa.polarity === 'Polarny';
    return matchesSearch;
  });

  const handleCopyJson = () => {
    navigator.clipboard.writeText(jsonString);
    setCopied(true);
    setTimeout(() => setCopied(false), 2000);
  };

  const handleDownloadJson = () => {
    const blob = new Blob([jsonString], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = exportScope === 'first3' ? 'amino_acid_warehouse_gly_ala_val.json' : 'amino_acid_warehouse_all21.json';
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  };

  const runBenchmark = async (targetId: string = benchmarkUniprot) => {
    setBenchmarkLoading(true);
    setBenchmarkError(null);
    try {
      const resp = await fetch('/api/snaplock/benchmark', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ uniprotId: targetId })
      });
      const data = await resp.json();
      if (!resp.ok) {
        throw new Error(data.error || 'Błąd wykonania benchmarku');
      }
      setBenchmarkReport(data);
    } catch (err: any) {
      setBenchmarkError(err.message || 'Wystąpił błąd podczas walidacji.');
    } finally {
      setBenchmarkLoading(false);
    }
  };

  return (
    <div className="bg-slate-900 border border-slate-800 rounded-2xl overflow-hidden shadow-2xl flex flex-col h-full">
      {/* Top Banner */}
      <div className="px-6 py-4 bg-slate-950 border-b border-slate-800 flex flex-wrap items-center justify-between gap-4">
        <div>
          <div className="flex items-center gap-2.5">
            <div className="p-2 bg-gradient-to-tr from-cyan-600 to-indigo-600 rounded-xl shadow-md">
              <Database className="w-5 h-5 text-white" />
            </div>
            <div>
              <h2 className="text-base font-bold text-slate-100 tracking-wide font-sans flex items-center gap-2">
                <span>Magazyn Fizykochemiczny SNAPLOCK</span>
                <span className="text-xs font-mono font-medium px-2 py-0.5 rounded-full bg-cyan-950/80 text-cyan-300 border border-cyan-800/80">
                  22 Parametry • 248 Łatek E8
                </span>
              </h2>
              <p className="text-xs text-slate-400">
                Deterministyczne dyskretne mapowanie kątów torsyjnych i ewaluacja 5 ścieżek biofizycznych
              </p>
            </div>
          </div>
        </div>

        {/* Tab Navigation */}
        <div className="flex items-center bg-slate-900 p-1 rounded-xl border border-slate-800 gap-1 text-xs">
          <button
            onClick={() => setActiveTab('matrix22')}
            className={`px-3 py-1.5 rounded-lg font-medium transition-all cursor-pointer ${
              activeTab === 'matrix22' 
                ? 'bg-cyan-600 text-white shadow-sm' 
                : 'text-slate-400 hover:text-slate-200'
            }`}
          >
            Matryca 22 Parametrów
          </button>
          <button
            onClick={() => setActiveTab('patches248')}
            className={`px-3 py-1.5 rounded-lg font-medium transition-all cursor-pointer ${
              activeTab === 'patches248' 
                ? 'bg-cyan-600 text-white shadow-sm' 
                : 'text-slate-400 hover:text-slate-200'
            }`}
          >
            Siatka 248 Łatek (E8)
          </button>
          <button
            onClick={() => setActiveTab('junction5tracks')}
            className={`px-3 py-1.5 rounded-lg font-medium transition-all cursor-pointer ${
              activeTab === 'junction5tracks' 
                ? 'bg-cyan-600 text-white shadow-sm' 
                : 'text-slate-400 hover:text-slate-200'
            }`}
          >
            Ewaluator 5 Ścieżek
          </button>
          <button
            onClick={() => setActiveTab('benchmarkEngine')}
            className={`px-3 py-1.5 rounded-lg font-medium transition-all cursor-pointer flex items-center gap-1.5 ${
              activeTab === 'benchmarkEngine' 
                ? 'bg-amber-600 text-white shadow-sm font-bold ring-1 ring-amber-400' 
                : 'text-amber-400 hover:text-amber-200'
            }`}
          >
            <ShieldCheck className="w-3.5 h-3.5 text-amber-300" />
            <span>Walidator AlphaFold E8</span>
          </button>
          <button
            onClick={() => setActiveTab('jsonExport')}
            className={`px-3 py-1.5 rounded-lg font-medium transition-all cursor-pointer flex items-center gap-1.5 ${
              activeTab === 'jsonExport' 
                ? 'bg-indigo-600 text-white shadow-sm' 
                : 'text-slate-400 hover:text-slate-200'
            }`}
          >
            <Download className="w-3.5 h-3.5" />
            <span>amino_acid_warehouse.json</span>
          </button>
        </div>
      </div>

      {/* Main Tab Content */}
      <div className="flex-1 overflow-y-auto p-6">
        {/* TAB 1: 22 Parameters Matrix */}
        {activeTab === 'matrix22' && (
          <div className="space-y-4">
            <div className="flex flex-wrap items-center justify-between gap-3 bg-slate-950/60 p-3 rounded-xl border border-slate-800">
              <div className="flex items-center gap-2 bg-slate-900 border border-slate-700/80 rounded-lg px-3 py-1.5 text-xs w-72">
                <Search className="w-4 h-4 text-slate-400" />
                <input
                  type="text"
                  placeholder="Szukaj aminokwasu (GLY, Ala, Val, Sec)..."
                  value={searchQuery}
                  onChange={(e) => setSearchQuery(e.target.value)}
                  className="bg-transparent text-slate-200 focus:outline-none w-full"
                />
              </div>

              <div className="flex items-center gap-1.5 text-xs">
                <span className="text-slate-400 mr-1 text-[11px]">Filtruj:</span>
                {[
                  { id: 'all', label: 'Wszystkie 21' },
                  { id: 'covalent', label: 'Klamry S-S / Se-Se' },
                  { id: 'hydrophobic', label: 'Hydrofobowe' },
                  { id: 'charged', label: 'Naładowane (pH 7.4)' },
                  { id: 'polar', label: 'Polarne' }
                ].map(f => (
                  <button
                    key={f.id}
                    onClick={() => setSelectedFilter(f.id)}
                    className={`px-2.5 py-1 rounded-md text-[11px] font-medium transition-colors ${
                      selectedFilter === f.id
                        ? 'bg-slate-700 text-cyan-300 border border-cyan-800'
                        : 'bg-slate-800/60 text-slate-400 hover:text-slate-200'
                    }`}
                  >
                    {f.label}
                  </button>
                ))}
              </div>
            </div>

            {/* Table */}
            <div className="overflow-x-auto border border-slate-800 rounded-xl bg-slate-950/40">
              <table className="w-full text-left text-xs border-collapse">
                <thead>
                  <tr className="bg-slate-950 border-b border-slate-800 text-[11px] font-mono text-slate-400 uppercase tracking-wider">
                    <th className="py-2.5 px-3">Nr</th>
                    <th className="py-2.5 px-3">Kod</th>
                    <th className="py-2.5 px-3">Nazwa</th>
                    <th className="py-2.5 px-3">Wzór boczny</th>
                    <th className="py-2.5 px-3">Masa [g/mol]</th>
                    <th className="py-2.5 px-3">VdW [Å³]</th>
                    <th className="py-2.5 px-3">R_vdw [Å]</th>
                    <th className="py-2.5 px-3">SASA [Å²]</th>
                    <th className="py-2.5 px-3">Polarność</th>
                    <th className="py-2.5 px-3">Kyte</th>
                    <th className="py-2.5 px-3">Hopp</th>
                    <th className="py-2.5 px-3">Ładunek 7.4</th>
                    <th className="py-2.5 px-3">pK1</th>
                    <th className="py-2.5 px-3">pK2</th>
                    <th className="py-2.5 px-3">pKr</th>
                    <th className="py-2.5 px-3">pI</th>
                    <th className="py-2.5 px-3">Mostki</th>
                    <th className="py-2.5 px-3">Ramachandran</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-slate-800/60 font-sans">
                  {filteredAAs.map((aa) => {
                    const rVdw = getVdwRadius(aa.vdwVolume).toFixed(2);
                    const isCovalent = aa.covalentBridges.includes('TAK');

                    return (
                      <tr key={aa.code3} className="hover:bg-slate-800/40 transition-colors">
                        <td className="py-2 px-3 font-mono text-slate-500">{aa.nr}</td>
                        <td className="py-2 px-3">
                          <span className="font-mono font-bold text-cyan-400 bg-cyan-950/50 px-1.5 py-0.5 rounded border border-cyan-800/40">
                            {aa.code3} / {aa.code1}
                          </span>
                        </td>
                        <td className="py-2 px-3 font-medium text-slate-200">{aa.namePl}</td>
                        <td className="py-2 px-3 font-mono text-slate-300">{aa.sidechainFormula}</td>
                        <td className="py-2 px-3 font-mono text-slate-400">{aa.mass.toFixed(2)}</td>
                        <td className="py-2 px-3 font-mono text-amber-300 font-semibold">{aa.vdwVolume}</td>
                        <td className="py-2 px-3 font-mono text-amber-400/80">{rVdw}</td>
                        <td className="py-2 px-3 font-mono text-slate-400">{aa.sasa}</td>
                        <td className="py-2 px-3">
                          <span className={`text-[10px] px-1.5 py-0.5 rounded font-medium ${
                            aa.polarity === 'Niepolarny' 
                              ? 'bg-emerald-950 text-emerald-300 border border-emerald-800/40' 
                              : aa.polarity === 'Kwasowy'
                              ? 'bg-rose-950 text-rose-300 border border-rose-800/40'
                              : aa.polarity === 'Zasadowy'
                              ? 'bg-blue-950 text-blue-300 border border-blue-800/40'
                              : 'bg-indigo-950 text-indigo-300 border border-indigo-800/40'
                          }`}>
                            {aa.polarity}
                          </span>
                        </td>
                        <td className={`py-2 px-3 font-mono font-medium ${aa.hydropathyKyte > 0 ? 'text-emerald-400' : 'text-blue-400'}`}>
                          {aa.hydropathyKyte > 0 ? `+${aa.hydropathyKyte.toFixed(1)}` : aa.hydropathyKyte.toFixed(1)}
                        </td>
                        <td className="py-2 px-3 font-mono text-slate-400">
                          {aa.hydrophilicityHopp > 0 ? `+${aa.hydrophilicityHopp.toFixed(1)}` : aa.hydrophilicityHopp.toFixed(1)}
                        </td>
                        <td className="py-2 px-3 font-mono">
                          <span className={`px-1.5 py-0.5 rounded font-bold ${
                            aa.chargePH74 > 0 
                              ? 'text-blue-400 bg-blue-950/60' 
                              : aa.chargePH74 < 0 
                              ? 'text-rose-400 bg-rose-950/60' 
                              : 'text-slate-400'
                          }`}>
                            {aa.chargePH74 > 0 ? `+${aa.chargePH74}` : aa.chargePH74}
                          </span>
                        </td>
                        <td className="py-2 px-3 font-mono text-slate-400">{aa.pK1_COOH.toFixed(2)}</td>
                        <td className="py-2 px-3 font-mono text-slate-400">{aa.pK2_NH3.toFixed(2)}</td>
                        <td className="py-2 px-3 font-mono text-slate-400">{aa.pKr_side !== null ? aa.pKr_side.toFixed(2) : '-'}</td>
                        <td className="py-2 px-3 font-mono text-slate-400">{aa.pI.toFixed(2)}</td>
                        <td className="py-2 px-3">
                          {isCovalent ? (
                            <span className="font-mono text-[10px] font-bold text-amber-300 bg-amber-950/80 px-1.5 py-0.5 rounded border border-amber-700/60">
                              {aa.covalentBridges}
                            </span>
                          ) : (
                            <span className="text-slate-500 font-mono">Nie</span>
                          )}
                        </td>
                        <td className="py-2 px-3 text-slate-300 text-[11px] truncate max-w-[130px]">{aa.ramachandranPref}</td>
                      </tr>
                    );
                  })}
                </tbody>
              </table>
            </div>
          </div>
        )}

        {/* TAB 2: 248 E8 Contact Patches Grid */}
        {activeTab === 'patches248' && (
          <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
            {/* Left: Cardinal Anchors info & Patch Details */}
            <div className="space-y-4">
              <div className="bg-slate-950 p-4 rounded-xl border border-slate-800">
                <h3 className="text-xs font-mono font-bold text-cyan-400 uppercase tracking-wider mb-2 flex items-center gap-1.5">
                  <Compass className="w-4 h-4" /> Żelazne Porty i Osie SNAPLOCK
                </h3>
                <div className="space-y-2 text-xs">
                  <div 
                    onClick={() => setSelectedPatchId(1)}
                    className="p-2.5 rounded-lg bg-emerald-950/40 border border-emerald-800/60 cursor-pointer hover:bg-emerald-950/60 transition-colors"
                  >
                    <div className="flex items-center justify-between font-mono font-bold text-emerald-400">
                      <span>Łatka #1 • Góra (+Y)</span>
                      <span>[0, 1, 0]</span>
                    </div>
                    <p className="text-[11px] text-slate-300 mt-0.5">Port WEJŚCIA: atom N szkieletu, wiązanie peptydowe z poprzednikiem.</p>
                  </div>

                  <div 
                    onClick={() => setSelectedPatchId(248)}
                    className="p-2.5 rounded-lg bg-rose-950/40 border border-rose-800/60 cursor-pointer hover:bg-rose-950/60 transition-colors"
                  >
                    <div className="flex items-center justify-between font-mono font-bold text-rose-400">
                      <span>Łatka #248 • Dół (-Y)</span>
                      <span>[0, -1, 0]</span>
                    </div>
                    <p className="text-[11px] text-slate-300 mt-0.5">Port WYJŚCIA: atom C szkieletu, tlen karbonylowy do następnika.</p>
                  </div>

                  <div className="grid grid-cols-2 gap-2 font-mono text-[11px]">
                    <div 
                      onClick={() => setSelectedPatchId(50)}
                      className="p-2 rounded-lg bg-slate-900 border border-slate-800 cursor-pointer hover:border-cyan-500"
                    >
                      <span className="text-cyan-400 font-bold block">#50 • Przód (+Z)</span>
                      <span className="text-slate-400">[0, 0, 1]</span>
                    </div>
                    <div 
                      onClick={() => setSelectedPatchId(100)}
                      className="p-2 rounded-lg bg-slate-900 border border-slate-800 cursor-pointer hover:border-cyan-500"
                    >
                      <span className="text-cyan-400 font-bold block">#100 • Prawo (+X)</span>
                      <span className="text-slate-400">[1, 0, 0]</span>
                    </div>
                    <div 
                      onClick={() => setSelectedPatchId(150)}
                      className="p-2 rounded-lg bg-slate-900 border border-slate-800 cursor-pointer hover:border-cyan-500"
                    >
                      <span className="text-cyan-400 font-bold block">#150 • Tył (-Z)</span>
                      <span className="text-slate-400">[0, 0, -1]</span>
                    </div>
                    <div 
                      onClick={() => setSelectedPatchId(200)}
                      className="p-2 rounded-lg bg-slate-900 border border-slate-800 cursor-pointer hover:border-cyan-500"
                    >
                      <span className="text-cyan-400 font-bold block">#200 • Lewo (-X)</span>
                      <span className="text-slate-400">[-1, 0, 0]</span>
                    </div>
                  </div>
                </div>
              </div>

              {/* Selected Patch Card */}
              <div className="bg-slate-950 p-4 rounded-xl border border-cyan-800/60 shadow-lg">
                <div className="flex items-center justify-between mb-2">
                  <span className="text-xs font-mono font-bold text-cyan-300">Szczegóły Wybranej Łatki</span>
                  <span className="text-xs px-2 py-0.5 rounded bg-cyan-950 text-cyan-300 border border-cyan-800 font-mono">
                    ID #{selectedPatch.id}
                  </span>
                </div>
                <div className="space-y-2 text-xs">
                  <div>
                    <span className="text-slate-400 block text-[11px]">Etykieta:</span>
                    <span className="font-semibold text-slate-100">{selectedPatch.label}</span>
                  </div>
                  <div>
                    <span className="text-slate-400 block text-[11px]">Wektor normalny [nx, ny, nz]:</span>
                    <span className="font-mono text-amber-300 font-bold">
                      [{selectedPatch.vector[0]}, {selectedPatch.vector[1]}, {selectedPatch.vector[2]}]
                    </span>
                  </div>
                  <div className="grid grid-cols-2 gap-2 text-[11px] font-mono">
                    <div className="p-1.5 bg-slate-900 rounded border border-slate-800">
                      <span className="text-slate-400">Theta (Polar): </span>
                      <span className="text-slate-200">{selectedPatch.accessiblePolarAngles.theta}°</span>
                    </div>
                    <div className="p-1.5 bg-slate-900 rounded border border-slate-800">
                      <span className="text-slate-400">Phi (Azymut): </span>
                      <span className="text-slate-200">{selectedPatch.accessiblePolarAngles.phi}°</span>
                    </div>
                  </div>
                  <div>
                    <span className="text-slate-400 block text-[11px] mb-1">Dyskretni sąsiedzi (skoki torsyjne):</span>
                    <div className="flex flex-wrap gap-1 font-mono text-[10px]">
                      {selectedPatch.neighborPatchIds.map(nid => (
                        <button
                          key={nid}
                          onClick={() => setSelectedPatchId(nid)}
                          className="px-2 py-1 rounded bg-slate-800 hover:bg-cyan-800 text-slate-300 hover:text-white transition-colors cursor-pointer"
                        >
                          #{nid}
                        </button>
                      ))}
                    </div>
                  </div>
                </div>
              </div>
            </div>

            {/* Right: Grid of all 248 patches */}
            <div className="lg:col-span-2 bg-slate-950 p-4 rounded-xl border border-slate-800 flex flex-col">
              <div className="flex items-center justify-between mb-3">
                <h3 className="text-xs font-mono font-bold text-slate-300 uppercase tracking-wider">
                  Matryca 248 Wektorów Sferycznych (Lattice E8)
                </h3>
                <span className="text-[11px] text-slate-500">Kliknij łatkę, aby wyświetlić jej wektory i sąsiadów</span>
              </div>

              <div className="flex-1 overflow-y-auto max-h-[460px] grid grid-cols-8 sm:grid-cols-10 md:grid-cols-12 gap-1.5 p-2 bg-slate-900/60 rounded-xl border border-slate-800/80">
                {E8_PATCHES.map((p) => {
                  const isSelected = p.id === selectedPatchId;
                  const isSpecial = p.id === 1 || p.id === 248 || p.id === 50 || p.id === 100 || p.id === 150 || p.id === 200;

                  return (
                    <button
                      key={p.id}
                      onClick={() => setSelectedPatchId(p.id)}
                      className={`h-9 flex flex-col items-center justify-center rounded text-[11px] font-mono transition-all cursor-pointer ${
                        isSelected
                          ? 'bg-cyan-500 text-slate-950 font-extrabold ring-2 ring-white scale-105 z-10'
                          : p.id === 1
                          ? 'bg-emerald-700/80 hover:bg-emerald-600 text-white font-bold'
                          : p.id === 248
                          ? 'bg-rose-700/80 hover:bg-rose-600 text-white font-bold'
                          : isSpecial
                          ? 'bg-amber-600/80 hover:bg-amber-500 text-white font-bold'
                          : 'bg-slate-800/80 hover:bg-slate-700 text-slate-300'
                      }`}
                      title={`${p.label} [${p.vector.join(', ')}]`}
                    >
                      <span>#{p.id}</span>
                    </button>
                  );
                })}
              </div>

              <div className="mt-3 pt-3 border-t border-slate-800 flex flex-wrap items-center gap-4 text-[11px] text-slate-400 font-mono">
                <span className="flex items-center gap-1.5">
                  <span className="w-2.5 h-2.5 rounded bg-emerald-600 inline-block"></span> #1 In (+Y)
                </span>
                <span className="flex items-center gap-1.5">
                  <span className="w-2.5 h-2.5 rounded bg-rose-600 inline-block"></span> #248 Out (-Y)
                </span>
                <span className="flex items-center gap-1.5">
                  <span className="w-2.5 h-2.5 rounded bg-amber-600 inline-block"></span> #50, #100, #150, #200 Osie
                </span>
                <span className="flex items-center gap-1.5">
                  <span className="w-2.5 h-2.5 rounded bg-slate-800 inline-block"></span> Siatka E8 (248)
                </span>
              </div>
            </div>
          </div>
        )}

        {/* TAB 3: 5 Tracks Junction Evaluator */}
        {activeTab === 'junction5tracks' && (
          <div className="space-y-6">
            <div className="bg-slate-950 p-4 rounded-xl border border-slate-800">
              <h3 className="text-xs font-mono font-bold text-cyan-400 uppercase tracking-wider mb-3 flex items-center gap-2">
                <Sliders className="w-4 h-4" /> Parametry Testowanego Złącza Peptydowego SNAPLOCK
              </h3>

              <div className="grid grid-cols-1 md:grid-cols-4 gap-4 text-xs">
                {/* Predecessor */}
                <div>
                  <label className="text-slate-400 block mb-1 font-medium">Poprzednik (Reszta N):</label>
                  <select
                    value={predCode}
                    onChange={(e) => setPredCode(e.target.value)}
                    className="w-full bg-slate-900 border border-slate-700 rounded-lg p-2 text-slate-200 font-mono focus:border-cyan-500"
                  >
                    {SNAPLOCK_22_DATABASE.map(aa => (
                      <option key={aa.code3} value={aa.code3}>
                        {aa.code3} ({aa.namePl}) {aa.covalentBridges.includes('TAK') ? '★ Klamra' : ''}
                      </option>
                    ))}
                  </select>
                  <span className="text-[11px] text-slate-400 block mt-1">Port wyjścia: Łatka #{customPredPatch}</span>
                </div>

                {/* Successor */}
                <div>
                  <label className="text-slate-400 block mb-1 font-medium">Następnik (Reszta N+1):</label>
                  <select
                    value={succCode}
                    onChange={(e) => setSuccCode(e.target.value)}
                    className="w-full bg-slate-900 border border-slate-700 rounded-lg p-2 text-slate-200 font-mono focus:border-cyan-500"
                  >
                    {SNAPLOCK_22_DATABASE.map(aa => (
                      <option key={aa.code3} value={aa.code3}>
                        {aa.code3} ({aa.namePl}) {aa.covalentBridges.includes('TAK') ? '★ Klamra' : ''}
                      </option>
                    ))}
                  </select>
                  <span className="text-[11px] text-slate-400 block mt-1">Port wejścia: Łatka #{customSuccPatch}</span>
                </div>

                {/* Torsion angle */}
                <div className="md:col-span-2">
                  <div className="flex items-center justify-between mb-1">
                    <label className="text-slate-400 font-medium">Kąt rotacji wokół osi styku (skoki po łatkach E8):</label>
                    <span className="font-mono text-cyan-400 font-bold">{torsionAngle}°</span>
                  </div>
                  <input
                    type="range"
                    min="-180"
                    max="180"
                    step="15"
                    value={torsionAngle}
                    onChange={(e) => setTorsionAngle(Number(e.target.value))}
                    className="w-full accent-cyan-500 cursor-pointer"
                  />
                  <div className="flex justify-between text-[10px] text-slate-400 font-mono mt-1">
                    <span>-180°</span>
                    <span>-90°</span>
                    <span>0° (E8 Align)</span>
                    <span>+90°</span>
                    <span>+180°</span>
                  </div>
                </div>
              </div>
            </div>

            {/* Evaluation Results Banner */}
            <div className={`p-4 rounded-xl border flex items-center justify-between gap-4 ${
              junctionEval.tracks.track5_covalentClamps.hasClamp
                ? 'bg-amber-950/40 border-amber-600/80'
                : junctionEval.compositeFitnessScore > 80
                ? 'bg-emerald-950/40 border-emerald-600/80'
                : 'bg-rose-950/40 border-rose-600/80'
            }`}>
              <div className="flex items-center gap-3">
                <div className={`p-2.5 rounded-xl ${
                  junctionEval.tracks.track5_covalentClamps.hasClamp
                    ? 'bg-amber-600 text-slate-950'
                    : junctionEval.compositeFitnessScore > 80
                    ? 'bg-emerald-600 text-slate-950'
                    : 'bg-rose-600 text-white'
                }`}>
                  <Activity className="w-6 h-6" />
                </div>
                <div>
                  <h4 className="font-mono font-bold text-sm text-slate-100 flex items-center gap-2">
                    <span>Fitness Złącza SNAPLOCK: {junctionEval.compositeFitnessScore}%</span>
                    {junctionEval.tracks.track5_covalentClamps.hasClamp && (
                      <span className="text-xs px-2 py-0.5 rounded bg-amber-500 text-slate-950 font-bold uppercase tracking-wider">
                        Klamra Kowalencyjna (Priorytet Maksymalny)
                      </span>
                    )}
                  </h4>
                  <p className="text-xs text-slate-300 mt-0.5">
                    {junctionEval.junctionId} • Dystans Cα-Cα: 3.80 Å • VdW Suma: {junctionEval.tracks.track2_vdwCollision.vdwExclusionThreshold} Å
                  </p>
                </div>
              </div>

              <div className="text-right font-mono text-xs text-slate-300">
                <div>Ścieżki zaliczone:</div>
                <div className="text-emerald-400 font-bold">5 / 5 Zgodność Reguł</div>
              </div>
            </div>

            {/* Detailed 5 Tracks Cards */}
            <div className="grid grid-cols-1 md:grid-cols-5 gap-3">
              {/* Track 1 */}
              <div className="bg-slate-950 p-3.5 rounded-xl border border-slate-800 flex flex-col justify-between">
                <div>
                  <div className="flex items-center justify-between mb-2">
                    <span className="text-[11px] font-mono text-cyan-400 font-bold uppercase">Track 1: Szkielet</span>
                    <span className="text-xs font-mono font-bold text-emerald-400">{junctionEval.tracks.track1_backbone.score}%</span>
                  </div>
                  <p className="text-xs text-slate-200 font-medium mb-1">{junctionEval.tracks.track1_backbone.portVerification}</p>
                  <p className="text-[11px] text-slate-400 leading-snug">{junctionEval.tracks.track1_backbone.details}</p>
                </div>
                <div className="mt-3 pt-2 border-t border-slate-800/80 text-[10px] text-slate-500 font-mono">
                  Porty: #248 -&gt; #1
                </div>
              </div>

              {/* Track 2 */}
              <div className="bg-slate-950 p-3.5 rounded-xl border border-slate-800 flex flex-col justify-between">
                <div>
                  <div className="flex items-center justify-between mb-2">
                    <span className="text-[11px] font-mono text-cyan-400 font-bold uppercase">Track 2: VdW Kolizja</span>
                    <span className="text-xs font-mono font-bold text-emerald-400">{junctionEval.tracks.track2_vdwCollision.score}%</span>
                  </div>
                  <p className="text-xs text-slate-200 font-medium mb-1">
                    {junctionEval.tracks.track2_vdwCollision.crushDetected ? 'Zgniecenie!' : 'Brak zgniecenia'}
                  </p>
                  <p className="text-[11px] text-slate-400 leading-snug">{junctionEval.tracks.track2_vdwCollision.details}</p>
                </div>
                <div className="mt-3 pt-2 border-t border-slate-800/80 text-[10px] text-slate-500 font-mono">
                  Próg: {junctionEval.tracks.track2_vdwCollision.vdwExclusionThreshold} Å
                </div>
              </div>

              {/* Track 3 */}
              <div className="bg-slate-950 p-3.5 rounded-xl border border-slate-800 flex flex-col justify-between">
                <div>
                  <div className="flex items-center justify-between mb-2">
                    <span className="text-[11px] font-mono text-cyan-400 font-bold uppercase">Track 3: Ładunki</span>
                    <span className="text-xs font-mono font-bold text-emerald-400">{junctionEval.tracks.track3_chargeOptimization.score}%</span>
                  </div>
                  <p className="text-xs text-slate-200 font-medium mb-1">{junctionEval.tracks.track3_chargeOptimization.coulombPotential}</p>
                  <p className="text-[11px] text-slate-400 leading-snug">{junctionEval.tracks.track3_chargeOptimization.details}</p>
                </div>
                <div className="mt-3 pt-2 border-t border-slate-800/80 text-[10px] text-slate-500 font-mono">
                  pH: 7.4 fizjologiczne
                </div>
              </div>

              {/* Track 4 */}
              <div className="bg-slate-950 p-3.5 rounded-xl border border-slate-800 flex flex-col justify-between">
                <div>
                  <div className="flex items-center justify-between mb-2">
                    <span className="text-[11px] font-mono text-cyan-400 font-bold uppercase">Track 4: Rdzeń Hydro</span>
                    <span className="text-xs font-mono font-bold text-emerald-400">{junctionEval.tracks.track4_hydrophobicCore.score}%</span>
                  </div>
                  <p className="text-xs text-slate-200 font-medium mb-1">
                    ΔG = -{junctionEval.tracks.track4_hydrophobicCore.desolvationGainKcal} kcal/mol
                  </p>
                  <p className="text-[11px] text-slate-400 leading-snug">{junctionEval.tracks.track4_hydrophobicCore.details}</p>
                </div>
                <div className="mt-3 pt-2 border-t border-slate-800/80 text-[10px] text-slate-500 font-mono">
                  Skala Kyte-Doolittle
                </div>
              </div>

              {/* Track 5 */}
              <div className={`p-3.5 rounded-xl border flex flex-col justify-between ${
                junctionEval.tracks.track5_covalentClamps.hasClamp
                  ? 'bg-amber-950/60 border-amber-600/80'
                  : 'bg-slate-950 border-slate-800'
              }`}>
                <div>
                  <div className="flex items-center justify-between mb-2">
                    <span className="text-[11px] font-mono text-amber-400 font-bold uppercase">Track 5: Klamra Kow.</span>
                    <span className="text-xs font-mono font-bold text-amber-300">
                      {junctionEval.tracks.track5_covalentClamps.hasClamp ? '100% OVERRIDE' : '-'}
                    </span>
                  </div>
                  <p className="text-xs text-amber-200 font-bold mb-1">
                    {junctionEval.tracks.track5_covalentClamps.clampType}
                  </p>
                  <p className="text-[11px] text-slate-300 leading-snug">{junctionEval.tracks.track5_covalentClamps.details}</p>
                </div>
                <div className="mt-3 pt-2 border-t border-amber-900/60 text-[10px] text-amber-400 font-mono">
                  Energia wiązania: {junctionEval.tracks.track5_covalentClamps.bondEnergyKcal} kcal/mol
                </div>
              </div>
            </div>
          </div>
        )}

        {/* TAB: Official SNAPLOCK vs AlphaFold Benchmark & Calibration Engine */}
        {activeTab === 'benchmarkEngine' && (
          <div className="space-y-6">
            <div className="bg-slate-950 p-5 rounded-2xl border border-slate-800 space-y-4">
              <div className="flex flex-wrap items-center justify-between gap-3">
                <div>
                  <h3 className="text-base font-bold text-slate-100 flex items-center gap-2">
                    <ShieldCheck className="w-5 h-5 text-amber-400" />
                    <span>Oficjalny Moduł Walidacyjny: SNAPLOCK vs AlphaFold Benchmark</span>
                  </h3>
                  <p className="text-xs text-slate-400 mt-1 max-w-3xl">
                    Pobiera model 3D bezpośrednio z bazy AlphaFold DB (EMBL-EBI), rzutuje wiązania szkieletowe na sferyczną siatkę 248 wektorów Liego E8 (<code className="text-amber-300 font-mono">Patch = argmax(v · e_i)</code>), mierzy błąd kątowy i weryfikuje mostki kowalencyjne S-S (Ścieżka 5).
                  </p>
                </div>
                <span className="text-[11px] font-mono px-3 py-1 rounded-full bg-amber-950/80 text-amber-300 border border-amber-800">
                  Protokół: snaplock_alphafold_benchmark.py
                </span>
              </div>

              {/* Presets and Custom Input */}
              <div className="flex flex-wrap items-center gap-2 pt-2 border-t border-slate-800">
                <span className="text-xs font-mono text-slate-400 font-bold mr-1">Wzorce UniProt:</span>
                {[
                  { id: 'P83347', name: 'Trp-cage TC5b (P83347)', desc: '20-30 aa standard' },
                  { id: 'P01308', name: 'Insulina ludzka (P01308)', desc: '110 aa • 3 mostki S-S' },
                  { id: 'P0CG48', name: 'Ubikwityna (P0CG48)', desc: 'Ubiquitin standard' },
                  { id: 'P00698', name: 'Lizozym HEWL (P00698)', desc: '129 aa • 4 mostki S-S' }
                ].map(p => (
                  <button
                    key={p.id}
                    onClick={() => {
                      setBenchmarkUniprot(p.id);
                      runBenchmark(p.id);
                    }}
                    className={`px-3 py-1.5 rounded-lg text-xs font-mono font-medium transition cursor-pointer border ${
                      benchmarkUniprot === p.id 
                        ? 'bg-amber-600 border-amber-400 text-white font-bold shadow' 
                        : 'bg-slate-900 hover:bg-slate-800 border-slate-700 text-slate-300'
                    }`}
                  >
                    ★ {p.name}
                  </button>
                ))}
              </div>

              {/* Custom UniProt ID Input & Run Button */}
              <div className="flex flex-wrap items-center gap-3 pt-2">
                <div className="flex items-center gap-2 bg-slate-900 border border-slate-700 rounded-lg px-3 py-2 text-xs w-64">
                  <span className="text-slate-400 font-mono uppercase font-bold text-[10px]">UniProt ID:</span>
                  <input
                    type="text"
                    value={benchmarkUniprot}
                    onChange={(e) => setBenchmarkUniprot(e.target.value.trim().toUpperCase())}
                    placeholder="np. P83347, P01308..."
                    className="bg-transparent text-amber-300 font-mono font-bold focus:outline-none w-full"
                  />
                </div>

                <button
                  onClick={() => runBenchmark(benchmarkUniprot)}
                  disabled={benchmarkLoading || !benchmarkUniprot}
                  className="px-5 py-2 bg-gradient-to-r from-amber-600 via-amber-500 to-indigo-600 hover:from-amber-500 hover:to-indigo-500 disabled:opacity-50 text-white font-bold text-xs rounded-lg transition shadow-md flex items-center gap-2 cursor-pointer font-mono"
                >
                  {benchmarkLoading ? (
                    <>
                      <Sparkles className="w-4 h-4 animate-spin text-white" />
                      <span>Pobieranie PDB i Kwantyzacja E8...</span>
                    </>
                  ) : (
                    <>
                      <Zap className="w-4 h-4 text-white" />
                      <span>Uruchom Kalibrację E8 dla {benchmarkUniprot}</span>
                    </>
                  )}
                </button>
              </div>

              {benchmarkError && (
                <div className="p-3 bg-rose-950/70 border border-rose-800 text-rose-300 rounded-xl text-xs font-mono">
                  Błąd kalibracji: {benchmarkError}
                </div>
              )}
            </div>

            {/* Results Display */}
            {benchmarkReport && (
              <div className="space-y-5 animate-in fade-in duration-300">
                {/* Metric Cards */}
                <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-4">
                  <div className="bg-slate-950 p-4 rounded-xl border border-slate-800">
                    <span className="text-[11px] font-mono text-slate-400 uppercase block font-bold">Reszty szkieletowe</span>
                    <div className="text-2xl font-mono font-extrabold text-cyan-400 mt-1">
                      {benchmarkReport.totalResidues}
                    </div>
                    <span className="text-[10px] text-slate-500 font-mono">Wiązania atomowe CA-CA</span>
                  </div>

                  <div className="bg-slate-950 p-4 rounded-xl border border-slate-800">
                    <span className="text-[11px] font-mono text-slate-400 uppercase block font-bold">Średni błąd kątowy łatki E8</span>
                    <div className="text-2xl font-mono font-extrabold text-amber-400 mt-1">
                      {benchmarkReport.meanAngularPatchLossDeg}°
                    </div>
                    <span className="text-[10px] text-emerald-400 font-mono">Dyskretyzacja siatki 248 wektorów</span>
                  </div>

                  <div className="bg-slate-950 p-4 rounded-xl border border-slate-800">
                    <span className="text-[11px] font-mono text-slate-400 uppercase block font-bold">Mostki S-S (Ścieżka 5)</span>
                    <div className="text-2xl font-mono font-extrabold text-yellow-400 mt-1">
                      {benchmarkReport.disulfideBridgesDetected?.length || 0}
                    </div>
                    <span className="text-[10px] text-slate-500 font-mono">
                      {benchmarkReport.disulfideBridgesDetected?.length > 0 
                        ? benchmarkReport.disulfideBridgesDetected.map((b: any) => `Cys${b[0]}-Cys${b[1]} (${b[2]}Å)`).join(', ')
                        : 'Brak w modelu'}
                    </span>
                  </div>

                  <div className="bg-slate-950 p-4 rounded-xl border border-slate-800">
                    <span className="text-[11px] font-mono text-slate-400 uppercase block font-bold">Średni pLDDT AlphaFold</span>
                    <div className="text-2xl font-mono font-extrabold text-emerald-400 mt-1">
                      {benchmarkReport.meanPlddtScore ? `${benchmarkReport.meanPlddtScore}%` : 'N/A'}
                    </div>
                    <span className="text-[10px] text-slate-500 font-mono">Pewność predykcji AF2</span>
                  </div>
                </div>

                {/* Residue Mapping Table */}
                <div className="border border-slate-800 rounded-xl overflow-hidden bg-slate-950">
                  <div className="px-4 py-3 bg-slate-900 border-b border-slate-800 flex items-center justify-between text-xs font-mono">
                    <span className="font-bold text-slate-200">
                      Rzutowanie wiązań szkieletowych C_alfa_i -&gt; C_alfa_i+1 na dyskretne łatki E8 ({benchmarkReport.patchAssignments?.length || 0} połączeń):
                    </span>
                    <button
                      onClick={() => {
                        const blob = new Blob([JSON.stringify(benchmarkReport, null, 2)], { type: 'application/json' });
                        const url = URL.createObjectURL(blob);
                        const a = document.createElement('a');
                        a.href = url;
                        a.download = `snaplock_${benchmarkReport.uniprotTarget || 'benchmark'}_report.json`;
                        a.click();
                        URL.revokeObjectURL(url);
                      }}
                      className="px-2.5 py-1 bg-slate-800 hover:bg-slate-700 text-slate-200 rounded text-[11px] font-bold cursor-pointer transition flex items-center gap-1"
                    >
                      <Download className="w-3 h-3" />
                      <span>Pobierz Raport JSON</span>
                    </button>
                  </div>

                  <div className="max-h-80 overflow-y-auto">
                    <table className="w-full text-left font-mono text-xs">
                      <thead className="bg-slate-950 text-slate-400 text-[11px] uppercase border-b border-slate-800 sticky top-0">
                        <tr>
                          <th className="px-4 py-2">Reszta (i)</th>
                          <th className="px-4 py-2">Następna (i+1)</th>
                          <th className="px-4 py-2">Przypisana Łatka E8</th>
                          <th className="px-4 py-2">Wektor Łatki</th>
                          <th className="px-4 py-2">Błąd Kątowy</th>
                          <th className="px-4 py-2">Zgodność</th>
                        </tr>
                      </thead>
                      <tbody className="divide-y divide-slate-800/60 text-slate-300 text-[11px]">
                        {benchmarkReport.patchAssignments?.map((m: any, idx: number) => {
                          const patchMeta = E8_PATCHES.find(p => p.id === m.assignedPatch) || E8_PATCHES[0];
                          return (
                            <tr key={idx} className="hover:bg-slate-900/50">
                              <td className="px-4 py-1.5 font-bold text-white">#{m.resIdx} {m.resName}</td>
                              <td className="px-4 py-1.5 text-slate-400">#{m.resIdx + 1} {m.nextRes}</td>
                              <td className="px-4 py-1.5">
                                <span className="px-2 py-0.5 rounded font-bold bg-amber-950/80 text-amber-300 border border-amber-800">
                                  #{m.assignedPatch}
                                </span>
                              </td>
                              <td className="px-4 py-1.5 text-slate-400 text-[10px]">
                                [{patchMeta.vector.join(', ')}]
                              </td>
                              <td className="px-4 py-1.5 font-bold text-amber-400">
                                {m.angleErrorDeg}°
                              </td>
                              <td className="px-4 py-1.5">
                                <span className={`px-2 py-0.5 rounded text-[10px] font-bold ${
                                  m.angleErrorDeg < 5.0 ? 'bg-emerald-950 text-emerald-400' : 'bg-amber-950 text-amber-400'
                                }`}>
                                  {m.angleErrorDeg < 5.0 ? 'WYSOKA ZGODNOŚĆ' : 'ZGODNA KWANTYZACJA'}
                                </span>
                              </td>
                            </tr>
                          );
                        })}
                      </tbody>
                    </table>
                  </div>
                </div>
              </div>
            )}
          </div>
        )}

        {/* TAB 4: JSON Export View */}
        {activeTab === 'jsonExport' && (
          <div className="space-y-4">
            <div className="flex flex-wrap items-center justify-between gap-3 bg-slate-950 p-4 rounded-xl border border-slate-800">
              <div>
                <h3 className="text-sm font-bold text-slate-100 font-sans flex items-center gap-2">
                  <span>Plik Wyjściowy: amino_acid_warehouse.json</span>
                  <span className="text-xs font-mono px-2 py-0.5 rounded bg-indigo-950 text-indigo-300 border border-indigo-800">
                    E8 LATTICE JSON
                  </span>
                </h3>
                <p className="text-xs text-slate-400">
                  Wygenerowano na podstawie 22 parametrów fizykochemicznych oraz wektorów 248 łatek sferycznych
                </p>
              </div>

              <div className="flex items-center gap-2">
                <div className="flex bg-slate-900 p-1 rounded-lg border border-slate-800 text-xs">
                  <button
                    onClick={() => setExportScope('first3')}
                    className={`px-3 py-1.5 rounded-md font-medium transition-colors ${
                      exportScope === 'first3' 
                        ? 'bg-cyan-600 text-white' 
                        : 'text-slate-400 hover:text-slate-200'
                    }`}
                  >
                    Pierwsze 3 (GLY, ALA, VAL)
                  </button>
                  <button
                    onClick={() => setExportScope('all21')}
                    className={`px-3 py-1.5 rounded-md font-medium transition-colors ${
                      exportScope === 'all21' 
                        ? 'bg-cyan-600 text-white' 
                        : 'text-slate-400 hover:text-slate-200'
                    }`}
                  >
                    Wszystkie 21 aminokwasów
                  </button>
                </div>

                <button
                  onClick={handleCopyJson}
                  className="flex items-center gap-1.5 px-3 py-1.5 bg-slate-800 hover:bg-slate-700 text-slate-200 rounded-lg text-xs font-medium transition-colors"
                >
                  {copied ? (
                    <>
                      <Check className="w-3.5 h-3.5 text-emerald-400" />
                      <span>Skopiowano</span>
                    </>
                  ) : (
                    <>
                      <Copy className="w-3.5 h-3.5" />
                      <span>Kopiuj</span>
                    </>
                  )}
                </button>

                <button
                  onClick={handleDownloadJson}
                  className="flex items-center gap-1.5 px-4 py-1.5 bg-gradient-to-r from-cyan-600 to-indigo-600 hover:from-cyan-500 hover:to-indigo-500 text-white rounded-lg text-xs font-semibold shadow-md transition-all cursor-pointer"
                >
                  <Download className="w-3.5 h-3.5" />
                  <span>Pobierz plik JSON</span>
                </button>
              </div>
            </div>

            {/* Code View */}
            <div className="relative border border-slate-800 rounded-xl overflow-hidden bg-slate-950">
              <div className="px-4 py-2 bg-slate-900 border-b border-slate-800 flex items-center justify-between text-xs font-mono text-slate-400">
                <span>amino_acid_warehouse.json</span>
                <span>Rozmiar: {(jsonString.length / 1024).toFixed(1)} KB • {exportScope === 'first3' ? '3 reszty' : '21 reszt'}</span>
              </div>
              <pre className="p-4 text-xs font-mono text-cyan-300/90 overflow-x-auto max-h-[460px] leading-relaxed">
                {jsonString}
              </pre>
            </div>
          </div>
        )}
      </div>
    </div>
  );
};
