import React, { useState } from 'react';
import { 
  Atom, 
  Search, 
  Plus, 
  ArrowRight, 
  Sliders, 
  ShieldCheck, 
  Flame, 
  Activity, 
  ExternalLink,
  Layers,
  Sparkles,
  Info
} from 'lucide-react';
import { SNAPLOCK_22_DATABASE, AminoAcidPhysicochemicalData } from '../snaplockEngine';

interface AminoAcidCatalogViewProps {
  onSelectAminoAcidForConstructor?: (code3: string) => void;
  onSwitchToConstructor?: () => void;
}

export const AminoAcidCatalogView: React.FC<AminoAcidCatalogViewProps> = ({
  onSelectAminoAcidForConstructor,
  onSwitchToConstructor
}) => {
  const [selectedCode, setSelectedCode] = useState<string>('ALA');
  const [filterType, setFilterType] = useState<'ALL' | 'HYDROPHOBIC' | 'POLAR' | 'CHARGED' | 'SPECIAL'>('ALL');
  const [search, setSearch] = useState<string>('');
  const [includeSpecial, setIncludeSpecial] = useState<boolean>(true); // SEC (21) & PYL (22)

  const standard20Codes = [
    'GLY', 'ALA', 'VAL', 'LEU', 'ILE', 
    'PRO', 'PHE', 'TYR', 'TRP', 'SER', 
    'THR', 'CYS', 'MET', 'ASN', 'GLN', 
    'ASP', 'GLU', 'LYS', 'ARG', 'HIS'
  ];

  const filteredAAs = SNAPLOCK_22_DATABASE.filter(aa => {
    const isSpecial = aa.code3 === 'SEC' || aa.code3 === 'PYL';
    if (!includeSpecial && isSpecial) return false;

    const matchesSearch = 
      aa.code3.toLowerCase().includes(search.toLowerCase()) ||
      aa.namePl.toLowerCase().includes(search.toLowerCase()) ||
      aa.code1.toLowerCase().includes(search.toLowerCase());
    if (!matchesSearch) return false;

    if (filterType === 'ALL') return true;
    if (filterType === 'HYDROPHOBIC') return aa.hydropathyKyte > 0;
    if (filterType === 'POLAR') return aa.polarity.toLowerCase().includes('polarny');
    if (filterType === 'CHARGED') return aa.chargePH74 !== 0;
    if (filterType === 'SPECIAL') return isSpecial || aa.covalentBridges.includes('TAK');
    return true;
  });

  const activeAA: AminoAcidPhysicochemicalData = 
    SNAPLOCK_22_DATABASE.find(a => a.code3 === selectedCode) || SNAPLOCK_22_DATABASE[0];

  return (
    <div className="flex-1 w-full h-full min-h-0 bg-[#05060b] text-zinc-100 flex flex-col font-sans overflow-hidden">
      
      {/* Header bar */}
      <div className="bg-[#0b0c16] border-b border-zinc-800 px-5 py-3 flex flex-wrap items-center justify-between gap-3 shrink-0">
        <div className="flex items-center gap-3">
          <div className="w-8 h-8 rounded-lg bg-cyan-950/70 border border-cyan-500/30 flex items-center justify-center">
            <Atom className="w-5 h-5 text-cyan-400" />
          </div>
          <div>
            <h1 className="text-sm font-bold text-white tracking-wide flex items-center gap-2">
              Katalog 20 Standardowych Aminokwasów
              <span className="text-[10px] px-2 py-0.5 rounded bg-zinc-800 text-cyan-300 font-mono border border-zinc-700">
                Jednostki: Ångström (Å) &amp; kcal/mol
              </span>
            </h1>
            <p className="text-[11px] text-zinc-400">
              Czysta baza fizykochemiczna bez siatki E8. 20 aminokwasów standardowych + opcjonalne SEC (21) i PYL (22).
            </p>
          </div>
        </div>

        <div className="flex items-center gap-2">
          <label className="flex items-center gap-2 text-xs text-zinc-300 cursor-pointer bg-black/60 px-2.5 py-1.5 rounded border border-zinc-800">
            <input 
              type="checkbox" 
              checked={includeSpecial}
              onChange={(e) => setIncludeSpecial(e.target.checked)}
              className="accent-cyan-500 rounded"
            />
            <span>Pokaż SEC (U) &amp; PYL (O)</span>
          </label>

          {onSwitchToConstructor && (
            <button
              onClick={onSwitchToConstructor}
              className="px-3 py-1.5 bg-indigo-600 hover:bg-indigo-500 text-white rounded text-xs font-bold flex items-center gap-1.5 transition cursor-pointer"
            >
              <span>Przejdź do Konstruktora 3D</span>
              <ArrowRight className="w-3.5 h-3.5" />
            </button>
          )}
        </div>
      </div>

      {/* Main split view: Catalog grid (left) + Residue physicochemical sheet (right) */}
      <div className="flex-1 w-full min-h-0 flex flex-col lg:flex-row overflow-hidden">
        
        {/* Left: AA Grid & Filter controls */}
        <div className="w-full lg:w-7/12 flex flex-col border-b lg:border-b-0 lg:border-r border-zinc-800 min-h-0 bg-[#070810]">
          
          {/* Filter Bar */}
          <div className="p-3 border-b border-zinc-800/80 flex flex-wrap items-center justify-between gap-2 shrink-0 bg-[#090b14]">
            <div className="flex items-center gap-1 flex-wrap">
              {(['ALL', 'HYDROPHOBIC', 'POLAR', 'CHARGED', 'SPECIAL'] as const).map(f => (
                <button
                  key={f}
                  onClick={() => setFilterType(f)}
                  className={`px-2.5 py-1 rounded text-[11px] font-mono font-bold transition cursor-pointer ${
                    filterType === f 
                      ? 'bg-cyan-500 text-black shadow-md' 
                      : 'bg-zinc-900 hover:bg-zinc-800 text-zinc-400'
                  }`}
                >
                  {f === 'ALL' && 'Wszystkie'}
                  {f === 'HYDROPHOBIC' && 'Hydrofobowe'}
                  {f === 'POLAR' && 'Polarne'}
                  {f === 'CHARGED' && 'Naładowane'}
                  {f === 'SPECIAL' && 'Specjalne/S-S'}
                </button>
              ))}
            </div>

            <div className="relative w-44">
              <Search className="w-3.5 h-3.5 absolute left-2.5 top-2 text-zinc-500" />
              <input
                type="text"
                value={search}
                onChange={(e) => setSearch(e.target.value)}
                placeholder="Szukaj (np. CYS, C, Alanina)..."
                className="w-full pl-8 pr-2 py-1 text-xs bg-black border border-zinc-800 rounded text-zinc-200 focus:outline-none focus:border-cyan-500"
              />
            </div>
          </div>

          {/* Cards Grid */}
          <div className="flex-1 overflow-y-auto p-4 grid grid-cols-2 sm:grid-cols-3 md:grid-cols-4 gap-2.5 scrollbar-thin">
            {filteredAAs.map(aa => {
              const isSelected = aa.code3 === selectedCode;
              const isSpecial = aa.code3 === 'SEC' || aa.code3 === 'PYL';
              const hydropathyColor = 
                aa.hydropathyKyte > 2.0 ? 'text-amber-400 border-amber-500/30' :
                aa.hydropathyKyte < -2.0 ? 'text-cyan-400 border-cyan-500/30' :
                'text-zinc-300 border-zinc-700';

              return (
                <div
                  key={aa.code3}
                  onClick={() => setSelectedCode(aa.code3)}
                  className={`p-3 rounded-xl border text-left cursor-pointer transition flex flex-col justify-between relative group ${
                    isSelected 
                      ? 'bg-cyan-950/40 border-cyan-500 shadow-lg shadow-cyan-500/10 ring-1 ring-cyan-500/40' 
                      : 'bg-[#0d0f1a] hover:bg-[#131626] border-zinc-800/90'
                  }`}
                >
                  {isSpecial && (
                    <span className="absolute top-2 right-2 text-[9px] font-mono px-1.5 py-0.2 rounded bg-purple-900/60 text-purple-300 border border-purple-700/50">
                      #{aa.nr}
                    </span>
                  )}
                  <div>
                    <div className="flex items-baseline justify-between mb-1">
                      <span className="font-mono text-base font-extrabold text-white group-hover:text-cyan-300 transition">
                        {aa.code3}
                      </span>
                      <span className="font-mono text-xs px-1.5 py-0.5 rounded bg-black/60 text-zinc-400 border border-zinc-800">
                        {aa.code1}
                      </span>
                    </div>
                    <div className="text-xs font-medium text-zinc-300 truncate">
                      {aa.namePl}
                    </div>
                    <div className="text-[10px] font-mono text-zinc-500 truncate mt-0.5">
                      R: {aa.sidechainFormula}
                    </div>
                  </div>

                  <div className="mt-3 pt-2 border-t border-zinc-800/70 grid grid-cols-2 gap-1 text-[10px] font-mono">
                    <div>
                      <span className="text-zinc-500">KD: </span>
                      <span className={aa.hydropathyKyte >= 0 ? 'text-amber-400' : 'text-blue-400'}>
                        {aa.hydropathyKyte > 0 ? `+${aa.hydropathyKyte}` : aa.hydropathyKyte}
                      </span>
                    </div>
                    <div>
                      <span className="text-zinc-500">VdW: </span>
                      <span className="text-zinc-300">{aa.vdwVolume} Å³</span>
                    </div>
                  </div>
                </div>
              );
            })}
          </div>

          <div className="p-2.5 bg-[#090b14] border-t border-zinc-800/80 text-[11px] text-zinc-400 flex items-center justify-between font-mono">
            <span>Widoczne aminokwasy: {filteredAAs.length}</span>
            <span className="text-zinc-500">1 Å = 0.1 nm = 100 pm</span>
          </div>
        </div>

        {/* Right: Detailed Physicochemical Profile & Physics Metrics in Å & kcal/mol */}
        <div className="w-full lg:w-5/12 flex flex-col bg-[#05060b] min-h-0 overflow-y-auto p-5 scrollbar-thin">
          <div className="bg-[#0b0c16] border border-zinc-800 rounded-2xl p-5 mb-4 shadow-xl">
            <div className="flex items-start justify-between">
              <div>
                <div className="flex items-center gap-3">
                  <span className="text-3xl font-extrabold font-mono text-white tracking-wider">
                    {activeAA.code3}
                  </span>
                  <span className="text-xl font-mono px-2.5 py-1 rounded bg-cyan-950/80 text-cyan-300 border border-cyan-500/40">
                    {activeAA.code1}
                  </span>
                  <span className="text-xs font-mono text-zinc-500">
                    Pozycja w tabeli: #{activeAA.nr}
                  </span>
                </div>
                <h2 className="text-base font-semibold text-zinc-200 mt-1">
                  {activeAA.namePl}
                </h2>
              </div>

              {onSelectAminoAcidForConstructor && (
                <button
                  onClick={() => onSelectAminoAcidForConstructor(activeAA.code3)}
                  className="px-3 py-1.5 bg-cyan-600 hover:bg-cyan-500 text-black font-bold rounded-lg text-xs flex items-center gap-1.5 transition cursor-pointer"
                  title="Wstaw ten aminokwas na koniec sekwencji w Konstruktorze 3D"
                >
                  <Plus className="w-3.5 h-3.5" />
                  <span>Dodaj do Konstruktora</span>
                </button>
              )}
            </div>

            {/* Formula & Core properties */}
            <div className="mt-4 p-3 rounded-xl bg-black/60 border border-zinc-800 grid grid-cols-2 gap-3 text-xs font-mono">
              <div>
                <span className="text-zinc-500 block text-[10px]">ŁAŃCUCH BOCZNY (GRUPA R):</span>
                <span className="text-cyan-300 font-bold text-sm">{activeAA.sidechainFormula}</span>
              </div>
              <div>
                <span className="text-zinc-500 block text-[10px]">MASA CZĄSTECZKOWA:</span>
                <span className="text-zinc-200 font-bold text-sm">{activeAA.mass.toFixed(2)} g/mol</span>
              </div>
            </div>
          </div>

          {/* Physical Parameters Table in real units: Å and kcal/mol */}
          <div className="space-y-3">
            <h3 className="text-xs font-bold text-zinc-400 uppercase tracking-wider font-mono flex items-center gap-2">
              <Sliders className="w-3.5 h-3.5 text-cyan-400" />
              <span>Parametry Fizykochemiczne w Prawdziwych Jednostkach</span>
            </h3>

            <div className="grid grid-cols-1 sm:grid-cols-2 gap-2 text-xs font-mono">
              
              <div className="p-3 rounded-xl bg-[#090b14] border border-zinc-800/80">
                <span className="text-zinc-500 text-[10px] block">OBJĘTOŚĆ VAN DER WAALSA:</span>
                <span className="text-base font-bold text-amber-300">{activeAA.vdwVolume} Å³</span>
                <span className="text-[10px] text-zinc-500 block mt-0.5">
                  R_vdw ≈ {(Math.cbrt((3 * activeAA.vdwVolume) / (4 * Math.PI))).toFixed(2)} Å
                </span>
              </div>

              <div className="p-3 rounded-xl bg-[#090b14] border border-zinc-800/80">
                <span className="text-zinc-500 text-[10px] block">DOSTĘPNA POWIERZCHNIA (SASA):</span>
                <span className="text-base font-bold text-emerald-300">{activeAA.sasa} Å²</span>
                <span className="text-[10px] text-zinc-500 block mt-0.5">Solvent-Accessible Surface Area</span>
              </div>

              <div className="p-3 rounded-xl bg-[#090b14] border border-zinc-800/80">
                <span className="text-zinc-500 text-[10px] block">HYDROPATIA (KYTE-DOOLITTLE):</span>
                <span className={`text-base font-bold ${activeAA.hydropathyKyte >= 0 ? 'text-amber-400' : 'text-blue-400'}`}>
                  {activeAA.hydropathyKyte > 0 ? `+${activeAA.hydropathyKyte}` : activeAA.hydropathyKyte}
                </span>
                <span className="text-[10px] text-zinc-500 block mt-0.5">
                  {activeAA.hydropathyKyte > 0 ? 'Rdzeń hydrofobowy' : 'Powierzchnia hydrofilowa'}
                </span>
              </div>

              <div className="p-3 rounded-xl bg-[#090b14] border border-zinc-800/80">
                <span className="text-zinc-500 text-[10px] block">HYDROFILOWOŚĆ (HOPP-WOODS):</span>
                <span className="text-base font-bold text-blue-300">
                  {activeAA.hydrophilicityHopp > 0 ? `+${activeAA.hydrophilicityHopp}` : activeAA.hydrophilicityHopp}
                </span>
                <span className="text-[10px] text-zinc-500 block mt-0.5">Eksponowanie na rozpuszczalnik</span>
              </div>

              <div className="p-3 rounded-xl bg-[#090b14] border border-zinc-800/80">
                <span className="text-zinc-500 text-[10px] block">ŁADUNEK W pH 7.4:</span>
                <span className={`text-base font-bold ${
                  activeAA.chargePH74 > 0 ? 'text-cyan-400' : activeAA.chargePH74 < 0 ? 'text-rose-400' : 'text-zinc-300'
                }`}>
                  {activeAA.chargePH74 > 0 ? `+${activeAA.chargePH74}` : activeAA.chargePH74} e
                </span>
                <span className="text-[10px] text-zinc-500 block mt-0.5">
                  pI = {activeAA.pI.toFixed(2)}
                </span>
              </div>

              <div className="p-3 rounded-xl bg-[#090b14] border border-zinc-800/80">
                <span className="text-zinc-500 text-[10px] block">STAŁE DYSOCJACJI (pKa):</span>
                <span className="text-xs text-zinc-300 block">
                  pK₁ (COOH): <strong className="text-white">{activeAA.pK1_COOH.toFixed(2)}</strong>
                </span>
                <span className="text-xs text-zinc-300 block">
                  pK₂ (NH₃⁺): <strong className="text-white">{activeAA.pK2_NH3.toFixed(2)}</strong>
                </span>
                {activeAA.pKr_side !== null && (
                  <span className="text-xs text-amber-300 block">
                    pKᵣ (grupa R): <strong>{activeAA.pKr_side.toFixed(2)}</strong>
                  </span>
                )}
              </div>

            </div>

            {/* Structural Preferences (Ramachandran, Covalent Bridges) */}
            <div className="p-4 rounded-xl bg-[#090b14] border border-zinc-800/80 space-y-2 text-xs font-mono">
              <div className="flex items-center justify-between">
                <span className="text-zinc-500">PREFERENCJA RAMACHANDRANA:</span>
                <span className="text-cyan-300 font-bold">{activeAA.ramachandranPref}</span>
              </div>
              <div className="flex items-center justify-between">
                <span className="text-zinc-500">ELASTYCZNOŚĆ ŁAŃCUCHA:</span>
                <span className="text-zinc-200">{activeAA.chainFlexibility}</span>
              </div>
              <div className="flex items-center justify-between">
                <span className="text-zinc-500">PREFERENCJE ROTAMEROWE:</span>
                <span className="text-zinc-200">{activeAA.rotamerPreferences}</span>
              </div>
              <div className="flex items-center justify-between">
                <span className="text-zinc-500">MOSTKI KOWALENCYJNE (S-S):</span>
                <span className={`font-bold ${activeAA.covalentBridges.includes('TAK') ? 'text-amber-400' : 'text-zinc-400'}`}>
                  {activeAA.covalentBridges}
                </span>
              </div>
            </div>

            <div className="p-3 rounded-xl bg-cyan-950/20 border border-cyan-800/30 text-[11px] text-zinc-300 flex items-start gap-2.5">
              <Info className="w-4 h-4 text-cyan-400 shrink-0 mt-0.5" />
              <span>
                Wszystkie parametry w tym katalogu są bezpośrednio sprzężone z silnikiem fizyki w ångströmach (Å) i energii w kcal/mol. Wyeliminowano wszelkie dyskretne zniekształcenia siatki E8.
              </span>
            </div>

          </div>
        </div>
      </div>
    </div>
  );
};
