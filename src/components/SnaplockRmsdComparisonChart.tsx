import React, { useMemo, useState } from 'react';
import {
  ResponsiveContainer,
  ComposedChart,
  Area,
  Bar,
  Line,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  Legend,
  ReferenceLine,
  Brush,
  Cell
} from 'recharts';
import {
  Activity,
  AlertTriangle,
  Atom,
  BarChart3,
  CheckCircle2,
  GitCompare,
  Layers,
  ShieldCheck
} from 'lucide-react';
import {
  AllAtomStructureResult,
  AlphaFoldRmsdComparisonResult,
  PerResidueRmsdEntry
} from '../snaplockAllAtomBuilder';
import { SNAPLOCK_22_DATABASE } from '../snaplockEngine';

export type RmsdChartViewMode = 'multilayer' | 'sidechain_vs_backbone' | 'threshold_bars';

export interface SnaplockRmsdComparisonChartProps {
  rmsdComparison: AlphaFoldRmsdComparisonResult;
  allAtomStructure?: AllAtomStructureResult | null;
  secondaryTypes?: ('HELIX' | 'SHEET' | 'TURN' | 'COIL')[];
  disulfideResidueSeqs?: number[]; // 1-based model residue indices involved in S-S / Se-Se clamps
  activeResidueIdx?: number;       // 0-based selected residue index in SNAPLOCK model
  onSelectResidue?: (zeroBasedIdx: number) => void;
  compact?: boolean;
}

interface ChartDatum extends PerResidueRmsdEntry {
  xLabel: string;
  aaNamePl: string;
  rotamerLabel: string;
  secondary: 'HELIX' | 'SHEET' | 'TURN' | 'COIL';
  isDisulfideClamp: boolean;
  isSelected: boolean;
}

function getPlddtTierLabel(plddt: number): { label: string; colorClass: string } {
  if (plddt >= 90) return { label: 'Bardzo wysoka (>90)', colorClass: 'text-sky-300' };
  if (plddt >= 70) return { label: 'Wysoka (70–90)', colorClass: 'text-cyan-300' };
  if (plddt >= 50) return { label: 'Niska (50–70)', colorClass: 'text-amber-300' };
  return { label: 'Nieuporządkowana (<50)', colorClass: 'text-rose-400' };
}

export const SnaplockRmsdComparisonChart: React.FC<SnaplockRmsdComparisonChartProps> = ({
  rmsdComparison,
  allAtomStructure,
  secondaryTypes = [],
  disulfideResidueSeqs = [],
  activeResidueIdx,
  onSelectResidue,
  compact = false
}) => {
  const [chartMode, setChartMode] = useState<RmsdChartViewMode>('multilayer');
  const [showPlddtOverlay, setShowPlddtOverlay] = useState<boolean>(true);
  const [thresholdMetric, setThresholdMetric] = useState<'backbone' | 'allAtom'>('allAtom');

  const disulfideSet = useMemo(() => new Set(disulfideResidueSeqs), [disulfideResidueSeqs]);

  const chartData: ChartDatum[] = useMemo(() => {
    return rmsdComparison.perResidueProfile.map((entry) => {
      const mIdx = entry.modelResSeq - 1;
      const aaMeta = SNAPLOCK_22_DATABASE.find((a) => a.code3 === entry.resName);
      const rot = allAtomStructure?.residueRotamers?.[mIdx];
      const sec = secondaryTypes[mIdx] || 'COIL';
      return {
        ...entry,
        xLabel: `${entry.code1}${entry.modelResSeq}`,
        aaNamePl: aaMeta?.namePl || entry.resName,
        rotamerLabel: rot?.rotamerName || aaMeta?.rotamerPreferences || 'Standard',
        secondary: sec,
        isDisulfideClamp:
          disulfideSet.has(entry.modelResSeq) ||
          entry.resName === 'CYS' ||
          entry.resName === 'SEC',
        isSelected: activeResidueIdx === mIdx
      };
    });
  }, [rmsdComparison.perResidueProfile, allAtomStructure, secondaryTypes, disulfideSet, activeResidueIdx]);

  // Biophysical & AlphaFold Confidence Analysis (per EMBL-EBI / AlphaFold DB guidelines)
  const analytics = useMemo(() => {
    const n = chartData.length;
    if (n === 0) {
      return {
        subAngstromBackboneCount: 0,
        subAngstromBackbonePct: 0,
        subAngstromAllAtomCount: 0,
        subAngstromAllAtomPct: 0,
        bestResidue: null as ChartDatum | null,
        worstResidue: null as ChartDatum | null,
        meanSidechainRmsd: 0,
        fracVeryHighPlddt: 0,
        fracConfidentPlddt: 0,
        fracLowPlddt: 0,
        fracDisorderedPlddt: 0,
        refStartSeq: 1,
        refEndSeq: 1,
        isSubsegmentAligned: false,
        isHighlyDisordered: false
      };
    }

    let subBb = 0;
    let subAll = 0;
    let scSum = 0;
    let scCnt = 0;
    let vh = 0;
    let conf = 0;
    let low = 0;
    let vlow = 0;

    let best = chartData[0];
    let worst = chartData[0];

    for (const d of chartData) {
      if (d.backboneRmsdAngstrom < 1.0) subBb++;
      if (d.allAtomRmsdAngstrom < 1.0) subAll++;
      if (d.matchedSidechainAtomCount > 0) {
        scSum += d.sidechainRmsdAngstrom;
        scCnt++;
      }

      if (d.allAtomRmsdAngstrom < best.allAtomRmsdAngstrom) best = d;
      if (d.allAtomRmsdAngstrom > worst.allAtomRmsdAngstrom) worst = d;

      if (d.refPlddt >= 90) vh++;
      else if (d.refPlddt >= 70) conf++;
      else if (d.refPlddt >= 50) low++;
      else vlow++;
    }

    const refStartSeq = chartData[0].refResSeq;
    const refEndSeq = chartData[chartData.length - 1].refResSeq;
    const isSubsegmentAligned =
      Math.abs(rmsdComparison.referenceLength - rmsdComparison.modelLength) >= 5;
    const fracDisordered = vlow / n;

    return {
      subAngstromBackboneCount: subBb,
      subAngstromBackbonePct: parseFloat(((subBb / n) * 100).toFixed(1)),
      subAngstromAllAtomCount: subAll,
      subAngstromAllAtomPct: parseFloat(((subAll / n) * 100).toFixed(1)),
      bestResidue: best,
      worstResidue: worst,
      meanSidechainRmsd: scCnt > 0 ? parseFloat((scSum / scCnt).toFixed(3)) : 0,
      fracVeryHighPlddt: parseFloat(((vh / n) * 100).toFixed(1)),
      fracConfidentPlddt: parseFloat(((conf / n) * 100).toFixed(1)),
      fracLowPlddt: parseFloat(((low / n) * 100).toFixed(1)),
      fracDisorderedPlddt: parseFloat((fracDisordered * 100).toFixed(1)),
      refStartSeq,
      refEndSeq,
      isSubsegmentAligned,
      isHighlyDisordered: fracDisordered >= 0.25 || rmsdComparison.meanRefPlddt < 65
    };
  }, [chartData, rmsdComparison]);

  const yMaxRmsd = useMemo(() => {
    let maxVal = 1.4;
    for (const d of chartData) {
      maxVal = Math.max(
        maxVal,
        d.allAtomRmsdAngstrom,
        d.backboneRmsdAngstrom,
        d.caDistAngstrom,
        chartMode === 'sidechain_vs_backbone' ? d.maxAtomDistAngstrom : 0
      );
    }
    return parseFloat((Math.ceil(maxVal * 1.18 * 10) / 10).toFixed(1));
  }, [chartData, chartMode]);

  const handleChartClick = (state: any) => {
    if (!onSelectResidue || !state) return;
    const payload = state?.activePayload?.[0]?.payload as ChartDatum | undefined;
    if (payload && typeof payload.modelResSeq === 'number') {
      onSelectResidue(payload.modelResSeq - 1);
    }
  };

  // Custom Scientific Biophysics Tooltip
  const renderCustomTooltip = ({ active, payload }: any) => {
    if (!active || !payload || !payload.length) return null;
    const datum = payload[0]?.payload as ChartDatum | undefined;
    if (!datum) return null;

    const plddtTier = getPlddtTierLabel(datum.refPlddt);

    return (
      <div className="bg-[#050711]/95 border border-zinc-700/90 rounded-lg p-3 shadow-2xl font-mono text-[11px] text-slate-200 min-w-[270px] space-y-2 tabular-nums">
        <div className="flex items-center justify-between gap-3 border-b border-zinc-800 pb-1.5">
          <div>
            <div className="font-extrabold text-white text-xs flex items-center gap-1.5">
              <span className="text-cyan-300">#{datum.modelResSeq}</span>
              <span>{datum.resName}</span>
              <span className="text-zinc-400 font-normal">({datum.aaNamePl})</span>
            </div>
            <div className="text-[10px] text-zinc-400">
              Pozycja AlphaFold/Ref: <strong className="text-purple-300">#{datum.refResSeq}</strong> • Rotamer:{' '}
              <strong className="text-indigo-300">{datum.rotamerLabel}</strong>
            </div>
          </div>
          {datum.isDisulfideClamp && (
            <span className="px-1.5 py-0.5 rounded bg-amber-950/90 border border-amber-600/70 text-amber-300 text-[9.5px] font-bold shrink-0">
              CYS/SEC
            </span>
          )}
        </div>

        <div className="space-y-1">
          <div className="flex items-center justify-between">
            <span className="text-zinc-400 flex items-center gap-1.5">
              <span className="w-2 h-2 rounded-xs bg-emerald-400 inline-block" />
              Backbone RMSD (N,CA,C,O):
            </span>
            <span
              className={`font-extrabold ${
                datum.backboneRmsdAngstrom < 1.0 ? 'text-emerald-300' : 'text-amber-300'
              }`}
            >
              {datum.backboneRmsdAngstrom.toFixed(2)} Å
            </span>
          </div>

          <div className="flex items-center justify-between">
            <span className="text-zinc-400 flex items-center gap-1.5">
              <span className="w-2 h-2 rounded-xs bg-sky-400 inline-block" />
              Odchylenie atomu Cα:
            </span>
            <span className="font-bold text-sky-300">{datum.caDistAngstrom.toFixed(2)} Å</span>
          </div>

          <div className="flex items-center justify-between">
            <span className="text-zinc-400 flex items-center gap-1.5">
              <span className="w-2 h-2 rounded-xs bg-indigo-400 inline-block" />
              Pełne All-Atom RMSD ({datum.matchedAtomCount} at.):
            </span>
            <span className="font-extrabold text-indigo-300">{datum.allAtomRmsdAngstrom.toFixed(2)} Å</span>
          </div>

          <div className="flex items-center justify-between">
            <span className="text-zinc-400 flex items-center gap-1.5">
              <span className="w-2 h-2 rounded-xs bg-fuchsia-400 inline-block" />
              Łańcuch boczny RMSD:
            </span>
            <span className="font-bold text-fuchsia-300">
              {datum.matchedSidechainAtomCount > 0
                ? `${datum.sidechainRmsdAngstrom.toFixed(2)} Å (${datum.matchedSidechainAtomCount} at.)`
                : 'Brak (-H)'}
            </span>
          </div>

          <div className="flex items-center justify-between pt-1 border-t border-zinc-800/80">
            <span className="text-zinc-400">Max rozbieżność atomu:</span>
            <span className="font-bold text-rose-300">
              {datum.maxAtomName} ({datum.maxAtomDistAngstrom.toFixed(2)} Å)
            </span>
          </div>

          <div className="flex items-center justify-between">
            <span className="text-zinc-400">AlphaFold B-factor (pLDDT):</span>
            <span className={`font-bold ${plddtTier.colorClass}`}>
              {datum.refPlddt.toFixed(1)} • {plddtTier.label}
            </span>
          </div>
        </div>
      </div>
    );
  };

  const chartHeight = compact ? 275 : 340;
  const showBrush = chartData.length > 28;

  return (
    <div className="bg-[#05060e] border border-zinc-800/90 rounded-xl overflow-hidden font-mono select-none">
      {/* TOP BAR: TITLE + VIEW MODE TABS + PLDDT TOGGLE */}
      <div className="px-4 py-3 bg-[#070915] border-b border-zinc-800/90 flex flex-wrap items-center justify-between gap-3">
        <div className="flex items-center gap-2.5">
          <div className="w-8 h-8 rounded-lg bg-indigo-950/90 border border-indigo-500/50 flex items-center justify-center text-indigo-300 shrink-0">
            <BarChart3 className="w-4 h-4" />
          </div>
          <div>
            <div className="flex items-center gap-2 flex-wrap">
              <h3 className="text-xs sm:text-sm font-bold text-white tracking-tight">
                Wykres Porównawczy RMSD Per-Aminokwas vs AlphaFold ({rmsdComparison.referenceId})
              </h3>
              <span className="text-[10px] px-2 py-0.5 rounded bg-emerald-950/90 border border-emerald-700/70 text-emerald-300 font-bold tabular-nums">
                {analytics.subAngstromBackbonePct}% reszt &lt; 1.0 Å (Szkielet)
              </span>
            </div>
            <p className="text-[10.5px] text-zinc-400 mt-0.5">
              Rozbieżności atomowe po superpozycji Kabscha-Horna SE(3) • pLDDT ekstrahowane z kolumny B-factor modelu referencyjnego
            </p>
          </div>
        </div>

        {/* Interactive Mode Switcher */}
        <div className="flex items-center gap-2 flex-wrap">
          <div className="flex items-center bg-black p-0.5 rounded-lg border border-zinc-800 text-[10.5px]">
            <button
              type="button"
              onClick={() => setChartMode('multilayer')}
              className={`px-2.5 py-1 rounded font-bold transition cursor-pointer flex items-center gap-1 ${
                chartMode === 'multilayer'
                  ? 'bg-indigo-950 border border-indigo-500/70 text-indigo-200'
                  : 'text-zinc-400 hover:text-white'
              }`}
            >
              <Layers className="w-3 h-3" />
              <span>Wielowarstwowy</span>
            </button>

            <button
              type="button"
              onClick={() => setChartMode('sidechain_vs_backbone')}
              className={`px-2.5 py-1 rounded font-bold transition cursor-pointer flex items-center gap-1 ${
                chartMode === 'sidechain_vs_backbone'
                  ? 'bg-fuchsia-950 border border-fuchsia-500/70 text-fuchsia-200'
                  : 'text-zinc-400 hover:text-white'
              }`}
            >
              <Atom className="w-3 h-3" />
              <span>Rotamery vs Szkielet</span>
            </button>

            <button
              type="button"
              onClick={() => setChartMode('threshold_bars')}
              className={`px-2.5 py-1 rounded font-bold transition cursor-pointer flex items-center gap-1 ${
                chartMode === 'threshold_bars'
                  ? 'bg-emerald-950 border border-emerald-500/70 text-emerald-200'
                  : 'text-zinc-400 hover:text-white'
              }`}
            >
              <GitCompare className="w-3 h-3" />
              <span>Klasyfikacja Progowa</span>
            </button>
          </div>

          {chartMode === 'threshold_bars' ? (
            <div className="flex items-center bg-black p-0.5 rounded-lg border border-zinc-800 text-[10px]">
              <button
                type="button"
                onClick={() => setThresholdMetric('allAtom')}
                className={`px-2 py-1 rounded font-bold cursor-pointer ${
                  thresholdMetric === 'allAtom' ? 'bg-zinc-800 text-white' : 'text-zinc-400 hover:text-white'
                }`}
              >
                All-Atom
              </button>
              <button
                type="button"
                onClick={() => setThresholdMetric('backbone')}
                className={`px-2 py-1 rounded font-bold cursor-pointer ${
                  thresholdMetric === 'backbone' ? 'bg-zinc-800 text-white' : 'text-zinc-400 hover:text-white'
                }`}
              >
                Szkielet
              </button>
            </div>
          ) : (
            <label className="flex items-center gap-1.5 px-2.5 py-1 rounded-lg bg-black border border-zinc-800 text-[10.5px] text-amber-200 cursor-pointer">
              <input
                type="checkbox"
                checked={showPlddtOverlay}
                onChange={(e) => setShowPlddtOverlay(e.target.checked)}
                className="accent-amber-500 rounded"
              />
              <span>Krzywa pLDDT</span>
            </label>
          )}
        </div>
      </div>

      {/* ALPHAFOLD SUB-DOMAIN / ISOFORM & DISORDER WARNINGS (EMBL-EBI PROTOCOL) */}
      {(analytics.isSubsegmentAligned || analytics.isHighlyDisordered) && (
        <div className="px-4 py-2 bg-amber-950/25 border-b border-amber-800/40 flex flex-wrap items-center justify-between gap-2 text-[10.5px]">
          {analytics.isSubsegmentAligned && (
            <div className="flex items-center gap-2 text-amber-200">
              <ShieldCheck className="w-3.5 h-3.5 text-amber-400 shrink-0" />
              <span>
                Dopasowanie podzbioru domenowego (Needleman-Wunsch): łańcuch modelu (
                <strong>{rmsdComparison.modelLength} aa</strong>) wyrównano do reszt{' '}
                <strong>
                  #{analytics.refStartSeq}–#{analytics.refEndSeq}
                </strong>{' '}
                pełnej sekwencji referencyjnej {rmsdComparison.referenceId} ({rmsdComparison.referenceLength} aa).
              </span>
            </div>
          )}
          {analytics.isHighlyDisordered && (
            <div className="flex items-center gap-2 text-rose-300">
              <AlertTriangle className="w-3.5 h-3.5 text-rose-400 shrink-0" />
              <span>
                Uwaga dot. nieuporządkowania (IDR): {analytics.fracDisorderedPlddt}% reszt posiada pLDDT &lt; 50 (niska
                pewność AlphaFold w kolumnie B-factor — zalecana analiza wyłącznie sztywnych domen rdzeniowych).
              </span>
            </div>
          )}
        </div>
      )}

      {/* RECHARTS CANVAS */}
      <div className="px-3 pt-3 pb-2 bg-[#04050c]">
        <div style={{ width: '100%', height: showBrush ? chartHeight + 34 : chartHeight }}>
          <ResponsiveContainer width="100%" height="100%">
            <ComposedChart
              data={chartData}
              margin={{ top: 10, right: showPlddtOverlay && chartMode !== 'threshold_bars' ? 14 : 18, left: 0, bottom: 4 }}
              onClick={handleChartClick}
            >
              <defs>
                <linearGradient id="snaplockAllAtomFill" x1="0" y1="0" x2="0" y2="1">
                  <stop offset="5%" stopColor="#818cf8" stopOpacity={0.38} />
                  <stop offset="95%" stopColor="#818cf8" stopOpacity={0.03} />
                </linearGradient>
              </defs>

              <CartesianGrid stroke="#181826" strokeDasharray="3 3" vertical={false} />

              <XAxis
                dataKey="xLabel"
                stroke="#71717a"
                tick={{ fill: '#a1a1aa', fontSize: 10, fontFamily: 'monospace' }}
                tickLine={{ stroke: '#27272a' }}
                interval={chartData.length > 90 ? Math.floor(chartData.length / 22) : chartData.length > 40 ? 2 : 0}
              />

              <YAxis
                yAxisId="rmsd"
                domain={[0, yMaxRmsd]}
                stroke="#94a3b8"
                tick={{ fill: '#cbd5e1', fontSize: 10, fontFamily: 'monospace' }}
                tickLine={{ stroke: '#27272a' }}
                width={46}
                tickFormatter={(v: number) => `${v.toFixed(1)}Å`}
              />

              {showPlddtOverlay && chartMode !== 'threshold_bars' && (
                <YAxis
                  yAxisId="plddt"
                  orientation="right"
                  domain={[0, 100]}
                  stroke="#fbbf24"
                  tick={{ fill: '#fcd34d', fontSize: 10, fontFamily: 'monospace' }}
                  tickLine={{ stroke: '#451a03' }}
                  width={44}
                  tickFormatter={(v: number) => `${v}`}
                />
              )}

              <Tooltip content={renderCustomTooltip} cursor={{ fill: 'rgba(56, 189, 248, 0.08)' }} />

              <Legend
                wrapperStyle={{ fontSize: '10.5px', fontFamily: 'monospace', paddingTop: '4px' }}
              />

              {/* Sub-Angstrom 1.0 Å Target Criterion */}
              <ReferenceLine
                yAxisId="rmsd"
                y={1.0}
                stroke="#10b981"
                strokeDasharray="5 4"
                strokeWidth={1.4}
                label={{
                  value: 'Próg 1.00 Å',
                  position: 'insideTopLeft',
                  fill: '#34d399',
                  fontSize: 10
                }}
              />

              {/* Global Mean Backbone RMSD Reference Line */}
              <ReferenceLine
                yAxisId="rmsd"
                y={rmsdComparison.backboneRmsdAngstrom}
                stroke="#22d3ee"
                strokeDasharray="2 3"
                strokeWidth={1.1}
                label={{
                  value: `Śr. BB: ${rmsdComparison.backboneRmsdAngstrom} Å`,
                  position: 'insideTopRight',
                  fill: '#67e8f9',
                  fontSize: 10
                }}
              />

              {chartMode === 'multilayer' && (
                <>
                  <Area
                    yAxisId="rmsd"
                    type="monotone"
                    dataKey="allAtomRmsdAngstrom"
                    name="All-Atom RMSD [Å]"
                    stroke="#818cf8"
                    strokeWidth={2}
                    fillOpacity={1}
                    fill="url(#snaplockAllAtomFill)"
                    activeDot={{ r: 4, stroke: '#ffffff', strokeWidth: 1.5 }}
                  />
                  <Bar
                    yAxisId="rmsd"
                    dataKey="backboneRmsdAngstrom"
                    name="Backbone N·CA·C·O [Å]"
                    fill="#10b981"
                    fillOpacity={0.82}
                    maxBarSize={14}
                    radius={[2, 2, 0, 0]}
                  >
                    {chartData.map((entry) => (
                      <Cell
                        key={`bb-${entry.modelResSeq}`}
                        fill={entry.isSelected ? '#22d3ee' : '#10b981'}
                        stroke={entry.isSelected ? '#ffffff' : 'none'}
                        strokeWidth={entry.isSelected ? 1.5 : 0}
                      />
                    ))}
                  </Bar>
                  <Line
                    yAxisId="rmsd"
                    type="monotone"
                    dataKey="caDistAngstrom"
                    name="Odchylenie Cα [Å]"
                    stroke="#38bdf8"
                    strokeWidth={1.8}
                    dot={chartData.length <= 55 ? { r: 2.3, fill: '#38bdf8' } : false}
                    activeDot={{ r: 4.5, fill: '#38bdf8', stroke: '#ffffff' }}
                  />
                  {showPlddtOverlay && (
                    <Line
                      yAxisId="plddt"
                      type="monotone"
                      dataKey="refPlddt"
                      name="AlphaFold pLDDT [0–100]"
                      stroke="#fbbf24"
                      strokeWidth={1.4}
                      strokeDasharray="4 3"
                      dot={false}
                    />
                  )}
                </>
              )}

              {chartMode === 'sidechain_vs_backbone' && (
                <>
                  <Bar
                    yAxisId="rmsd"
                    dataKey="backboneRmsdAngstrom"
                    name="Szkielet N·CA·C·O [Å]"
                    fill="#22d3ee"
                    maxBarSize={12}
                    radius={[2, 2, 0, 0]}
                  />
                  <Bar
                    yAxisId="rmsd"
                    dataKey="sidechainRmsdAngstrom"
                    name="Łańcuch Boczny R [Å]"
                    fill="#c084fc"
                    maxBarSize={12}
                    radius={[2, 2, 0, 0]}
                  />
                  <Line
                    yAxisId="rmsd"
                    type="monotone"
                    dataKey="maxAtomDistAngstrom"
                    name="Max Rozbieżność Atomu [Å]"
                    stroke="#fb7185"
                    strokeWidth={1.8}
                    dot={chartData.length <= 55 ? { r: 2.3, fill: '#fb7185' } : false}
                  />
                  {showPlddtOverlay && (
                    <Line
                      yAxisId="plddt"
                      type="monotone"
                      dataKey="refPlddt"
                      name="AlphaFold pLDDT [0–100]"
                      stroke="#fbbf24"
                      strokeWidth={1.3}
                      strokeDasharray="4 3"
                      dot={false}
                    />
                  )}
                </>
              )}

              {chartMode === 'threshold_bars' && (
                <Bar
                  yAxisId="rmsd"
                  dataKey={thresholdMetric === 'allAtom' ? 'allAtomRmsdAngstrom' : 'backboneRmsdAngstrom'}
                  name={
                    thresholdMetric === 'allAtom'
                      ? 'All-Atom RMSD [Å] (<1Å Zielony • 1–2Å Bursztynowy • >2Å Koralowy)'
                      : 'Backbone RMSD [Å] (<1Å Zielony • 1–2Å Bursztynowy • >2Å Koralowy)'
                  }
                  maxBarSize={18}
                  radius={[2, 2, 0, 0]}
                >
                  {chartData.map((entry) => {
                    const val =
                      thresholdMetric === 'allAtom' ? entry.allAtomRmsdAngstrom : entry.backboneRmsdAngstrom;
                    const color = val < 1.0 ? '#10b981' : val < 2.0 ? '#f59e0b' : '#f43f5e';
                    return (
                      <Cell
                        key={`th-${entry.modelResSeq}`}
                        fill={color}
                        stroke={entry.isSelected ? '#ffffff' : 'none'}
                        strokeWidth={entry.isSelected ? 1.8 : 0}
                      />
                    );
                  })}
                </Bar>
              )}

              {showBrush && (
                <Brush
                  dataKey="xLabel"
                  height={22}
                  stroke="#38bdf8"
                  fill="#090b16"
                  travellerWidth={8}
                  startIndex={0}
                  endIndex={Math.min(chartData.length - 1, 59)}
                />
              )}
            </ComposedChart>
          </ResponsiveContainer>
        </div>
      </div>

      {/* BOTTOM BIOPHYSICAL DIAGNOSTIC STRIP */}
      <div className="px-4 py-2.5 bg-[#070914] border-t border-zinc-800/90 grid grid-cols-2 lg:grid-cols-5 gap-2.5 text-[10.5px] tabular-nums">
        <div className="bg-black/60 border border-zinc-800/80 rounded px-2.5 py-1.5">
          <div className="text-zinc-500 uppercase text-[9.5px]">Zgodność Szkieletu &lt; 1.0 Å</div>
          <div className="font-extrabold text-emerald-300 text-xs mt-0.5 flex items-center gap-1">
            <CheckCircle2 className="w-3.5 h-3.5 text-emerald-400" />
            <span>
              {analytics.subAngstromBackboneCount}/{chartData.length} reszt ({analytics.subAngstromBackbonePct}%)
            </span>
          </div>
        </div>

        <div className="bg-black/60 border border-zinc-800/80 rounded px-2.5 py-1.5">
          <div className="text-zinc-500 uppercase text-[9.5px]">Śr. RMSD Łańcuchów Bocznych</div>
          <div className="font-extrabold text-fuchsia-300 text-xs mt-0.5">
            {analytics.meanSidechainRmsd} Å <span className="text-zinc-500 font-normal">(All: {rmsdComparison.allHeavyAtomRmsdAngstrom} Å)</span>
          </div>
        </div>

        <div className="bg-black/60 border border-zinc-800/80 rounded px-2.5 py-1.5">
          <div className="text-zinc-500 uppercase text-[9.5px]">Najlepsze Dopasowanie</div>
          <div className="font-bold text-cyan-300 text-xs mt-0.5 truncate">
            {analytics.bestResidue
              ? `#${analytics.bestResidue.modelResSeq} ${analytics.bestResidue.resName}: ${analytics.bestResidue.allAtomRmsdAngstrom.toFixed(2)} Å`
              : '—'}
          </div>
        </div>

        <div className="bg-black/60 border border-zinc-800/80 rounded px-2.5 py-1.5">
          <div className="text-zinc-500 uppercase text-[9.5px]">Max Odchylenie Atomowe</div>
          <div className="font-bold text-rose-300 text-xs mt-0.5 truncate">
            {analytics.worstResidue
              ? `#${analytics.worstResidue.modelResSeq} ${analytics.worstResidue.resName} (${analytics.worstResidue.maxAtomName}: ${analytics.worstResidue.maxAtomDistAngstrom.toFixed(2)} Å)`
              : '—'}
          </div>
        </div>

        <div className="bg-black/60 border border-zinc-800/80 rounded px-2.5 py-1.5 col-span-2 lg:col-span-1">
          <div className="text-zinc-500 uppercase text-[9.5px] flex items-center justify-between">
            <span>Rozkład pLDDT (AF)</span>
            <Activity className="w-3 h-3 text-amber-400" />
          </div>
          <div className="font-bold text-zinc-300 text-[10px] mt-0.5 truncate" title="Udział przedziałów pewności pLDDT z kolumny B-factor: >90 Bardzo wysoka / 70-90 Wysoka / <50 Nieuporządkowana">
            <span className="text-sky-300">&gt;90:{analytics.fracVeryHighPlddt}%</span>{' '}
            <span className="text-cyan-300">70+:{analytics.fracConfidentPlddt}%</span>{' '}
            <span className="text-rose-400">&lt;50:{analytics.fracDisorderedPlddt}%</span>
          </div>
        </div>
      </div>
    </div>
  );
};
