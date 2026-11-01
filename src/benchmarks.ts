import { BenchmarkProtein, BenchmarkProteinId } from './types';

export const BENCHMARK_PROTEINS: Record<BenchmarkProteinId, BenchmarkProtein> = {
  '7nhu': {
    id: '7nhu',
    name: 'Insulina desB30',
    fullName: 'Human desB30 Insulin (PDB: 7NHU)',
    pdbId: '7NHU',
    afId: 'AF-P01308-F1',
    length: 50,
    resolution: '1.40 Å',
    organism: 'Homo sapiens',
    description: 'Hormon peptydowy regulujący metabolizm glukozy. Złoty standard o ultrawysokiej rozdzielczości krystalograficznej 1.40 Å z 3 mostkami dwusiarczkowymi.',
    sequence3Letter: 'GLY ILE VAL GLU GLN CYS CYS THR SER ILE CYS SER LEU TYR GLN LEU GLU ASN TYR CYS ASN PHE VAL ASN GLN HIS LEU CYS GLY SER HIS LEU VAL GLU ALA LEU TYR LEU VAL CYS GLY GLU ARG GLY PHE PHE TYR THR PRO LYS',
    afPlddt: 91.4,
    disulfideBonds: [
      [6, 11],   // A6-A11 intra-chain
      [7, 28],   // A7-B7 inter-chain
      [20, 40]   // A20-B19 inter-chain
    ],
    functionalSites: [
      'Chain A: reszty 1-21',
      'Chain B: reszty 22-50 (desB30)',
      'Mostek S-S A6-A11 (wewnątrzłańcuchowy)',
      'Mostki S-S A7-B7 i A20-B19 (międzyłańcuchowe)',
      'Kieszeń wiążąca receptor insulinowy: Tyr19, Phe24, Phe25'
    ],
    optimalDielectric: 78.4,
    optimalHbondWeight: 1.05,
    optimalHydrophobicWeight: 1.12,
    baseRmsd: 1.78,
    calibratedRmsd: 1.08,
    baseTmScore: 0.80,
    calibratedTmScore: 0.95
  },
  'ubq': {
    id: 'ubq',
    name: 'Ubikwityna',
    fullName: 'Human Ubiquitin (PDB: 1UBQ • AlphaFold AF-P0CG48)',
    pdbId: '1UBQ',
    afId: 'AF-P0CG48-F1',
    length: 76,
    resolution: '1.80 Å',
    organism: 'Homo sapiens',
    description: 'Kluczowy białkowy znacznik degradacji proteasomalnej i sygnalizacji komórkowej. Klasyczny motyw "beta-grasp" z 5-niciową harmonijką beta i helisą alfa.',
    sequence3Letter: 'MET GLN ILE PHE VAL LYS THR LEU THR GLY LYS THR ILE THR LEU GLU VAL GLU PRO SER ASP THR ILE GLU ASN VAL LYS ALA LYS ILE GLN ASP LYS GLU GLY ILE PRO PRO ASP GLN GLN ARG LEU ILE PHE ALA GLY LYS GLN LEU GLU ASP GLY ARG THR LEU SER ASP TYR ASN ILE GLN LYS GLU SER THR LEU HIS LEU VAL LEU ARG LEU ARG GLY GLY',
    afPlddt: 96.2,
    disulfideBonds: [], // Ubiquitin naturally has no cysteines
    functionalSites: [
      'Łata hydrofobowa wiążąca UBD: Leu8, Ile44, Val70',
      'K48: Rozgałęzienie poliubikwitynacji degradacyjnej 26S',
      'K63: Rozgałęzienie poliubikwitynacji sygnalizacyjnej i naprawy DNA',
      'C-końcowy ogon Gly75-Gly76: Miejsce wiązania izopeptydowego z substratem'
    ],
    optimalDielectric: 76.5,
    optimalHbondWeight: 1.15,
    optimalHydrophobicWeight: 1.25,
    baseRmsd: 1.62,
    calibratedRmsd: 0.94,
    baseTmScore: 0.84,
    calibratedTmScore: 0.97
  },
  'lyz': {
    id: 'lyz',
    name: 'Lizozym',
    fullName: 'Hen Egg-White Lysozyme / HEWL (PDB: 1LYZ • AlphaFold AF-P00698)',
    pdbId: '1LYZ',
    afId: 'AF-P00698-F1',
    length: 129,
    resolution: '2.00 Å',
    organism: 'Gallus gallus',
    description: 'Modelowy enzym hydrolityczny przecinający wiązania beta(1->4) w peptydoglikanie ścian bakterii. Posiada 4 kluczowe mostki dwusiarczkowe stabilizujące szczelinę aktywną.',
    sequence3Letter: 'LYS VAL PHE GLY ARG CYS GLU LEU ALA ALA ALA MET LYS ARG HIS GLY LEU ASP ASN TYR ARG GLY TYR SER LEU GLY ASN TRP VAL CYS ALA ALA LYS PHE GLU SER ASN PHE ASN THR GLN ALA THR ASN ARG ASN THR ASP GLY SER THR ASP TYR GLY ILE LEU GLN ILE ASN SER ARG TRP TRP CYS ASN ASP GLY ARG THR PRO GLY SER ARG ASN LEU CYS ASN ILE PRO CYS SER ALA LEU LEU SER SER ASP ILE THR ALA SER VAL ASN CYS ALA LYS LYS ILE VAL SER ASP GLY ASN GLY MET ASN ALA TRP VAL ALA TRP ARG ASN ARG CYS LYS GLY THR ASP VAL GLN ALA TRP ILE ARG GLY CYS ARG LEU',
    afPlddt: 94.8,
    disulfideBonds: [
      [6, 127],  // Cys6 - Cys127
      [30, 115], // Cys30 - Cys115
      [64, 80],  // Cys64 - Cys80
      [76, 94]   // Cys76 - Cys94
    ],
    functionalSites: [
      'Centrum katalityczne: Glu35 (donor protonu, kwas ogólny)',
      'Centrum katalityczne: Asp52 (nukleofil stabilizujący stan przejściowy)',
      'Szczelina wiążąca heksasacharyd: pod-miejsca A-F',
      '4 mostki S-S: C6-C127, C30-C115, C64-C80, C76-C94'
    ],
    optimalDielectric: 74.2,
    optimalHbondWeight: 1.18,
    optimalHydrophobicWeight: 1.22,
    baseRmsd: 1.88,
    calibratedRmsd: 1.14,
    baseTmScore: 0.79,
    calibratedTmScore: 0.94
  }
};
