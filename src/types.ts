export interface Atom {
  id: string;
  element: 'C' | 'H' | 'O' | 'N' | 'S' | 'Se';
  x: number; // 3D coordinates relative to center of gravity
  y: number;
  z: number;
  mass: number;
}

export interface Bond {
  from: string;
  to: string;
  type: 'single' | 'double';
}

export interface AminoAcid {
  id: string; // e.g. ALA, VAL
  namePl: string;
  nameEn: string;
  code3: string;
  code1: string;
  type: 'hydrophobic' | 'polar' | 'acidic' | 'basic' | 'special';
  formula: string;
  mass: number; // Molar mass in g/mol
  description: string;
  color: string; // Aesthetic color matching its chemical category
  atoms: Atom[];
  bonds: Bond[];
  vanDerWaals: number;     // kJ/mol interactive volume
  hydrophobicity: number;  // scale -4.5 to 4.5
  ramachandranPhi?: number;
  ramachandranPsi?: number;
}

export type BenchmarkProteinId = '7nhu' | 'ubq' | 'lyz';

export interface BenchmarkProtein {
  id: BenchmarkProteinId;
  name: string;
  fullName: string;
  pdbId: string;
  afId: string;
  length: number;
  resolution: string;
  organism: string;
  description: string;
  sequence3Letter: string;
  afPlddt: number;
  disulfideBonds: [number, number][]; // 1-indexed residue pairs
  functionalSites: string[];
  optimalDielectric: number;
  optimalHbondWeight: number;
  optimalHydrophobicWeight: number;
  baseRmsd: number;
  calibratedRmsd: number;
  baseTmScore: number;
  calibratedTmScore: number;
}

export interface NhuCalibrationState {
  isCalibrated: boolean;
  calibrationLevel: number; // 0 to 100%
  dielectricConstant: number; // e.g. 78.4
  hydrogenBondWeight: number; // e.g. 1.05
  hydrophobicCoreWeight: number; // e.g. 1.12
  rmsdVs7nhu: number;
  tmScoreVs7nhu: number;
  rmsdVsAlphaFold: number;
  activeBenchmark: BenchmarkProteinId;
}

export interface ConnectedAminoAcid {
  id: string;
  aminoAcid: AminoAcid;
  fromPatch: number; // 1-248 of preceeding node
  toPatch: number;   // 1-248 of this node
  torsionAngle: number; // angle of twist in degrees (-180 to 180)
}

export interface SoundKey {
  digit: number | string;
  frequency: number;
  note: string;
  chord: string; // e.g., C-dur, F-dur
  description: string;
}
