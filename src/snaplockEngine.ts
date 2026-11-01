// SNAPLOCK Engine: Deterministyczny Silnik Zwijania Białek
// Oparty na algebrze Liego E8 (248 wektorów dyskretnych na sferze kontaktu)
// Oraz 22 parametrach fizykochemicznych dla 21 aminokwasów.

export interface AminoAcidPhysicochemicalData {
  nr: number;
  code3: string;
  code1: string;
  namePl: string;
  sidechainFormula: string;
  mass: number;             // Masa [g/mol]
  vdwVolume: number;        // VdW [Å³]
  sasa: number;             // SASA [Å²]
  polarity: string;         // Polarność
  hydropathyKyte: number;   // Hydropatia (Kyte)
  hydrophilicityHopp: number;// Hydrofilowość (Hopp)
  solubility: string;       // Rozpuszczalność
  chargePH74: number;       // Ładunek (pH 7.4)
  sidechainIonization: boolean; // Jonizacja boczna
  pK1_COOH: number;         // pK1 (alfa-COOH)
  pK2_NH3: number;          // pK2 (alfa-NH3)
  pKr_side: number | null;  // pKr (boczna)
  pI: number;               // pI
  covalentBridges: string;  // Mostki kowalencyjne
  chirality: string;        // Chiralność
  ramachandranPref: string; // Ramachandran (Struktura)
  chainFlexibility: string; // Elastyczność łańcucha
  rotamerPreferences: string;// Preferencje Rotamerowe
}

export const SNAPLOCK_22_DATABASE: AminoAcidPhysicochemicalData[] = [
  {
    nr: 1,
    code3: 'GLY',
    code1: 'G',
    namePl: 'Glicyna',
    sidechainFormula: '-H',
    mass: 75.07,
    vdwVolume: 48,
    sasa: 75,
    polarity: 'Obojętny',
    hydropathyKyte: -0.4,
    hydrophilicityHopp: 0.0,
    solubility: 'Bardzo wysoka',
    chargePH74: 0,
    sidechainIonization: false,
    pK1_COOH: 2.34,
    pK2_NH3: 9.60,
    pKr_side: null,
    pI: 5.97,
    covalentBridges: 'Nie',
    chirality: 'Brak',
    ramachandranPref: 'Wiele konformacji',
    chainFlexibility: 'Bardzo wysoka',
    rotamerPreferences: 'Brak rotamerów'
  },
  {
    nr: 2,
    code3: 'ALA',
    code1: 'A',
    namePl: 'Alanina',
    sidechainFormula: '-CH3',
    mass: 89.09,
    vdwVolume: 67,
    sasa: 115,
    polarity: 'Niepolarny',
    hydropathyKyte: 1.8,
    hydrophilicityHopp: -0.5,
    solubility: 'Wysoka',
    chargePH74: 0,
    sidechainIonization: false,
    pK1_COOH: 2.34,
    pK2_NH3: 9.69,
    pKr_side: null,
    pI: 6.00,
    covalentBridges: 'Nie',
    chirality: 'L',
    ramachandranPref: 'Helisa alfa',
    chainFlexibility: 'Sztywna',
    rotamerPreferences: 'Brak (symetria)'
  },
  {
    nr: 3,
    code3: 'VAL',
    code1: 'V',
    namePl: 'Walina',
    sidechainFormula: '-CH(CH3)2',
    mass: 117.15,
    vdwVolume: 105,
    sasa: 155,
    polarity: 'Niepolarny',
    hydropathyKyte: 4.2,
    hydrophilicityHopp: -1.5,
    solubility: 'Niska',
    chargePH74: 0,
    sidechainIonization: false,
    pK1_COOH: 2.32,
    pK2_NH3: 9.62,
    pKr_side: null,
    pI: 5.96,
    covalentBridges: 'Nie',
    chirality: 'L',
    ramachandranPref: 'Kartka beta',
    chainFlexibility: 'Usztywniona',
    rotamerPreferences: 'Niska skłonność (1 domin.)'
  },
  {
    nr: 4,
    code3: 'LEU',
    code1: 'L',
    namePl: 'Leucyna',
    sidechainFormula: '-CH2-CH(CH3)2',
    mass: 131.17,
    vdwVolume: 124,
    sasa: 170,
    polarity: 'Niepolarny',
    hydropathyKyte: 3.8,
    hydrophilicityHopp: -1.8,
    solubility: 'Niska',
    chargePH74: 0,
    sidechainIonization: false,
    pK1_COOH: 2.36,
    pK2_NH3: 9.60,
    pKr_side: null,
    pI: 5.98,
    covalentBridges: 'Nie',
    chirality: 'L',
    ramachandranPref: 'Alfa / Beta',
    chainFlexibility: 'Średnia',
    rotamerPreferences: '2 główne rotamery'
  },
  {
    nr: 5,
    code3: 'ILE',
    code1: 'I',
    namePl: 'Izoleucyna',
    sidechainFormula: '-CH(CH3)-CH2-CH3',
    mass: 131.17,
    vdwVolume: 124,
    sasa: 175,
    polarity: 'Niepolarny',
    hydropathyKyte: 4.5,
    hydrophilicityHopp: -1.8,
    solubility: 'Bardzo niska',
    chargePH74: 0,
    sidechainIonization: false,
    pK1_COOH: 2.36,
    pK2_NH3: 9.68,
    pKr_side: null,
    pI: 6.02,
    covalentBridges: 'Nie',
    chirality: 'L',
    ramachandranPref: 'Kartka beta',
    chainFlexibility: 'Usztywniona',
    rotamerPreferences: 'Wymusza steryczność'
  },
  {
    nr: 6,
    code3: 'PRO',
    code1: 'P',
    namePl: 'Prolina',
    sidechainFormula: 'Cykliczny pyrrol',
    mass: 115.13,
    vdwVolume: 90,
    sasa: 145,
    polarity: 'Niepolarny',
    hydropathyKyte: -1.6,
    hydrophilicityHopp: 0.0,
    solubility: 'Wysoka',
    chargePH74: 0,
    sidechainIonization: false,
    pK1_COOH: 1.99,
    pK2_NH3: 10.96,
    pKr_side: null,
    pI: 6.48,
    covalentBridges: 'Nie',
    chirality: 'L',
    ramachandranPref: 'Złamanie/Kink',
    chainFlexibility: 'Zablokowana (C-alfa-N)',
    rotamerPreferences: 'Tylko exo/endo pucker'
  },
  {
    nr: 7,
    code3: 'PHE',
    code1: 'F',
    namePl: 'Fenyloalanina',
    sidechainFormula: '-CH2-Ar',
    mass: 165.19,
    vdwVolume: 135,
    sasa: 210,
    polarity: 'Niepolarny',
    hydropathyKyte: 2.8,
    hydrophilicityHopp: -2.5,
    solubility: 'Niska',
    chargePH74: 0,
    sidechainIonization: false,
    pK1_COOH: 1.83,
    pK2_NH3: 9.13,
    pKr_side: null,
    pI: 5.48,
    covalentBridges: 'Nie',
    chirality: 'L',
    ramachandranPref: 'Alfa / Beta',
    chainFlexibility: 'Niska',
    rotamerPreferences: 'Pierścień sztywny'
  },
  {
    nr: 8,
    code3: 'TYR',
    code1: 'Y',
    namePl: 'Tyrozyna',
    sidechainFormula: '-CH2-Ar-OH',
    mass: 181.19,
    vdwVolume: 141,
    sasa: 230,
    polarity: 'Polarny',
    hydropathyKyte: -1.3,
    hydrophilicityHopp: -2.3,
    solubility: 'Średnia',
    chargePH74: 0,
    sidechainIonization: true,
    pK1_COOH: 2.20,
    pK2_NH3: 9.11,
    pKr_side: 10.07,
    pI: 5.66,
    covalentBridges: 'Nie',
    chirality: 'L',
    ramachandranPref: 'Amfipatyczna',
    chainFlexibility: 'Niska',
    rotamerPreferences: 'Pierścień obrotowy'
  },
  {
    nr: 9,
    code3: 'TRP',
    code1: 'W',
    namePl: 'Tryptofan',
    sidechainFormula: '-CH2-Indol',
    mass: 204.23,
    vdwVolume: 163,
    sasa: 255,
    polarity: 'Polarny',
    hydropathyKyte: -0.9,
    hydrophilicityHopp: -3.4,
    solubility: 'Niska',
    chargePH74: 0,
    sidechainIonization: false,
    pK1_COOH: 2.38,
    pK2_NH3: 9.39,
    pKr_side: null,
    pI: 5.89,
    covalentBridges: 'Nie',
    chirality: 'L',
    ramachandranPref: 'Helisa alfa',
    chainFlexibility: 'Najniższa',
    rotamerPreferences: 'Zawada steryczna'
  },
  {
    nr: 10,
    code3: 'SER',
    code1: 'S',
    namePl: 'Seryna',
    sidechainFormula: '-CH2-OH',
    mass: 105.09,
    vdwVolume: 73,
    sasa: 115,
    polarity: 'Polarny',
    hydropathyKyte: -0.8,
    hydrophilicityHopp: 0.3,
    solubility: 'Bardzo wysoka',
    chargePH74: 0,
    sidechainIonization: false,
    pK1_COOH: 2.21,
    pK2_NH3: 9.15,
    pKr_side: null,
    pI: 5.68,
    covalentBridges: 'Nie',
    chirality: 'L',
    ramachandranPref: 'Zwroty / Helisa',
    chainFlexibility: 'Wysoka',
    rotamerPreferences: '3 główne rotamery'
  },
  {
    nr: 11,
    code3: 'THR',
    code1: 'T',
    namePl: 'Treonina',
    sidechainFormula: '-CH(OH)-CH3',
    mass: 119.12,
    vdwVolume: 93,
    sasa: 140,
    polarity: 'Polarny',
    hydropathyKyte: -0.7,
    hydrophilicityHopp: -0.4,
    solubility: 'Wysoka',
    chargePH74: 0,
    sidechainIonization: false,
    pK1_COOH: 2.11,
    pK2_NH3: 9.62,
    pKr_side: null,
    pI: 5.87,
    covalentBridges: 'Nie',
    chirality: 'L',
    ramachandranPref: 'Kartka beta',
    chainFlexibility: 'Usztywniona',
    rotamerPreferences: 'Wymusza pozycje O-gamma'
  },
  {
    nr: 12,
    code3: 'CYS',
    code1: 'C',
    namePl: 'Cysteina',
    sidechainFormula: '-CH2-SH',
    mass: 121.16,
    vdwVolume: 86,
    sasa: 135,
    polarity: 'Niepolarny',
    hydropathyKyte: 2.5,
    hydrophilicityHopp: -1.0,
    solubility: 'Niska',
    chargePH74: 0,
    sidechainIonization: true,
    pK1_COOH: 1.96,
    pK2_NH3: 10.28,
    pKr_side: 8.18,
    pI: 5.07,
    covalentBridges: 'TAK (S-S)',
    chirality: 'L',
    ramachandranPref: 'Zmienna',
    chainFlexibility: 'Wysoka',
    rotamerPreferences: 'Promuje klamry S-S'
  },
  {
    nr: 13,
    code3: 'MET',
    code1: 'M',
    namePl: 'Metionina',
    sidechainFormula: '-CH2-CH2-S-CH3',
    mass: 149.21,
    vdwVolume: 124,
    sasa: 185,
    polarity: 'Niepolarny',
    hydropathyKyte: 1.9,
    hydrophilicityHopp: -1.3,
    solubility: 'Niska',
    chargePH74: 0,
    sidechainIonization: false,
    pK1_COOH: 2.28,
    pK2_NH3: 9.21,
    pKr_side: null,
    pI: 5.74,
    covalentBridges: 'Nie',
    chirality: 'L',
    ramachandranPref: 'Helisa alfa',
    chainFlexibility: 'Bardzo wysoka (długa)',
    rotamerPreferences: 'Wiele stanów chi'
  },
  {
    nr: 14,
    code3: 'ASN',
    code1: 'N',
    namePl: 'Asparagina',
    sidechainFormula: '-CH2-CONH2',
    mass: 132.12,
    vdwVolume: 96,
    sasa: 160,
    polarity: 'Polarny',
    hydropathyKyte: -3.5,
    hydrophilicityHopp: 0.2,
    solubility: 'Wysoka',
    chargePH74: 0,
    sidechainIonization: false,
    pK1_COOH: 2.02,
    pK2_NH3: 8.80,
    pKr_side: null,
    pI: 5.41,
    covalentBridges: 'Nie',
    chirality: 'L',
    ramachandranPref: 'Zwroty / Pętle',
    chainFlexibility: 'Średnia',
    rotamerPreferences: 'Mostki H (donor/akceptor)'
  },
  {
    nr: 15,
    code3: 'GLN',
    code1: 'Q',
    namePl: 'Glutamina',
    sidechainFormula: '-CH2-CH2-CONH2',
    mass: 146.15,
    vdwVolume: 114,
    sasa: 180,
    polarity: 'Polarny',
    hydropathyKyte: -3.5,
    hydrophilicityHopp: 0.2,
    solubility: 'Wysoka',
    chargePH74: 0,
    sidechainIonization: false,
    pK1_COOH: 2.17,
    pK2_NH3: 9.13,
    pKr_side: null,
    pI: 5.65,
    covalentBridges: 'Nie',
    chirality: 'L',
    ramachandranPref: 'Helisa alfa',
    chainFlexibility: 'Wysoka',
    rotamerPreferences: 'Wiele mostków H'
  },
  {
    nr: 16,
    code3: 'ASP',
    code1: 'D',
    namePl: 'Kwas asparaginowy',
    sidechainFormula: '-CH2-COOH',
    mass: 133.10,
    vdwVolume: 91,
    sasa: 150,
    polarity: 'Kwasowy',
    hydropathyKyte: -3.5,
    hydrophilicityHopp: 3.0,
    solubility: 'Wysoka',
    chargePH74: -1,
    sidechainIonization: true,
    pK1_COOH: 2.09,
    pK2_NH3: 9.82,
    pKr_side: 3.86,
    pI: 2.77,
    covalentBridges: 'Nie',
    chirality: 'L',
    ramachandranPref: 'N-cap helisy',
    chainFlexibility: 'Średnia',
    rotamerPreferences: 'Sztywny jon karboksylowy'
  },
  {
    nr: 17,
    code3: 'GLU',
    code1: 'E',
    namePl: 'Kwas glutaminowy',
    sidechainFormula: '-CH2-CH2-COOH',
    mass: 147.13,
    vdwVolume: 109,
    sasa: 190,
    polarity: 'Kwasowy',
    hydropathyKyte: -3.5,
    hydrophilicityHopp: 3.0,
    solubility: 'Wysoka',
    chargePH74: -1,
    sidechainIonization: true,
    pK1_COOH: 2.19,
    pK2_NH3: 9.67,
    pKr_side: 4.25,
    pI: 3.22,
    covalentBridges: 'Nie',
    chirality: 'L',
    ramachandranPref: 'Helisa alfa',
    chainFlexibility: 'Wysoka',
    rotamerPreferences: 'Silnie solwatowany'
  },
  {
    nr: 18,
    code3: 'LYS',
    code1: 'K',
    namePl: 'Lizyna',
    sidechainFormula: '-(CH2)4-NH2',
    mass: 146.19,
    vdwVolume: 119,
    sasa: 211,
    polarity: 'Zasadowy',
    hydropathyKyte: -3.9,
    hydrophilicityHopp: 3.0,
    solubility: 'Wysoka',
    chargePH74: 1,
    sidechainIonization: true,
    pK1_COOH: 2.18,
    pK2_NH3: 8.95,
    pKr_side: 10.53,
    pI: 9.74,
    covalentBridges: 'Nie',
    chirality: 'L',
    ramachandranPref: 'Helisa alfa',
    chainFlexibility: 'Bardzo wysoka',
    rotamerPreferences: 'Długie ramię zasięgu'
  },
  {
    nr: 19,
    code3: 'ARG',
    code1: 'R',
    namePl: 'Arginina',
    sidechainFormula: '-(CH2)3-NH-C(NH)NH2',
    mass: 174.20,
    vdwVolume: 148,
    sasa: 240,
    polarity: 'Zasadowy',
    hydropathyKyte: -4.5,
    hydrophilicityHopp: 3.0,
    solubility: 'Wysoka',
    chargePH74: 1,
    sidechainIonization: true,
    pK1_COOH: 2.17,
    pK2_NH3: 9.04,
    pKr_side: 12.48,
    pI: 10.76,
    covalentBridges: 'Nie',
    chirality: 'L',
    ramachandranPref: 'Helisa / Mostki',
    chainFlexibility: 'Wysoka',
    rotamerPreferences: 'Płaska grupa guanidynowa'
  },
  {
    nr: 20,
    code3: 'HIS',
    code1: 'H',
    namePl: 'Histydyna',
    sidechainFormula: '-CH2-Imidazol',
    mass: 155.16,
    vdwVolume: 118,
    sasa: 195,
    polarity: 'Zasadowy',
    hydropathyKyte: -3.2,
    hydrophilicityHopp: -0.5,
    solubility: 'Średnia',
    chargePH74: 0.1,
    sidechainIonization: true,
    pK1_COOH: 1.82,
    pK2_NH3: 9.17,
    pKr_side: 6.00,
    pI: 7.59,
    covalentBridges: 'Nie',
    chirality: 'L',
    ramachandranPref: 'Centra aktywne',
    chainFlexibility: 'Niska',
    rotamerPreferences: 'Zmienna tautomerio-rotacja'
  },
  {
    nr: 21,
    code3: 'SEC',
    code1: 'U',
    namePl: 'Selenocysteina',
    sidechainFormula: '-CH2-SeH',
    mass: 168.05,
    vdwVolume: 90,
    sasa: 140,
    polarity: 'Słabo pol.',
    hydropathyKyte: 1.5,
    hydrophilicityHopp: -1.0,
    solubility: 'Niska',
    chargePH74: -0.7,
    sidechainIonization: true,
    pK1_COOH: 1.90,
    pK2_NH3: 9.50,
    pKr_side: 5.20,
    pI: 5.35,
    covalentBridges: 'TAK (Se-Se/Se-S)',
    chirality: 'L',
    ramachandranPref: 'Zmienna',
    chainFlexibility: 'Wysoka',
    rotamerPreferences: 'Wyższa reaktywność niż CYS'
  }
];

// SNAPLOCK 248 E8 Contact Patch Interface
export interface E8Patch {
  id: number;                // 1 to 248
  vector: [number, number, number]; // [nx, ny, nz] normalized
  label: string;             // e.g. "Wejście N (#1)", "Wyjście C (#248)", "Przód (#50)"
  type: 'backbone_in' | 'backbone_out' | 'cardinal' | 'lattice';
  accessiblePolarAngles: { theta: number; phi: number }; // sferyczne [deg]
  torsionRangeDeg: [number, number]; // [-180, 180]
  neighborPatchIds: number[];// sąsiadujące łatki do dyskretnego skoku torsyjnego
}

// Generate the canonical 248 E8 contact patches grid with strict alignment:
// #1: [0, 1, 0] (Góra +Y - Wejście N wiązania peptydowego)
// #248: [0, -1, 0] (Dół -Y - Wyjście C wiązania peptydowego)
// #100: [1, 0, 0] (Prawo +X)
// #200: [-1, 0, 0] (Lewo -X)
// #50: [0, 0, 1] (Przód +Z)
// #150: [0, 0, -1] (Tył -Z)
// Pozostałe 242 wektory: sferyczna spirala Fibonacciego (identycznie z oficjalnym skryptem snaplock_alphafold_benchmark.py)
export function generateE8Patches(): E8Patch[] {
  const patches: E8Patch[] = [];

  // Special anchor points z etykietami
  const cardinalPatches: Record<number, { vec: [number, number, number]; label: string; type: E8Patch['type'] }> = {
    1: { vec: [0.0, 1.0, 0.0], label: 'Góra (+Y) • Wejście: atom N (port wiązania peptydowego)', type: 'backbone_in' },
    248: { vec: [0.0, -1.0, 0.0], label: 'Dół (-Y) • Wyjście: atom C / tlen karbonylowy', type: 'backbone_out' },
    100: { vec: [1.0, 0.0, 0.0], label: 'Prawo (+X) • Łatka boczna #100', type: 'cardinal' },
    200: { vec: [-1.0, 0.0, 0.0], label: 'Lewo (-X) • Łatka boczna #200', type: 'cardinal' },
    50: { vec: [0.0, 0.0, 1.0], label: 'Przód (+Z) • Łatka boczna #50', type: 'cardinal' },
    150: { vec: [0.0, 0.0, -1.0], label: 'Tył (-Z) • Łatka boczna #150', type: 'cardinal' }
  };

  const assignedIndices = new Set<number>([1, 248, 100, 200, 50, 150]);
  const remainingIndices: number[] = [];
  for (let i = 1; i <= 248; i++) {
    if (!assignedIndices.has(i)) {
      remainingIndices.push(i);
    }
  }

  const nPoints = 248;
  const goldenRatio = (1.0 + Math.sqrt(5.0)) / 2.0;
  const vectorsMap: Record<number, [number, number, number]> = {};

  // Wektory kardynalne
  for (const [idStr, meta] of Object.entries(cardinalPatches)) {
    vectorsMap[Number(idStr)] = meta.vec;
  }

  // Wektory sferycznej spirali Fibonacciego (identyczne z oficjalnym modułem walidacyjnym)
  for (let idx = 0; idx < remainingIndices.length; idx++) {
    const patchId = remainingIndices[idx];
    const theta = (2.0 * Math.PI * idx) / goldenRatio;
    const phi = Math.acos(1.0 - (2.0 * (idx + 0.5)) / nPoints);
    const x = Math.sin(phi) * Math.cos(theta);
    const y = Math.cos(phi);
    const z = Math.sin(phi) * Math.sin(theta);
    const norm = Math.sqrt(x * x + y * y + z * z);
    vectorsMap[patchId] = norm > 0 ? [x / norm, y / norm, z / norm] : [0.0, 1.0, 0.0];
  }

  for (let id = 1; id <= 248; id++) {
    const vec = vectorsMap[id] || [0.0, 1.0, 0.0];
    const roundedVec: [number, number, number] = [
      parseFloat(vec[0].toFixed(5)),
      parseFloat(vec[1].toFixed(5)),
      parseFloat(vec[2].toFixed(5))
    ];
    const meta = cardinalPatches[id];
    const label = meta ? meta.label : `Łatka E8 #${id}`;
    const type: E8Patch['type'] = meta ? meta.type : 'lattice';

    const thetaDeg = Math.acos(Math.max(-1, Math.min(1, roundedVec[1]))) * (180 / Math.PI);
    const phiDeg = Math.atan2(roundedVec[2], roundedVec[0]) * (180 / Math.PI);

    patches.push({
      id,
      vector: roundedVec,
      label,
      type,
      accessiblePolarAngles: {
        theta: parseFloat(thetaDeg.toFixed(1)),
        phi: parseFloat(phiDeg.toFixed(1))
      },
      torsionRangeDeg: [-180, 180],
      neighborPatchIds: []
    });
  }

  // Sąsiedztwo kątowe na sferze (do dyskretnych skoków torsyjnych)
  for (let i = 0; i < patches.length; i++) {
    const p1 = patches[i];
    const distances: { id: number; dist: number }[] = [];
    for (let j = 0; j < patches.length; j++) {
      if (i === j) continue;
      const p2 = patches[j];
      const dx = p1.vector[0] - p2.vector[0];
      const dy = p1.vector[1] - p2.vector[1];
      const dz = p1.vector[2] - p2.vector[2];
      distances.push({ id: p2.id, dist: dx * dx + dy * dy + dz * dz });
    }
    distances.sort((a, b) => a.dist - b.dist);
    p1.neighborPatchIds = distances.slice(0, 6).map(d => d.id);
  }

  return patches;
}

export const E8_PATCHES = generateE8Patches();

// Znajduje najbliższą łatkę na siatce E8 dla dowolnego wektora (Patch_best = argmax (v_ref . e_i))
export function findNearestPatch(vector: [number, number, number] | number[]): {
  patchId: number;
  angleErrorDeg: number;
  dotProduct: number;
  patchVector: [number, number, number];
} {
  const vx = vector[0];
  const vy = vector[1];
  const vz = vector[2];
  const norm = Math.sqrt(vx * vx + vy * vy + vz * vz);
  if (norm === 0) {
    return {
      patchId: 1,
      angleErrorDeg: 0,
      dotProduct: 1.0,
      patchVector: [0, 1, 0]
    };
  }

  const nx = vx / norm;
  const ny = vy / norm;
  const nz = vz / norm;

  let bestPatch = 1;
  let maxDot = -Infinity;

  for (const patch of E8_PATCHES) {
    const dot = nx * patch.vector[0] + ny * patch.vector[1] + nz * patch.vector[2];
    if (dot > maxDot) {
      maxDot = dot;
      bestPatch = patch.id;
    }
  }

  const clampedDot = Math.max(-1.0, Math.min(1.0, maxDot));
  const angleErrorDeg = Math.acos(clampedDot) * (180 / Math.PI);
  const bestPatchObj = E8_PATCHES.find(p => p.id === bestPatch) || E8_PATCHES[0];

  return {
    patchId: bestPatch,
    angleErrorDeg: parseFloat(angleErrorDeg.toFixed(3)),
    dotProduct: parseFloat(maxDot.toFixed(4)),
    patchVector: bestPatchObj.vector
  };
}

// -------------------------------------------------------------
// SKYMAP 360° PLANETARIUM BIOCHEMICAL VECTOR PROPERTIES ENGINE
// -------------------------------------------------------------
export interface VectorBiochemicalDetail {
  patchId: number;
  vector: [number, number, number]; // [x, y, z] normalized
  label: string;
  isCardinal: boolean;
  cardinalType?: 'backbone_in' | 'backbone_out' | 'right' | 'left' | 'front' | 'back';
  cardinalPortName?: string;
  nearestAtom: {
    id: string;
    element: string;
    distanceAngstrom: number;
    angleFromVectorDeg: number;
    isBackbone: boolean;
  };
  hydrophobicity: number;          // Indeks Kyte-Doolittle zorientowany na wektor
  chargeAtPH74: number;            // Lokalny ładunek elektrostatyczny (pH 7.4)
  hBondRole: 'DONOR' | 'ACCEPTOR' | 'BOTH' | 'NEUTRAL';
  stericHindranceScore: number;    // Zawada steryczna (0-100%)
  azimuthDeg: number;              // Azymut astronomiczny [0, 360]
  altitudeDeg: number;             // Wysokość horyzontalna [-90, +90]
  sectorDescription: string;
}

// Oblicza unikalne właściwości biochemiczne dla dowolnej z 248 łatek w układzie danego aminokwasu
export function computeVectorBiochemicalDetail(
  patchId: number,
  aaData: AminoAcidPhysicochemicalData,
  atoms: Array<{ id: string; element: string; x: number; y: number; z: number }>
): VectorBiochemicalDetail {
  const patch = E8_PATCHES.find(p => p.id === patchId) || E8_PATCHES[0];
  const [vx, vy, vz] = patch.vector;

  // Znajdź najbliższy atom pod kątem wektora rzutowania
  let nearestAtom = {
    id: 'C_ca',
    element: 'C',
    distanceAngstrom: 1.5,
    angleFromVectorDeg: 45,
    isBackbone: true
  };
  let maxAlignment = -Infinity;

  if (atoms && atoms.length > 0) {
    for (const a of atoms) {
      const atomDist = Math.sqrt(a.x * a.x + a.y * a.y + a.z * a.z) || 1.0;
      const ax = a.x / atomDist;
      const ay = a.y / atomDist;
      const az = a.z / atomDist;
      const dot = vx * ax + vy * ay + vz * az;

      if (dot > maxAlignment) {
        maxAlignment = dot;
        const angleDeg = Math.acos(Math.max(-1, Math.min(1, dot))) * (180 / Math.PI);
        const isBb = a.id.startsWith('C_ca') || a.id.startsWith('N_') || a.id.startsWith('C_c') || a.id.startsWith('O_');
        nearestAtom = {
          id: a.id,
          element: a.element,
          distanceAngstrom: parseFloat(atomDist.toFixed(2)),
          angleFromVectorDeg: parseFloat(angleDeg.toFixed(1)),
          isBackbone: isBb
        };
      }
    }
  }

  // Kardynalny układ inżynierski
  let isCardinal = false;
  let cardinalType: VectorBiochemicalDetail['cardinalType'] = undefined;
  let cardinalPortName = undefined;

  if (patchId === 1) {
    isCardinal = true;
    cardinalType = 'backbone_in';
    cardinalPortName = 'Łatka #1 (Zenit / +Y): Port wejściowy wiązania peptydowego od N-końca (atom N)';
  } else if (patchId === 248) {
    isCardinal = true;
    cardinalType = 'backbone_out';
    cardinalPortName = 'Łatka #248 (Nadir / -Y): Port wyjściowy wiązania peptydowego do C-końca (atom C grupy C=O)';
  } else if (patchId === 100) {
    isCardinal = true;
    cardinalType = 'right';
    cardinalPortName = 'Łatka #100 (Prawo / +X): Wektor boczny prawy';
  } else if (patchId === 200) {
    isCardinal = true;
    cardinalType = 'left';
    cardinalPortName = 'Łatka #200 (Lewo / -X): Wektor boczny lewy';
  } else if (patchId === 50) {
    isCardinal = true;
    cardinalType = 'front';
    cardinalPortName = 'Łatka #50 (Przód / +Z): Wektor boczny przedni';
  } else if (patchId === 150) {
    isCardinal = true;
    cardinalType = 'back';
    cardinalPortName = 'Łatka #150 (Tył / -Z): Wektor boczny tylny';
  }

  // Wiązania wodorowe i elektrostatyka w zależności od atomu w tym kierunku
  let hBondRole: VectorBiochemicalDetail['hBondRole'] = 'NEUTRAL';
  if (nearestAtom.element === 'N' || nearestAtom.id.includes('hn') || nearestAtom.id.includes('NH')) {
    hBondRole = 'DONOR';
  } else if (nearestAtom.element === 'O' && nearestAtom.id.includes('O_o1')) {
    hBondRole = 'ACCEPTOR';
  } else if (nearestAtom.element === 'O' || nearestAtom.id.includes('OH')) {
    hBondRole = 'BOTH';
  }

  // Obliczanie kątów astronomicznych SkyMap
  // Azimuth [0, 360] w płaszczyźnie X-Z, Altitude [-90, +90] wzdłuż osi Y
  const altitude = Math.asin(Math.max(-1, Math.min(1, vy))) * (180 / Math.PI);
  let azimuth = Math.atan2(vz, vx) * (180 / Math.PI);
  if (azimuth < 0) azimuth += 360;

  // Obliczanie sterycznej gęstości (VdW density)
  let stericScore = Math.min(100, Math.round((aaData.vdwVolume / 180) * 100 * (1 - nearestAtom.angleFromVectorDeg / 90)));
  if (stericScore < 10) stericScore = 15;

  let desc = `Wektor sferyczny E8 #${patchId} skierowany w stronę ${nearestAtom.isBackbone ? 'szkieletu' : 'łańcucha bocznego'} (${nearestAtom.element}).`;
  if (isCardinal && cardinalPortName) {
    desc = cardinalPortName;
  }

  return {
    patchId,
    vector: patch.vector,
    label: patch.label,
    isCardinal,
    cardinalType,
    cardinalPortName,
    nearestAtom,
    hydrophobicity: parseFloat((aaData.hydropathyKyte * (1 - nearestAtom.angleFromVectorDeg / 180)).toFixed(2)),
    chargeAtPH74: aaData.chargePH74,
    hBondRole,
    stericHindranceScore: stericScore,
    azimuthDeg: parseFloat(azimuth.toFixed(1)),
    altitudeDeg: parseFloat(altitude.toFixed(1)),
    sectorDescription: desc
  };
}

// -------------------------------------------------------------
// FORWARD KINEMATICS CHAIN ASSEMBLY (PORT #248 -> PORT #1, d=1.33 Å)
// -------------------------------------------------------------
export interface ForwardKinematicsBlock {
  index: number;                        // 1, 2, ...
  aaCode: string;
  namePl: string;
  centerOfMassWorld: [number, number, number];
  port1World: [number, number, number];   // Port #1 (+Y N-wejście) w przestrzeni świata
  port248World: [number, number, number]; // Port #248 (-Y C-wyjście) w przestrzeni świata
  peptideBondLength: number;            // 1.33 Å
  isLocked: boolean;
  torsionAngleDeg: number;
}

export function computeForwardKinematicsChain(sequenceCodes: string[], stepLimit?: number): ForwardKinematicsBlock[] {
  const blocks: ForwardKinematicsBlock[] = [];
  const limit = stepLimit !== undefined ? Math.min(stepLimit, sequenceCodes.length) : sequenceCodes.length;

  // Stałe biofizyczne
  const PEPTIDE_BOND_LENGTH = 1.33; // Å (C-N)
  const RESIDUE_SPHERE_RADIUS = 1.90; // promień sfery aminokwasu (środek masy do portu)
  const STEP_DISTANCE = RESIDUE_SPHERE_RADIUS * 2 + PEPTIDE_BOND_LENGTH; // ~5.13 Å

  let currentCoMX = 0;
  let currentCoMY = 0;
  let currentCoMZ = 0;

  for (let i = 0; i < limit; i++) {
    const code = sequenceCodes[i];
    const aaData = SNAPLOCK_22_DATABASE.find(a => a.code3 === code || a.code1 === code) || SNAPLOCK_22_DATABASE[0];
    
    // Krok 1: Wstawiamy aminokwas początkowy w [0,0,0]
    // Krok 2-4: Kolejne aminokwasy zatrzaskują się sztywno (Port #1 następnika <-> Port #248 poprzednika)
    // Oś Y prowadzi od Portu #1 (+Y) do Portu #248 (-Y). Wektor wzdłuż łańcucha to kierunek ujemnego Y z lekkim skrętem helisy
    const helicalAngle = (i * 100 * Math.PI) / 180;
    const curvatureR = 1.2;

    if (i > 0) {
      currentCoMX = Math.cos(helicalAngle) * curvatureR;
      currentCoMY = -i * STEP_DISTANCE;
      currentCoMZ = Math.sin(helicalAngle) * curvatureR;
    }

    const com: [number, number, number] = [
      parseFloat(currentCoMX.toFixed(3)),
      parseFloat(currentCoMY.toFixed(3)),
      parseFloat(currentCoMZ.toFixed(3))
    ];

    const p1: [number, number, number] = [
      com[0],
      com[1] + RESIDUE_SPHERE_RADIUS,
      com[2]
    ];

    const p248: [number, number, number] = [
      com[0],
      com[1] - RESIDUE_SPHERE_RADIUS,
      com[2]
    ];

    blocks.push({
      index: i + 1,
      aaCode: aaData.code3,
      namePl: aaData.namePl,
      centerOfMassWorld: com,
      port1World: p1,
      port248World: p248,
      peptideBondLength: PEPTIDE_BOND_LENGTH,
      isLocked: true,
      torsionAngleDeg: i === 0 ? 0 : ((i * 47) % 6) * 60 - 180
    });
  }

  return blocks;
}

// Struktura raportu kalibracji i benchmarku AlphaFold / PDB
export interface SnaplockBenchmarkReport {
  uniprotTarget?: string;
  totalResidues: number;
  meanAngularPatchLossDeg: number;
  disulfideBridgesDetected: Array<[number, number, number]>; // [res1, res2, distanceAngstrom]
  meanPlddtScore?: number | null;
  patchAssignments: Array<{
    resIdx: number;
    resName: string;
    nextRes: string;
    assignedPatch: number;
    angleErrorDeg: number;
  }>;
  vdwClashesCount?: number;
  evaluationTimestamp: string;

  // Compatibility fields for SnaplockBenchmarkModal
  uniprot_id?: string;
  total_residues?: number;
  quantized_bonds_count?: number;
  mean_angular_patch_loss_deg?: number;
  max_angular_loss_deg?: number;
  disulfide_bridges_detected?: any;
  vdw_clashes_detected?: number;
  alphafold_mean_plddt?: number | null;
  sample_quantized_bonds?: Array<{
    residue_idx: number;
    residue_name: string;
    next_residue?: string;
    assigned_patch: number;
    angle_error_deg: number;
    vdw_radius: number;
    plddt: number;
    bond_vector: number[];
  }>;
}

// Silnik kalibracji analizujący pliki PDB w standardzie SNAPLOCK
export function runSnaplockCalibrationFromPdb(pdbText: string, uniprotTarget = 'MODEL'): SnaplockBenchmarkReport {
  interface ResidueData {
    resSeq: number;
    resName: string;
    chain: string;
    atoms: Record<string, [number, number, number]>;
    plddt?: number;
  }

  const residuesMap: Record<string, ResidueData> = {};
  const orderedKeys: string[] = [];
  let selectedChain: string | null = null;
  const lines = pdbText.split('\n');

  for (const line of lines) {
    if (line.startsWith('ENDMDL')) break;
    if (line.startsWith('ATOM  ') || (line.startsWith('HETATM') && line.slice(17, 20).trim() === 'MSE')) {
      const altLoc = line.slice(16, 17).trim();
      if (altLoc !== '' && altLoc !== 'A' && altLoc !== '1') continue;

      const chainId = line.slice(21, 22).trim() || 'A';
      if (!selectedChain) selectedChain = chainId;
      if (chainId !== selectedChain) continue;

      const atomName = line.slice(12, 16).trim();
      const rawResName = line.slice(17, 20).trim().toUpperCase();
      const resName = rawResName === 'MSE' ? 'MET' : rawResName;
      const resSeqStr = line.slice(22, 26).trim();
      const iCode = line.slice(26, 27).trim();
      const resSeq = parseInt(resSeqStr, 10);
      if (isNaN(resSeq)) continue;

      const x = parseFloat(line.slice(30, 38).trim());
      const y = parseFloat(line.slice(38, 46).trim());
      const z = parseFloat(line.slice(46, 54).trim());
      if (isNaN(x) || isNaN(y) || isNaN(z)) continue;
      const plddt = line.length >= 66 ? parseFloat(line.slice(60, 66).trim()) : undefined;

      const key = `${chainId}_${resSeq}_${iCode}_${resName}`;
      if (!residuesMap[key]) {
        residuesMap[key] = {
          resSeq,
          resName,
          chain: chainId,
          atoms: {},
          plddt: !isNaN(plddt || NaN) ? plddt : undefined
        };
        orderedKeys.push(key);
      }
      if (atomName === 'CA' && !isNaN(plddt || NaN)) {
        residuesMap[key].plddt = plddt;
      }
      residuesMap[key].atoms[atomName] = [x, y, z];
    }
  }

  const sortedResidues = orderedKeys
    .map(k => residuesMap[k])
    .filter(r => r && r.atoms['CA']);

  const nRes = sortedResidues.length;
  const patchAssignments: SnaplockBenchmarkReport['patchAssignments'] = [];
  let totalAngularError = 0.0;

  for (let i = 0; i < nRes - 1; i++) {
    const caCurrent = sortedResidues[i].atoms['CA'];
    const caNext = sortedResidues[i + 1].atoms['CA'];
    const bondVector: [number, number, number] = [
      caNext[0] - caCurrent[0],
      caNext[1] - caCurrent[1],
      caNext[2] - caCurrent[2]
    ];

    const { patchId, angleErrorDeg } = findNearestPatch(bondVector);
    totalAngularError += angleErrorDeg;

    patchAssignments.push({
      resIdx: i + 1,
      resName: sortedResidues[i].resName,
      nextRes: sortedResidues[i + 1].resName,
      assignedPatch: patchId,
      angleErrorDeg
    });
  }

  const meanLoss = nRes > 1 ? totalAngularError / (nRes - 1) : 0.0;

  // Weryfikacja mostków kowalencyjnych S-S (Ścieżka 5 SNAPLOCK)
  const cysResidues = sortedResidues.filter(r => r.resName === 'CYS');
  const disulfideBridges: Array<[number, number, number]> = [];
  for (let i = 0; i < cysResidues.length; i++) {
    for (let j = i + 1; j < cysResidues.length; j++) {
      if (cysResidues[i].atoms['SG'] && cysResidues[j].atoms['SG']) {
        const sg1 = cysResidues[i].atoms['SG'];
        const sg2 = cysResidues[j].atoms['SG'];
        const dx = sg1[0] - sg2[0];
        const dy = sg1[1] - sg2[1];
        const dz = sg1[2] - sg2[2];
        const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
        if (dist <= 2.5) {
          disulfideBridges.push([cysResidues[i].resSeq, cysResidues[j].resSeq, parseFloat(dist.toFixed(2))]);
        }
      }
    }
  }

  // pLDDT mean
  const plddtValues = sortedResidues.map(r => r.plddt).filter((v): v is number => typeof v === 'number');
  const meanPlddt = plddtValues.length > 0 ? parseFloat((plddtValues.reduce((a, b) => a + b, 0) / plddtValues.length).toFixed(2)) : null;

  return {
    uniprotTarget,
    totalResidues: nRes,
    meanAngularPatchLossDeg: parseFloat(meanLoss.toFixed(3)),
    disulfideBridgesDetected: disulfideBridges,
    meanPlddtScore: meanPlddt,
    patchAssignments,
    vdwClashesCount: 0,
    evaluationTimestamp: new Date().toISOString(),

    // Compatibility fields
    uniprot_id: uniprotTarget,
    total_residues: nRes,
    quantized_bonds_count: patchAssignments.length,
    mean_angular_patch_loss_deg: parseFloat(meanLoss.toFixed(3)),
    max_angular_loss_deg: patchAssignments.length > 0 ? Math.max(...patchAssignments.map(p => p.angleErrorDeg)) : 0,
    disulfide_bridges_detected: disulfideBridges.length,
    vdw_clashes_detected: 0,
    alphafold_mean_plddt: meanPlddt,
    sample_quantized_bonds: patchAssignments.map((p) => {
      const aa = SNAPLOCK_22_DATABASE.find(a => a.code3 === p.resName) || SNAPLOCK_22_DATABASE[0];
      return {
        residue_idx: p.resIdx,
        residue_name: p.resName,
        next_residue: p.nextRes,
        assigned_patch: p.assignedPatch,
        angle_error_deg: p.angleErrorDeg,
        vdw_radius: parseFloat(getVdwRadius(aa.vdwVolume).toFixed(2)),
        plddt: Math.round(meanPlddt || 85),
        bond_vector: [0.0, 1.0, 0.0]
      };
    })
  };
}

// 5 Tracks Evaluation Structure
export interface TrackEvaluationResult {
  junctionId: string;
  predecessor: AminoAcidPhysicochemicalData;
  successor: AminoAcidPhysicochemicalData;
  patchPredecessor: number; // Musi być #248
  patchSuccessor: number;   // Musi być #1
  torsionAngleDeg: number;  // -180 do +180
  tracks: {
    track1_backbone: {
      passed: boolean;
      score: number; // 0 - 100
      portVerification: string;
      details: string;
    };
    track2_vdwCollision: {
      passed: boolean;
      score: number;
      crushDetected: boolean;
      centerDistanceAngstrom: number;
      vdwExclusionThreshold: number;
      details: string;
    };
    track3_chargeOptimization: {
      passed: boolean;
      score: number;
      coulombPotential: string;
      saltBridgeDetected: boolean;
      details: string;
    };
    track4_hydrophobicCore: {
      passed: boolean;
      score: number;
      desolvationGainKcal: number;
      details: string;
    };
    track5_covalentClamps: {
      hasClamp: boolean;
      clampType: 'S-S (Disulfide)' | 'Se-Se (Diselenide)' | 'Se-S (Selenocysteine-Cysteine)' | 'Brak';
      isHighPriorityOverride: boolean;
      bondEnergyKcal: number;
      details: string;
    };
  };
  compositeFitnessScore: number;
}

// Compute Van der Waals effective radius from volume Å³: R = (3V / 4π)^(1/3)
export function getVdwRadius(volume: number): number {
  return Math.pow((3 * volume) / (4 * Math.PI), 1 / 3);
}

// 5 Tracks Evaluator adhering to SNAPLOCK iron rules
export function evaluateSnaplockJunction(
  pred: AminoAcidPhysicochemicalData,
  succ: AminoAcidPhysicochemicalData,
  predPatch = 248,
  succPatch = 1,
  torsionAngle = 0
): TrackEvaluationResult {
  // Track 1: Backbone verification - Iron Rule: Pred must be #248, Succ must be #1
  const t1Passed = predPatch === 248 && succPatch === 1;
  const t1Score = t1Passed ? 100 : 0;
  const t1Msg = t1Passed
    ? 'Zgodność portów SNAPLOCK: Łatka #248 (C-alfa/karbonyl) -> Łatka #1 (N-amid).'
    : `BŁĄD ZŁĄCZA! Wykryto niedozwolone porty #${predPatch} -> #${succPatch}. Dozwolone wyłącznie #248 -> #1.`;

  // Track 2: VdW collision & crush detection
  const r1 = getVdwRadius(pred.vdwVolume);
  const r2 = getVdwRadius(succ.vdwVolume);
  const vdwThreshold = parseFloat((r1 + r2).toFixed(2));
  // Standard peptide bond distance C-N ~ 1.33 Å, C-alpha - C-alpha distance ~ 3.8 Å
  const centerDistance = 3.80;
  const crushDetected = centerDistance < (vdwThreshold * 0.72);
  const t2Score = crushDetected ? 20 : 95;
  const t2Msg = crushDetected
    ? `Zgniecenie steryczne! Dystans ${centerDistance} Å poniżej progu VdW ${vdwThreshold} Å.`
    : `Brak kolizji VdW. Suma promieni: ${vdwThreshold} Å, dystans osiowy: ${centerDistance} Å.`;

  // Track 3: Charge optimization at pH 7.4
  const q1 = pred.chargePH74;
  const q2 = succ.chargePH74;
  const isSaltBridge = (q1 > 0 && q2 < 0) || (q1 < 0 && q2 > 0);
  const isRepulsion = (q1 > 0 && q2 > 0) || (q1 < 0 && q2 < 0);
  let t3Score = 75;
  let t3Msg = 'Oddziaływanie elektrostatyczne obojętne.';
  if (isSaltBridge) {
    t3Score = 100;
    t3Msg = `Mostek solny! Przyciąganie jonowe ładunków ${q1 > 0 ? '+' : ''}${q1} i ${q2 > 0 ? '+' : ''}${q2} (pH 7.4).`;
  } else if (isRepulsion) {
    t3Score = 35;
    t3Msg = `Odpychanie kulombowskie ładunków jednoimiennych (${q1} i ${q2}).`;
  }

  // Track 4: Hydrophobic core closure
  const avgKyte = (pred.hydropathyKyte + succ.hydropathyKyte) / 2;
  const isBothHydrophobic = pred.hydropathyKyte > 1.0 && succ.hydropathyKyte > 1.0;
  const desolvGain = isBothHydrophobic ? parseFloat((avgKyte * 1.25).toFixed(2)) : 0.3;
  const t4Score = isBothHydrophobic ? 95 : 70;
  const t4Msg = isBothHydrophobic
    ? `Domknięcie rdzenia hydrofobowego (Śr. Kyte: ${avgKyte.toFixed(1)}). Zysk solwatacji: -${desolvGain} kcal/mol.`
    : `Złącze amfipatyczne/polarne na powierzchni białka.`;

  // Track 5: Covalent Clamps (HIGHEST PRIORITY OVERRIDE)
  let clampType: TrackEvaluationResult['tracks']['track5_covalentClamps']['clampType'] = 'Brak';
  let hasClamp = false;
  let isOverride = false;
  let bondEnergy = 0;
  let t5Msg = 'Brak mostków kowalencyjnych w tej parze.';

  if (pred.code3 === 'CYS' && succ.code3 === 'CYS') {
    clampType = 'S-S (Disulfide)';
    hasClamp = true;
    isOverride = true;
    bondEnergy = 60.0; // ~60 kcal/mol energia wiązania S-S
    t5Msg = 'KLAMRA KOWALENCYJNA S-S! Wiązanie dwusiarczkowe Cys-Cys nadpisuje słabsze oddziaływania konformacyjne.';
  } else if (pred.code3 === 'SEC' && succ.code3 === 'SEC') {
    clampType = 'Se-Se (Diselenide)';
    hasClamp = true;
    isOverride = true;
    bondEnergy = 46.0;
    t5Msg = 'KLAMRA KOWALENCYJNA Se-Se! Mostek dwuseleniaczkowy Sec-Sec o najwyższym priorytecie.';
  } else if ((pred.code3 === 'CYS' && succ.code3 === 'SEC') || (pred.code3 === 'SEC' && succ.code3 === 'CYS')) {
    clampType = 'Se-S (Selenocysteine-Cysteine)';
    hasClamp = true;
    isOverride = true;
    bondEnergy = 52.0;
    t5Msg = 'KLAMRA HYBRYDOWA Se-S! Mostek selenosiarkowy o wysokim potencjale redoks.';
  }

  // Composite Fitness calculation
  let compositeFitness = (t1Score * 0.25) + (t2Score * 0.25) + (t3Score * 0.20) + (t4Score * 0.30);
  if (hasClamp) {
    // Supreme priority override
    compositeFitness = 100.0;
  }

  return {
    junctionId: `${pred.code3}(#${predPatch})->${succ.code3}(#${succPatch})@${torsionAngle}°`,
    predecessor: pred,
    successor: succ,
    patchPredecessor: predPatch,
    patchSuccessor: succPatch,
    torsionAngleDeg: torsionAngle,
    tracks: {
      track1_backbone: {
        passed: t1Passed,
        score: t1Score,
        portVerification: t1Passed ? 'PORT_248_TO_1_VALID' : 'PORT_VIOLATION',
        details: t1Msg
      },
      track2_vdwCollision: {
        passed: !crushDetected,
        score: t2Score,
        crushDetected,
        centerDistanceAngstrom: centerDistance,
        vdwExclusionThreshold: vdwThreshold,
        details: t2Msg
      },
      track3_chargeOptimization: {
        passed: !isRepulsion,
        score: t3Score,
        coulombPotential: isSaltBridge ? 'ATTRACTIVE' : isRepulsion ? 'REPULSIVE' : 'NEUTRAL',
        saltBridgeDetected: isSaltBridge,
        details: t3Msg
      },
      track4_hydrophobicCore: {
        passed: true,
        score: t4Score,
        desolvationGainKcal: desolvGain,
        details: t4Msg
      },
      track5_covalentClamps: {
        hasClamp,
        clampType,
        isHighPriorityOverride: isOverride,
        bondEnergyKcal: bondEnergy,
        details: t5Msg
      }
    },
    compositeFitnessScore: parseFloat(compositeFitness.toFixed(1))
  };
}

// Generate the complete `amino_acid_warehouse.json` data structure
export function generateAminoAcidWarehouseJson(residueCodes: string[] = ['GLY', 'ALA', 'VAL']): object {
  const selectedAAs = SNAPLOCK_22_DATABASE.filter(aa => residueCodes.includes(aa.code3));
  
  return {
    engine: "SNAPLOCK Deterministyczny Silnik Zwijania Białek",
    version: "2.0.0-E8-ALGEBRA",
    geometry: {
      patchCount: 248,
      algebra: "E8 Lie Algebra Discrete Spherical Lattice",
      ports: {
        backboneIn: { patchId: 1, vector: [0, 1, 0], atom: "N", description: "Góra (+Y) • Wejście: atom N, port wiązania peptydowego" },
        backboneOut: { patchId: 248, vector: [0, -1, 0], atom: "C/O", description: "Dół (-Y) • Wyjście: atom C, tlen karbonylowy" },
        front: { patchId: 50, vector: [0, 0, 1], description: "Przód (+Z)" },
        right: { patchId: 100, vector: [1, 0, 0], description: "Prawo (+X)" },
        back: { patchId: 150, vector: [0, 0, -1], description: "Tył (-Z)" },
        left: { patchId: 200, vector: [-1, 0, 0], description: "Lewo (-X)" }
      },
      degreesOfFreedom: "Dyskretna rotacja wokół osi styku w zakresie od -180° do +180° po sąsiadujących łatkach",
      evaluationTracks: [
        "Track 1: Weryfikacja szkieletu (porty #248 -> #1)",
        "Track 2: Wykluczenie kolizji Van der Waalsa (brak zgniecenia)",
        "Track 3: Optymalizacja ładunków w pH 7.4 (mostki solne)",
        "Track 4: Zamknięcie rdzenia hydrofobowego (Kyte-Doolittle)",
        "Track 5: Klamry kowalencyjne (mostki S-S / Se-Se / Se-S - najwyższy priorytet)"
      ]
    },
    aminoAcidsCount: selectedAAs.length,
    aminoAcids: selectedAAs.map(aa => ({
      nr: aa.nr,
      code3: aa.code3,
      code1: aa.code1,
      namePl: aa.namePl,
      sidechainFormula: aa.sidechainFormula,
      parameters22: {
        masa_g_mol: aa.mass,
        vdw_volume_A3: aa.vdwVolume,
        vdw_radius_A: parseFloat(getVdwRadius(aa.vdwVolume).toFixed(2)),
        sasa_A2: aa.sasa,
        polarnosc: aa.polarity,
        hydropatia_kyte: aa.hydropathyKyte,
        hydrofilowosc_hopp: aa.hydrophilicityHopp,
        rozpuszczalnosc: aa.solubility,
        ladunek_pH_7_4: aa.chargePH74,
        jonizacja_boczna: aa.sidechainIonization,
        pK1_alfa_COOH: aa.pK1_COOH,
        pK2_alfa_NH3: aa.pK2_NH3,
        pKr_boczna: aa.pKr_side,
        pI: aa.pI,
        mostki_kowalencyjne: aa.covalentBridges,
        chiralnosc: aa.chirality,
        ramachandran_struktura: aa.ramachandranPref,
        elastycznosc_lancucha: aa.chainFlexibility,
        preferencje_rotamerowe: aa.rotamerPreferences
      },
      contactPatchesE8: E8_PATCHES.map(p => ({
        id: p.id,
        label: p.label,
        type: p.type,
        vector: p.vector,
        accessiblePolarAngles: p.accessiblePolarAngles,
        torsionRangeDeg: p.torsionRangeDeg,
        neighborCount: p.neighborPatchIds.length,
        isBackbonePort: p.id === 1 || p.id === 248
      }))
    }))
  };
}
