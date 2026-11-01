import { AminoAcid, SoundKey, Atom, Bond } from './types';

// Helper to create atomic masses
const ATOMIC_MASS: Record<string, number> = {
  C: 12.011,
  H: 1.008,
  O: 15.999,
  N: 14.007,
  S: 32.06,
  Se: 78.96
};

// Helper function to build standard amino acid backbone atoms and bonds
function makeChemicalStructure(
  sidechainId: string,
  sideAtoms: Omit<Atom, 'id' | 'mass'>[],
  sideBonds: Bond[]
): { atoms: Atom[]; bonds: Bond[] } {
  // Backbone atoms common to standard amino acids
  const backboneAtoms: Omit<Atom, 'id' | 'mass'>[] = [
    { element: 'C', x: 0.0, y: 0.0, z: 0.0 },       // C-alpha [0]
    { element: 'H', x: 0.0, y: 0.0, z: 1.0 },       // H-alpha [1]
    { element: 'N', x: -1.4, y: -0.4, z: -0.3 },    // Amino N [2]
    { element: 'H', x: -1.5, y: -1.4, z: -0.3 },    // Amino H1 [3]
    { element: 'H', x: -2.1, y: 0.1, z: 0.1 },      // Amino H2 [4]
    { element: 'C', x: 1.2, y: -0.8, z: -0.4 },     // Carboxyl C [5]
    { element: 'O', x: 1.3, y: -2.0, z: -0.2 },     // Carboxyl double-oxygen [6]
    { element: 'O', x: 2.2, y: -0.1, z: -0.9 },     // Carboxyl hydroxyl O [7]
    { element: 'H', x: 3.0, y: -0.6, z: -1.1 }      // Carboxyl Hydroxyl H [8]
  ];

  const backboneBonds: Bond[] = [
    { from: 'C_ca', to: 'H_ha', type: 'single' },
    { from: 'C_ca', to: 'N_n', type: 'single' },
    { from: 'N_n', to: 'H_hn1', type: 'single' },
    { from: 'N_n', to: 'H_hn2', type: 'single' },
    { from: 'C_ca', to: 'C_c', type: 'single' },
    { from: 'C_c', to: 'O_o1', type: 'double' },
    { from: 'C_c', to: 'O_o2', type: 'single' },
    { from: 'O_o2', to: 'H_ho2', type: 'single' }
  ];

  // Map backbone names
  const bNames = ['C_ca', 'H_ha', 'N_n', 'H_hn1', 'H_hn2', 'C_c', 'O_o1', 'O_o2', 'H_ho2'];
  const finalAtoms: Atom[] = backboneAtoms.map((ba, idx) => ({
    id: bNames[idx],
    element: ba.element as any,
    x: ba.x,
    y: ba.y,
    z: ba.z,
    mass: ATOMIC_MASS[ba.element] || 12
  }));

  // Append Sidechain Atoms
  sideAtoms.forEach((sa, idx) => {
    finalAtoms.push({
      id: `S_${idx}`,
      element: sa.element as any,
      x: sa.x,
      y: sa.y,
      z: sa.z,
      mass: ATOMIC_MASS[sa.element] || 12
    });
  });

  // Combine bonds
  const finalBonds = [...backboneBonds, ...sideBonds];

  // NOW - we satisfy the user request:
  // "środek ciężkości tej cząsteczki na dole a najlżejsza na górze"
  // Let's compute Center of Gravity (CoG)
  let totalMass = 0;
  let sX = 0, sY = 0, sZ = 0;
  finalAtoms.forEach(a => {
    totalMass += a.mass;
    sX += a.x * a.mass;
    sY += a.y * a.mass;
    sZ += a.z * a.mass;
  });

  const cogX = sX / totalMass;
  const cogY = sY / totalMass;
  const cogZ = sZ / totalMass;

  // Center the atoms around (0,0,0)
  finalAtoms.forEach(a => {
    a.x -= cogX;
    a.y -= cogY;
    a.z -= cogZ;
  });

  // Separate heavy atoms (C, N, O, S, Se) and light atoms (H) to orient weight-wise
  // Light atoms pointing UP (positive Y), heavy atoms pointing DOWN (negative Y)
  let heavyX = 0, heavyY = 0, heavyZ = 0, heavyCount = 0;
  let lightX = 0, lightY = 0, lightZ = 0, lightCount = 0;

  finalAtoms.forEach(a => {
    if (a.element === 'H') {
      lightX += a.x;
      lightY += a.y;
      lightZ += a.z;
      lightCount++;
    } else {
      heavyX += a.x;
      heavyY += a.y;
      heavyZ += a.z;
      heavyCount++;
    }
  });

  const hAvg = lightCount > 0 ? { x: lightX / lightCount, y: lightY / lightCount, z: lightZ / lightCount } : { x: 0, y: 1, z: 0 };
  const hvAvg = heavyCount > 0 ? { x: heavyX / heavyCount, y: heavyY / heavyCount, z: heavyZ / heavyCount } : { x: 0, y: -1, z: 0 };

  // Alignment vector: heavy -> light
  const alignX = hAvg.x - hvAvg.x;
  const alignY = hAvg.y - hvAvg.y;
  const alignZ = hAvg.z - hvAvg.z;
  const alignLen = Math.sqrt(alignX * alignX + alignY * alignY + alignZ * alignZ) || 1;

  // Let's rotate so align vector matches the positive Y axis (0, 1, 0)
  // Axis of rotation = Align x (0, 1, 0)
  const rotAxisX = alignZ / alignLen;
  const rotAxisY = 0;
  const rotAxisZ = -alignX / alignLen;
  const rotAxisLen = Math.sqrt(rotAxisX * rotAxisX + rotAxisZ * rotAxisZ) || 1;

  // Angle to rotate: cos(theta) = alignY / alignLen
  const cosTheta = alignY / alignLen;
  const sinTheta = Math.sqrt(1 - cosTheta * cosTheta);

  const ux = rotAxisX / rotAxisLen;
  const uz = rotAxisZ / rotAxisLen;

  // Rotate each atom
  finalAtoms.forEach(a => {
    // Rotation around arbitrary axis (ux, 0, uz) by theta
    const px = a.x;
    const py = a.y;
    const pz = a.z;

    // Standard 3D Rodrigues rotation matrix
    const c = cosTheta;
    const s = sinTheta;
    const t = 1 - c;

    const nx = (c + ux * ux * t) * px + (-uz * s) * py + (ux * uz * t) * pz;
    const ny = (uz * s) * px + c * py + (-ux * s) * pz;
    const nz = (ux * uz * t) * px + (ux * s) * py + (c + uz * uz * t) * pz;

    a.x = nx;
    a.y = ny;
    a.z = nz;
  });

  return { atoms: finalAtoms, bonds: finalBonds };
}

// 21 amino acids with distinct 3D atomic shapes and chemical bonds
export const AMINO_ACIDS: AminoAcid[] = [
  {
    id: 'ALA',
    namePl: 'Alanina',
    nameEn: 'Alanine',
    code3: 'Ala',
    code1: 'A',
    type: 'hydrophobic',
    formula: '-CH₃',
    mass: 89.1,
    description: 'Jeden z najprostszych i najpowszechniejszych aminokwasów. Stabilizuje helisy alfa w białku.',
    color: '#f59e0b',
    vanDerWaals: 88,
    hydrophobicity: 1.8,
    ...makeChemicalStructure(
      'ALA',
      [
        { element: 'C', x: 0.0, y: 1.5, z: 0.0 },  // S_0 Beta Carbon
        { element: 'H', x: -0.9, y: 1.9, z: 0.4 },  // S_1
        { element: 'H', x: 0.9, y: 1.9, z: -0.4 },  // S_2
        { element: 'H', x: 0.0, y: 1.6, z: -1.0 }   // S_3
      ],
      [
        { from: 'C_ca', to: 'S_0', type: 'single' },
        { from: 'S_0', to: 'S_1', type: 'single' },
        { from: 'S_0', to: 'S_2', type: 'single' },
        { from: 'S_0', to: 'S_3', type: 'single' }
      ]
    )
  },
  {
    id: 'ARG',
    namePl: 'Arginina',
    nameEn: 'Arginine',
    code3: 'Arg',
    code1: 'R',
    type: 'basic',
    formula: '-(CH₂)₃NH-C(NH)NH₂',
    mass: 174.2,
    description: 'Silnie zasadowy, naładowany dodatnio aminokwas. Tworzy stabilne mostki solne z resztami ujemnymi.',
    color: '#8b5cf6',
    vanDerWaals: 173,
    hydrophobicity: -4.5,
    ...makeChemicalStructure(
      'ARG',
      [
        { element: 'C', x: 0.0, y: 1.5, z: 0.0 },   // S_0 (Beta C)
        { element: 'C', x: 0.1, y: 2.7, z: -0.9 },  // S_1 (Gamma C)
        { element: 'C', x: -0.2, y: 3.9, z: -0.1 }, // S_2 (Delta C)
        { element: 'N', x: -0.1, y: 5.1, z: -0.9 }, // S_3 (Epsilon N)
        { element: 'C', x: 0.7, y: 6.1, z: -0.8 },  // S_4 (Zeta C)
        { element: 'N', x: 1.6, y: 6.1, z: -1.7 },  // S_5 (Imine N)
        { element: 'N', x: 0.6, y: 7.2, z: 0.0 }    // S_6 (Amino N)
      ],
      [
        { from: 'C_ca', to: 'S_0', type: 'single' },
        { from: 'S_0', to: 'S_1', type: 'single' },
        { from: 'S_1', to: 'S_2', type: 'single' },
        { from: 'S_2', to: 'S_3', type: 'single' },
        { from: 'S_3', to: 'S_4', type: 'single' },
        { from: 'S_4', to: 'S_5', type: 'double' },
        { from: 'S_4', to: 'S_6', type: 'single' }
      ]
    )
  },
  {
    id: 'ASN',
    namePl: 'Asparagina',
    nameEn: 'Asparagine',
    code3: 'Asn',
    code1: 'N',
    type: 'polar',
    formula: '-CH₂CONH₂',
    mass: 132.1,
    description: 'Polarny hydrofilowy o strukturze amidowej. Łatwo bierze udział w wiązaniach wodorowych na zewnątrz białka.',
    color: '#06b6d4',
    vanDerWaals: 114,
    hydrophobicity: -3.5,
    ...makeChemicalStructure(
      'ASN',
      [
        { element: 'C', x: 0.0, y: 1.5, z: 0.0 },   // S_0 (Beta C)
        { element: 'C', x: 1.2, y: 2.3, z: 0.3 },   // S_1 (Amide C)
        { element: 'O', x: 2.2, y: 1.9, z: -0.3 },  // S_2 (Amide O)
        { element: 'N', x: 1.1, y: 3.4, z: 1.1 }    // S_3 (Amino N)
      ],
      [
        { from: 'C_ca', to: 'S_0', type: 'single' },
        { from: 'S_0', to: 'S_1', type: 'single' },
        { from: 'S_1', to: 'S_2', type: 'double' },
        { from: 'S_1', to: 'S_3', type: 'single' }
      ]
    )
  },
  {
    id: 'ASP',
    namePl: 'Kwas asparaginowy',
    nameEn: 'Aspartic Acid',
    code3: 'Asp',
    code1: 'D',
    type: 'acidic',
    formula: '-CH₂COOH',
    mass: 133.1,
    description: 'Kwas dakarboksylowy, naładowany ujemnie w pH 7. Służy do przewodnictwa i mostków jonowych.',
    color: '#ef4444',
    vanDerWaals: 111,
    hydrophobicity: -3.5,
    ...makeChemicalStructure(
      'ASP',
      [
        { element: 'C', x: 0.0, y: 1.5, z: 0.0 },   // S_0
        { element: 'C', x: 1.2, y: 2.3, z: 0.4 },   // S_1
        { element: 'O', x: 2.3, y: 2.1, z: -0.2 },  // S_2 (carbonyl)
        { element: 'O', x: 1.0, y: 3.3, z: 1.3 }    // S_3 (hydroxyl)
      ],
      [
        { from: 'C_ca', to: 'S_0', type: 'single' },
        { from: 'S_0', to: 'S_1', type: 'single' },
        { from: 'S_1', to: 'S_2', type: 'double' },
        { from: 'S_1', to: 'S_3', type: 'single' }
      ]
    )
  },
  {
    id: 'CYS',
    namePl: 'Cysteina',
    nameEn: 'Cysteine',
    code3: 'Cys',
    code1: 'C',
    type: 'special',
    formula: '-CH₂SH',
    mass: 121.2,
    description: 'Kluczowa dla trzeciorzędowej formy białka. Posiada grupę tiolową, która paruje się tworząc mostki dwusiarczkowe S-S.',
    color: '#10b981',
    vanDerWaals: 104,
    hydrophobicity: 2.5,
    ...makeChemicalStructure(
      'CYS',
      [
        { element: 'C', x: 0.0, y: 1.5, z: 0.0 },  // S_0 C-beta
        { element: 'S', x: 1.2, y: 2.7, z: -0.5 }  // S_1 Siarka
      ],
      [
        { from: 'C_ca', to: 'S_0', type: 'single' },
        { from: 'S_0', to: 'S_1', type: 'single' }
      ]
    )
  },
  {
    id: 'GLU',
    namePl: 'Kwas glutaminowy',
    nameEn: 'Glutamic Acid',
    code3: 'Glu',
    code1: 'E',
    type: 'acidic',
    formula: '-(CH₂)₂COOH',
    mass: 147.1,
    description: 'Kwasowy hydrofilowy aminokwas o dłuższym łańcuchu. Bardzo ważny neurotransmiter w organizmie.',
    color: '#f43f5e',
    vanDerWaals: 138,
    hydrophobicity: -3.5,
    ...makeChemicalStructure(
      'GLU',
      [
        { element: 'C', x: 0.0, y: 1.5, z: 0.0 },   // Beta C
        { element: 'C', x: 0.1, y: 2.7, z: -0.9 },  // Gamma C
        { element: 'C', x: 1.3, y: 3.5, z: -0.5 },  // Delta C
        { element: 'O', x: 2.4, y: 3.1, z: -1.0 },  // O Double
        { element: 'O', x: 1.2, y: 4.6, z: 0.3 }    // O Single
      ],
      [
        { from: 'C_ca', to: 'S_0', type: 'single' },
        { from: 'S_0', to: 'S_1', type: 'single' },
        { from: 'S_1', to: 'S_2', type: 'single' },
        { from: 'S_2', to: 'S_3', type: 'double' },
        { from: 'S_2', to: 'S_4', type: 'single' }
      ]
    )
  },
  {
    id: 'GLN',
    namePl: 'Glutamina',
    nameEn: 'Glutamine',
    code3: 'Gln',
    code1: 'Q',
    type: 'polar',
    formula: '-(CH₂)₂CONH₂',
    mass: 146.1,
    description: 'Polarny aminokwas obficie występujący w mięśniach. Odpowiada za transport wolnego amoniaku.',
    color: '#3b82f6',
    vanDerWaals: 144,
    hydrophobicity: -3.5,
    ...makeChemicalStructure(
      'GLN',
      [
        { element: 'C', x: 0.0, y: 1.5, z: 0.0 },   // S_0
        { element: 'C', x: 0.1, y: 2.7, z: -0.9 },  // S_1
        { element: 'C', x: 1.3, y: 3.5, z: -0.4 },  // S_2
        { element: 'O', x: 2.4, y: 3.1, z: -0.9 },  // S_3
        { element: 'N', x: 1.1, y: 4.6, z: 0.4 }    // S_4
      ],
      [
        { from: 'C_ca', to: 'S_0', type: 'single' },
        { from: 'S_0', to: 'S_1', type: 'single' },
        { from: 'S_1', to: 'S_2', type: 'single' },
        { from: 'S_2', to: 'S_3', type: 'double' },
        { from: 'S_2', to: 'S_4', type: 'single' }
      ]
    )
  },
  {
    id: 'GLY',
    namePl: 'Glicyna',
    nameEn: 'Glycine',
    code3: 'Gly',
    code1: 'G',
    type: 'special',
    formula: '-H',
    mass: 75.1,
    description: 'Najmniejszy, achiralny i najbardziej gibki aminokwas. Pozwala białku na skręty o bardzo ostrych kątach.',
    color: '#14b8a6',
    vanDerWaals: 60,
    hydrophobicity: -0.4,
    ...makeChemicalStructure(
      'GLY',
      [
        { element: 'H', x: 0.0, y: 1.2, z: 0.0 }  // S_0
      ],
      [
        { from: 'C_ca', to: 'S_0', type: 'single' }
      ]
    )
  },
  {
    id: 'HIS',
    namePl: 'Histydyna',
    nameEn: 'Histidine',
    code3: 'His',
    code1: 'H',
    type: 'basic',
    formula: '-CH₂-C₃H₃N₂',
    mass: 155.2,
    description: 'Posiada dodatni ładunek o pH fizjologicznym. Tworzy pierścień imidazolowy zaangażowany w katalizowanie reakcji enzymatycznych.',
    color: '#6366f1',
    vanDerWaals: 153,
    hydrophobicity: -3.2,
    ...makeChemicalStructure(
      'HIS',
      [
        { element: 'C', x: 0.0, y: 1.5, z: 0.0 },   // S_0 Beta-C
        { element: 'C', x: -0.5, y: 2.7, z: -0.7 },  // S_1 ring Carbon
        { element: 'N', x: -1.7, y: 3.1, z: -0.5 },  // S_2 ring Nitrogen
        { element: 'C', x: -1.8, y: 4.1, z: -1.5 },  // S_3 ring Carbon
        { element: 'N', x: -0.7, y: 4.4, z: -2.3 },  // S_4 ring Nitrogen
        { element: 'C', x: 0.1, y: 3.5, z: -1.8 }    // S_5 ring Carbon
      ],
      [
        { from: 'C_ca', to: 'S_0', type: 'single' },
        { from: 'S_0', to: 'S_1', type: 'single' },
        { from: 'S_1', to: 'S_2', type: 'single' },
        { from: 'S_2', to: 'S_3', type: 'double' },
        { from: 'S_3', to: 'S_4', type: 'single' },
        { from: 'S_4', to: 'S_5', type: 'single' },
        { from: 'S_5', to: 'S_1', type: 'double' }
      ]
    )
  },
  {
    id: 'ILE',
    namePl: 'Izoleucyna',
    nameEn: 'Isoleucine',
    code3: 'Ile',
    code1: 'I',
    type: 'hydrophobic',
    formula: '-CH(CH₃)CH₂CH₃',
    mass: 131.2,
    description: 'Doskonałe geometryczne wypełnienie wnętrza hydrofobowego rdzenia białek.',
    color: '#f59e0b',
    vanDerWaals: 147,
    hydrophobicity: 4.5,
    ...makeChemicalStructure(
      'ILE',
      [
        { element: 'C', x: 0.0, y: 1.5, z: 0.0 },   // Beta C
        { element: 'C', x: -1.1, y: 2.3, z: -0.6 },  // Gamma C1
        { element: 'C', x: 1.3, y: 2.3, z: -0.1 },   // Gamma C2
        { element: 'C', x: 1.3, y: 3.6, z: -0.9 }    // Delta C
      ],
      [
        { from: 'C_ca', to: 'S_0', type: 'single' },
        { from: 'S_0', to: 'S_1', type: 'single' },
        { from: 'S_0', to: 'S_2', type: 'single' },
        { from: 'S_2', to: 'S_3', type: 'single' }
      ]
    )
  },
  {
    id: 'LEU',
    namePl: 'Leucyna',
    nameEn: 'Leucine',
    code3: 'Leu',
    code1: 'L',
    type: 'hydrophobic',
    formula: '-CH₂CH(CH₃)₂',
    mass: 131.2,
    description: 'Bogato występuje w strukturach białkowych mięśni gładkich. Stabilizuje hydrofobowy rdzeń.',
    color: '#eab308',
    vanDerWaals: 147,
    hydrophobicity: 3.8,
    ...makeChemicalStructure(
      'LEU',
      [
        { element: 'C', x: 0.0, y: 1.5, z: 0.0 },   // Beta C
        { element: 'C', x: 0.1, y: 2.8, z: -0.8 },  // Gamma C
        { element: 'C', x: -1.1, y: 3.6, z: -0.7 }, // Delta C1
        { element: 'C', x: 1.3, y: 3.6, z: -0.3 }   // Delta C2
      ],
      [
        { from: 'C_ca', to: 'S_0', type: 'single' },
        { from: 'S_0', to: 'S_1', type: 'single' },
        { from: 'S_1', to: 'S_2', type: 'single' },
        { from: 'S_1', to: 'S_3', type: 'single' }
      ]
    )
  },
  {
    id: 'LYS',
    namePl: 'Lizyna',
    nameEn: 'Lysine',
    code3: 'Lys',
    code1: 'K',
    type: 'basic',
    formula: '-(CH₂)₄NH₂',
    mass: 146.2,
    description: 'Długi, giętki, naładowany dodatnio łańcuch. Stabilizuje histonowe pętle DNA.',
    color: '#a855f7',
    vanDerWaals: 162,
    hydrophobicity: -3.9,
    ...makeChemicalStructure(
      'LYS',
      [
        { element: 'C', x: 0.0, y: 1.5, z: 0.0 },   // Beta C
        { element: 'C', x: 0.1, y: 2.7, z: -0.9 },  // Gamma C
        { element: 'C', x: -0.2, y: 3.9, z: -0.1 }, // Delta C
        { element: 'C', x: -0.1, y: 5.1, z: -0.9 }, // Epsilon C
        { element: 'N', x: -0.4, y: 6.3, z: -0.1 }  // Amine Nitrogen
      ],
      [
        { from: 'C_ca', to: 'S_0', type: 'single' },
        { from: 'S_0', to: 'S_1', type: 'single' },
        { from: 'S_1', to: 'S_2', type: 'single' },
        { from: 'S_2', to: 'S_3', type: 'single' },
        { from: 'S_3', to: 'S_4', type: 'single' }
      ]
    )
  },
  {
    id: 'MET',
    namePl: 'Metionina',
    nameEn: 'Methionine',
    code3: 'Met',
    code1: 'M',
    type: 'hydrophobic',
    formula: '-(CH₂)₂SCH₃',
    mass: 149.2,
    description: 'Pierwszy aminokwas inicjujący syntezę każdego łańcucha białkowego (kodon start AUG).',
    color: '#d97706',
    vanDerWaals: 160,
    hydrophobicity: 1.9,
    ...makeChemicalStructure(
      'MET',
      [
        { element: 'C', x: 0.0, y: 1.5, z: 0.0 },   // Beta C
        { element: 'C', x: 0.1, y: 2.7, z: -0.9 },  // Gamma C
        { element: 'S', x: -1.3, y: 3.7, z: -0.5 }, // Siarka Ether
        { element: 'C', x: -1.1, y: 4.8, z: 0.9 }   // Methyl Carbon terminal
      ],
      [
        { from: 'C_ca', to: 'S_0', type: 'single' },
        { from: 'S_0', to: 'S_1', type: 'single' },
        { from: 'S_1', to: 'S_2', type: 'single' },
        { from: 'S_2', to: 'S_3', type: 'single' }
      ]
    )
  },
  {
    id: 'PHE',
    namePl: 'Fenyloalanina',
    nameEn: 'Phenylalanine',
    code3: 'Phe',
    code1: 'F',
    type: 'hydrophobic',
    formula: '-CH₂-C₆H₅',
    mass: 165.2,
    description: 'Duża hydrofobowa grupa aromatyczna tworząca stabilne ułożenia pierścieniowe w rdzeniu białka.',
    color: '#ca8a04',
    vanDerWaals: 190,
    hydrophobicity: 2.8,
    ...makeChemicalStructure(
      'PHE',
      [
        { element: 'C', x: 0.0, y: 1.5, z: 0.0 },   // Beta C
        { element: 'C', x: -0.3, y: 2.7, z: -0.8 },  // Ring C1
        { element: 'C', x: -1.5, y: 3.4, z: -0.7 },  // Ring C2
        { element: 'C', x: -1.8, y: 4.5, z: -1.5 },  // Ring C3
        { element: 'C', x: -0.8, y: 4.9, z: -2.3 },  // Ring C4
        { element: 'C', x: 0.4, y: 4.2, z: -2.4 },   // Ring C5
        { element: 'C', x: 0.7, y: 3.1, z: -1.6 }    // Ring C6
      ],
      [
        { from: 'C_ca', to: 'S_0', type: 'single' },
        { from: 'S_0', to: 'S_1', type: 'single' },
        { from: 'S_1', to: 'S_2', type: 'double' },
        { from: 'S_2', to: 'S_3', type: 'single' },
        { from: 'S_3', to: 'S_4', type: 'double' },
        { from: 'S_4', to: 'S_5', type: 'single' },
        { from: 'S_5', to: 'S_6', type: 'double' },
        { from: 'S_6', to: 'S_1', type: 'single' }
      ]
    )
  },
  {
    id: 'PRO',
    namePl: 'Prolina',
    nameEn: 'Proline',
    code3: 'Pro',
    code1: 'P',
    type: 'special',
    formula: '-C₃H₆ cyclic',
    mass: 115.1,
    description: 'Ma sztywny pierścień łączący się z azotem szkieletu ("imonokwas"). Wymusza ostre zagięcia łańcucha.',
    color: '#ec4899',
    vanDerWaals: 112,
    hydrophobicity: -1.6,
    ...makeChemicalStructure(
      'PRO',
      [
        { element: 'C', x: 0.0, y: 1.5, z: 0.0 },   // Beta C
        { element: 'C', x: -1.1, y: 2.1, z: -0.7 },  // Gamma C
        { element: 'C', x: -2.1, y: 0.9, z: -0.6 }   // Delta C (connecting to amine N)
      ],
      [
        { from: 'C_ca', to: 'S_0', type: 'single' },
        { from: 'S_0', to: 'S_1', type: 'single' },
        { from: 'S_1', to: 'S_2', type: 'single' },
        { from: 'S_2', to: 'N_n', type: 'single' } // Close custom covalent ring!
      ]
    )
  },
  {
    id: 'SER',
    namePl: 'Seryna',
    nameEn: 'Serine',
    code3: 'Ser',
    code1: 'S',
    type: 'polar',
    formula: '-CH₂OH',
    mass: 105.1,
    description: 'Polarny hydroksyaminokwas. Częsty substrat do fosforylacji, regulujący funkcje życiowe białka.',
    color: '#0ea5e9',
    vanDerWaals: 89,
    hydrophobicity: -0.8,
    ...makeChemicalStructure(
      'SER',
      [
        { element: 'C', x: 0.0, y: 1.5, z: 0.0 },  // Beta C
        { element: 'O', x: 1.1, y: 2.2, z: -0.4 }  // Hydroxyl O
      ],
      [
        { from: 'C_ca', to: 'S_0', type: 'single' },
        { from: 'S_0', to: 'S_1', type: 'single' }
      ]
    )
  },
  {
    id: 'THR',
    namePl: 'Treonina',
    nameEn: 'Threonine',
    code3: 'Thr',
    code1: 'T',
    type: 'polar',
    formula: '-CH(OH)CH₃',
    mass: 119.1,
    description: 'Uczestniczy w potranslacyjnych modyfikacjach i buduje silne wiązania wodorowe.',
    color: '#2563eb',
    vanDerWaals: 116,
    hydrophobicity: -0.7,
    ...makeChemicalStructure(
      'THR',
      [
        { element: 'C', x: 0.0, y: 1.5, z: 0.0 },   // Beta C
        { element: 'O', x: 1.1, y: 2.2, z: -0.3 },  // O
        { element: 'C', x: -1.1, y: 2.4, z: -0.4 }  // C-methyl
      ],
      [
        { from: 'C_ca', to: 'S_0', type: 'single' },
        { from: 'S_0', to: 'S_1', type: 'single' },
        { from: 'S_0', to: 'S_2', type: 'single' }
      ]
    )
  },
  {
    id: 'TRP',
    namePl: 'Tryptofan',
    nameEn: 'Tryptophan',
    code3: 'Trp',
    code1: 'W',
    type: 'hydrophobic',
    formula: '-CH₂-C₈H₆N',
    mass: 204.2,
    description: 'Największy ze standardowych aminokwasów. Posiada potężny pierścień indolowy absorbujący fale światła UV.',
    color: '#b45309',
    vanDerWaals: 227,
    hydrophobicity: 1.9,
    ...makeChemicalStructure(
      'TRP',
      [
        { element: 'C', x: 0.0, y: 1.5, z: 0.0 },   // Beta-C
        { element: 'C', x: -0.4, y: 2.7, z: -0.7 },  // Indole C1
        { element: 'N', x: -1.7, y: 3.1, z: -0.5 },  // Indole N
        { element: 'C', x: -1.8, y: 4.2, z: -1.4 },  // Indole C2
        { element: 'C', x: -0.6, y: 4.6, z: -2.3 },  // Indole C3
        { element: 'C', x: -0.1, y: 5.8, z: -3.2 },  // Benzene C1
        { element: 'C', x: 1.1, y: 5.9, z: -3.9 },   // Benzene C2
        { element: 'C', x: 1.9, y: 4.8, z: -3.7 },   // Benzene C3
        { element: 'C', x: 1.4, y: 3.6, z: -2.8 }    // Benzene C4
      ],
      [
        { from: 'C_ca', to: 'S_0', type: 'single' },
        { from: 'S_0', to: 'S_1', type: 'single' },
        { from: 'S_1', to: 'S_2', type: 'single' },
        { from: 'S_2', to: 'S_3', type: 'double' },
        { from: 'S_3', to: 'S_4', type: 'single' },
        { from: 'S_4', to: 'S_5', type: 'double' },
        { from: 'S_5', to: 'S_6', type: 'single' },
        { from: 'S_6', to: 'S_7', type: 'double' },
        { from: 'S_7', to: 'S_8', type: 'single' },
        { from: 'S_8', to: 'S_1', type: 'double' }
      ]
    )
  },
  {
    id: 'TYR',
    namePl: 'Tyrozyna',
    nameEn: 'Tyrosine',
    code3: 'Tyr',
    code1: 'Y',
    type: 'polar',
    formula: '-CH₂-C₆H₄OH',
    mass: 181.2,
    description: 'Polarny aromatyczny aminokwas. Posiada czynną grupę hydroksylową na wierzchołku pierścienia, kluczową w enzymach.',
    color: '#38bdf8',
    vanDerWaals: 191,
    hydrophobicity: -1.3,
    ...makeChemicalStructure(
      'TYR',
      [
        { element: 'C', x: 0.0, y: 1.5, z: 0.0 },   // Beta C
        { element: 'C', x: -0.3, y: 2.7, z: -0.8 },  // Ring C1
        { element: 'C', x: -1.5, y: 3.4, z: -0.7 },  // Ring C2
        { element: 'C', x: -1.8, y: 4.5, z: -1.5 },  // Ring C3
        { element: 'C', x: -0.8, y: 4.9, z: -2.3 },  // Ring C4 (hydroxyl attach)
        { element: 'O', x: -1.1, y: 6.0, z: -3.1 },  // Phenol O
        { element: 'C', x: 0.4, y: 4.2, z: -2.4 },   // Ring C5
        { element: 'C', x: 0.7, y: 3.1, z: -1.6 }    // Ring C6
      ],
      [
        { from: 'C_ca', to: 'S_0', type: 'single' },
        { from: 'S_0', to: 'S_1', type: 'single' },
        { from: 'S_1', to: 'S_2', type: 'double' },
        { from: 'S_2', to: 'S_3', type: 'single' },
        { from: 'S_3', to: 'S_4', type: 'double' },
        { from: 'S_4', to: 'S_5', type: 'single' }, // Connect to Oxygen
        { from: 'S_4', to: 'S_6', type: 'single' },
        { from: 'S_6', to: 'S_7', type: 'double' },
        { from: 'S_7', to: 'S_1', type: 'single' }
      ]
    )
  },
  {
    id: 'VAL',
    namePl: 'Walina',
    nameEn: 'Valine',
    code3: 'Val',
    code1: 'V',
    type: 'hydrophobic',
    formula: '-CH(CH₃)₂',
    mass: 117.1,
    description: 'Rozgałęziony aminokwas niezbędny. Trwałość i hydrofobowość na wysokim poziomie.',
    color: '#ca8a04',
    vanDerWaals: 105,
    hydrophobicity: 4.2,
    ...makeChemicalStructure(
      'VAL',
      [
        { element: 'C', x: 0.0, y: 1.5, z: 0.0 },   // S_0 Beta Carbon
        { element: 'C', x: -1.1, y: 2.3, z: -0.7 },  // S_1 Gamma Carbon 1
        { element: 'C', x: 1.1, y: 2.3, z: 0.7 }    // S_2 Gamma Carbon 2
      ],
      [
        { from: 'C_ca', to: 'S_0', type: 'single' },
        { from: 'S_0', to: 'S_1', type: 'single' },
        { from: 'S_0', to: 'S_2', type: 'single' }
      ]
    )
  },
  {
    id: 'SEC',
    namePl: 'Selenocysteina',
    nameEn: 'Selenocysteine',
    code3: 'Sec',
    code1: 'U',
    type: 'special',
    formula: '-CH₂SeH',
    mass: 168.1,
    description: 'Niezwykle fascynujący 21. AMINOKWAS! Zawiera atom Selektywny (Se-Selen) zamiast siarki, o wysokim potencjale przeciwutleniającym.',
    color: '#10b981',
    vanDerWaals: 118,
    hydrophobicity: 1.2,
    ...makeChemicalStructure(
      'SEC',
      [
        { element: 'C', x: 0.0, y: 1.5, z: 0.0 },  // Beta Carbon
        { element: 'Se', x: 1.3, y: 2.8, z: -0.6 } // Selen atom!
      ],
      [
        { from: 'C_ca', to: 'S_0', type: 'single' },
        { from: 'S_0', to: 'S_1', type: 'single' }
      ]
    )
  }
];

// Sound Keys mapping specific chords for left & right marginal columns in English and Polish
export const SOUND_KEYS: SoundKey[] = [
  { digit: 0, frequency: 261.63, note: 'C4', chord: 'C-dur (C Major)', description: 'Pure, joyful triad of safety and warmth' },
  { digit: 1, frequency: 261.63, note: 'C4 (Min)', chord: 'c-moll (C Minor)', description: 'Grave, serious and sorrowful mool tone' },
  { digit: 2, frequency: 392.00, note: 'G4', chord: 'G-dur (G Major)', description: 'Bright, open and radiant primary dominant' },
  { digit: 3, frequency: 392.00, note: 'G4 (Min)', chord: 'g-moll (G Minor)', description: 'Melancholic, sentimental and poetic state' },
  { digit: 4, frequency: 293.66, note: 'D4', chord: 'D-dur (D Major)', description: 'Triumphant, clear, proud and uplifting' },
  { digit: 5, frequency: 293.66, note: 'D4 (Min)', chord: 'd-moll (D Minor)', description: 'Deep, dramatic and cinematic dark gravity' },
  { digit: 6, frequency: 440.00, note: 'A4', chord: 'A-dur (A Major)', description: 'Luminous, glowing and bright spiritual key' },
  { digit: 7, frequency: 440.00, note: 'A4 (Min)', chord: 'a-moll (A Minor)', description: 'Nostalgic, gentle and classical focus' },
  { digit: 8, frequency: 329.63, note: 'E4', chord: 'E-dur (E Major)', description: 'Energetic, powerful, lively and electric' },
  { digit: 9, frequency: 329.63, note: 'E4 (Min)', chord: 'e-moll (E Minor)', description: 'Unsettling, tense and movie-like suspense' },
  { digit: 'X', frequency: 466.16, note: 'Bb4', chord: 'B-dur (Bb Major)', description: 'Noble, royal, stable and comforting' },
  { digit: 'Y', frequency: 466.16, note: 'Bb4 (Min)', chord: 'b-moll (Bb Minor)', description: 'Dark, gloomy, mool and mysterious black' },
  { digit: 'Z', frequency: 349.23, note: 'F4', chord: 'F-dur (F Major)', description: 'Peaceful, calm, pastoral and organic breeze' },
  { digit: '-', frequency: 349.23, note: 'F4 (Min)', chord: 'f-moll (F Minor)', description: 'Passionate, deep, intense and emotional' }
];

// Coordinate helper function to generate 3D spatial coordinate for patch IDs [1-248]
// satisfying user specified conditions:
// - patch 1: absolute top (0, 1, 0)
// - patch 248: absolute bottom (0, -1, 0)
// - patch 50: positive X direction (1, 0, 0)
// - patch 100: negative X direction (-1, 0, 0)
// - patch 150: positive Z direction (0, 0, 1)
// - patch 200: negative Z direction (0, 0, -1)
export function getPatchCoordinate(patchId: number): { x: number; y: number; z: number } {
  // Let's protect bounds
  const id = Math.min(248, Math.max(1, Math.round(patchId)));

  if (id === 1) return { x: 0, y: 1, z: 0 };       // Góra (+Y - Port wejściowy N)
  if (id === 248) return { x: 0, y: -1, z: 0 };     // Dół (-Y - Port wyjściowy C)
  if (id === 50) return { x: 0, y: 0, z: 1 };       // Przód (+Z)
  if (id === 100) return { x: 1, y: 0, z: 0 };      // Prawo (+X)
  if (id === 150) return { x: 0, y: 0, z: -1 };     // Tył (-Z)
  if (id === 200) return { x: -1, y: 0, z: 0 };     // Lewo (-X)

  // For other patch IDs, let's distribute them deterministically on a sphere
  // using a Fibonacci spiral, but with slight radial offsets or distortions
  // to perfectly preserve smooth spherical interpolation from top to bottom
  const count = 248;
  const i = id - 0.5;
  const phi = Math.acos(1 - 2 * i / count); // latitude [0, PI]
  const theta = Math.PI * (1 + 5 ** 0.5) * i; // longitude based on golden ratio

  return {
    x: Math.sin(phi) * Math.cos(theta),
    y: Math.sin(phi) * Math.sin(theta),
    z: Math.cos(phi)
  };
}
