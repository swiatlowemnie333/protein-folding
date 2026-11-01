// SNAPLOCK Universal Loop Closure & Adaptive Disulfide Restraint Solver
// Implements:
// 1. Universal Disulfide Topology Inference & Explicit SSBOND Hard Boundary Constraints
// 2. Two-Way Inverse Kinematics (Bidirectional Cyclic Coordinate Descent - CCD IK) on flexible hinge residues
// 3. Adaptive Exponential Spring Constant: E_ssbond = w_ss * max(0, d(Sg, Sg) - d_target)^2
//    with automatic exponential growth w_ss^(k) = w_ss^(0) * exp(beta * k) when initial deviation > 2.5 Å
// 4. Stereochemical CB & SG Reconstruction (CA-CB = 1.53 Å, CB-SG = 1.81 Å, target SG-SG = 2.05 Å, CA-CA = 5.50 Å)
// 5. Drift-Compensated Micro-Facet Angular Quantization (Level 2: 2232 facets, < 1.0°/res without breaking SSBOND < 0.5 Å)

import * as THREE from 'three';
import { SNAPLOCK_21_BLOCKS, generate248SphericalFacets } from './snaplock21Blocks';

export type SolverPipelineMode = 'STANDARD_L2_PIPELINE' | 'TWO_WAY_IK_ADAPTIVE_SSBOND';

export interface DisulfideBridgeDiagnostic {
  bridgeIndex: number;
  idx1: number;                   // 0-indexed residue A
  idx2: number;                   // 0-indexed residue B
  res1: number;                   // 1-indexed PDB residue A
  res2: number;                   // 1-indexed PDB residue B
  initialCaDist: number;          // [Å] Pre-closure Cα-Cα distance
  finalCaDist: number;            // [Å] Post-closure Cα-Cα distance (target 5.50 Å)
  caDeviation: number;            // [Å] |finalCaDist - 5.50|
  initialSgDist: number;          // [Å] Pre-closure Sγ-Sγ distance
  finalSgDist: number;            // [Å] Post-closure Sγ-Sγ distance (target 2.05 Å)
  sgDeviation: number;            // [Å] |finalSgDist - 2.05|
  energySsbond: number;           // [kcal/mol] w_ss * max(0, d(Sg,Sg) - 2.05)^2
  hingeResidues: number[];        // 1-indexed flexible hinge residues used in CCD IK
  passedSubHalfAngstrom: boolean; // Both Sγ-Sγ and Cα-Cα deviations < 0.5 Å
}

export interface AdaptiveIterationTrace {
  iteration: number;
  wSs: number;                    // Adaptive spring weight w_ss^(k)
  maxCaDeviation: number;         // [Å]
  maxSgDeviation: number;         // [Å]
  totalSsbondEnergy: number;      // [kcal/mol] sum of w_ss * max(0, d_sg - 2.05)^2
  meanAngularErrorDeg: number;    // [deg]
  vdwClashes: number;
  dMaxActual: number;             // [Å]
}

export interface UniversalLoopClosureResult {
  pipelineMode: SolverPipelineMode;
  positions: [number, number, number][];      // Closed & quantized Cα coordinates [Å]
  continuousClosedPositions: [number, number, number][]; // Continuous IK-closed Cα target skeleton [Å]
  cbPositions: ([number, number, number] | null)[]; // Reconstructed Cβ coordinates [Å]
  sgPositions: ([number, number, number] | null)[]; // Reconstructed Sγ coordinates [Å]
  ssPairs: [number, number][];
  bridgeDiagnostics: DisulfideBridgeDiagnostic[];
  initialMaxDeviationAngstrom: number;
  finalMaxSgDeviationAngstrom: number;
  finalMaxCaDeviationAngstrom: number;
  exponentialBoostTriggered: boolean;         // True if initial deviation > 2.5 Å
  initialSpringWeight: number;
  finalSpringWeight: number;
  totalSsbondPenaltyEnergy: number;           // Exact sum of w_ss * max(0, d(Sg,Sg) - d_target)^2
  meanMicroFacetErrorDeg: number;             // Global angular error on micro-facets (< 1.0° on L2)
  maxMicroFacetErrorDeg: number;
  perBondPatchIds: number[];
  perBondNormals: [number, number, number][];
  perBondAngularErrorsDeg: number[];
  vdwClashesCount: number;
  iterationTraces: AdaptiveIterationTrace[];
}

export const TARGET_SG_SG_DIST = 2.05; // [Å] Covalent disulfide Sγ-Sγ bond length
export const TARGET_CA_CA_DIST = 5.50; // [Å] Canonical Cα-Cα disulfide bridge span
export const PEPTIDE_BOND_CA_DIST = 3.80; // [Å] Standard trans-peptide Cα-Cα distance
export const EXPONENTIAL_DEVIATION_THRESHOLD = 2.5; // [Å] Threshold triggering exponential w_ss growth

// -----------------------------------------------------------------------------
// 1. PARSE EXPLICIT SSBOND RECORDS OR INFER UNIVERSAL BIOPHYSICAL TOPOLOGY
// -----------------------------------------------------------------------------

/**
 * Extracts explicit SSBOND pairs from PDB or FASTA header remarks if present.
 * Supports standard PDB `SSBOND` lines as well as inline tags like `SSBOND: 3-40, 4-32, 16-26`.
 */
export function parseExplicitSsbondFromText(rawText: string, maxLength: number): [number, number][] {
  if (!rawText || maxLength < 2) return [];
  const pairs: [number, number][] = [];
  const used = new Set<number>();

  const lines = rawText.split('\n');
  for (const line of lines) {
    const trimmed = line.trim();
    if (trimmed.startsWith('SSBOND')) {
      // Standard PDB SSBOND format or free-form numbers
      const nums = trimmed
        .replace(/^SSBOND\s+\d+/, '')
        .match(/\b\d+\b/g);
      if (nums && nums.length >= 2) {
        const r1 = parseInt(nums[0], 10) - 1;
        const r2 = parseInt(nums[1], 10) - 1;
        const a = Math.min(r1, r2);
        const b = Math.max(r1, r2);
        if (a >= 0 && b < maxLength && a < b && !used.has(a) && !used.has(b)) {
          pairs.push([a, b]);
          used.add(a);
          used.add(b);
        }
      }
    } else if (/SSBOND\s*[:=]/i.test(trimmed)) {
      const matchPart = trimmed.split(/SSBOND\s*[:=]/i)[1] || '';
      const pairTokens = matchPart.match(/(\d+)\s*[-–:]\s*(\d+)/g) || [];
      for (const tok of pairTokens) {
        const parts = tok.split(/[-–:]/).map(s => parseInt(s.trim(), 10) - 1);
        if (parts.length === 2 && !isNaN(parts[0]) && !isNaN(parts[1])) {
          const a = Math.min(parts[0], parts[1]);
          const b = Math.max(parts[0], parts[1]);
          if (a >= 0 && b < maxLength && a < b && !used.has(a) && !used.has(b)) {
            pairs.push([a, b]);
            used.add(a);
            used.add(b);
          }
        }
      }
    }
  }

  return pairs.sort((p, q) => p[0] - q[0]);
}

/**
 * Universal biophysical disulfide topology matcher (zero hardcoded sequence strings).
 * Evaluates steric vicinal barriers (|j - i| >= 4), nested globular β/α cradles,
 * canonical CC-{3,4}-C helical staples, and Inhibitor Cystine Knot (ICK) motifs.
 */
export function inferUniversalDisulfideTopology(
  sequenceCodes: string[],
  explicitPairs?: [number, number][]
): [number, number][] {
  const n = sequenceCodes?.length ?? 0;
  if (n < 2) return [];

  if (explicitPairs && explicitPairs.length > 0) {
    const used = new Set<number>();
    const valid: [number, number][] = [];
    for (const pair of explicitPairs) {
      if (!pair || pair.length < 2) continue;
      const a = Math.min(pair[0], pair[1]);
      const b = Math.max(pair[0], pair[1]);
      if (a >= 0 && b < n && a < b && !used.has(a) && !used.has(b)) {
        valid.push([a, b]);
        used.add(a);
        used.add(b);
      }
    }
    if (valid.length > 0) return valid.sort((x, y) => x[0] - y[0]);
  }

  const cysIndices: number[] = [];
  for (let i = 0; i < n; i++) {
    const code = sequenceCodes[i];
    if (code === 'CYS' || code === 'SEC') {
      cysIndices.push(i);
    }
  }

  const m = cysIndices.length;
  if (m < 2) return [];

  // If odd number of cysteines > 2, pick the (m - 1) subset with optimal pairing
  if (m <= 10) {
    const candidateLists: number[][] = [];
    if (m % 2 === 0) {
      candidateLists.push(cysIndices);
    } else {
      for (let skip = 0; skip < m; skip++) {
        candidateLists.push(cysIndices.filter((_, idx) => idx !== skip));
      }
    }

    let bestPairs: [number, number][] = [];
    let bestScore = Infinity;

    const scoreMatching = (pairs: [number, number][], activeCys: number[]): number => {
      let score = 0;
      const spans = pairs.map(([a, b]) => b - a);
      const meanSpan = spans.reduce((acc, s) => acc + s, 0) / Math.max(1, spans.length);

      // Detect if there is an internal (non-N-terminal) vicinal CC doublet (characteristic of ICK toxins)
      let hasInternalVicinalCC = false;
      for (let k = 1; k < activeCys.length - 2; k++) {
        if (activeCys[k + 1] - activeCys[k] === 1) {
          hasInternalVicinalCC = true;
        }
      }

      for (const [a, b] of pairs) {
        const sep = b - a;
        if (sep < 4) {
          // Strong steric barrier against strained vicinal disulfides
          score += 5000 * (4 - sep) * (4 - sep);
        }
        // Loop entropy cost
        score += 0.25 * Math.log(sep + 1);
        if (hasInternalVicinalCC) {
          score += 0.08 * (sep - meanSpan) * (sep - meanSpan);
        }
      }

      for (let i = 0; i < pairs.length; i++) {
        for (let j = i + 1; j < pairs.length; j++) {
          const [a, b] = pairs[i];
          const [c, d] = pairs[j];
          const isNested = a < c && d < b;
          const isCrossing = a < c && c < b && b < d;

          if (isNested) {
            // Nested globular cradle stabilization (e.g., Crambin CYS3-40 > CYS4-32 > CYS16-26)
            score -= hasInternalVicinalCC ? 1.5 : 8.5;
          } else if (isCrossing) {
            // Check canonical CC-{3,4}-C helical N-cap staple (e.g., CYS(a), CYS(a+1) with CYS(a+5/6))
            const isHelicalCCStaple = c - a === 1 && (b - a === 5 || b - a === 6);
            if (isHelicalCCStaple) {
              score -= 24.0;
            } else if (hasInternalVicinalCC) {
              // Canonical Inhibitor Cystine Knot (ICK) [I-IV, II-V, III-VI] topology stabilization
              score -= 9.5;
            } else {
              // Topological knot penalty in standard globular folds
              score += 6.5;
            }
          }
        }
      }

      return score;
    };

    const enumerateMatchings = (remaining: number[], acc: [number, number][], fullList: number[]) => {
      if (remaining.length === 0) {
        const sc = scoreMatching(acc, fullList);
        if (sc < bestScore) {
          bestScore = sc;
          bestPairs = acc.map(p => [p[0], p[1]]);
        }
        return;
      }
      const first = remaining[0];
      for (let k = 1; k < remaining.length; k++) {
        const partner = remaining[k];
        const nextRem = remaining.filter((_, idx) => idx !== 0 && idx !== k);
        acc.push([first, partner]);
        enumerateMatchings(nextRem, acc, fullList);
        acc.pop();
      }
    };

    for (const subList of candidateLists) {
      enumerateMatchings(subList, [], subList);
    }

    return bestPairs.sort((x, y) => x[0] - y[0]);
  }

  // For very large cysteine-rich proteins (>10 CYS), pair non-vicinal cradle intervals greedily
  const used = new Set<number>();
  const largePairs: [number, number][] = [];
  for (let i = 0; i < m; i++) {
    const c1 = cysIndices[i];
    if (used.has(c1)) continue;
    let chosen = -1;
    for (let j = i + 1; j < m; j++) {
      const c2 = cysIndices[j];
      if (!used.has(c2) && c2 - c1 >= 4) {
        chosen = c2;
        break;
      }
    }
    if (chosen !== -1) {
      used.add(c1);
      used.add(chosen);
      largePairs.push([c1, chosen]);
    }
  }
  return largePairs.sort((x, y) => x[0] - y[0]);
}

// -----------------------------------------------------------------------------
// 2. KINEMATIC HINGE FLEXIBILITY FROM 22-PARAMETER SNAPLOCK MATRIX
// -----------------------------------------------------------------------------

export function getResidueKinematicFlexibility(
  code3: string,
  secondaryType?: 'HELIX' | 'SHEET' | 'TURN' | 'COIL'
): number {
  const block = SNAPLOCK_21_BLOCKS[code3] || SNAPLOCK_21_BLOCKS['ALA'];
  const flexLabel = block.physico.chainFlexibility;

  let base = 0.55;
  if (code3 === 'PRO' || flexLabel.includes('Zablokowana')) {
    base = 0.08;
  } else if (flexLabel.includes('Najniższa')) {
    base = 0.20;
  } else if (flexLabel.includes('Niska')) {
    base = 0.35;
  } else if (flexLabel.includes('Sztywna') || flexLabel.includes('Usztywniona')) {
    base = 0.45;
  } else if (flexLabel.includes('Średnia')) {
    base = 0.65;
  } else if (flexLabel.includes('Bardzo wysoka')) {
    base = 1.00;
  } else if (flexLabel.includes('Wysoka')) {
    base = 0.85;
  }

  if (secondaryType === 'TURN' || secondaryType === 'COIL') {
    base = Math.min(1.0, base * 1.28);
  } else if (secondaryType === 'HELIX' || secondaryType === 'SHEET') {
    base *= 0.78;
  }

  return base;
}

// -----------------------------------------------------------------------------
// 3. STEREOCHEMICAL CB & SG RECONSTRUCTION & DISULFIDE ENERGY EVALUATION
// -----------------------------------------------------------------------------

/**
 * Reconstructs realistic 3D Cβ (1.53 Å from Cα) and Sγ (1.81 Å from Cβ, ∠CA-CB-SG = 114°)
 * atomic coordinates for all residues, optimizing χ1 rotamers across SSBOND pairs so
 * d(Sγ, Sγ) approaches 2.05 Å whenever Cα-Cα is within kinematic reach.
 */
export function reconstructSidechainCbAndSg(
  sequenceCodes: string[],
  caPositions: [number, number, number][],
  ssPairs: [number, number][]
): {
  cbPositions: ([number, number, number] | null)[];
  sgPositions: ([number, number, number] | null)[];
} {
  const n = Math.min(sequenceCodes?.length ?? 0, caPositions?.length ?? 0);
  const cbPositions: ([number, number, number] | null)[] = new Array(n).fill(null);
  const sgPositions: ([number, number, number] | null)[] = new Array(n).fill(null);
  if (n === 0) return { cbPositions, sgPositions };

  const partnerMap = new Map<number, number>();
  for (const [i1, i2] of ssPairs) {
    if (i1 >= 0 && i2 < n) {
      partnerMap.set(i1, i2);
      partnerMap.set(i2, i1);
    }
  }

  const CA_CB_LEN = 1.53;
  const CB_SG_LEN = 1.81;
  const CONE_ANGLE_RAD = (66.0 * Math.PI) / 180.0; // Supplement of 114° tetrahedral angle
  const cosCone = Math.cos(CONE_ANGLE_RAD);
  const sinCone = Math.sin(CONE_ANGLE_RAD);

  // First compute CB for all non-GLY residues
  const cbAxisVecs: (THREE.Vector3 | null)[] = new Array(n).fill(null);

  for (let i = 0; i < n; i++) {
    const code = sequenceCodes[i] || 'ALA';
    if (code === 'GLY') continue;

    const ca = new THREE.Vector3(
      caPositions[i]?.[0] ?? 0,
      caPositions[i]?.[1] ?? 0,
      caPositions[i]?.[2] ?? 0
    );
    const prev = i > 0
      ? new THREE.Vector3(caPositions[i - 1]?.[0] ?? 0, caPositions[i - 1]?.[1] ?? 0, caPositions[i - 1]?.[2] ?? 0)
      : ca.clone().add(new THREE.Vector3(0, -1, 0));
    const next = i < n - 1
      ? new THREE.Vector3(caPositions[i + 1]?.[0] ?? 0, caPositions[i + 1]?.[1] ?? 0, caPositions[i + 1]?.[2] ?? 0)
      : ca.clone().add(new THREE.Vector3(0, 1, 0));

    const tangent = next.clone().sub(prev);
    if (tangent.lengthSq() < 1e-6) tangent.set(0, 1, 0);
    tangent.normalize();

    let cbDir: THREE.Vector3;
    const partnerIdx = partnerMap.get(i);
    if (partnerIdx !== undefined && caPositions[partnerIdx]) {
      const partnerCa = new THREE.Vector3(
        caPositions[partnerIdx]?.[0] ?? 0,
        caPositions[partnerIdx]?.[1] ?? 0,
        caPositions[partnerIdx]?.[2] ?? 0
      );
      const bridgeAxis = partnerCa.clone().sub(ca);
      if (bridgeAxis.lengthSq() < 1e-6) bridgeAxis.set(1, 0, 0);
      bridgeAxis.normalize();

      // Build perpendicular basis to bridgeAxis for tetrahedral 35° tilt
      const refUp = Math.abs(bridgeAxis.y) < 0.9 ? new THREE.Vector3(0, 1, 0) : new THREE.Vector3(1, 0, 0);
      const perp = new THREE.Vector3().crossVectors(bridgeAxis, refUp).normalize();
      const tiltRad = (34.0 * Math.PI) / 180.0;
      cbDir = bridgeAxis.clone().multiplyScalar(Math.cos(tiltRad)).addScaledVector(perp, Math.sin(tiltRad)).normalize();
    } else {
      const bisector = prev.clone().sub(ca).normalize().add(next.clone().sub(ca).normalize());
      if (bisector.lengthSq() > 1e-4) {
        cbDir = bisector.multiplyScalar(-1).normalize();
      } else {
        const ref = Math.abs(tangent.y) < 0.9 ? new THREE.Vector3(0, 1, 0) : new THREE.Vector3(1, 0, 0);
        cbDir = new THREE.Vector3().crossVectors(tangent, ref).normalize();
      }
    }

    cbAxisVecs[i] = cbDir;
    const cb = ca.clone().addScaledVector(cbDir, CA_CB_LEN);
    cbPositions[i] = [
      parseFloat(cb.x.toFixed(3)),
      parseFloat(cb.y.toFixed(3)),
      parseFloat(cb.z.toFixed(3))
    ];
  }

  // Reconstruct SG for paired cysteines via discrete χ1 rotamer optimization
  const handled = new Set<number>();
  for (const [i1, i2] of ssPairs) {
    if (i1 < 0 || i2 >= n || !cbPositions[i1] || !cbPositions[i2] || !cbAxisVecs[i1] || !cbAxisVecs[i2]) continue;
    handled.add(i1);
    handled.add(i2);

    const cb1 = new THREE.Vector3(...cbPositions[i1]!);
    const cb2 = new THREE.Vector3(...cbPositions[i2]!);
    const axis1 = cbAxisVecs[i1]!;
    const axis2 = cbAxisVecs[i2]!;

    const buildOrtho = (ax: THREE.Vector3): [THREE.Vector3, THREE.Vector3] => {
      const ref = Math.abs(ax.y) < 0.9 ? new THREE.Vector3(0, 1, 0) : new THREE.Vector3(1, 0, 0);
      const u = new THREE.Vector3().crossVectors(ax, ref).normalize();
      const v = new THREE.Vector3().crossVectors(ax, u).normalize();
      return [u, v];
    };

    const [u1, v1] = buildOrtho(axis1);
    const [u2, v2] = buildOrtho(axis2);

    const STEPS = 24;
    let bestSg1 = cb1.clone().addScaledVector(axis1, CB_SG_LEN);
    let bestSg2 = cb2.clone().addScaledVector(axis2, CB_SG_LEN);
    let bestDiff = Infinity;

    for (let s1 = 0; s1 < STEPS; s1++) {
      const chi1 = (s1 * 2 * Math.PI) / STEPS;
      const radial1 = u1.clone().multiplyScalar(Math.cos(chi1)).addScaledVector(v1, Math.sin(chi1));
      const cand1 = cb1
        .clone()
        .addScaledVector(axis1, CB_SG_LEN * cosCone)
        .addScaledVector(radial1, CB_SG_LEN * sinCone);

      for (let s2 = 0; s2 < STEPS; s2++) {
        const chi2 = (s2 * 2 * Math.PI) / STEPS;
        const radial2 = u2.clone().multiplyScalar(Math.cos(chi2)).addScaledVector(v2, Math.sin(chi2));
        const cand2 = cb2
          .clone()
          .addScaledVector(axis2, CB_SG_LEN * cosCone)
          .addScaledVector(radial2, CB_SG_LEN * sinCone);

        const d = cand1.distanceTo(cand2);
        const diff = Math.abs(d - TARGET_SG_SG_DIST);
        if (diff < bestDiff) {
          bestDiff = diff;
          bestSg1 = cand1;
          bestSg2 = cand2;
        }
      }
    }

    sgPositions[i1] = [
      parseFloat(bestSg1.x.toFixed(3)),
      parseFloat(bestSg1.y.toFixed(3)),
      parseFloat(bestSg1.z.toFixed(3))
    ];
    sgPositions[i2] = [
      parseFloat(bestSg2.x.toFixed(3)),
      parseFloat(bestSg2.y.toFixed(3)),
      parseFloat(bestSg2.z.toFixed(3))
    ];
  }

  // Unpaired CYS / SEC free thiols
  for (let i = 0; i < n; i++) {
    const code = sequenceCodes[i];
    if ((code === 'CYS' || code === 'SEC') && !handled.has(i) && cbPositions[i] && cbAxisVecs[i]) {
      const cb = new THREE.Vector3(...cbPositions[i]!);
      const ax = cbAxisVecs[i]!;
      const ref = Math.abs(ax.y) < 0.9 ? new THREE.Vector3(0, 1, 0) : new THREE.Vector3(1, 0, 0);
      const u = new THREE.Vector3().crossVectors(ax, ref).normalize();
      const sg = cb
        .clone()
        .addScaledVector(ax, CB_SG_LEN * cosCone)
        .addScaledVector(u, CB_SG_LEN * sinCone);
      sgPositions[i] = [
        parseFloat(sg.x.toFixed(3)),
        parseFloat(sg.y.toFixed(3)),
        parseFloat(sg.z.toFixed(3))
      ];
    }
  }

  return { cbPositions, sgPositions };
}

/**
 * Evaluates the exact adaptive harmonic disulfide penalty formula:
 * E_ssbond = w_ss * max(0, d(Sγ, Sγ) - d_target)^2
 */
export function computeAdaptiveSsbondEnergy(
  dSgSg: number,
  wSs: number,
  dTarget: number = TARGET_SG_SG_DIST
): number {
  const stretch = Math.max(0, dSgSg - dTarget);
  return wSs * stretch * stretch;
}

// -----------------------------------------------------------------------------
// 4. DENSE SUBDIVIDED MICRO-FACET + TORSION RING QUANTIZER
// -----------------------------------------------------------------------------

export interface CachedPolarMicroLattice {
  level: 0 | 1 | 2;
  facetCount: number;
  normals: [number, number, number][];
  polarThetas: number[];       // Polar angle from +Y [0, PI] for each facet
  torsionSteps: number;        // Discrete azimuthal contact steps around #248 -> #1 axis
}

const POLAR_LATTICE_CACHE: Partial<Record<0 | 1 | 2, CachedPolarMicroLattice>> = {};

function subdivideTri1to3(
  tri: [[number, number, number], [number, number, number], [number, number, number]]
): [[number, number, number], [number, number, number], [number, number, number]][] {
  const [a, b, c] = tri;
  const cx = a[0] + b[0] + c[0];
  const cy = a[1] + b[1] + c[1];
  const cz = a[2] + b[2] + c[2];
  const len = Math.sqrt(cx * cx + cy * cy + cz * cz) || 1.0;
  const m: [number, number, number] = [cx / len, cy / len, cz / len];
  return [
    [a, b, m],
    [b, c, m],
    [c, a, m]
  ];
}

function triNormal(
  tri: [[number, number, number], [number, number, number], [number, number, number]]
): [number, number, number] {
  const [a, b, c] = tri;
  const x = a[0] + b[0] + c[0];
  const y = a[1] + b[1] + c[1];
  const z = a[2] + b[2] + c[2];
  const len = Math.sqrt(x * x + y * y + z * z) || 1.0;
  return [x / len, y / len, z / len];
}

export function getPolarMicroLattice(level: 0 | 1 | 2): CachedPolarMicroLattice {
  if (POLAR_LATTICE_CACHE[level]) return POLAR_LATTICE_CACHE[level]!;

  const base = generate248SphericalFacets();
  let tris = base.vertices;
  let norms = base.normals;

  if (level >= 1) {
    const t1: typeof tris = [];
    const n1: typeof norms = [];
    for (const t of tris) {
      for (const st of subdivideTri1to3(t)) {
        t1.push(st);
        n1.push(triNormal(st));
      }
    }
    tris = t1;
    norms = n1;
  }

  if (level >= 2) {
    const t2: typeof tris = [];
    const n2: typeof norms = [];
    for (const t of tris) {
      for (const st of subdivideTri1to3(t)) {
        t2.push(st);
        n2.push(triNormal(st));
      }
    }
    tris = t2;
    norms = n2;
  }

  const polarThetas = norms.map(n => Math.acos(Math.max(-1, Math.min(1, n[1]))));
  const torsionSteps = level === 2 ? 360 : level === 1 ? 96 : 32;

  const built: CachedPolarMicroLattice = {
    level,
    facetCount: norms.length,
    normals: norms,
    polarThetas,
    torsionSteps
  };
  POLAR_LATTICE_CACHE[level] = built;
  return built;
}

/**
 * Snaps a target bond vector onto the discrete SNAPLOCK sphere state
 * (micro-facet polar index + discrete torsion jump around contact axis #248 -> #1).
 */
export function snapVectorToDiscreteMicroState(
  targetVec: THREE.Vector3,
  incomingAxis: THREE.Vector3,
  lattice: CachedPolarMicroLattice
): {
  patchId: number;
  snappedUnitDir: THREE.Vector3;
  angleErrorDeg: number;
} {
  const targetDir = targetVec.clone();
  if (targetDir.lengthSq() < 1e-8) {
    return {
      patchId: 1,
      snappedUnitDir: new THREE.Vector3(0, 1, 0),
      angleErrorDeg: 0
    };
  }
  targetDir.normalize();

  const axis = incomingAxis.clone();
  if (axis.lengthSq() < 1e-8) axis.set(0, 1, 0);
  axis.normalize();

  // Build local frame where +Y aligns with incomingAxis
  const frameQuat = new THREE.Quaternion().setFromUnitVectors(new THREE.Vector3(0, 1, 0), axis);
  const invFrameQuat = frameQuat.clone().invert();
  const localTarget = targetDir.clone().applyQuaternion(invFrameQuat);

  const targetTheta = Math.acos(Math.max(-1, Math.min(1, localTarget.y)));
  const targetPhi = Math.atan2(localTarget.z, localTarget.x);

  // 1. Find micro-facet with closest polar angle theta (with secondary preference for azimuth)
  let bestIdx = 0;
  let bestPolarDiff = Infinity;
  for (let i = 0; i < lattice.polarThetas.length; i++) {
    const diff = Math.abs(lattice.polarThetas[i] - targetTheta);
    if (diff < bestPolarDiff) {
      bestPolarDiff = diff;
      bestIdx = i;
    }
  }

  const chosenNormal = lattice.normals[bestIdx];
  const chosenPhi = Math.atan2(chosenNormal[2], chosenNormal[0]);

  // 2. Quantize discrete torsion jump around contact axis (+Y)
  // Note: In right-handed (X, Y, Z), rotation around +Y by +alpha maps atan2(z, x) -> phi - alpha,
  // so rotating by -quantizedDeltaPhi aligns localSnapped azimuth from chosenPhi to targetPhi.
  const stepRad = (2 * Math.PI) / lattice.torsionSteps;
  let rawDeltaPhi = targetPhi - chosenPhi;
  while (rawDeltaPhi > Math.PI) rawDeltaPhi -= 2 * Math.PI;
  while (rawDeltaPhi < -Math.PI) rawDeltaPhi += 2 * Math.PI;
  const discreteK = Math.round(rawDeltaPhi / stepRad);
  const quantizedDeltaPhi = discreteK * stepRad;

  const localSnapped = new THREE.Vector3(chosenNormal[0], chosenNormal[1], chosenNormal[2]);
  localSnapped.applyAxisAngle(new THREE.Vector3(0, 1, 0), -quantizedDeltaPhi);

  const worldSnapped = localSnapped.applyQuaternion(frameQuat).normalize();
  const dot = Math.max(-1, Math.min(1, worldSnapped.dot(targetDir)));
  const errDeg = (Math.acos(dot) * 180.0) / Math.PI;

  return {
    patchId: bestIdx + 1,
    snappedUnitDir: worldSnapped,
    angleErrorDeg: parseFloat(errDeg.toFixed(3))
  };
}

// -----------------------------------------------------------------------------
// 5. MAIN TWO-WAY IK LOOP CLOSURE + ADAPTIVE EXPONENTIAL SPRING SOLVER
// -----------------------------------------------------------------------------

export interface LoopClosureSolverOptions {
  sequenceCodes: string[];
  initialCaPositions: [number, number, number][];
  secondaryTypes?: ('HELIX' | 'SHEET' | 'TURN' | 'COIL')[];
  explicitDisulfidePairs?: [number, number][];
  subdivisionLevel?: 0 | 1 | 2;
  dMaxTarget?: number;
  enforceDMax?: boolean;
  baseSpringWeight?: number;   // Initial w_ss^(0)
  kappaBend?: number;          // Backbone bending stiffness
  maxIterations?: number;
  centerAtOrigin?: boolean;
}

export function runUniversalLoopClosureSolver(
  options: LoopClosureSolverOptions
): UniversalLoopClosureResult {
  const {
    sequenceCodes,
    initialCaPositions,
    secondaryTypes = [],
    explicitDisulfidePairs,
    subdivisionLevel = 2,
    dMaxTarget = 0,
    enforceDMax = false,
    baseSpringWeight = 3.0,
    kappaBend = 42.5,
    maxIterations = 24,
    centerAtOrigin = true
  } = options;

  const n = Math.min(sequenceCodes?.length ?? 0, initialCaPositions?.length ?? 0);
  if (n === 0) {
    return {
      pipelineMode: 'STANDARD_L2_PIPELINE',
      positions: [],
      continuousClosedPositions: [],
      cbPositions: [],
      sgPositions: [],
      ssPairs: [],
      bridgeDiagnostics: [],
      initialMaxDeviationAngstrom: 0,
      finalMaxSgDeviationAngstrom: 0,
      finalMaxCaDeviationAngstrom: 0,
      exponentialBoostTriggered: false,
      initialSpringWeight: baseSpringWeight,
      finalSpringWeight: baseSpringWeight,
      totalSsbondPenaltyEnergy: 0,
      meanMicroFacetErrorDeg: 0,
      maxMicroFacetErrorDeg: 0,
      perBondPatchIds: [],
      perBondNormals: [],
      perBondAngularErrorsDeg: [],
      vdwClashesCount: 0,
      iterationTraces: []
    };
  }

  // Clone input positions safely
  const pts: THREE.Vector3[] = [];
  for (let i = 0; i < n; i++) {
    const p = initialCaPositions[i];
    pts.push(new THREE.Vector3(p?.[0] ?? 0, p?.[1] ?? 0, p?.[2] ?? 0));
  }

  // Reference peptide bond lengths (clamped around canonical 3.80 Å)
  const bondLengths: number[] = [];
  for (let i = 0; i < n - 1; i++) {
    const rawLen = pts[i].distanceTo(pts[i + 1]);
    bondLengths.push(
      rawLen > 2.8 && rawLen < 4.4 ? rawLen : PEPTIDE_BOND_CA_DIST
    );
  }

  // Residue hinge compliance weights [0.08 .. 1.0]
  const flexWeights: number[] = [];
  for (let i = 0; i < n; i++) {
    flexWeights.push(getResidueKinematicFlexibility(sequenceCodes[i] || 'ALA', secondaryTypes[i]));
  }

  // Infer or validate SSBOND pairs
  const ssPairs = inferUniversalDisulfideTopology(sequenceCodes, explicitDisulfidePairs);
  const hasDisulfides = ssPairs.length > 0;
  const pipelineMode: SolverPipelineMode = hasDisulfides
    ? 'TWO_WAY_IK_ADAPTIVE_SSBOND'
    : 'STANDARD_L2_PIPELINE';

  // Measure initial pre-closure disulfide distances & deviations
  const rawInitialCoords: [number, number, number][] = pts.map(v => [v.x, v.y, v.z]);
  const initialSidechains = reconstructSidechainCbAndSg(sequenceCodes, rawInitialCoords, ssPairs);

  const initialBridgeMetrics = ssPairs.map(([i1, i2]) => {
    const caDist = pts[i1].distanceTo(pts[i2]);
    const sg1 = initialSidechains.sgPositions[i1];
    const sg2 = initialSidechains.sgPositions[i2];
    const sgDist =
      sg1 && sg2
        ? Math.sqrt((sg1[0] - sg2[0]) ** 2 + (sg1[1] - sg2[1]) ** 2 + (sg1[2] - sg2[2]) ** 2)
        : Math.max(TARGET_SG_SG_DIST, caDist - 3.45);
    return {
      caDist,
      sgDist,
      caDev: Math.abs(caDist - TARGET_CA_CA_DIST),
      sgDev: Math.abs(sgDist - TARGET_SG_SG_DIST)
    };
  });

  let initialMaxDev = 0;
  for (const m of initialBridgeMetrics) {
    initialMaxDev = Math.max(initialMaxDev, m.caDev, m.sgDev);
  }

  // Trigger automatic exponential w_ss growth if initial deviation > 2.5 Å
  const exponentialBoostTriggered = hasDisulfides && initialMaxDev > EXPONENTIAL_DEVIATION_THRESHOLD;

  const anchorSet = new Set<number>();
  for (const [i1, i2] of ssPairs) {
    anchorSet.add(i1);
    anchorSet.add(i2);
  }

  const lattice = getPolarMicroLattice(subdivisionLevel);
  const iterationTraces: AdaptiveIterationTrace[] = [];
  let currentWSs = Math.max(0.5, baseSpringWeight);

  // Helper: center points around [0, 0, 0]
  const recenterPoints = () => {
    const com = new THREE.Vector3(0, 0, 0);
    for (let i = 0; i < n; i++) com.add(pts[i]);
    com.divideScalar(Math.max(1, n));
    for (let i = 0; i < n; i++) pts[i].sub(com);
  };

  // Helper: compute D_max
  const computeCurrentDMax = (): number => {
    let maxD = 0;
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        const d = pts[i].distanceTo(pts[j]);
        if (d > maxD) maxD = d;
      }
    }
    return maxD;
  };

  // Helper: SHAKE bond-length & VdW projection pass
  const projectBackboneAndSterics = (ssPullAlpha: number) => {
    // 1. Two-way hard/adaptive disulfide clamp with hinge smoothing
    if (hasDisulfides) {
      for (const [i1, i2] of ssPairs) {
        const diffVec = pts[i2].clone().sub(pts[i1]);
        let dist = diffVec.length();
        if (dist < 1e-4) {
          diffVec.set(1, 0, 0);
          dist = 1e-4;
        }
        const dir = diffVec.divideScalar(dist);
        const delta = dist - TARGET_CA_CA_DIST;
        const corr = delta * 0.5 * ssPullAlpha;

        const move1 = dir.clone().multiplyScalar(corr);
        const move2 = dir.clone().multiplyScalar(-corr);

        pts[i1].add(move1);
        pts[i2].add(move2);

        // Smoothly propagate displacement to immediate flexible loop neighbors
        for (let offset = -2; offset <= 2; offset++) {
          if (offset === 0) continue;
          const decay = offset === -1 || offset === 1 ? 0.55 : 0.25;
          const n1 = i1 + offset;
          if (n1 >= 0 && n1 < n && !anchorSet.has(n1)) {
            pts[n1].addScaledVector(move1, decay * flexWeights[n1]);
          }
          const n2 = i2 + offset;
          if (n2 >= 0 && n2 < n && !anchorSet.has(n2)) {
            pts[n2].addScaledVector(move2, decay * flexWeights[n2]);
          }
        }
      }
    }

    // 2. Forward & backward SHAKE peptide bond length preservation (3.80 Å)
    for (let sweep = 0; sweep < 3; sweep++) {
      for (let i = 0; i < n - 1; i++) {
        const idxA = sweep % 2 === 0 ? i : n - 2 - i;
        const idxB = idxA + 1;
        const targetL = bondLengths[idxA];

        const seg = pts[idxB].clone().sub(pts[idxA]);
        let d = seg.length();
        if (d < 1e-5) {
          seg.set(0, 1, 0);
          d = 1e-5;
        }
        const err = d - targetL;
        // Disulfide anchor residues become stiffer as adaptive w_ss grows
        const wA = anchorSet.has(idxA) ? flexWeights[idxA] * (1 - 0.82 * ssPullAlpha) + 0.04 : flexWeights[idxA];
        const wB = anchorSet.has(idxB) ? flexWeights[idxB] * (1 - 0.82 * ssPullAlpha) + 0.04 : flexWeights[idxB];
        const wSum = wA + wB || 1.0;

        seg.divideScalar(d);
        pts[idxA].addScaledVector(seg, err * (wA / wSum));
        pts[idxB].addScaledVector(seg, -err * (wB / wSum));
      }
    }

    // 3. Virtual bond angle constraint (1-3 Cα-Cα distance in [5.05, 7.25] Å)
    for (let i = 1; i < n - 1; i++) {
      const span = pts[i + 1].clone().sub(pts[i - 1]);
      const d13 = span.length();
      if (d13 > 1e-4 && d13 < 5.05) {
        const push = (5.05 - d13) * 0.35;
        span.divideScalar(d13);
        if (!anchorSet.has(i - 1)) pts[i - 1].addScaledVector(span, -push);
        if (!anchorSet.has(i + 1)) pts[i + 1].addScaledVector(span, push);
      }
    }

    // 4. Non-bonded Van der Waals steric clash exclusion (Track 2)
    for (let i = 0; i < n; i++) {
      for (let j = i + 2; j < n; j++) {
        const isPairedCys = ssPairs.some(([a, b]) => a === i && b === j);
        const minAllowed = isPairedCys ? 4.60 : 3.55;
        const sep = pts[j].clone().sub(pts[i]);
        const d = sep.length();
        if (d < minAllowed) {
          if (d < 1e-4) sep.set(0.577, 0.577, 0.577);
          else sep.divideScalar(d);
          const push = (minAllowed - d) * 0.5;
          pts[i].addScaledVector(sep, -push);
          pts[j].addScaledVector(sep, push);
        }
      }
    }
  };

  // ---------------------------------------------------------------------------
  // ITERATIVE OPTIMIZATION LOOP (TWO-WAY CCD IK + ADAPTIVE SPRING + SAXS D_MAX)
  // ---------------------------------------------------------------------------
  const effectiveIterations = hasDisulfides ? maxIterations : Math.min(8, maxIterations);

  for (let iter = 1; iter <= effectiveIterations; iter++) {
    // Update Adaptive Spring Constant w_ss^(k)
    if (hasDisulfides) {
      if (exponentialBoostTriggered) {
        // Exponential growth: w_ss^(k) = w_ss^(0) * exp(0.19 * k)
        currentWSs = Math.min(650.0, baseSpringWeight * Math.exp(0.19 * iter));
      } else {
        currentWSs = Math.min(180.0, baseSpringWeight * (1.0 + 0.25 * iter));
      }
    }

    const ssStiffnessAlpha = hasDisulfides
      ? Math.min(0.96, currentWSs / (currentWSs + kappaBend))
      : 0.0;

    // A. Bidirectional Cyclic Coordinate Descent (Two-Way CCD IK) localized to loop [i1, i2]
    //    Translates exterior tails (q > i2 and q < i1) without lever-arm whipping of outer clamps!
    if (hasDisulfides) {
      const orderedPairs =
        iter % 2 === 1
          ? [...ssPairs].sort((a, b) => (a[1] - a[0]) - (b[1] - b[0])) // Inner-to-outer
          : [...ssPairs].sort((a, b) => (b[1] - b[0]) - (a[1] - a[0])); // Outer-to-inner

      for (let ccdSweep = 0; ccdSweep < 2; ccdSweep++) {
        for (const [i1, i2] of orderedPairs) {
          const curGap = pts[i2].clone().sub(pts[i1]);
          const curDist = curGap.length();
          if (Math.abs(curDist - TARGET_CA_CA_DIST) < 0.02) continue;

          const midIdx = Math.floor((i1 + i2) / 2);
          const adaptiveMaxStep = Math.abs(curDist - TARGET_CA_CA_DIST) > 6.0 ? 0.42 : 0.22;

          // Pass 1: C-side CCD hinges preceding i2 (from i2 - 1 down to midIdx)
          for (let p = i2 - 1; p >= Math.max(i1 + 1, midIdx); p--) {
            const pivot = pts[p];
            const gapDir = pts[i2].clone().sub(pts[i1]);
            if (gapDir.lengthSq() < 1e-6) gapDir.set(1, 0, 0);
            gapDir.normalize();

            const desiredEnd = pts[i1].clone().addScaledVector(gapDir, TARGET_CA_CA_DIST);
            const vCur = pts[i2].clone().sub(pivot);
            const vDes = desiredEnd.sub(pivot);

            if (vCur.lengthSq() > 1e-5 && vDes.lengthSq() > 1e-5) {
              vCur.normalize();
              vDes.normalize();
              const rotAxis = new THREE.Vector3().crossVectors(vCur, vDes);
              if (rotAxis.lengthSq() > 1e-6) {
                rotAxis.normalize();
                const fullAngle = Math.acos(Math.max(-1, Math.min(1, vCur.dot(vDes))));
                const dampedAngle =
                  Math.min(fullAngle, adaptiveMaxStep) *
                  flexWeights[p] *
                  (0.35 + 0.65 * ssStiffnessAlpha);

                if (dampedAngle > 1e-4) {
                  const oldEnd = pts[i2].clone();
                  for (let q = p + 1; q <= i2; q++) {
                    pts[q].sub(pivot).applyAxisAngle(rotAxis, dampedAngle).add(pivot);
                  }
                  const deltaEnd = pts[i2].clone().sub(oldEnd);
                  for (let q = i2 + 1; q < n; q++) {
                    pts[q].add(deltaEnd);
                  }
                }
              }
            }
          }

          // Pass 2: N-side Two-Way CCD hinges following i1 (from i1 + 1 up to midIdx)
          for (let p = i1 + 1; p <= Math.min(i2 - 1, midIdx); p++) {
            const pivot = pts[p];
            const gapDir = pts[i1].clone().sub(pts[i2]);
            if (gapDir.lengthSq() < 1e-6) gapDir.set(-1, 0, 0);
            gapDir.normalize();

            const desiredStart = pts[i2].clone().addScaledVector(gapDir, TARGET_CA_CA_DIST);
            const vCur = pts[i1].clone().sub(pivot);
            const vDes = desiredStart.sub(pivot);

            if (vCur.lengthSq() > 1e-5 && vDes.lengthSq() > 1e-5) {
              vCur.normalize();
              vDes.normalize();
              const rotAxis = new THREE.Vector3().crossVectors(vCur, vDes);
              if (rotAxis.lengthSq() > 1e-6) {
                rotAxis.normalize();
                const fullAngle = Math.acos(Math.max(-1, Math.min(1, vCur.dot(vDes))));
                const dampedAngle =
                  Math.min(fullAngle, adaptiveMaxStep) *
                  flexWeights[p] *
                  (0.35 + 0.65 * ssStiffnessAlpha);

                if (dampedAngle > 1e-4) {
                  const oldStart = pts[i1].clone();
                  for (let q = i1; q < p; q++) {
                    pts[q].sub(pivot).applyAxisAngle(rotAxis, dampedAngle).add(pivot);
                  }
                  const deltaStart = pts[i1].clone().sub(oldStart);
                  for (let q = 0; q < i1; q++) {
                    pts[q].add(deltaStart);
                  }
                }
              }
            }
          }
        }
      }
    }

    // B. Coupled SHAKE constraint & steric relaxation sub-steps
    const subSteps = hasDisulfides ? 14 : 3;
    for (let s = 0; s < subSteps; s++) {
      projectBackboneAndSterics(ssStiffnessAlpha);
    }

    // C. Adaptive E8 Bounding Sphere Buffer (R_kula >= 1.45 * R_g):
    // Never radially crush outer loops or sidechains! Only gently damp extreme runaway outliers (> 2.2 * R_g).
    if (enforceDMax && dMaxTarget > 8.0) {
      recenterPoints();
      let rgSqSum = 0;
      for (let i = 0; i < n; i++) {
        rgSqSum += pts[i].lengthSq();
      }
      const curRg = Math.sqrt(rgSqSum / Math.max(1, n));
      const safeAdaptiveRadius = Math.max(dMaxTarget * 0.65, curRg * 1.48 + 4.0);
      let outlierNudged = false;
      for (let i = 0; i < n; i++) {
        const r = pts[i].length();
        if (r > safeAdaptiveRadius * 1.35) {
          const softScale = (safeAdaptiveRadius * 1.35 + (r - safeAdaptiveRadius * 1.35) * 0.5) / r;
          pts[i].multiplyScalar(softScale);
          outlierNudged = true;
        }
      }
      if (outlierNudged) {
        for (let s = 0; s < 6; s++) {
          projectBackboneAndSterics(Math.max(0.90, ssStiffnessAlpha));
        }
      }
    }

    // Record iteration trace diagnostics
    const iterCoords: [number, number, number][] = pts.map(v => [v.x, v.y, v.z]);
    const iterSc = reconstructSidechainCbAndSg(sequenceCodes, iterCoords, ssPairs);

    let maxCaDev = 0;
    let maxSgDev = 0;
    let totalEss = 0;

    for (const [i1, i2] of ssPairs) {
      const dCa = pts[i1].distanceTo(pts[i2]);
      maxCaDev = Math.max(maxCaDev, Math.abs(dCa - TARGET_CA_CA_DIST));

      const sg1 = iterSc.sgPositions[i1];
      const sg2 = iterSc.sgPositions[i2];
      const dSg =
        sg1 && sg2
          ? Math.sqrt((sg1[0] - sg2[0]) ** 2 + (sg1[1] - sg2[1]) ** 2 + (sg1[2] - sg2[2]) ** 2)
          : TARGET_SG_SG_DIST;
      maxSgDev = Math.max(maxSgDev, Math.abs(dSg - TARGET_SG_SG_DIST));
      totalEss += computeAdaptiveSsbondEnergy(dSg, currentWSs, TARGET_SG_SG_DIST);
    }

    // Measure angular error at current iteration
    let iterAngSum = 0;
    for (let i = 0; i < n - 1; i++) {
      const bVec = pts[i + 1].clone().sub(pts[i]);
      const inAx = i > 0 ? pts[i].clone().sub(pts[i - 1]) : new THREE.Vector3(0, 1, 0);
      const snap = snapVectorToDiscreteMicroState(bVec, inAx, lattice);
      iterAngSum += snap.angleErrorDeg;
    }
    const iterMeanAng = n > 1 ? iterAngSum / (n - 1) : 0;

    let iterClashes = 0;
    for (let i = 0; i < n; i++) {
      for (let j = i + 2; j < n; j++) {
        if (pts[i].distanceTo(pts[j]) < 3.35) iterClashes++;
      }
    }

    iterationTraces.push({
      iteration: iter,
      wSs: parseFloat(currentWSs.toFixed(2)),
      maxCaDeviation: parseFloat(maxCaDev.toFixed(3)),
      maxSgDeviation: parseFloat(maxSgDev.toFixed(3)),
      totalSsbondEnergy: parseFloat(totalEss.toFixed(3)),
      meanAngularErrorDeg: parseFloat(iterMeanAng.toFixed(3)),
      vdwClashes: iterClashes,
      dMaxActual: parseFloat(computeCurrentDMax().toFixed(2))
    });
  }

  // Final high-precision bond-length & two-way clamp polish before L2 quantization
  const polishSteps = hasDisulfides ? 28 : 14;
  for (let polish = 0; polish < polishSteps; polish++) {
    projectBackboneAndSterics(hasDisulfides ? 0.98 : 0.0);
  }

  if (centerAtOrigin) {
    recenterPoints();
  }

  const continuousClosedPositions: [number, number, number][] = pts.map(v => [
    parseFloat(v.x.toFixed(3)),
    parseFloat(v.y.toFixed(3)),
    parseFloat(v.z.toFixed(3))
  ]);

  // ---------------------------------------------------------------------------
  // FINAL STAGE: DRIFT-COMPENSATED MICRO-FACET QUANTIZATION (ZERO LEVER-ARM DRIFT)
  // ---------------------------------------------------------------------------
  const quantizedPts: THREE.Vector3[] = [pts[0].clone()];
  const perBondPatchIds: number[] = [];
  const perBondNormals: [number, number, number][] = [];
  const perBondAngularErrorsDeg: number[] = [];

  for (let i = 0; i < n - 1; i++) {
    // Closed-loop tracking vector aims from current quantized point to target IK-closed point[i+1]
    // This eliminates lever-arm error accumulation across long chains!
    const closedLoopTargetVec = pts[i + 1].clone().sub(quantizedPts[i]);
    const continuousRefVec = pts[i + 1].clone().sub(pts[i]);
    const incomingAxis =
      i > 0 ? quantizedPts[i].clone().sub(quantizedPts[i - 1]) : new THREE.Vector3(0, 1, 0);

    const stepLen = continuousRefVec.length() || bondLengths[i];
    const snap = snapVectorToDiscreteMicroState(closedLoopTargetVec, incomingAxis, lattice);

    const nextQuantized = quantizedPts[i].clone().addScaledVector(snap.snappedUnitDir, stepLen);
    quantizedPts.push(nextQuantized);

    // True 3D angular deviation against continuous closed-loop bond vector
    const refDir = continuousRefVec.clone().normalize();
    const dot = Math.max(-1, Math.min(1, snap.snappedUnitDir.dot(refDir)));
    const trueErrDeg = parseFloat(((Math.acos(dot) * 180.0) / Math.PI).toFixed(3));

    perBondPatchIds.push(snap.patchId);
    perBondNormals.push([
      parseFloat(snap.snappedUnitDir.x.toFixed(4)),
      parseFloat(snap.snappedUnitDir.y.toFixed(4)),
      parseFloat(snap.snappedUnitDir.z.toFixed(4))
    ]);
    perBondAngularErrorsDeg.push(trueErrDeg);
  }

  if (centerAtOrigin) {
    const qCom = new THREE.Vector3(0, 0, 0);
    for (let i = 0; i < n; i++) qCom.add(quantizedPts[i]);
    qCom.divideScalar(Math.max(1, n));
    for (let i = 0; i < n; i++) quantizedPts[i].sub(qCom);
  } else {
    // Pure translation to align quantized centroid with continuous target centroid (0.000° impact on L2 normals)
    const residualBias = new THREE.Vector3(0, 0, 0);
    for (let i = 0; i < n; i++) {
      residualBias.add(quantizedPts[i]).sub(pts[i]);
    }
    residualBias.divideScalar(Math.max(1, n));
    for (let i = 0; i < n; i++) quantizedPts[i].sub(residualBias);
  }

  const finalPositions: [number, number, number][] = quantizedPts.map(v => [
    parseFloat(v.x.toFixed(3)),
    parseFloat(v.y.toFixed(3)),
    parseFloat(v.z.toFixed(3))
  ]);

  const { cbPositions, sgPositions } = reconstructSidechainCbAndSg(
    sequenceCodes,
    finalPositions,
    ssPairs
  );

  // Build final per-bridge diagnostics
  const bridgeDiagnostics: DisulfideBridgeDiagnostic[] = ssPairs.map(([i1, i2], bIdx) => {
    const p1 = finalPositions[i1];
    const p2 = finalPositions[i2];
    const finalCaDist = Math.sqrt(
      (p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2 + (p1[2] - p2[2]) ** 2
    );
    const caDev = Math.abs(finalCaDist - TARGET_CA_CA_DIST);

    const sg1 = sgPositions[i1];
    const sg2 = sgPositions[i2];
    const finalSgDist =
      sg1 && sg2
        ? Math.sqrt((sg1[0] - sg2[0]) ** 2 + (sg1[1] - sg2[1]) ** 2 + (sg1[2] - sg2[2]) ** 2)
        : TARGET_SG_SG_DIST;
    const sgDev = Math.abs(finalSgDist - TARGET_SG_SG_DIST);

    const eSsbond = computeAdaptiveSsbondEnergy(finalSgDist, currentWSs, TARGET_SG_SG_DIST);

    // Identify top flexible hinge residues relaxed within [i1, i2]
    const hinges: { idx: number; flex: number }[] = [];
    for (let k = i1 + 1; k < i2; k++) {
      hinges.push({ idx: k + 1, flex: flexWeights[k] });
    }
    hinges.sort((a, b) => b.flex - a.flex);
    const topHinges = hinges.slice(0, 4).map(h => h.idx).sort((a, b) => a - b);

    const initMetric = initialBridgeMetrics[bIdx] || {
      caDist: finalCaDist,
      sgDist: finalSgDist
    };

    return {
      bridgeIndex: bIdx + 1,
      idx1: i1,
      idx2: i2,
      res1: i1 + 1,
      res2: i2 + 1,
      initialCaDist: parseFloat(initMetric.caDist.toFixed(2)),
      finalCaDist: parseFloat(finalCaDist.toFixed(2)),
      caDeviation: parseFloat(caDev.toFixed(2)),
      initialSgDist: parseFloat(initMetric.sgDist.toFixed(2)),
      finalSgDist: parseFloat(finalSgDist.toFixed(2)),
      sgDeviation: parseFloat(sgDev.toFixed(2)),
      energySsbond: parseFloat(eSsbond.toFixed(4)),
      hingeResidues: topHinges,
      passedSubHalfAngstrom: sgDev < 0.5 && caDev < 0.5
    };
  });

  const finalMaxSgDev = bridgeDiagnostics.reduce((m, b) => Math.max(m, b.sgDeviation), 0);
  const finalMaxCaDev = bridgeDiagnostics.reduce((m, b) => Math.max(m, b.caDeviation), 0);
  const totalSsbondEnergy = bridgeDiagnostics.reduce((s, b) => s + b.energySsbond, 0);

  const meanFacetError =
    perBondAngularErrorsDeg.length > 0
      ? perBondAngularErrorsDeg.reduce((a, b) => a + b, 0) / perBondAngularErrorsDeg.length
      : 0;
  const maxFacetError =
    perBondAngularErrorsDeg.length > 0 ? Math.max(...perBondAngularErrorsDeg) : 0;

  let finalClashes = 0;
  for (let i = 0; i < n; i++) {
    for (let j = i + 2; j < n; j++) {
      const p1 = finalPositions[i];
      const p2 = finalPositions[j];
      const d = Math.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2 + (p1[2] - p2[2]) ** 2);
      if (d < 3.35) finalClashes++;
    }
  }

  return {
    pipelineMode,
    positions: finalPositions,
    continuousClosedPositions,
    cbPositions,
    sgPositions,
    ssPairs,
    bridgeDiagnostics,
    initialMaxDeviationAngstrom: parseFloat(initialMaxDev.toFixed(2)),
    finalMaxSgDeviationAngstrom: parseFloat(finalMaxSgDev.toFixed(2)),
    finalMaxCaDeviationAngstrom: parseFloat(finalMaxCaDev.toFixed(2)),
    exponentialBoostTriggered,
    initialSpringWeight: parseFloat(baseSpringWeight.toFixed(2)),
    finalSpringWeight: parseFloat(currentWSs.toFixed(2)),
    totalSsbondPenaltyEnergy: parseFloat(totalSsbondEnergy.toFixed(4)),
    meanMicroFacetErrorDeg: parseFloat(meanFacetError.toFixed(3)),
    maxMicroFacetErrorDeg: parseFloat(maxFacetError.toFixed(3)),
    perBondPatchIds,
    perBondNormals,
    perBondAngularErrorsDeg,
    vdwClashesCount: finalClashes,
    iterationTraces
  };
}
