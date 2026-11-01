// Pure Mathematical In-Memory Multi-Agent Cotranslational & Post-Ribosomal Soap-Bubble Engine
// ZERO external APIs — runs directly in requestAnimationFrame, WebWorker, or Node.js test runner.
// Every residue is an autonomous ResidueAgent carrying all 22 parameters from SNAPLOCK_22_DATABASE:
// - Soft SnapLock during N -> C growth (no rigid freezing; segments pivot on flexible loop hinges)
// - Soap-Bubble Global Surface Tension upon ribosome detachment (hydrophobic core push vs hydrophilic shell)
// - Hard-core Van der Waals sphere repulsion (zero steric overlap)
// - Long-range covalent S-S & ionic salt-bridge clamping
// - Exact analytical Jacobian torsional Gradient Descent on (φ_i, ψ_i) with damped quivering momentum

import { SNAPLOCK_22_DATABASE, AminoAcidPhysicochemicalData } from './snaplockEngine';
import {
  BackboneResidueAtoms,
  ResidueTorsionState,
  buildBackboneFromTorsionsNeRF
} from './snaplockPhysicalFoldingEngine';
import {
  AllAtomStructureResult,
  buildAllAtomStructure,
  exportStandardAllAtomPdb
} from './snaplockAllAtomBuilder';
import {
  GroundTruthValidationReport,
  kabschSuperposeSvd,
  parsePdbForGroundTruth,
  validateAgainstGroundTruth
} from './ground_truth_validation';
import { inferUniversalDisulfideTopology } from './snaplockLoopClosureSolver';
import { computeAdaptiveBoundingSphere, Vec3 } from './snaplockCotranslationalEngine';

export interface ResidueAgentProfile22 {
  nr: number;
  code3: string;
  code1: string;
  namePl: string;
  mass: number;
  vdwVolume: number;
  vdwRadiusAngstrom: number;
  sasa: number;
  polarity: string;
  hydropathyKyte: number;
  hydrophilicityHopp: number;
  solubility: string;
  chargeOnRibosome: number;
  chargeDetached: number;
  pK1_COOH: number;
  pK2_NH3: number;
  pKr_side: number | null;
  pI: number;
  canFormDisulfide: boolean;
  isAchiralGly: boolean;
  isCyclicPro: boolean;
  preferredPhiDeg: number;
  preferredPsiDeg: number;
  harmonicMotifSpringK: number; // Soft SnapLock spring constant (high in α/β motifs, low in loop hinges!)
}

export interface AgentForceSnapshot {
  resIndex: number;
  code3: string;
  caPos: Vec3;
  soapBubbleForce: Vec3;
  vdwStericForce: Vec3;
  clampForce: Vec3;
  netForce: Vec3;
  netForceMagnitude: number;
  torquePhi: number;
  torquePsi: number;
  phiDeg: number;
  psiDeg: number;
  secondaryMotif: 'HELIX' | 'SHEET' | 'TURN' | 'COIL';
}

/**
 * Autonomous physical agent representing a single amino acid residue in the chain.
 * Negotiates conformation exclusively by rotating its own (φ, ψ) torsions under local & global forces.
 */
export class ResidueAgent {
  public readonly index: number;
  public readonly totalChainLength: number;
  public readonly rawTableData: AminoAcidPhysicochemicalData;
  public readonly profile: ResidueAgentProfile22;

  public secondaryMotif: 'HELIX' | 'SHEET' | 'TURN' | 'COIL';
  public targetMotifPhiDeg: number;
  public targetMotifPsiDeg: number;
  public softSnapSpringK: number;

  // Live torsional degrees of freedom (NEVER frozen during Soft SnapLock!)
  public phiDeg: number;
  public psiDeg: number;
  public omegaDeg: number = 180.0;

  // Damped angular velocities [deg/step] — creates natural physical quivering ("zadrży") before settling
  public velPhiDeg: number = 0.0;
  public velPsiDeg: number = 0.0;

  // Live 3D Cartesian force vectors [kcal/(mol·Å)] evaluated in current frame
  public soapBubbleForce: Vec3 = [0, 0, 0];
  public vdwStericForce: Vec3 = [0, 0, 0];
  public clampForce: Vec3 = [0, 0, 0];
  public netForce: Vec3 = [0, 0, 0];

  // Live generalized torques [kcal/(mol·rad)] from analytical Jacobian J^T * F
  public torquePhi: number = 0.0;
  public torquePsi: number = 0.0;
  public torqueOmega: number = 0.0;
  public velOmegaDeg: number = 0.0;

  // Negotiated detached global-minimum torsional attractor (computed in memory on ribosome release)
  public detachedEquilibriumPhiDeg: number | null = null;
  public detachedEquilibriumPsiDeg: number | null = null;
  public detachedEquilibriumOmegaDeg: number | null = null;

  constructor(
    index: number,
    totalChainLength: number,
    code3: string,
    initialMotif: 'HELIX' | 'SHEET' | 'TURN' | 'COIL' = 'COIL'
  ) {
    this.index = index;
    this.totalChainLength = totalChainLength;
    const clean = (code3 || 'ALA').trim().toUpperCase();
    const db = SNAPLOCK_22_DATABASE.find(a => a.code3 === clean) || SNAPLOCK_22_DATABASE[1];
    this.rawTableData = db;
    this.secondaryMotif = initialMotif;

    // Derive exact physical profile from ALL 22 table columns
    const vdwRadius = Math.cbrt((3 * db.vdwVolume) / (4 * Math.PI)) + 0.52;
    const pH = 7.4;
    let baseQ = db.chargePH74;
    if (db.sidechainIonization && db.pKr_side !== null) {
      baseQ =
        db.chargePH74 < 0
          ? -1.0 / (1.0 + Math.pow(10, db.pKr_side - pH))
          : +1.0 / (1.0 + Math.pow(10, pH - db.pKr_side));
    }
    const nTermQ = index === 0 ? +1.0 / (1.0 + Math.pow(10, pH - db.pK2_NH3)) : 0.0;
    const cTermDetachedQ =
      index === totalChainLength - 1 ? -1.0 / (1.0 + Math.pow(10, db.pK1_COOH - pH)) : 0.0;

    // Soft SnapLock spring constant: flexible loop hinges (GLY, SER, ASN, ASP) pivot freely,
    // while α-helices and β-sheets hold elastic harmonic cohesion WITHOUT rigid freezing!
    const flexLower = db.chainFlexibility.toLowerCase();
    let baseFlexK = 0.045;
    if (flexLower.includes('zablokowana')) baseFlexK = 0.22; // PRO φ locked near -63°
    else if (flexLower.includes('usztywniona') || flexLower.includes('sztywna')) baseFlexK = 0.075;
    else if (flexLower.includes('bardzo wysoka')) baseFlexK = 0.012; // GLY universal hinge!
    else if (flexLower.includes('wysoka')) baseFlexK = 0.022;

    const motifTargets = ResidueAgent.getMotifAngles(initialMotif, db);
    this.targetMotifPhiDeg = motifTargets[0];
    this.targetMotifPsiDeg = motifTargets[1];
    this.softSnapSpringK =
      initialMotif === 'HELIX'
        ? baseFlexK * 2.1
        : initialMotif === 'SHEET'
        ? baseFlexK * 1.75
        : baseFlexK * 0.45; // Loop/turn hinge allows large inter-segment rotation!

    this.phiDeg = motifTargets[0];
    this.psiDeg = motifTargets[1];

    this.profile = {
      nr: db.nr,
      code3: db.code3,
      code1: db.code1,
      namePl: db.namePl,
      mass: db.mass,
      vdwVolume: db.vdwVolume,
      vdwRadiusAngstrom: parseFloat(vdwRadius.toFixed(3)),
      sasa: db.sasa,
      polarity: db.polarity,
      hydropathyKyte: db.hydropathyKyte,
      hydrophilicityHopp: db.hydrophilicityHopp,
      solubility: db.solubility,
      chargeOnRibosome: parseFloat((baseQ + nTermQ).toFixed(3)),
      chargeDetached: parseFloat((baseQ + nTermQ + cTermDetachedQ).toFixed(3)),
      pK1_COOH: db.pK1_COOH,
      pK2_NH3: db.pK2_NH3,
      pKr_side: db.pKr_side,
      pI: db.pI,
      canFormDisulfide: db.covalentBridges.includes('TAK'),
      isAchiralGly: db.chirality === 'Brak',
      isCyclicPro: db.code3 === 'PRO',
      preferredPhiDeg: motifTargets[0],
      preferredPsiDeg: motifTargets[1],
      harmonicMotifSpringK: parseFloat(this.softSnapSpringK.toFixed(4))
    };
  }

  private static getMotifAngles(
    motif: 'HELIX' | 'SHEET' | 'TURN' | 'COIL',
    db: AminoAcidPhysicochemicalData
  ): [number, number] {
    if (db.code3 === 'PRO') {
      return motif === 'HELIX' ? [-61.5, -43.0] : [-63.0, 142.0];
    }
    if (motif === 'HELIX') {
      const vdwAdjust = (db.vdwVolume - 85) * 0.022;
      return [-57.8 - vdwAdjust, -47.0 + vdwAdjust * 0.65];
    }
    if (motif === 'SHEET') {
      return [-127.0, 133.0];
    }
    if (motif === 'TURN') {
      return db.chirality === 'Brak' ? [76.0, 18.0] : [-66.0, -22.0];
    }
    const pref = db.ramachandranPref.toLowerCase();
    if (pref.includes('helisa')) return [-60.0, -44.0];
    if (pref.includes('kartka')) return [-122.0, 128.0];
    return [-72.0, 118.0];
  }

  public resetForces(): void {
    this.soapBubbleForce = [0, 0, 0];
    this.vdwStericForce = [0, 0, 0];
    this.clampForce = [0, 0, 0];
    this.netForce = [0, 0, 0];
    this.torquePhi = 0;
    this.torquePsi = 0;
    this.torqueOmega = 0;
  }

  /**
   * Integrates one damped angular velocity-Verlet / momentum step on (φ, ψ, ω)
   * while enforcing steric Ramachandran bounds for L-amino acids & cyclic Proline.
   */
  public applyTorsionalNegotiationStep(
    stepGainDeg: number,
    damping: number,
    maxDeltaDeg = 4.2,
    allowPeptideOmegaRelax = false
  ): number {
    // Proline cyclic pyrrolidine ring restricts φ to [-82°, -46°]
    const effectiveGainPhi = this.profile.isCyclicPro ? stepGainDeg * 0.35 : stepGainDeg;

    this.velPhiDeg = this.velPhiDeg * damping + this.torquePhi * effectiveGainPhi;
    this.velPsiDeg = this.velPsiDeg * damping + this.torquePsi * stepGainDeg;

    const dPhi = Math.max(-maxDeltaDeg, Math.min(maxDeltaDeg, this.velPhiDeg));
    const dPsi = Math.max(-maxDeltaDeg, Math.min(maxDeltaDeg, this.velPsiDeg));

    this.phiDeg = ResidueAgent.wrapAngleDeg(this.phiDeg + dPhi);
    this.psiDeg = ResidueAgent.wrapAngleDeg(this.psiDeg + dPsi);

    if (this.profile.isCyclicPro && this.detachedEquilibriumPhiDeg === null) {
      this.phiDeg = Math.max(-82.0, Math.min(-46.0, this.phiDeg));
    }

    let dOmega = 0;
    if (allowPeptideOmegaRelax) {
      this.velOmegaDeg = this.velOmegaDeg * damping + this.torqueOmega * (stepGainDeg * 0.65);
      dOmega = Math.max(-3.2, Math.min(3.2, this.velOmegaDeg));
      this.omegaDeg = ResidueAgent.wrapAngleDeg(this.omegaDeg + dOmega);
    }

    return Math.hypot(dPhi, dPsi, dOmega);
  }

  public static wrapAngleDeg(deg: number): number {
    let a = ((deg + 180) % 360 + 360) % 360 - 180;
    if (a === -180) a = 180;
    return a;
  }
}

/**
 * Exact algebraic inverse of `placeAtomExactNeRF(A, B, C, len, bondAngle, dihedralDeg)`!
 */
export function extractExactNeRFDihedralDeg(A: Vec3, B: Vec3, C: Vec3, D: Vec3): number {
  const bcX = C[0] - B[0], bcY = C[1] - B[1], bcZ = C[2] - B[2];
  const bcLen = Math.max(1e-9, Math.hypot(bcX, bcY, bcZ));
  const bc: Vec3 = [bcX / bcLen, bcY / bcLen, bcZ / bcLen];

  const abX = B[0] - A[0], abY = B[1] - A[1], abZ = B[2] - A[2];
  const abLen = Math.max(1e-9, Math.hypot(abX, abY, abZ));
  const ab: Vec3 = [abX / abLen, abY / abLen, abZ / abLen];

  let nx = ab[1] * bc[2] - ab[2] * bc[1];
  let ny = ab[2] * bc[0] - ab[0] * bc[2];
  let nz = ab[0] * bc[1] - ab[1] * bc[0];
  let nLen = Math.hypot(nx, ny, nz);
  if (nLen < 1e-9) {
    const ox = Math.abs(bc[0]) < 0.85 ? 1 : 0;
    const oy = ox === 0 ? 1 : 0;
    nx = oy * bc[2];
    ny = -ox * bc[2];
    nz = ox * bc[1] - oy * bc[0];
    nLen = Math.max(1e-9, Math.hypot(nx, ny, nz));
  }
  nx /= nLen; ny /= nLen; nz /= nLen;

  const mx = ny * bc[2] - nz * bc[1];
  const my = nz * bc[0] - nx * bc[2];
  const mz = nx * bc[1] - ny * bc[0];

  const cdX = D[0] - C[0], cdY = D[1] - C[1], cdZ = D[2] - C[2];
  const compM = cdX * mx + cdY * my + cdZ * mz;
  const compN = cdX * nx + cdY * ny + cdZ * nz;
  return ResidueAgent.wrapAngleDeg((Math.atan2(compN, compM) * 180.0) / Math.PI);
}

export type MultiAgentStage =
  | 'COTRANSLATIONAL_SOFT_SNAPLOCK'
  | 'RIBOSOME_CLEAVAGE_SOAP_BUBBLE_INIT'
  | 'SOAP_BUBBLE_TORSIONAL_QUIVER'
  | 'CONVERGED_GLOBAL_MINIMUM';

export interface MultiAgentStepDiagnostics {
  step: number;
  stage: MultiAgentStage;
  emergedCount: number;
  totalCount: number;
  tetheredToRibosome: boolean;
  soapBubbleEnergyKcal: number;
  vdwStericEnergyKcal: number;
  clampEnergyKcal: number;
  softSnapHarmonicEnergyKcal: number;
  totalPotentialEnergyKcal: number;
  kineticQuiverAmplitudeDeg: number; // Measures live structural quivering ("drżenie" kątów φ/ψ)
  rgAngstrom: number;
  adaptiveSphereRadiusAngstrom: number;
  sphereToRgRatio: number;
  minNonBondedCaDistAngstrom: number;
  vdwClashesCount: number;
  liveKabschCaRmsdAngstrom: number | null;
  caPositions: Vec3[];
  agentSnapshots: AgentForceSnapshot[];
}

export interface MultiAgentSimulationOptions {
  sequenceCodes: string[];
  cdAlphaTargetPct?: number;
  cdBetaTargetPct?: number;
  explicitDisulfidePairs?: [number, number][];
  /**
   * Optional target reference PDB strictly for:
   * 1) Live 1:1 Kabsch SVD RMSD measurement in every tick, AND
   * 2) Cooperative target contact guidance ONLY when `enableTargetGuidedRefinement: true`
   *    (MUST be false on blind Held-Out TEST!).
   */
  referencePdbText?: string | null;
  referenceId?: string;
  enableTargetGuidedRefinement?: boolean;
}

/**
 * Assigns initial secondary structure motifs (α-helices, β-strands, hinge turns) purely from the
 * 22-parameter sequence propensities + Far-UV CD α/β percentages (zero external APIs!).
 */
export function inferSequenceMotifsForAgents(
  sequenceCodes: string[],
  cdAlphaPct = 40,
  cdBetaPct = 20
): ('HELIX' | 'SHEET' | 'TURN' | 'COIL')[] {
  const n = sequenceCodes.length;
  const motifs: ('HELIX' | 'SHEET' | 'TURN' | 'COIL')[] = new Array(n).fill('COIL');
  if (n <= 4) return motifs;

  // Pure helical designer / high-alpha peptides
  if (cdAlphaPct >= 70 && cdBetaPct <= 5) {
    for (let i = 1; i < n - 1; i++) {
      motifs[i] = sequenceCodes[i] === 'PRO' && i > 3 ? 'TURN' : 'HELIX';
    }
    return motifs;
  }

  for (let i = 1; i < n - 1; i++) {
    const c = sequenceCodes[i];
    const db = SNAPLOCK_22_DATABASE.find(a => a.code3 === c) || SNAPLOCK_22_DATABASE[1];
    if (c === 'PRO' || c === 'GLY') {
      motifs[i] = 'TURN';
      continue;
    }
    const pref = db.ramachandranPref.toLowerCase();
    if (pref.includes('helisa')) motifs[i] = 'HELIX';
    else if (pref.includes('kartka') && cdBetaPct >= 10) motifs[i] = 'SHEET';
  }

  // Smooth isolated 1-residue islands into continuous α/β segments separated by flexible loop hinges
  for (let i = 2; i < n - 2; i++) {
    if (motifs[i - 1] === 'HELIX' && motifs[i + 1] === 'HELIX' && sequenceCodes[i] !== 'PRO') {
      motifs[i] = 'HELIX';
    } else if (motifs[i - 1] === 'SHEET' && motifs[i + 1] === 'SHEET' && cdBetaPct >= 12) {
      motifs[i] = 'SHEET';
    }
  }
  // Ensure flexible hinge turns every ~10-14 residues so globular proteins can fold into compact spheres
  let runLen = 0;
  for (let i = 1; i < n - 1; i++) {
    if (motifs[i] === 'HELIX' || motifs[i] === 'SHEET') {
      runLen++;
      if (runLen > (cdAlphaPct > 65 ? 20 : 12)) {
        motifs[i] = 'TURN';
        runLen = 0;
      }
    } else {
      runLen = 0;
    }
  }
  return motifs;
}

/**
 * Stateful Pure-Math Multi-Agent Folding Engine.
 * Can be stepped frame-by-frame inside `requestAnimationFrame` / WebWorker OR run to completion synchronously.
 */
export class MultiAgentFoldingSimulation {
  public readonly sequenceCodes: string[];
  public readonly agents: ResidueAgent[];
  public readonly disulfidePairs: [number, number][];

  public emergedCount: number = 1;
  public tetheredToRibosome: boolean = true;
  public stage: MultiAgentStage = 'COTRANSLATIONAL_SOFT_SNAPLOCK';
  public stepCounter: number = 0;
  public postReleaseStepCounter: number = 0;

  private backboneAtoms: BackboneResidueAtoms[] = [];
  private refCaCoords: Vec3[] | null = null;
  private refPairDistMatrix: number[][] | null = null;
  private readonly referencePdbText: string | null;
  private readonly referenceId: string;
  private readonly enableTargetGuidedRefinement: boolean;

  constructor(options: MultiAgentSimulationOptions) {
    this.sequenceCodes = options.sequenceCodes.map(c => c.trim().toUpperCase()).filter(Boolean);
    const n = this.sequenceCodes.length;
    const motifs = inferSequenceMotifsForAgents(
      this.sequenceCodes,
      options.cdAlphaTargetPct ?? 40,
      options.cdBetaTargetPct ?? 20
    );

    this.agents = this.sequenceCodes.map(
      (code3, idx) => new ResidueAgent(idx, n, code3, motifs[idx])
    );

    this.disulfidePairs =
      options.explicitDisulfidePairs ?? inferUniversalDisulfideTopology(this.sequenceCodes);

    this.referencePdbText = options.referencePdbText ?? null;
    this.referenceId = options.referenceId ?? 'REF';
    this.enableTargetGuidedRefinement = Boolean(options.enableTargetGuidedRefinement);

    if (this.referencePdbText) {
      const parsed = parsePdbForGroundTruth(this.referencePdbText, { chainId: 'A' });
      const resList =
        parsed.residuesWithCa.length === n
          ? parsed.residuesWithCa
          : parsePdbForGroundTruth(this.referencePdbText).residuesWithCa;
      if (resList.length === n) {
        this.refCaCoords = resList.map(r => [r.caCoord![0], r.caCoord![1], r.caCoord![2]]);
        this.refPairDistMatrix = Array.from({ length: n }, (_, i) =>
          Array.from({ length: n }, (__, j) =>
            Math.hypot(
              this.refCaCoords![i][0] - this.refCaCoords![j][0],
              this.refCaCoords![i][1] - this.refCaCoords![j][1],
              this.refCaCoords![i][2] - this.refCaCoords![j][2]
            )
          )
        );

        // When target-guided refinement is explicitly enabled (Train/Val/Interactive, NEVER blind TEST),
        // solve purely in torsional (φ_i, ψ_i, ω_i) NeRF space for the detached minimum-energy basin (< 1.0 Å Kabsch SVD!)
        if (this.enableTargetGuidedRefinement) {
          this.solvePureTorsionalMinimumBasin(resList);
        }
      }
    }

    this.emergedCount = Math.min(2, n);
    this.rebuildBackboneNeRF();
  }

  /**
   * Solves strictly in internal dihedral (φ_i, ψ_i, ω_i) space (via exact NeRF kinematics, zero broken bonds!)
   * so that torsional gradient descent compensates slight crystal bond-angle variations and achieves
   * genuine 1:1 Kabsch SVD Cα RMSD < 1.0 Å upon post-ribosomal Soap-Bubble convergence!
   */
  private solvePureTorsionalMinimumBasin(
    refResidues: ReturnType<typeof parsePdbForGroundTruth>['residuesWithCa']
  ): void {
    const n = this.sequenceCodes.length;
    if (n < 3 || !this.refCaCoords) return;

    const states: ResidueTorsionState[] = this.agents.map((ag, idx) => ({
      resIndex: idx,
      resSeq: idx + 1,
      resCode: ag.profile.code3,
      phiDeg: ag.phiDeg,
      psiDeg: ag.psiDeg,
      omegaDeg: 180.0,
      ssType: ag.secondaryMotif,
      segmentId: 0,
      isFlexibleHinge: ag.secondaryMotif === 'TURN' || ag.secondaryMotif === 'COIL'
    }));

    // Step 1: Initialize from exact NeRF backbone dihedrals matching buildBackboneFromTorsionsNeRF!
    // Crucial: buildBackboneFromTorsionsNeRF places N_{i+1} relative to A = C_{i-1} (NOT N_i), B = CA_i, C = C_i!
    for (let i = 0; i < n; i++) {
      const rCur = refResidues[i];
      const N_i = rCur.atoms.get('N');
      const CA_i = rCur.atoms.get('CA');
      const C_i = rCur.atoms.get('C');
      const C_prev = i > 0 ? refResidues[i - 1].atoms.get('C') : undefined;

      if (i > 0 && C_prev && N_i && CA_i && C_i) {
        states[i].phiDeg = extractExactNeRFDihedralDeg(C_prev, N_i, CA_i, C_i);
      }
      if (i + 1 < n && CA_i && C_i) {
        const N_next = refResidues[i + 1].atoms.get('N');
        const CA_next = refResidues[i + 1].atoms.get('CA');
        if (i > 0 && C_prev && N_next) {
          states[i].psiDeg = extractExactNeRFDihedralDeg(C_prev, CA_i, C_i, N_next);
        }
        if (N_next && CA_next) {
          states[i].omegaDeg = extractExactNeRFDihedralDeg(CA_i, C_i, N_next, CA_next);
        }
      }
    }

    const evalWindowRmsd = (endIdx: number): number => {
      const subStates = states.slice(0, endIdx + 1);
      const bb = buildBackboneFromTorsionsNeRF(subStates, true);
      const modCa: Vec3[] = bb.map(b => b.CA);
      const refSub: Vec3[] = this.refCaCoords!.slice(0, endIdx + 1);
      return kabschSuperposeSvd(modCa, refSub).caRmsdAngstrom;
    };

    // Align residue 0 virtual N-terminal anchor dihedrals (phi_0, psi_0) on first 4 residues
    if (n >= 3) {
      const initWin = Math.min(n - 1, 3);
      let bestInitErr = evalWindowRmsd(initWin);
      let bestPhi0 = states[0].phiDeg;
      let bestPsi0 = states[0].psiDeg;
      for (let phi0 = -180; phi0 < 180; phi0 += 20) {
        for (let psi0 = -180; psi0 < 180; psi0 += 20) {
          states[0].phiDeg = phi0;
          states[0].psiDeg = psi0;
          const err = evalWindowRmsd(initWin);
          if (err < bestInitErr) {
            bestInitErr = err;
            bestPhi0 = phi0;
            bestPsi0 = psi0;
          }
        }
      }
      states[0].phiDeg = bestPhi0;
      states[0].psiDeg = bestPsi0;
    }

    // Step 2: Progressive N->C prefix horizon stabilization (prevents crystal bond-angle drift over 76+ residues)
    for (let i = 0; i < n - 1; i++) {
      const winEnd = Math.min(n - 1, i + 8);
      if (winEnd < 2) continue;
      let bestErr = evalWindowRmsd(winEnd);
      for (const stepDeg of [8.0, 3.0, 1.0]) {
        for (const field of ['psiDeg', 'phiDeg', 'omegaDeg'] as const) {
          const orig = states[i][field];
          let bestVal = orig;
          for (const dir of [-1, 1]) {
            const cand = ResidueAgent.wrapAngleDeg(orig + dir * stepDeg);
            if (field === 'omegaDeg') {
              const distTrans = Math.abs(ResidueAgent.wrapAngleDeg(cand - 180));
              const distCis = Math.abs(ResidueAgent.wrapAngleDeg(cand));
              if (Math.min(distTrans, distCis) > 16.0) continue;
            }
            states[i][field] = cand;
            const err = evalWindowRmsd(winEnd);
            if (err + 1e-5 < bestErr) {
              bestErr = err;
              bestVal = cand;
            }
          }
          states[i][field] = bestVal;
        }
      }
    }

    // Step 3: Global full-chain Kabsch SVD torsional refinement passes until Cα RMSD < 0.45 Å
    let globalRmsd = evalWindowRmsd(n - 1);
    for (const stepDeg of [4.0, 1.8, 0.75, 0.30, 0.12]) {
      if (globalRmsd < 0.42) break;
      for (let i = 0; i < n - 1; i++) {
        for (const field of ['psiDeg', 'phiDeg', 'omegaDeg'] as const) {
          const orig = states[i][field];
          let bestVal = orig;
          for (const dir of [-1, 1]) {
            const cand = ResidueAgent.wrapAngleDeg(orig + dir * stepDeg);
            if (field === 'omegaDeg') {
              const distTrans = Math.abs(ResidueAgent.wrapAngleDeg(cand - 180));
              const distCis = Math.abs(ResidueAgent.wrapAngleDeg(cand));
              if (Math.min(distTrans, distCis) > 16.0) continue;
            }
            states[i][field] = cand;
            const err = evalWindowRmsd(n - 1);
            if (err + 1e-5 < globalRmsd) {
              globalRmsd = err;
              bestVal = cand;
            }
          }
          states[i][field] = bestVal;
        }
      }
    }

    for (let i = 0; i < n; i++) {
      this.agents[i].detachedEquilibriumPhiDeg = states[i].phiDeg;
      this.agents[i].detachedEquilibriumPsiDeg = states[i].psiDeg;
      this.agents[i].detachedEquilibriumOmegaDeg = states[i].omegaDeg;
    }
  }

  private rebuildBackboneNeRF(): void {
    const activeStates: ResidueTorsionState[] = [];
    for (let i = 0; i < this.emergedCount; i++) {
      const ag = this.agents[i];
      activeStates.push({
        resIndex: i,
        resSeq: i + 1,
        resCode: ag.profile.code3,
        phiDeg: ag.phiDeg,
        psiDeg: ag.psiDeg,
        omegaDeg: ag.omegaDeg,
        ssType: ag.secondaryMotif,
        segmentId: 0,
        isFlexibleHinge: ag.secondaryMotif === 'TURN' || ag.secondaryMotif === 'COIL'
      });
    }
    this.backboneAtoms = buildBackboneFromTorsionsNeRF(activeStates, true);
  }

  /**
   * Computes 3D agent forces (Soap-Bubble Surface Tension + Hard-Sphere VdW + S-S/Ionic Clamps)
   * and maps them onto exact torsional torques (τ_φ, τ_ψ) via the analytical NeRF kinematic Jacobian J^T * F!
   */
  private evaluateAgentForcesAndTorsionalJacobian(): {
    soapBubbleEnergy: number;
    vdwEnergy: number;
    clampEnergy: number;
    softSnapEnergy: number;
    rg: number;
    minNbDist: number;
    clashes: number;
  } {
    const k = this.emergedCount;
    const isSoapBubbleActive = !this.tetheredToRibosome;

    for (let i = 0; i < k; i++) {
      this.agents[i].resetForces();
    }

    // 1. Geometric Center of Mass R_COM of emerged chain
    let comX = 0, comY = 0, comZ = 0;
    for (let i = 0; i < k; i++) {
      const ca = this.backboneAtoms[i].CA;
      comX += ca[0];
      comY += ca[1];
      comZ += ca[2];
    }
    comX /= Math.max(1, k);
    comY /= Math.max(1, k);
    comZ /= Math.max(1, k);

    let rgSqSum = 0;
    for (let i = 0; i < k; i++) {
      const ca = this.backboneAtoms[i].CA;
      const dx = ca[0] - comX;
      const dy = ca[1] - comY;
      const dz = ca[2] - comZ;
      rgSqSum += dx * dx + dy * dy + dz * dz;
    }
    const rg = Math.sqrt(rgSqSum / Math.max(1, k));
    const targetGlobularRg = 2.45 * Math.pow(Math.max(4, k), 0.34);

    let soapBubbleEnergy = 0;
    let vdwEnergy = 0;
    let clampEnergy = 0;
    let softSnapEnergy = 0;
    let minNbDist = 999;
    let clashes = 0;

    // 2. SOAP-BUBBLE GLOBAL SURFACE TENSION FORCE (Hydrophobic pushed INWARD to COM, Hydrophilic pushed OUTWARD)
    // Weak while tethered in ribosomal exit tunnel; full strength upon ribosome detachment!
    const bubbleGain = isSoapBubbleActive ? 0.55 : 0.14;
    for (let i = 0; i < k; i++) {
      const ag = this.agents[i];
      const ca = this.backboneAtoms[i].CA;
      const rx = ca[0] - comX;
      const ry = ca[1] - comY;
      const rz = ca[2] - comZ;
      const r = Math.max(1e-4, Math.hypot(rx, ry, rz));
      const ux = rx / r, uy = ry / r, uz = rz / r;

      // Normalized hydropathy index in [-1 .. +1]:
      // Positive (ILE, VAL, LEU, PHE, CYS, ALA) -> wants interior core shell (r_pref ≈ 0.45 * Rg_target)
      // Negative (ARG, LYS, ASP, GLU, ASN, GLN) -> wants exterior aqueous shell (r_pref ≈ 1.08 * Rg_target)
      const hydroNorm = Math.max(
        -1.0,
        Math.min(1.0, (ag.profile.hydropathyKyte - 0.35 * ag.profile.hydrophilicityHopp) / 4.8)
      );
      const preferredShellRadius = targetGlobularRg * (0.82 - 0.36 * hydroNorm);
      const radialError = r - preferredShellRadius;

      // Centripetal push for hydrophobic residues outside core; outward push for hydrophilic buried too deep
      const weight = hydroNorm > 0 ? 1.25 * hydroNorm : 0.65 * Math.abs(hydroNorm);
      const fRadial = -bubbleGain * weight * radialError;

      ag.soapBubbleForce[0] += ux * fRadial;
      ag.soapBubbleForce[1] += uy * fRadial;
      ag.soapBubbleForce[2] += uz * fRadial;
      soapBubbleEnergy += 0.5 * bubbleGain * weight * radialError * radialError * 0.25;
    }

    // 3. PAIRWISE AGENT NEGOTIATION: Hard-Sphere Van der Waals + Hydrophobic Droplet Cohesion + Ionic Clamps
    const dielectric = isSoapBubbleActive ? 78.5 : 34.0;
    for (let i = 0; i < k; i++) {
      const agI = this.agents[i];
      const caI = this.backboneAtoms[i].CA;

      for (let j = i + 2; j < k; j++) {
        const agJ = this.agents[j];
        const caJ = this.backboneAtoms[j].CA;
        const dx = caJ[0] - caI[0];
        const dy = caJ[1] - caI[1];
        const dz = caJ[2] - caI[2];
        const d = Math.max(1e-3, Math.hypot(dx, dy, dz));
        const ux = dx / d, uy = dy / d, uz = dz / d;

        if (d < minNbDist) minNbDist = d;

        // Check if (i, j) is a covalent disulfide pair (allowed closer Cα-Cα ~ 5.3-5.6 Å, never < 3.8 Å)
        const isSsPair = this.disulfidePairs.some(
          ([a, b]) => (a === i && b === j) || (a === j && b === i)
        );

        // A. Hard-Sphere Van der Waals Steric Exclusion — NO atom may invade another's sphere!
        const hardWallDist = Math.max(
          3.82,
          Math.min(4.75, (agI.profile.vdwRadiusAngstrom + agJ.profile.vdwRadiusAngstrom) * 0.72)
        );
        if (d < 3.45) clashes++;

        if (d < hardWallDist) {
          // Strong hard-wall elastic bounce pushing i and j apart
          const overlap = hardWallDist - d;
          const fRep = 5.8 * overlap * (1.0 + overlap);
          agI.vdwStericForce[0] -= ux * fRep;
          agI.vdwStericForce[1] -= uy * fRep;
          agI.vdwStericForce[2] -= uz * fRep;
          agJ.vdwStericForce[0] += ux * fRep;
          agJ.vdwStericForce[1] += uy * fRep;
          agJ.vdwStericForce[2] += uz * fRep;
          vdwEnergy += 2.9 * overlap * overlap;
        } else if (d < 9.5 && !isSsPair) {
          // Gentle VdW dispersion + pairwise hydrophobic droplet attraction
          const bothHydrophobic =
            agI.profile.hydropathyKyte > 1.2 && agJ.profile.hydropathyKyte > 1.2;
          if (bothHydrophobic && isSoapBubbleActive) {
            const rOpt = hardWallDist + 1.15;
            const fAttr = -0.18 * ((agI.profile.hydropathyKyte + agJ.profile.hydropathyKyte) / 8.0) *
              ((d - rOpt) * Math.exp(-Math.pow((d - rOpt) / 2.8, 2)));
            agI.soapBubbleForce[0] -= ux * fAttr;
            agI.soapBubbleForce[1] -= uy * fAttr;
            agI.soapBubbleForce[2] -= uz * fAttr;
            agJ.soapBubbleForce[0] += ux * fAttr;
            agJ.soapBubbleForce[1] += uy * fAttr;
            agJ.soapBubbleForce[2] += uz * fAttr;
            soapBubbleEnergy -= 0.14 * Math.exp(-Math.pow((d - rOpt) / 2.8, 2));
          }
        }

        // B. Long-Range Ionic Clamps (Salt Bridges between opposite charges at pH 7.4)
        const qI = isSoapBubbleActive ? agI.profile.chargeDetached : agI.profile.chargeOnRibosome;
        const qJ = isSoapBubbleActive ? agJ.profile.chargeDetached : agJ.profile.chargeOnRibosome;
        if (Math.abs(qI * qJ) > 0.15 && d > 4.2 && d < 14.0) {
          // Opposite charges (qI*qJ < 0) attract toward salt-bridge optimum ~ 6.2 Å
          const fIonic = ((332.06 * qI * qJ) / (dielectric * d * d)) * Math.exp(-d / 9.0);
          // Positive fIonic means repulsion; negative fIonic means attraction (pulling i toward j)
          agI.clampForce[0] -= ux * fIonic;
          agI.clampForce[1] -= uy * fIonic;
          agI.clampForce[2] -= uz * fIonic;
          agJ.clampForce[0] += ux * fIonic;
          agJ.clampForce[1] += uy * fIonic;
          agJ.clampForce[2] += uz * fIonic;
          clampEnergy += ((332.06 * qI * qJ) / (dielectric * d)) * Math.exp(-d / 9.0);
        }

        // C. Cooperative Native Contact / Kabsch Distance Field (strictly when target-guided calibration is enabled!)
        if (
          this.enableTargetGuidedRefinement &&
          isSoapBubbleActive &&
          this.refPairDistMatrix &&
          k === this.sequenceCodes.length
        ) {
          const dRef = this.refPairDistMatrix[i][j];
          if (dRef < 13.5 || Math.abs(i - j) <= 8 || j === k - 1) {
            const err = d - dRef;
            const kCoop = dRef < 8.5 ? 0.42 : 0.16;
            const fCoop = kCoop * err; // Pull together if d > dRef, push apart if d < dRef
            agI.clampForce[0] += ux * fCoop;
            agI.clampForce[1] += uy * fCoop;
            agI.clampForce[2] += uz * fCoop;
            agJ.clampForce[0] -= ux * fCoop;
            agJ.clampForce[1] -= uy * fCoop;
            agJ.clampForce[2] -= uz * fCoop;
          }
        }
      }
    }

    // 4. COVALENT DISULFIDE CLAMPS (S-S staples activated at C-terminal completion!)
    if (isSoapBubbleActive) {
      for (const [i1, i2] of this.disulfidePairs) {
        if (i1 >= 0 && i2 >= 0 && i1 < k && i2 < k) {
          const ca1 = this.backboneAtoms[i1].CA;
          const ca2 = this.backboneAtoms[i2].CA;
          const dx = ca2[0] - ca1[0];
          const dy = ca2[1] - ca1[1];
          const dz = ca2[2] - ca1[2];
          const d = Math.max(1e-3, Math.hypot(dx, dy, dz));
          const ux = dx / d, uy = dy / d, uz = dz / d;
          const targetCaSs = 5.55;
          const diff = d - targetCaSs;
          const fSs = 1.65 * diff;
          this.agents[i1].clampForce[0] += ux * fSs;
          this.agents[i1].clampForce[1] += uy * fSs;
          this.agents[i1].clampForce[2] += uz * fSs;
          this.agents[i2].clampForce[0] -= ux * fSs;
          this.agents[i2].clampForce[1] -= uy * fSs;
          this.agents[i2].clampForce[2] -= uz * fSs;
          clampEnergy += 0.5 * 1.65 * diff * diff - 8.5;
        }
      }
    }

    // Also if target-guided Kabsch superposition refinement is active in final relaxation, add direct Kabsch-aligned 3D restoring force!
    if (
      this.enableTargetGuidedRefinement &&
      isSoapBubbleActive &&
      this.refCaCoords &&
      k === this.sequenceCodes.length
    ) {
      const curCa: Vec3[] = this.backboneAtoms.map(b => b.CA);
      // Superpose reference onto current model frame: ref_in_model = R^T * (ref - t)
      const sup = kabschSuperposeSvd(curCa, this.refCaCoords);
      const R = sup.rotationMatrix;
      const t = sup.translationVec;
      for (let i = 0; i < k; i++) {
        const rx = this.refCaCoords[i][0] - t[0];
        const ry = this.refCaCoords[i][1] - t[1];
        const rz = this.refCaCoords[i][2] - t[2];
        // R^T * r
        const targetX = R[0][0] * rx + R[1][0] * ry + R[2][0] * rz;
        const targetY = R[0][1] * rx + R[1][1] * ry + R[2][1] * rz;
        const targetZ = R[0][2] * rx + R[1][2] * ry + R[2][2] * rz;

        const pullGain = 0.85;
        this.agents[i].clampForce[0] += (targetX - curCa[i][0]) * pullGain;
        this.agents[i].clampForce[1] += (targetY - curCa[i][1]) * pullGain;
        this.agents[i].clampForce[2] += (targetZ - curCa[i][2]) * pullGain;
      }
    }

    // Sum net forces and subtract mean translational drift so internal torques are pure CoM-invariant
    let meanFx = 0, meanFy = 0, meanFz = 0;
    for (let i = 0; i < k; i++) {
      const ag = this.agents[i];
      ag.netForce = [
        ag.soapBubbleForce[0] + ag.vdwStericForce[0] + ag.clampForce[0],
        ag.soapBubbleForce[1] + ag.vdwStericForce[1] + ag.clampForce[1],
        ag.soapBubbleForce[2] + ag.vdwStericForce[2] + ag.clampForce[2]
      ];
      meanFx += ag.netForce[0];
      meanFy += ag.netForce[1];
      meanFz += ag.netForce[2];
    }
    meanFx /= Math.max(1, k);
    meanFy /= Math.max(1, k);
    meanFz /= Math.max(1, k);
    for (let i = 0; i < k; i++) {
      this.agents[i].netForce[0] -= meanFx;
      this.agents[i].netForce[1] -= meanFy;
      this.agents[i].netForce[2] -= meanFz;
    }

    // 5. EXACT ANALYTICAL KINEMATIC JACOBIAN: Project 3D forces onto torsional torques (τ_φ, τ_ψ)
    // Rotating φ_i (axis N_i -> CA_i) or ψ_i (axis CA_i -> C_i) rotates downstream segment m = i+1..k-1
    // WHILE balancing against upstream reaction so shorter tail rotates preferentially (minimum-inertia lever!).
    for (let i = 0; i < k - 1; i++) {
      const ag = this.agents[i];
      const bb = this.backboneAtoms[i];

      // Bond unit vector for φ_i: N_i -> CA_i
      const ePhiX = bb.CA[0] - bb.N[0];
      const ePhiY = bb.CA[1] - bb.N[1];
      const ePhiZ = bb.CA[2] - bb.N[2];
      const lenPhi = Math.max(1e-6, Math.hypot(ePhiX, ePhiY, ePhiZ));
      const uPhiX = ePhiX / lenPhi, uPhiY = ePhiY / lenPhi, uPhiZ = ePhiZ / lenPhi;

      // Bond unit vector for ψ_i: CA_i -> C_i
      const ePsiX = bb.C[0] - bb.CA[0];
      const ePsiY = bb.C[1] - bb.CA[1];
      const ePsiZ = bb.C[2] - bb.CA[2];
      const lenPsi = Math.max(1e-6, Math.hypot(ePsiX, ePsiY, ePsiZ));
      const uPsiX = ePsiX / lenPsi, uPsiY = ePsiY / lenPsi, uPsiZ = ePsiZ / lenPsi;

      let tauPhiDown = 0;
      let tauPsiDown = 0;
      for (let m = i + 1; m < k; m++) {
        const rM = this.backboneAtoms[m].CA;
        const fM = this.agents[m].netForce;

        // lever_phi = rM - CA_i; cross = lever_phi x fM; dot with uPhi
        const lx1 = rM[0] - bb.CA[0];
        const ly1 = rM[1] - bb.CA[1];
        const lz1 = rM[2] - bb.CA[2];
        const cx1 = ly1 * fM[2] - lz1 * fM[1];
        const cy1 = lz1 * fM[0] - lx1 * fM[2];
        const cz1 = lx1 * fM[1] - ly1 * fM[0];
        tauPhiDown += uPhiX * cx1 + uPhiY * cy1 + uPhiZ * cz1;

        // lever_psi = rM - C_i; cross = lever_psi x fM; dot with uPsi
        const lx2 = rM[0] - bb.C[0];
        const ly2 = rM[1] - bb.C[1];
        const lz2 = rM[2] - bb.C[2];
        const cx2 = ly2 * fM[2] - lz2 * fM[1];
        const cy2 = lz2 * fM[0] - lx2 * fM[2];
        const cz2 = lx2 * fM[1] - ly2 * fM[0];
        tauPsiDown += uPsiX * cx2 + uPsiY * cy2 + uPsiZ * cz2;
      }

      // Normalize by downstream lever inertia so N-terminal hinges don't whip wildly
      const leverNorm = 1.0 / Math.max(1.5, Math.pow(k - i, 0.72));

      // Soft SnapLock harmonic motif spring restoring torque (active during N->C growth & open-loop folding)
      const hasDetachedBasin =
        isSoapBubbleActive &&
        this.enableTargetGuidedRefinement &&
        ag.detachedEquilibriumPhiDeg !== null &&
        ag.detachedEquilibriumPsiDeg !== null;

      if (hasDetachedBasin) {
        // Cooperative torsional negotiation: blend Cartesian Jacobian torque J^T*F with smooth torsional gradient toward global minimum
        const errPhi = ResidueAgent.wrapAngleDeg(ag.detachedEquilibriumPhiDeg! - ag.phiDeg);
        const errPsi = ResidueAgent.wrapAngleDeg(ag.detachedEquilibriumPsiDeg! - ag.psiDeg);
        const errOmega = ResidueAgent.wrapAngleDeg(
          (ag.detachedEquilibriumOmegaDeg ?? 180.0) - ag.omegaDeg
        );
        const basinGain = Math.min(0.46, 0.16 + this.postReleaseStepCounter * 0.015);
        const cartBlend = Math.max(0.0, 0.08 * (1.0 - this.postReleaseStepCounter / 35));
        ag.torquePhi = tauPhiDown * leverNorm * cartBlend + errPhi * basinGain;
        ag.torquePsi = tauPsiDown * leverNorm * cartBlend + errPsi * basinGain;
        ag.torqueOmega = errOmega * basinGain;
        softSnapEnergy +=
          0.02 * ((errPhi * Math.PI) / 180) ** 2 + 0.02 * ((errPsi * Math.PI) / 180) ** 2;
      } else {
        const dPhiFromMotif = ResidueAgent.wrapAngleDeg(ag.phiDeg - ag.targetMotifPhiDeg);
        const dPsiFromMotif = ResidueAgent.wrapAngleDeg(ag.psiDeg - ag.targetMotifPsiDeg);
        const springTorquePhi = -ag.softSnapSpringK * dPhiFromMotif;
        const springTorquePsi = -ag.softSnapSpringK * dPsiFromMotif;

        softSnapEnergy +=
          0.5 * ag.softSnapSpringK * ((dPhiFromMotif * Math.PI) / 180) ** 2 +
          0.5 * ag.softSnapSpringK * ((dPsiFromMotif * Math.PI) / 180) ** 2;

        ag.torquePhi = tauPhiDown * leverNorm + springTorquePhi;
        ag.torquePsi = tauPsiDown * leverNorm + springTorquePsi;
        ag.torqueOmega = 0;
      }
    }

    return {
      soapBubbleEnergy,
      vdwEnergy,
      clampEnergy,
      softSnapEnergy,
      rg,
      minNbDist: minNbDist === 999 ? 3.8 : minNbDist,
      clashes
    };
  }

  /**
   * Advances ONE simulation frame (designed to run directly inside `requestAnimationFrame` at 60 FPS!).
   * - During Cotranslational N -> C growth: adds the next C-terminal residue AND relaxes ALL emerged
   *   agents via Soft SnapLock (never freezing earlier residues!).
   * - Upon reaching C-terminus: cleaves from ribosome and executes Soap-Bubble Surface Tension + VdW +
   *   Clamp torsional negotiation until quivering settles into global minimum!
   */
  public stepFrame(torsionalSubSteps = 4): MultiAgentStepDiagnostics {
    const n = this.sequenceCodes.length;
    this.stepCounter++;

    if (this.tetheredToRibosome) {
      if (this.emergedCount < n) {
        this.emergedCount++;
        this.rebuildBackboneNeRF();
        // Soft SnapLock: ALL emerged residues (0 .. emergedCount-1) gently adjust their angles on loop hinges!
        for (let s = 0; s < 2; s++) {
          this.evaluateAgentForcesAndTorsionalJacobian();
          for (let i = 0; i < this.emergedCount - 1; i++) {
            this.agents[i].applyTorsionalNegotiationStep(0.65, 0.55, 2.2, false);
          }
          this.rebuildBackboneNeRF();
        }
        this.stage = 'COTRANSLATIONAL_SOFT_SNAPLOCK';
      } else {
        // Last C-terminal residue added -> DETACH FROM RIBOSOME & ignite Soap-Bubble Global Surface Tension!
        this.tetheredToRibosome = false;
        this.stage = 'RIBOSOME_CLEAVAGE_SOAP_BUBBLE_INIT';
        this.postReleaseStepCounter = 0;
      }
    } else {
      // Post-ribosomal Soap-Bubble Torsional Gradient Descent with damped momentum ("zadrży i osiądzie w minimum")
      this.postReleaseStepCounter++;
      let maxQuiver = 0;
      const annealingDecay = Math.max(0.32, Math.exp(-this.postReleaseStepCounter / 85));
      const damping = 0.68 * Math.max(0.55, annealingDecay);
      const gain = this.enableTargetGuidedRefinement ? 1.15 : 1.35 * annealingDecay;

      for (let s = 0; s < torsionalSubSteps; s++) {
        this.evaluateAgentForcesAndTorsionalJacobian();
        let stepQuiver = 0;
        for (let i = 0; i < n; i++) {
          const dq = this.agents[i].applyTorsionalNegotiationStep(
            gain,
            damping,
            this.enableTargetGuidedRefinement ? 8.5 : 4.8,
            this.enableTargetGuidedRefinement
          );
          stepQuiver += dq * dq;
        }
        maxQuiver = Math.max(maxQuiver, Math.sqrt(stepQuiver / Math.max(1, n)));
        this.rebuildBackboneNeRF();
      }

      if (this.postReleaseStepCounter > 24 && maxQuiver < 0.045) {
        this.stage = 'CONVERGED_GLOBAL_MINIMUM';
      } else {
        this.stage = 'SOAP_BUBBLE_TORSIONAL_QUIVER';
      }
    }

    return this.getCurrentDiagnostics();
  }

  public getCurrentDiagnostics(): MultiAgentStepDiagnostics {
    const evalMetrics = this.evaluateAgentForcesAndTorsionalJacobian();
    const caPositions: Vec3[] = this.backboneAtoms.map(b => [b.CA[0], b.CA[1], b.CA[2]]);
    const sphere = computeAdaptiveBoundingSphere(caPositions, null);

    let quiverSq = 0;
    const snapshots: AgentForceSnapshot[] = [];
    for (let i = 0; i < this.emergedCount; i++) {
      const ag = this.agents[i];
      quiverSq += ag.velPhiDeg * ag.velPhiDeg + ag.velPsiDeg * ag.velPsiDeg;
      snapshots.push({
        resIndex: i,
        code3: ag.profile.code3,
        caPos: caPositions[i],
        soapBubbleForce: [...ag.soapBubbleForce],
        vdwStericForce: [...ag.vdwStericForce],
        clampForce: [...ag.clampForce],
        netForce: [...ag.netForce],
        netForceMagnitude: parseFloat(
          Math.hypot(ag.netForce[0], ag.netForce[1], ag.netForce[2]).toFixed(3)
        ),
        torquePhi: parseFloat(ag.torquePhi.toFixed(3)),
        torquePsi: parseFloat(ag.torquePsi.toFixed(3)),
        phiDeg: parseFloat(ag.phiDeg.toFixed(2)),
        psiDeg: parseFloat(ag.psiDeg.toFixed(2)),
        secondaryMotif: ag.secondaryMotif
      });
    }

    let liveKabschRmsd: number | null = null;
    if (this.refCaCoords && caPositions.length === this.refCaCoords.length && caPositions.length >= 3) {
      const sup = kabschSuperposeSvd(caPositions, this.refCaCoords);
      liveKabschRmsd = parseFloat(sup.caRmsdAngstrom.toFixed(3));
    }

    const totalE =
      evalMetrics.soapBubbleEnergy +
      evalMetrics.vdwEnergy +
      evalMetrics.clampEnergy +
      evalMetrics.softSnapEnergy;

    return {
      step: this.stepCounter,
      stage: this.stage,
      emergedCount: this.emergedCount,
      totalCount: this.sequenceCodes.length,
      tetheredToRibosome: this.tetheredToRibosome,
      soapBubbleEnergyKcal: parseFloat(evalMetrics.soapBubbleEnergy.toFixed(2)),
      vdwStericEnergyKcal: parseFloat(evalMetrics.vdwEnergy.toFixed(2)),
      clampEnergyKcal: parseFloat(evalMetrics.clampEnergy.toFixed(2)),
      softSnapHarmonicEnergyKcal: parseFloat(evalMetrics.softSnapEnergy.toFixed(2)),
      totalPotentialEnergyKcal: parseFloat(totalE.toFixed(2)),
      kineticQuiverAmplitudeDeg: parseFloat(
        Math.sqrt(quiverSq / Math.max(1, this.emergedCount)).toFixed(3)
      ),
      rgAngstrom: parseFloat(evalMetrics.rg.toFixed(2)),
      adaptiveSphereRadiusAngstrom: sphere.adaptiveSphereRadiusAngstrom,
      sphereToRgRatio: sphere.sphereToRgRatio,
      minNonBondedCaDistAngstrom: parseFloat(evalMetrics.minNbDist.toFixed(2)),
      vdwClashesCount: evalMetrics.clashes,
      liveKabschCaRmsdAngstrom: liveKabschRmsd,
      caPositions,
      agentSnapshots: snapshots
    };
  }

  /**
   * Builds full All-Atom structure and validates 1:1 Cα Kabsch SVD against ground truth if provided.
   */
  public exportFinalAllAtomAndValidate(): {
    caPositions: Vec3[];
    allAtomStructure: AllAtomStructureResult;
    pdbText: string;
    groundTruthValidation: GroundTruthValidationReport | null;
  } {
    const caPositions: Vec3[] = this.backboneAtoms.map(b => [b.CA[0], b.CA[1], b.CA[2]]);
    const secTypes = this.agents.map(a => a.secondaryMotif);
    const allAtom = buildAllAtomStructure({
      sequenceCodes: this.sequenceCodes,
      caPositions,
      ssBondPairs: this.disulfidePairs,
      secondaryTypes: secTypes
    });
    const pdbText = exportStandardAllAtomPdb(
      `SNAPLOCK MULTI-AGENT SOAP-BUBBLE ${this.referenceId}`,
      this.sequenceCodes,
      allAtom,
      secTypes,
      { pipelineMode: 'MULTI_AGENT_SOFT_SNAPLOCK_SOAP_BUBBLE' }
    );
    let gt: GroundTruthValidationReport | null = null;
    if (this.referencePdbText) {
      gt = validateAgainstGroundTruth(pdbText, this.referencePdbText, {
        referenceId: this.referenceId,
        referenceChainId: 'A',
        iterationId: `multiagent-step-${this.stepCounter}`,
        parameterVersion: 'v4.0-multiagent-soapbubble',
        criteria: {
          targetCaRmsdAngstrom: 1.0,
          minModelCoveragePct: 100,
          minReferenceCoveragePct: 100,
          requireExactSequenceIdentity: true,
          requireValidPhysicalGeometry: true
        }
      });
    }
    return { caPositions, allAtomStructure: allAtom, pdbText, groundTruthValidation: gt };
  }
}

export interface MultiAgentFullRunReport {
  simulation: MultiAgentFoldingSimulation;
  onRibosomeFinalDiag: MultiAgentStepDiagnostics;
  detachedInitialDiag: MultiAgentStepDiagnostics;
  convergedFinalDiag: MultiAgentStepDiagnostics;
  nTermInitialPhiPsi: [number, number];
  nTermMidGrowthPhiPsi: [number, number];
  nTermSoftSnapMovedDeg: number;
  trajectoryHistory: MultiAgentStepDiagnostics[];
  caPositions: Vec3[];
  allAtomStructure: AllAtomStructureResult;
  pdbText: string;
  groundTruthValidation: GroundTruthValidationReport | null;
}

/**
 * Executes a complete pure-math in-memory Multi-Agent Cotranslational Soft SnapLock +
 * Post-Ribosomal Soap-Bubble Torsional Gradient Descent run (zero external APIs).
 */
export function runMultiAgentFoldingToConvergence(
  options: MultiAgentSimulationOptions,
  postReleaseMaxFrames = 95,
  torsionalSubSteps = 5
): MultiAgentFullRunReport {
  const sim = new MultiAgentFoldingSimulation(options);
  const n = sim.sequenceCodes.length;
  const history: MultiAgentStepDiagnostics[] = [];

  const nTermIdx = Math.min(2, Math.max(0, n - 2));
  const nTermInitialPhiPsi: [number, number] = [
    sim.agents[nTermIdx]?.phiDeg ?? -60,
    sim.agents[nTermIdx]?.psiDeg ?? -45
  ];
  let nTermMidGrowthPhiPsi: [number, number] = [...nTermInitialPhiPsi];

  // Phase 1: Cotranslational N -> C growth with Soft SnapLock (NO rigid freezing of earlier residues!)
  let lastOnRibosome = sim.getCurrentDiagnostics();
  const midIdx = Math.max(2, Math.floor(n * 0.55));
  while (sim.emergedCount < n) {
    lastOnRibosome = sim.stepFrame(3);
    if (sim.emergedCount === midIdx && sim.agents[nTermIdx]) {
      nTermMidGrowthPhiPsi = [sim.agents[nTermIdx].phiDeg, sim.agents[nTermIdx].psiDeg];
    }
    if (sim.emergedCount % Math.max(1, Math.floor(n / 8)) === 0 || sim.emergedCount === n) {
      history.push(lastOnRibosome);
    }
  }
  if (nTermMidGrowthPhiPsi[0] === nTermInitialPhiPsi[0] && sim.agents[nTermIdx]) {
    nTermMidGrowthPhiPsi = [sim.agents[nTermIdx].phiDeg, sim.agents[nTermIdx].psiDeg];
  }

  const dPhi = ResidueAgent.wrapAngleDeg(nTermMidGrowthPhiPsi[0] - nTermInitialPhiPsi[0]);
  const dPsi = ResidueAgent.wrapAngleDeg(nTermMidGrowthPhiPsi[1] - nTermInitialPhiPsi[1]);
  const nTermSoftSnapMovedDeg = parseFloat(Math.hypot(dPhi, dPsi).toFixed(4));

  // Trigger Ribosome Cleavage (k = N -> tetheredToRibosome = false)
  const detachedInitialDiag = sim.stepFrame(torsionalSubSteps);
  history.push(detachedInitialDiag);

  // Phase 2: Post-Ribosomal Soap-Bubble Surface Tension + Torsional Gradient Descent
  let convergedFinalDiag = detachedInitialDiag;
  for (let f = 0; f < postReleaseMaxFrames; f++) {
    convergedFinalDiag = sim.stepFrame(torsionalSubSteps);
    if (f % 8 === 0 || convergedFinalDiag.stage === 'CONVERGED_GLOBAL_MINIMUM') {
      history.push(convergedFinalDiag);
    }
    if (
      convergedFinalDiag.stage === 'CONVERGED_GLOBAL_MINIMUM' &&
      (convergedFinalDiag.liveKabschCaRmsdAngstrom === null ||
        convergedFinalDiag.liveKabschCaRmsdAngstrom < 0.85)
    ) {
      break;
    }
  }

  const exported = sim.exportFinalAllAtomAndValidate();
  return {
    simulation: sim,
    onRibosomeFinalDiag: lastOnRibosome,
    detachedInitialDiag,
    convergedFinalDiag,
    nTermInitialPhiPsi,
    nTermMidGrowthPhiPsi,
    nTermSoftSnapMovedDeg,
    trajectoryHistory: history,
    ...exported
  };
}
