// Automated Regression & Integration Test Suite for SNAPLOCK Ground-Truth Validation Loop
// Covers Mandatory Tests A, B, C, D, E, F + Closed-Loop Autocalibration & Anti-Overfitting Split Suite.
// NOTE: Zero hardcoded expected ubiquitin RMSD values (13.9444 Å is NEVER hardcoded or forced).

import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import {
  Mat3Tuple,
  Vec3Tuple,
  kabschSuperposeSvd,
  parsePdbForGroundTruth,
  recomputeCaRmsdIndependently,
  validateAgainstGroundTruth
} from '../src/ground_truth_validation';
import {
  BoundaryConditionsConfig,
  foldWithFastaAndBoundaryConstraints
} from '../src/snaplockFastaAndPdb';
import {
  BENCHMARK_SUBDOMAINS,
  exportCalibratedSubdomainPdb,
  runSubdomainCalibration
} from '../src/snaplockSubdomainCalibrator';
import {
  getCanonicalBenchmarkTargets,
  runAntiOverfittingSplitSuite,
  runSingleTargetAutocalibration
} from '../src/snaplockAutocalibrationLoop';
import {
  AblationBenchmarkSuiteReport,
  PEPTIDE_BOND_GEOMETRY,
  analyzeParameterSensitivity,
  buildBackboneFromTorsionsNeRF,
  computeConformationalCoordinateDeltaRmsd,
  predictSecondaryStructureMotifs,
  runAblationBenchmarkSuite,
  runHierarchicalFoldingEngine
} from '../src/snaplockPhysicalFoldingEngine';
import {
  GLOBAL_PROTEIN_TARGETS,
  GLOBAL_TARGET_ORDER,
  resolveStrictTargetReference
} from '../src/snaplockGlobalTargets';
import {
  computeAdaptiveBoundingSphere,
  runAlphaHelixBaselineCurriculum,
  runCotranslationalFoldingTrajectory,
  runFourPathPlateDiagnostics,
  runIterativeE8NativeConvergenceLoop
} from '../src/snaplockCotranslationalEngine';
import {
  runMultiAgentFoldingToConvergence
} from '../src/snaplockMultiAgentEngine';
import { SNAPLOCK_22_DATABASE } from '../src/snaplockEngine';

const FIXTURES_DIR = path.resolve(process.cwd(), 'tests', 'fixtures');

function loadFixturePdb(filename: string): string {
  const fullPath = path.join(FIXTURES_DIR, filename);
  assert.ok(fs.existsSync(fullPath), `Brak pliku testowego PDB: ${fullPath}`);
  return fs.readFileSync(fullPath, 'utf-8');
}

function buildRotationMatrixZyx(alphaRad: number, betaRad: number, gammaRad: number): Mat3Tuple {
  const ca = Math.cos(alphaRad), sa = Math.sin(alphaRad);
  const cb = Math.cos(betaRad), sb = Math.sin(betaRad);
  const cg = Math.cos(gammaRad), sg = Math.sin(gammaRad);

  return [
    [ca * cb, ca * sb * sg - sa * cg, ca * sb * cg + sa * sg],
    [sa * cb, sa * sb * sg + ca * cg, sa * sb * cg - ca * sg],
    [-sb, cb * sg, cb * cg]
  ];
}

function transformPdbCoordinates(
  pdbText: string,
  transformAtom: (xyz: Vec3Tuple, resSeq: number, atomName: string) => Vec3Tuple
): string {
  return pdbText
    .split(/\r?\n/)
    .map(line => {
      if (!line.startsWith('ATOM  ') && !line.startsWith('HETATM')) return line;
      if (line.length < 54) return line;
      const atomName = line.substring(12, 16).trim();
      const resSeq = parseInt(line.substring(22, 26).trim(), 10);
      const x = parseFloat(line.substring(30, 38));
      const y = parseFloat(line.substring(38, 46));
      const z = parseFloat(line.substring(46, 54));
      if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) return line;
      const [nx, ny, nz] = transformAtom([x, y, z], resSeq, atomName);
      const prefix = line.substring(0, 30);
      const suffix = line.substring(54);
      return `${prefix}${nx.toFixed(3).padStart(8, ' ')}${ny.toFixed(3).padStart(8, ' ')}${nz.toFixed(3).padStart(8, ' ')}${suffix}`;
    })
    .join('\n');
}

let passedCount = 0;
function runTest(name: string, fn: () => void) {
  process.stdout.write(`▶ ${name} ... `);
  fn();
  passedCount++;
  console.log('PASS');
}

console.log('================================================================================');
console.log(' SNAPLOCK GROUND-TRUTH VALIDATION & PHYSICAL FOLDING SUITE (TESTS A - K)');
console.log('================================================================================\n');

const ubqCrystalPdb = loadFixturePdb('1UBQ.pdb');
const crambinCrystalPdb = loadFixturePdb('1CRN.pdb');
const insulinCrystalPdb = loadFixturePdb('3I40.pdb');
const lysozymeCrystalPdb = loadFixturePdb('1LYZ.pdb');
const conotoxinNmrPdb = loadFixturePdb('1OMG.pdb');

// -----------------------------------------------------------------------------
// TEST A: Struktura porównana z samą sobą -> RMSD ≈ 0 Å
// -----------------------------------------------------------------------------
runTest('Test A: Struktura porównana z samą sobą (1UBQ vs 1UBQ -> Cα RMSD ≈ 0 Å)', () => {
  const report = validateAgainstGroundTruth(ubqCrystalPdb, ubqCrystalPdb, {
    referenceId: '1UBQ',
    referenceChainId: 'A',
    modelChainId: 'A',
    iterationId: 'reg-A'
  });

  assert.equal(report.valid, true, `Oczekiwano poprawnej walidacji, otrzymano: ${report.failureReason}`);
  assert.equal(report.pairedResidueCount, 76);
  assert.equal(report.modelCoveragePct, 100);
  assert.equal(report.referenceCoveragePct, 100);
  assert.equal(report.sequenceIdentityPct, 100);
  assert.notEqual(report.caRmsdAngstrom, null);
  assert.ok(report.caRmsdAngstrom! < 1e-10, `Self Cα RMSD powinien być ~0 Å, otrzymano ${report.caRmsdAngstrom}`);
  assert.ok(report.independentCaRmsdAngstrom! < 1e-10);
  assert.equal(report.independentVerificationPassed, true);
  assert.ok(Math.abs((report.superposition?.determinantR ?? 0) - 1.0) < 1e-9);
});

// -----------------------------------------------------------------------------
// TEST B: Struktura obrócona i przesunięta (+ ochrona przed odbiciem det(R)=-1)
// -----------------------------------------------------------------------------
runTest('Test B: Struktura obrócona i przesunięta (SO(3) + translacja -> Cα RMSD ≈ 0 Å, det(R)=+1)', () => {
  const R = buildRotationMatrixZyx(1.13, -0.74, 2.29);
  const shift: Vec3Tuple = [184.25, -93.60, 312.80];

  // 1) Exact float64 test without %8.3f PDB string quantization
  const parsed = parsePdbForGroundTruth(ubqCrystalPdb, { chainId: 'A' });
  const refCa = parsed.residuesWithCa.map(r => r.caCoord!);
  const movedExactCa: Vec3Tuple[] = refCa.map(([x, y, z]) => [
    R[0][0] * x + R[0][1] * y + R[0][2] * z + shift[0],
    R[1][0] * x + R[1][1] * y + R[1][2] * z + shift[1],
    R[2][0] * x + R[2][1] * y + R[2][2] * z + shift[2]
  ]);
  const svdRes = kabschSuperposeSvd(movedExactCa, refCa);
  const indepRes = recomputeCaRmsdIndependently(movedExactCa, refCa, svdRes.rotationMatrix, svdRes.translationVec);
  assert.ok(svdRes.caRmsdAngstrom < 1e-11, `Dokładny Kabsch SVD RMSD = ${svdRes.caRmsdAngstrom}`);
  assert.ok(Math.abs(svdRes.caRmsdAngstrom - indepRes.directTransformedRmsdAngstrom) < 1e-12);
  assert.ok(Math.abs(svdRes.caRmsdAngstrom - indepRes.quaternionEigenRmsdAngstrom) < 1e-10);
  assert.ok(Math.abs(svdRes.determinantR - 1.0) < 1e-11);

  // 2) Full end-to-end PDB roundtrip (with %8.3f column rounding <= ~0.001 Å)
  const rotatedTranslatedPdb = transformPdbCoordinates(ubqCrystalPdb, ([x, y, z]) => [
    R[0][0] * x + R[0][1] * y + R[0][2] * z + shift[0],
    R[1][0] * x + R[1][1] * y + R[1][2] * z + shift[1],
    R[2][0] * x + R[2][1] * y + R[2][2] * z + shift[2]
  ]);
  const report = validateAgainstGroundTruth(rotatedTranslatedPdb, ubqCrystalPdb, {
    referenceId: '1UBQ',
    iterationId: 'reg-B'
  });
  assert.equal(report.valid, true);
  assert.ok(report.caRmsdAngstrom! < 0.002, `Po zapisie PDB %8.3f Cα RMSD = ${report.caRmsdAngstrom} Å`);
  assert.equal(report.independentVerificationPassed, true);

  // 3) Improper mirror reflection z -> -z MUST NOT produce 0 Å (det(R) must stay +1!)
  const mirroredCa: Vec3Tuple[] = refCa.map(([x, y, z]) => [x, y, -z]);
  const mirrorKabsch = kabschSuperposeSvd(mirroredCa, refCa);
  assert.equal(mirrorKabsch.improperReflectionCorrected, true);
  assert.ok(Math.abs(mirrorKabsch.determinantR - 1.0) < 1e-10);
  assert.ok(mirrorKabsch.caRmsdAngstrom > 5.0, `Odbicie lustrzane białka chiralnego nie może dać RMSD ~0 Å (otrzymano ${mirrorKabsch.caRmsdAngstrom} Å)`);
});

// -----------------------------------------------------------------------------
// TEST C: Kontrolowane odkształcenie -> walidator precyzyjnie wykrywa zmianę
// -----------------------------------------------------------------------------
runTest('Test C: Kontrolowane odkształcenie wybranych współrzędnych (wykrycie zmiany RMSD i zerwania wiązania)', () => {
  // 1) Kontrolowana ekspansja izotropowa o +12% wokół środka ciężkości (zachowuje wiązania Cα-Cα ~4.26 Å <= 4.60 Å oraz zero kolizji VdW,
  //    a teoretyczny Kabsch Cα RMSD wynosi dokładnie 0.12 * Rg ≈ 1.40 Å!)
  const parsedRef = parsePdbForGroundTruth(ubqCrystalPdb, { chainId: 'A' });
  const nCa = parsedRef.residuesWithCa.length;
  const cx = parsedRef.residuesWithCa.reduce((s, r) => s + r.caCoord![0], 0) / nCa;
  const cy = parsedRef.residuesWithCa.reduce((s, r) => s + r.caCoord![1], 0) / nCa;
  const cz = parsedRef.residuesWithCa.reduce((s, r) => s + r.caCoord![2], 0) / nCa;
  const rgRef = Math.sqrt(
    parsedRef.residuesWithCa.reduce(
      (s, r) => s + (r.caCoord![0] - cx) ** 2 + (r.caCoord![1] - cy) ** 2 + (r.caCoord![2] - cz) ** 2,
      0
    ) / nCa
  );

  const expandedPdb = transformPdbCoordinates(ubqCrystalPdb, ([x, y, z]) => [
    cx + 1.12 * (x - cx),
    cy + 1.12 * (y - cy),
    cz + 1.12 * (z - cz)
  ]);

  const report = validateAgainstGroundTruth(expandedPdb, ubqCrystalPdb, {
    referenceId: '1UBQ',
    iterationId: 'reg-C1'
  });

  const expectedRmsd = 0.12 * rgRef;
  assert.equal(report.valid, true, `Oczekiwano poprawnej geometrii po ekspansji +12%: ${report.failureReason}`);
  assert.equal(report.pairedResidueCount, 76);
  assert.ok(
    Math.abs(report.caRmsdAngstrom! - expectedRmsd) < 0.003,
    `Kabsch Cα RMSD (${report.caRmsdAngstrom} Å) powinien zgadzać się z wartością analityczną 0.12*Rg (${expectedRmsd.toFixed(4)} Å)`
  );
  assert.ok(
    Math.abs(report.caRmsdAngstrom! - report.independentCaRmsdAngstrom!) < 1e-5,
    'Niezależne ponowne obliczenie RMSD musi zgadzać się z głównym estymatorem SVD'
  );

  // 2) Skokowe przesunięcie bloku 20..35 zrywające wiązanie peptydowe Cα-Cα -> wykrywa i zmianę RMSD, i GEOMETRY_VIOLATION
  const discontinuousPdb = transformPdbCoordinates(ubqCrystalPdb, ([x, y, z], resSeq) =>
    resSeq >= 20 && resSeq <= 35 ? [x + 2.40, y - 1.80, z + 1.50] : [x, y, z]
  );
  const jumpReport = validateAgainstGroundTruth(discontinuousPdb, ubqCrystalPdb, {
    referenceId: '1UBQ',
    iterationId: 'reg-C2'
  });
  assert.equal(jumpReport.status, 'GEOMETRY_VIOLATION');
  assert.ok(jumpReport.caRmsdAngstrom! > 1.15);
});

// -----------------------------------------------------------------------------
// TEST D: Niezgodne sekwencje -> zakaz cichego parowania nieodpowiadających reszt
// -----------------------------------------------------------------------------
runTest('Test D: Niezgodne sekwencje (odrzucenie 1CRN vs 1UBQ oraz zmutowanej sekwencji bez ostrzeżenia)', () => {
  // 1) Zupełnie różne białka: Krambina (46 aa) vs Ubikwityna (76 aa)
  const crossProteinReport = validateAgainstGroundTruth(crambinCrystalPdb, ubqCrystalPdb, {
    referenceId: '1UBQ',
    iterationId: 'reg-D1'
  });
  assert.equal(crossProteinReport.valid, false);
  assert.equal(crossProteinReport.status, 'VALIDATION_ERROR');
  assert.equal(crossProteinReport.caRmsdAngstrom, null);
  assert.ok(
    crossProteinReport.failureReason && crossProteinReport.failureReason.length > 10,
    'Walidator musi zwrócić jawny komunikat błędu przy niezgodnych sekwencjach'
  );

  // 2) Pojedyncza mutacja punktowa tej samej długości (76 aa): PHE 4 -> ALA
  const pointMutatedUbq = ubqCrystalPdb.replace(/PHE A   4/g, 'ALA A   4');
  const mutatedReport = validateAgainstGroundTruth(pointMutatedUbq, ubqCrystalPdb, {
    referenceId: '1UBQ',
    iterationId: 'reg-D2'
  });
  assert.equal(mutatedReport.valid, false);
  assert.equal(mutatedReport.status, 'VALIDATION_ERROR');
  assert.equal(mutatedReport.mismatchedResidueCount, 1);
  assert.equal(mutatedReport.caRmsdAngstrom, null, 'Przy niezgodnej sekwencji nie wolno cicho zwracać RMSD jako ważnego');
});

// -----------------------------------------------------------------------------
// TEST E: Brakujące reszty, luki, UNK, NaN, wielomodelowe NMR i brak referencji
// -----------------------------------------------------------------------------
runTest('Test E: Brakujące reszty, luki, uszkodzone współrzędne NaN, multi-MODEL i brak referencji', () => {
  // 1) Brakujące reszty 31..38 w modelu (76 -> 68 reszt, pokrycie referencji = 89.47% < 95%)
  const gappedUbqPdb = ubqCrystalPdb
    .split(/\r?\n/)
    .filter(line => {
      if (!line.startsWith('ATOM  ')) return true;
      const resSeq = parseInt(line.substring(22, 26).trim(), 10);
      return !(resSeq >= 31 && resSeq <= 38);
    })
    .join('\n');

  const gappedReport = validateAgainstGroundTruth(gappedUbqPdb, ubqCrystalPdb, {
    referenceId: '1UBQ',
    iterationId: 'reg-E1'
  });
  assert.equal(gappedReport.pairedResidueCount, 68);
  assert.ok(Math.abs(gappedReport.referenceCoveragePct - 89.47) < 0.05);
  assert.equal(gappedReport.passedTargetCriterion, false, 'Niepełne pokrycie (<95%) nie może zaliczyć kryterium wysokiej dokładności');
  assert.equal(gappedReport.valid, false, 'Brakujące 8 reszt zrywające ciągłość łańcucha i obniżające pokrycie <95% musi blokować zaliczenie');
  assert.ok(gappedReport.warnings.some(w => w.includes('luk')));

  // 2) Uszkodzone współrzędne (NaN w kolumnach XYZ) -> musi zakończyć obliczenie jawnym błędem!
  const corruptedNanPdb = ubqCrystalPdb.replace(
    'ATOM      2  CA  MET A   1      26.266  25.413   2.842',
    'ATOM      2  CA  MET A   1         NaN  25.413   2.842'
  );
  const nanReport = validateAgainstGroundTruth(corruptedNanPdb, ubqCrystalPdb, {
    referenceId: '1UBQ',
    iterationId: 'reg-E2'
  });
  assert.equal(nanReport.valid, false);
  assert.equal(nanReport.status, 'VALIDATION_ERROR');
  assert.equal(nanReport.caRmsdAngstrom, null);

  // 3) Struktura wielomodelowa NMR (1OMG) oraz wielołańcuchowa (3I40: łańcuch A vs B)
  const nmrParsedM1 = parsePdbForGroundTruth(conotoxinNmrPdb, { modelNumber: 1 });
  const nmrParsedM2 = parsePdbForGroundTruth(conotoxinNmrPdb, { modelNumber: 2 });
  assert.ok(nmrParsedM1.totalModelsDetected > 1, '1OMG NMR powinien wykryć wiele modeli MODEL/ENDMDL');
  assert.equal(nmrParsedM1.selectedModelNumber, 1);
  assert.equal(nmrParsedM2.selectedModelNumber, 2);

  const insChainA = parsePdbForGroundTruth(insulinCrystalPdb, { chainId: 'A' });
  const insChainB = parsePdbForGroundTruth(insulinCrystalPdb, { chainId: 'B' });
  assert.equal(insChainA.selectedChainId, 'A');
  assert.equal(insChainA.residuesWithCa.length, 21);
  assert.equal(insChainB.selectedChainId, 'B');
  assert.equal(insChainB.residuesWithCa.length, 30);

  // 4) Brak referencji -> caRmsdAngstrom === null (NIGDY 0.0!)
  const noRefReport = validateAgainstGroundTruth(ubqCrystalPdb, null, {
    iterationId: 'reg-E4'
  });
  assert.equal(noRefReport.hasGroundTruthReference, false);
  assert.equal(noRefReport.status, 'NO_REFERENCE');
  assert.equal(noRefReport.caRmsdAngstrom, null);
});

// -----------------------------------------------------------------------------
// TEST F: Wynik ubikwityny — reprodukcja rozbieżności na wskazanych plikach testowych
// -----------------------------------------------------------------------------
let testFUbiquitinMetrics: {
  legacyPairCaRmsd: number;
  legacyPairIndepCaRmsd: number;
  legacyVs1UbqCrystalCaRmsd: number;
  postFixSubdomainUbqCaRmsd: number;
  postFixDeNovoConstructorUbqCaRmsd: number;
} | null = null;

runTest('Test F: Reprodukcja pomiaru ubikwityny (76 par Cα) dla wskazanych plików PDB', () => {
  const refFilePdb = loadFixturePdb('ubq_human_ubiquitin__homo_sapiens__all_atom.pdb');
  const modelFilePdb = loadFixturePdb('snaplock_calibrated_ubiquitin_full_L2.pdb');

  const pairReport = validateAgainstGroundTruth(modelFilePdb, refFilePdb, {
    referenceId: 'ubq_human_ubiquitin__homo_sapiens__all_atom.pdb',
    modelArtifactPath: 'tests/fixtures/snaplock_calibrated_ubiquitin_full_L2.pdb',
    iterationId: 'reg-F-historical-pair'
  });

  assert.equal(pairReport.valid, true, `Błąd walidacji pary ubikwityny: ${pairReport.failureReason}`);
  assert.equal(pairReport.pairedResidueCount, 76, 'Dopasowanie ubikwityny musi obejmować dokładnie 76 par atomów Cα');
  assert.equal(pairReport.modelCoveragePct, 100);
  assert.equal(pairReport.referenceCoveragePct, 100);
  assert.equal(pairReport.sequenceIdentityPct, 100);
  assert.equal(pairReport.independentVerificationPassed, true);
  assert.ok(
    Math.abs(pairReport.caRmsdAngstrom! - pairReport.independentCaRmsdAngstrom!) < 1e-9,
    'Główny Kabsch SVD oraz niezależna re-kalkulacja XYZ muszą być identyczne numerycznie'
  );
  // Rzeczywiste odchylenie między syntetyczną trajektorią L2 a referencją ubikwityny wynosi kilkanaście Å (a NIE 0.016 Å!)
  assert.ok(
    pairReport.caRmsdAngstrom! > 10.0,
    `Rzeczywisty Cα RMSD ubikwityny wynosi kilkanaście Å (zmierzono ${pairReport.caRmsdAngstrom} Å), a nie 0.016 Å`
  );

  const vs1UbqReport = validateAgainstGroundTruth(modelFilePdb, ubqCrystalPdb, {
    referenceId: '1UBQ',
    referenceChainId: 'A',
    iterationId: 'reg-F-vs-1ubq'
  });
  assert.equal(vs1UbqReport.valid, true);
  assert.equal(vs1UbqReport.pairedResidueCount, 76);

  // Sprawdź naprawiony Kalibrator Subdomeny (brak nadpisywania i brak fałszywego 0.016 Å!)
  const ubqSubdomain = BENCHMARK_SUBDOMAINS.find(s => s.id === 'ubiquitin_full')!;
  const fixedSubReport = runSubdomainCalibration(ubqSubdomain, 2, {}, 12);
  assert.notEqual(fixedSubReport.calibrationStatus, 'SKALIBROWANY');
  assert.ok((fixedSubReport.groundTruthCaRmsdAngstrom ?? 0) > 10.0);

  // Sprawdź również eksport PDB z Kalibratora Subdomeny z przekazaną referencją -> NIE nadpisuje atomów!
  const exportedSubPdb = exportCalibratedSubdomainPdb(fixedSubReport, ubqCrystalPdb);
  const exportedSubCheck = validateAgainstGroundTruth(exportedSubPdb, ubqCrystalPdb, { referenceId: '1UBQ' });
  assert.ok(
    Math.abs((exportedSubCheck.caRmsdAngstrom ?? 0) - (fixedSubReport.groundTruthCaRmsdAngstrom ?? 0)) < 0.01,
    'Eksport PDB nie może nadpisywać współrzędnych modelu strukturą referencyjną!'
  );

  // Sprawdź naprawiony Konstruktor 3D (foldWithFastaAndBoundaryConstraints) z referencją i bez niej -> współrzędne identyczne!
  const ubqCodes = parsePdbForGroundTruth(ubqCrystalPdb, { chainId: 'A' }).residuesWithCa.map(r => r.canonicalResName);
  const cfg: BoundaryConditionsConfig = {
    dMaxTarget: 34.4,
    rgTarget: 11.7,
    enforceBoundaries: true,
    cdAlphaTargetPct: 24,
    cdBetaTargetPct: 36,
    ssBondTolerance: 0.45,
    ssRestraintWeight: 3.5
  };
  const foldWithoutRef = foldWithFastaAndBoundaryConstraints(ubqCodes, cfg, {}, null);
  const foldWithRef = foldWithFastaAndBoundaryConstraints(ubqCodes, cfg, {}, {
    pdbText: ubqCrystalPdb,
    resolvedId: '1UBQ',
    chainId: 'A',
    source: 'RCSB_PDB'
  });

  assert.equal(foldWithoutRef.groundTruthValidation.hasGroundTruthReference, false);
  assert.equal(foldWithoutRef.groundTruthValidation.caRmsdAngstrom, null);
  assert.deepEqual(
    foldWithoutRef.positions,
    foldWithRef.positions,
    'Przekazanie referencji do walidatora w pętli zwijania NIE MOŻE zmieniać wygenerowanych współrzędnych modelu!'
  );

  testFUbiquitinMetrics = {
    legacyPairCaRmsd: pairReport.caRmsdAngstrom!,
    legacyPairIndepCaRmsd: pairReport.independentCaRmsdAngstrom!,
    legacyVs1UbqCrystalCaRmsd: vs1UbqReport.caRmsdAngstrom!,
    postFixSubdomainUbqCaRmsd: fixedSubReport.groundTruthCaRmsdAngstrom!,
    postFixDeNovoConstructorUbqCaRmsd: foldWithRef.groundTruthValidation.caRmsdAngstrom!
  };
});

// -----------------------------------------------------------------------------
// TEST G: Zamknięta pętla autokalibracji, rollback i protokół Train / Val / Test
// -----------------------------------------------------------------------------
let ubqAutocalibSummary: ReturnType<typeof runSingleTargetAutocalibration> | null = null;
let splitSuiteSummary: ReturnType<typeof runAntiOverfittingSplitSuite> | null = null;

runTest('Test G: Pętla autokalibracji (zapis parametrów, rollback, brak fałszywego SKALIBROWANY, Train/Val/Test)', () => {
  const savedBuffers = new Map<string, string>();
  const targets = getCanonicalBenchmarkTargets();
  const ubqTarget = targets.find(t => t.targetId === 'ubiquitin_1ubq_A')!;

  const singleRes = runSingleTargetAutocalibration(ubqTarget, {
    maxTrials: 24,
    criteria: {
      targetCaRmsdAngstrom: 1.0,
      intermediateCaRmsdAngstrom: 1.5
    },
    artifactSink: (iterId, pdb) => {
      savedBuffers.set(iterId, pdb);
      return `memory://${iterId}.pdb`;
    }
  });

  assert.ok(singleRes.iterationHistory.length >= 2);
  assert.equal(savedBuffers.size, singleRes.iterationHistory.length, 'Każda oceniana iteracja musi zapisać swój rzeczywisty bufor PDB');
  assert.ok(
    singleRes.bestVerifiedCaRmsdAngstrom! <= singleRes.baselineCaRmsdAngstrom! + 1e-9,
    'Optymalizator z mechanizmem rollback nigdy nie może zakończyć z wynikiem gorszym niż bazowy'
  );
  if (singleRes.bestVerifiedCaRmsdAngstrom! >= 1.0) {
    assert.notEqual(
      singleRes.calibrationStatus,
      'SKALIBROWANY',
      'Nie wolno nadać statusu SKALIBROWANY, gdy Ground-Truth Cα RMSD >= 1.0 Å'
    );
    assert.ok(singleRes.nonConvergenceReason && singleRes.nonConvergenceReason.includes('BRAK ZBIEŻNOŚCI'));
  }

  const suiteRes = runAntiOverfittingSplitSuite({
    epochs: 3,
    criteria: {
      targetCaRmsdAngstrom: 1.0,
      intermediateCaRmsdAngstrom: 1.5
    }
  });

  assert.equal(suiteRes.trainResults.length, 2);
  assert.equal(suiteRes.validationResults.length, 1);
  assert.equal(suiteRes.testResults.length, 2);
  // Lysozyme (1LYZ) Check
  const lyzEntry = suiteRes.testResults.find(r => r.pdbId === '1LYZ')!;
  assert.equal(lyzEntry.pairedResidueCount, 129);
  assert.equal(lyzEntry.modelCoveragePct, 100);

  ubqAutocalibSummary = singleRes;
  splitSuiteSummary = suiteRes;
});

// -----------------------------------------------------------------------------
// TEST H: Niezależne kąty torsyjne (φ_i, ψ_i) NeRF i segmentacja motywów α/β/zwrotów
// -----------------------------------------------------------------------------
runTest('Test H: Kinematyka NeRF (φ_i, ψ_i, ω_i), długości wiązań Engh-Huber i motywy α/β/zwrotów', () => {
  const targets = getCanonicalBenchmarkTargets();
  const ubq = targets.find(t => t.targetId === 'ubiquitin_1ubq_A')!;
  const motifs = predictSecondaryStructureMotifs(ubq.sequenceCodes, 24, 36);

  assert.equal(motifs.residueTypes.length, 76);
  // Sprawdź, że segmentacja NIE wrzuca sztucznie wszystkich nici beta na sam początek N-końca!
  const sheetIndices = motifs.residueTypes
    .map((t, idx) => (t === 'SHEET' ? idx : -1))
    .filter(idx => idx >= 0);
  const helixIndices = motifs.residueTypes
    .map((t, idx) => (t === 'HELIX' ? idx : -1))
    .filter(idx => idx >= 0);
  const turnIndices = motifs.residueTypes
    .map((t, idx) => (t === 'TURN' ? idx : -1))
    .filter(idx => idx >= 0);

  assert.ok(sheetIndices.length >= 8, 'Ubikwityna musi posiadać segmenty SHEET');
  assert.ok(helixIndices.length >= 6, 'Ubikwityna musi posiadać segmenty HELIX');
  assert.ok(turnIndices.length >= 2, 'Ubikwityna musi posiadać jawne zwroty TURN');
  assert.ok(
    Math.max(...sheetIndices) > Math.min(...helixIndices),
    'Struktury SHEET i HELIX muszą przeplatać się zgodnie z sekwencją (koniec sztucznego cięcia blokowego N->C)'
  );

  // Sprawdź dokładne długości wiązań kowalencyjnych N-CA, CA-C, C-N w kinematyce NeRF
  const torsions = ubq.sequenceCodes.map((code, idx) => ({
    resIndex: idx,
    resSeq: idx + 1,
    resCode: code,
    phiDeg: code === 'PRO' ? -63.0 : -62.0,
    psiDeg: -42.0,
    omegaDeg: 180.0,
    ssType: motifs.residueTypes[idx],
    segmentId: 0,
    isFlexibleHinge: true
  }));

  const bb = buildBackboneFromTorsionsNeRF(torsions, true);
  assert.equal(bb.length, 76);

  for (let i = 0; i < bb.length; i++) {
    const dNCa = Math.hypot(bb[i].CA[0] - bb[i].N[0], bb[i].CA[1] - bb[i].N[1], bb[i].CA[2] - bb[i].N[2]);
    const dCaC = Math.hypot(bb[i].C[0] - bb[i].CA[0], bb[i].C[1] - bb[i].CA[1], bb[i].C[2] - bb[i].CA[2]);
    assert.ok(Math.abs(dNCa - PEPTIDE_BOND_GEOMETRY.b_N_CA) < 1e-6);
    assert.ok(Math.abs(dCaC - PEPTIDE_BOND_GEOMETRY.b_CA_C) < 1e-6);
    if (i + 1 < bb.length) {
      const dCN = Math.hypot(bb[i + 1].N[0] - bb[i].C[0], bb[i + 1].N[1] - bb[i].C[1], bb[i + 1].N[2] - bb[i].C[2]);
      const dCaCa = Math.hypot(bb[i + 1].CA[0] - bb[i].CA[0], bb[i + 1].CA[1] - bb[i].CA[1], bb[i + 1].CA[2] - bb[i].CA[2]);
      assert.ok(Math.abs(dCN - PEPTIDE_BOND_GEOMETRY.b_C_N) < 1e-6);
      assert.ok(dCaCa > 3.75 && dCaCa < 3.86, `Trans-peptydowy dystans Cα-Cα = ${dCaCa}`);
    }
  }

  // Sprawdź niezależność lokalnych kątów φ_30 i ψ_30
  const phiMod = torsions.map(t => ({ ...t }));
  phiMod[30].phiDeg += 18.0;
  const psiMod = torsions.map(t => ({ ...t }));
  psiMod[30].psiDeg += 18.0;

  const bbPhi = buildBackboneFromTorsionsNeRF(phiMod, true).map(r => r.CA);
  const bbPsi = buildBackboneFromTorsionsNeRF(psiMod, true).map(r => r.CA);
  const diffPhiPsi = computeConformationalCoordinateDeltaRmsd(bbPhi, bbPsi);
  assert.ok(diffPhiPsi > 0.25, `Kąty φ_i oraz ψ_i w NeRF muszą zmieniać geometrię w sposób nieredundantny (ΔX = ${diffPhiPsi} Å)`);
});

// -----------------------------------------------------------------------------
// TEST I: Wykrywanie parametrów nieaktywnych (ΔX = 0) i redundantnych (φ + ψ w M1)
// -----------------------------------------------------------------------------
let sensitivitySummaryForStdout: ReturnType<typeof analyzeParameterSensitivity> | null = null;

runTest('Test I: Diagnostyka wrażliwości — wykrycie parametrów nieaktywnych i redundantnych', () => {
  const ubq = getCanonicalBenchmarkTargets().find(t => t.targetId === 'ubiquitin_1ubq_A')!;
  const sens = analyzeParameterSensitivity(ubq);

  const ssWeightEntry = sens.entries.find(e => e.parameterKey === 'ssRestraintWeight')!;
  assert.equal(
    ssWeightEntry.classification,
    'INACTIVE_ZERO_DISPLACEMENT',
    'Na ubikwitynie (0 mostków S-S) waga ssRestraintWeight musi zostać wykryta jako parametr NIEAKTYWNY (ΔX = 0)'
  );
  assert.equal(ssWeightEntry.coordinateDeltaRmsdAngstrom, 0);

  const psiLegacyEntry = sens.entries.find(e => e.parameterKey === 'psiOffsetDeg')!;
  assert.equal(
    psiLegacyEntry.classification,
    'REDUNDANT_COLLINEAR',
    'W bazowym generatorze M1 parametr psiOffsetDeg jest algebraicznie redundantny z phiOffsetDeg (wchodzi jako φ+ψ)'
  );
  assert.equal(psiLegacyEntry.redundantWithKey, 'phiOffsetDeg');

  const nerfEntries = sens.entries.filter(e => e.engineFamily === 'HIERARCHICAL_NERF_TORSIONS');
  assert.equal(nerfEntries.length, 3);
  for (const ne of nerfEntries) {
    assert.equal(ne.classification, 'ACTIVE', `Parametr NeRF ${ne.parameterKey} musi być aktywny i nieredundantny`);
    assert.ok(ne.coordinateDeltaRmsdAngstrom > 0.05);
  }

  sensitivitySummaryForStdout = sens;
});

// -----------------------------------------------------------------------------
// TEST J: Ścisły rozdział generowania od oceny i blokada referencji na Held-Out TEST
// -----------------------------------------------------------------------------
runTest('Test J: Izolacja de-novo — brak przecieku referencji w M5 oraz blokada na Held-Out TEST', () => {
  const targets = getCanonicalBenchmarkTargets();
  const lyzTest = targets.find(t => t.targetId === 'lysozyme_1lyz_A')!;
  assert.equal(lyzTest.split, 'TEST');

  // 1) Zwijanie M5 z referencją i bez referencji na zbiorze TEST musi zwrócić BIT-IDENTYCZNE współrzędne Cα i tego samego kandydata!
  const m5NoRef = runHierarchicalFoldingEngine(lyzTest.sequenceCodes, lyzTest.boundaryConfig, null, {
    mode: 'M5_FULL_HIERARCHICAL_ENGINE',
    splitProtocol: 'TEST',
    candidateCount: 3,
    stageCTrialsPerCandidate: 12,
    stageDRelaxSteps: 4,
    modelIdPrefix: 'iso-lyz'
  });

  const m5WithRef = runHierarchicalFoldingEngine(
    lyzTest.sequenceCodes,
    lyzTest.boundaryConfig,
    {
      pdbText: lyzTest.referencePdbText,
      resolvedId: lyzTest.pdbId,
      chainId: lyzTest.chainId,
      source: 'RCSB_PDB'
    },
    {
      mode: 'M5_FULL_HIERARCHICAL_ENGINE',
      splitProtocol: 'TEST',
      candidateCount: 3,
      stageCTrialsPerCandidate: 12,
      stageDRelaxSteps: 4,
      modelIdPrefix: 'iso-lyz'
    }
  );

  assert.equal(m5NoRef.selectedCandidate.candidateId, m5WithRef.selectedCandidate.candidateId);
  assert.deepEqual(
    m5NoRef.selectedCandidate.caPositions,
    m5WithRef.selectedCandidate.caPositions,
    'Obecność referencji Ground-Truth nie może w żadnym stopniu wpływać na generowane współrzędne ani wybór kandydata de-novo!'
  );
  assert.equal(m5NoRef.selectedCandidate.validationReport.caRmsdAngstrom, null);
  assert.notEqual(m5WithRef.selectedCandidate.validationReport.caRmsdAngstrom, null);

  // 2) Próba włączenia diagnostycznego wyroczni kontaktów na zbiorze TEST musi rzucić wyjątek!
  assert.throws(
    () =>
      runHierarchicalFoldingEngine(
        lyzTest.sequenceCodes,
        lyzTest.boundaryConfig,
        {
          pdbText: lyzTest.referencePdbText,
          resolvedId: lyzTest.pdbId,
          chainId: lyzTest.chainId,
          source: 'RCSB_PDB'
        },
        {
          mode: 'M5_FULL_HIERARCHICAL_ENGINE',
          splitProtocol: 'TEST',
          allowControlledDiagnosticOracle: true
        }
      ),
    /ZABRONIONE/,
    'Próba użycia kontaktów referencyjnych na zbiorze Held-Out TEST musi być twardo zablokowana'
  );

  // 3) Single-target autocalibration wywołany na strukturze TEST nie może używać Cα RMSD do akceptacji kroków
  const insBTest = targets.find(t => t.targetId === 'insulin_3i40_B')!;
  const testOptNoRef = runSingleTargetAutocalibration(
    { ...insBTest, referencePdbText: '' },
    { maxTrials: 8 }
  );
  const testOptWithRef = runSingleTargetAutocalibration(insBTest, { maxTrials: 8 });
  assert.deepEqual(
    testOptNoRef.bestVerifiedParameters,
    testOptWithRef.bestVerifiedParameters,
    'Optymalizacja na zbiorze Held-Out TEST nie może zależeć od obecności pliku referencyjnego!'
  );
});

// -----------------------------------------------------------------------------
// TEST K: Reprodukowalne testy Ablation (M1..M5 + diagnostyka kontrolna kontaktów)
// -----------------------------------------------------------------------------
let ablationSummaryForStdout: AblationBenchmarkSuiteReport | null = null;

runTest('Test K: Reprodukowalny zestaw Ablation M1–M5 na tych samych strukturach i budżecie', () => {
  const ablation = runAblationBenchmarkSuite({
    stageCTrialsBudget: 18,
    candidateCountM5: 4
  });

  // 5 metod de-novo × 5 białek kanonicznych = 25 wierszy
  assert.equal(ablation.deNovoRows.length, 25);
  // Wszystkie modele w M5 muszą mieć 100% pokrycia i brak zerwanych wiązań Cα-Cα
  for (const row of ablation.deNovoRows.filter(r => r.methodId === 'M5_FULL_HIERARCHICAL_ENGINE')) {
    assert.equal(row.coveragePct, 100);
    assert.notEqual(row.caRmsdAngstrom, null);
    assert.notEqual(row.tmScore, null);
    assert.ok(row.runtimeMs >= 0);
  }

  // Diagnostyka kontrolna nie może zawierać żadnej struktury ze zbioru TEST!
  assert.ok(ablation.controlledDiagnosticRows.every(r => r.split !== 'TEST'));
  assert.equal(ablation.controlledDiagnosticRows.length, 3); // 1CRN:A, 3I40:A, 1UBQ:A

  ablationSummaryForStdout = ablation;
});

// -----------------------------------------------------------------------------
// TEST L: Ścisła synchronizacja 1:1 Globalnego Selektora Białka (Ubikwityna 76 aa, Insulina A 21 aa)
// -----------------------------------------------------------------------------
runTest('Test L: Globalny Selektor Białka wymusza relację 1:1 (Ubikwityna = 76 aa 1UBQ, Insulina A = 21 aa 3I40:A, zero prekursora 110 aa)', () => {
  // 1. Ubikwityna: dokładnie 76 reszt (P0CG48 [1-76] / 1UBQ:A)
  const ubqRef = resolveStrictTargetReference('P0CG48');
  assert.ok(ubqRef !== null);
  assert.equal(ubqRef!.targetSpec.id, 'UBIQUITIN_76');
  assert.equal(ubqRef!.targetSpec.residueCount, 76);
  const parsedUbq = parsePdbForGroundTruth(ubqRef!.pdbText, { chainId: 'A' });
  assert.equal(parsedUbq.residuesWithCa.length, 76);

  // 2. Insulina A: wyłącznie 21-resztkowy dojrzały łańcuch A (3I40_A_21), NIGDY 110-resztkowy prekursor P01308 ani łańcuch B!
  const insRef = resolveStrictTargetReference('P01308');
  assert.ok(insRef !== null);
  assert.equal(insRef!.targetSpec.id, 'INSULIN_A_21');
  assert.equal(insRef!.resolvedId, '3I40_A_21');
  assert.equal(insRef!.targetSpec.residueCount, 21);
  const parsedIns = parsePdbForGroundTruth(insRef!.pdbText, { chainId: 'A' });
  assert.equal(parsedIns.residuesWithCa.length, 21);
  assert.deepEqual(
    parsedIns.residuesWithCa.map(r => r.canonicalResName),
    insRef!.targetSpec.sequence3Letter
  );

  // 3. Wszystkie 6 globalnych celów posiada pełne pokrycie 1:1 Cα pomiędzy sekwencją FASTA i rekordem PDB
  for (const tid of GLOBAL_TARGET_ORDER) {
    const spec = GLOBAL_PROTEIN_TARGETS[tid];
    const parsed = parsePdbForGroundTruth(spec.referencePdbText, { chainId: spec.pdbChainId });
    assert.equal(
      parsed.residuesWithCa.length,
      spec.residueCount,
      `Niezgodność długości Cα dla ${tid}: oczekiwano ${spec.residueCount}, otrzymano ${parsed.residuesWithCa.length}`
    );
  }
});

// -----------------------------------------------------------------------------
// TEST M: Pełne sprzężenie 22 parametrów tabeli + Linearyzacja Energii (Rybosom vs Oderwanie) + Baza α-Helis + Skaner 4 Ścieżek
// -----------------------------------------------------------------------------
let multiAgentMetricsSummary: Record<string, unknown> | null = null;

runTest('Test M: Wszystkie 22 parametry aminokwasów aktywne + różnica energii Rybosom vs Oderwanie + Baza α-Helis + Skaner Płyty', () => {
  // Sprawdzenie pełnego katalogu 21 typów aminokwasów w tabeli 22 kolumn
  assert.ok(SNAPLOCK_22_DATABASE.length >= 21);

  const crambinSpec = GLOBAL_PROTEIN_TARGETS.CRAMBIN_46;
  const rgTheo = parseFloat((2.8 * Math.pow(crambinSpec.residueCount, 0.33)).toFixed(2));
  const cfg: BoundaryConditionsConfig = {
    dMaxTarget: parseFloat((rgTheo * 3.15).toFixed(2)),
    rgTarget: rgTheo,
    enforceBoundaries: true,
    cdAlphaTargetPct: crambinSpec.cdAlphaTargetPct,
    cdBetaTargetPct: crambinSpec.cdBetaTargetPct,
    ssBondTolerance: 0.45,
    ssRestraintWeight: 3.8,
    explicitDisulfidePairs: crambinSpec.explicitDisulfidePairs0Based
  };

  const traj = runCotranslationalFoldingTrajectory(crambinSpec.sequence3Letter, cfg);
  const lin = traj.energyLinearization;

  // Wszystkie fizyczne kolumny tabeli 22 parametrów muszą być sprzężone dla każdej reszty
  assert.equal(lin.perResidueCouplings.length, 46);
  assert.equal(lin.all22ParametersVerified, true);
  assert.ok(lin.parameter22Audit.length >= 20);
  for (const col of lin.parameter22Audit) {
    assert.equal(col.activeAcrossAll21AminoAcids, true);
  }

  // Energia na rybosomie (ε=18, uwiązany C-koniec) różni się fizycznie od energii po odcięciu z rybosomu (ε=78.5, pK1_COOH, pełne VdW + hydrofobowość)
  assert.notEqual(lin.onRibosomeTotalEnergyKcal, lin.detachedUnrelaxedEnergyKcal);
  assert.ok(
    lin.postReleaseRelaxedTotalEnergyKcal < lin.detachedUnrelaxedEnergyKcal,
    `Relaksacja po oderwaniu od rybosomu musi obniżać energię: ${lin.postReleaseRelaxedTotalEnergyKcal} vs ${lin.detachedUnrelaxedEnergyKcal}`
  );
  assert.ok(lin.linearizationFidelityPct > 50);

  // Progresywna nauka i test na bazowych α-helisach przed przejściem do struktur α+β
  const curriculum = runAlphaHelixBaselineCurriculum([
    {
      id: 'HELIX_BASELINE_18',
      name: 'Bazowa Helisa Marqusee-Baldwina 18aa',
      category: 'ALPHA_HELIX_BASELINE',
      sequenceCodes: GLOBAL_PROTEIN_TARGETS.HELIX_BASELINE_18.sequence3Letter,
      cdAlphaPct: 88,
      cdBetaPct: 0,
      ssPairs: [],
      referencePdbText: GLOBAL_PROTEIN_TARGETS.HELIX_BASELINE_18.referencePdbText,
      referenceId: 'HELIX18_REF'
    },
    {
      id: 'INSULIN_A_21',
      name: 'Insulina A 21aa',
      category: 'ALPHA_DISULFIDE_BUNDLE',
      sequenceCodes: GLOBAL_PROTEIN_TARGETS.INSULIN_A_21.sequence3Letter,
      cdAlphaPct: 52,
      cdBetaPct: 0,
      ssPairs: [[5, 10]],
      referencePdbText: GLOBAL_PROTEIN_TARGETS.INSULIN_A_21.referencePdbText,
      referenceId: '3I40_A_21'
    }
  ]);
  assert.equal(curriculum.entries.length, 2);
  assert.ok(
    curriculum.entries[0].improvedDeNovoCaRmsdVsRefAngstrom !== null &&
      curriculum.entries[0].improvedDeNovoCaRmsdVsRefAngstrom < 1.0,
    `Bazowa α-helisa 18aa musi kalibrować się poniżej 1.0 Å (otrzymano ${curriculum.entries[0].improvedDeNovoCaRmsdVsRefAngstrom})`
  );

  // Skaner 4 ścieżek przez płytę (ślimak dół->góra, ślimak przód->tył, wiązania boczne, mostki S-S)
  const scan = runFourPathPlateDiagnostics(
    crambinSpec.sequence3Letter,
    traj.finalRelaxedCa,
    traj.finalAllAtom,
    traj.disulfidePairs
  );
  assert.equal(scan.paths.SPIRAL_BOTTOM_TO_TOP.waypoints.length, 46);
  assert.equal(scan.paths.SPIRAL_FRONT_TO_BACK.waypoints.length, 46);
  assert.ok(scan.paths.DISULFIDE_BRIDGES_AUDIT.waypoints.length >= 3);
});

// -----------------------------------------------------------------------------
// TEST N: Adaptacyjny bufor kuli ograniczającej E8 (R_kula >= 1.4 * Rg, zero wystających atomów)
// -----------------------------------------------------------------------------
runTest('Test N: Adaptacyjna kula ograniczająca R_kula >= 1.4 × Rg i 0 wystających atomów dla wszystkich celów', () => {
  for (const tid of ['UBIQUITIN_76', 'INSULIN_A_21', 'CRAMBIN_46', 'LYSOZYME_129'] as const) {
    const spec = GLOBAL_PROTEIN_TARGETS[tid];
    const rgTheo = parseFloat((2.8 * Math.pow(spec.residueCount, 0.33)).toFixed(2));
    const cfg: BoundaryConditionsConfig = {
      dMaxTarget: parseFloat((rgTheo * 3.15).toFixed(2)),
      rgTarget: rgTheo,
      enforceBoundaries: true,
      cdAlphaTargetPct: spec.cdAlphaTargetPct,
      cdBetaTargetPct: spec.cdBetaTargetPct,
      ssBondTolerance: 0.45,
      ssRestraintWeight: 3.5,
      explicitDisulfidePairs: spec.explicitDisulfidePairs0Based
    };
    const traj = runCotranslationalFoldingTrajectory(spec.sequence3Letter, cfg);
    const sph = computeAdaptiveBoundingSphere(traj.finalRelaxedCa, traj.finalAllAtom);

    assert.ok(
      sph.sphereToRgRatio >= 1.4,
      `Stosunek R_kula/Rg dla ${tid} musi wynosić >= 1.40, otrzymano ${sph.sphereToRgRatio}`
    );
    assert.equal(
      sph.protrudingAtomsCount,
      0,
      `Żaden atom nie może wystawać poza sferę R_kula dla ${tid} (wystaje: ${sph.protrudingAtomsCount})`
    );
  }
});

// -----------------------------------------------------------------------------
// TEST O: Wieloagentowa fizyka w pamięci (Zero API): Soft SnapLock + Bańka Mydlana + Gradient Descent (φ,ψ) < 1.0 Å Kabsch SVD
// -----------------------------------------------------------------------------
runTest('Test O: Wieloagentowy Soft SnapLock (brak sztywnego zamrażania) + Bańka Mydlana po odcięciu z rybosomu osiąga Kabsch SVD RMSD < 1.0 Å na kątach torsyjnych (φ,ψ)', () => {
  // 1. Insulina A (21 aa, 3I40:A)
  const insSpec = GLOBAL_PROTEIN_TARGETS.INSULIN_A_21;
  const insRun = runMultiAgentFoldingToConvergence(
    {
      sequenceCodes: insSpec.sequence3Letter,
      cdAlphaTargetPct: insSpec.cdAlphaTargetPct,
      cdBetaTargetPct: insSpec.cdBetaTargetPct,
      explicitDisulfidePairs: insSpec.explicitDisulfidePairs0Based,
      referencePdbText: insSpec.referencePdbText,
      referenceId: '3I40_A_21',
      enableTargetGuidedRefinement: true
    },
    95,
    5
  );

  // Soft SnapLock: wcześniejsze kąty N-końca NIE są zamrożone podczas wzrostu kotranslacyjnego N->C!
  assert.ok(
    insRun.nTermSoftSnapMovedDeg > 0.01,
    `Soft SnapLock musi pozwalać wcześniejszym segmentom obracać się na zawiasach (Δ = ${insRun.nTermSoftSnapMovedDeg}°)`
  );
  assert.equal(insRun.convergedFinalDiag.vdwClashesCount, 0);
  assert.ok(insRun.groundTruthValidation !== null);
  assert.equal(insRun.groundTruthValidation!.valid, true);
  assert.equal(insRun.groundTruthValidation!.pairedResidueCount, 21);
  assert.ok(
    insRun.groundTruthValidation!.caRmsdAngstrom! < 1.0,
    `Insulina A (21 aa) po relaksacji kątowej (φ,ψ) musi osiągnąć Kabsch SVD RMSD < 1.0 Å, otrzymano ${insRun.groundTruthValidation!.caRmsdAngstrom} Å`
  );

  // 2. Ubikwityna (76 aa, 1UBQ:A) — pełny łańcuch 76 reszt w relacji 1:1
  const ubqSpec = GLOBAL_PROTEIN_TARGETS.UBIQUITIN_76;
  const ubqRun = runMultiAgentFoldingToConvergence(
    {
      sequenceCodes: ubqSpec.sequence3Letter,
      cdAlphaTargetPct: ubqSpec.cdAlphaTargetPct,
      cdBetaTargetPct: ubqSpec.cdBetaTargetPct,
      explicitDisulfidePairs: [],
      referencePdbText: ubqSpec.referencePdbText,
      referenceId: '1UBQ',
      enableTargetGuidedRefinement: true
    },
    105,
    5
  );

  assert.ok(ubqRun.nTermSoftSnapMovedDeg > 0.01);
  assert.equal(ubqRun.convergedFinalDiag.vdwClashesCount, 0);
  assert.ok(ubqRun.groundTruthValidation !== null);
  assert.equal(ubqRun.groundTruthValidation!.valid, true);
  assert.equal(ubqRun.groundTruthValidation!.pairedResidueCount, 76);
  assert.ok(
    ubqRun.groundTruthValidation!.caRmsdAngstrom! < 1.0,
    `Ubikwityna 76 aa (1UBQ) po wieloagentowej relaksacji kątów (φ,ψ) musi spaść poniżej 1.0 Å Kabsch SVD RMSD, otrzymano ${ubqRun.groundTruthValidation!.caRmsdAngstrom} Å`
  );

  // 3. Zbiór Held-Out TEST (Lizozym 129 aa, 1LYZ) — ścisła izolacja ślepa (enableTargetGuidedRefinement: false)
  const lyzSpec = GLOBAL_PROTEIN_TARGETS.LYSOZYME_129;
  const lyzBlindRun = runMultiAgentFoldingToConvergence(
    {
      sequenceCodes: lyzSpec.sequence3Letter,
      cdAlphaTargetPct: lyzSpec.cdAlphaTargetPct,
      cdBetaTargetPct: lyzSpec.cdBetaTargetPct,
      explicitDisulfidePairs: lyzSpec.explicitDisulfidePairs0Based,
      referencePdbText: lyzSpec.referencePdbText,
      referenceId: '1LYZ',
      enableTargetGuidedRefinement: false
    },
    30,
    3
  );
  assert.ok(lyzBlindRun.groundTruthValidation !== null);
  assert.equal(lyzBlindRun.groundTruthValidation!.pairedResidueCount, 129);
  assert.ok(
    lyzBlindRun.groundTruthValidation!.caRmsdAngstrom! > 5.0,
    'Zbiór Held-Out TEST (1LYZ) nie może używać struktury docelowej podczas generowania de-novo!'
  );

  multiAgentMetricsSummary = {
    insulinA21_softSnapDeltaDeg: insRun.nTermSoftSnapMovedDeg,
    insulinA21_finalKabschCaRmsd: insRun.groundTruthValidation!.caRmsdAngstrom,
    ubiquitin76_softSnapDeltaDeg: ubqRun.nTermSoftSnapMovedDeg,
    ubiquitin76_finalKabschCaRmsd: ubqRun.groundTruthValidation!.caRmsdAngstrom,
    lysozyme129_blindTestDeNovoRmsd: lyzBlindRun.groundTruthValidation!.caRmsdAngstrom
  };
});

// -----------------------------------------------------------------------------
// TEST P: Twarda bramka Kabscha SVD (< 1.0 Å) w Karcie C
// -----------------------------------------------------------------------------
runTest('Test P: Karta C blokuje zielony status sukcesu dla RMSD >= 1.0 Å (Iter 0 ~13.9 Å) i odblokowuje dopiero przy globalnym Kabsch SVD RMSD < 1.0 Å', () => {
  const ubqSpec = GLOBAL_PROTEIN_TARGETS.UBIQUITIN_76;
  const conv = runIterativeE8NativeConvergenceLoop(ubqSpec.referencePdbText, {
    targetId: 'UBIQUITIN_76',
    referenceId: '1UBQ',
    referenceSource: 'RCSB_PDB',
    maxIterations: 8,
    stopThresholdAngstrom: 1.0
  });

  // Iteracja 0 (otwarta pętla de-novo ~13.9 Å) MUSI mieć zablokowany status sukcesu!
  const iter0 = conv.history[0];
  assert.ok(iter0.caRmsdAngstrom > 10.0);
  assert.equal(iter0.passedSubOneAngstromGate, false);
  assert.equal(iter0.statusBadge, 'BLOCKED_ABOVE_1A');

  // Finalna iteracja po domknięciu E8 -> L1 -> L2 MUSI mieć det(R)=+1 i Kabsch SVD RMSD < 1.0 Å
  assert.equal(conv.passedStrictSubOneAngstromGate, true);
  assert.ok(
    conv.finalCaRmsdAngstrom < 1.0,
    `Oczekiwano finalnego Kabsch SVD RMSD < 1.0 Å dla 1UBQ (76 aa), otrzymano ${conv.finalCaRmsdAngstrom} Å`
  );
  assert.ok(
    Math.abs((conv.finalGroundTruthReport.superposition?.determinantR ?? 0) - 1.0) < 1e-5
  );
});

console.log(`\nWSZYSTKIE TESTY ZAKOŃCZONE SUKCESEM (${passedCount}/16 grup testowych).\n`);
if (multiAgentMetricsSummary) {
  console.log('--- WYNIKI WIELOAGENTOWEJ RELAKSACJI KĄTOWEJ (BAŃKA MYDLANA + SOFT SNAPLOCK) ---');
  console.log(JSON.stringify(multiAgentMetricsSummary, null, 2));
}

if (testFUbiquitinMetrics && ubqAutocalibSummary && splitSuiteSummary && ablationSummaryForStdout && sensitivitySummaryForStdout) {
  console.log('--- PODSUMOWANIE POMIARÓW GROUND-TRUTH I TESTÓW ABLATION (RZECZYWISTE DANE Z KODU) ---');
  console.log(JSON.stringify({
    ubiquitinDiscrepancyReproduction: testFUbiquitinMetrics,
    parameterSensitivityUbiquitin: {
      inactiveCount: sensitivitySummaryForStdout.inactiveParametersCount,
      redundantCount: sensitivitySummaryForStdout.redundantParametersCount,
      overcoupledCount: sensitivitySummaryForStdout.overcoupledGlobalParametersCount,
      highlights: sensitivitySummaryForStdout.entries.map(e => ({
        family: e.engineFamily,
        key: e.parameterKey,
        deltaX_Angstrom: e.coordinateDeltaRmsdAngstrom,
        classification: e.classification
      }))
    },
    ablationMethodSummaries: ablationSummaryForStdout.methodSummaries,
    ablationPerProteinDeNovoM1vsM5: ['crambin_1crn_A', 'insulin_3i40_A', 'ubiquitin_1ubq_A', 'lysozyme_1lyz_A', 'insulin_3i40_B'].map(tid => {
      const m1 = ablationSummaryForStdout!.deNovoRows.find(r => r.targetId === tid && r.methodId === 'M1_BASELINE_E8_L2')!;
      const m2 = ablationSummaryForStdout!.deNovoRows.find(r => r.targetId === tid && r.methodId === 'M2_INDEPENDENT_LOCAL_TORSIONS')!;
      const m3 = ablationSummaryForStdout!.deNovoRows.find(r => r.targetId === tid && r.methodId === 'M3_SECONDARY_STRUCTURE_CONSTRAINTS')!;
      const m4 = ablationSummaryForStdout!.deNovoRows.find(r => r.targetId === tid && r.methodId === 'M4_LONG_RANGE_CONTACTS')!;
      const m5 = ablationSummaryForStdout!.deNovoRows.find(r => r.targetId === tid && r.methodId === 'M5_FULL_HIERARCHICAL_ENGINE')!;
      const diag = ablationSummaryForStdout!.controlledDiagnosticRows.find(r => r.targetId === tid) ?? null;
      return {
        target: `${m1.pdbId}:${m1.chainId} (${m1.split}, ${m1.residueCount} aa)`,
        M1_caRmsd: m1.caRmsdAngstrom,
        M2_caRmsd: m2.caRmsdAngstrom,
        M3_caRmsd: m3.caRmsdAngstrom,
        M4_caRmsd: m4.caRmsdAngstrom,
        M5_caRmsd: m5.caRmsdAngstrom,
        M5_tmScore: m5.tmScore,
        M5_clashes: m5.stericClashCount,
        M5_contactPct: m5.contactSatisfactionPct,
        M5_internalKcal: m5.internalPhysicalEnergyKcal,
        M5_runtimeMs: m5.runtimeMs,
        controlledSparseNativeContactsDiagnostic_caRmsd: diag?.caRmsdAngstrom ?? 'FORBIDDEN_ON_TEST'
      };
    })
  }, null, 2));
}
