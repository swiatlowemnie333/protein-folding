# SNAPLOCK - protein structure engine with an honest validation harness

**Author:** Michał Zazuniuk (concepts, design, validation). Code written with AI coding agents, reviewed and run by the author.

SNAPLOCK is a TypeScript/React app for building, visualising and **rigorously validating** 3D protein backbones.
Its most valuable part is the evaluation harness:

- Kabsch (SVD) superposition with independent quaternion cross-check, Ca / backbone / all-heavy-atom RMSD, TM-score
- NeRF backbone construction, all-atom builder, PDB import/export, steric-clash checks
- **TRAIN / VALIDATION / held-out TEST split**; on TEST any use of native coordinates is blocked (`FORBIDDEN_ON_TEST`)
- **Ablation study M1-M5** isolating each engine component

## Result: de-novo prediction does not work (yet) - reported honestly

Ca RMSD in Angstrom (lower is better; < 2 A is a good prediction). TM-score < 0.17 = random.

| Target | M1 baseline | M2 local torsions | M3 + sec. structure | M4 + long-range | **M5 full** | TM (M5) |
|---|---|---|---|---|---|---|
| Insulin A, 21 aa (TRAIN) | 9.0 | 7.4 | 7.9 | 6.0 | **5.8** | 0.02 |
| Ubiquitin 1UBQ, 76 aa (VALIDATION) | 10.3 | 22.0 | 34.6 | 14.1 | **23.3** | 0.04 |
| **Lysozyme 1LYZ, 129 aa (TEST)** | 15.5 | 24.5 | 44.1 | 12.2 | **27.4** | 0.05 |
| Insulin B, 30 aa (TEST) | 6.1 | 9.4 | 14.8 | 7.6 | **11.6** | 0.03 |

1. All TM-scores are in the random range: the de-novo engine does not predict the native fold.
2. The full engine (M5) beats the baseline only on its TRAIN target and is worse on VALIDATION/TEST: overfitting, exposed by the split and the ablation.
3. Reproduce: `npm install && npm test`

## What the "< 1 A" numbers in the UI are
Some panels show Ca RMSD < 1 A. They come from **reconstruction loops that compare against the known reference at every iteration** (e.g. `runIterativeE8NativeConvergenceLoop`, open-loop start ~13.9 A). This is re-encoding a known structure, **not prediction**. `snaplock_alphafold_benchmark.py` likewise measures discretisation error of AlphaFold structures.

## Naming note
The "248 vectors" are a **248-direction spherical (Fibonacci) discretisation of S2**, named after the dimension of the E8 Lie algebra - not the E8 root system (240 roots in 8D).

## Next step: real predictor in the same harness
`tests/esmfold_benchmark.ts` sends the same sequences to **ESMFold** (Meta, public API) and scores them with the same Kabsch/TM-score code. Run: `npx tsx tests/esmfold_benchmark.ts`

## Run locally
    npm install
    npm run dev     # app on http://localhost:3000
    npm test        # validation + ablation suite

## ESMFold in the same harness (Nov 2026)
Same targets, same Kabsch / TM-score code (`tests/esmfold_benchmark.ts`, raw PDBs in `results/esmfold/`).

| Target | ESMFold Ca RMSD | SNAPLOCK M5 | ESMFold TM-score | mean pLDDT |
|---|---|---|---|---|
| **Lysozyme 1LYZ, 129 aa (TEST)** | **0.55 A** | 27.4 A | **0.983** | 95.1 |
| Ubiquitin 1UBQ, 76 aa (VALIDATION) | **0.83 A** | 23.3 A | 0.954 | 90.5 |
| Insulin B, 30 aa (TEST) | 2.30 A | 11.6 A | 0.459 | 77.9 |
| Insulin A, 21 aa (TRAIN) | 0.87 A | 5.8 A | 0.330 | 85.7 |
| Crambin 1CRN, 46 aa (TRAIN) | 7.86 A | - | 0.138 | 50.3 |

- The harness reproduces near-atomic accuracy for a state-of-the-art predictor, which **validates the harness itself**.
- ESMFold confidence (pLDDT) tracks real error: the one miss (crambin) is flagged by pLDDT ~50.
- Isolated insulin chains score low on TM: in vivo A and B fold together via disulfide bonds; single-chain input lacks the partner.
- Conclusion: SNAPLOCK's value is the **validation/visualisation layer**; prediction should come from established models.
