
# Protein Folding Cosmic Synthesizer
### Spatial Bio-Acoustic Probe & Folding Engine (WIP)

> ⚠️ **Project Status & Engineering Notice (Current Sprint)**  
> **UI, 3D Spatial Rendering & Multi-track Sonification Engine:** Fully functional and interactive.  
> **Core Biophysical Folding Engine:** Currently undergoing active refactoring.  
> 
> The initial prototype utilized heuristic/geometric anchors (SNAPLOCK model) to visualize target conformations (e.g. native 7NHU insulin, ubiquitin). However, benchmark verification revealed that the ab-initio folding minimization algorithm requires a fundamental overhaul to prevent artificial energy fitting and unphysical backbone rigidity. We are actively decoupling hardcoded PDB constraints and rebuilding the physics-driven conformational search engine from the ground up (implementing true Ramachandran torsion optimization, realistic steric clash avoidance, and accurate solvent-accessible surface area potentials).

---

## 🛠️ Roadmap & Current Work (Engine Overhaul)

- [x] **Interactive WebGL / 3D Canvas:** Real-time residue manipulation, center-of-mass tracking, and spatial axis projections.
- [x] **Harmonic Sonification Matrix:** Web Audio API pipeline converting structural coordinates into dynamic tonal palettes and chord progressions.
- [ ] **[IN PROGRESS] Physics Engine Overhaul:** Transitioning from discrete geometric locking to a continuous gradient/torsion angle optimization engine.
- [ ] **[IN PROGRESS] Unbiased Energy Scoring:** Replacing preset cohesion values with empirical force-field potentials (AMBER/CHARMM-compatible terms).
- [ ] **[PLANNED] Rigorous PDB/AlphaFold2 Validation:** Publishing genuine, unconstrained RMSD metrics across benchmark sets (7NHU, 1UBQ, 1LYZ).
