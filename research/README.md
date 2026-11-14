# research/ — polowanie na reguły (uczciwe wyniki)

Stan 2026-11-14. Dane: struktury z GitHuba (RCSB niedostępne z sandboxa), 142 niepowtarzalne białka, 20 635 reszt.
- `scripts/gly_scan.py`, `gly_motif.py` — glicyna (phi>0 u 63% Gly vs 2% reszty; kontekst: C-cap helisy 93%).
- `scripts/scan20.py` → `an20_seq.py` — wszystkie 20 aminokwasów, test leave-protein-out.
- `scripts/error_scan.ts` (+ kopia w korzeniu, bo importuje ./src) — komparator błędu modelu vs struktura natywna.
- `results/20aa_sequence_only.txt` — UCZCIWA wersja (tylko sekwencja + zakopanie).
- `results/20aa_with_SS_circular.txt` — wersja z etykietami SS z geometrii CA: CZĘŚCIOWO KOŁOWA (SS i region phi/psi z tej samej geometrii), nie używać jako dowodu.
Uwaga: SS liczone z geometrii CA, nie DSSP; warunek "bur<99" = zakopanie ≥20 sąsiadów CA w 10 Å.
