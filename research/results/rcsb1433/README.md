# Wyniki na 1433 niepowtarzalnych bialkach z RCSB (14.11.2026)
Zrodlo: RCSB, X-ray <=2.0 A, 50-300 aa, klastry 30% identycznosci (skrypt: research/scripts/download_pdb.ps1), po usunieciu >80% identycznych.
205 998 reszt. Test na bialkach NIEWIDZIANYCH (20% bialek, podzial po bialkach, nie po resztach).
- comparator_learning_curve.txt: krzywa uczenia (logloss 0,766 -> 0,749 -> 0,736 -> 0,731 dla 160/320/640/1147 bialek: splaszcza sie)
- 20aa_sequence_only_1433.txt: tabele dla 20 aminokwasow (tylko sekwencja + zakopanie, bez etykiet SS)
Sufit sekwencji: ok. 0,71-0,73 logloss; dodanie zakopania (info ze struktury) daje tylko -0,025.
