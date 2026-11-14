// Wygenerowane przez research/scripts (opiekun reguł) z 1433 białek RCSB (~206 tys. reszt). NIE edytować ręcznie.
// Frakcje (%) regionów Ramachandrana: A=alfa, B=beta, L=lewoskrętna (phi>0), O=inne. 'dominuje' = z danych, nie z opisu tekstowego.
export const RAMACHANDRAN_DATA_1433: Record<string, { A: number; B: number; L: number; O: number; n: number; dominuje: string }> = {
 "L": {
  "A": 56.4,
  "B": 41.8,
  "L": 0.7,
  "O": 1.1,
  "n": 18912,
  "dominuje": "alfa"
 },
 "F": {
  "A": 44.3,
  "B": 52.2,
  "L": 1.6,
  "O": 1.9,
  "n": 8494,
  "dominuje": "mieszana alfa/beta"
 },
 "I": {
  "A": 42.7,
  "B": 56.5,
  "L": 0.1,
  "O": 0.7,
  "n": 11934,
  "dominuje": "beta"
 },
 "A": {
  "A": 62.2,
  "B": 34.9,
  "L": 1.4,
  "O": 1.5,
  "n": 16252,
  "dominuje": "alfa"
 },
 "Y": {
  "A": 42.1,
  "B": 54.1,
  "L": 1.8,
  "O": 2.1,
  "n": 7369,
  "dominuje": "beta"
 },
 "D": {
  "A": 52.8,
  "B": 37.6,
  "L": 5.7,
  "O": 3.9,
  "n": 12084,
  "dominuje": "alfa"
 },
 "E": {
  "A": 63.5,
  "B": 33.2,
  "L": 2.2,
  "O": 1.1,
  "n": 14463,
  "dominuje": "alfa"
 },
 "R": {
  "A": 55.4,
  "B": 40.2,
  "L": 2.9,
  "O": 1.5,
  "n": 10687,
  "dominuje": "alfa"
 },
 "T": {
  "A": 43.7,
  "B": 55.1,
  "L": 0.4,
  "O": 0.7,
  "n": 11387,
  "dominuje": "beta"
 },
 "K": {
  "A": 56.7,
  "B": 38.4,
  "L": 3.4,
  "O": 1.4,
  "n": 12211,
  "dominuje": "alfa"
 },
 "G": {
  "A": 20.5,
  "B": 21.3,
  "L": 56.5,
  "O": 1.6,
  "n": 14724,
  "dominuje": "L (lewoskretna)"
 },
 "H": {
  "A": 47.1,
  "B": 45.6,
  "L": 4.2,
  "O": 3.2,
  "n": 4914,
  "dominuje": "mieszana alfa/beta"
 },
 "N": {
  "A": 46.5,
  "B": 36.0,
  "L": 12.8,
  "O": 4.8,
  "n": 8770,
  "dominuje": "alfa"
 },
 "W": {
  "A": 46.1,
  "B": 50.5,
  "L": 1.4,
  "O": 2.0,
  "n": 2923,
  "dominuje": "mieszana alfa/beta"
 },
 "S": {
  "A": 48.0,
  "B": 48.3,
  "L": 2.3,
  "O": 1.3,
  "n": 12234,
  "dominuje": "mieszana alfa/beta"
 },
 "C": {
  "A": 40.7,
  "B": 54.3,
  "L": 2.3,
  "O": 2.7,
  "n": 3152,
  "dominuje": "beta"
 },
 "P": {
  "A": 40.7,
  "B": 57.5,
  "L": 0.0,
  "O": 1.8,
  "n": 9265,
  "dominuje": "beta"
 },
 "V": {
  "A": 36.3,
  "B": 62.9,
  "L": 0.2,
  "O": 0.6,
  "n": 15234,
  "dominuje": "beta"
 },
 "Q": {
  "A": 60.3,
  "B": 35.4,
  "L": 2.7,
  "O": 1.6,
  "n": 7982,
  "dominuje": "alfa"
 },
 "M": {
  "A": 53.7,
  "B": 43.4,
  "L": 1.2,
  "O": 1.7,
  "n": 3007,
  "dominuje": "alfa"
 }
};
// Piki z danych (siatka 5°), do porównań z bazenami silnika:
export const RAMACHANDRAN_PEAKS_1433 = { alfa: [-62,-42], beta: [-118,128], ppii: [-62,138], gly_lewy: [82,2], gly_beta: [-78,172], nieGly_zwrot: [-62,-28] } as const;
