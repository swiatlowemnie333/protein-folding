"""OPIEKUN REGUŁ (rules keeper): audyt stałych w 'magazynie' SNAPLOCK.
Warstwa 1: spójność wewnętrzna tabeli 22 aminokwasów (pI z pK, ładunek z Henderson-Hasselbalch, jonizacja, null-e).
Warstwa 2: kolumna 'preferencja Ramachandrana' vs dane (1433 białka), bazeny kątów silnika vs piki danych.
Uruchom: python3 rules_keeper.py db22.json rows.pkl
"""
import json, pickle, sys, math
import numpy as np
db = {r['code1']: r for r in json.load(open(sys.argv[1]))}
rows = pickle.load(open(sys.argv[2], 'rb'))
out = []; bad = 0
def flag(level, msg):
    global bad
    if level != 'OK': bad += 1
    out.append(f'[{level}] {msg}')

print('=== WARSTWA 1: spójność wewnętrzna ===')
for c, r in db.items():
    pk1, pk2, pkr, pi, q = r['pK1_COOH'], r['pK2_NH3'], r['pKr_side'], r['pI'], r['chargePH74']
    # pI
    if pkr is None: pi_calc = (pk1+pk2)/2
    elif r['chargePH74'] < 0 or c in 'DEC': pi_calc = (pk1+pkr)/2
    else: pi_calc = (pkr+pk2)/2
    if c == 'Y': pi_calc = (pk1+pk2)/2  # konwencja: pKr(OH)=10.07 > pK2
    if abs(pi_calc-pi) > 0.1: flag('PODEJRZANE', f'{c} pI={pi} vs z pK={pi_calc:.2f}')
    # ładunek
    if pkr is not None:
        if c in 'DEUC' or (c=='Y'):  q_calc = -1/(1+10**(pkr-7.4))
        else: q_calc = 1/(1+10**(7.4-pkr))
        if abs(q_calc-q) > 0.15 and c != 'Y': flag('PODEJRZANE', f'{c} ładunek pH7.4={q} vs Henderson-Hasselbalch={q_calc:.2f} (pKr={pkr})')
    elif q != 0: flag('PODEJRZANE', f'{c} ładunek {q} bez pKr')
    if r['sidechainIonization'] != (pkr is not None and abs(q) > 0.3) and c not in 'CYH':
        flag('UWAGA', f'{c} flaga jonizacji={r["sidechainIonization"]} vs pKr/ładunek')
# monotoniczność wielkości: objętość vdW Cys<Sec, Ser<Cys itd.
if db['U']['vdwVolume'] <= db['C']['vdwVolume']: flag('PODEJRZANE', f"SEC vdW {db['U']['vdwVolume']} <= CYS {db['C']['vdwVolume']} (Se > S, powinno być większe)")
if db['U']['mass'] < 168 or db['U']['mass'] > 169: flag('PODEJRZANE', 'masa SEC')

print('=== WARSTWA 2: kolumna Ramachandran vs dane ===')
def region(phi, psi):
    if phi > 0: return 'L'
    if -160 <= phi <= -20 and -120 <= psi <= 50: return 'A'
    if phi < -45 and (psi >= 90 or psi <= -150): return 'B'
    return 'O'
cnt = {}
for x in rows:
    cnt.setdefault(x['aa'], {'A':0,'B':0,'L':0,'O':0})[region(x['phi'], x['psi'])] += 1
txt_class = {  # co mówi tabela -> oczekiwana dominująca klasa
 'Helisa alfa':'A','Kartka beta':'B','Alfa / Beta':'AB','Złamanie/Kink':'P','Wiele konformacji':'G',
 'Zwroty / Helisa':'AO','Zwroty / Pętle':'O','N-cap helisy':'A','Zmienna':'?','Amfipatyczna':'(nie klasa Rama)',
 'Centra aktywne':'(nie klasa Rama)','Helisa / Mostki':'A','Helisa / Pętla':'A'}
print(f"{'aa':3} {'A%':>5} {'B%':>5} {'L%':>5} {'O%':>5}  tabela-mówi")
for c in 'GAVLIPFYWSTCMNQDEKRH':
    n = sum(cnt[c].values()); f = {k: 100*v/n for k, v in cnt[c].items()}
    t = db[c]['ramachandranPref']
    note = ''
    exp = txt_class.get(t, '?')
    if exp.startswith('('): note = '<- opis nie jest klasą Ramachandrana (do przepisania na dane)'
    elif exp == 'A' and f['A'] < f['B']: note = '<- DANE: więcej beta niż alfa'
    elif exp == 'B' and f['B'] < f['A']: note = '<- DANE: więcej alfa niż beta'
    print(f"{c:3} {f['A']:5.1f} {f['B']:5.1f} {f['L']:5.1f} {f['O']:5.1f}  {t:18} {note}")
    if note: flag('PODEJRZANE', f'{c} "{t}" {note}')

print('=== BAZENY KĄTÓW SILNIKA vs piki danych (5° siatka) ===')
def peak(sel, lo_phi, hi_phi, lo_psi, hi_psi, bin=5):
    pts = np.array([(x['phi'], x['psi']) for x in sel if lo_phi<=x['phi']<=hi_phi and lo_psi<=x['psi']<=hi_psi])
    H, xe, ye = np.histogram2d(pts[:,0], pts[:,1], bins=[np.arange(lo_phi,hi_phi+bin,bin), np.arange(lo_psi,hi_psi+bin,bin)])
    i, j = np.unravel_index(H.argmax(), H.shape)
    return xe[i]+bin/2, ye[j]+bin/2, len(pts)
gly = [x for x in rows if x['aa']=='G']; non = [x for x in rows if x['aa'] not in 'GP']; pro = [x for x in rows if x['aa']=='P']
tests = [('alfa (nie-Gly) -63,-42', non, -120,-20,-80,10, (-63,-42)),
         ('beta (nie-Gly) -120,130', non, -180,-60,90,180, (-120,130)),
         ('PPII (nie-Gly) -75,145', non, -100,-50,110,180, (-75,145)),
         ('baz (-80,0) nie-Gly', non, -120,-50,-30,30, (-80,0)),
         ('Gly lewy baz 60,45', gly, 30,110,-30,90, (60,45)),
         ('Gly prawy alfa', gly, -120,-40,-70,20, (-63,-42)),
         ('Gly beta', gly, -180,-60,90,180, (-120,130)),
         ('Pro (-65,145)', pro, -100,-40,100,180, (-75,145))]
for name, sel, a,b,c_,d, eng in tests:
    px, py, n = peak(sel, a,b,c_,d)
    dist = math.hypot(px-eng[0], py-eng[1])
    lvl = 'OK' if dist < 15 else 'PODEJRZANE'
    flag(lvl, f'bazen {name}: silnik {eng}, pik danych ({px:.0f},{py:.0f}) n={n}, odległość {dist:.0f}°')
    print(out[-1])
print(f'\nRAZEM: {bad} pozycji do sprawdzenia')
open('/home/claude/zamys-/research/results/opiekun/audyt_regul.txt','w').write('\n'.join(out)+f'\nRAZEM {bad}\n')
