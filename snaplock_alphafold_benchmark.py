#!/usr/bin/env python3
"""
SNAPLOCK vs AlphaFold Benchmark & Calibration Engine v1.0
Deterministyczne mapowanie struktury 3D na sferyczną siatkę 248 wektorów E8.
Oficjalny moduł walidacyjny silnika SNAPLOCK.
"""

import json
import math
import sys
import os
import urllib.request
import urllib.error

# -----------------------------------------------------------------------------
# 1. GENERATOR DYSKRETNEJ SIATKI 248 WEKTORÓW (E8 SFERYCZNE MAPOWANIE)
# -----------------------------------------------------------------------------
def generate_e8_248_grid():
    """
    Generuje 248 znormalizowanych wektorów na sferze jednostkowej S^2
    z rygorystycznym uwzględnieniem bazowych osi SNAPLOCK:
    #1: [0, 1, 0] (Góra +Y - Wejście N)
    #248: [0, -1, 0] (Dół -Y - Wyjście C)
    #100: [1, 0, 0] (Prawo +X)
    #200: [-1, 0, 0] (Lewo -X)
    #50: [0, 0, 1] (Przód +Z)
    #150: [0, 0, -1] (Tył -Z)
    """
    vectors = {}
    
    # Punkty kardynalne
    vectors[1] = [0.0, 1.0, 0.0]
    vectors[248] = [0.0, -1.0, 0.0]
    vectors[100] = [1.0, 0.0, 0.0]
    vectors[200] = [-1.0, 0.0, 0.0]
    vectors[50] = [0.0, 0.0, 1.0]
    vectors[150] = [0.0, 0.0, -1.0]
    
    # Deterministyczny rozkład pozostałych 242 wektorów (spirala sferyczna Fibonacciego)
    assigned_indices = set(vectors.keys())
    remaining_indices = [i for i in range(1, 249) if i not in assigned_indices]
    
    n_points = 248
    golden_ratio = (1.0 + math.sqrt(5.0)) / 2.0
    
    for idx, patch_id in enumerate(remaining_indices):
        theta = 2.0 * math.pi * idx / golden_ratio
        phi = math.acos(1.0 - 2.0 * (idx + 0.5) / n_points)
        x = math.sin(phi) * math.cos(theta)
        y = math.cos(phi)
        z = math.sin(phi) * math.sin(theta)
        
        norm = math.sqrt(x*x + y*y + z*z)
        if norm > 0:
            vectors[patch_id] = [x / norm, y / norm, z / norm]
        else:
            vectors[patch_id] = [0.0, 1.0, 0.0]
        
    return vectors

# -----------------------------------------------------------------------------
# 2. POBIERANIE STRUKTURY Z DARMOWEJ BAZY ALPHAFOLD DB (EBI API)
# -----------------------------------------------------------------------------
def fetch_alphafold_pdb(uniprot_id: str) -> str:
    """
    Pobiera strukturę PDB bezpośrednio z AlphaFold Database (EMBL-EBI).
    """
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    req = urllib.request.Request(api_url, headers={'User-Agent': 'SNAPLOCK-Engine-Validator/1.0'})
    
    try:
        with urllib.request.urlopen(req, timeout=20) as resp:
            data = json.loads(resp.read().decode('utf-8'))
    except Exception as e:
        raise ValueError(f"Błąd pobierania metadanych AlphaFold dla UniProt ID {uniprot_id}: {e}")
    
    if not data or "pdbUrl" not in data[0]:
        raise ValueError(f"Brak pliku PDB dla ID: {uniprot_id}")
    
    pdb_url = data[0]["pdbUrl"]
    pdb_req = urllib.request.Request(pdb_url, headers={'User-Agent': 'SNAPLOCK-Engine-Validator/1.0'})
    
    try:
        with urllib.request.urlopen(pdb_req, timeout=30) as pdb_resp:
            pdb_content = pdb_resp.read().decode('utf-8')
    except Exception as e:
        raise ValueError(f"Nie udało się pobrać pliku PDB z adresu {pdb_url}: {e}")
        
    local_filename = f"{uniprot_id}_alphafold.pdb"
    with open(local_filename, "w", encoding="utf-8") as f:
        f.write(pdb_content)
        
    print(f"[*] Pomyślnie pobrano model AlphaFold: {local_filename}")
    return local_filename

# -----------------------------------------------------------------------------
# 3. SILNIK KWANTYZACJI I KALIBRACJI SNAPLOCK
# -----------------------------------------------------------------------------
class SnaplockCalibrator:
    def __init__(self):
        self.grid = generate_e8_248_grid()

    def find_nearest_patch(self, vector):
        """
        Znajduje indeks łatki o najwyższym iloczynie skalarnym (najmniejszy kąt).
        Patch_best = argmax (v_ref . e_i)
        """
        vx, vy, vz = vector
        norm = math.sqrt(vx*vx + vy*vy + vz*vz)
        if norm == 0:
            return 1, 0.0
        vx, vy, vz = vx / norm, vy / norm, vz / norm
        
        best_patch = None
        max_dot = -float('inf')
        
        for patch_id, (ex, ey, ez) in self.grid.items():
            dot = vx * ex + vy * ey + vz * ez
            if dot > max_dot:
                max_dot = dot
                best_patch = patch_id
                
        # Obliczenie błędu kątowego (w stopniach)
        clamped_dot = max(-1.0, min(1.0, max_dot))
        angle_error_deg = math.degrees(math.acos(clamped_dot))
        return best_patch, angle_error_deg

    def run_benchmark(self, pdb_file: str):
        # Deterministyczny parser PDB
        residues_dict = {}
        with open(pdb_file, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith("ATOM  "):
                    atom_name = line[12:16].strip()
                    res_name = line[17:20].strip()
                    chain_id = line[21].strip()
                    res_seq = int(line[22:26].strip())
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    plddt = float(line[60:66].strip()) if len(line) >= 66 else 0.0
                    
                    key = (chain_id, res_seq, res_name)
                    if key not in residues_dict:
                        residues_dict[key] = {"res_seq": res_seq, "res_name": res_name, "chain": chain_id, "atoms": {}, "plddt": plddt}
                    residues_dict[key]["atoms"][atom_name] = [x, y, z]

        residues = []
        for key in sorted(residues_dict.keys(), key=lambda k: k[1]):
            res = residues_dict[key]
            if "CA" in res["atoms"]:
                residues.append(res)
                        
        n_res = len(residues)
        print(f"[*] Liczba zanalizowanych reszt szkieletowych: {n_res}")
        
        # Test wektorów translacji szkieletu C_alfa_i -> C_alfa_{i+1}
        patch_assignments = []
        total_angular_error = 0.0
        
        for i in range(n_res - 1):
            ca_current = residues[i]["atoms"]["CA"]
            ca_next = residues[i+1]["atoms"]["CA"]
            bond_vector = [
                ca_next[0] - ca_current[0],
                ca_next[1] - ca_current[1],
                ca_next[2] - ca_current[2]
            ]
            
            patch_id, err_deg = self.find_nearest_patch(bond_vector)
            total_angular_error += err_deg
            
            patch_assignments.append({
                "res_idx": i + 1,
                "res_name": residues[i]["res_name"],
                "next_res": residues[i+1]["res_name"],
                "assigned_patch": patch_id,
                "angle_error_deg": round(err_deg, 3)
            })

        mean_angular_loss = total_angular_error / (n_res - 1) if n_res > 1 else 0.0
        
        # Weryfikacja mostków kowalencyjnych S-S (Ścieżka 5 SNAPLOCK)
        cys_residues = [r for r in residues if r["res_name"] == "CYS"]
        disulfide_bridges = []
        for i in range(len(cys_residues)):
            for j in range(i + 1, len(cys_residues)):
                if "SG" in cys_residues[i]["atoms"] and "SG" in cys_residues[j]["atoms"]:
                    sg1 = cys_residues[i]["atoms"]["SG"]
                    sg2 = cys_residues[j]["atoms"]["SG"]
                    dx = sg1[0] - sg2[0]
                    dy = sg1[1] - sg2[1]
                    dz = sg1[2] - sg2[2]
                    dist = math.sqrt(dx*dx + dy*dy + dz*dz)
                    if dist <= 2.5: # Standardowy mostek S-S ~2.05 Å
                        disulfide_bridges.append((cys_residues[i]["res_seq"], cys_residues[j]["res_seq"], round(dist, 2)))

        plddt_list = [r["plddt"] for r in residues if r.get("plddt")]
        plddt_mean = round(sum(plddt_list) / len(plddt_list), 2) if plddt_list else None

        report = {
            "total_residues": n_res,
            "mean_angular_patch_loss_deg": round(mean_angular_loss, 3),
            "disulfide_bridges_detected": disulfide_bridges,
            "mean_plddt_score": plddt_mean,
            "sample_mappings": patch_assignments[:5],
            "patch_assignments_count": len(patch_assignments)
        }
        
        return report

# -----------------------------------------------------------------------------
# 4. PUNKT STARTOWY
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    # Domyślny test: TC5b (Trp-cage, 20 reszt, UniProt: P83347) lub Insulina (P01308)
    uniprot_target = "P83347"
    if len(sys.argv) > 1:
        uniprot_target = sys.argv[1]
        
    print(f"=== ROZPOCZYNAM KALIBRACJĘ SNAPLOCK DLA {uniprot_target} ===")
    pdb_path = fetch_alphafold_pdb(uniprot_target)
    
    calibrator = SnaplockCalibrator()
    results = calibrator.run_benchmark(pdb_path)
    
    print("\n--- WYNIKI KALIBRACJI Z ALPHAFOLD DB ---")
    print(json.dumps(results, indent=2))
