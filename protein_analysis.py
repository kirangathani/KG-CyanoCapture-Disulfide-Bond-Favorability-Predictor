import os
import subprocess
import sys
import venv
import platform
import pandas as pd
from Bio.PDB import PDBParser
import requests

def setup_environment():
    venv_dir = 'protein_analysis_env'
    
    if not os.path.exists(venv_dir):
        print("Creating virtual environment...")
        venv.create(venv_dir, with_pip=True)
    
    if platform.system() == "Windows":
        activate_this = os.path.join(venv_dir, 'Scripts', 'activate_this.py')
    else:
        activate_this = os.path.join(venv_dir, 'bin', 'activate_this.py')
    
    if os.path.exists(activate_this):
        with open(activate_this) as file:
            exec(file.read(), {'__file__': activate_this})
    else:
        print(f"Warning: Could not find {activate_this}")
        print("Attempting to continue without activating the virtual environment...")
    
    print("Installing required packages...")
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "pandas", "biopython", "requests"])
    except subprocess.CalledProcessError as e:
        print(f"Error installing packages: {e}")
        sys.exit(1)
    
    try:
        subprocess.check_call(["propka3", "--version"])
    except FileNotFoundError:
        print("Installing PROPKA3...")
        try:
            subprocess.check_call([sys.executable, "-m", "pip", "install", "propka"])
        except subprocess.CalledProcessError as e:
            print(f"Error installing PROPKA3: {e}")
            sys.exit(1)

def download_pdb(pdb_id):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        with open(f"{pdb_id}.pdb", 'wb') as f:
            f.write(response.content)
        print(f"Downloaded {pdb_id}.pdb")
        return True
    else:
        print(f"Failed to download {pdb_id}.pdb")
        return False

def run_propka(pdb_file, ph_value):
    command = f"propka3 {pdb_file} --pH {ph_value}"
    try:
        result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
        print(f"PROPKA stdout: {result.stdout}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running PROPKA: {e}")
        print(f"PROPKA stdout: {e.stdout}")
        print(f"PROPKA stderr: {e.stderr}")
        return False

def parse_propka_output(pdb_id):
    pka_file = f"{pdb_id}.pka"
    if not os.path.exists(pka_file):
        print(f"PROPKA output file {pka_file} not found")
        return {}, set()
    cys_pkas = {}
    disulfide_cys = set()
    print(f"Parsing PROPKA output for {pdb_id}")
    with open(pka_file, 'r') as f:
        for line in f:
            if line.startswith('CYS'):
                parts = line.split()
                print(f"Found CYS line: {line.strip()}")
                residue_number = int(parts[1])
                chain = parts[2]
                pka = float(parts[3])
                if pka == 99.99:
                    disulfide_cys.add((chain, residue_number))
                    print(f"CYS {chain}:{residue_number} is in a disulfide bond")
                else:
                    cys_pkas[(chain, residue_number)] = pka
                    print(f"Stored pKa {pka} for CYS {chain}:{residue_number}")
    if not cys_pkas and not disulfide_cys:
        print(f"No valid CYS pKa values found for {pdb_id}")
    return cys_pkas, disulfide_cys

def get_existing_disulfides(pdb_file):
    existing_disulfides = set()
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('SSBOND'):
                parts = line.split()
                cys1 = (parts[3], int(parts[4]))
                cys2 = (parts[6], int(parts[7]))
                existing_disulfides.add((cys1, cys2))
    return existing_disulfides

def calculate_total_cysteines(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)
    total_cysteines = 0
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.resname == 'CYS':
                    total_cysteines += 1
    return total_cysteines

def predict_protonation_state(pka, ph):
    return "Deprotonated" if pka < ph else "Protonated"

def count_protonation_states(cys_pkas, ph):
    protonated = sum(1 for pka in cys_pkas.values() if pka >= ph)
    deprotonated = len(cys_pkas) - protonated
    print(f"pH: {ph}, Protonated: {protonated}, Deprotonated: {deprotonated}")
    return protonated, deprotonated

def analyze_protein(pdb_id, pdb_file, ph_values):
    print(f"\nAnalyzing {pdb_id}")
    
    if not os.path.exists(pdb_file):
        if not download_pdb(pdb_id):
            print(f"Failed to download PDB file for {pdb_id}")
            return None, None, None, None, None, None, None

    if not run_propka(pdb_file, 7.0):  # Run PROPKA at neutral pH
        print(f"PROPKA failed for {pdb_id}")
        return None, None, None, None, None, None, None

    cys_pkas, disulfide_cys = parse_propka_output(pdb_id)
    print(f"CYS pKas for {pdb_id}: {cys_pkas}")
    print(f"Disulfide cysteines for {pdb_id}: {disulfide_cys}")

    existing_disulfides = get_existing_disulfides(pdb_file)
    print(f"Existing disulfide bonds for {pdb_id}: {existing_disulfides}")

    total_cysteines = calculate_total_cysteines(pdb_file)
    print(f"Total cysteines in {pdb_id}: {total_cysteines}")

    protonation_counts = {}
    for env, ph in ph_values.items():
        protonated, deprotonated = count_protonation_states(cys_pkas, ph)
        protonation_counts[env] = {
            'Protonated': protonated, 
            'Deprotonated': deprotonated,
            'In Disulfide': len(disulfide_cys)
        }

    return existing_disulfides, total_cysteines, protonation_counts, cys_pkas, disulfide_cys

def main():
    setup_environment()
    
    df = pd.read_csv("protein_data.csv")
    ph_values = {'E. coli': 7.2, 'Cyanobacteria': 7.5}
    
    summary_data = []

    for _, row in df.iterrows():
        protein_name = row['protein_name']
        pdb_id = row['pdb_id']
        pdb_file = f"{pdb_id}.pdb"
        
        print(f"\nAnalyzing {protein_name} (PDB: {pdb_id}):")
        
        results = analyze_protein(pdb_id, pdb_file, ph_values)
        
        if results is None:
            print(f"Analysis failed for {protein_name} (PDB: {pdb_id})")
            continue

        existing_disulfides, total_cysteines, protonation_counts, cys_pkas, disulfide_cys = results

        protein_data = {
            'Protein Name': protein_name,
            'PDB ID': pdb_id,
            'Existing Disulfides': len(existing_disulfides),
            'Total Cysteines': total_cysteines,
            'Cysteines with pKa': len(cys_pkas),
            'Cysteines in Disulfide Bonds': len(disulfide_cys)
        }

        for env in ph_values.keys():
            protein_data.update({
                f'Protonated Cysteines ({env})': protonation_counts[env]['Protonated'],
                f'Deprotonated Cysteines ({env})': protonation_counts[env]['Deprotonated'],
                f'Disulfide Cysteines ({env})': protonation_counts[env]['In Disulfide']
            })

        summary_data.append(protein_data)

    summary_df = pd.DataFrame(summary_data)
    print("\nSummary Table:")
    print(summary_df.to_string(index=False))
    
    summary_df.to_csv("cysteine_analysis_summary.csv", index=False)
    print("\nSummary table saved to cysteine_analysis_summary.csv")

if __name__ == "__main__":
    main()