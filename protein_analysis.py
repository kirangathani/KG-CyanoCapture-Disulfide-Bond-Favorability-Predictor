import os
import subprocess
import sys
import venv
import platform
import pandas as pd
from Bio.PDB import PDBParser, Selection
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
        subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
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
        return {}
    cys_pkas = {}
    with open(pka_file, 'r') as f:
        for line in f:
            if line.startswith('CYS'):
                parts = line.split()
                residue_number = int(parts[1])
                chain = parts[2]
                pka = float(parts[3])
                if pka != 99.99:
                    cys_pkas[(chain, residue_number)] = pka
    return cys_pkas

def get_existing_disulfides(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)
    existing_disulfides = set()
    ssbond_records = set()

    # First, read SSBOND records from the PDB file
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('SSBOND'):
                chain1 = line[15]
                resnum1 = int(line[17:21])
                chain2 = line[29]
                resnum2 = int(line[31:35])
                ssbond_records.add(tuple(sorted([(chain1, resnum1), (chain2, resnum2)])))

    for model in structure:
        for chain in model:
            cys_residues = [res for res in chain if res.resname == 'CYS' and 'SG' in res]
            
            for i, cys1 in enumerate(cys_residues):
                for cys2 in cys_residues[i+1:]:
                    distance = cys1['SG'] - cys2['SG']
                    if distance < 2.1:  # Slightly relaxed threshold
                        bond = tuple(sorted([(chain.id, cys1.id[1]), (chain.id, cys2.id[1])]))
                        existing_disulfides.add(bond)
                        print(f"Detected disulfide bond: {bond}, distance: {distance:.2f} Å")

    # Cross-reference with SSBOND records
    for bond in ssbond_records:
        if bond not in existing_disulfides:
            print(f"WARNING: SSBOND record {bond} not detected by distance method")
            existing_disulfides.add(bond)

    return existing_disulfides

def calculate_cys_distances(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)
    cys_atoms = {}
    total_cysteines = 0
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.resname == 'CYS':
                    total_cysteines += 1
                    if 'SG' in residue:
                        cys_atoms[(chain.id, residue.id[1])] = residue['SG']

    distances = {}
    for cys1, atom1 in cys_atoms.items():
        for cys2, atom2 in cys_atoms.items():
            if cys1 < cys2:
                distance = atom1 - atom2
                distances[(cys1, cys2)] = distance
    return distances, cys_atoms, total_cysteines

def create_reduced_pdb(pdb_file):
    reduced_pdb = f"{os.path.splitext(pdb_file)[0]}_reduced.pdb"
    with open(pdb_file, 'r') as original, open(reduced_pdb, 'w') as reduced:
        for line in original:
            if not line.startswith('SSBOND'):
                reduced.write(line)
    return reduced_pdb
