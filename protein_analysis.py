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
    
    # Create virtual environment if it doesn't exist
    if not os.path.exists(venv_dir):
        print("Creating virtual environment...")
        venv.create(venv_dir, with_pip=True)
    
    # Determine the correct activation script based on the OS
    if platform.system() == "Windows":
        activate_this = os.path.join(venv_dir, 'Scripts', 'activate_this.py')
    else:  # macOS and Linux
        activate_this = os.path.join(venv_dir, 'bin', 'activate_this.py')
    
    # Activate the virtual environment
    if os.path.exists(activate_this):
        with open(activate_this) as file:
            exec(file.read(), {'__file__': activate_this})
    else:
        print(f"Warning: Could not find {activate_this}")
        print("Attempting to continue without activating the virtual environment...")
    
    # Install required packages
    print("Installing required packages...")
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "pandas", "biopython", "requests"])
    except subprocess.CalledProcessError as e:
        print(f"Error installing packages: {e}")
        sys.exit(1)
    
    # Install PROPKA3 if not already installed
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
    else:
        print(f"Failed to download {pdb_id}.pdb")
        return False
    return True

def run_propka(pdb_file, ph_value):
    command = f"propka3 {pdb_file} --pH {ph_value}"
    try:
        subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running PROPKA: {e}")
        print(f"PROPKA stdout: {e.stdout}")
        print(f"PROPKA stderr: {e.stderr}")
        return False
    return True

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
                if pka != 99.99:  # Only store valid pKa values
                    cys_pkas[(chain, residue_number)] = pka
    return cys_pkas

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

def estimate_disulfide_potential(cys_pkas, distances, existing_disulfides, ph):
    potential_disulfides = []
    for (cys1, cys2), distance in distances.items():
        if (cys1, cys2) not in existing_disulfides and (cys2, cys1) not in existing_disulfides:
            if distance <= 8.0:
                pka1 = cys_pkas.get(cys1, 8.3)  # Use 8.3 as default pKa
                pka2 = cys_pkas.get(cys2, 8.3)
                
                reactivity1 = 1 / (1 + 10**(pka1 - ph))
                reactivity2 = 1 / (1 + 10**(pka2 - ph))
                avg_reactivity = (reactivity1 + reactivity2) / 2
                
                # Modified distance factor calculation
                distance_factor = 1 - (distance / 8.0)
                
                # Calculate favorability based on reactivity and distance
                favorability = (avg_reactivity + distance_factor) / 2
                
                potential_disulfides.append((cys1, cys2, distance, avg_reactivity, favorability))
                
                # Print intermediate values for debugging
                print(f"  Pair: {cys1} - {cys2}")
                print(f"  pKa values: {pka1:.2f}, {pka2:.2f}")
                print(f"  Reactivities: {reactivity1:.4f}, {reactivity2:.4f}")
                print(f"  Distance: {distance:.2f}, Distance factor: {distance_factor:.4f}")
                print(f"  Favorability: {favorability:.4f}")
                print()
                
    return sorted(potential_disulfides, key=lambda x: x[4], reverse=True)

def calculate_environment_favorability(results):
    return sum(bond[4] for bond in results)

def analyze_protein(pdb_id, pdb_file, ph_values):
    if not os.path.exists(pdb_file):
        if not download_pdb(pdb_id):
            return None, None, None, None, None, None

    if not run_propka(pdb_file, 7.0):  # Run PROPKA at neutral pH
        return None, None, None, None, None, None

    cys_pkas = parse_propka_output(pdb_id)
    existing_disulfides = get_existing_disulfides(pdb_file)
    distances, cys_atoms, total_cysteines = calculate_cys_distances(pdb_file)

    results = {}
    favorability_scores = {}
    for ph in ph_values:
        potential_disulfides = estimate_disulfide_potential(cys_pkas, distances, existing_disulfides, ph)
        results[ph] = potential_disulfides
        favorability_scores[ph] = calculate_environment_favorability(potential_disulfides)

    return results, existing_disulfides, cys_pkas, distances, favorability_scores, total_cysteines

def main():
    # Set up the environment
    setup_environment()
    
    df = pd.read_csv("protein_data.csv")
    ph_values = {'E. coli': 7.2, 'Cyanobacteria': 7.5}
    
    summary_data = []

    for _, row in df.iterrows():
        protein_name = row['protein_name']
        pdb_id = row['pdb_id']
        pdb_file = f"{pdb_id}.pdb"
        
        print(f"\nAnalyzing {protein_name} (PDB: {pdb_id}):")
        
        if not os.path.exists(pdb_file):
            if not download_pdb(pdb_id):
                print(f"Skipping analysis for {protein_name} due to download failure.")
                continue
        
        results, existing_disulfides, cys_pkas, distances, favorability_scores, total_cysteines = analyze_protein(pdb_id, pdb_file, ph_values.values())
        
        if results is None:
            print(f"Analysis failed for {protein_name} (PDB: {pdb_id})")
            continue

        print("Existing disulfide bonds:")
        for bond in existing_disulfides:
            print(f"  {bond[0]} - {bond[1]}")
        
        print("\nCysteine pKa values:")
        for cys, pka in cys_pkas.items():
            print(f"  {cys}: {pka:.2f}")
        
        print("\nDistances between cysteine pairs:")
        for (cys1, cys2), distance in distances.items():
            print(f"  {cys1} - {cys2}: {distance:.2f} Å")
        
        potential_disulfides = {}
        for env, ph in ph_values.items():
            print(f"\nPotential disulfide bonds in {env} (pH {ph}):")
            potential_disulfides[env] = len(results[ph])
            for cys1, cys2, distance, reactivity, favorability in results[ph]:
                print(f"  {cys1} - {cys2}: Distance = {distance:.2f} Å, Reactivity = {reactivity:.2f}, Favorability = {favorability:.4f}")
            print(f"Overall favorability score for {env}: {favorability_scores[ph]:.4f}")
        
        if favorability_scores[ph_values['E. coli']] > favorability_scores[ph_values['Cyanobacteria']]:
            preferred_environment = 'E. coli'
        elif favorability_scores[ph_values['E. coli']] < favorability_scores[ph_values['Cyanobacteria']]:
            preferred_environment = 'Cyanobacteria'
        else:
            preferred_environment = 'Equal'
        
        summary_data.append({
            'Protein Name': protein_name,
            'PDB ID': pdb_id,
            'Existing Disulfides': len(existing_disulfides),
            'Cysteines': total_cysteines,
            'Potential Disulfides (E. coli)': potential_disulfides['E. coli'],
            'Potential Disulfides (Cyanobacteria)': potential_disulfides['Cyanobacteria'],
            'Favorability Score (E. coli)': favorability_scores[ph_values['E. coli']],
            'Favorability Score (Cyanobacteria)': favorability_scores[ph_values['Cyanobacteria']],
            'Preferred Environment': preferred_environment
        })

    # Create summary table
    summary_df = pd.DataFrame(summary_data)
    print("\nSummary Table:")
    print(summary_df.to_string(index=False))
    
    # Save summary table to CSV
    summary_df.to_csv("disulfide_analysis_summary.csv", index=False)
    print("\nSummary table saved to disulfide_analysis_summary.csv")

if __name__ == "__main__":
    main()