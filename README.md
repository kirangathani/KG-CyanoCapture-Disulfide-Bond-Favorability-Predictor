# Protein Cysteine Analysis

This Python script analyzes protein structures to determine the protonation states of cysteine residues under different pH conditions. It uses PROPKA for pKa predictions and provides insights into disulfide bonds and cysteine behavior in various environments.

## Features

- Downloads PDB files for specified proteins
- Runs PROPKA to predict pKa values for cysteine residues
- Analyzes existing disulfide bonds
- Calculates protonation states for different pH environments (E. coli and Cyanobacteria)
- Generates a summary CSV file with analysis results

## Requirements

- Python 3.6+
- Packages: pandas, biopython, requests
- PROPKA3 (installed via pip)

## Setup

1. Ensure you have Python 3.6 or higher installed.
2. Clone this repository or download the script.
3. The script will create a virtual environment and install required packages automatically.

## Usage

1. Prepare a CSV file named `protein_data.csv` with columns:
   - protein_name
   - pdb_id
   - uniprot_id
   - sequence

2. Run the script:
python protein_analysis.py
Copy
3. The script will:
- Set up a virtual environment
- Install required packages
- Download PDB files
- Run PROPKA analysis
- Generate a summary CSV file

## Output

The script generates a CSV file named `cysteine_analysis_summary.csv` with the following columns:

- Protein Name
- PDB ID
- Existing Disulfides
- Total Cysteines
- Cysteines with pKa
- Cysteines in Disulfide Bonds
- Protonated Cysteines (E. coli)
- Deprotonated Cysteines (E. coli)
- Disulfide Cysteines (E. coli)
- Protonated Cysteines (Cyanobacteria)
- Deprotonated Cysteines (Cyanobacteria)
- Disulfide Cysteines (Cyanobacteria)

## Notes

- The script uses pH 7.2 for E. coli and pH 7.5 for Cyanobacteria environments.
- PROPKA analysis is run at neutral pH (7.0).
- Ensure you have a stable internet connection for downloading PDB files.

## Troubleshooting

If you encounter issues:
- Check your internet connection
- Ensure you have the correct permissions to create directories and files
- Verify that the `protein_data.csv` file is in the correct format and location

For any other issues, please open an issue in the repository.
