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

- PROPKA analysis is run at a stable, neutral pH (7.0) for all proteins.
- The script uses the pKa values from this neutral pH analysis to predict protonation states at pH 7.2 (E. coli environment) and pH 7.5 (Cyanobacteria environment).
- This approach assumes that pKa values don't significantly change over this small pH range.
- PropKa is a tool for predicting pKa values of titratable groups in proteins, including cysteines. This README explains how PropKa calculates pKa values for cysteine residues.
Calculation Method
- PropKa uses a semi-empirical approach that considers various structural and environmental factors:
1. Base pKa: Starts with a reference pKa of ~9.0 for cysteines in solution.
2. Desolvation Effects: Calculates burial depth within the protein structure.
3. Hydrogen Bonding: Identifies potential H-bonds involving the thiol group.
4. Electrostatic Interactions: Considers nearby charged groups.
5. Coulombic Interactions: Calculates longer-range electrostatic effects.
6. Empirical Rules: Applies rules for specific structural contexts.
7. Iterative Calculation: Converges on consistent pKa predictions for all titratable groups.
8. Temperature Dependence: Adjusts predictions based on temperature.

## Troubleshooting

If you encounter issues:
- Check your internet connection
- Ensure you have the correct permissions to create directories and files
- Verify that the `protein_data.csv` file is in the correct format and location

For any other issues, please open an issue in the repository.
