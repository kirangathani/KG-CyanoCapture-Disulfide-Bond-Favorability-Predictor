# CyanoCapture-Disulfide-Bond-Favorability-Predictor

# Protein Disulfide Bond Analysis Script

This script analyzes protein structures to predict the likelihood of disulfide bond formation in different cellular environments (E. coli and Cyanobacteria). It uses PROPKA3 for pKa predictions and considers factors such as cysteine distances, existing disulfide bonds, and environmental pH.

## Prerequisites

- Python 3.6 or higher
- PROPKA3 installed and accessible from the command line
- Required Python packages: pandas, biopython, requests

You can install the required packages using:
pip install pandas biopython requests

## Input

The script requires a CSV file named "protein_data.csv" in the same directory, containing the following columns:
- protein_name
- pdb_id
- uniprot_id
- sequence

## Usage

1. Ensure "protein_data.csv" is in the same directory as the script.
2. Run the script using: `python protein_analysis.py`

## Output

The script provides:
1. Detailed analysis for each protein, including:
   - Existing disulfide bonds
   - Cysteine pKa values
   - Distances between cysteine pairs
   - Potential disulfide bonds in E. coli and Cyanobacteria environments
   - Favorability scores for disulfide bond formation in each environment
2. A summary table (printed to console and saved as "disulfide_analysis_summary.csv") containing:
   - Protein Name
   - PDB ID
   - Number of existing disulfide bonds
   - Total number of cysteines
   - Number of potential disulfide bonds in E. coli and Cyanobacteria
   - Favorability scores for E. coli and Cyanobacteria
   - Preferred environment for disulfide bond formation

## Interpreting the Favorability Score

The favorability score is a measure of how likely a disulfide bond is to form, based on two main factors:
1. The reactivity of the cysteine residues
2. The distance between the cysteine residues

The score ranges from 0 to 1, where:
- 0 indicates very unfavorable conditions for disulfide bond formation
- 1 indicates highly favorable conditions for disulfide bond formation

The score is calculated as follows:
1. Reactivity: Based on the pKa of the cysteine and the pH of the environment. Higher reactivity occurs when the pH is above the pKa.
2. Distance factor: Based on the distance between cysteine residues. Shorter distances are more favorable.
3. The final score is the average of the reactivity and distance factors.

## How Disulfide Bonds Influence the Score

Existing disulfide bonds in a protein structure are taken into account in the following ways:
1. They are excluded from the potential new disulfide bond calculations.
2. They may influence the pKa values of nearby cysteines, potentially affecting their reactivity.
3. They contribute to the overall count of disulfide bonds in the protein, which is reported in the summary.

The presence of existing disulfide bonds can indirectly affect the favorability scores of potential new bonds by altering the protein's structure and the local environment of other cysteines.

## Limitations

- The analysis is based on static protein structures and does not account for protein dynamics.
- It does not consider other factors that might affect protein folding and stability in different cellular environments.
- The prediction model is simplified and may not capture all complexities of cellular environments.

## Troubleshooting

If you encounter issues:
1. Ensure PROPKA3 is correctly installed and accessible from the command line.
2. Check that all required Python packages are installed.
3. Verify that "protein_data.csv" is correctly formatted and in the same directory as the script.
4. Check console output for any error messages during execution.
