# PeptideScan

A Python package for peptide sequence analysis and biochemical characterization.

PeptideScan is an open-source Python package designed to analyze peptide sequences provided in FASTA files. It provides tools for peptide characterization, amino acid composition analysis, biochemical group classification, customizable filtering, and data export.

The package is intended to support researchers working with peptides, food proteins, proteomics, and bioactive peptide research.

## Features

PeptideScan currently provides the following functionalities:

**FASTA input**: read peptide sequences from FASTA files while preserving sequence identifiers and descriptions.

**Peptide properties**: calculate peptide length and monoisotopic molecular mass.

**Amino acid composition**: determine the count and relative frequency of each of the 20 standard amino acids.

**Biochemical group analysis**: calculate the composition and frequency of amino acid groups, including hydrophobic, polar uncharged, acidic, basic, and aromatic residues.

**Customizable criteria**: identify peptides enriched in selected amino acid groups using user-defined frequency thresholds.

**Peptide filtering**: select sequences according to length, molecular mass, amino acid occurrence, and biochemical group criteria.

**Data export**: save analysis results in CSV format or export selected peptide sequences to FASTA.

## Installation

Clone the repository:

[git clone] (https://github.com/jana-melo/PeptideScan.git)

Navigate to the project directory:

cd PeptideScan

### Install the package:

python -m pip install .

To install the development dependencies, including testing tools:

python -m pip install -e ".[dev]"

Quick start

The following example demonstrates how to analyze peptide sequences from a FASTA file and export the results to CSV.

from peptidescan import analyze_and_export

results = analyze_and_export(
    file_path="peptides.fasta",
    output_path="peptide_analysis.csv"
)

print(f"Analyzed peptides: {len(results['records'])}")
print(f"Results saved to: {results['csv_path']}")

## Available functions

The main functions available in PeptideScan include:

Module

Main functions

fasta.py

read_fasta()

properties.py

calculate_length(), calculate_molecular_mass()

composition.py

amino_acid_count(), amino_acid_composition(), contains_amino_acids()

amino_acid_groups.py

amino_acid_group_composition(), identify_rich_groups(), evaluate_group_criteria()

filtering.py

filter_peptides()

analysis.py

analyze_fasta()

export.py

export_to_csv(), export_to_fasta()

workflow.py

analyze_and_export()

Input format

PeptideScan accepts FASTA files containing peptide sequences.

### Example:

>P001 Peptide example 1
AGKPKY
>P002 Peptide example 2
FLWPK
>P003 Peptide example 3
VILMGA

The current implementation validates sequences against the 20 standard amino acid residue codes.

## Scientific considerations

Molecular mass calculations refer to monoisotopic masses of linear, unmodified peptide sequences.

Amino acid biochemical groups are operational classifications. Some groups overlap; for example, aromatic residues can also be classified as hydrophobic.

Frequency thresholds are user-defined screening criteria and should not be interpreted as universal thresholds for biological activity.

Amino acid composition and biochemical properties can help prioritize peptide candidates, but do not independently demonstrate biological activity.

## Testing

PeptideScan uses pytest for automated testing.

Run the test suite from the project root directory:

python -m pytest -v

### Project structure

```
PeptideScan/
│
├── pyproject.toml
├── README.md
├── LICENSE
│
├── src/
│   └── peptidescan/
│       ├── __init__.py
│       ├── fasta.py
│       ├── properties.py
│       ├── composition.py
│       ├── amino_acid_groups.py
│       ├── filtering.py
│       ├── export.py
│       ├── analysis.py
│       └── workflow.py
│
├── tests/
│
└── examples/
``` 


## Dependencies

PeptideScan uses the following Python libraries:

Biopython — FASTA sequence handling.

Pyteomics — molecular mass calculations.

pytest — automated testing during development.

## Contributing

Contributions, suggestions, and bug reports are welcome.

If you would like to contribute, please open an issue or submit a pull request through GitHub.

## License

This project is distributed under the MIT License. See the LICENSE file for details.

## Authors

Janaina Melo and Renata Menezes

GitHub: https://github.com/jana-melo

PeptideScan — peptide sequence analysis and biochemical characterization.