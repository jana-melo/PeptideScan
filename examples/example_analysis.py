from pathlib import Path
import pandas as pd

from peptidescan.fasta import read_fasta

from peptidescan.properties import (
    calculate_length,
    calculate_molecular_mass
)

from peptidescan.amino_acid_groups import (
    amino_acid_group_composition,
    identify_rich_groups
)


# --------------------------------------------------
# 1. Locate FASTA file
# --------------------------------------------------

project_root = Path(__file__).resolve().parent.parent

fasta_file = project_root / "examples" / "examples_proteins.fasta"


# --------------------------------------------------
# 2. Read FASTA sequences
# --------------------------------------------------

records = read_fasta(fasta_file)

print(f"\nSequences found: {len(records)}")


# --------------------------------------------------
# 3. Analyze peptide sequences
# --------------------------------------------------

results = []

for record in records:

    sequence = record["sequence"]

    # Basic properties
    length = calculate_length(sequence)

    molecular_mass = calculate_molecular_mass(sequence)

    # Amino acid group composition
    group_composition = amino_acid_group_composition(sequence)

    # Identify enriched groups (threshold = 30%)
    rich_groups = identify_rich_groups(
        sequence,
        threshold=30
    )

    # Organize results
    results.append({

        "ID": record["id"],

        "Sequence": sequence,

        "Length": length,

        "Molecular_Mass_Da": molecular_mass,

        "Hydrophobic (%)":
            group_composition["hydrophobic"]["frequency"],

        "Polar_uncharged (%)":
            group_composition["polar_uncharged"]["frequency"],

        "Acidic (%)":
            group_composition["acidic"]["frequency"],

        "Basic (%)":
            group_composition["basic"]["frequency"],

        "Aromatic (%)":
            group_composition["aromatic"]["frequency"],

        "Rich_groups": ", ".join(rich_groups)

    })


# --------------------------------------------------
# 4. Create results table
# --------------------------------------------------

df = pd.DataFrame(results)

print("\nPeptideScan analysis results:\n")

print(df.to_string(index=False))


# --------------------------------------------------
# 5. Export results to CSV
# --------------------------------------------------

output_file = (
    project_root /
    "examples" /
    "results_peptidescan.csv"
)

df.to_csv(
    output_file,
    index=False,
    encoding="utf-8-sig"
)

print("\nAnalysis completed successfully!")