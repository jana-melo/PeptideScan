from .fasta import read_fasta
from .properties import calculate_length, calculate_molecular_mass
from .composition import validate_sequence, amino_acid_composition
from .amino_acid_groups import amino_acid_group_composition


def analyze_fasta(file_path):
    """
    Analyze all peptide sequences in a FASTA file.

    This function calculates peptide properties and composition
    without applying any filtering criteria.

    Parameters
    ----------
    file_path : str or Path
        Path to the input FASTA file.

    Returns
    -------
    list of dict
        A list containing the characterization of each peptide.
    """
    # Read peptide sequences from the FASTA file.
    records = read_fasta(file_path)

    analyzed_records = []

    # Analyze each peptide independently.
    for record in records:
        sequence = validate_sequence(record["sequence"])

        # Calculate basic peptide properties.
        length = calculate_length(sequence)
        molecular_mass = calculate_molecular_mass(sequence)

        # Calculate individual amino acid composition.
        composition = amino_acid_composition(sequence)

        # Calculate biochemical group composition.
        group_composition = amino_acid_group_composition(sequence)

        # Store all results while preserving the original FASTA data.
        result = {
            "id": record["id"],
            "description": record["description"],
            "sequence": sequence,
            "length": length,
            "molecular_mass": molecular_mass,
            "amino_acid_composition": composition,
            "group_composition": group_composition,
            "group_frequencies": {
                group: values["frequency"]
                for group, values in group_composition.items()
            }
        }

        analyzed_records.append(result)

    return analyzed_records
