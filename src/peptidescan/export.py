import csv
from pathlib import Path


def _validate_records(records):
    """
    Validate the list of peptide records before export.
    """
    if not isinstance(records, list):
        raise TypeError("Records must be provided as a list.")

    for index, record in enumerate(records):
        if not isinstance(record, dict):
            raise TypeError(f"Record {index} must be a dictionary.")

        if "sequence" not in record:
            raise ValueError(f"Record {index} is missing the sequence field.")

        if not isinstance(record["sequence"], str) or not record["sequence"]:
            raise ValueError(f"Record {index} has an invalid sequence.")


def export_to_csv(records, file_path):
    """
    Export peptide analysis results to a CSV file.

    Includes basic peptide properties, individual amino acid
    composition, and biochemical group composition.
    """
    import csv
    from pathlib import Path

    from .amino_acid_groups import AMINO_ACID_GROUPS

    _validate_records(records)

    file_path = Path(file_path)

    # Standard amino acid codes, in a fixed order.
    amino_acid_codes = "ACDEFGHIKLMNPQRSTVWY"

    # Basic information about each peptide.
    fieldnames = [
        "id",
        "description",
        "sequence",
        "length",
        "molecular_mass",
        "matched_amino_acids",
        "matched_groups",
        "criteria_mode"
    ]

    # Add count and frequency columns for each amino acid.
    for amino_acid in amino_acid_codes:
        fieldnames.extend([
            f"{amino_acid}_count",
            f"{amino_acid}_frequency"
        ])

    # Add count and frequency columns for each biochemical group.
    for group in AMINO_ACID_GROUPS:
        fieldnames.extend([
            f"{group}_count",
            f"{group}_frequency"
        ])

    with file_path.open("w", newline="", encoding="utf-8-sig") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()

        for record in records:
            row = {
                "id": record.get("id", ""),
                "description": record.get("description", ""),
                "sequence": record.get("sequence", ""),
                "length": record.get("length", ""),
                "molecular_mass": record.get("molecular_mass", ""),
                "matched_amino_acids": ", ".join(
                    record.get("matched_amino_acids") or []
                ),
                "matched_groups": ", ".join(
                    record.get("matched_groups") or []
                ),
                "criteria_mode": record.get("criteria_mode", "")
            }

            # Retrieve individual amino acid composition.
            amino_acid_composition = record.get(
                "amino_acid_composition", {}
            )

            for amino_acid in amino_acid_codes:
                values = amino_acid_composition.get(amino_acid, {})
                row[f"{amino_acid}_count"] = values.get("count", "")
                row[f"{amino_acid}_frequency"] = values.get(
                    "frequency", ""
                )

            # Retrieve biochemical group composition.
            group_composition = record.get("group_composition", {})
            group_frequencies = record.get("group_frequencies", {})

            for group in AMINO_ACID_GROUPS:
                values = group_composition.get(group, {})

                row[f"{group}_count"] = values.get("count", "")
                row[f"{group}_frequency"] = values.get(
                    "frequency",
                    group_frequencies.get(group, "")
                )

            writer.writerow(row)

    return file_path

def export_to_fasta(records, file_path):
    """
    Export peptide sequences to a FASTA file.

    Parameters
    ----------
    records : list of dict
        Peptide records, usually returned by filter_peptides().
    file_path : str or Path
        Destination path for the FASTA file.

    Returns
    -------
    Path
        Path to the generated FASTA file.
    """
    _validate_records(records)

    file_path = Path(file_path)

    with file_path.open("w", encoding="utf-8") as fasta_file:
        for index, record in enumerate(records):
            record_id = record.get("id")

            if not record_id:
                raise ValueError(
                    f"Record {index} is missing an identifier (id)."
                )

            # Preserve the original FASTA description when available.
            header = record.get("description") or record_id

            fasta_file.write(f">{header}\n")

            sequence = record["sequence"]

            # Write sequences in lines of up to 60 residues.
            for start in range(0, len(sequence), 60):
                fasta_file.write(sequence[start:start + 60] + "\n")

    # Return the path of the generated file.
    return file_path