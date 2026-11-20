from pathlib import Path

from .analysis import analyze_fasta
from .export import export_to_csv


def analyze_and_export(file_path, output_path):
    """
    Analyze peptides from a FASTA file and export the results to CSV.

    Parameters
    ----------
    file_path : str or Path
        Path to the input FASTA file.

    output_path : str or Path
        Destination path for the CSV file.

    Returns
    -------
    dict
        Dictionary containing the analysis results and CSV file path.
    """

    # Step 1: Analyze the FASTA file.
    analyzed_records = analyze_fasta(file_path)

    # Step 2: Export the results to CSV.
    csv_path = export_to_csv(analyzed_records, output_path)

    # Step 3: Return the analysis results and output location.
    return {
        "records": analyzed_records,
        "csv_path": Path(csv_path)
    }
    