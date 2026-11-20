import csv

from peptidescan import analyze_and_export

def test_analyze_and_export_creates_csv(tmp_path):
    # Create an example FASTA file.
    fasta_file = tmp_path / "peptides.fasta"

    fasta_file.write_text(
        ">P001 Peptide example 1\n"
        "AGKPKY\n"
        ">P002 Peptide example 2\n"
        "FLWPK\n",
        encoding="utf-8"
    )

    # Define the CSV output path.
    csv_file = tmp_path / "results.csv"

    # Run the complete workflow.
    result = analyze_and_export(fasta_file, csv_file)

    # Check whether the CSV file was created.
    assert csv_file.exists()

    # Check the returned analysis results.
    assert len(result["records"]) == 2
    assert result["csv_path"] == csv_file

    # Read the generated CSV file.
    with csv_file.open(
        "r", newline="", encoding="utf-8-sig"
    ) as file:
        reader = csv.DictReader(file)
        rows = list(reader)

    # Check the number of exported peptides.
    assert len(rows) == 2

    # Check peptide identification and sequence.
    assert rows[0]["id"] == "P001"
    assert rows[0]["sequence"] == "AGKPKY"

    # Check individual amino acid composition.
    assert int(rows[0]["A_count"]) == 1
    assert int(rows[0]["K_count"]) == 2

    # Check biochemical group information.
    assert rows[0]["hydrophobic_frequency"] != ""
