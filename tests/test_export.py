import csv
import pytest

from peptidescan.export import export_to_csv, export_to_fasta


@pytest.fixture
def peptide_records():
    return [
        {
            "id": "P001",
            "description": "P001 Peptide example 1",
            "sequence": "AGKPKY",
            "length": 6,
            "molecular_mass": 662.37,
            "matched_amino_acids": ["K", "Y"],
            "matched_groups": ["hydrophobic"],
            "criteria_mode": "AND",
            "group_frequencies": {
                "hydrophobic": 33.33,
                "aromatic": 16.67
            }
        },
        {
            "id": "P002",
            "description": "P002 Peptide example 2",
            "sequence": "FLWPK",
            "length": 5,
            "molecular_mass": 656.36,
            "matched_amino_acids": ["F"],
            "matched_groups": ["hydrophobic", "aromatic"],
            "criteria_mode": "OR",
            "group_frequencies": {
                "hydrophobic": 80.0,
                "aromatic": 40.0
            }
        }
    ]


def test_export_to_csv(peptide_records, tmp_path):
    output_file = tmp_path / "results.csv"

    result = export_to_csv(peptide_records, output_file)

    assert result.exists()

    with output_file.open(encoding="utf-8-sig", newline="") as file:
        reader = csv.DictReader(file)
        rows = list(reader)

    assert len(rows) == 2
    assert rows[0]["id"] == "P001"
    assert rows[0]["sequence"] == "AGKPKY"
    assert rows[0]["hydrophobic_frequency"] == "33.33"


def test_export_to_fasta(peptide_records, tmp_path):
    output_file = tmp_path / "results.fasta"

    result = export_to_fasta(peptide_records, output_file)

    assert result.exists()

    content = output_file.read_text(encoding="utf-8")

    assert ">P001 Peptide example 1" in content
    assert "AGKPKY" in content
    assert ">P002 Peptide example 2" in content
    assert "FLWPK" in content


def test_export_empty_records_to_csv(tmp_path):
    output_file = tmp_path / "empty.csv"

    export_to_csv([], output_file)

    assert output_file.exists()

    with output_file.open(encoding="utf-8-sig", newline="") as file:
        reader = csv.DictReader(file)
        assert list(reader) == []


def test_export_empty_records_to_fasta(tmp_path):
    output_file = tmp_path / "empty.fasta"

    export_to_fasta([], output_file)

    assert output_file.exists()
    assert output_file.read_text(encoding="utf-8") == ""


def test_export_invalid_records():
    with pytest.raises(TypeError):
        export_to_csv("invalid", "results.csv")


def test_export_record_without_sequence():
    records = [{"id": "P001"}]

    with pytest.raises(ValueError):
        export_to_csv(records, "results.csv")


def test_export_fasta_without_id(tmp_path):
    records = [{"sequence": "AGKPKY"}]
    output_file = tmp_path / "results.fasta"

    with pytest.raises(ValueError):
        export_to_fasta(records, output_file)