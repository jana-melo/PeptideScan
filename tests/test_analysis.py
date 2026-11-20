import pytest

from peptidescan.analysis import analyze_fasta


@pytest.fixture
def fasta_file(tmp_path):
    file_path = tmp_path / "peptides.fasta"

    content = (
        ">P001 Peptide example 1\n"
        "AGKPKY\n"
        ">P002 Peptide example 2\n"
        "FLWPK\n"
        ">P003 Peptide example 3\n"
        "VILMGA\n"
    )

    file_path.write_text(content, encoding="utf-8")

    return file_path


def test_analyze_fasta_number_of_records(fasta_file):
    results = analyze_fasta(fasta_file)

    assert len(results) == 3


def test_analyze_fasta_preserves_identifiers(fasta_file):
    results = analyze_fasta(fasta_file)

    assert results[0]["id"] == "P001"
    assert results[0]["description"] == "P001 Peptide example 1"


def test_analyze_fasta_calculates_length(fasta_file):
    results = analyze_fasta(fasta_file)

    assert results[0]["length"] == 6
    assert results[1]["length"] == 5
    assert results[2]["length"] == 6


def test_analyze_fasta_calculates_mass(fasta_file):
    results = analyze_fasta(fasta_file)

    assert all(
        record["molecular_mass"] > 0
        for record in results
    )


def test_analyze_fasta_calculates_composition(fasta_file):
    results = analyze_fasta(fasta_file)

    composition = results[0]["amino_acid_composition"]

    assert composition["A"]["count"] == 1
    assert composition["G"]["count"] == 1
    assert composition["K"]["count"] == 2
    assert composition["P"]["count"] == 1
    assert composition["Y"]["count"] == 1


def test_analyze_fasta_calculates_group_composition(fasta_file):
    results = analyze_fasta(fasta_file)

    group_composition = results[0]["group_composition"]

    assert "hydrophobic" in group_composition
    assert "aromatic" in group_composition

    assert group_composition["hydrophobic"]["frequency"] > 0


def test_analyze_fasta_invalid_sequence(tmp_path):
    file_path = tmp_path / "invalid.fasta"

    file_path.write_text(
        ">P001 Invalid peptide\nAGKX\n",
        encoding="utf-8"
    )

    with pytest.raises(ValueError):
        analyze_fasta(file_path)


def test_analyze_empty_fasta(tmp_path):
    file_path = tmp_path / "empty.fasta"
    file_path.write_text("", encoding="utf-8")

    results = analyze_fasta(file_path)

    assert results == []