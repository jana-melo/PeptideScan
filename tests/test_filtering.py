import pytest

from peptidescan.filtering import filter_peptides
from peptidescan.properties import calculate_molecular_mass


@pytest.fixture
def peptide_records():
    return [
        {"id": "P001", "sequence": "AGKPKY"},
        {"id": "P002", "sequence": "FLWPK"},
        {"id": "P003", "sequence": "VILMGA"},
        {"id": "P004", "sequence": "DEDEK"},
    ]


def test_filter_by_min_length(peptide_records):
    result = filter_peptides(
        peptide_records,
        min_length=6
    )

    assert all(record["length"] >= 6 for record in result)
    assert len(result) == 2


def test_filter_by_length_range(peptide_records):
    result = filter_peptides(
        peptide_records,
        min_length=5,
        max_length=6
    )

    assert all(5 <= record["length"] <= 6 for record in result)


def test_filter_by_mass(peptide_records):
    sequence = "AGKPKY"
    expected_mass = calculate_molecular_mass(sequence)

    result = filter_peptides(
        peptide_records,
        min_mass=expected_mass,
        max_mass=expected_mass
    )

    assert len(result) == 1
    assert result[0]["id"] == "P001"


def test_filter_by_amino_acid_or(peptide_records):
    result = filter_peptides(
        peptide_records,
        required_amino_acids=["F", "D"],
        amino_acid_mode="OR"
    )

    assert {record["id"] for record in result} == {"P002", "P004"}


def test_filter_by_amino_acid_and(peptide_records):
    result = filter_peptides(
        peptide_records,
        required_amino_acids=["K", "Y"],
        amino_acid_mode="AND"
    )

    assert {record["id"] for record in result} == {"P001"}


def test_filter_by_biochemical_groups(peptide_records):
    result = filter_peptides(
        peptide_records,
        group_thresholds={"hydrophobic": 30},
        group_mode="AND"
    )

    assert all(
        record["group_frequencies"]["hydrophobic"] >= 30
        for record in result
    )


def test_combined_criteria_and(peptide_records):
    result = filter_peptides(
        peptide_records,
        min_length=5,
        required_amino_acids=["K"],
        criteria_mode="AND"
    )

    assert all(
        record["length"] >= 5 and "K" in record["sequence"]
        for record in result
    )


def test_combined_criteria_or(peptide_records):
    result = filter_peptides(
        peptide_records,
        min_length=6,
        required_amino_acids=["F"],
        criteria_mode="OR"
    )

    assert {record["id"] for record in result} == {
        "P001", "P002", "P003"
    }


def test_invalid_length_range(peptide_records):
    with pytest.raises(ValueError):
        filter_peptides(
            peptide_records,
            min_length=10,
            max_length=5
        )


def test_invalid_criteria_mode(peptide_records):
    with pytest.raises(ValueError):
        filter_peptides(
            peptide_records,
            min_length=5,
            criteria_mode="XOR"
        )


def test_invalid_sequence():
    records = [{"id": "P005", "sequence": "AGKX"}]

    with pytest.raises(ValueError):
        filter_peptides(records)


def test_empty_criteria(peptide_records):
    result = filter_peptides(peptide_records)

    assert len(result) == len(peptide_records)