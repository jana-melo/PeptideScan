import pytest

from peptidescan.composition import (
    amino_acid_count,
    amino_acid_composition,
    contains_amino_acids,
    validate_sequence
)


def test_amino_acid_count():
    assert amino_acid_count("AGKPKY", "K") == 2


def test_amino_acid_composition():

    result = amino_acid_composition("AGKPKY")

    assert result["K"]["count"] == 2
    assert result["A"]["count"] == 1
    assert result["Y"]["count"] == 1


def test_amino_acid_frequency():

    result = amino_acid_composition("AGKPKY")

    assert result["K"]["frequency"] == pytest.approx(33.3333, rel=1e-3)


def test_contains_amino_acids():

    assert contains_amino_acids("AGKPKY", ["K", "R"]) is True

    assert contains_amino_acids("AGKPKY", ["D", "E"]) is False


def test_invalid_amino_acid():

    with pytest.raises(ValueError):
        amino_acid_count("AGKPKY", "B")


def test_invalid_sequence():

    with pytest.raises(ValueError):
        validate_sequence("AGKX")
