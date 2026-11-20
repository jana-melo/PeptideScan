import pytest

from peptidescan.amino_acid_groups import (
    amino_acid_group_composition,
    identify_rich_groups
)


def test_amino_acid_group_composition():

    result = amino_acid_group_composition("AGKPKY")

    assert result["hydrophobic"]["count"] == 2
    assert result["basic"]["count"] == 2
    assert result["aromatic"]["count"] == 1


def test_group_frequency():

    result = amino_acid_group_composition("AGKPKY")

    assert result["hydrophobic"]["frequency"] == pytest.approx(33.3333, rel=1e-3)


def test_identify_rich_groups():

    result = identify_rich_groups("AGKPKY", threshold=30)

    assert "hydrophobic" in result
    assert "basic" in result
    assert "aromatic" not in result


def test_custom_threshold():

    result = identify_rich_groups("AGKPKY", threshold=15)

    assert "aromatic" in result


def test_invalid_threshold():

    with pytest.raises(ValueError):
        identify_rich_groups("AGKPKY", threshold=110)


def test_invalid_sequence():

    with pytest.raises(ValueError):
        amino_acid_group_composition("AGKX")
