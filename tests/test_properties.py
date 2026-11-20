from peptidescan.properties import (
    calculate_length,
    calculate_molecular_mass
)


def test_calculate_length():
    sequence = "AGK"

    result = calculate_length(sequence)

    assert result == 3


def test_calculate_molecular_mass():
    sequence = "AGK"

    result = calculate_molecular_mass(sequence)

    assert isinstance(result, float)
    assert result > 0
