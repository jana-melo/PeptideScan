from collections import Counter


# Aminoácidos padrão e seus códigos de uma letra
AMINO_ACIDS = {
    "A": "Alanina",
    "R": "Arginina",
    "N": "Asparagina",
    "D": "Ácido aspártico",
    "C": "Cisteína",
    "E": "Ácido glutâmico",
    "Q": "Glutamina",
    "G": "Glicina",
    "H": "Histidina",
    "I": "Isoleucina",
    "L": "Leucina",
    "K": "Lisina",
    "M": "Metionina",
    "F": "Fenilalanina",
    "P": "Prolina",
    "S": "Serina",
    "T": "Treonina",
    "W": "Triptofano",
    "Y": "Tirosina",
    "V": "Valina"
}


def validate_sequence(sequence):
    """
    Validate a protein or peptide sequence
    using the 20 standard amino acids.
    """

    sequence = sequence.upper()

    invalid = set(sequence) - set(AMINO_ACIDS)

    if not sequence:
        raise ValueError("Sequence cannot be empty.")

    if invalid:
        raise ValueError(
            f"Invalid amino acid codes: {sorted(invalid)}"
        )

    return sequence


def amino_acid_count(sequence, amino_acid):
    """
    Count occurrences of a specific amino acid.
    """

    sequence = validate_sequence(sequence)
    amino_acid = amino_acid.upper()

    if amino_acid not in AMINO_ACIDS:
        raise ValueError("Invalid amino acid code.")

    return sequence.count(amino_acid)


def amino_acid_composition(sequence):
    """
    Calculate amino acid counts and frequencies.
    """

    sequence = validate_sequence(sequence)

    counts = Counter(sequence)

    results = {}

    for code, name in AMINO_ACIDS.items():

        count = counts.get(code, 0)

        frequency = (count / len(sequence)) * 100

        results[code] = {
            "name": name,
            "count": count,
            "frequency": frequency
        }

    return results


def contains_amino_acids(sequence, amino_acids):
    """
    Check whether a sequence contains
    at least one of the specified amino acids.
    """

    sequence = validate_sequence(sequence)

    amino_acids = [aa.upper() for aa in amino_acids]

    invalid = set(amino_acids) - set(AMINO_ACIDS)

    if invalid:
        raise ValueError(
            f"Invalid amino acid codes: {sorted(invalid)}"
        )

    return any(aa in sequence for aa in amino_acids)