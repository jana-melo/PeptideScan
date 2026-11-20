from .composition import validate_sequence


# Operational classification of amino acids.
# Some groups overlap, such as hydrophobic and aromatic.
#
# These categories represent an operational classification
# adopted for this package and are not mutually exclusive.

AMINO_ACID_GROUPS = {
    "hydrophobic": {
        "A", "V", "I", "L", "M", "F", "W", "P"
    },

    "polar_uncharged": {
        "S", "T", "N", "Q", "C", "G", "Y"
    },

    "acidic": {
        "D", "E"
    },

    "basic": {
        "K", "R", "H"
    },

    "aromatic": {
        "F", "Y", "W"
    }
}


def amino_acid_group_composition(sequence):
    """
    Calculate the count and frequency of amino acid groups.

    Parameters
    ----------
    sequence : str
        Peptide amino acid sequence.

    Returns
    -------
    dict
        Counts and frequencies (%) for each amino acid group.
    """

    sequence = validate_sequence(sequence)

    results = {}

    for group, amino_acids in AMINO_ACID_GROUPS.items():

        count = sum(
            sequence.count(amino_acid)
            for amino_acid in amino_acids
        )

        frequency = (count / len(sequence)) * 100

        results[group] = {
            "count": count,
            "frequency": frequency
        }

    return results


def identify_rich_groups(sequence, threshold=30.0):
    """
    Identify amino acid groups with frequencies
    equal to or greater than a specified threshold.

    Parameters
    ----------
    sequence : str
        Peptide amino acid sequence.

    threshold : float
        Minimum frequency (%) to classify a group
        as enriched. Default is 30%.

    Returns
    -------
    list
        Amino acid groups meeting the threshold.
    """

    if isinstance(threshold, bool) or not isinstance(
        threshold, (int, float)
    ):
        raise TypeError("Threshold must be a number.")

    if not 0 <= threshold <= 100:
        raise ValueError("Threshold must be between 0 and 100.")

    composition = amino_acid_group_composition(sequence)

    rich_groups = [
        group
        for group, values in composition.items()
        if values["frequency"] >= threshold
    ]

    return rich_groups


def evaluate_group_criteria(sequence, thresholds, mode="AND"):
    """
    Evaluate a peptide against user-defined amino acid
    group frequency thresholds.

    Parameters
    ----------
    sequence : str
        Peptide amino acid sequence.

    thresholds : dict
        Dictionary containing biochemical groups and
        their minimum frequency thresholds.

        Example:
        {"hydrophobic": 30, "aromatic": 15}

    mode : str
        "AND": all selected groups must meet their thresholds.
        "OR": at least one selected group must meet its threshold.

    Returns
    -------
    dict
        Evaluation results, including:
        - passes: whether the peptide meets the criteria.
        - matched_groups: groups meeting their thresholds.
        - group_frequencies: calculated frequencies.
        - mode: selected evaluation mode.
    """

    # Validate thresholds
    if not isinstance(thresholds, dict) or not thresholds:
        raise ValueError(
            "Provide at least one group and its frequency threshold."
        )

    # Validate evaluation mode
    if not isinstance(mode, str):
        raise TypeError("Mode must be a string.")

    mode = mode.upper()

    if mode not in ("AND", "OR"):
        raise ValueError("Mode must be 'AND' or 'OR'.")

    # Validate group names
    invalid_groups = set(thresholds) - set(AMINO_ACID_GROUPS)

    if invalid_groups:
        raise ValueError(
            f"Invalid amino acid groups: {sorted(invalid_groups)}"
        )

    # Validate threshold values
    for group, threshold in thresholds.items():

        if isinstance(threshold, bool) or not isinstance(
            threshold, (int, float)
        ):
            raise TypeError(
                f"Threshold for {group} must be a number."
            )

        if not 0 <= threshold <= 100:
            raise ValueError(
                f"Threshold for {group} must be between 0 and 100."
            )

    # Calculate group composition using the existing function
    composition = amino_acid_group_composition(sequence)

    frequencies = {
        group: composition[group]["frequency"]
        for group in thresholds
    }

    # Identify groups that meet their respective thresholds
    matched_groups = [
        group
        for group, threshold in thresholds.items()
        if frequencies[group] >= threshold
    ]

    # Evaluate the selected groups using AND or OR logic
    if mode == "AND":
        passes = len(matched_groups) == len(thresholds)

    else:
        passes = len(matched_groups) > 0
    
    return {
        "passes": passes,
        "matched_groups": matched_groups,
        "group_frequencies": frequencies,
        "mode": mode
    }
