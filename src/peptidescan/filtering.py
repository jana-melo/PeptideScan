from .properties import (
    calculate_length,
    calculate_molecular_mass
)

from .composition import validate_sequence

from .amino_acid_groups import (
    evaluate_group_criteria
)


def filter_peptides(
    records,
    min_length=None,
    max_length=None,
    min_mass=None,
    max_mass=None,
    required_amino_acids=None,
    amino_acid_mode="OR",
    group_thresholds=None,
    group_mode="AND",
    criteria_mode="AND"
):
    """
    Filter peptide records using user-defined criteria.

    Parameters
    ----------
    records : list
        List of dictionaries returned by read_fasta().
        Each dictionary must contain a "sequence" key.

    min_length, max_length : int, optional
        Minimum and maximum peptide length.

    min_mass, max_mass : float, optional
        Minimum and maximum monoisotopic molecular mass in Da.

    required_amino_acids : list, optional
        Amino acid codes required in the peptide.

    amino_acid_mode : str
        "OR": at least one specified amino acid must be present.
        "AND": all specified amino acids must be present.

    group_thresholds : dict, optional
        Dictionary with biochemical groups and minimum
        frequency thresholds.
        Example: {"hydrophobic": 30, "aromatic": 15}

    group_mode : str
        "AND" or "OR" for the selected biochemical groups.

    criteria_mode : str
        "AND": all active filter categories must be satisfied.
        "OR": at least one active filter category must be satisfied.

    Returns
    -------
    list
        Filtered peptide records with calculated properties.
    """

    # 1. Validate general input

    if not isinstance(records, list):
        raise TypeError("Records must be provided as a list.")

    if not isinstance(criteria_mode, str):
        raise TypeError("Criteria mode must be a string.")

    if not isinstance(amino_acid_mode, str):
        raise TypeError("Amino acid mode must be a string.")

    criteria_mode = criteria_mode.upper()
    amino_acid_mode = amino_acid_mode.upper()

    if criteria_mode not in ("AND", "OR"):
        raise ValueError("Criteria mode must be 'AND' or 'OR'.")

    if amino_acid_mode not in ("AND", "OR"):
        raise ValueError("Amino acid mode must be 'AND' or 'OR'.")

    # 2. Validate length limits

    for name, value in [
        ("min_length", min_length),
        ("max_length", max_length)
    ]:

        if value is not None:

            if isinstance(value, bool) or not isinstance(value, int):
                raise TypeError(f"{name} must be an integer.")

            if value < 1:
                raise ValueError(f"{name} must be at least 1.")

    if (
        min_length is not None
        and max_length is not None
        and min_length > max_length
    ):
        raise ValueError("Minimum length cannot exceed maximum length.")

    # 3. Validate molecular mass limits

    for name, value in [
        ("min_mass", min_mass),
        ("max_mass", max_mass)
    ]:

        if value is not None:

            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise TypeError(f"{name} must be a number.")

            if value <= 0:
                raise ValueError(f"{name} must be greater than zero.")

    if (
        min_mass is not None
        and max_mass is not None
        and min_mass > max_mass
    ):
        raise ValueError("Minimum mass cannot exceed maximum mass.")

    # 4. Validate required amino acids

    if required_amino_acids is not None:

        if not isinstance(required_amino_acids, (list, tuple)):
            raise TypeError(
                "Required amino acids must be a list or tuple."
            )

        normalized_amino_acids = []

        for aa in required_amino_acids:

            if not isinstance(aa, str):
                raise TypeError("Amino acid codes must be strings.")

            aa = aa.upper()

            if len(aa) != 1:
                raise ValueError(
                    "Each required amino acid must have one letter."
                )

            validate_sequence(aa)

            normalized_amino_acids.append(aa)

        # Remove repeated amino acid codes
        required_amino_acids = list(
            dict.fromkeys(normalized_amino_acids)
        )

    # 5. Prepare the result list

    filtered_records = []

    # 6. Analyze each peptide

    for record in records:

        if not isinstance(record, dict):
            raise TypeError("Each record must be a dictionary.")

        if "sequence" not in record:
            raise ValueError("Each record must contain a sequence.")

        sequence = validate_sequence(record["sequence"])

        # Reuse existing functions
        length = calculate_length(sequence)

        molecular_mass = calculate_molecular_mass(sequence)

        conditions = []

        # 7. Length criterion

        if min_length is not None or max_length is not None:

            length_passes = True

            if min_length is not None:
                length_passes = length_passes and length >= min_length

            if max_length is not None:
                length_passes = length_passes and length <= max_length

            conditions.append(length_passes)

        # 8. Molecular mass criterion

        if min_mass is not None or max_mass is not None:

            mass_passes = True

            if min_mass is not None:
                mass_passes = mass_passes and molecular_mass >= min_mass

            if max_mass is not None:
                mass_passes = mass_passes and molecular_mass <= max_mass

            conditions.append(mass_passes)

        # 9. Required amino acid criterion

        matched_amino_acids = []

        if required_amino_acids:

            matched_amino_acids = [
                aa for aa in required_amino_acids
                if aa in sequence
            ]

            if amino_acid_mode == "AND":

                amino_acids_pass = (
                    len(matched_amino_acids) == len(required_amino_acids)
                )

            else:

                amino_acids_pass = len(matched_amino_acids) > 0

            conditions.append(amino_acids_pass)

        # 10. Biochemical group criterion

        group_evaluation = {
            "passes": True,
            "matched_groups": [],
            "group_frequencies": {},
            "mode": group_mode.upper() if isinstance(group_mode, str) else group_mode
        }

        if group_thresholds is not None:

            group_evaluation = evaluate_group_criteria(
                sequence,
                group_thresholds,
                mode=group_mode
            )

            conditions.append(group_evaluation["passes"])

        # 11. Combine all active criteria

        if not conditions:

            passes = True

        elif criteria_mode == "AND":

            passes = all(conditions)

        else:

            passes = any(conditions)

        # 12. Save selected peptide and calculated properties

        if passes:

            result = record.copy()

            result["sequence"] = sequence
            result["length"] = length
            result["molecular_mass"] = molecular_mass

            result["matched_amino_acids"] = matched_amino_acids

            result["matched_groups"] = (
                group_evaluation["matched_groups"]
            )

            result["group_frequencies"] = (
                group_evaluation["group_frequencies"]
            )

            result["criteria_mode"] = criteria_mode

            filtered_records.append(result)

    # 13. Return the filtered records

    return filtered_records
