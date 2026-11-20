from pyteomics import mass

def calculate_length(sequence):
    """
    Calculate the number of amino acid residues
    in a peptide sequence.

    Parameters
    ----------
    sequence : str
        Amino acid sequence.

    Returns
    -------
    int
        Number of amino acid residues.
    """
    return len(sequence)

def calculate_molecular_mass(sequence):
    """
    Calculate the monoisotopic molecular mass
    of a linear peptide in Daltons (Da).

    Parameters
    ----------
    sequence : str
        Amino acid sequence.

    Returns
    -------
    float
        Molecular mass in Daltons.
    """

    sequence = sequence.upper()

    return mass.calculate_mass(sequence=sequence)
