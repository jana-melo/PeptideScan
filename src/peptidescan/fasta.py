from Bio import SeqIO


def read_fasta(file_path):
    """
    Read protein or peptide sequences from a FASTA file.

    Parameters
    ----------
    file_path : str
        Path to the FASTA file.

    Returns
    -------
    list
        List of dictionaries containing sequence IDs,
        descriptions and amino acid sequences.
    """

    records = []

    for record in SeqIO.parse(file_path, "fasta"):

        sequence = str(record.seq).upper()

        records.append({
            "id": record.id,
            "description": record.description,
            "sequence": sequence
        })

    return records
