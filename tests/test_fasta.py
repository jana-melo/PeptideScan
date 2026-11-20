from peptidescan.fasta import read_fasta

def test_read_fasta():
    records = read_fasta("examples/examples_proteins.fasta")
    # Verify the number of sequences
    assert len(records) == 50
    # Verify that sequences were read
    assert all(record["sequence"] for record in records)
# Print results
    for record in records:
        print(
            record["id"],
            record["sequence"],
            len(record["sequence"])
        )
