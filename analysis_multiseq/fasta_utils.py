def read_fa(multi_fa: str) -> list[str]:
    """
    Read the lines of a multi-FASTA file.

    Args:
        multi_fa (str): Path to the FASTA file.

    Returns:
        list[str]: Lines of the file.
    """
    with open(multi_fa, "r") as file:
        return file.readlines()


def count_seq(multi_fa: str) -> int:
    """
    Count the number of sequences in a multi-FASTA file.

    Args:
        multi_fa (str): Path to the FASTA file.

    Returns:
        int: Number of sequences.
    """
    count = 0
    lines = read_fa(multi_fa)
    for line in lines:
        if line.startswith(">"):
            count += 1
    return count