#__all__ = ["seq_dict_from_fasta", "count_sequences", "count_seq_from_file"]

def seq_dict_from_fasta(multi_fa: str) -> dict[str, str]:
    """
    Create a dictionary containing sequence IDs and their sequences.

    Args:
        multi_fa (str): Path to the FASTA file.

    Returns:
        dict: {sequence_id: sequence}
    """
    seq_dict = {}
    seq_id = None

    with open(multi_fa, "r") as file:
        for line in file:
            line = line.strip()

            if line.startswith(">"):
                seq_id = line[1:]  # remove ">" from ID
                seq_dict[seq_id] = ""
            else:
                if seq_id is None:
                    raise ValueError("FASTA file does not start with a sequence ID")
                seq_dict[seq_id] += line

    return seq_dict


def count_sequences(seq_dictionary: dict) -> int:
    """
    Count the number of sequences in a sequence dictionary.

    Args:
        seq_dictionary (dict): Dictionary of sequences.

    Returns:
        int: Number of sequences.
    """
    return len(seq_dictionary)