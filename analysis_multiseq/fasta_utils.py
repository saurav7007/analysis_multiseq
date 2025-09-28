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


def sequences_lens(seq_dictionary: dict) -> dict[str, int]:
    """
    Count the length of every sequence in a sequence dictionary.

    Args:
        seq_dictionary (dict): Dictionary of sequences.

    Returns:
        dict: {sequence_id: sequence length}.
    """
    return {seq_id: len(seq) for seq_id, seq in seq_dictionary.items()}


def find_extrema(sequences_lens: dict, extrema: str = "min") -> tuple[int, list[str]]:
    """
    Find the sequences with maximum or minimun length.

    Args:
        seq_lens (dict): {seq_id: len}
        extreme (str): "min" or "max"

    Returns:
        tuple: (extreme_len, [seq_ids])
    """
    if extrema == "min":
        extreme_len = min(sequences_lens.values())
    elif extrema == 'max':
        extreme_len = max(sequences_lens.values())
    else:
        raise ValueError("extrema must be 'min' or 'max'")
    
    seq_ids = [seq_id for seq_id, length in sequences_lens.items() if length == extreme_len]

    return (extreme_len, seq_ids)