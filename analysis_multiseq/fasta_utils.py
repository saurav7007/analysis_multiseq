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


def orf_finder(sequence: str = "ATGCCCATGTAGAATTCAATGTTATAG", start_codons: list = ["ATG"], stop_codons: list = ["TAA","TAG","TGA"], frame: int = 1) -> list:
    """
    Find all the ORFs present in the sequence.

    Args:
        A string of sequence.

    Returns:
        list: A list of ORFs per sequence.
    """
    if frame not in {1, 2, 3}:
        raise ValueError("Frame must be 1, 2, or 3")
    updated_seq = sequence[frame - 1:]

    codons = [updated_seq[i:i+3] for i in range(0, len(updated_seq) - len(updated_seq) % 3, 3)]

    start_codon_pos = [pos * 3 + frame for pos, codon in enumerate(codons) if codon in start_codons]
    stop_codon_pos = [pos * 3 + frame for pos, codon in enumerate(codons) if codon in stop_codons]

    valid_orfs_pos = []

    for start in start_codon_pos:
        for stop in stop_codon_pos:
            if start < stop:
                valid_orfs_pos.append((start, stop))
                break
                
    return [(orf_pos[0], sequence[orf_pos[0] - 1 : orf_pos[1] + 2]) for orf_pos in valid_orfs_pos]


def orf_dict(seq_dictionary: dict, frames: list = [1, 2, 3]) -> dict[str, list]:
        """
        Converts seq_dictionary into orf dictionary containing sequence_id and list of ORFs.

        Args:
            seq_dictionary (dict): Dictionary of sequences {seq_id: seq}.
            frames (list): List of frames to scan (1, 2, 3).

        Returns:
            dict: {sequence_id: [(orf1_pos, orf1_seq), (orf2_pos, orf2_seq), ...]}
        """
        orf_dict = {}

        for seq_id, sequence in seq_dictionary.items():
            all_orf = []
            for frame in frames:
                orf_in_frame = orf_finder(sequence = sequence, frame = frame)
                all_orf.extend(orf_in_frame)
            orf_dict[seq_id] = all_orf
        
        return orf_dict