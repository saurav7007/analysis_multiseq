from collections import Counter

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
                seq_id = line[1:].split(" ")[0]  # remove ">" from ID
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


def len_extrema(sequences_lens: dict, extrema: str = "min") -> tuple[int, list[str]]:
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
        list: [(orf1_pos, orf1_seq, orf1_len)]
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
                
    return [(orf_pos[0], sequence[orf_pos[0] - 1 : orf_pos[1] + 2], orf_pos[1] - orf_pos[0]) for orf_pos in valid_orfs_pos]


def orf_dict(seq_dictionary: dict, frames: list = [1, 2, 3]) -> dict[str, list[tuple[int, str, int]]]:
    """
    Converts seq_dictionary into orf dictionary containing sequence_id and list of ORFs.

    Args:
        seq_dictionary (dict): Dictionary of sequences {seq_id: seq}.
        frames (list): List of frames to scan (1, 2, 3).

    Returns:
        dict: {sequence_id: [(orf1_pos, orf1_seq, orf1_len), (orf2_pos, orf2_seq, orf2_len), ...]}
    """
    orf_dict = {}

    for seq_id, sequence in seq_dictionary.items():
        all_orf = []
        for frame in frames:
            orf_in_frame = orf_finder(sequence = sequence, frame = frame)
            all_orf.extend(orf_in_frame)
        orf_dict[seq_id] = all_orf
        
    return orf_dict


def orf_extremas(orf_dictionary: dict, extrema: str = "min") -> dict[str, list[tuple[int, str, int]]]:
    """
    Find ORFs with maximum and minimum lengths for a sequence_id.

    Args:
        orf_dictionary (dict): Dictionary of ORFS {sequence_id: [(orf1_pos, orf1_seq, orf1_len), (orf2_pos, orf2_seq, orf2_len), ...]}
        extreme (str): "min" or "max"

    Returns:
        dict: {sequence_id: [(orf1_pos, orf1_seq, orf1_len), orf2_pos, orf2_seq, orf2_len, ...]}
    """
    extremas = {}

    for seq_id, orf_list in orf_dictionary.items():
        if orf_list == []:
            extremas[seq_id] = []
            continue
        elif extrema == "min":
            extrema_len = min(orf[2] for orf in orf_list)
        elif extrema == "max":
            extrema_len = max(orf[2] for orf in orf_list)
        else:
            raise ValueError("extrema must be 'min' or 'max'")
        
        orf_extrema = [orf for orf in orf_list if orf[2] == extrema_len]

        extremas[seq_id] = orf_extrema
    
    return extremas


def longest_orf(orf_dictionary: dict) -> dict[str, list[tuple[int, str, int]]]:
    """
    Get the longest ORF(s) across all sequences.

    Args:
        orf_dictionary (dict): Dictionary of ORFs
        {sequence_id: [(orf1_pos, orf1_seq, orf1_len), ...]}

    Returns:
        dict: {sequence_id: [(longest_orf_pos, longest_orf_seq, longest_orf_len), ...]}
              Includes all sequences that contain the globally longest ORF(s).
    """
    seq_longest_orfs = orf_extremas(orf_dictionary, "max")

    global_max_len = max(orfs[0][2] for orfs in seq_longest_orfs.values() if orfs)

    # Keep only those sequences whose ORFs have the global max length
    global_longest_orfs = {seq_id: orfs for seq_id, orfs in seq_longest_orfs.items() if orfs and orfs[0][2] == global_max_len}

    return global_longest_orfs


def shortest_orf(orf_dictionary: dict) -> dict[str, list[tuple[int, str, int]]]:
    """
    Get the shortest orf from the orf_dictionary.

    Args:
        orf_dictionary (dict): Dictionary of ORFS {sequence_id: [(orf1_pos, orf1_seq, orf1_len), (orf2_pos, orf2_seq, orf2_len), ...]}

    Returns:
        dict: {sequence_id: [(shortest_orf_pos, shortest_orf_seq, shortest_orf_len)]}
    """
    seq_shortest_orfs = orf_extremas(orf_dictionary, "min")

    global_min_len = min(orfs[0][2] for orfs in seq_shortest_orfs.values() if orfs)

    # Keep only those sequences whose ORFs have the global max length
    global_shortest_orfs = {seq_id: orfs for seq_id, orfs in seq_shortest_orfs.items() if orfs and orfs[0][2] == global_min_len}

    return global_shortest_orfs


def scan_motif(seq_dictionary: dict, size: int) -> dict[str , list[tuple[int, str]]]:
    """
    Find all the repeating region in a given sequence dictionary along with their count.

    Args:
        seq_dictionary (dict): Dictionary of sequences: {seq_id: seq}.
        size (int): The size of the motif to search.

    Returns:
        dict: {sequence_id: [(motif1_count, motif1_seq), (motif2_count, motif2_seq), ...]}
    """
    motifs = {}

    for seq_id, seq in seq_dictionary.items():
        k_mers = [seq[i:i+size] for i in range(0, (len(seq))) if len(seq) - i >= size]
        counts = Counter(k_mers)
        motif_list = [(count, k_mer) for k_mer, count in counts.items() if count > 1]
        motifs[seq_id] = motif_list
    
    return size, motifs


def find_common_motif(motif_dictionary: dict, size: int) -> dict[str , list[tuple[int, str]]]:
    """
    Find most common motif per sequence in sequence dictionary.

    Args:
        seq_dictionary (dict): Dictionary of sequences: {seq_id: seq}.
        size (int): The size of the motif to search.
    
    Returns:
        dict: {seq1: [(seq1_common_motif_count, seq1_common_motif_seq)], seq2: [(seq2_common_motif_count, seq1_common_motif_seq)], ...}
    """
    size, motifs = scan_motif(motif_dictionary, size)

    common_motif = {}

    for seq_id, motif_list in motifs.items():
        max_freq = max(motif[0] for motif in motif_list)
        common_motif[seq_id] = [motif for motif in motif_list if motif[0] == max_freq]

    return size, common_motif


def most_common_motif(motif_dictionary: dict, size: int) -> dict[str , list[tuple[int, str]]]:
    """
    Find most common motif across all the sequence in sequence dictionary.

    Args:
        orf_dictionary (dict): Dictionary of ORFS {sequence_id: [(orf1_pos, orf1_seq, orf1_len), (orf2_pos, orf2_seq, orf2_len), ...]}

    Returns:
        dict: {sequence_id: [(most_common_motif_count, most_common_motif_seq), ...]}
    """
    size, common_motifs = find_common_motif(motif_dictionary, size)

    global_max_count = max(motif[0][0] for motif in common_motifs.values())

    # Keep only those motifs with global max count
    most_common_motif = {seq_id: motif_list for seq_id, motif_list in common_motifs.items() if motif_list[0][0] == global_max_count}

    return size, most_common_motif