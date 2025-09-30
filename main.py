import sys
from analysis_multiseq import logger, fasta_utils as fs

log = logger.get_logger()

def main():
    if len(sys.argv) != 2:
        logger.error("Usage: python main.py <multi_fasta_file>")
        sys.exit(1)

    fasta_file = sys.argv[1]

    # Create sequence dictionary
    sequences = fs.seq_dict_from_fasta(fasta_file)

    # Count sequences
    seq_count = fs.count_sequences(sequences)

    log.info(f"There are {seq_count} sequences in {fasta_file}")

    # Measure the length of every sequence
    seq_len = fs.sequences_lens(sequences)
    min_len, shortest_seqs = fs.len_extrema(seq_len, "min")
    max_len, longest_seqs = fs.len_extrema(seq_len, "max")

    log.info(f"There are {len(shortest_seqs)} shortest sequences with ids {shortest_seqs} and length {min_len}.")
    log.info(f"There are {len(shortest_seqs)} longest sequences with ids {longest_seqs} and length {max_len}.")

    # Measure all the ORFs in sequence dictionary
    orfs = fs.orf_dict(sequences, frames = [1, 2, 3])
    long_orfs = fs.orf_extremas(orfs, "max")

    seq_id = "gi|142022655|gb|EQ086233.1|16"
    log.info(f"The longest forward ORF in sequence {seq_id} is {long_orfs[seq_id]}.")

    longest_orf = fs.longest_orf(orfs)

    log.info(f"{len(longest_orf)} sequences {list(longest_orf.keys())} are found in {fasta_file} with longest ORFs with length {list(longest_orf.values())[0][0][2]} at position {list(longest_orf.values())[0][0][0]}.")

    # Detect the motif frequency
    all_motifs = fs.scan_motif(sequences, size=7)

    common_motifs = fs.find_common_motif(sequences, size=7)

    motif_size, most_common_motifs = fs.most_common_motif(sequences, size=7)

    print(most_common_motifs)

    log.info(f"The most frequent motif of length {motif_size} is {list(most_common_motifs.values())[0][0][1]}. It occurred {list(most_common_motifs.values())[0][0][0]} times in sequence {list(longest_orf.keys())}")

if __name__ == "__main__":
    main()