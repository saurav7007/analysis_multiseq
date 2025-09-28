import sys
from analysis_multiseq import logger as log, fasta_utils as fs

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
    min_len, shortest_seqs = fs.find_extrema(seq_len, "min")
    max_len, longest_seqs = fs.find_extrema(seq_len, "max")

    log.info(f"There are {len(shortest_seqs)} shortest sequences with ids {shortest_seqs} and length {min_len}.")
    log.info(f"There are {len(shortest_seqs)} longest sequences with ids {longest_seqs} and length {max_len}.")


if __name__ == "__main__":
    main()