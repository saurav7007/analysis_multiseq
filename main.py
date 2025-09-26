import sys
from analysis_multiseq.fasta_utils import *
from analysis_multiseq.logger import get_logger

logger = get_logger()

def main():
    if len(sys.argv) != 2:
        logger.error("Usage: python main.py <multi_fasta_file>")
        sys.exit(1)

    fasta_file = sys.argv[1]

    # Create sequence dictionary
    sequences = seq_dict_from_fasta(fasta_file)

    # Count sequences
    count = count_sequences(sequences)
    logger.info("There are %d sequences in %s", count, fasta_file)


if __name__ == "__main__":
    main()