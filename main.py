import sys
from analysis_multiseq.fasta_utils import count_seq
from analysis_multiseq.logger import get_logger

logger = get_logger()

def main():
    if len(sys.argv) != 2:
        logger.error("Usage: python main.py <multi_fasta_file>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    count = count_seq(fasta_file)
    logger.info("There are %d sequences in %s", count, fasta_file)

if __name__ == "__main__":
    main()