
# analysis_multiseq

A Python package for analyzing DNA sequence data, with utilities for sequence length, ORF, and motif analysis—optimized for multi-sequence FASTA files.

## Features

- Read and process multi-FASTA files
- Analyze DNA sequence length, open reading frames (ORFs), and motifs
- Logging support for analysis workflows

## Project Structure

```
analysis_multiseq/
    __init__.py
    fasta_utils.py   # FASTA file utilities
    logger.py        # Logging utilities
main.py             # Example usage or entry point
dna.example.fasta   # Example FASTA file
dna2.fasta          # Another example FASTA file
```

git clone https://github.com/saurav7007/analysis_multiseq.git
## Installation

Clone the repository:

```bash
git clone https://github.com/saurav7007/analysis_multiseq.git
cd analysis_multiseq
```

(Optional) Create and activate a virtual environment for isolation.

## Usage

To analyze a FASTA file, run:

```bash
python3 main.py dna2.fasta
```

See `main.py` for more usage examples.

## Requirements

- Python 3.7 or higher

## License

This project is licensed under the MIT License.

## Author

[Saurav](https://github.com/saurav7007)
