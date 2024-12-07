# Sequence Alignment Tool

## Overview
This project provides a Python-based command-line tool to compute optimal alignments and scores for pairs of amino acid sequences using either the Needleman-Wunsch (global alignment) or Smith-Waterman (local alignment) algorithms. The tool supports affine gap penalties for local alignment and linear gap penalties for global alignment.

## Features
- **Global alignment** using the Needleman-Wunsch algorithm.
- **Local alignment** using the Smith-Waterman algorithm.

## Requirements
- Python 3.6 or later.

## Installation
1. Clone this repository:
   ```bash
   git clone <repository_url>
   cd <repository_folder>
   ```
2. Ensure Python and the required dependencies are installed.

## Usage
The tool is run via the command line and requires specific input parameters. Below is the basic syntax:

```bash
python main.py <sequence_one> <sequence_two> <scoring_matrix_file> <gap_penalty> --mode <alignment_mode> [--e <extension_penalty>]
```

### Positional Arguments
1. **`sequence_one`**: The first sequence for alignment.
2. **`sequence_two`**: The second sequence for alignment.
3. **`scoring_matrix_file`**: Path to the scoring matrix file in BLAST format.
4. **`gap_penalty`**: The gap penalty for alignment (integer).

### Optional Arguments
1. **`--mode`** *(required)*:
   - `global`: For global alignment (Needleman-Wunsch).
   - `local`: For local alignment (Smith-Waterman).
2. **`--e`** *(required only in local mode)*: The gap extension penalty (integer).

### Example Commands
#### Global Alignment
```bash
python main.py ACGTACGT TGCATGCA scoring_matrix.txt 5 --mode global
```

#### Local Alignment
```bash
python main.py ACGTACGT TGCATGCA scoring_matrix.txt 5 --mode local --e 2
```

## Input Requirements
### Scoring Matrix
- Scoring matrices must be provided in **BLAST format**.
- The matrix defines match/mismatch scores for amino acids.

### Sequences
- Input sequences must be strings composed of valid amino acid or nucleotide symbols.

## Output
The program outputs:
1. The optimal alignment score.
2. The number of optimal alignments found.
3. Each optimal alignment displayed in a two-line format.

### Example Output
```plaintext
Score: 13
Found 4 optimal global alignments:
-ACGTACGT
TGCATGC-A

A-CGTACGT
TGCATGC-A

-ACGTACGT
TGCATGCA-

A-CGTACGT
TGCATGCA-
```

## Contributing
Contributions are welcome! Feel free to submit issues or pull requests to enhance functionality or fix bugs.

## License
This project is licensed under the MIT License. See the LICENSE file for details.

---
