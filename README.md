# Global-Alignment
## Introduction
A python implementation of the Needleman–Wunsch algorithm. This program will find the optimal global alignments for two DNA, RNA or amino acid sequences.
## How to get the code
Type the following command into the terminal

`$ git clone https://github.com/Brandonag1997/Global-Alignment.git`
## Usage
The global alignment program takes 2 texts files containing either nucleotide or amino-acid sequences. To run navigate to ./Global-Alignment and make globalalign.py executable by running

`$ chmod +x globalalign.py`

The program takes at least 2 arguments, the locations of the sequences you want to align. For example

`$ ./globalalign.py ./examples/exampleSequence1.txt ./examples/exampleSequence2.txt`

This will output the best alignments for the 2 sequences and save the results to Output.txt

## Options
- `--similarity_matrix`
  - This optional parameter specifies the similarity matrix if no matrix is specied a simple match/mismatch/gap scoring matrix will be used
  - Available similarity matrices include
    - BLOSUM62 
    - BLOSUM30 increased by 5 up to BLOSUM90
    - BLUSOM100
    - PAM10 increased by 10 up to PAM500
    - NUC.4.2 and NUC.4.4
    - DAYHOFF
    - GONNET
    - IDENTITY
    - MATCH (gives the same results as the default)
- `--match`
  - The score that is assigned to 2 nucleotides or amino acids that match. This value defaults to 1. This is only used when a similarity matrix is not specified.
- `--mismatch`
  - The score that is assigned to 2 nucleotides or amino acids that do not match. This value defaults to -1. This is only used when a similarity matrix is not specified.
- `--gap`
  - The score that is assigned to a nucleotide or amino acid that is aligned to a gap. This value defaults to -1.
