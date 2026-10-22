# Codon Table Generator

An educational Python command-line tool that displays all 64 standard DNA
codons and translates DNA/RNA input in reading frame 1.

## Requirements

Python 3 with f-string support (3.6 or newer); verified with Python 3.12.8.
No third-party dependencies are required for running the program or tests.

## Run

```bash
python3 codon_table_generator.py
```

Enter `ATGGCTTTTTAA` when prompted.
Expected: Methionine - Alanine - Phenylalanine - Stop.
Run again for each new input.

## Features

- Full standard-code dictionary, displayed alphabetically
- Case-insensitive DNA/RNA input with whitespace removal
- Invalid-character, empty-input and mixed-T/U checks
- Triplet splitting and codon-by-codon translation
- ATG flags and all three standard stop-codon flags
- Incomplete terminal-codon reporting
- DNA/RNA output labels and RNA-spelling display

## Assumptions

Input DNA represents the coding strand written 5' to 3'. The program starts
at the first base, does not search for an ORF and does not stop at the first
stop. `Stop` is a label, not an amino acid. All in-frame ATGs are marked;
these are not experimentally established translation initiation sites.
ACG-only input is labelled DNA by convention because its alphabet alone
cannot distinguish DNA from RNA. Ambiguous IUPAC symbols and FASTA headers
are not accepted.

## Run tests

```bash
python3 -m unittest -v
```

The suite contains 26 input cases plus 3 structural/behavior tests.
`test_results.csv` records the 26 executed input cases;
`codon_verification.csv` records all 64 codon assignments.

## Reference and GitHub comparison

NCBI standard code:
https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter

Biotite CodonTable:
https://github.com/biotite-dev/biotite/blob/main/src/biotite/sequence/codon.py

Its dictionary export stores codon-to-single-letter mappings; internally,
Biotite uses encoded numerical arrays. Its default table uses the standard
assignments with ATG as the sole marked start. The 64 project mappings were
compared with NCBI table 1 and Biotite's standard-table data after converting
full amino-acid names to symbols; all mappings matched. These references
share an NCBI origin and are not independent biological determinations.

Python Codon Tables:
https://github.com/Edinburgh-Genome-Foundry/python_codon_tables

This repository stores nested amino-acid/codon/relative-frequency dictionaries.
It concerns codon usage and is not equivalent to a simple translation table.
No organism-specific frequency data are claimed for this project.

## Verification scope

The accompanying report and logs concern the bundled implementation tested
in an isolated Python environment. They do not establish that a modified
local copy has passed or that the project has been uploaded to GitHub.
