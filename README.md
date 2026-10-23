# Codon Table Generator

Educational Python program displaying all 64 standard codon assignments and mapping a DNA/RNA input in frame 1.

## Run

```bash
python codon_table_generator.py
```

Example input: `ATGGCTTTTTAA`.
Expected mapping: Methionine - Alanine - Phenylalanine - Stop.

## Input policy

- Coding-strand DNA or RNA written 5′ to 3′.
- Lowercase and whitespace normalized; one input line per interactive run.
- T and U together rejected; unsupported characters, including N, rejected.
- A/C/G-only input labeled DNA by convention; alphabet alone is inconclusive.
- Only ATG flagged as a potential start; not an exhaustive initiation model.
- All complete codons mapped even after Stop; not an ORF finder.
- One or two terminal bases reported without translation.
- Positions are one-based codon indices, not nucleotide coordinates.

## Tests

```bash
python -m unittest -v
```

The bundled verification recorded 29 passing unit tests, covering 26 function-level input cases, plus 24 successful single-line interactive checks. Two multiline strings were function-level tests only. Retest your edited local version.

## Optional reference comparison

```bash
python -m pip install biopython
python compare_codon_table.py
```

The checker writes comparison_results.csv. Bundled standard assignments matched 64/64 codons using Biopython 1.88. The generator itself does not need Biopython.

## GitHub submission

Open this assignment folder in VS Code. Preserve existing repository history. Review changes, stage intended files, commit with an accurate message, then push or publish to GitHub. Do not upload a virtual environment, credentials, or unrelated work. No GitHub repository or commit history was created by generating this package.

## Reference links

Standard assignments: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi

Comparison: https://github.com/biopython/biopython/blob/master/Bio/Data/CodonTable.py

Structural comparison: https://github.com/biotite-dev/biotite/blob/v1.8.0/src/biotite/sequence/codon.py

See COMPARISON.md and the complete report for policy distinctions, all ten assignment steps, comments versus commits, and the submission checklist.
