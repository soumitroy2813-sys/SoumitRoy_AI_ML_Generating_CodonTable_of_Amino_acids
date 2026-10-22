# Development of a Codon Table Generator Using Python

**Academic project report — 20 October 2026**

## Abstract

A Python command-line application was developed to display the standard codon table and map DNA or RNA input to codon-level amino-acid labels. Its explicit dictionary contains 64 codons. Modular functions clean and validate input, normalize RNA spelling, extract complete triplets, translate codons, report start/stop occurrences, and identify incomplete terminal fragments. The bundled implementation passed 26 input cases and three additional structural or behavioral tests under Python 3.12.8. Every assignment matched NCBI standard-code data and Biotite's standard-table representation. These results concern the supplied implementation, not an unobserved local copy; GitHub publication has not been performed.

## Aim and requirements

**Aim:** Implement an explainable codon-table generator and nucleotide-sequence translator using fundamental Python programming constructs.

The assignment requires a 64-entry dictionary, DNA input validation, triplet splitting, amino-acid translation, start/stop identification, incomplete-codon handling, tabular output, comparison with GitHub implementations, and repository submission; RNA support is an additional challenge.[^1]

| Requirement | Implementation |
|---|---|
| Complete codon dictionary | `CODON_TABLE`, with 64 DNA keys |
| Sequence input | One interactive input per run |
| Validation | Empty-input, invalid-character, and mixed-T/U checks |
| Codon extraction | Non-overlapping triplets beginning at nucleotide 1 |
| Translation | Dictionary lookup for every complete codon |
| Start/stop reporting | ATG-only marking and all three standard stop labels |
| Incomplete codons | Report the terminal one- or two-base fragment |
| Display | Full reference table and input-specific translation |
| RNA extension | U-to-T normalization internally; RNA spelling in output |
| GitHub comparison | Biotite and Python Codon Tables |
| GitHub publication | Instructions supplied; upload requires the user |

## Biological background

The standard genetic code assigns 61 codons to 20 amino acids and three to termination. The DNA-alphabet representation uses T instead of U; ATG maps to methionine, while TAA, TAG, and TGA are termination codons. The program stores termination as the label `Stop`, not as another amino acid.[^2][^3]

Alternative genetic codes and initiation rules exist. NCBI table 1 records alternative initiation possibilities, including TTG and CTG, whereas this educational implementation flags only ATG. Its initiation policy is therefore narrower than the full biological rules and should be described separately from its standard amino-acid assignments.[^2]

DNA input is assumed to represent the coding strand written 5′ to 3′. The program does not reverse-complement a template strand, remove introns, discover genes, or establish genuine initiation sites. It processes only reading frame 1 of the supplied orientation.

## Design and methodology

### Dictionary representation

Each DNA triplet is a unique dictionary key, and the value is a full amino-acid name or `Stop`. For example, `CODON_TABLE["ATG"]` gives `"Methionine"`. Dictionaries support retrieval by keys, while sets support membership and alphabet checks.[^4]

The full table is displayed alphabetically before the input prompt. The program does not infer the genetic code from the entered sequence; it uses an explicit reference mapping.

### Input processing

The sequence is cleaned using:

```python
sequence = "".join(raw_sequence.split()).upper()
```

Whitespace is removed and lowercase letters become uppercase. The relevant normalization operations are documented Python string methods.[^5]

Validation proceeds in this order:

1. Reject empty input after cleaning.
2. Reject simultaneous T and U.
3. Reject characters outside A/T/G/C/U.
4. Label input containing U as RNA; otherwise label it DNA.
5. Replace U with T internally for dictionary lookup.

Ambiguous nucleotide symbols such as N and B are deliberately unsupported. A sequence containing only A, C, and G is labelled DNA by convention, because these letters alone cannot establish whether the input represents DNA or RNA.

### Triplets and incomplete input

```python
complete_length = len(dna_sequence) - (len(dna_sequence) % 3)
codons = [dna_sequence[i:i + 3]
          for i in range(0, complete_length, 3)]
remainder = dna_sequence[complete_length:]
amino_acids = [CODON_TABLE[codon] for codon in codons]
```

Only complete triplets undergo lookup. Any final fragment containing one or two bases is reported separately, preventing an incomplete key from reaching the dictionary.

### Start and stop policy

Every in-frame ATG receives a start flag. The flag does not prove that translation begins there in a cell. Start/stop positions are 1-based **codon numbers**, not nucleotide coordinates.

Stop codons are annotated without terminating the mapping process. Thus, `ATGTAAGCT` yields `Methionine - Stop - Alanine`. This is a codon-annotation result, not a prediction that a single ribosome continues through ordinary termination.

### Function organization

| Function | Purpose |
|---|---|
| `clean_and_validate()` | Clean input, validate alphabet, determine label, normalize RNA |
| `translate_sequence()` | Split validated DNA-spelling input and translate complete codons |
| `display_codon_table()` | Print all 64 reference codons |
| `display_translation()` | Print codon rows, translation, positions, and remainder |
| `main()` | Coordinate interaction and validation-error display |

`translate_sequence()` is a helper expecting validated input. Validation occurs before it is called in the main application.

## Verification

All 64 dictionary assignments were checked against NCBI table 1 and Biotite's standard-table data after converting full amino-acid names to one-letter symbols. All assignments matched. Biotite's representation is NCBI-derived, so this is a consistency check rather than two independent biological determinations.



Twenty-six deliberately constructed teaching/test inputs were executed. They cover ordinary DNA, all three stop codons, absent starts, repeated ATGs, lowercase input, spaces, RNA, short sequences, remainders, continuation after stops, empty input, invalid letters, numbers, and mixed T/U notation. These strings are software test cases, not claimed naturally occurring genes.



A reproducible standard-library unit-test suite additionally checks dictionary completeness, DNA/RNA equivalence, and start/stop positions. Together, the suite contains **29 tests, all of which passed**. The supplied project package contains the execution log and verification records; these results do not establish that a locally modified copy has passed.

## Sample output

Input:

```text
ATGGCTTTTTAA
```

Observed codon-specific output:

| Position | Codon | Label |
|---:|---|---|
| 1 | ATG | Methionine (Start) |
| 2 | GCT | Alanine |
| 3 | TTT | Phenylalanine |
| 4 | TAA | Stop |

```text
Amino acid sequence: Methionine - Alanine - Phenylalanine - Stop
Start codon position(s): [^1]
Stop codon(s): [(4, 'TAA')]
Incomplete codon: None
```

The complete console output includes the 64-row reference table before these results.



## Test results

The following selected records illustrate the executed behavior. Here `*` represents a stop label and single letters represent amino-acid assignments.[^2]

| Input | Translation or response | Remainder |
|---|---|---|
| `ATGAAAGGGTAG` | `MKG*` | None |
| `ATGTGGTGTTGA` | `MWC*` | None |
| `ATGAAAGGGTGGTGTTAA` | `MKGWC*` | None |
| `GCTTTTGGG` | `AFG`; no start or stop required | None |
| `ATGATGTAA` | `MM*`; ATG flags at positions 1 and 2 | None |
| `ATGGCTTAAC` | `MA*` | C |
| `ATGGCTTAACC` | `MA*` | CC |
| `atgaaagggtag` | `MKG*` | None |
| `ATG AAA GGG TAG` | `MKG*` | None |
| `AUGAAAGGGUAG` | `MKG*`; RNA spelling retained in display | None |
| `AT` | No complete codon | AT |
| `ATGTAAGCT` | `M*A`; mapping continues after stop | None |
| `AUGU` | `M` | U |
| `ATGBCT` | Invalid nucleotide B | Not translated |
| `ATGNAA` | Unsupported ambiguous nucleotide N | Not translated |
| `ATG123` | Invalid characters 1, 2, and 3 | Not translated |
| `ATUG` | Mixed T/U input rejected | Not translated |
| Empty input | Empty-sequence error | Not translated |

The full report and CSV include all 26 cases and pass/fail outcomes.



## GitHub comparison

Biotite provides a reusable `CodonTable` class with codon mappings, configurable starts, and NCBI table loading. Python Codon Tables provides organism-specific codon-usage data rather than a simple amino-acid translation dictionary.[^6][^7][^8]

| Dimension | Project implementation | Biotite | Python Codon Tables |
|---|---|---|---|
| Main purpose | Educational sequence annotation | Codon mapping in a sequence library.[^8] | Codon-usage tables.[^6] |
| Representation | Codon → full amino-acid name | Dictionary export uses symbols or numerical codes.[^8] | Amino-acid symbol → codon → relative frequency.[^7] |
| Standard ATG lookup | `Methionine` | `M`.[^8] | A lookup returns usage frequency, not a translation label.[^7] |
| Code-table flexibility | Standard assignments only | NCBI tables load by name or ID.[^8] | Multiple organism-specific usage tables.[^6] |
| Start marking | ATG only | Default table also marks only ATG; configurable.[^8] | Not the role of a usage-frequency table.[^6] |
| Translation handling | All complete frame-1 codons, remainder reported | Translation belongs to `NucleotideSequence`; complete mode requires a multiple-of-three length.[^9] | Not the principal operation of the usage API.[^7] |

The 64 project assignments agree with Biotite's standard-table data after representation normalization. Translation identity and usage frequency answer different questions, so frequency values must not replace amino-acid labels in the student's dictionary.[^7]

## Execution instructions

Run from the extracted project folder:

```bash
python3 codon_table_generator.py
```

Enter the nucleotide sequence when prompted. Run the command again for each new sequence. To verify the supplied code locally:

```bash
python3 -m unittest -v
```

For the previously supplied Mac paths:

```bash
/Users/soumitroy/SR_AI_ML/Python_AI_ML_venv/bin/python /Users/soumitroy/SR_AI_ML/codon_table_generator.py
```

The demonstrated earlier failure was caused by a bare `codon_table_generator` identifier on line 1. A filename does not define that Python variable; remove the stray line or make it a comment. Python raises `NameError` when an unqualified name is not found.[^10]

## GitHub submission

Create an empty repository on GitHub, then stage and commit the project files locally. GitHub documents configuring a remote and pushing the branch to upload locally hosted code.[^11]

```bash
git init
git add .
git commit -m "Add codon generator, tests and documentation"
git branch -M main
```

Copy the actual remote URL from GitHub and use it after `git remote add origin` on the same command line. Then run:

```bash
git remote -v
git push -u origin main
```

If an instructor has supplied an existing repository, use its clone/fork and submission instructions instead. No authenticated GitHub upload was carried out for this report; the resulting repository link must be supplied by the user.

## Limitations and extensions

- Only reading frame 1 of the supplied orientation is processed.
- DNA input is assumed to be the coding strand, not the template strand.
- ATG marking is an educational motif rule, not start-site prediction.
- Stop labels do not terminate codon mapping.
- Alternative codes, ambiguity handling, intron removal, and recoding are absent.
- FASTA files, batch processing, and repeated-input menus are not implemented.
- A/C/G-only input receives a conventional DNA label.
- Test results concern the bundled implementation, not unobserved local edits.

Future extensions could add explicit input-type selection, alternate reading frames, reverse-complement processing, genetic-code selection, separate ORF and codon-annotation modes, FASTA support, and CSV/JSON export.

For an input of length n, the design performs input-length-dependent cleaning, validation, splitting, and mapping passes. With its fixed 64-entry table, ordinary dictionary lookup gives linear overall scaling and linear stored-sequence/output space. This is an algorithmic assessment, not a measured performance benchmark.

## Deliverables

The downloadable report is provided in PDF and editable Word formats with the complete source-code appendix and full test matrix. The project archive contains the runnable Python script, 29-test suite, README, ignore rules, sample output, input-test CSV, codon-verification CSV, and execution records. GitHub publication is documented as a remaining user action rather than reported as completed.

---

## References

1. [Python_Project.pdf](https://ppl-ai-file-upload.s3.amazonaws.com/web/direct-files/attachments/82685109/069e93de-e4db-490b-9898-f836fabf3cfe/Python_Project.pdf?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Checksum-Mode=ENABLED&X-Amz-Credential=ASIA2F3EMEYE7R3JLGGV%2F20261022%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20261022T045747Z&X-Amz-Expires=2029&X-Amz-Security-Token=IQoJb3JpZ2luX2VjEKX%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJHMEUCIBDC67QIudDAcYFoW1c8C%2B9I8fUv%2BKSAP9pGX%2FpYGZevAiEA2q5J6F9t18P8Ei7dUVwwYWZIHkn7iIoWlXxsLnotAl8qgAUIbRABGgw2OTk3NTMzMDk3MDUiDFY9KElBt2GZygactSrdBOILGDlq7kyZ01Wxx1B4sVIca7d3ACM8H8dcCu40Cz%2FmN3ALU37s6tI%2FmQVJgGH8LeRkwJ7TfsddkgSxEpDoxGRXJfo60Yz88BJPTR7QI6EupmRMnLfbQYcBmnmJWzCnp1RBt0fCBd%2F7kfpXKrl5ooFaDfY84orayzcuiSUkf3z6PSEMEkYGT%2F9CgkjLS13iSCpAdO87H47wGOnAL19e8mAs6uSMvVGIF9ASD82I2%2BodkyV%2BGyi5da%2BWBG7sTVFbWpgErLvoYU4J7j9tGdJ4dwWf%2BtsExN9i%2FZfVs61tLFc%2F8gzHgeHvepXjVH1gZHTtVpOd01nBlP1VqkXbMBVntrRDkK9BwuAu3I9ka8BM1Za3iMIkwnbEIO105O7k2g5KMZxRi%2FaDasMgA8du7JJjbXPYaqDhXKFCBv7Qwh5pihIAyi1ZDZ2KUhhac3csOcC4Q9s9xi7S%2B%2Bqjmsxh%2BwCcGsvRcJ%2BxRP2RPVNN8Fnwiv5q3Gg3jl4%2FgJq2byWH2Dczbw75qWke%2FjLb94gNKPu0bDil%2FxbWejfCYZNv6iFuyQBIT3vU7jV5jdzcE6pNOzurul8LFoosWSKN0OTKFB8qSBN2YEuTgmcc2%2BUpsba4n78wYkUGxom%2BdpEhyP%2BW8Kg1KUMxXjUvBktZsnPvvCP9BMQeQdPpL6E0J0MyxuD3PTDrkDMMy%2FB9Kw38Dou5V1BeaBiy1GYlk%2BdtKVdiZib%2F%2FP50wAyJfrVkR9lsszhSHD7FA914A45QZPgkRS1KZOI91%2BV0sNl3BgCq8a9jZnUj%2B5nURSVutdFWDPT3w7zcMOWj5tYGOpgBvg8UNPnYY%2Btq8nAZI%2Fg%2FarbOWgDwyBGi7N5LU5KijVfaZYNHHvpS8LEfGD%2FYc9Dzqn02tfOns06MpanLW0AKj09wNLX7xCdvmdQzaZ0WycEiqPAGaJ3xRkQtaslQVLAPHdC2uy9Rgb9Q8RrI5YLuqpnWvx6uOGrwRQeMk40ksG0EontLQS%2B7ZHvQJ9Cuttw22mmp63ugw6s%3D&X-Amz-SignedHeaders=host&x-id=GetObject&X-Amz-Signature=d742506b07494138b54c8551d5a4dec9ee9cd0d235529a7726a5743cfca8069e) - page-1 Python Project Development of a Codon Table GeneratorObjective Develop a Python program that ...

2. [genetic code (translation table) - NCBI Taxonomy Homepage - NIH](https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes) - Thus, for the convenience of people reading GenBank records, the genetic code tables shown here use ...

3. [Genetic Codes - NCBI Taxonomy Homepage - NIH](https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes) - ## 1. The Standard Code (transl_table=1) By default all transl_table in GenBank flatfiles are equal ...

4. [5. Data Structures — Python 3.15.0 documentation](https://docs.python.org/3/tutorial/datastructures.html) - 5.5. Dictionaries¶ · It is best to think of a dictionary as a set of key: value pairs, with the requ...

5. [Built-in Types — Python 3.12.15 documentation](https://docs.python.org/3.12/library/stdtypes.html) - Return a copy of the string with all the cased characters [4] converted to uppercase.

6. [Edinburgh-Genome-Foundry/python_codon_tables - GitHub](https://github.com/Edinburgh-Genome-Foundry/python_codon_tables) - Provides codon usage tables as dictionaries, for Python.
...
All tables are from kazusa.or.jp (codon...

7. [python_codon_tables/python_codon_tables/python_codon_tables.py at master · Edinburgh-Genome-Foundry/python_codon_tables](https://github.com/Edinburgh-Genome-Foundry/python_codon_tables/blob/master/python_codon_tables/python_codon_tables.py) - Returns a dict {"*": {'TAA': 0.64...}, 'K': {'AAA': 0.76...}, ...}

8. [CodonTable — Biotite](https://www.biotite-python.org/latest/apidoc/biotite.sequence.CodonTable.html) - The table is equal to the NCBI “Standard” codon table, with the difference that only “ATG” is a star...

9. [NucleotideSequence - biotite.sequence](https://www.biotite-python.org/latest/apidoc/biotite.sequence.NucleotideSequence.html) - codon_tableCodonTable, optional. The codon table to be used. By default the default table will be us...

10. [Built-in Exceptions](https://docs.python.org/3/library/exceptions.html?highlight=attributeerror) - Raised when a local or global name is not found.
This applies only to unqualified names.

11. [docs.github.com · enterprise-cloud@latestAdding locally hosted code to GitHub](https://docs.github.com/en/enterprise-cloud@latest/migrations/importing-source-code/using-the-command-line-to-import-source-code/adding-locally-hosted-code-to-github) - 7. To push the changes in your local repository to GitHub, run the following command.
...
git push -...

