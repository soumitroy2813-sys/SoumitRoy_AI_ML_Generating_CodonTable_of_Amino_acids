# Codon-table comparison

## References inspected

- NCBI standard table 1: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
- Biopython source: https://github.com/biopython/biopython/blob/master/Bio/Data/CodonTable.py
- Biotite versioned source: https://github.com/biotite-dev/biotite/blob/v1.8.0/src/biotite/sequence/codon.py
- Biotite tutorial: https://www.biotite-python.org/latest/tutorial/sequence/types.html

## Representation

The local implementation has 64 DNA-spelling keys and full amino-acid names, including Stop. Biopython's standard table has 61 one-letter sense mappings and a separate list of three stops. The checker normalizes names and stop labels and verifies missing and unexpected keys.

## Recorded bundled result

- Biopython version 1.88; standard DNA table 1.
- 64/64 amino-acid/stop assignments matched.
- No missing or unexpected codon keys.
- Direct NCBI assignments also matched 64/64; evidence in ncbi_verification.csv.
- Biopython derives tables from NCBI, so agreement is not independent experimental validation.

## Start policy

The educational script flags ATG only. Biopython standard table 1 lists TTG, CTG, and ATG as possible initiators. Biotite's documented default uses ATG only. A start-site policy difference is not an amino-acid assignment mismatch. Internal ATG is methionine and is not by itself proof of initiation.

## Scope

Biopython comparison was executed. Biotite comparison was structural/documentation-based; no Biotite runtime benchmark was performed. Results apply to the bundled implementation, not automatically to an edited local file. Run `python compare_codon_table.py` before making a final local-verification claim.
