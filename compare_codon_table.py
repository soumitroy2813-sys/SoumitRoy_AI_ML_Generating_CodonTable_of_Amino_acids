"""Verify the local dictionary using Biopython's standard DNA table 1.

Optional dependency: python -m pip install biopython
Run: python compare_codon_table.py
"""
import csv
from itertools import product
from pathlib import Path
from codon_table_generator import CODON_TABLE

NAME_TO_SYMBOL = dict(zip(
    ["Alanine", "Arginine", "Asparagine", "Aspartic acid", "Cysteine",
     "Glutamine", "Glutamic acid", "Glycine", "Histidine", "Isoleucine",
     "Leucine", "Lysine", "Methionine", "Phenylalanine", "Proline", "Serine",
     "Threonine", "Tryptophan", "Tyrosine", "Valine", "Stop"],
    "ARNDCQEGHILKMFPSTWYV*"
))


def main():
    try:
        from Bio import __version__
        from Bio.Data import CodonTable
    except ImportError:
        print("Optional checker dependency: python -m pip install biopython")
        return 2
    reference = CodonTable.unambiguous_dna_by_id[1]
    expected_keys = {"".join(bases) for bases in product("ATGC", repeat=3)}
    missing = sorted(expected_keys - set(CODON_TABLE))
    extra = sorted(repr(k) for k in set(CODON_TABLE) - expected_keys)
    rows = []
    for codon in sorted(expected_keys):
        label = CODON_TABLE.get(codon)
        symbol = NAME_TO_SYMBOL.get(label)
        # Stops are separate from Biopython's 61 sense-codon forward mappings.
        expected = "*" if codon in reference.stop_codons else reference.forward_table[codon]
        rows.append([codon, label, symbol, expected,
                     "PASS" if symbol == expected and codon in CODON_TABLE else "FAIL"])
    target = Path(__file__).resolve().with_name("comparison_results.csv")
    with target.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["codon", "local_label", "local_symbol", "reference_symbol", "result"])
        writer.writerows(rows)
    matches = sum(row[-1] == "PASS" for row in rows)
    passed = matches == 64 and not missing and not extra
    print("Biopython version:", __version__)
    print("Reference: standard DNA table 1")
    print("Dictionary entries:", len(CODON_TABLE))
    print("Missing codons:", missing)
    print("Unexpected keys:", extra)
    print(f"Matching standard assignments: {matches}/64")
    print("Reference starts:", reference.start_codons)
    print("Start-site policy is separate from amino-acid mapping.")
    print("Overall:", "PASS" if passed else "FAIL")
    print("Saved:", target.name)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
