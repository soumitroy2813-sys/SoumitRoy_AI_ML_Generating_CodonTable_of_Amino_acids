"""Run using python -m unittest -v; no third-party dependency is required."""
import contextlib
import io
import unittest
from itertools import product
from codon_table_generator import (
    CODON_TABLE, clean_and_validate, translate_sequence, find_positions,
    display_translation,
)
from compare_codon_table import NAME_TO_SYMBOL

# Constructed software-test inputs, not samples claimed to be natural genes.
CASES = [
    ("original", "ATGGCTTTTTAA", "MAF*", "", "DNA"),
    ("tag_stop", "ATGAAAGGGTAG", "MKG*", "", "DNA"),
    ("tga_stop", "ATGTGGTGTTGA", "MWC*", "", "DNA"),
    ("no_start", "GCTTTTGGG", "AFG", "", "DNA"),
    ("two_starts", "ATGATGTAA", "MM*", "", "DNA"),
    ("one_leftover", "ATGGCTTAAC", "MA*", "C", "DNA"),
    ("two_leftovers", "ATGGCTTAACC", "MA*", "CC", "DNA"),
    ("lowercase", "atgaaagggtag", "MKG*", "", "DNA"),
    ("spaces", "ATG AAA GGG TAG", "MKG*", "", "DNA"),
    ("rna", "AUGAAAGGGUAG", "MKG*", "", "RNA"),
    ("rna_lowercase", "augucuuga", "MS*", "", "RNA"),
    ("rna_whitespace", "AUG\tGCU\nUAA", "MA*", "", "RNA"),
    ("one_base", "A", "", "A", "DNA"),
    ("two_bases", "AT", "", "AT", "DNA"),
    ("only_stops", "TAATAGTGA", "***", "", "DNA"),
    ("no_stop", "ATGAAA", "MK", "", "DNA"),
    ("after_stop", "ATGTAAGCT", "M*A", "", "DNA"),
    ("rna_one_leftover", "AUGU", "M", "U", "RNA"),
    ("rna_two_leftovers", "AUGUC", "M", "UC", "RNA"),
    ("acg_convention", "ACG", "T", "", "DNA"),
    ("empty", "", "ERROR", "The sequence cannot be empty.", None),
    ("whitespace_only", " \t\n", "ERROR", "The sequence cannot be empty.", None),
    ("invalid_b", "ATGBCT", "ERROR", "Invalid nucleotide(s): B", None),
    ("ambiguous_n", "ATGNAA", "ERROR", "Invalid nucleotide(s): N", None),
    ("numbers", "ATG123", "ERROR", "Invalid nucleotide(s): 1, 2, 3", None),
    ("mixed_tu", "ATUG", "ERROR", "both T and U", None),
]


class CodonTests(unittest.TestCase):
    def test_dictionary_complete(self):
        expected = {"".join(parts) for parts in product("ATGC", repeat=3)}
        self.assertEqual(set(CODON_TABLE), expected)
        self.assertEqual(sum(value == "Stop" for value in CODON_TABLE.values()), 3)
        self.assertEqual(len(set(CODON_TABLE.values()) - {"Stop"}), 20)

    def test_positions_are_one_based_codons(self):
        self.assertEqual(find_positions(["ATG", "ATG", "TAA"]),
                         ([1, 2], [(3, "TAA")]))

    def test_rna_output_and_remainder(self):
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            display_translation("AUGU", "RNA", ["ATG"], ["Methionine"], "T")
        self.assertIn("AUG", output.getvalue())
        self.assertIn("Incomplete codon ignored: U (1 nucleotide(s))", output.getvalue())


def make_test(case):
    def test(self):
        _, raw, expected, remainder, kind = case
        if expected == "ERROR":
            with self.assertRaises(ValueError) as caught:
                clean_and_validate(raw)
            self.assertIn(remainder, str(caught.exception))
            return
        original, dna, observed_kind = clean_and_validate(raw)
        codons, labels, rem = translate_sequence(dna)
        self.assertEqual("".join(NAME_TO_SYMBOL[label] for label in labels), expected)
        self.assertEqual(rem.replace("T", "U") if kind == "RNA" else rem, remainder)
        self.assertEqual(observed_kind, kind)
        self.assertEqual("".join(codons) + rem, dna)
    return test


for number, case in enumerate(CASES, 1):
    setattr(CodonTests, f"test_case_{number:02d}_{case[0]}", make_test(case))


if __name__ == "__main__":
    unittest.main()
