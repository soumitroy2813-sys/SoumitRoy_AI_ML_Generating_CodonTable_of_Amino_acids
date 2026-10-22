"""Run with: python -m unittest -v"""
import itertools
import json
from pathlib import Path
import unittest
from codon_table_generator import (
    CODON_TABLE, START_CODON, STOP_CODONS,
    clean_and_validate, translate_sequence,
)

SYMBOLS = {
    "Phenylalanine": "F", "Leucine": "L", "Serine": "S",
    "Tyrosine": "Y", "Stop": "*", "Cysteine": "C", "Tryptophan": "W",
    "Proline": "P", "Histidine": "H", "Glutamine": "Q", "Arginine": "R",
    "Isoleucine": "I", "Methionine": "M", "Threonine": "T",
    "Asparagine": "N", "Lysine": "K", "Valine": "V", "Alanine": "A",
    "Aspartic acid": "D", "Glutamic acid": "E", "Glycine": "G",
}


class CodonTableTests(unittest.TestCase):
    def test_complete_dictionary(self):
        expected = {"".join(c) for c in itertools.product("ATGC", repeat=3)}
        self.assertEqual(set(CODON_TABLE), expected)
        self.assertEqual({c for c, a in CODON_TABLE.items() if a == "Stop"}, STOP_CODONS)
        self.assertEqual(CODON_TABLE[START_CODON], "Methionine")

    def test_dna_rna_equivalence(self):
        dna = clean_and_validate("ATGGCTTTTTAA")[1]
        rna = clean_and_validate("AUGGCUUUUUAA")[1]
        self.assertEqual(translate_sequence(dna), translate_sequence(rna))

    def test_start_stop_positions(self):
        codons, _, _ = translate_sequence("ATGATGTAATAGTGA")
        self.assertEqual([i + 1 for i, c in enumerate(codons) if c == START_CODON], [1, 2])
        self.assertEqual([i + 1 for i, c in enumerate(codons) if c in STOP_CODONS], [3, 4, 5])


def build_case_test(case):
    def test(self):
        _, _, raw, kind, expected_codons, expected_amino, expected_remainder = case
        if kind == "error":
            with self.assertRaises(ValueError) as context:
                clean_and_validate(raw)
            self.assertEqual(str(context.exception), expected_codons)
        else:
            _, dna, actual_kind = clean_and_validate(raw)
            codons, amino, remainder = translate_sequence(dna)
            self.assertEqual(actual_kind, kind)
            self.assertEqual(" ".join(codons), expected_codons)
            self.assertEqual("".join(SYMBOLS[a] for a in amino), expected_amino)
            self.assertEqual(remainder, expected_remainder)
    return test


cases_path = Path(__file__).with_name("test_cases.json")
for case in json.loads(cases_path.read_text()):
    setattr(CodonTableTests, "test_" + case[0].lower(), build_case_test(case))

if __name__ == "__main__":
    unittest.main()
