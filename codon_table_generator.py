"""Display the standard codon table and annotate a DNA/RNA sequence.

Frame 1 only; flag ATG; map all codons after stop labels.
The input must be a coding-strand DNA or RNA sequence written 5' to 3'.
"""

# 1. Keep all 64 standard DNA codons in a readable lookup dictionary.
CODON_TABLE = {
    "TTT": "Phenylalanine", "TTC": "Phenylalanine",
    "TTA": "Leucine", "TTG": "Leucine",
    "TCT": "Serine", "TCC": "Serine", "TCA": "Serine", "TCG": "Serine",
    "TAT": "Tyrosine", "TAC": "Tyrosine", "TAA": "Stop", "TAG": "Stop",
    "TGT": "Cysteine", "TGC": "Cysteine", "TGA": "Stop", "TGG": "Tryptophan",
    "CTT": "Leucine", "CTC": "Leucine", "CTA": "Leucine", "CTG": "Leucine",
    "CCT": "Proline", "CCC": "Proline", "CCA": "Proline", "CCG": "Proline",
    "CAT": "Histidine", "CAC": "Histidine", "CAA": "Glutamine", "CAG": "Glutamine",
    "CGT": "Arginine", "CGC": "Arginine", "CGA": "Arginine", "CGG": "Arginine",
    "ATT": "Isoleucine", "ATC": "Isoleucine", "ATA": "Isoleucine", "ATG": "Methionine",
    "ACT": "Threonine", "ACC": "Threonine", "ACA": "Threonine", "ACG": "Threonine",
    "AAT": "Asparagine", "AAC": "Asparagine", "AAA": "Lysine", "AAG": "Lysine",
    "AGT": "Serine", "AGC": "Serine", "AGA": "Arginine", "AGG": "Arginine",
    "GTT": "Valine", "GTC": "Valine", "GTA": "Valine", "GTG": "Valine",
    "GCT": "Alanine", "GCC": "Alanine", "GCA": "Alanine", "GCG": "Alanine",
    "GAT": "Aspartic acid", "GAC": "Aspartic acid",
    "GAA": "Glutamic acid", "GAG": "Glutamic acid",
    "GGT": "Glycine", "GGC": "Glycine", "GGA": "Glycine", "GGG": "Glycine",
}
START_CODONS = {"ATG"}
STOP_CODONS = {"TAA", "TAG", "TGA"}


def clean_and_validate(raw_sequence):
    """Return (cleaned input, DNA spelling, input label), or raise ValueError."""
    # 2. Ignore whitespace and letter case before checking the alphabet.
    sequence = "".join(raw_sequence.split()).upper()
    # 3. Reject empty input instead of producing an empty result silently.
    if not sequence:
        raise ValueError("The sequence cannot be empty.")
    # 4. Require either DNA or RNA notation, never a mixture of T and U.
    if "T" in sequence and "U" in sequence:
        raise ValueError("The sequence contains both T and U; enter DNA or RNA, not both.")
    # 5. Reject unsupported characters before any codon dictionary lookup.
    invalid = sorted(set(sequence) - set("ATGCU"))
    if invalid:
        raise ValueError("Invalid nucleotide(s): " + ", ".join(invalid))
    sequence_type = "RNA" if "U" in sequence else "DNA"
    # 6. Normalize RNA to DNA spelling so a single dictionary supports both.
    dna_sequence = sequence.replace("U", "T")
    return sequence, dna_sequence, sequence_type


def translate_sequence(dna_sequence):
    """Map complete codons in a validated DNA-spelling string; return remainder."""
    # 7. Exclude the terminal one or two bases from complete-codon extraction.
    complete_length = len(dna_sequence) - len(dna_sequence) % 3
    # 8. Read non-overlapping triplets starting at the first nucleotide.
    codons = [dna_sequence[i:i + 3] for i in range(0, complete_length, 3)]
    # 9. Map every complete codon, including those after a Stop label.
    amino_acids = [CODON_TABLE[codon] for codon in codons]
    return codons, amino_acids, dna_sequence[complete_length:]


def find_positions(codons):
    """Return one-based codon positions for starts and (position, codon) stops."""
    # 10. Use one-based codon positions rather than nucleotide coordinates.
    starts = [position for position, codon in enumerate(codons, 1)
              if codon in START_CODONS]
    stops = [(position, codon) for position, codon in enumerate(codons, 1)
             if codon in STOP_CODONS]
    return starts, stops


def display_codon_table():
    """Print the entire reference dictionary, not just the input's codons."""
    # 11. Sort codon keys to keep the reference display deterministic.
    print("STANDARD DNA CODON TABLE (64 CODONS)")
    print(f"{'Codon':<8}Amino acid")
    print("-" * 32)
    for codon in sorted(CODON_TABLE):
        print(f"{codon:<8}{CODON_TABLE[codon]}")


def display_translation(original, sequence_type, codons, amino_acids, remainder):
    """Print codon mapping, start/stop annotations, and any incomplete fragment."""
    # 12. Restore RNA spelling for output while retaining DNA keys internally.
    def shown(codon):
        return codon.replace("T", "U") if sequence_type == "RNA" else codon

    print("\nSEQUENCE ANALYSIS")
    print("Input type:", sequence_type)
    print("Cleaned sequence:", original)
    print(f"{'Position':<10}{'Codon':<8}Amino acid")
    print("-" * 48)
    for position, (codon, amino_acid) in enumerate(zip(codons, amino_acids), 1):
        label = amino_acid + " (Start)" if codon in START_CODONS else amino_acid
        print(f"{position:<10}{shown(codon):<8}{label}")
    if not codons:
        print("No complete codon is available.")
    print("Amino acid sequence:", " - ".join(amino_acids) or "None")
    starts, stops = find_positions(codons)
    print("Start codon position(s):", starts if starts else "None")
    stop_info = [(position, shown(codon)) for position, codon in stops]
    print("Stop codon(s):", stop_info if stop_info else "None")
    # 13. Report incomplete bases without assigning them an amino acid.
    if remainder:
        print(f"Incomplete codon ignored: {shown(remainder)} "
              f"({len(remainder)} nucleotide(s))")
    else:
        print("Incomplete codon: None")


def main():
    """Run one interactive analysis; return a nonzero status for invalid input."""
    display_codon_table()
    raw_sequence = input("\nEnter a DNA or RNA sequence: ")
    # 14. Catch validation errors and display an understandable message.
    try:
        original, dna_sequence, sequence_type = clean_and_validate(raw_sequence)
        codons, amino_acids, remainder = translate_sequence(dna_sequence)
        display_translation(original, sequence_type, codons, amino_acids, remainder)
    except ValueError as error:
        print("Error:", error)
        return 1
    return 0


# 15. Prevent interactive input when another module imports the dictionary.
if __name__ == "__main__":
    raise SystemExit(main())
