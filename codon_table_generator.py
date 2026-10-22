"""Codon table generator and DNA/RNA translator."""

CODON_TABLE = {
    "TTT": "Phenylalanine", "TTC": "Phenylalanine", "TTA": "Leucine", "TTG": "Leucine",
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
    "GAT": "Aspartic acid", "GAC": "Aspartic acid", "GAA": "Glutamic acid", "GAG": "Glutamic acid",
    "GGT": "Glycine", "GGC": "Glycine", "GGA": "Glycine", "GGG": "Glycine"
}

START_CODON = "ATG"
STOP_CODONS = {"TAA", "TAG", "TGA"}


def clean_and_validate(raw_sequence):
    sequence = "".join(raw_sequence.split()).upper()
    if not sequence:
        raise ValueError("The sequence cannot be empty.")

    alphabet = set(sequence)
    if "T" in alphabet and "U" in alphabet:
        raise ValueError("Enter DNA or RNA, not a mixture containing both T and U.")
    if not alphabet <= set("ATGCU"):
        invalid = sorted(alphabet - set("ATGCU"))
        raise ValueError("Invalid nucleotide(s): " + ", ".join(invalid))

    sequence_type = "RNA" if "U" in alphabet else "DNA"
    dna_sequence = sequence.replace("U", "T")
    return sequence, dna_sequence, sequence_type


def translate_sequence(dna_sequence):
    complete_length = len(dna_sequence) - (len(dna_sequence) % 3)
    codons = [dna_sequence[i:i + 3] for i in range(0, complete_length, 3)]
    amino_acids = [CODON_TABLE[codon] for codon in codons]
    remainder = dna_sequence[complete_length:]
    return codons, amino_acids, remainder


def display_codon_table():
    print("\nSTANDARD DNA CODON TABLE")
    print("-" * 32)
    print(f"{'Codon':<8}{'Amino acid'}")
    print("-" * 32)
    for codon in sorted(CODON_TABLE):
        print(f"{codon:<8}{CODON_TABLE[codon]}")


def display_translation(original_sequence, sequence_type, codons, amino_acids, remainder):
    shown_codons = [c.replace("T", "U") for c in codons] if sequence_type == "RNA" else codons
    start_positions = [i + 1 for i, codon in enumerate(codons) if codon == START_CODON]
    stop_info = [(i + 1, shown_codons[i]) for i, codon in enumerate(codons) if codon in STOP_CODONS]

    print(f"\nInput type: {sequence_type}")
    print(f"Clean sequence: {original_sequence}")
    print("\nTRANSLATION")
    print("-" * 35)
    print(f"{'Position':<10}{'Codon':<10}{'Amino acid'}")
    print("-" * 35)
    for position, (codon, amino_acid) in enumerate(zip(shown_codons, amino_acids), start=1):
        label = amino_acid
        if codons[position - 1] == START_CODON:
            label += " (Start)"
        print(f"{position:<10}{codon:<10}{label}")

    print("\nAmino acid sequence:", " - ".join(amino_acids) if amino_acids else "No complete codon")
    print("Start codon position(s):", start_positions if start_positions else "None")
    print("Stop codon(s):", stop_info if stop_info else "None")
    if remainder:
        shown_remainder = remainder.replace("T", "U") if sequence_type == "RNA" else remainder
        print(f"Incomplete codon ignored: {shown_remainder} ({len(remainder)} nucleotide(s))")
    else:
        print("Incomplete codon: None")


def main():
    display_codon_table()
    raw_sequence = input("\nEnter a DNA or RNA sequence: ")
    try:
        original, dna_sequence, sequence_type = clean_and_validate(raw_sequence)
        codons, amino_acids, remainder = translate_sequence(dna_sequence)
        display_translation(original, sequence_type, codons, amino_acids, remainder)
    except ValueError as error:
        print("Error:", error)


if __name__ == "__main__":
    main()
