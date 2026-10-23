# Codon Table Generator — Complete Code Explanation



## 1. Program map

| Lines | Block | Job |
|---|---|---|
| 1–5 | Module docstring | States purpose and input assumptions |
| 7–29 | `CODON_TABLE`, start/stop sets | Reference data |
| 32–49 | `clean_and_validate()` | Clean, validate, convert U to T |
| 52–60 | `translate_sequence()` | Split codons, map labels, keep remainder |
| 63–70 | `find_positions()` | One-based start/stop positions |
| 73–80 | `display_codon_table()` | Print all 64 codons |
| 83–109 | `display_translation()` | Print analysis |
| 112–124 | `main()` | Run the pipeline, catch errors |
| 127–129 | Entry guard | Run `main()` only when executed directly |

**Execution order:** definitions load first, then the guard calls `main()` → `display_codon_table()` → `input()` → `clean_and_validate()` → `translate_sequence()` → `display_translation()` (which calls `find_positions()`) → exit status.

## 2. Docstring, dictionary, sets (lines 1–29)

```python
  1  """Display the standard codon table and annotate a DNA/RNA sequence.
  2  
  3  Educational policy: frame 1 only; flag ATG; map all codons after stop labels.
  4  The input must be a coding-strand DNA or RNA sequence written 5' to 3'.
  5  """
  6  
  7  # 1. Keep all 64 standard DNA codons in a readable lookup dictionary.
  8  CODON_TABLE = {
```

- **Lines 1–5:** a triple-quoted string at the top documents the module. It states frame 1 only, ATG flagged, all codons mapped even after a stop, and coding-strand 5'→3' input.
- **Line 8–27:** `CODON_TABLE` is a dictionary of `"codon": "label"` pairs. Keys are unique; many codons share one amino acid. `"Stop"` is a label, not an amino acid.
- **Lines 28–29:** `START_CODONS` and `STOP_CODONS` are *sets* used for fast `codon in SET` tests. Only `ATG` is treated as a start.
- All keys use `T`; RNA is converted to this spelling before lookup.


### The 64 entries, grouped by label

| Label | Codons | Count |
|---|---|---|
| Alanine | GCT, GCC, GCA, GCG | 4 |
| Arginine | CGT, CGC, CGA, CGG, AGA, AGG | 6 |
| Asparagine | AAT, AAC | 2 |
| Aspartic acid | GAT, GAC | 2 |
| Cysteine | TGT, TGC | 2 |
| Glutamic acid | GAA, GAG | 2 |
| Glutamine | CAA, CAG | 2 |
| Glycine | GGT, GGC, GGA, GGG | 4 |
| Histidine | CAT, CAC | 2 |
| Isoleucine | ATT, ATC, ATA | 3 |
| Leucine | TTA, TTG, CTT, CTC, CTA, CTG | 6 |
| Lysine | AAA, AAG | 2 |
| Methionine | ATG | 1 |
| Phenylalanine | TTT, TTC | 2 |
| Proline | CCT, CCC, CCA, CCG | 4 |
| Serine | TCT, TCC, TCA, TCG, AGT, AGC | 6 |
| Stop | TAA, TAG, TGA | 3 |
| Threonine | ACT, ACC, ACA, ACG | 4 |
| Tryptophan | TGG | 1 |
| Tyrosine | TAT, TAC | 2 |
| Valine | GTT, GTC, GTA, GTG | 4 |

Total: 64 entries (61 amino-acid codons across 20 names, plus 3 stops). [file:211]

## 3. clean_and_validate() (lines 32–49)

```python
 32  def clean_and_validate(raw_sequence):
 33      """Return (cleaned input, DNA spelling, input label), or raise ValueError."""
 34      # 2. Ignore whitespace and letter case before checking the alphabet.
 35      sequence = "".join(raw_sequence.split()).upper()
 36      # 3. Reject empty input instead of producing an empty result silently.
 37      if not sequence:
 38          raise ValueError("The sequence cannot be empty.")
 39      # 4. Require either DNA or RNA notation, never a mixture of T and U.
 40      if "T" in sequence and "U" in sequence:
 41          raise ValueError("The sequence contains both T and U; enter DNA or RNA, not both.")
 42      # 5. Reject unsupported characters before any codon dictionary lookup.
 43      invalid = sorted(set(sequence) - set("ATGCU"))
 44      if invalid:
 45          raise ValueError("Invalid nucleotide(s): " + ", ".join(invalid))
 46      sequence_type = "RNA" if "U" in sequence else "DNA"
 47      # 6. Normalize RNA to DNA spelling so a single dictionary supports both.
 48      dna_sequence = sequence.replace("U", "T")
 49      return sequence, dna_sequence, sequence_type
```

| Lines | Explanation |
|---|---|
| 35 | `split()` breaks on whitespace, `"".join()` rejoins without spaces, `upper()` fixes case |
| 37–38 | Empty or whitespace-only input raises `ValueError` |
| 40–41 | Containing both `T` and `U` is rejected (checked before invalid characters) |
| 43–45 | `set(sequence) - set("ATGCU")` leaves illegal characters; `sorted` gives stable messages |
| 46 | `"RNA" if "U" in sequence else "DNA"`; DNA is the default label |
| 48 | `replace("U","T")` makes the lookup string |
| 49 | Returns a 3-value tuple: cleaned input, DNA spelling, type |

Example: `"  aug gcu uaa  "` returns `("AUGGCUUAA", "ATGGCTTAA", "RNA")`.


## 4. translate_sequence() (lines 52–60)

```python
 52  def translate_sequence(dna_sequence):
 53      """Map complete codons in a validated DNA-spelling string; return remainder."""
 54      # 7. Exclude the terminal one or two bases from complete-codon extraction.
 55      complete_length = len(dna_sequence) - len(dna_sequence) % 3
 56      # 8. Read non-overlapping triplets starting at the first nucleotide.
 57      codons = [dna_sequence[i:i + 3] for i in range(0, complete_length, 3)]
 58      # 9. Map every complete codon, including those after a Stop label.
 59      amino_acids = [CODON_TABLE[codon] for codon in codons]
 60      return codons, amino_acids, dna_sequence[complete_length:]
```

- **55:** `len % 3` is the leftover count; subtracting it gives the length of complete codons.
- **57:** `range(0, complete_length, 3)` yields offsets 0, 3, 6…; `dna_sequence[i:i+3]` slices three bases (end index excluded).
- **59:** each codon is looked up in `CODON_TABLE`. There is no break at a stop.
- **60:** returns codons, labels, and the 0–2 leftover bases.

| Length | `% 3` | `complete_length` | Leftover |
|---|---|---|---|
| 12 | 0 | 12 | 0 |
| 13 | 1 | 12 | 1 |
| 14 | 2 | 12 | 2 |
| 2 | 2 | 0 | 2 |


## 5. find_positions() (lines 63–70)

```python
 63  def find_positions(codons):
 64      """Return one-based codon positions for starts and (position, codon) stops."""
 65      # 10. Use one-based codon positions rather than nucleotide coordinates.
 66      starts = [position for position, codon in enumerate(codons, 1)
 67                if codon in START_CODONS]
 68      stops = [(position, codon) for position, codon in enumerate(codons, 1)
 69               if codon in STOP_CODONS]
 70      return starts, stops
```

- `enumerate(codons, 1)` pairs each codon with a one-based position.
- Two list comprehensions keep starts (positions only) and stops (`(position, codon)` tuples).
- For `["ATG","GCT","ATG","TAA"]` the result is `([1, 3], [(4, "TAA")])`.
- Codon position *p* begins at nucleotide 3(p − 1) + 1 (one-based). Every in-frame `ATG` is flagged, which does not prove a gene starts there.


## 6. display_codon_table() (lines 73–80)

```python
 73  def display_codon_table():
 74      """Print the entire reference dictionary, not just the input's codons."""
 75      # 11. Sort codon keys to keep the reference display deterministic.
 76      print("STANDARD DNA CODON TABLE (64 CODONS)")
 77      print(f"{'Codon':<8}Amino acid")
 78      print("-" * 32)
 79      for codon in sorted(CODON_TABLE):
 80          print(f"{codon:<8}{CODON_TABLE[codon]}")
```

- Line 77: an f-string; `{'Codon':<8}` left-aligns in 8 columns.
- Line 78: `"-" * 32` repeats the dash 32 times.
- Line 79: `sorted(CODON_TABLE)` sorts the keys without changing the dictionary.
- The reference always prints in DNA spelling.


## 7. display_translation() (lines 83–109)

```python
 83  def display_translation(original, sequence_type, codons, amino_acids, remainder):
 84      """Print codon mapping, start/stop annotations, and any incomplete fragment."""
 85      # 12. Restore RNA spelling for output while retaining DNA keys internally.
 86      def shown(codon):
 87          return codon.replace("T", "U") if sequence_type == "RNA" else codon
 88  
 89      print("\nSEQUENCE ANALYSIS")
 90      print("Input type:", sequence_type)
 91      print("Cleaned sequence:", original)
 92      print(f"{'Position':<10}{'Codon':<8}Amino acid")
 93      print("-" * 48)
 94      for position, (codon, amino_acid) in enumerate(zip(codons, amino_acids), 1):
 95          label = amino_acid + " (Start)" if codon in START_CODONS else amino_acid
 96          print(f"{position:<10}{shown(codon):<8}{label}")
 97      if not codons:
 98          print("No complete codon is available.")
 99      print("Amino acid sequence:", " - ".join(amino_acids) or "None")
100      starts, stops = find_positions(codons)
101      print("Start codon position(s):", starts if starts else "None")
102      stop_info = [(position, shown(codon)) for position, codon in stops]
103      print("Stop codon(s):", stop_info if stop_info else "None")
104      # 13. Report incomplete bases without assigning them an amino acid.
105      if remainder:
106          print(f"Incomplete codon ignored: {shown(remainder)} "
107                f"({len(remainder)} nucleotide(s))")
108      else:
109          print("Incomplete codon: None")
```

- **86–87:** nested `shown()` converts `T` to `U` for display when the input was RNA; internal lookup stays DNA.
- **94:** `zip(codons, amino_acids)` pairs each codon with its label, and `enumerate(..., 1)` adds the position; the target `position, (codon, amino_acid)` unpacks both levels.
- **95:** appends `" (Start)"` to the label when the codon is `ATG`.
- **97–98:** prints a message if there are no complete codons.
- **99:** `" - ".join(amino_acids) or "None"` falls back to `None` when the list is empty.
- **100–103:** calls `find_positions()` and prints starts and stops (stops shown in RNA spelling if needed).
- **105–109:** prints the incomplete fragment and its length, or `None`.


## 8. main() and entry point (lines 112–129)

```python
112  def main():
113      """Run one interactive analysis; return a nonzero status for invalid input."""
114      display_codon_table()
115      raw_sequence = input("\nEnter a DNA or RNA sequence: ")
116      # 14. Catch validation errors and display an understandable message.
117      try:
118          original, dna_sequence, sequence_type = clean_and_validate(raw_sequence)
119          codons, amino_acids, remainder = translate_sequence(dna_sequence)
120          display_translation(original, sequence_type, codons, amino_acids, remainder)
121      except ValueError as error:
122          print("Error:", error)
123          return 1
124      return 0
125  
126  
127  # 15. Prevent interactive input when another module imports the dictionary.
128  if __name__ == "__main__":
129      raise SystemExit(main())
```

- `main()` prints the table, reads one line with `input()`, then runs the pipeline inside `try`.
- `except ValueError as error` prints `Error: ...` and returns `1`; success returns `0`.
- `if __name__ == "__main__":` runs only when the file is executed directly, not when imported.
- `raise SystemExit(main())` exits with that status.
- Only `ValueError` is caught; `KeyboardInterrupt` and `EOFError` are not.


## 9. Worked examples (actual output of this file)

**Input:** `'ATGGCTTTTTAA'`  → exit status 0

```text
SEQUENCE ANALYSIS
Input type: DNA
Cleaned sequence: ATGGCTTTTTAA
Position  Codon   Amino acid
------------------------------------------------
1         ATG     Methionine (Start)
2         GCT     Alanine
3         TTT     Phenylalanine
4         TAA     Stop
Amino acid sequence: Methionine - Alanine - Phenylalanine - Stop
Start codon position(s): [1]
Stop codon(s): [(4, 'TAA')]
Incomplete codon: None
```

**Input:** `'AUGGCUUUUUAAC'`  → exit status 0

```text
SEQUENCE ANALYSIS
Input type: RNA
Cleaned sequence: AUGGCUUUUUAAC
Position  Codon   Amino acid
------------------------------------------------
1         AUG     Methionine (Start)
2         GCU     Alanine
3         UUU     Phenylalanine
4         UAA     Stop
Amino acid sequence: Methionine - Alanine - Phenylalanine - Stop
Start codon position(s): [1]
Stop codon(s): [(4, 'UAA')]
Incomplete codon ignored: C (1 nucleotide(s))
```

**Input:** `'ATGTAAGCT'`  → exit status 0

```text
SEQUENCE ANALYSIS
Input type: DNA
Cleaned sequence: ATGTAAGCT
Position  Codon   Amino acid
------------------------------------------------
1         ATG     Methionine (Start)
2         TAA     Stop
3         GCT     Alanine
Amino acid sequence: Methionine - Stop - Alanine
Start codon position(s): [1]
Stop codon(s): [(2, 'TAA')]
Incomplete codon: None
```

**Input:** `'ATGATGTAA'`  → exit status 0

```text
SEQUENCE ANALYSIS
Input type: DNA
Cleaned sequence: ATGATGTAA
Position  Codon   Amino acid
------------------------------------------------
1         ATG     Methionine (Start)
2         ATG     Methionine (Start)
3         TAA     Stop
Amino acid sequence: Methionine - Methionine - Stop
Start codon position(s): [1, 2]
Stop codon(s): [(3, 'TAA')]
Incomplete codon: None
```

**Input:** `'AT'`  → exit status 0

```text
SEQUENCE ANALYSIS
Input type: DNA
Cleaned sequence: AT
Position  Codon   Amino acid
------------------------------------------------
No complete codon is available.
Amino acid sequence: None
Start codon position(s): None
Stop codon(s): None
Incomplete codon ignored: AT (2 nucleotide(s))
```

**Input:** `'ATGBCT'`  → exit status 1

```text
Error: Invalid nucleotide(s): B
```

**Input:** `'ATUG'`  → exit status 1

```text
Error: The sequence contains both T and U; enter DNA or RNA, not both.
```

**Input:** `''`  → exit status 1

```text
Error: The sequence cannot be empty.
```

## 10. The 15 existing comments

| # | Line | Purpose |
|---|---|---|
| 1 | 7 | Dictionary of 64 codons |
| 2 | 34 | Ignore whitespace and case |
| 3 | 36 | Reject empty input |
| 4 | 39 | Reject mixed T/U |
| 5 | 42 | Reject unsupported characters |
| 6 | 47 | Normalize RNA to DNA spelling |
| 7 | 54 | Exclude terminal 1–2 bases |
| 8 | 56 | Non-overlapping triplets from base 1 |
| 9 | 58 | Map codons even after Stop |
| 10 | 65 | One-based codon positions |
| 11 | 75 | Sorted reference display |
| 12 | 85 | Restore RNA spelling for output |
| 13 | 104 | Report incomplete bases |
| 14 | 116 | Catch validation errors |
| 15 | 127 | Main guard |

Comments are separate from Git commits; saving code does not create a commit.

## 11. Syntax cheat sheet

| Syntax | Meaning here |
|---|---|
| `def f(x):` | Define a function |
| `return a, b, c` | Return a tuple |
| `s[i:i+3]` | Slice, end excluded |
| `% 3` | Remainder after division by 3 |
| `[x for ... if ...]` | List comprehension with filter |
| `enumerate(xs, 1)` | Counter starting at 1 |
| `zip(a, b)` | Pair elements by position |
| `a or "None"` | Fallback if `a` is empty |
| `f"{x:<8}"` | Left-align in 8 columns |
| `raise ValueError(...)` | Signal bad input |
| `try / except ... as e` | Catch an error |
| `raise SystemExit(n)` | Exit with status n |

## 12. Limits and viva answers

- Reads frame 1 only; no ORF search, reverse complement, or alternative genetic codes.
- Flags every in-frame `ATG` and maps codons after stops.
- Rejects ambiguity codes (N, R, Y), digits, and FASTA headers.
- Defaults to the label DNA when neither T nor U appears.

| Question | Answer |
|---|---|
| Why a dictionary? | Direct codon-to-label lookup |
| Why convert U to T? | One dictionary serves DNA and RNA |
| Why `len % 3`? | Separates complete codons from leftovers |
| Does Stop end the loop? | No; it is only a label |
| Why the main guard? | Allows safe import |
| What does `return 1` mean? | Handled validation failure |

## Appendix: full source

```python
  1  """Display the standard codon table and annotate a DNA/RNA sequence.
  2  
  3  Educational policy: frame 1 only; flag ATG; map all codons after stop labels.
  4  The input must be a coding-strand DNA or RNA sequence written 5' to 3'.
  5  """
  6  
  7  # 1. Keep all 64 standard DNA codons in a readable lookup dictionary.
  8  CODON_TABLE = {
  9      "TTT": "Phenylalanine", "TTC": "Phenylalanine",
 10      "TTA": "Leucine", "TTG": "Leucine",
 11      "TCT": "Serine", "TCC": "Serine", "TCA": "Serine", "TCG": "Serine",
 12      "TAT": "Tyrosine", "TAC": "Tyrosine", "TAA": "Stop", "TAG": "Stop",
 13      "TGT": "Cysteine", "TGC": "Cysteine", "TGA": "Stop", "TGG": "Tryptophan",
 14      "CTT": "Leucine", "CTC": "Leucine", "CTA": "Leucine", "CTG": "Leucine",
 15      "CCT": "Proline", "CCC": "Proline", "CCA": "Proline", "CCG": "Proline",
 16      "CAT": "Histidine", "CAC": "Histidine", "CAA": "Glutamine", "CAG": "Glutamine",
 17      "CGT": "Arginine", "CGC": "Arginine", "CGA": "Arginine", "CGG": "Arginine",
 18      "ATT": "Isoleucine", "ATC": "Isoleucine", "ATA": "Isoleucine", "ATG": "Methionine",
 19      "ACT": "Threonine", "ACC": "Threonine", "ACA": "Threonine", "ACG": "Threonine",
 20      "AAT": "Asparagine", "AAC": "Asparagine", "AAA": "Lysine", "AAG": "Lysine",
 21      "AGT": "Serine", "AGC": "Serine", "AGA": "Arginine", "AGG": "Arginine",
 22      "GTT": "Valine", "GTC": "Valine", "GTA": "Valine", "GTG": "Valine",
 23      "GCT": "Alanine", "GCC": "Alanine", "GCA": "Alanine", "GCG": "Alanine",
 24      "GAT": "Aspartic acid", "GAC": "Aspartic acid",
 25      "GAA": "Glutamic acid", "GAG": "Glutamic acid",
 26      "GGT": "Glycine", "GGC": "Glycine", "GGA": "Glycine", "GGG": "Glycine",
 27  }
 28  START_CODONS = {"ATG"}
 29  STOP_CODONS = {"TAA", "TAG", "TGA"}
 30  
 31  
 32  def clean_and_validate(raw_sequence):
 33      """Return (cleaned input, DNA spelling, input label), or raise ValueError."""
 34      # 2. Ignore whitespace and letter case before checking the alphabet.
 35      sequence = "".join(raw_sequence.split()).upper()
 36      # 3. Reject empty input instead of producing an empty result silently.
 37      if not sequence:
 38          raise ValueError("The sequence cannot be empty.")
 39      # 4. Require either DNA or RNA notation, never a mixture of T and U.
 40      if "T" in sequence and "U" in sequence:
 41          raise ValueError("The sequence contains both T and U; enter DNA or RNA, not both.")
 42      # 5. Reject unsupported characters before any codon dictionary lookup.
 43      invalid = sorted(set(sequence) - set("ATGCU"))
 44      if invalid:
 45          raise ValueError("Invalid nucleotide(s): " + ", ".join(invalid))
 46      sequence_type = "RNA" if "U" in sequence else "DNA"
 47      # 6. Normalize RNA to DNA spelling so a single dictionary supports both.
 48      dna_sequence = sequence.replace("U", "T")
 49      return sequence, dna_sequence, sequence_type
 50  
 51  
 52  def translate_sequence(dna_sequence):
 53      """Map complete codons in a validated DNA-spelling string; return remainder."""
 54      # 7. Exclude the terminal one or two bases from complete-codon extraction.
 55      complete_length = len(dna_sequence) - len(dna_sequence) % 3
 56      # 8. Read non-overlapping triplets starting at the first nucleotide.
 57      codons = [dna_sequence[i:i + 3] for i in range(0, complete_length, 3)]
 58      # 9. Map every complete codon, including those after a Stop label.
 59      amino_acids = [CODON_TABLE[codon] for codon in codons]
 60      return codons, amino_acids, dna_sequence[complete_length:]
 61  
 62  
 63  def find_positions(codons):
 64      """Return one-based codon positions for starts and (position, codon) stops."""
 65      # 10. Use one-based codon positions rather than nucleotide coordinates.
 66      starts = [position for position, codon in enumerate(codons, 1)
 67                if codon in START_CODONS]
 68      stops = [(position, codon) for position, codon in enumerate(codons, 1)
 69               if codon in STOP_CODONS]
 70      return starts, stops
 71  
 72  
 73  def display_codon_table():
 74      """Print the entire reference dictionary, not just the input's codons."""
 75      # 11. Sort codon keys to keep the reference display deterministic.
 76      print("STANDARD DNA CODON TABLE (64 CODONS)")
 77      print(f"{'Codon':<8}Amino acid")
 78      print("-" * 32)
 79      for codon in sorted(CODON_TABLE):
 80          print(f"{codon:<8}{CODON_TABLE[codon]}")
 81  
 82  
 83  def display_translation(original, sequence_type, codons, amino_acids, remainder):
 84      """Print codon mapping, start/stop annotations, and any incomplete fragment."""
 85      # 12. Restore RNA spelling for output while retaining DNA keys internally.
 86      def shown(codon):
 87          return codon.replace("T", "U") if sequence_type == "RNA" else codon
 88  
 89      print("\nSEQUENCE ANALYSIS")
 90      print("Input type:", sequence_type)
 91      print("Cleaned sequence:", original)
 92      print(f"{'Position':<10}{'Codon':<8}Amino acid")
 93      print("-" * 48)
 94      for position, (codon, amino_acid) in enumerate(zip(codons, amino_acids), 1):
 95          label = amino_acid + " (Start)" if codon in START_CODONS else amino_acid
 96          print(f"{position:<10}{shown(codon):<8}{label}")
 97      if not codons:
 98          print("No complete codon is available.")
 99      print("Amino acid sequence:", " - ".join(amino_acids) or "None")
100      starts, stops = find_positions(codons)
101      print("Start codon position(s):", starts if starts else "None")
102      stop_info = [(position, shown(codon)) for position, codon in stops]
103      print("Stop codon(s):", stop_info if stop_info else "None")
104      # 13. Report incomplete bases without assigning them an amino acid.
105      if remainder:
106          print(f"Incomplete codon ignored: {shown(remainder)} "
107                f"({len(remainder)} nucleotide(s))")
108      else:
109          print("Incomplete codon: None")
110  
111  
112  def main():
113      """Run one interactive analysis; return a nonzero status for invalid input."""
114      display_codon_table()
115      raw_sequence = input("\nEnter a DNA or RNA sequence: ")
116      # 14. Catch validation errors and display an understandable message.
117      try:
118          original, dna_sequence, sequence_type = clean_and_validate(raw_sequence)
119          codons, amino_acids, remainder = translate_sequence(dna_sequence)
120          display_translation(original, sequence_type, codons, amino_acids, remainder)
121      except ValueError as error:
122          print("Error:", error)
123          return 1
124      return 0
125  
126  
127  # 15. Prevent interactive input when another module imports the dictionary.
128  if __name__ == "__main__":
129      raise SystemExit(main())
```
