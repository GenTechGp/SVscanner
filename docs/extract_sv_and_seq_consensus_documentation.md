
# Construction of `seq_consensus` for Structural Variant (SV) Types

This document summarizes the method used to construct the `seq_consensus` for each SV type in a VCF file. The consensus sequence is intended to represent the local genomic context around a variant and is output as a FASTA sequence.

---

## Common Parameters and Definitions

- **`f_len`**: Flanking length calculated as:

  `f_len = min(args.flen, svlen * args.ffac)`

  `[args.flen - "The maximum detectable period size supported by TRF to determine the length of flanking sequences", optional, default=2000]`
  
  `[args.ffac - "Multiplication factor for SVLEN to determine the length of flanking sequences", optional, default=10]`

- **`seq_fl`**: Sequence to the left of the SV (left flank).
- **`seq_fr`**: Sequence to the right of the SV (right flank).
- **`seq`**: Reference sequence or ALT sequence involved in the SV.

---

## SV Types and `seq_consensus` Construction

### DEL (Deletion)

**Description**: A segment of the reference genome is deleted.

**Construction**:  
`seq_consensus = left_flank + deleted_sequence_from_reference + right_flank`

**Code Behavior**:  
- The deleted segment is included in the consensus.
- Flanking regions are fetched and concatenated with the deleted region:

  `seq = fasta.fetch(region=f"{chrom}:{start_fl}-{end_fr}")`  
  `seq_consensus = seq`

---

### INS (Insertion)

**Description**: A new sequence is inserted at a given position.

**Construction**:  
`seq_consensus = left_flank + inserted_sequence_from_ALT + right_flank`

**Code Behavior**:  
- ALT field of VCF provides the inserted sequence.
- The sequence is not present in the reference.
- Inserted into context using:

  `seq_consensus = seq_fl + seq + seq_fr`

---

### DUP (Duplication)

**Description**: A segment of the genome is duplicated.

**Construction**:  
`seq_consensus = left_flank + duplicated_sequence + duplicated_sequence + right_flank`

**Code Behavior**:  
- Assumes tandem duplication (immediate repeat).
- Duplication inserted between flanks:

  `seq_consensus = seq_fl + seq + seq + seq_fr`

---

### INV (Inversion)

**Description**: A segment of the genome is inverted (reverse-complemented).

**Construction**:  
`seq_consensus = left_flank + reverse_complement(variant_sequence) + right_flank`

**Code Behavior**:  
- Inversion represented by reverse complement of the sequence:

  `seq_rc = Seq(seq).reverse_complement()`  
  `seq_consensus = seq_fl + seq_rc + seq_fr`

---

### BND (Breakend)/ TRA (Translocation)

**Description**: A complex breakpoint (translocation, etc.).

**Construction**:  
`seq_consensus = flanking_sequence_around_breakpoint`

**Code Behavior**:  
- Two flanking sequences are extracted (no insertion or rearranged sequence) for each genomic position described in the vcf record.

  `seq1 = fasta.fetch(region=f"{chrom1}:{start_f}-{end_f}")`  
  `seq2 = fasta.fetch(region=f"{chrom2}:{start_f}-{end_f}")`  
  `seq_consensus1 = seq1`
  `seq_consensus2 = seq2`

---

## Notes

- For INS, the actual length is determined from `len(ALT[0])` rather than from `SVLEN` due to known inaccuracies in some callers.
- For DUP, INV, and DEL, sequence integrity is verified against expected `svlen`.
- `seq_consensus` is written to a `.fa` FASTA file named by SV type and index (`<SVTYPE>.<i>.fa`).
