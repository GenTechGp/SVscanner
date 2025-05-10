# Understanding RepeatMasker Output Formats

RepeatMasker (RM) is a tool for identifying and classifying repetitive DNA elements. It produces several output formats that differ in detail and intended use. Two of the most commonly used formats are `.out` and `.cat`.

## 1. `.out` Format and Explanation of Columns

The `.out` file is a summary report listing each repetitive element found in the input sequence. Here's an example record:

```
 140   14.6  3.5  1.9  BND.1      2002  2204  (1796) + L1M3         LINE/L1            3726   3931 (2209)    1
```

**Column Explanation**:

| Column | Description |
|--------|-------------|
| 1      | Smith-Waterman score of the match |
| 2–4    | % divergence, % deletions, % insertions from the consensus |
| 5      | Name of the query sequence |
| 6      | starting position of match in query sequence |
| 7      | ending position of match in query sequence |
| 8      | no. of bases in query sequence past the ending position of match |
| 9      | Strand (`+` or `C` for complement) |
| 10     | Name of the matching repeat |
| 11     | Repeat class/family |
| 12     | starting position of match in repeat sequence |
| 13     | ending position of match in repeat sequence |
| 14     | no. of bases in repeat ***consensus** sequence past the ending position of match |
| 15     | ID number of the repeat match (for internal tracking) |

Note - The above explanation is for `+`. For `C` strand the description vary specially for 8th and 14th columns.

## 2. Corresponding `.cat` Record

The `.cat` file contains more detailed alignment information. The corresponding `.cat` record for the same repeat match above might look like:

```
140 14.57 3.45 1.94 BND.1 2002 2204 (1796) L1M3_orf2#LINE/L1 1617 1822 (1466) m_b1s452i0
```

### Key Differences:

- The `.cat` file specifies the exact **subfamily** (`L1M3_orf2`) whereas `.out` uses the **collapsed family name** (`L1M3`).
- The `.cat` file shows correct alignment coordinates based on the actual consensus sequence (e.g., `1617–1822` for a 3.3k model), unlike `.out` which may show inflated virtual coordinates (e.g., `3726–3931` in a ~6k coordinate system).

## 3. Visualization and Explanation

Below is a screenshot of Dfam entries showing consensus lengths for different L1M3 subfamilies:

![Dfam L1M3 Entries](/images/dfam_records.png)

### Explanation

- None of the individual L1M3 subfamilies shown in Dfam have consensus sequences over 6,000 base pairs.
- The `.out` file uses a **composite coordinate space** for repeat families. That’s why `L1M3` can appear to extend past 6,000 bp — it’s not a bug, but an internal linear coordinate space combining fragments like `L1M3a_5end`, `L1M3_orf2`, etc.
- The `.cat` file and Dfam use the **true length** of the specific subfamily model.

According to the [RepeatMasker documentation](http://www.repeatmasker.org), these inconsistencies arise due to:

- **Fragment stitching:** RM may merge fragments from similar subfamilies or attempt to assign a unified consensus coordinate range even across fragmented elements.
- **LINE structure normalization:** For LINEs (like L1M3), consensus sequences are **not full-length** in the database. Instead, LINEs are represented in parts (5' UTRs, 3' UTRs, ORFs). The `.out` coordinates are mapped relative to a normalized full-length LINE (e.g., L1PA2), hence appearing longer than any actual fragment.
- **Annotation legibility and biological relevance:** The `.cat` file and alignment outputs reflect more biologically plausible matches and correct subfamily lengths.

### From RepeatMasker Docs:
> *"...all position numbers in older LINE1 subfamilies are adjusted to the position of ORF2 (the conserved part of LINE1) in a complete L1PA2 element...this sometimes results in the assignment of negative position numbers for the 5' end of LINEs..."*
