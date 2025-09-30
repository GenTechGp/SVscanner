# Repeat Annotation VCF tags

| CALLER_ID              | Description                                                                                                                                                                                                 | Type         |
|-------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------|
| RM_CLASSIFICATION       | Repeat class(es) overlapping the SV: SINE, LINE, LTR, DNA, Retroposon, NON-REPETITIVE.                                                                                                                     | list[string] |
| RM_SV_COVERAGE          | For each overlap, fraction of SV length covered (overlap_len / sv_len).                                                                                                                                     | list[float]  |
| RM_ELEMENTS_COVERAGE    | For each overlap, fraction of the repeat element covered by the SV (overlap_len / element_len).                                                                                                             | list[float]  |
| RM_ELEMENT_PROPORTION   | For each Repeat Masker hit, proportion of the query sequence (SV + flanks) aligning to that element.                                                                                                        | list[float]  |
| RM_RECIPROCAL           | Overlap class for the dominant RepeatMasker class. Full if at least ~75% of the SV is covered and each overlapped repeat element is mostly covered (≥75%). Partial if ≥75% SV coverage but some <75%.       | string       |
| RM_TOTAL_SV_COVERAGE    | Proportion of the SV covered by RepeatMasker hits, with overlaps merged so regions aren’t counted twice [0–1].                                                                                              | float        |
| TRF_CLASSIFICATION      | TRF repeat class(es) overlapping the SV: HOMO (homopolymer), STR, VNTR, TR, or NON-REPETITIVE.                                                                                                              | list[string] |
| TRF_SV_COVERAGE         | For each TRF hit, fraction of SV length covered.                                                                                                                                                            | list[float]  |
| TRF_PERIOD_SIZE         | Period size(s) in bp for each TRF hit.                                                                                                                                                                      | list[integer]|
| TRF_COPY_NUMBER         | Copy number(s) for each TRF hit.                                                                                                                                                                            | list[float]  |
| TRF_TOTAL_SV_COVERAGE   | Proportion of the SV covered by TRF hits, with overlaps merged so regions aren’t counted twice [0–1].                                                                                                      | float        |
| CONSENSUS_REPEAT        | Dominant repeat motif from TRF across the SV.                                                                                                                                                               | string       |
| FINAL_CLASSIFICATION    | Final repeat class for the SV, combining TRF and RepeatMasker. If only one method reports a repeat, use that class; if both do, choose the stronger/cleaner call; if neither, NON-REPETITIVE/NA.             | string       |
| DISEASE_GENE*           | STR-associated disease gene (from STRchive).                                                                                                                                                                | string       |
| STRCHIVE_MOTIF*         | Pathogenic motif(s) matched in STRchive (normalised for rotation/complement).                                                                                                                               | string       |
| PATHOGENIC_MIN*         | Minimum pathogenic repeat count for the matched gene/motif (STRchive).                                                                                                                                      | integer      |

*Added as INFO attribute only when SV intersects with position of gene

# Examples Records

| CALLER_ID             | Record 1                | Record 2 | Record 3         | Record 4          | Note                                                                 |
|-----------------------|-------------------------|----------|------------------|------------------|----------------------------------------------------------------------|
| RM_CLASSIFICATION     | LINE,LINE,SINE          | LTR      | NON_REPETITIVE   | NA               | Always available (value or NA)                                       |
| RM_SV_COVERAGE        | 0.53,0.35,0.1           | 1        |                  |                  | Not available if RM_CLASSIFICATION is NON_REPETITIVE or NA            |
| RM_ELEMENTS_COVERAGE  | 0.09,0.06,0.34          | 0.07     |                  |                  | Not available if RM_CLASSIFICATION is NON_REPETITIVE or NA            |
| RM_ELEMENT_PROPORTION | 0.09,0.06,0.98          | 0.15     |                  |                  | Not available if RM_CLASSIFICATION is NON_REPETITIVE or NA            |
| RM_RECIPROCAL         | Partial                 | Partial  | NA               | NA               | Always available (value or NA)                                       |
| RM_TOTAL_SV_COVERAGE  | 0.99                    | 1        | 0.07             |                  | Not available if RM_CLASSIFICATION is NA                             |
| TRF_CLASSIFICATION    | NON_REPETITIVE          | NA       | VNTR             | VNTR             | Always available (value or NA)                                       |
| TRF_SV_COVERAGE       |                         |          | 0.77             | 1                | Not available if TRF_CLASSIFICATION is NON_REPETITIVE or NA           |
| TRF_PERIOD_SIZE       |                         |          | 74               | 66               | Not available if TRF_CLASSIFICATION is NON_REPETITIVE or NA           |
| TRF_COPY_NUMBER       |                         |          | 8.8              | 44               | Not available if TRF_CLASSIFICATION is NON_REPETITIVE or NA           |
| TRF_TOTAL_SV_COVERAGE | 0                       |          | 0.77             | 1                | Not available if TRF_CLASSIFICATION is NA                             |
| CONSENSUS_REPEAT      | NA                      | NA       | ATAGGTGTTGGC     | GGAACAGTCGAGTG   | Always available (value or NA)                                       |
| FINAL_CLASSIFICATION  | LINE                    | LTR      | VNTR             | VNTR             | Always available                                                     |
