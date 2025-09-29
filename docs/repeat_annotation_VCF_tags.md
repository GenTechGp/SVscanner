# Repeat Annotation VCF tags

| CALLER_ID | Description |
|-----------|-------------|
| RM_CLASSIFICATION | Repeat class(es) overlapping the SV: SINE, LINE, LTR, DNA, Retroposon, NON-REPETITIVE. |
| RM_SV_COVERAGE | For each overlap, fraction of SV length covered (overlap_len / sv_len). |
| RM_ELEMENTS_COVERAGE | For each overlap, fraction of the repeat element covered by the SV (overlap_len / element_len). |
| RM_ELEMENT_PROPORTION | For each Repeat Masker hit, proportion of the query sequence (SV + flanks) aligning to that element. |
| RM_RECIPROCAL | Overlap class for the dominant RepeatMasker class. Full if at least ~75% of the SV is covered and each overlapped repeat element is mostly covered by the SV (≥75% of the element’s length). Partial if ≥75% of the SV is covered but one or more overlapped elements are not mostly covered. Minimal if <75% of the SV is covered. |
| RM_TOTAL_SV_COVERAGE | Proportion of the SV covered by RepeatMasker hits, with overlaps merged so regions aren’t counted twice [0–1]. |
| TRF_CLASSIFICATION | TRF repeat class(es) overlapping the SV: HOMO (homopolymer), STR, VNTR, TR, or NON-REPETITIVE. |
| TRF_SV_COVERAGE | For each TRF hit, fraction of SV length covered. |
| TRF_PERIOD_SIZE | Period size(s) in bp for each TRF hit. |
| TRF_COPY_NUMBER | Copy number(s) for each TRF hit. |
| TRF_TOTAL_SV_COVERAGE | Proportion of the SV covered by TRF hits, with overlaps merged so regions aren’t counted twice [0–1]. |
| CONSENSUS_REPEAT | Dominant repeat motif from TRF across the SV. |
| FINAL_CLASSIFICATION | Final repeat class for the SV, combining TRF and RepeatMasker. If only one method reports a repeat, use that class; if both do, choose the stronger/cleaner call as defined in the methods; if neither, NON-REPETITIVE/NA. RECIPROCAL is reported for mobile-element classes and NA for TRF-only classes. |
| DISEASE_GENE* | STR-associated disease gene (from STRchive). |
| STRCHIVE_MOTIF* | Pathogenic motif(s) matched in STRchive (normalised for rotation/complement). |
| PATHOGENIC_MIN* | Minimum pathogenic repeat count for the matched gene/motif (STRchive). |

*Added as INFO attribute only when SV intersects with position of gene