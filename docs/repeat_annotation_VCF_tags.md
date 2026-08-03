# Repeat Annotation VCF tags

| CALLER_ID               | Description                                                                                                                                                                                                 | Type         |
|-------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------|
| RM_CLASSIFICATION       | Repeat class(es) overlapping the SV: SINE, LINE, LTR, DNA, Retroposon, NON-REPETITIVE.                                                                                                                      | list[string] |
| RM_FAMILY               | RepeatMasker family for each hit, in the same order as RM_CLASSIFICATION (e.g. L1, Alu, ERVL-MaLR). '.' if no family sub-classification exists for that element.                                            | list[string] |
| RM_SUBFAMILY            | RepeatMasker subfamily (repeat model name) for each hit, in the same order as RM_CLASSIFICATION (e.g. L1M4, AluYa5, MLT1D).                                                                                 | list[string] |
| RM_SV_COVERAGE          | For each overlap, fraction of SV length covered (overlap_len / sv_len).                                                                                                                                     | list[float]  |
| RM_ELEMENT_PROPORTION   | For each RepeatMasker hit, fraction of the element's consensus model covered by the alignment (consensus_aligned / consensus_length). Always in [0, 1]; values near 1 indicate a complete element insertion. | list[float]  |
| RM_RECIPROCAL           | Completeness of the dominant RepeatMasker class. Full if the dominant class covers ≥75% of the SV and every element is ≥75% of its consensus model; Partial if ≥75% SV coverage but at least one element is <75% of its consensus; NA otherwise. Values: {Full, Partial, NA}. | string       |
| RM_TOTAL_SV_COVERAGE    | Proportion of the SV covered by RepeatMasker hits, with overlaps merged so regions aren’t counted twice [0–1].                                                                                              | float        |
| TRF_CLASSIFICATION      | TRF repeat class(es) overlapping the SV: HOMO (homopolymer), STR, VNTR, TR, or NON-REPETITIVE.                                                                                                              | list[string] |
| TRF_SV_COVERAGE         | For each TRF hit, fraction of SV length covered.                                                                                                                                                            | list[float]  |
| TRF_PERIOD_SIZE         | Period size(s) in bp for each TRF hit.                                                                                                                                                                      | list[integer]|
| TRF_COPY_NUMBER         | Copy number(s) for each TRF hit.                                                                                                                                                                            | list[float]  |
| TRF_TOTAL_SV_COVERAGE   | Proportion of the SV covered by TRF hits, with overlaps merged so regions aren’t counted twice [0–1].                                                                                                       | float        |
| CONSENSUS_REPEAT        | Dominant repeat motif from TRF across the SV.                                                                                                                                                               | list[string] |
| FINAL_CLASSIFICATION    | Final repeat class for the SV, combining TRF and RepeatMasker. If only one method reports a repeat, use that class; if both do, choose the stronger/cleaner call; if neither, NON-REPETITIVE/NA.            | string       |
| DISEASE_GENE*           | STR-associated disease gene (from STRchive).                                                                                                                                                                | string       |
| STRCHIVE_MOTIF*         | Pathogenic motif(s) matched in STRchive (normalised for rotation/complement).                                                                                                                               | string       |
| PATHOGENIC_MIN*         | Minimum pathogenic repeat count for the matched gene/motif (STRchive).                                                                                                                                      | integer      |
| BND_MATE_INFO           | Annotations for the second breakend (mate) query sequence in BND/TRA records. Pipe-delimited key=value pairs mirroring primary tags                                                                         | string      |

*Added as INFO attribute only when SV intersects with position of gene

# Examples Records

| CALLER_ID             | Record 1                | Record 2 | Record 3         | Record 4          | Note                                                                 |
|-----------------------|-------------------------|----------|------------------|------------------|----------------------------------------------------------------------|
| RM_CLASSIFICATION     | LINE,LINE,SINE          | LTR      | NON_REPETITIVE   | NA               | Always available (value or NA)                                       |
| RM_FAMILY             | L1,L1,Alu               | ERVL-MaLR|                  |                  | Not available if RM_CLASSIFICATION is NON_REPETITIVE or NA; '.' per element if no family exists |
| RM_SUBFAMILY          | L1M4,L1M4b,AluYa5       | MLT1D    |                  |                  | Not available if RM_CLASSIFICATION is NON_REPETITIVE or NA            |
| RM_SV_COVERAGE        | 0.53,0.35,0.1           | 1        |                  |                  | Not available if RM_CLASSIFICATION is NON_REPETITIVE or NA            |
| RM_ELEMENT_PROPORTION | 0.09,0.06,0.98          | 0.15     |                  |                  | Not available if RM_CLASSIFICATION is NON_REPETITIVE or NA            |
| RM_RECIPROCAL         | Partial                 | Partial  | NA               | NA               | Always available (value or NA)                                       |
| RM_TOTAL_SV_COVERAGE  | 0.99                    | 1        | 0.07             |                  | Not available if RM_CLASSIFICATION is NA                             |
| TRF_CLASSIFICATION    | NON_REPETITIVE          | NA       | VNTR             | VNTR             | Always available (value or NA)                                       |
| TRF_SV_COVERAGE       |                         |          | 0.77             | 1                | Not available if TRF_CLASSIFICATION is NON_REPETITIVE or NA           |
| TRF_PERIOD_SIZE       |                         |          | 74               | 66               | Not available if TRF_CLASSIFICATION is NON_REPETITIVE or NA           |
| TRF_COPY_NUMBER       |                         |          | 8.8              | 44               | Not available if TRF_CLASSIFICATION is NON_REPETITIVE or NA           |
| TRF_TOTAL_SV_COVERAGE | 0                       |          | 0.77             | 1                | Not available if TRF_CLASSIFICATION is NA                             |
| CONSENSUS_REPEAT      | NA                      | NA       | ATAGGTGTTGGC     | GGAACAGTCGAGTG   | Always available (value or NA)                                       |
| FINAL_CLASSIFICATION  | Repetitive/Mobile/LINE  | Repetitive/Mobile/LTR    | Repetitive/Tandem/VNTR             | Repetitive/Tandem/VNTR             | Always available                                                     |


### Example record with BND_MATE_INFO tag
```
chr21	5216247	Sniffles2.BND.A1CS14	G	G]chr20:30346078]	59	GT	PRECISE;SVTYPE=BND;SUPPORT=8;COVERAGE=0,0,55,62,67;STRAND=+-;CHR2=chr20;PHASE=2,5216322,7,8,FAIL,PASS;STDEV_POS=0.5;VAF=0.138;RM_CLASSIFICATION=NON_REPETITIVE;RM_RECIPROCAL=NA;RM_TOTAL_SV_COVERAGE=0.47;TRF_CLASSIFICATION=NA;CONSENSUS_REPEAT=NA;FINAL_CLASSIFICATION=NON_REPETITIVE;BND_MATE_INFO=RM_CLASSIFICATION=NON_REPETITIVE|RM_RECIPROCAL=NA|RM_TOTAL_SV_COVERAGE=0.5|TRF_CLASSIFICATION=NON_REPETITIVE|TRF_TOTAL_SV_COVERAGE=0.0|CONSENSUS_REPEAT=NA|FINAL_CLASSIFICATION=Repetitive/Mobile	GT:GQ:DR:DV:PS	0/0:59:50:8:5216322
```


# FINAL_CLASSIFICATION

`FINAL_CLASSSIFICATION` can be one of the following values:
1. Repetitive/Tandem
2. Repetitive/Tandem/HOMO
3. Repetitive/Tandem/STR
4. Repetitive/Tandem/VNTR
5. Repetitive/Tandem/TR
6. Repetitive/Mobile
7. Repetitive/Mobile/SINE
8. Repetitive/Mobile/LINE
9. Repetitive/Mobile/LTR
10. Repetitive/Mobile/DNA
11. Repetitive/Mobile/Retroposon
12. Repetitive/Mixed
13. Repetitive/Mixed/HOMO
14. Repetitive/Mixed/STR
15. Repetitive/Mixed/VNTR
16. Repetitive/Mixed/TR
17. Repetitive/Mixed/SINE
18. Repetitive/Mixed/LINE
19. Repetitive/Mixed/LTR
20. Repetitive/Mixed/DNA
21. Repetitive/Mixed/Retroposon
22. NON_REPETITIVE

The process of determing the `FINAL_CLASSIFICATION` is explained in [repeat annotation process](repeat_annotation_steps.md)