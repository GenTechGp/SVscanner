# User Guide (Algorithms)

## SVChecker
To identify supporting reads for inversions and duplications, we locate alignments overlapping the left and right breakpoints of the SV and assess whether each read supported the SV. The two-letter code indicates whether the start or end of the read’s alignment overlaps a breakpoint
* LE (Left breakpoint, End of read alignment)  
* LS (Left breakpoint, Start of read alignment)  
* RE (Right breakpoint, End of read alignment)  
* RS (Right breakpoint, Start of read alignment)

### Inversions
A read is **valid** when the start and end of the alignments have opposite orientations at the left and right SV breakpoints. 
* BOTH: Alignments that span the entire inversion must intersect both the left and right breakpoints with the start and end of an alignment to be valid.   
* LEFT/RIGHT: Reads that span one SV breakpoint in the sample must have the start or end of alignments support the left and right SV breakpoints in the reference.
![GS](/images/inversion_passed_dark.svg#gh-dark-mode-only)
![GS](/images/inversion_passed.svg#gh-light-mode-only)


Additional FLAGS are implemented for edge cases 
* MULTI: multiple alignments intersect an SV breakpoints  
* EXTENSION: an alignment extended to or beyond the left or right breakpoint 
![GS](/images/inversion_flags_dark.svg#gh-dark-mode-only)
![GS](/images/inversion_flags.svg#gh-light-mode-only)
A read is rejected when it does not pass quality checks and/or one or both of the breakpoints are not supported.
![GS](/images/inversion_rejected_dark.svg#gh-dark-mode-only)
![GS](/images/inversion_rejected.svg#gh-light-mode-only)

### Duplications
A read is **valid** when the intersection between the start of an alignment with the left breakpoint (LS) and the end of another alignment with the right breakpoint (RE) can be paired. A valid pairing requires consecutive LS and RE alignments in the same orientation.  
![GS](/images/duplication_pair_dark.svg#gh-dark-mode-only)
![GS](/images/duplication_pair.svg#gh-light-mode-only)

The read order and number of alignments is assessed to classify alignments as tandem (consecutive) or interspersed (gaps) and duplication or repeat. 
* TANDEM\_DUPLICATION  
* TANDEM\_REPEAT  
* INTERSPERSED\_DUPLICATION  
* INTERSPERSED\_REPEAT

If any alignment intersects are unpaired, these flags are appended with the maximal number of pairs that should exist. 
![GS](/images/duplication_passed_dark.svg#gh-dark-mode-only)
![GS](/images/duplication_passed.svg#gh-light-mode-only)

### Extended Details

| Common Flags |  |
| :---- | :---- |
| NUM\_SPLIT | Rejected: The number of alignments exceeds the max number of splits  |

| Inversion Flags | Description  |
| :---- | :---- |
| PASSED\_BOTH | Passed: Start and end of alignments support the SV breakpoints (BOTH in count results) |
| PASSED\_RIGHT | Passed: End of alignments support the SV breakpoints (RIGHT\_BP in count results) |
| PASSED\_LEFT | Passed: Start of alignments support the SV breakpoints (LEFT\_BP in count results) |
| MUTLI\_\[LE,LS,RE,RS\] | Flag: Multiple alignments with the same intersect are present (MULTI in count results ) |
| PASSED\_EXTENSION\_\[LS,RE\] | Passed: An alignment extends beyond the start or end of another alignment (EXTENSION in count results) |
| MISSING\_BOTH | Rejected: Neither breakpoint is supported \[(LS,RE), (LE, RS)\] |
| MISSING\_LEFT\_BP | Rejected: The left breakpoint of the SV is not supported \[(LE, RE, RS), (RE,RS)\] |
| MISSING\_RIGHT\_BP | Rejected: The right breakpoint of the SV is not supported \[LE, LS, RS), (LE,LS)\] |
| ORIENTATION | Rejected: The orientation of the alignments do not support an inversion (ORIENTATION in count results) |
| PRIMARY\_CHR | Rejected: The primary chromosome is located on another chromosome OR has MAPQ \< 25 |

| Duplication Flags | Description  |
| :---- | :---- |
| TANDEM\_DUPLICATION | Passed: Read is made up of two consecutive alignments |
| TANDEM\_REPEAT | Passed: Read is made up of more than two consecutive alignments |
| INTERSPERSED\_DUPLICATION | Passed: Read contains two alignments with a break  |
| INTERSPERSED\_REPEAT | Passed: Read contains multiple alignments that are not consecutive |
| ___.PAIRS	 | Flag: Read contains unpaired intersects \-\> a count of the maximal number of pairs that should exist. |
| DUPLEX | Passed: The read is potentially a duplex \- contains two pairs with opposite orientation  |
| MISSING\_INTERSECT | Rejected: The read is missing alignments with intersect LS or RE |
| MISSING\_PAIRS | Rejected: No valid pairs created from the alignments  |
| MAPQ | Rejected: The read is made up of majority alignments of MAPQ is below \<25 |

## SVClassifier
For **Tandem Repeat Finder**, entries are determined by prioritising maximal intersection between SV and repeat and minimal period size in intervals (e.g. 0.05). Non-overlapping entries are selected in order of priority. Based on the period size of the repeat, each filtered entry is identified as HOMO (1bp), STR (2-12bp), TR (\>12bp).   
![](/images/TRF_workflow.svg) 
![GS](/images/TRF_workflow_dark.svg#gh-dark-mode-only)
![GS](/images/TRF_workflow.svg#gh-light-mode-only)

For **RepeatMasker**, entries are determined by prioritising those with maximal intersection between the SV and repeat entry. Entries are grouped into their repeat class (e.g. SINE, LINE) and in order of priority non-overlapping entries are selected within each class. Based on the element coverage and SV coverage, the SV is classified as a complete transposition when there is a reciprocal overlap of \>=75% or a fragment.
![GS](/images/RM_workflow_dark.svg.svg#gh-dark-mode-only)
![GS](/images/RM_workflow.svg.svg#gh-light-mode-only)
The ‘Repetitive’ classification is given to SVs where the intersecting repeat elements cover a threshold (e.g. 0.5). For SVs containing both tandem repeats and mobile elements, the final classification is determined by the type with the highest total coverage. 

### Annotations

| CALLER\_ID | Caller ID for the SV |
| :---- | :---- |
| RM\_CLASSIFICATION | (list) Classification of repeat class covering the SV \[SINE,LINE,LTR,DNA,Retroposon or NON-REPETITIVE\] |
| RM\_ELEMENTS\_COVERAGE | (list) Fraction of the mobile element covered by the SV |
| RM\_ELEMENT\_PROPORTION | (list) Proportion of the query sequence (includes flanking region) found in the mobile element |
| RM\_TRANSPOSITION | Type of transposition \[COMPLETE/FRAGMENT\] |
| RM\_SV\_COVERAGE | (list) Fraction of the SV covered by the mobile element intersection/SV\_length |
| RM\_TOTAL\_SV\_COVERAGE | (float) Total coverage of SV covered by mobile elements |
| TRF\_CLASSIFICATION | Classification(s) of tandem repeat class covering the SV \[HOMO,STR,TR or NON-REPETITIVE\] |
| TRF\_SV\_COVERAGE | Fraction of the SV covered by the tandem repeat |
| TRF\_PERIOD\_SIZE | Period size of the repeat(s) |
| TRF\_COPY\_NUMBER | Copy number of the repeat(s) |
| TRF\_TOTAL\_SV\_COVERAGE | Total coverage of SV covered by tandem repeats |
| CONSENSUS\_REPEAT | Motif of repeat(s) found by Tandem Repeat Finder |
| FINAL\_CLASSIFICATION | Classification of SV as repetitive element based on TRF and RepeatMasker results |
| DISEASE_GENE* | STR disease associated with gene (annotated by STRchive) |
| STRCHIVE_MOTIF* | Consensus Repeat is a version (rotation/complement) of pathogenic motif(s) annotated in STRchive |
| PATHOGENIC_MIN* | Minimum pathogenic number annotated in STRchive |

*Added as INFO attribute only when SV intersects with position of gene

![Illustration](/images/annotation_tag_illustration.png)

Note - `ri_l` is not necessarily equal to `ri_sv_l`. It can be less than, equal or greater than `ri_sv_l`.

* RM\_ELEMENTS\_COVERAGE = `ri_sv_l/ri` e.g. [`r1_sv_l/r1`, `r2_sv_l/r2`]
* RM\_ELEMENT\_PROPORTION = `ri_l/ri` e.g. [`r1_l/r1`, `r2_l/r2`]
* RM\_SV\_COVERAGE = `ri_sv_l/sv_l` e.g. [`r1_sv_l/sv_l`, `r2_sv_l/sv_l`]
* RM\_TOTAL\_SV\_COVERAGE = `sum(ri_sv_l)/sv_l` e.g. `(r1_sv_l+r2_sv_l)/sv_l`