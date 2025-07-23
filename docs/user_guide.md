# User Guide (Algorithms)

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

* RM\_ELEMENTS\_COVERAGE = `ri_sv_l/ri_L` e.g. [`r1_sv_l/r1_L`, `r2_sv_l/r2_L`]
* RM\_ELEMENT\_PROPORTION = `ri_l/ri_L` e.g. [`r1_l/r1_L`, `r2_l/r2_L`]
* RM\_SV\_COVERAGE = `ri_sv_l/sv_L` e.g. [`r1_sv_l/sv_L`, `r2_sv_l/sv_L`]
* RM\_TOTAL\_SV\_COVERAGE = `sum(ri_sv_l)/sv_L` e.g. `(r1_sv_l+r2_sv_l)/sv_L`