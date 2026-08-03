# Repeat/SV Diagram

This document describes the format used to represent structural variants (SVs) with annotated **TRF (Tandem Repeats)** and **RM (RepeatMasker)** elements, highlighting their overlap with the variant region and flanking sequences. The terms "repeat," "element," "repeat element," and "mobile element" are used interchangeably in this context.

Each record starts with a line beginning with a `>` symbol.
- **`SV_ID`**: A unique identifier for the structural variant.
- **`SOURCE_ID`**: Original SV ID from caller, e.g., Sniffles.
- **`GENOMIC_REGION`**: Region of the genome spanned by the SV including flanks (e.g., `chr1:8692-12927`).
- **`SV_REGION`**: The precise SV site or inserted sequence location (e.g., `chr1:10692-10927`).
- **`SV_LENGTH`**: The length of the SV in base pairs.
- **`FINAL_CLASSIFICATION`**: The final repeat classification for this SV (e.g. `Repetitive/Mobile/LINE`, `NON_REPETITIVE`).

Note that diagrams show only the parts of the repeat elements that align to the SV.

```

>INS.112	(Sniffles2.INS.5AS0)	chr1:34257940-34263488	chr1:34259940-34259940	1549	Repetitive/Mobile
 SV & flanking      ------------------------------------[##########################]------------------------------------ 
 -RM-
  Retroposon                                            [..................]                                             	67.79%[96.97%]	chr1:34259940-34259940	rm_362*
  SINE                                                                     [..]                                          	9.1%[42.44%]	chr1:34259940-34259940	rm_363*
  SINE                                                                        [...]                                      	16.66%[84.89%]	chr1:34259940-34259940	rm_364*
 -TRF-
  49                                                              [......]                                               	27.44%	chr1:34259940-34259940	trf_122*
  40                                                                     [..]                                            	10.91%	chr1:34259940-34259940	trf_123

 SV Diagram         [##################################################################################################] 
 -RM-
  Retroposon          [...................................................................]                              	SVA_A
  SINE                                                                                    [........]                     	AluY
  SINE                                                                                             [................]    	AluYk11
 -TRF-
  49                                                    [...........................]                                    	TCCTCACTTCCCAGACGGGATGGCGGCCGGGAAGGCCGGGCAGAGACGC


>INS.113	(Sniffles2.INS.5BS0)	chr1:34769095-34770501	chr1:34769765-34769765	67	Repetitive/Tandem/VNTR
 SV & flanking      -----------------------------------------------[####]----------------------------------------------- 
 -RM-
  Simple_repeat                                                     [...]                                                	91.04%[98.48%]	chr1:34769765-34769765	rm_4423
 -TRF-
  14                                                             [......]                                                	100.0%	chr1:34769737-34769765	trf_1484*
  3                                                                 [...]                                                	91.04%	chr1:34769765-34769765	trf_1485

 SV Diagram         [##################################################################################################] 
 -TRF-
  14                .................................................................................................... 	CTCCCCACCACCAC


>INS.114	(Sniffles2.INS.5CS0)	chr1:34989130-34995645	chr1:34991130-34991130	2516	Repetitive/Mobile/LINE
 SV & flanking      ------------------------------[######################################]------------------------------ 
 -RM-
  LINE                                            [.....................................]                                	96.74%[39.55%]	chr1:34991130-34991130	rm_5551*

 SV Diagram         [##################################################################################################] 
 -RM-
  LINE              [................................................................................................]   	L1HS


>INS.115	(Sniffles2.INS.5DS0)	chr1:35079307-35080860	chr1:35080047-35080047	74	Repetitive/Mixed/VNTR
 SV & flanking      -----------------------------------------------[####]----------------------------------------------- 
 -RM-
  LINE                                                   [.............]                                                 	85.14%[2.75%]	chr1:35079887-35080047	rm_4612*
  SINE                                                                 [.................]                               	13.51%[87.18%]	chr1:35080047-35080309	rm_4613*
 -TRF-
  11                                                             [......]                                                	100.0%	chr1:35080019-35080048	trf_1566*

 SV Diagram         [##################################################################################################] 
 -RM-
  LINE              .....................................................................................]               	L1MC3
  SINE                                                                                                    [............. 	AluSx
 -TRF-
  11                .................................................................................................... 	TATTAAATCAT



```
Immediately following the header is a **flanking diagram** that visualizes the SV relative to its surrounding sequence:
- `-`: Flanking genomic sequence  
- `[####]`: The SV itself (e.g., deletion, insertion, etc.)

Annotations from TRF and RM tools are listed below the flanking diagram. They indicate repetitive elements that **intersect with the SV** or flanking sequence.

- **`PERIOD`**: Repeat unit period (TRF) or class (RM like LINE, SINE)
- **`[REPEAT_ALIGNMNET_SHAPE]`**: ASCII diagram of the repeat’s alignment position relative to the flanking/SV diagram
- **`PERCENT`**:  
  - For TRF: `sv_coverage%`  
  - For RM: `sv_coverage%[element_proportion%]`, e.g., `84.11%[39.73%]`. `element_proportion` is the fraction of the element's consensus model covered by the alignment (always ≤ 100%).
- **`REPEAT_ALIGNMNET_GENOMIC_RANGE`**: Start–end coordinates of the repeat's alignment
- **`ID`**: Unique identifier (e.g., `rm_668`, `trf_216`)  
  - A trailing `*` means this repeat **was selected at the pre-classification stage**

A second visualization shows a zoomed-in, **SV-only** view (no flanks) containing only the elements selected at the pre-classification stage (those marked `*` in the flanking section above):

- For **RM**, the last column shows the **subfamily name** (e.g., `L1P1`, `AluY`).
- For **TRF**, the motif sequence may be displayed directly or replaced with `*80plus` if it exceeds 80 bp.

The information in the `diagram.txt` is available in `rm_diagram.tsv` and `trf_diagram.tsv` for easy parsing if necessary.

## Tracing an ID back to the raw tool output

The `rm_` and `trf_` IDs are sequential counters assigned as SVscanner reads the raw tool outputs, so they can be used to locate the original annotation line.

- **RM (`rm_N`)**: the counter increments once per line of the raw RepeatMasker output (the file passed to `--rm`, e.g. `rm.tab`), starting at `0` and with no header skipped. So `rm_N` corresponds to the **`N+1`-th line (1-based)**: e.g. `rm_19` is line 20, retrievable with `sed -n '20p' rm.tab` (or `awk 'NR==20'`).
- **TRF (`trf_N`)**: the counter increments per repeat line **only within the SV blocks SVscanner actually processed** (the `@`-delimited entries in the `--trf` file), skipping headers and unmatched blocks. It is therefore *not* a direct line number of the raw `.dat`; find the entry by its SV block and coordinates/motif instead.

These lookups are only valid against the exact `--rm`/`--trf` files used for that run — regenerating or re-sorting them changes the numbering.
