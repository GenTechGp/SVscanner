# Repeat/SV Diagram

This document describes the format used to represent structural variants (SVs) with annotated **TRF (Tandem Repeats)** and **RM (RepeatMasker)** elements, highlighting their overlap with the variant region and flanking sequences. The terms "repeat," "element," "repeat element," and "mobile element" are used interchangeably in this context.

Each record starts with a line beginning with a `>` symbol.
- **`SV_ID`**: A unique identifier for the structural variant.
- **`SOURCE_ID`**: Original SV ID from caller, e.g., Sniffles.
- **`GENOMIC_REGION`**: Region of the genome spanned by the SV including flanks (e.g., `chr1:8692-12927`).
- **`SV_REGION`**: The precise SV site or inserted sequence location (e.g., `chr1:10692-10927`).
- **`SV_LENGTH`**: The length of the SV in base pairs.



```

>DEL.0	(0_Sniffles2.DEL.5624S0)	chr1:8692-12927	chr1:10692-10927	235	
 SV & flanking      -----------------------------------------------[####]----------------------------------------------- 
 -TRF-
  29                                                             [...]                                                   	46.38%	chr1:10629-10801	trf_216
  76                                                                [.....]                                              	71.49%	chr1:10759-10998	trf_217

 SV Diagram         [##################################################################################################] 
 -TRF-
  29                ..............................................]                                                      	AGGCGCGCCGCGCCGGCGCAGGCGCAGAG
  76                                            [....................................................................... 	GGCGCAGGCGCAGAGAGGCGCGCCGCGCCGGCGCAGGCGCAGAGACACATGCTAGCGCGTCCAGGGGTGGAGGCGT

>INS.16	(0_Sniffles2.INS.33S0)	chr1:430868-435559	chr1:432868-432868	692	
 SV & flanking      ------------------------------------------[##############]------------------------------------------ 
 -RM-
  LINE                                                         []                                                        	11.42%[1.28%]	chr1:432868-432868	rm_294
  LINE                                                          [.]                                                      	7.95%[0.89%]	chr1:432868-432868	rm_295
  LINE                                                           [.]                                                     	10.55%[1.18%]	chr1:432868-432868	rm_296
  LINE                                                             [.]                                                   	11.42%[1.28%]	chr1:432868-432868	rm_297
  LINE                                                               []                                                  	11.42%[1.28%]	chr1:432868-432868	rm_298
  LINE                                                                [.]                                                	11.27%[1.26%]	chr1:432868-432868	rm_299
  LINE                                                                  [.]                                              	11.56%[1.3%]	chr1:432868-432868	rm_300
  LINE                                                                   [.]                                             	14.02%[1.57%]	chr1:432868-432868	rm_301
  LINE                                                                     [.]                                           	11.71%[1.31%]	chr1:432868-432895	rm_302
 -TRF-
  77                                                    [........]                                                       	21.97%	chr1:432584-432868	trf_75*
  123                                                          [....]                                                    	34.68%	chr1:432868-432868	trf_76
  77                                                             [...........]                                           	77.46%	chr1:432868-432868	trf_77*
  105                                                                      [...]                                         	15.46%	chr1:432868-432973	trf_80
  182                                                                    [.......]                                       	26.59%	chr1:432868-433050	trf_81

 SV Diagram         [##################################################################################################] 
 -RM-
  LINE                 [...........]                                                                                     	L1P4
  LINE                            [.......]                                                                              	L1P3
  LINE                                    [.........]                                                                    	L1P4
  LINE                                              [...........]                                                        	L1P4
  LINE                                                         [...........]                                             	L1P4
  LINE                                                                    [...........]                                  	L1P4
  LINE                                                                                [..........]                       	L1P4
  LINE                                                                                        [.............]            	L1P3
  LINE                                                                                                      [........... 	L1P3
 -TRF-
  77                .....................]                                                                               	CAAACACGTGGATACATGGAGGGGAACAACACACACCAGGGCCTCTCAGGGGGACAGGGGGTAGGAGACCATCAGGA
  123                    [.................................]                                                             	*80plus
  77                                      [............................................................................] 	TGGGTACATGGAGGGGAACAACACACACCAGGGCCTCTCAGGGGGACAGGGGGTAGGAGACCATCAGGACAAACACG
  105                                                                                                   [............... 	*80plus
  182                                                                                        [.......................... 	*80plus


```
Immediately following the header is a **flanking diagram** that visualizes the SV relative to its surrounding sequence:
- `-`: Flanking genomic sequence  
- `[####]`: The SV itself (e.g., deletion, insertion, etc.)

Annotations from TRF and RM tools are listed below the flanking diagram. They indicate repetitive elements that **intersect with the SV** or flanking sequence.

- **`PERIOD`**: Repeat unit period (TRF) or class (RM like LINE, SINE)
- **`[REPEAT_ALIGNMNET_SHAPE]`**: ASCII diagram of the repeat’s alignment position relative to the flanking/SV diagram
- **`PERCENT`**:  
  - For TRF: `sv_coverage%`  
  - For RM: `sv_coverage%[element_coverag%]`, e.g., `84.11%[39.73%]`
- **`REPEAT_ALIGNMNET_GENOMIC_RANGE`**: Start–end coordinates of the repeat's alignment
- **`ID`**: Unique identifier (e.g., `rm_668`, `trf_216`)  
  - A trailing `*` means this repeat **was selected** during final classification

A second visualization shows a zoomed-in, **SV-only** view (no flanks):

- For **RM**, the last column shows the **family name** (e.g., `L1P1`, `AluY`).
- For **TRF**, the motif sequence may be displayed directly or replaced with `*80plus` if it exceeds 80 bp.

The information in the `diagram.txt` is available in `rm_diagram.tsv` and `trf_diagram.tsv` for easy parsing if necessary.
