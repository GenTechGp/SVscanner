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
  - For RM: `sv_coverage%[element_proportion%]`, e.g., `84.11%[39.73%]`. `element_proportion` is the fraction of the element's consensus model covered by the part of the alignment that overlaps the SV (always ≤ 100%).
- **`REPEAT_ALIGNMNET_GENOMIC_RANGE`**: Start–end coordinates of the repeat's alignment
- **`ID`**: Unique identifier (e.g., `rm_668`, `trf_216`)  
  - A trailing `*` means this repeat **was selected at the pre-classification stage**
- **`rm_id=N`** (RM rows): RepeatMasker's own ID (column 15 of `rm.tab`). Rows with the same `rm_id` are fragments RepeatMasker joined into one element, e.g. the two halves of a reference element split by an insertion.

### RM flank

The `-RM-` rows list only hits that overlap the SV by more than `--min_sv_coverage`. A `-RM flank-` section before them lists the RepeatMasker hits near the SV that are left out: hits at least partly in the flanks and within `--diagram_flank_dist` bp of the SV (off by default; e.g. `--diagram_flank_dist 100`). These rows are drawn only. They are not used for any annotation and are not written to `rm_diagram.tsv`.

```
>INS.2350	(Sniffles2.INS.17BS1)	chr2:49634763-49636127	chr2:49635413-49635413	65	Repetitive/Tandem/VNTR
 SV & flanking      -----------------------------------------------[####]----------------------------------------------- 
 -RM flank-
  LINE                                                                    [::::::::::::::::::::::::::::]                 	chr2:49635441-49635831	rm_114735	rm_id=318	L1PA15	ref_cov=6.42%	dist=28bp
 -RM-
  LINE                  [...........................................]                                                    	16.92%[0.18%]	chr2:49634822-49635413	rm_114734*	rm_id=318
```

- `:` instead of `.`: a flank hit.
- Columns: genomic range, `ID`, `rm_id`, the subfamily name (there is no `SV Diagram` row for it), `ref_cov` and `dist`.
- `ref_cov`: the fraction of the element's consensus model the whole hit spans. Unlike `element_proportion` in the `-RM-` rows, it is not scaled to the part inside the SV.
- `dist`: the distance from the hit to the SV (0 if it reaches into the SV but covers too little of it).

Here the VNTR insertion sits inside a reference L1PA15: the hit left of the SV and the one 28 bp right of it share `rm_id=318`, one element RepeatMasker joined across the insertion. Flank hits are also useful for `NON_REPETITIVE` SVs, BND breakpoints in particular, and for SVs flanked by homologous repeats (e.g. the same Alu subfamily on both sides), a pattern typical of non-allelic homologous recombination (NAHR).

A second visualization shows a zoomed-in, **SV-only** view (no flanks) containing only the elements selected at the pre-classification stage (those marked `*` in the flanking section above):

- For **RM**, the last column shows the **subfamily name** (e.g., `L1P1`, `AluY`).
- For **TRF**, the motif sequence may be displayed directly or replaced with `*80plus` if it exceeds 80 bp.

The information in the `diagram.txt` is available in `rm_diagram.tsv` and `trf_diagram.tsv` for easy parsing if necessary.

## RM chain block

A third block, `RM chain`, redraws the same selected RM elements from RepeatMasker's own base-level alignments (`.align`, written with RepeatMasker `-a`). It shows how much of each element's consensus the SV matches, and where sequence the consensus does not have sits in the SV.

```
>INS.32	(Sniffles2.INS.32CS2)	chr3:147888048-147894775	chr3:147890048-147890048	2728	Repetitive/Mobile/Retroposon
 ...
 SV Diagram         [##################################################################################################] 
 -RM-
  Retroposon        [..................................................................................................] 	SVA_F

 RM chain           [##################################################################################################] 
 Retroposon         [.......**********************************************************.................................] 	SVA_F	cov=40.4%	ident=96.1%	ins=1596	ref_cov=80.1%	blocks=18
```

Here RepeatMasker reports one `SVA_F` covering 99% of the SV, but only 40% of the SV matches the consensus: 1.6 kb is extra VNTR copies, an expansion of the SVA's VNTR, drawn where it sits.

**Why it is needed.** RepeatMasker reports one element as several alignment blocks under one ID, e.g. one block per VNTR copy of an expanded SVA, each aligned to the same stretch of consensus, or an L1 split over its `_5end`/`_orf2`/`_3end` models. The `.out` line only gives the outer range. SVscanner picks one path through the blocks: query and consensus both increase in the element's direction, as many bases as possible are matched, and each jump between blocks costs `--rm_chain_jump_cost` matched bases (default 30). Query bases inside the element that are not on the path are **inserted**.

- **Row**: one per selected element. Several selected `.out` lines of one joined element (same RepeatMasker ID) get one row.
- **`[...]`**: the element's extent over the SV, with the same scaling and bracket rules as the `SV Diagram` block (no bracket where the element continues past the SV).
- **`*`**: inserted stretches of ≥ 50 bp.
- **`cov`**: SV bases matched to the consensus / SV length.
- **`ident`**: matched SV bases identical to the consensus / matched SV bases.
- **`ins`**: SV bases inside the element that are not matched (inserted), in bp.
- **`ref_cov`**: consensus positions matched by SV bases / consensus length (of the model on the `.out` line).
- **`blocks`**: alignment blocks chained.

**Reading it**

- `[.....]` with no `*` and `ref_cov` near 100%: a full-length element that matches its consensus end to end.
- `*` inside the row: the element carries sequence its consensus does not have, e.g. a VNTR expansion (SVA), or another sequence it was split around.
- All `*`, no brackets (e.g. `cov=0.0%`, `ins` = SV length): the SV lies inside a larger element that RepeatMasker joined across it. For an INS this is typically extra copies of an existing element's tandem repeat, or a new insertion into an older element.
- `cov` well below RM's `sv_coverage` and `ref_cov` well below 100%: a 5'-truncated or partial element.

**Consensus coordinates.** Blocks are placed on the coordinates the `.out` line reports: RepeatMasker already puts an L1's `_5end`/`_orf2`/`_3end` models on one full-length model (e.g. `L1M5_orf2` → `L1M5`). Sibling subfamilies under one ID (SVA_A–F tiles under an `SVA_F` line) share coordinates, which differ by only a few bp. An ERV joined from an LTR and its `-int` model is placed LTR, then `-int`, then LTR. A tail of another class or strand joined to an element (an Alu's `(T)n`) is left out of the chain.

The values come from `rm_chain.tab` in the output directory, one row per SV and RepeatMasker ID whose alignments touch the SV:

| Column | Description |
|---|---|
| `sv_id` | SV ID (as in `info.tab`) |
| `rm_id` | RepeatMasker ID: the last column of the `.out` line, column 15 of `rm.tab` |
| `repeat`, `class_family` | Name and class/family of the element's longest `.out` line |
| `strand` | `+` or `-` |
| `models` | Consensus models of the chained blocks |
| `n_blocks`, `n_left_out` | Blocks chained; blocks left out (another class or strand) |
| `element_start`, `element_end` | Element extent in query coordinates (1-based, like `rm.tab`) |
| `consensus_len` | Length of the `.out` line's consensus model |
| `matched`, `identical`, `inserted` | SV bases matched, of those identical to the consensus, and inserted |
| `consensus_hit` | Consensus positions matched by SV bases |
| `ins_blocks` | Inserted stretches of ≥ 50 bp: `start-end@consensus_position` (query coordinates; the consensus position reached when the stretch starts), `;`-separated |

The block is on by default; `run_workflow.sh --no_rm_chain` turns it off (RepeatMasker then runs without `-a`). Runs made with `--no_rm_chain`, or before this block was added, have no `rm_chain.tab`; `--resume` on them draws diagrams without it.

## Tracing an ID back to the raw tool output

The `rm_` and `trf_` IDs are sequential counters assigned as SVscanner reads the raw tool outputs, so they can be used to locate the original annotation line.

- **RM (`rm_N`)**: the counter increments once per line of the raw RepeatMasker output (the file passed to `--rm`, e.g. `rm.tab`), starting at `0` and with no header skipped. So `rm_N` corresponds to the **`N+1`-th line (1-based)**: e.g. `rm_19` is line 20, retrievable with `sed -n '20p' rm.tab` (or `awk 'NR==20'`).
- **TRF (`trf_N`)**: the counter increments per repeat line **only within the SV blocks SVscanner actually processed** (the `@`-delimited entries in the `--trf` file), skipping headers and unmatched blocks. It is therefore *not* a direct line number of the raw `.dat`; find the entry by its SV block and coordinates/motif instead.

These lookups are only valid against the exact `--rm`/`--trf` files used for that run — regenerating or re-sorting them changes the numbering.

- **RM chain rows** are keyed by `sv_id` and `rm_id`, RepeatMasker's own ID (column 15 of `rm.tab`). That ID is numbered per RepeatMasker run, so it is only unique together with the SV ID: e.g. `awk '$5=="INS.32" && $15==7' rm.tab` lists the `.out` lines of the element in the `rm_chain.tab` row with `sv_id` `INS.32` and `rm_id` `7`.
