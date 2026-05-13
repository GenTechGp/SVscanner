# Repeat Annotation Steps

The steps followed in annotating the vcf records using TRF and RM information are given below. For VCF INFO tags refer [repeat_annotation_VCF_tags.md](repeat_annotation_VCF_tags.md)

## SV/Repeat overlap calculations

![Illustration](/images/annotation_tag_illustration.png)


* RM\_ELEMENTS\_COVERAGE = `ri_sv_l/ri_L` e.g. [`r1_sv_l/r1_L`, `r2_sv_l/r2_L`]
* RM\_ELEMENT\_PROPORTION = `ri_l/ri_L` e.g. [`r1_l/r1_L`, `r2_l/r2_L`]
* RM\_SV\_COVERAGE = `ri_sv_l/sv_L` e.g. [`r1_sv_l/sv_L`, `r2_sv_l/sv_L`]
  - TRF_SV_COVERAGE is also calculated similarly.
* RM\_TOTAL\_SV\_COVERAGE = `sum(ri_sv_l)/sv_L` e.g. `(r1_sv_l+r2_sv_l)/sv_L`
   - TRF_TOTAL_SV_COVERAGE is also calculated similarly.

### Note
1.  `ri_l` is not necessarily equal to `ri_sv_l`. It can be less than, equal or greater than `ri_sv_l`.
2. The dotted lines in the figure show the boundaries of the alignment.
3. Diagrams visualizing Repeat/SV annotations ([diagram file format](Repeat-SV_diagram.md)) show only the parts of the repeat elements that align to the SV.

## Important parameters with default values

| Parameter | Default | Description |
|----------|---------|-------------|
| min_sv_coverage | 0.05 | The minimum intersection between a repeat element and SV (aka sv_coverage) e.g. 0.05 (5%) (0 < min_sv_coverage < 1). |
| min_total_sv_coverage | 0.75 | The tight minimum total sv coverage by repeat elements to be considered repetitive. |
| sv_coverage_cutoff | 0.5 | The soft minimum total sv coverage by repeat elements to be considered repetitive. |
| min_class_sv_coverage | 0.25 | The minimum class sv coverage by repeat elements to be considered repetitive. |
| div | 0.05 | The chosen intervals to prioritise period size over intersection (0 < divisor < 1). |
| max_trf_overlap | 0.1 | The maximum TRF element overlap to be considered non-overlapping (0 < max_trf_overlap < 1). |


## Step 1 (TRF preprocessing)

For **Tandem Repeat Finder (TRF)**, entries are determined by prioritising maximal intersection between SV and repeat and minimal period size in intervals (e.g. 0.05). Non-overlapping entries are selected in order of priority.

![Illustration](/images/TRF_workflow.png)

1. Parse `sample_trf.tab` file:
   - The file can have multiple lines (elements) per `sv_id`.
   - For each line, save the element under the corresponding `sv_id`.

2. Determine each TRF element's classification using table below.

    | Condition               | Classification |
    |------------------------|----------------|
    | period_size == 1       | HOMO           |
    | 1 < period_size < 7    | STR            |
    | 6 < period_size < 101  | VNTR           |
    | 100 < period_size      | TR             |

4. Calculate each element's `sv_coverage`.

5. Filter out elements where `sv_coverage < args.min_sv_coverage`.

6. Store the remaining elements in an array for each `sv_id`:
   - Do not separate elements by classification; store all together under the `sv_id`.


## Step 2 (TRF filtering)
1. For each `sv_id`:
   - Get the `trf` elements list created above.
   - Sort the repeat elements in the list based on the element's repeat start index.

2. Calculate total SV coverage (TC).

3. Group the elements into `args.div%` bins based on the element's `sv_coverage`:
   - Within each bin, sort elements from smallest `period_size` to largest.

4. Sort the bins in decreasing order:
   - Elements with higher `sv_coverage` are prioritised.
   - Within the same bin, smaller `period_size` elements are prioritised.

5. Create a list of non-overlapping elements using the sorted list above:
   - Allow overlaps up to `args.max_trf_overlap`.

6. Sort the non-overlapping list using each element's repeat start index.
   - Call this list `TRF_LIST`.

7. Calculate total non-overlapping elements SV coverage (`TRF_TTC`) using `TRF_LIST`.


## Step 3 (RM preprocessing)
For **RepeatMasker (RM)**, entries are determined by prioritising those with maximal intersection between the SV and repeat entry. Entries are grouped into their repeat class (e.g. SINE, LINE) and in order of priority non-overlapping entries are selected within each class. Based on the element coverage and SV coverage, the SV is classified as a Full/Partial RECIPROCAL or Minimal.

![Illustration](/images/RM_workflow.png)

1. Go through each row (element) of the `SAMPLE_rm.tab` file.

2. Obtain the element's classification (also known as repeat class), e.g., SINE, LINE, LTR, DNA, Retroposon (col 11 of `.out`, part before `/`). Also extract:
   - `rm_family`: col 11 part after `/` (e.g. `L1` from `LINE/L1`); set to `.` if col 11 has no `/` (bare class).
   - `rm_subfamily`: col 10 repeat name (e.g. `L1M4`, `AluYa5`).

3. Filter out elements that are not one of the following:
   - SINE, LINE, LTR, DNA, Retroposon

4. Calculate the element's `sv_coverage`.

5. Filter out elements where `sv_coverage < args.min_sv_coverage`.

6. Calculate the element's `element_coverage`.

7. Categorise and store elements in a dictionary:
   - Key: classification (e.g., "SINE", "LINE", etc.)
   - Value: list of elements belonging to that classification.

## Step 4 (RM filtering)
1. Iterate through the populated dictionary key by key (each key is a classification).

2. For each classification:
   - Obtain the element list for that classification.
   - Sort the element list in decreasing order of `sv_coverage`.
   - Iterate through the sorted list and pick elements that do not overlap:
     - The element with the highest `sv_coverage` is always picked.
   - Calculate the total SV coverage (TC) using the selected non-overlapping elements for that classification.
   - If every non-overlapping element's `element_coverage` >= 0.75, then Element Coverage Complete (ECC) is `True`.
   - Determine `RECIPROCAL` for the classification using **Table 1**.
   - Store `RECIPROCAL` for the classification.

3. After processing all classifications:
   - Concatenate all selected non-overlapping elements (from all classifications) into a list called `RM_LIST`.
   - Note: `RM_LIST` may contain elements that overlap across different classifications.

4. Sort `RM_LIST` by each element's repeat start.

5. Calculate the total SV coverage `RM_TTC` using `RM_LIST` (classification is not considered for this total).

## Step 5 (Decide Final Classification)
1. Prepare `RM` tags and `TRF` tags as shown in **Table 2**.

2. Check if `RM_TTC < args.min_total_sv_coverage`:
   - If true, set `RM_CLASSIFICATION = NON_REPETITIVE`.
   - Set every `RM` tag to `NA`.

3. Check if `TRF_TTC < args.min_total_sv_coverage`:
   - If true, set `TRF_CLASSIFICATION = NON_REPETITIVE`.
   - Set every `TRF` tag to `NA`.

4. Set `RM_TOTAL_SV_COVERAGE = RM_TTC`.

5. Set `TRF_TOTAL_SV_COVERAGE = TRF_TTC`.

6. Use **Table 3** to determine:
   - `PRE_CLASSIFICATION`
   - `RECIPROCAL`

7. Use **Table 8** to determine:
  - `FINAL_CLASSIFICATION`

8. Print the following:
   - `RM` tags
   - `RM_TOTAL_SV_COVERAGE`
   - `TRF` tags
   - `TRF_TOTAL_SV_COVERAGE`
   - `RECIPROCAL`
   - `PRE_CLASSIFICATION`

9. In addition, if a STRchive bed file is provided and the SV intersects with the position of a gene then following tags are also added to INFO
   - `DISEASE_GENE` STR disease associated with gene
   - `STRCHIVE_MOTIF` Pathogenic motif
   - `PATHOGENIC_MIN` Minimum pathogenic number annotated in STRchive


### Table 1: Determine ECC
| Condition         | ECC = True | ECC = False |
|---------------------------|------------|-------------|
| TC ≥ 0.75                 | Full       | Partial     |
| TC < 0.75                 | Minimal    | Minimal     |

### Table 2: VCF INFO tags
| Key                     | Value (comma-separated string) |
|-------------------------|--------------------------------|
| RM_CLASSIFICATION       | Classification (repeat class) of each element in the RM_LIST |
| RM_FAMILY               | Family of each element in the RM_LIST (e.g. L1, Alu). `.` if no family sub-classification exists for that element. |
| RM_SUBFAMILY            | Subfamily (repeat model name) of each element in the RM_LIST (e.g. L1M4, AluYa5). |
| RM_SV_COVERAGE          | `sv_coverage` of each element in the RM_LIST (sum of these values can be > 1) |
| RM_ELEMENT_COVERAGE     | `element_coverage` of each element in the RM_LIST (sum of these values can be > 1) |
| RM_ELEMENT_PROPORTION   | `element_proportion` of each element in the RM_LIST (sum of these values can be > 1) |
| TRF_CLASSIFICATION      | Classification of each element in the TRF_LIST |
| TRF_SV_COVERAGE         | `sv_coverage` of each element in the TRF_LIST |
| TRF_PERIOD_SIZE         | Period size of each element in the TRF_LIST |
| TRF_COPY_NUMBER         | Copy number of each element in the TRF_LIST |
| CONSENSUS_REPEAT        | Motif of each element in the TRF_LIST |

### Table 3
 TRF_CLASSIFICATION | RM_CLASSIFICATION | PRE_CLASSIFICATION | RECIPROCAL |
|--------------------|-------------------|-----------------------|-----------------------|
| NON REPETITIVE/NA | NON REPETITIVE | NON REPETITIVE/NA | determine using **Table 6** |
| valid | NON REPETITIVE/NA | determine using **Table 4** | determine using **Table 6** |
| NON REPETITIVE/NA | valid | determine using **Table 4** | determine using **Table 6** |
| valid | valid | calculate for both TRF and RM separately using **Table 4** (let's call them TRF_PRE_CLASSIFICATION and RM_PRE_CLASSIFICATION respectively) and then select using **Table 7** | determine using **Table 6** |

### Table 4
| condition | PRE_CLASSIFICATION |
|-----------|----------------------|
| highest class TC < `args.min_class_sv_coverage` | NON REPETITIVE |
| else | determine using **Table 5** for the class with highest TC |

### Table 5
| Class | RECIPROCAL | PRE_CLASSIFICATION |
|-------|------------|----------------------|
| SINE/LINE/LTR/DNA/Retroposon | Full/Partial | Class |
| SINE/LINE/LTR/DNA/Retroposon | Minimal | NON REPETITIVE |
| HOMO/STR/VNTR/TR | NA | Class |

### Table 6
| PRE_CLASSIFICATION | RECIPROCAL |
|----------------------|------------|
| SINE/LINE/LTR/DNA/Retroposon | RECIPROCAL of the class |
| HOMO/STR/VNTR/TR | NA |

### Table 7
| TRF_PRE_CLASSIFICATION | RM_PRE_CLASSIFICATION | PRE_CLASSIFICATION |
|--------------------------|-------------------------|----------------------|
| NON REPETITIVE | NON REPETITIVE | NON REPETITIVE |
| NON REPETITIVE | valid | RM_PRE_CLASSIFICATION |
| valid | NON REPETITIVE | TRF_PRE_CLASSIFICATION |
| valid | valid | if RM TTC > TRF TTC and RECIPROCAL is Full then RM_PRE_CLASSIFICATION else TRF_PRE_CLASSIFICATION |


### Table 8
| PRE_CLASSIFICATION | `TRF_TTC` = x | `RM_TTC` = y | Region | Final classification |
|--------------------|---------------|---------------|--------|----------------------|
| NON_REPETITIVE     | 0 ≤ x < S     | 0 ≤ y < S     | A      | NON_REPETITIVE       |
|                    | S ≤ x ≤ 1     | 0 ≤ y < S     | B, C   | Repetitive / Tandem  |
|                    | 0 ≤ x < S     | S ≤ y < T     | D      | Repetitive / Mobile  |
|                    | S ≤ x < T     | S ≤ y < T     | E      | Repetitive / Mixed   |
|                    | T ≤ x ≤ 1     | S ≤ y < T     | F      | Repetitive / Tandem  |
|                    | 0 ≤ x < T     | T ≤ y ≤ 1     | G, H   | Repetitive / Mobile  |
|                    | T ≤ x ≤ 1     | T ≤ y ≤ 1     | I      | Repetitive / Mixed   |
| M                  | 0 ≤ x < T     | T ≤ y ≤ 1     | G, H   | Repetitive / Mobile / M |
| N                  | T ≤ x ≤ 1     | 0 ≤ y < T     | F, C   | Repetitive / Tandem / N |
| M or N             | T ≤ x ≤ 1     | T ≤ y ≤ 1     | I      | Repetitive / Mixed / M or N |
| M or N             | 0 ≤ x < T     | 0 ≤ y < T     | A, B, D, E | Error!            |


**Note:**  
- S = `args.sv_coverage_cutoff` (0.5)  
- T = `args.min_total_sv_coverage` (0.75)
- if `TRF_TTC` is missing then set it to 0
- if `RM_TTC` is missing then set it to 0
- M = if PRE_CLASSIFICATION in {DNA,Retroposon,SINE,LINE,LTR}
- N = if PRE_CLASSIFICATION in {TR,VNTR,STR,HOMO}


![Illustration](/images/final_classification_regions.png)
