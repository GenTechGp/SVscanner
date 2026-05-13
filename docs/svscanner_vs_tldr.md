# SVscanner vs TLDR — Annotation Comparison

## Overview

Tool-to-tool annotation comparison between TLDR v1.3.0 and SVscanner v0.4.0 on HG002 ONT long-read data aligned to T2T-CHM13 (hs1). There is no ground truth — this comparison characterises where the two tools agree and disagree at the repeat class, family, and subfamily levels.

**Key difference in approach:** TLDR makes its own MEI calls directly from the BAM; SVscanner annotates an existing Sniffles2 call set. The two tools can therefore produce non-overlapping call sets even when both are correct.

---

## Data

- **Sample:** HG002 ONT long reads
- **Reference:** T2T-CHM13 hs1 (`hs1.fa`)
- **TLDR:** `tldr.table.txt` — all calls, all families, PASS and FAIL
- **SVscanner v0.4.0:** `run2/annotated.vcf.gz`
- **Sniffles2 v2.6.0 command:** `--reference hs1.fa --input sorted.haplotagged.bam --minsvlen 50 --phase`
- **SVscanner v0.3.0:** `annotated.vcf.gz` — Sniffles2 VCF annotated with RepeatMasker (Dfam) and TRF; key tag: `FINAL_CLASSIFICATION`

---

## Matching approach

Script: `compare_svscanner_tldr.py`

SVscanner INS records (SVTYPE=INS only; 15,753 of 29,146 total Sniffles2 calls) are indexed by chromosome and position. For each TLDR call, the script finds the closest unmatched SVscanner INS on the same chromosome within ±50 bp of the TLDR `Start` coordinate, using greedy 1:1 matching.

Three output zones:

| Zone | Definition |
|------|-----------|
| matched | TLDR call paired with an SVscanner INS call |
| tldr_only | TLDR call with no SVscanner INS within ±50 bp |
| svscanner_only | SVscanner INS classified as Mobile, with no matched TLDR call |

---

## Results

### Overview

| | Count |
|---|---|
| TLDR calls (all families, all filters) | 3,596 |
| SVscanner INS calls | 15,753 |
| Matched | 2,147 (59.7% of TLDR) |
| TLDR-only | 1,449 |
| SVscanner MEI-only | 683 |

### TLDR call breakdown by family

| Family | All calls | PASS |
|--------|----------|------|
| NA (unclassified) | 1,787 | 0 |
| ALU | 1,463 | 1,189 |
| L1 | 201 | 116 |
| SVA | 119 | 32 |
| ERV | 26 | 3 |

All NA-family calls are non-PASS. Concordance evaluation below is restricted to the four classified families (ALU, L1, SVA, ERV).

---

### Annotation concordance — matched calls

Concordance is assessed at the family level: TLDR family maps to an expected SVscanner `FINAL_CLASSIFICATION` substring (ALU → `SINE`, L1 → `LINE`, SVA → `Retroposon`, ERV → `LTR`).

- **correct** — `FINAL_CLASSIFICATION` contains the expected family keyword
- **ambiguous** — `Mixed` classification or bare `Repetitive/Mobile` / `Repetitive/Tandem` (repetitive content detected but family not resolved)
- **wrong** — `NON_REPETITIVE`, missing (`.`), or a clearly different family

#### PASS calls only

| Family | n | Correct | Ambiguous | Wrong |
|--------|---|---------|-----------|-------|
| ALU | 1,123 | 1,088 (97%) | 35 (3%) | 0 (0%) |
| L1 | 113 | 107 (95%) | 5 (4%) | 1 (1%) |
| SVA | 22 | 3 (14%) | 16 (73%) | 3 (14%) |
| ERV | 3 | 3 (100%) | 0 | 0 |

#### All filters

| Family | n | Correct | Ambiguous | Wrong |
|--------|---|---------|-----------|-------|
| ALU | 1,290 | 1,123 (87%) | 101 (8%) | 66 (5%) |
| L1 | 152 | 117 (77%) | 22 (14%) | 13 (9%) |
| SVA | 57 | 4 (7%) | 40 (70%) | 13 (23%) |
| ERV | 5 | 5 (100%) | 0 | 0 |

The drop from PASS-only to all-filters concordance is substantial for ALU (97% → 87%) and L1 (95% → 77%). TLDR FAIL calls have more `NON_REPETITIVE` and missing SVscanner classifications than PASS calls: ALU gains 55 `NON_REPETITIVE` wrongs and 57 bare `Repetitive/Mobile` ambiguous entries; L1 gains 5 `NON_REPETITIVE` and 20 bare `Repetitive/Mobile`.

#### Full classification breakdown — PASS calls

**ALU (n=1,123):**
- 1,077 → `Repetitive/Mobile/SINE` ✓
- 29 → `Repetitive/Mixed/TR` ~
- 11 → `Repetitive/Mixed/SINE` ✓
- 6 → `Repetitive/Mobile` ~

**L1 (n=113):**
- 107 → `Repetitive/Mobile/LINE` ✓
- 4 → `Repetitive/Mobile` ~
- 1 → `Repetitive/Mixed/TR` ~
- 1 → `.` ✗

**SVA (n=22):**
- 5 → `Repetitive/Mixed` ~
- 4 → `Repetitive/Mixed/VNTR` ~
- 3 → `Repetitive/Tandem/STR` ✗
- 3 → `Repetitive/Mobile` ~
- 3 → `Repetitive/Mobile/Retroposon` ✓
- 2 → `Repetitive/Tandem` ~
- 2 → `Repetitive/Mixed/TR` ~

**ERV (n=3):**
- 3 → `Repetitive/Mobile/LTR` ✓

The SVA concordance pattern (14% correct, 73% ambiguous) is consistent with SVA's VNTR domain: SVA insertions contain tandemly repeated sequence that TRF classifies independently of the RepeatMasker mobile element signal, producing Mixed or Tandem labels rather than Retroposon.

---

### SVscanner classification for matched TLDR NA-family calls

TLDR assigns NA when it cannot classify the insertion. For the 643 NA-family TLDR calls that positionally matched an SVscanner INS, SVscanner classified:

| SVscanner classification | Count |
|--------------------------|-------|
| `Repetitive/Mobile` | 152 |
| `Repetitive/Tandem/VNTR` | 112 |
| `NON_REPETITIVE` | 98 |
| `Repetitive/Mobile/LINE` | 88 |
| `Repetitive/Tandem/STR` | 52 |
| `Repetitive/Mixed/TR` | 23 |
| `Repetitive/Tandem/TR` | 20 |
| `Repetitive/Mobile/Retroposon` | 19 |
| `Repetitive/Mixed` | 18 |
| `Repetitive/Mobile/SINE` | 15 |
| `Repetitive/Mixed/VNTR` | 15 |
| `Repetitive/Mobile/LTR` | 13 |
| other / missing | 18 |

88 of these calls receive `LINE` and 19 receive `Retroposon` from SVscanner — indicating L1 and SVA insertions that TLDR detected but could not assemble well enough to assign a subfamily.

---

### TLDR-only calls (no SVscanner INS within ±50 bp)

| Family | All | PASS |
|--------|-----|------|
| NA | 1,144 | 0 |
| ALU | 173 | 66 |
| SVA | 62 | 10 |
| L1 | 49 | 3 |
| ERV | 21 | 0 |

1,449 TLDR calls have no matching Sniffles2 INS within ±50 bp. 1,144 are NA-family (all non-PASS).

#### Positional displacement

Widening the matching window resolves a portion of unmatched calls:

| Wiggle | Matched | % of TLDR |
|--------|---------|-----------|
| ±50 bp | 2,147 | 59.7% |
| ±100 bp | 2,267 | 63.0% |
| ±200 bp | 2,358 | 65.6% |
| ±500 bp | 2,434 | 67.7% |

287 additional calls match when the window is expanded from ±50 bp to ±500 bp. 1,162 TLDR calls remain unmatched even at ±500 bp.

#### Spanning read support — PASS classified calls

For the 79 PASS classified TLDR-only calls (66 ALU, 10 SVA, 3 L1), SpanReads (reads spanning the insertion breakpoints as reported by TLDR) are substantially lower than for matched PASS calls of the same family:

| Family | TLDR-only PASS | | | Matched PASS | | |
|--------|---------------|---|---|-------------|---|---|
| | n | median SpanReads | median LengthIns | n | median SpanReads | median LengthIns |
| ALU | 66 | 4 | 302 | 1,123 | 20 | 305 |
| L1 | 3 | 4 | 646 | 113 | 17 | 4,112 |
| SVA | 10 | 2 | 252 | 22 | 20 | 633 |

For ALU, insert sizes are similar between the two groups (median 302 bp vs 305 bp), so size is not a factor in the ALU mismatch. For SVA, unmatched PASS calls are substantially shorter (median 252 bp) than matched PASS calls (median 633 bp).

---

### SVscanner MEI-only calls (no matched TLDR call)

683 SVscanner INS records with a Mobile classification were not matched to any TLDR call:

| SVscanner classification | Count |
|--------------------------|-------|
| `Repetitive/Mobile` | 206 |
| `Repetitive/Mobile/SINE` | 184 |
| `Repetitive/Mobile/LINE` | 157 |
| `Repetitive/Mobile/LTR` | 107 |
| `Repetitive/Mobile/DNA` | 17 |
| `Repetitive/Mobile/Retroposon` | 12 |

184 SINE and 157 LINE insertions were detected by Sniffles2 and annotated by SVscanner but not matched to a TLDR call within ±50 bp — these represent insertions that TLDR either missed or placed outside the matching window. The 107 LTR calls have no TLDR counterpart for two reasons: TLDR detects only HERVK (all 26 ERV calls in the TLDR table are HERVK subfamily), while SVscanner's `Repetitive/Mobile/LTR` covers all LTR element classes (HERVK, HERV1, HERVL, MaLR/THE1, etc.); non-HERVK LTR insertions are outside TLDR's detection scope. The 12 Retroposon calls represent SVA insertions Sniffles2 detected but TLDR missed or displaced. The 206 bare `Repetitive/Mobile` entries are calls where SVscanner detected mobile element content but did not resolve to a specific family.

---

## Family and subfamily level comparison (SVscanner v0.4.0)

SVscanner v0.4.0 adds `RM_FAMILY` and `RM_SUBFAMILY` INFO fields — comma-separated lists parallel to `RM_CLASSIFICATION`, giving the RepeatMasker family (e.g. `Alu`, `L1`, `SVA`, `ERVK`) and subfamily/model name (e.g. `AluYa5`, `L1HS`, `SVA_F`) for each hit. This enables comparison with TLDR's `Family` and `Subfamily` columns at a finer level than FINAL_CLASSIFICATION alone.

All results below are for **PASS matched calls** only (n=1,161 across ALU, L1, SVA, ERV).

**Methodological note on dominant hit selection:** `RM_CLASSIFICATION`, `RM_FAMILY`, and `RM_SUBFAMILY` are parallel comma-separated lists — one entry per RepeatMasker hit. SVscanner can report multiple hits of the same class for a single call (e.g. `SINE,SINE` → `Alu,Alu` → `AluSx,AluY`). TLDR always reports a single family and subfamily per call. To compare, the script extracts the first RM hit whose class matches the class resolved in FINAL_CLASSIFICATION. When multiple same-class hits exist, this first-hit selection is arbitrary and introduces uncertainty into the subfamily comparison.

### RM_FAMILY concordance

RM_FAMILY concordance rates are consistent with the class-level FINAL_CLASSIFICATION rates.

| Family | n | Expected RM_FAMILY | Match | Mismatch | Unresolved |
|--------|---|-------------------|-------|----------|------------|
| ALU | 1,123 | `Alu` | 1,085 (97%) | 3 (0%) | 35 (3%) |
| L1 | 113 | `L1` | 107 (95%) | 0 (0%) | 6 (5%) |
| SVA | 22 | `SVA` | 3 (14%) | 0 (0%) | 19 (86%) |
| ERV | 3 | `ERVK` | 3 (100%) | 0 (0%) | 0 (0%) |

Unresolved entries are calls where FINAL_CLASSIFICATION did not resolve to a specific mobile class (e.g. `Repetitive/Mobile` or `Repetitive/Mixed/*`), so no dominant RM_FAMILY can be extracted.

### RM_SUBFAMILY concordance

TLDR's legacy L1 name `L1Ta` is normalised to the Dfam name `L1HS` before comparison.

| Family | n | Match | Mismatch | Unresolved |
|--------|---|-------|----------|------------|
| ALU | 1,123 | 119 (11%) | 969 (86%) | 35 (3%) |
| L1 | 113 | 77 (68%) | 30 (27%) | 6 (5%) |
| SVA | 22 | 0 (0%) | 3 (14%) | 19 (86%) |
| ERV | 3 | 0 (0%) | 3 (100%) | 0 (0%) |

#### ALU

The 86% apparent mismatch rate reflects a naming granularity difference rather than a classification disagreement. TLDR reports specific modern-clade names (AluYa5, AluYb8); the dominant SVscanner RM_SUBFAMILY for most calls is the ancestral consensus model `AluY`. Top mismatch pairs:

| Count | TLDR | SVscanner |
|-------|------|-----------|
| 638 | AluYa5 | AluY |
| 233 | AluYb8 | AluY |
| 28 | AluYb9 | AluY |
| 21 | AluYa5 | AluSx |
| 10 | AluYa5 | AluYe5 |

Note that the first-hit selection issue above applies here: for calls where SVscanner reports multiple Alu hits (e.g. `AluSx,AluY`), the reported mismatch may depend on hit order rather than a genuine difference.

#### L1

77/113 (68%) match after normalising L1Ta to L1HS. The 30 mismatches are calls where SVscanner's first L1 hit resolves to an older subfamily (L1PA2, L1P1, L1P2, L1PA4):

| Count | TLDR | SVscanner |
|-------|------|-----------|
| 10 | L1Ta | L1PA2 |
| 5 | L1Ta | L1P1 |
| 2 | L1Ta | L1P2 |
| 2 | L1Ta | L1PA4 |
| 1 | L1preTa | L1P1 |

#### SVA

19/22 calls are unresolved (no dominant Retroposon hit extracted), consistent with the RM_FAMILY result. Of the 3 mismatches: TLDR=SVA_A vs SVscanner=SVA_D, TLDR=SVA_E vs SVscanner=SVA_C, TLDR=SVA_F vs SVscanner=SVA_B.

#### ERV

TLDR assigns the family-level label `HERVK`; SVscanner RM_SUBFAMILY provides model-level names (`LTR5_Hs`, `HERVK-int`). These refer to the same element — the mismatch is a naming convention difference, not a classification disagreement.
