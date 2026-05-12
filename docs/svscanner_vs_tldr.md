# SVscanner vs TLDR — Annotation Comparison

## Overview

Tool-to-tool annotation comparison between TLDR v1.3.0 and SVscanner v0.3.0 on HG002 ONT long-read data aligned to T2T-CHM13 (hs1). There is no ground truth — this comparison characterises where the two tools agree and disagree at the repeat family level, which is the finest common granularity.

**Key difference in approach:** TLDR makes its own MEI calls directly from the BAM; SVscanner annotates an existing Sniffles2 call set. The two tools can therefore produce non-overlapping call sets even when both are correct.

---

## Data

- **Sample:** HG002 ONT long reads
- **Reference:** T2T-CHM13 hs1 (`hs1.fa`)
- **TLDR:** `tldr.table.txt` — all calls, all families, PASS and FAIL
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
