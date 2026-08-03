# RM element selection: greedy vs minimum interval cover DP

When annotating a structural variant (SV) with RepeatMasker (RM) hits, multiple partially-overlapping hits of the same repeat class may cover the SV. A selection method is needed to choose a representative subset for two purposes:

1. **Reporting** — selected elements appear in VCF tags (`RM_CLASSIFICATION`, `RM_SV_COVERAGE`, etc.)
2. **Coverage metric** — `RM_TOTAL_SV_COVERAGE` is the union of selected elements' spans, computed via a sweep-line algorithm

The selection method determines both what is reported and whether coverage is accurately captured.

---

## Problem definition

Given a set of RM intervals on the SV, find the subset S such that:

1. **Primary**: union_coverage(S) is maximised — S accounts for as much of the SV as possible
2. **Secondary**: |S| is minimised — no redundant elements are reported
3. **Tie-break**: among all minimum-cardinality sets with maximum coverage, prefer the one with the highest total `sv_coverage` weight (sum of individual element sv_coverage values) — this selects the highest-quality alignments when alternatives exist

Note: the weight used for tie-breaking (sum of individual `sv_coverage` values) can exceed 1.0 when selected elements overlap. This is intentional — it measures selection quality, not coverage. The final `RM_TOTAL_SV_COVERAGE` is always computed via union (sweep-line), never from the weight sum.

---

## Illustrative example

An SV of length 1000 bp with three RM hits of the same repeat class:

```
SV:  |---------------------------------------------|  [0 : 1000]

A:   |=================================|                 [0 : 700]   sv_cov = 0.70
B:                         |=======================|  [500 : 1000]  sv_cov = 0.50
C:               |=======================|            [300 : 800]  sv_cov = 0.50
```

Pairwise overlaps:
- A ∩ B = [500 : 700] = 200 bp
- A ∩ C = [300 : 700] = 400 bp
- B ∩ C = [500 : 800] = 300 bp

Union of all three = [0 : 1000] = **100%**.

The minimum-cardinality subset achieving 100% is **{A, B}** (count = 2):
- {A, B}: [0:700] ∪ [500:1000] = [0:1000] = 100% ✓
- {A, C}: [0:800] = 80% — does not reach 100%
- {B, C}: [300:1000] = 70% — does not reach 100%
- {A, B, C}: 100% but count = 3 — not minimum

---

## Method 1: Greedy + unique coverage threshold

### Algorithm

1. Sort elements by `sv_coverage` descending.
2. Maintain a running set of already-covered SV positions.
3. For each candidate (in sorted order):
   - Compute its **unique contribution**: SV positions it covers not already covered.
   - Accept if `unique_contribution / sv_length ≥ min_sv_coverage` (default 5%).
4. Compute `RM_TOTAL_SV_COVERAGE` as the union of accepted elements.

### Applied to the example

Sort order: A (0.70), then B and C tied at (0.50). Order of B and C is arbitrary.

**Ordering 1 — B before C:**

| Step | Candidate | Already covered | Unique region | Unique / sv_len | Accepted? |
|------|-----------|----------------|---------------|----------------|-----------|
| 1 | A [0:700] | — | [0:700] = 700 bp | 70% | ✓ |
| 2 | B [500:1000] | [0:700] | [700:1000] = 300 bp | 30% | ✓ |
| 3 | C [300:800] | [0:1000] | 0 bp | 0% | ✗ |

**Selected: {A, B}** — 2 elements, 100% ✓ (optimal)

**Ordering 2 — C before B:**

| Step | Candidate | Already covered | Unique region | Unique / sv_len | Accepted? |
|------|-----------|----------------|---------------|----------------|-----------|
| 1 | A [0:700] | — | [0:700] = 700 bp | 70% | ✓ |
| 2 | C [300:800] | [0:700] | [700:800] = 100 bp | 10% | ✓ |
| 3 | B [500:1000] | [0:800] | [800:1000] = 200 bp | 20% | ✓ |

**Selected: {A, C, B}** — 3 elements, 100% ✗ (not minimum)

### The problem

B and C have identical `sv_coverage` (0.50). The greedy sort order between them is arbitrary. Depending on which is processed first, the result is either the optimal 2-element set or a suboptimal 3-element set. Greedy is **order-dependent**: the order of evaluation determines which element "claims" the contested [700:800] region, changing whether the other is needed.

---

## Method 2: Minimum interval cover DP

### Overview

The DP explores all possible selection paths simultaneously and finds the globally optimal set regardless of element ordering. It operates in two phases.

### Phase 1 — decompose into connected components

Compute the union of all elements to find maximal contiguous covered regions. Gaps (positions covered by no element) separate independent components. Run the DP on each component separately.

For our example: one component [0:1000] (no gaps, all positions covered by at least one element).

For the SVA2 real case (see below): two components — [19:1426] and [1681:2701] — separated by a gap at [1426:1681].

### Phase 2 — DP per component

For a component spanning `[a, b]`:

**State:** `dp[x] = (min_count, max_weight)` — minimum elements to cover `[a, x]`, and the maximum total `sv_coverage` weight achievable with that count.

**Initialisation:** `dp[a] = (0, 0)`

**Transition:** sort elements by right endpoint. For each element `e = [l, r]`:
- Find the best predecessor `dp[x']` for any defined `x'` in `[l, r)`.
- The condition `x' ≥ l` ensures no gap is left between existing coverage and where `e` starts.
- Among valid predecessors, pick the one with smallest count; break ties by largest weight.
- Update: `dp[r] = (best.count + 1, best.weight + e.sv_cov)` if this improves current `dp[r]`.

**Answer:** `dp[b]`. Reconstruct selected elements by tracing predecessor links back to `dp[a]`.

**Complexity:** O(n²) — adequate since n (elements per SV) is small in practice.

### Weight in the DP

The weight accumulates the sum of individual `sv_coverage` values of selected elements. It can exceed 1.0 when selected elements overlap — this is expected. The weight is used only for tie-breaking during selection; it is never reported. The final `RM_TOTAL_SV_COVERAGE` is always recomputed from the selected elements using the sweep-line union.

Using the actual union coverage as the weight does not work for tie-breaking: all minimum-count paths reaching `dp[r]` cover the same range `[a, r]` by construction, so their union coverage is identical and provides no differentiation.

### Applied to the example

Sort elements by right endpoint: A (700), C (800), B (1000).

`dp[0] = (0, 0)`

**Process A = [0, 700]:**
- Valid predecessors x' ∈ [0, 700): dp[0] = (0, 0) ✓
- dp[700] = **(1, 0.70)**

**Process C = [300, 800]:**
- Valid predecessors x' ∈ [300, 800): dp[700] = (1, 0.70) ✓
- dp[800] = **(2, 1.20)**

**Process B = [500, 1000]:**
- Valid predecessors x' ∈ [500, 1000): dp[700] = (1, 0.70) and dp[800] = (2, 1.20)
- Best: dp[700] = **(1, 0.70)** — count 1 < count 2
- dp[1000] = **(2, 1.20)**

**Answer: dp[1000] = (2, 1.20)**

Trace back: last element added to reach dp[1000] is B (predecessor dp[700]); element added to reach dp[700] is A (predecessor dp[0]).

**Selected: {A, B}** — 2 elements, 100% ✓ (optimal, regardless of element ordering)

The DP finds dp[700] as the best predecessor for B, meaning "A already covers [0:700]; adding B alone reaches 1000." It does not need C because C's contribution is fully subsumed once both A and B are selected.

### Comparison on the illustrative example

| Method | Selected elements | Count | RM_TOTAL_SV_COVERAGE | Optimal? |
|--------|-----------------|-------|---------------------|----------|
| Greedy (B before C) | A, B | 2 | 100% | ✓ |
| Greedy (C before B) | A, C, B | 3 | 100% | ✗ |
| **DP** | **A, B** | **2** | **100%** | **✓ always** |

---

## Real case: SVA VNTR insertion (SVscanner Bug 1)

For variant `Sniffles2.INS.D2S0` (chr1:73679867, SVLEN=2748), four Retroposon hits are produced. RepeatMasker generates a chain of partially-overlapping SVA2 hits across the expanded VNTR domain — adjacent hits overlap at the VNTR boundary because the aligner bridges the expansion (conserved flanking sequences on both sides anchor the alignment, causing extra tandem copies to appear as large query insertions rather than stopping the alignment).

```
SV:       |--------------------------------------------------|  [0 : 2748]

rm_5044:  |==========|                                        [19  : 642]   sv_cov = 0.227
rm_5045:           |===========|                              [506 : 1230]  sv_cov = 0.264
rm_5046:                       |====|                         [1094: 1426]  sv_cov = 0.121
                                    gap [1426:1681]
rm_5048:                               |================|     [1681: 2701]  sv_cov = 0.371
```

*(Positions shown relative to SV start. Gap [1426:1681] = 255 bp is covered by no RM hit.)*

### Phase 1 — connected components

- Component 1: [19 : 1426] — elements rm_5044, rm_5045, rm_5046
- Component 2: [1681 : 2701] — element rm_5048

### Phase 2 — DP on component 1 ([19, 1426])

Sort by right endpoint: rm_5044 (642), rm_5045 (1230), rm_5046 (1426).

`dp[19] = (0, 0)`

| Element | Valid predecessor x' | dp update |
|---------|---------------------|-----------|
| rm_5044 [19:642] | dp[19]=(0,0), x'=19 ≥ 19 ✓ | dp[642] = **(1, 0.227)** |
| rm_5045 [506:1230] | dp[642]=(1,0.227), x'=642 ≥ 506 ✓ | dp[1230] = **(2, 0.491)** |
| rm_5046 [1094:1426] | dp[1230]=(2,0.491), x'=1230 ≥ 1094 ✓ | dp[1426] = **(3, 0.612)** |

All three elements are essential — each uniquely covers a segment no other element reaches.

Trace back: rm_5046 ← rm_5045 ← rm_5044. **Selected: {rm_5044, rm_5045, rm_5046}**

### DP on component 2 ([1681, 2701])

Only rm_5048. dp[2701] = (1, 0.371). **Selected: {rm_5048}**

### Result

**Selected: {rm_5044, rm_5045, rm_5046, rm_5048}** — 4 elements (all essential, none removable).

Union: [19:1426] ∪ [1681:2701] = 1407 + 1020 = 2427 bp → **RM_TOTAL_SV_COVERAGE = 88.3%**

### Comparison with current code

| Method | Selected elements | RM_TOTAL_SV_COVERAGE | RM_CLASSIFICATION |
|--------|-----------------|---------------------|------------------|
| Strict (non-overlapping) | rm_5045, rm_5048 | 63.5% | NON_REPETITIVE ✗ |
| **DP (minimum interval cover)** | **all four** | **88.3%** | **Retroposon ✓** |

The strict method rejects rm_5044 and rm_5046 because they overlap with rm_5045, even though each contributes a unique tail (487 bp and 196 bp respectively) that rm_5045 does not cover.

---

## Summary

| | Greedy + unique threshold | Minimum interval cover DP |
|---|---|---|
| **Optimal (minimum set)?** | Not always — order-dependent when sv_cov values are tied | Always |
| **Order-dependent?** | Yes | No |
| **Handles VNTR chains?** | Yes (same result as DP when no ties) | Yes |
| **Weight tie-break?** | N/A (order determines result) | Additive sv_cov sum; higher = better alignment quality |
| **Final coverage metric** | `calculate_total_coverage` (union) | `calculate_total_coverage` (union) |
| **Complexity** | O(n log n) | O(n²) per component |
| **New parameters** | Reuses `min_sv_coverage` | None |
