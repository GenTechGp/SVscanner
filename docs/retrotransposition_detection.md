# Retrotransposition Detection

SVscanner optionally flags insertion SVs that are likely retrotransposon insertions — either *de novo* events or polymorphic insertions segregating in the population. This is implemented as a Level 3 post-processing step in `src/retrotransposition_detection.py`.

---

## Design philosophy

RETROTP is a **high-specificity annotation flag**, not a MEI discovery tool. Its role is to add confident biological interpretation to variants already identified by an SV caller, based on RepeatMasker evidence present in the VCF.

Dedicated MEI callers such as [TLDR](https://github.com/adamewing/tldr) achieve higher sensitivity by detecting sequence-level signatures — target site duplications (TSDs), poly-A tails, split-read mapping — that SVscanner cannot access from VCF tags alone. The intended division of labour is:

| Tool | Role | Evidence base |
|---|---|---|
| Dedicated MEI caller (TLDR, PALMER, etc.) | MEI discovery | TSD, poly-A tail, sequence-level mapping |
| SVscanner RETROTP | MEI annotation / confirmation | RM subfamily + SV body coverage + SVLEN |

**Tag absent ≠ not an MEI.** Absence means SVscanner's RM-derived evidence did not meet the specificity threshold — not that the insertion is not a mobile element. For comprehensive MEI discovery, cross-validation with a dedicated MEI caller is recommended.

**Tag present = strong RM evidence.** When `RETROTP_CANDIDATE=HIGH` is set, the variant has satisfied three independent RM-derived criteria: the RM hit matches a known active subfamily by name, SV body coverage meets the per-class threshold, and SVLEN falls within the expected size range for that element class.

This design choice means RETROTP will not match the recall of a dedicated MEI caller. That is intentional.

---

## Biological background

### Mechanism

Active retrotransposons mobilise via a "copy-and-paste" mechanism known as target-primed reverse transcription (TPRT):

1. The source element is transcribed into RNA.
2. The RNA is reverse-transcribed into cDNA at the target site.
3. The cDNA is integrated into the genome by nicking the target DNA strand.

The result is a new copy of the element inserted at a new genomic location. Unlike DNA transposons, the source element is not excised — hence "copy-and-paste."

### Currently active elements in humans

Not all copies of a retrotransposon class are capable of transposition. Activity is restricted to young, recently diverged subfamilies. The three active non-LTR retrotransposon lineages in humans are:

**Alu (SINE)**
- Only AluY subfamilies are currently active. AluYa5 and AluYb8 are the most prolific.
- AluS subfamilies (~30 Mya) are largely extinct but some are still polymorphic in human populations.
- AluJ subfamilies (~65 Mya) are ancient and not active.
- Full-length Alu: ~282 bp. Commonly 5'-truncated in older insertions, but new insertions are typically near full-length.
- References: Konkel & Batzer (2010) *Genome Res*; Stewart et al. (2011) *Genome Res*

**L1 (LINE)**
- Only L1Hs (human-specific) is robustly active. Approximately 80–100 copies per genome are estimated to be retrotransposition-competent ("hot L1s").
- L1PA2 is also relatively recent and may contribute polymorphic insertions.
- L1PA3 and older (L1PA4+, L1M*) are progressively less active to completely extinct.
- Full-length L1Hs: ~6019 bp. The vast majority of L1 insertions are 5'-truncated; sizes range from ~100 bp to full-length.
- References: Beck et al. (2010) *Cell*; Brouha et al. (2003) *PNAS*

**SVA (Retroposon)**
- SVA is a composite element (SINE-R, VNTR, and Alu-like domains).
- SVA_F and SVA_E are currently active; SVA_D shows lower activity.
- Size range: ~700–4500 bp (average ~2.5 kb).
- References: Wang et al. (2005) *Genome Res*; Hancks & Kazazian (2010) *Hum Mutat*

### Molecular signatures of a genuine insertion event

| Signature | Description | Detectable from SVscanner VCF tags? |
|---|---|---|
| Young/active subfamily | Only specific subfamilies transpose | **Yes** — `RM_SUBFAMILY` |
| High TE coverage of insertion | Most of the SV body is the element | **Yes** — `RM_SV_COVERAGE` |
| SVLEN consistent with element class | Full-length or 5'-truncated size | **Yes** — `SVLEN` |
| SVTYPE = INS | Retrotransposition produces an insertion | **Yes** — `SVTYPE` |
| Target site duplications (TSDs) | 7–20 bp direct repeats flanking the insertion | **No** — sequence-level, not in current VCF tags |
| Poly-A tail | 3' poly-A from the RNA intermediate | **No** — positional sequence feature, not tagged |

**Known limitation:** SVscanner cannot confirm TSDs or the poly-A tail position from existing VCF tags alone. Calls are therefore retrotransposition *candidates*, not definitive calls. This is a deliberate consequence of the high-specificity design: RETROTP only fires on evidence derivable from RM annotations; signatures requiring sequence-level analysis are outside its scope (see [Design philosophy](#design-philosophy)).

---

## Detection logic

The detection script (`src/retrotransposition_detection.py`) applies the following per-record algorithm:

```
For each SV record in the annotated VCF:

  1. SVTYPE == "INS"
       → else: skip (retrotransposition produces insertions, not deletions)

  2. SVLEN must be present
       → abs(SVLEN) used throughout (some callers report negative SVLEN for INS)

  3. For each i in zip(RM_SUBFAMILY, RM_SV_COVERAGE):
       a. RM_SUBFAMILY[i] must be in the active subfamily config (retrotp_params.tsv)
       b. RM_SV_COVERAGE[i] >= config row's min_rm_sv_coverage
       c. abs(SVLEN) in [config row's svlen_min, config row's svlen_max]
       If all conditions met: emit (RM_SUBFAMILY[i], confidence) at position i
       Otherwise:            emit ('.', '.') at position i

  4. If at least one position emitted a non-'.' value:
       RETROTP_ELEMENT   = comma-joined values (same length as RM_SUBFAMILY)
       RETROTP_CANDIDATE = comma-joined values (same length as RM_SUBFAMILY)
     Else: both tags absent (record omitted from output TSV entirely)
```

### Notes on list handling

`RETROTP_CANDIDATE` and `RETROTP_ELEMENT` are **parallel to `RM_SUBFAMILY`** — they always have the same number of entries as `RM_SUBFAMILY`. Positions where the RM hit did not qualify (inactive subfamily, coverage below threshold, or SVLEN outside range) receive `'.'` as a placeholder. This is consistent with the parallel-list convention of all other `RM_*` tags in the SVscanner VCF, and allows direct positional comparison between `RETROTP_ELEMENT[i]` and `RM_SUBFAMILY[i]`.

The tag is absent entirely when no position in `RM_SUBFAMILY` qualifies — i.e., when all entries would be `'.'`. This preserves the "absent means no candidate" semantics at the record level.

No deduplication is performed: if RepeatMasker reports two hits of the same active subfamily (e.g., fragmented annotation of a single Alu), both positions are evaluated independently and both report their result.

If `RM_SUBFAMILY` and `RM_SV_COVERAGE` have unequal lengths (should not occur in well-formed SVscanner output), the shorter list is padded — coverage defaults to 0.0, which fails the threshold and emits `'.'`.

---

## Configuration

All thresholds are defined in `config/retrotp_params.tsv`. Each row specifies one active subfamily and its per-subfamily detection criteria. There is no runtime CLI override — tune the file directly.

| Column | Description |
|---|---|
| `subfamily` | RepeatMasker repeat name (matches `RM_SUBFAMILY` values) |
| `class` | RM class: `SINE`, `LINE`, or `Retroposon` |
| `confidence` | `HIGH` (well-established active subfamily) or `LOW` (borderline/older) |
| `min_rm_sv_coverage` | Minimum `RM_SV_COVERAGE` fraction for this hit to count |
| `svlen_min` | Minimum `abs(SVLEN)` in bp |
| `svlen_max` | Maximum `abs(SVLEN)` in bp |

Lower `min_rm_sv_coverage` thresholds are used for LINE/L1 (0.50) and SVA (0.60) compared to Alu (0.70), reflecting that heavily 5'-truncated L1 insertions and composite SVA elements are commonly underannotated by RepeatMasker relative to the full SV length.

---

## Output

The script produces a TSV (`retrotp_annotations.tsv`) compatible with `src/annotate_vcf.py`:

```
ID    RETROTP_CANDIDATE    RETROTP_ELEMENT
SV001    HIGH                AluYa5
SV002    HIGH,HIGH           AluYa5,AluYa5
SV003    HIGH,LOW            L1Hs,L1PA3
SV004    .,HIGH              .,AluYa5
```

In SV004, the first RM hit (e.g. an inactive L1 fragment) did not qualify; the second hit (AluYa5) did. Both positions are emitted, keeping the lists parallel to `RM_SUBFAMILY`.

Only records with at least one candidate hit are written. Absent records are not written (annotate_vcf.py leaves their tags absent).

### New VCF INFO tags

| Tag | Type | Description |
|---|---|---|
| `RETROTP_CANDIDATE` | Number=., String | Confidence tier per matched RM hit: `HIGH` or `LOW`. Parallel to `RETROTP_ELEMENT`. |
| `RETROTP_ELEMENT` | Number=., String | Active retrotransposon subfamily per matched RM hit. Parallel to `RETROTP_CANDIDATE`. |

---

## Usage

**Standalone:**
```bash
python3 src/retrotransposition_detection.py \
    --vcf annotated.vcf.gz \
    --config config/retrotp_params.tsv \
    --out retrotp_annotations.tsv

python3 src/annotate_vcf.py \
    --vcf annotated.vcf.gz \
    --header config/retrotp_header.txt \
    --tsv retrotp_annotations.tsv \
    --output annotated_retrotp.vcf
```

**Via pipeline:**
```bash
bash scripts/run_workflow.sh --retrotp --out OUT_DIR --vcf INPUT.vcf.gz --ref REF.fa
```

---

## Known false positive scenarios

- **Inherited TE at a breakpoint:** A pre-existing Alu in the reference flanking a deletion breakpoint can cause RepeatMasker to annotate the extracted SV+flank sequence with high Alu coverage. If the SV caller mislabels such a deletion as an insertion, SVscanner may produce `RM_SV_COVERAGE` values that pass the threshold. Always check SVTYPE carefully.
- **Highly repetitive loci:** Regions with tandem Alu arrays may produce INS calls where the insertion body contains multiple Alu elements, all annotated at high coverage. These are more likely to be assembly/calling artefacts than genuine retrotransposition events.
- **SVA composite annotation:** RepeatMasker may split an SVA into its component parts (VNTR, HERV-K-derived). Only the SVA-annotated hit (RM_SUBFAMILY = SVA_D/E/F) passes; component hits annotated as other classes do not.

## Known false negative scenarios

- **RepeatMasker reporting generic `AluY` instead of a specific active subfamily.** The RM library contains both the generic `AluY` consensus and specific subfamily consensuses (AluYa5, AluYb8, etc.). For many genuine AluY-clade insertions, RM scores the hit to the generic consensus rather than the specific subfamily, particularly when post-insertion mutations obscure the diagnostic positions. `AluY` is intentionally absent from `retrotp_params.tsv`: it is the ancestral consensus of the clade, not a characterised active lineage in its own right, and adding it would weaken the specificity guarantee of the `HIGH` tier. These cases are best detected by a dedicated MEI caller. Comparison with TLDR on a real long-read dataset shows this is the dominant source of false negatives for the ALU class.
- **SVA composite annotation by RepeatMasker.** RepeatMasker routinely fragments SVA insertions into their component parts (VNTR domain, SINE-R/HERV-K-derived region, Alu-like domain), none of which carry `SVA_D`, `SVA_E`, or `SVA_F` as the RM repeat name. RETROTP therefore produces no call for SVA insertions even when coverage and SVLEN would otherwise pass. SVA detection is effectively a limitation of RM annotation rather than of the detection logic itself. Dedicated MEI callers handle SVA more reliably.
- **Heavily 5'-truncated L1 insertions** shorter than `svlen_min` (100 bp) will not be called, even if the RM hit clearly identifies L1Hs.
- **Young Alu subfamilies not yet in the config** will be missed. The `config/retrotp_params.tsv` should be updated as new active subfamilies are characterised in the literature.
- **Low-quality SV calls** where RepeatMasker did not annotate the SV body (resulting in no `RM_SUBFAMILY` entries) will produce no call regardless of true content.

**Note — ancient subfamily insertions are not false negatives.** A full-length, high-coverage AluJ or AluS insertion with no RETROTP tag is a correct non-call, not a missed event. Ancient elements cannot transpose: their internal RNA Pol III promoters have accumulated too many post-insertion mutations to support efficient transcription, and they are not competitive substrates for L1 ORF2p machinery. Such a variant most likely represents a segregating polymorphic locus (the insertion is inherited and population-level polymorphic, absent from the reference but present in this individual) or an Alu-mediated NAHR product (a genomic rearrangement, not retrotransposition). The exclusion of ancient subfamilies from `retrotp_params.tsv` is therefore a biological decision, not only a specificity trade-off.
