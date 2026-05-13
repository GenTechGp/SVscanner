# Retrotransposition Detection

SVscanner optionally flags insertion SVs that are likely retrotransposon insertions — either *de novo* events or polymorphic insertions segregating in the population. This is implemented as a Level 3 post-processing step in `src/retrotransposition_detection.py`.

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

**Known limitation:** SVscanner cannot confirm TSDs or the poly-A tail position from existing VCF tags alone. Calls are therefore retrotransposition *candidates*, not definitive calls. Confirmation requires sequence-level inspection (e.g., IGV, or dedicated MEI callers such as TLDR).

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
       If all conditions met: record (RM_SUBFAMILY[i], confidence) as a hit

  4. If any hits found:
       RETROTP_ELEMENT   = comma-joined hit subfamilies  (parallel list)
       RETROTP_CANDIDATE = comma-joined confidence tiers (parallel list)
     Else: both tags absent
```

### Notes on list handling

`RM_SUBFAMILY` and `RM_SV_COVERAGE` are parallel lists — one entry per RepeatMasker hit. The algorithm evaluates each hit independently. No deduplication is performed: if RepeatMasker reports two hits of the same subfamily (e.g., fragmented annotation of a single Alu), both are reported in the output lists. This is consistent with the parallel-list convention of other SVscanner RM tags.

If `RM_SUBFAMILY` and `RM_SV_COVERAGE` have unequal lengths (should not occur in well-formed SVscanner output), the shorter list is padded — coverage defaults to 0.0, which fails the threshold.

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
```

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

- **Heavily 5'-truncated L1 insertions** shorter than `svlen_min` (100 bp) will not be called, even if the RM hit clearly identifies L1Hs.
- **Young Alu subfamilies not yet in the config** will be missed. The `config/retrotp_params.tsv` should be updated as new active subfamilies are characterised in the literature.
- **Low-quality SV calls** where RepeatMasker did not annotate the SV body (resulting in no `RM_SUBFAMILY` entries) will produce no call regardless of true content.
