#!/bin/bash

set -euo pipefail

# module load bcftools
# module load bedtools

# Inputs
VCF="sniffles.vcf.gz"
BED="base_ref/visor_hack.bed"
OUTVCF="sniffles.filtered.vcf.gz"
BUFFER=20  # Buffer in bp for insertion proximity

# Temporary files
TMP_DIR="tmp_dir"
DELS_BED="$TMP_DIR/sv_dels.bed"
INS_BED="$TMP_DIR/sv_ins.bed"
MATCHED_IDS="$TMP_DIR/keep_ids.txt"
BED_MATCHED="matched_regions.bed"

mkdir -p ${TMP_DIR}

# Step 1: Filter for FILTER=PASS variants only
bcftools view -f PASS "$VCF" -Oz -o "$TMP_DIR/pass.vcf.gz"
bcftools index -t "$TMP_DIR/pass.vcf.gz"

# Step 2: Split into INS and other SVs (DEL/DUP/etc)
bcftools view -i 'INFO/SVTYPE!="INS" && INFO/END!="."' "$TMP_DIR/pass.vcf.gz" | \
  bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\n' | \
  awk '$2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/' > "$DELS_BED"
bcftools view -i 'INFO/SVTYPE="INS"' "$TMP_DIR/pass.vcf.gz" | bcftools query -f '%CHROM\t%POS\t%POS\t%ID\n' > "$INS_BED"

# Step 3a: Intersect DEL/DUP with BED using reciprocal overlap
bedtools intersect -a "$DELS_BED" -b "$BED" -f 0.20 -r -wa -wb > "$TMP_DIR/intersected_dels.tsv"
cut -f4 "$TMP_DIR/intersected_dels.tsv" > "$TMP_DIR/del_ids.txt"

# Step 3b: Intersect INS using window proximity
bedtools window -a "$INS_BED" -b "$BED" -w "$BUFFER" > "$TMP_DIR/intersected_ins.tsv"
cut -f4 "$TMP_DIR/intersected_ins.tsv" > "$TMP_DIR/ins_ids.txt"

# Step 4: Combine IDs to keep
cat "$TMP_DIR/del_ids.txt" "$TMP_DIR/ins_ids.txt" | sort -u > "$MATCHED_IDS"

# Step 5: Filter original VCF
#bcftools view -i "ID=@$MATCHED_IDS" "$TMP_DIR/pass.vcf.gz" -Oz -o "$OUTVCF"
#bcftools index -t "$OUTVCF"

# Step 6: Extract matched BED regions
cut -f4-7 "$TMP_DIR/intersected_dels.tsv" > "$TMP_DIR/dels_regions.bed"
cut -f4-7 "$TMP_DIR/intersected_ins.tsv" > "$TMP_DIR/ins_regions.bed"
cat "$TMP_DIR/dels_regions.bed" "$TMP_DIR/ins_regions.bed" | sort -u > "$BED_MATCHED"

# Step 7: Get unmatched BED regions
bedtools intersect -v -a "$BED" -b <( cat "$BED_MATCHED" | cut -f2-4) > unmatched_regions.bed

# Step 8: Annotate TSV
python /home/hirsam/storage/SVtoolkit/src/annotate_matched_regions.py

# Step 9 (Alternative): Filter VCF by IDs retained in matched_annotated.tsv
cut -f1 matched_annotated.tsv | sort -u > "$TMP_DIR/final_ids.txt"

bcftools view -i "ID=@$TMP_DIR/final_ids.txt" "$TMP_DIR/pass.vcf.gz" -Oz -o "$OUTVCF"
bcftools index -t "$OUTVCF"

# Cleanup
# rm -r "$TMP_DIR"

# Output info
echo "Filtered VCF: $OUTVCF"
echo "Matched BED regions: $BED_MATCHED"
echo "Unmatched BED regions: unmatched_regions.bed"
