#!/bin/bash

# set -x
RED='\033[0;31m' ; GREEN='\033[0;32m' ; NC='\033[0m' # No Color
die() { echo -e "${RED}$1${NC}" >&2 ; echo ; exit 1 ; } # terminate script
info() {  echo -e "${GREEN}$1${NC}" >&2 ; }
info "$(date)"

BCFTOOLS=$(realpath bcftools-1.21/bcftools)
BGZIP=$(realpath htslib-1.21/bgzip)
TABIX=$(realpath htslib-1.21/tabix)

OUTPUT_DIR="test/out_longtr_fix"
ANNOT_TABLE="test/out_longtr_fix/annotation.tsv"
VCF_HEADER="test/out_longtr_fix/header.txt"
INPUT_VCF="test/LongTR/norm.vcf"

set -x 

COL_LIST=$(head -n 1 ${ANNOT_TABLE} | cut -c2- | tr '\t' ',')
${BGZIP} ${ANNOT_TABLE} -c > ${OUTPUT_DIR}/anno_table.gz || die "bgzip failed"
${TABIX} -s1 -b2 -e2 ${OUTPUT_DIR}/anno_table.gz || die "tabix failed"
${BCFTOOLS} annotate -a ${OUTPUT_DIR}/anno_table.gz -c ${COL_LIST} -h ${VCF_HEADER} ${INPUT_VCF} -o ${OUTPUT_DIR}/annotated.vcf || die "annotate failed"

info "success"