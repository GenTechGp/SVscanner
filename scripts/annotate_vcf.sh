#!/bin/bash

# set -x
RED='\033[0;31m' ; GREEN='\033[0;32m' ; NC='\033[0m' # No Color
die() { echo -e "${RED}$1${NC}" >&2 ; echo ; exit 1 ; } # terminate script
info() {  echo -e "${GREEN}$1${NC}" >&2 ; }
info "$(date)"

DIR="test/HG002_subset_mini"
BCFTOOLS=$(realpath bcftools-1.21/bcftools)
BGZIP=$(realpath htslib-1.21/bgzip)
TABIX=$(realpath htslib-1.21/tabix)

TABLE="test/anno_table.txt"

set -x 

COL_LIST=$(head -n 1 ${TABLE} | cut -c2- | tr '\t' ',')
${BGZIP} ${TABLE} -c > ${DIR}/anno_table.gz || die "bgzip failed"
${TABIX} -s1 -b2 -e2 ${DIR}/anno_table.gz || die "tabix failed"
${BCFTOOLS} annotate -a ${DIR}/anno_table.gz -c ${COL_LIST} -h ${DIR}/vcf_header.txt ${DIR}/HG002_subset_mini.vcf -o ${DIR}/anno_HG002_subset_mini.vcf || die "annotate failed"

info "success"