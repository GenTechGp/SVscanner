#!/bin/bash

# set -x
RED='\033[0;31m' ; GREEN='\033[0;32m' ; NC='\033[0m' # No Color
die() { echo -e "${RED}$1${NC}" >&2 ; echo ; exit 1 ; } # terminate script
info() {  echo -e "${GREEN}$1${NC}" >&2 ; }
info "$(date)"

OUTPUT_DIR=$(realpath "test/custom_vcf")

REF="genome_B.fa"
READS="read_0.fasta"
SAM="mapped.sam"
BAM="mapped.bam"
# VCF="sniffles_207.vcf"
# VENV_PATH="/data/hiruna/sniffles2/sniffles_207"
VCF="sniffles_253.vcf"
VENV_PATH="/data/hiruna/sniffles2/sniffles_253"

cd ${OUTPUT_DIR}

align() {
    minimap2 -ax map-ont ${REF} -t32 --secondary=no ${READS} > ${SAM} || die "minimap2 failed"
    samtools sort ${SAM} -o ${BAM} && samtools index ${BAM} || die "samtools sort and index failed"
}

variant_call() {

    # Check if the virtual environment exists
    if [ ! -d "$VENV_PATH" ]; then
        echo "Error: Virtual environment not found at $VENV_PATH"
        exit 1
    fi

    # Activate the virtual environment
    source "$VENV_PATH/bin/activate"

    sniffles --reference ${REF} --input ${BAM} --vcf ${VCF} --allow-overwrite --min-alignment-length 50 --minsupport 1 --minsvlen 2 --no-qc || die "sniffles failed"

    # Deactivate the virtual environment
    deactivate
    rm -rf ${VCF}.gz && bgzip -c ${VCF} > ${VCF}.gz && tabix -p vcf ${VCF}.gz || die "bgzip tabix failed"
}

make_consensus() {
    bcftools consensus -f ${REF} ${VCF}.gz -o with_variant_${REF} || die "bcftools failed"
    samtools faidx with_variant_${REF} || die "samtools faidx failed"
}

align
variant_call
# make_consensus
