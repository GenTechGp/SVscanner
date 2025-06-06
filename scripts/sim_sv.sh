#!/bin/bash

# set -x
RED='\033[0;31m' ; GREEN='\033[0;32m' ; NC='\033[0m' # No Color
die() { echo -e "${RED}$1${NC}" >&2 ; echo ; exit 1 ; } # terminate script
info() {  echo -e "${GREEN}$1${NC}" >&2 ; }
info "$(date)"

SEED=42

OUTPUT_DIR=$(realpath "test/sim_ref_run")

BASE_REF=${OUTPUT_DIR}/base_ref/base_ref.fa
SV_TREATED_REF=${OUTPUT_DIR}/sv_treated_ref.fa

SV_COUNT=1000

READS="read_0.fasta"
SIM_READS="${OUTPUT_DIR}/sim_reads.fq.gz"
SAM="${OUTPUT_DIR}/mapped.sam"
BAM="${OUTPUT_DIR}/mapped.bam"
VCF="${OUTPUT_DIR}/sniffles.vcf"

SVTOOLS_VENV_PATH="/data/hiruna/SVtoolkit/svtools"
SNIFFLES_VENV_PATH="/data/install/sniffles_260"
# VISOR_VENV_PATH="/data/hiruna/SVtoolkit/VISOR-1.1.2.1/visor"
VISOR_VENV_PATH="/data/hiruna/VISOR/visor_dev"
PACBIO_CCS="/data/install/ccs_v6.4.0/ccs"
BCFTOOLS="bcftools-1.21/bcftools"
SAMTOOLS="samtools-1.21/samtools"
TABIX="htslib-1.21/tabix"
BGZIP="htslib-1.21/bgzip"

SIM_REF="src/simulate_sv.py"
PBSIM3="/data/install/pbsim3-3.0.5/src/pbsim"
PBSIM3_DATA="/data/install/pbsim3-3.0.5/data"

SV_CLASSIFIER="scripts/run_classifier.sh"

TRF_BED="/genome/hg38.trf.bed"
MOBILE_ELEMENTS="test/databases/dfam_selected_species.fasta"
REP_ELEMENTS="test/databases/trcatalog_bins.bed"

READ_DEPTH=30
PASS_NUM=20

create_output_dir() {
	test -d "${OUTPUT_DIR}" && rm -r "${OUTPUT_DIR}"
    mkdir ${OUTPUT_DIR}
}

create_visor_bed() {
    source "${SVTOOLS_VENV_PATH}/bin/activate"
    python ${SIM_REF} --frac --simple -n ${SV_COUNT} --mob ${MOBILE_ELEMENTS} --rep ${REP_ELEMENTS} --out ${OUTPUT_DIR}/base_ref --seed ${SEED} || die "sim ref failed"
    deactivate
}

simulate_sv_using_visor() {
    source "${VISOR_VENV_PATH}/bin/activate"
    VISOR HACk --seed ${SEED} -b ${OUTPUT_DIR}/base_ref/visor_hack.bed -g ${BASE_REF} -o ${OUTPUT_DIR}/visor_hack || die "visor hack failed"
    awk '/^>/ {if (seq) print seq; print; seq=""; next} {seq=seq $0} END {if (seq) print seq}' ${OUTPUT_DIR}/visor_hack/h1.fa > ${SV_TREATED_REF} || die "awk after visor failed"
    deactivate
}

# sim_reads_ont() {
#     ${PBSIM3} --strategy wgs \
#       --method qshmm \
#       --qshmm ${PBSIM3_DATA}/QSHMM-ONT-HQ.model \
#       --depth ${READ_DEPTH} \
#       --genome ${SV_TREATED_REF} \
#       --seed ${SEED} \
#       --prefix ${OUTPUT_DIR}/simulated_reads || die "pbsim3 failed"

#     mv ${OUTPUT_DIR}/simualted_reads*.fastq.gz ${SIM_READS}
# }

simulated_reads_pacbio_hifi() {
    ${PBSIM3} --strategy wgs \
      --method qshmm \
      --qshmm ${PBSIM3_DATA}/QSHMM-RSII.model \
      --depth ${READ_DEPTH} \
      --genome ${BASE_REF} \
      --seed ${SEED} \
      --pass-num ${PASS_NUM}} \
      --accuracy-mean 1 \
      --accuracy-min 1 \
      --difference-ratio 0:0:0 \
      --prefix ${OUTPUT_DIR}/simulated_reads || die "pbsim3 failed"

    ${PACBIO_CCS} ${OUTPUT_DIR}/simulated_reads_0001.bam ${OUTPUT_DIR}/base_reads.fq.gz || die "pacbio ccs failed"
    rm ${OUTPUT_DIR}/simulated_reads_0001.bam*

    ${PBSIM3} --strategy wgs \
      --method qshmm \
      --qshmm ${PBSIM3_DATA}/QSHMM-RSII.model \
      --depth ${READ_DEPTH} \
      --genome ${SV_TREATED_REF} \
      --seed ${SEED} \
      --pass-num ${PASS_NUM}} \
      --accuracy-mean 1 \
      --accuracy-min 1 \
      --difference-ratio 0:0:0 \
      --prefix ${OUTPUT_DIR}/simulated_reads || die "pbsim3 failed"

    ${PACBIO_CCS} ${OUTPUT_DIR}/simulated_reads_0001.bam ${OUTPUT_DIR}/sv_reads.fq.gz || die "pacbio ccs failed"

    zcat ${OUTPUT_DIR}/base_reads.fq.gz ${OUTPUT_DIR}/sv_reads.fq.gz | gzip > ${SIM_READS} || die "zcat failed"
}

# align_reads_ont() {
#     minimap2 -cx map-ont ${BASE_REF} -t32 --secondary=no ${SIM_READS} | awk '{print $10/$11}' | datamash mean 1 sstdev 1 q1 1 median 1 q3 1 count 1 || die "minimap2 accuracy metric failed"
#     minimap2 -ax map-ont ${BASE_REF} -t32 --secondary=no ${SIM_READS} > ${SAM} || die "minimap2 failed"
#     samtools sort ${SAM} -o ${BAM} && samtools index ${BAM} || die "samtools sort and index failed"
# }

align_reads_pacbio_hifi() {
    minimap2 -cx map-hifi ${SV_TREATED_REF} -t32 --secondary=no ${SIM_READS} | awk '{print $10/$11}' | datamash mean 1 sstdev 1 q1 1 median 1 q3 1 count 1 || die "minimap2 accuracy metric failed"
    minimap2 -cx map-hifi ${BASE_REF} -t32 --secondary=no ${SIM_READS} | awk '{print $10/$11}' | datamash mean 1 sstdev 1 q1 1 median 1 q3 1 count 1 || die "minimap2 accuracy metric failed"
    minimap2 -ax map-hifi ${BASE_REF} -t32 --secondary=no ${SIM_READS} > ${SAM} || die "minimap2 failed"
    samtools sort ${SAM} -o ${BAM} && samtools index ${BAM} || die "samtools sort and index failed"
}

variant_call_sniffles() {
    # Check if the virtual environment exists
    if [ ! -d "${SNIFFLES_VENV_PATH}" ]; then
        echo "Error: Virtual environment not found at ${SNIFFLES_VENV_PATH}"
        exit 1
    fi
    source "${SNIFFLES_VENV_PATH}/bin/activate"
    # sniffles --reference ${BASE_REF} --input ${BAM} --vcf ${VCF} --phase --minsvlen 50 --allow-overwrite  --no-qc || die "sniffles failed"
    sniffles --reference ${BASE_REF} --input ${BAM} --vcf ${VCF} --minsvlen 50 --allow-overwrite  --no-qc || die "sniffles failed"
    deactivate
    rm -rf ${VCF}.gz && ${BGZIP} -c ${VCF} > ${VCF}.gz && ${TABIX} -p vcf ${VCF}.gz || die "bgzip and tabix failed"
}

# variant_call_bcftools() {
#     ${BCFTOOLS} mpileup -f ${BASE_REF} ${BAM} | ${BCFTOOLS} call -mv -Ov -o ${OUTPUT_DIR}/bcftools_mpileup.vcf
# }

# make_consensus() {
#     ${BCFTOOLS} view -e 'ALT ~ "<"' ${VCF}.gz -O b -o ${OUTPUT_DIR}/sym_filtered.vcf.gz
#     ${TABIX} ${OUTPUT_DIR}/sym_filtered.vcf.gz
#     ${BCFTOOLS} consensus -f ${BASE_REF} ${OUTPUT_DIR}/sym_filtered.vcf.gz -o ${OUTPUT_DIR}/with_variant_ref.fa 2>${OUTPUT_DIR}/bcftools.stderr || die "bcftools failed"
#     ${SAMTOOLS} faidx ${OUTPUT_DIR}/with_variant_ref.fa || die "samtools faidx failed"
# }

run_svclassifier() {
    source "${SVTOOLS_VENV_PATH}/bin/activate"
    ${SV_CLASSIFIER} --output_dir ${OUTPUT_DIR}/svclassifier --sample simulated --sv_vcf ${VCF}.gz --ref_fasta ${BASE_REF} > ${OUTPUT_DIR}/svclass_stdout || die "sv classifier script failed"
    deactivate
}

run_func() {
    local func_name="$1"
    if ! declare -f "$func_name" > /dev/null; then
        echo "Function '$func_name' not found"
        return 1
    fi

    local start_time=$(date +%s.%N)
    "$func_name"
    local end_time=$(date +%s.%N)

    local duration=$(echo "$end_time - $start_time" | bc)
    echo "Function:$func_name took ${duration} seconds"
}


run_func create_output_dir
run_func create_visor_bed
run_func simulate_sv_using_visor
run_func simulated_reads_pacbio_hifi
run_func align_reads_pacbio_hifi
run_func variant_call_sniffles
run_func run_svclassifier

info "success"