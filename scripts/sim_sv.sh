#!/bin/bash

# set -x
RED='\033[0;31m' ; GREEN='\033[0;32m' ; NC='\033[0m' # No Color
die() { echo -e "${RED}$1${NC}" >&2 ; echo ; exit 1 ; } # terminate script
info() {  echo -e "${GREEN}$1${NC}" >&2 ; }
info "$(date)"

SEED=42

OUTPUT_DIR=$(realpath "test/sim_ref_run")

BASE_REF=${OUTPUT_DIR}/base_ref/base_ref.fa
BOTH_REF=${OUTPUT_DIR}/both_ref.fa

SV_COUNT=100

READS="read_0.fasta"
SIM_READS="${OUTPUT_DIR}/sim_reads.fq.gz"
SAM="${OUTPUT_DIR}/mapped.sam"
BAM="${OUTPUT_DIR}/mapped.bam"
VCF="${OUTPUT_DIR}/sniffles.vcf"

SVSCANNER_VENV_PATH="/data/hiruna/SVscanner/svscanner"
SNIFFLES_VENV_PATH="/data/install/sniffles_260"
# VISOR_VENV_PATH="/data/hiruna/SVtoolkit/VISOR-1.1.2.1/visor"
VISOR_VENV_PATH="/data/hiruna/VISOR/visor_dev"
PACBIO_CCS="/data/install/ccs_v6.4.0/ccs"
BCFTOOLS="bcftools"
TABIX="tabix"
BGZIP="bgzip"

SIMULATE_SV="src/simulate_sv.py"
PBSIM3="/data/install/pbsim3-3.0.5/src/pbsim"
PBSIM3_DATA="/data/install/pbsim3-3.0.5/data"

SVSCANNER="scripts/run_workflow.sh"

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
    split=$1
    source "${SVSCANNER_VENV_PATH}/bin/activate"
    python ${SIMULATE_SV} --split ${split} --frac --simple -n ${SV_COUNT} --mob ${MOBILE_ELEMENTS} --rep ${REP_ELEMENTS} --out ${OUTPUT_DIR}/base_ref --seed ${SEED} || die "sim ref failed"
    deactivate
}

# b_ref, sv_ref, bed_file, ref_name
simulate_sv_using_visor() {
    base_ref=$1
    sv_ref=$2
    bed_file=$3
    ref_name=$4
    source "${VISOR_VENV_PATH}/bin/activate"
    rm -rf ${OUTPUT_DIR}/visor_hack 
    VISOR HACk --seed ${SEED} -b ${bed_file} -g ${base_ref} -o ${OUTPUT_DIR}/visor_hack || die "visor hack failed"
    deactivate
    awk '/^>/ {if (seq) print seq; print; seq=""; next} {seq=seq $0} END {if (seq) print seq}' ${OUTPUT_DIR}/visor_hack/h1.fa > ${sv_ref} || die "awk after visor failed"
    sed -i "1s/^>ref0$/>${ref_name}/" ${sv_ref} || die "renaming ref failed"
}

simulated_reads_pacbio_hifi() {
    ref_0=$1
    ref_1=$2
    
    cat ${ref_0} ${ref_1} > ${BOTH_REF}

    ${PBSIM3} --strategy wgs \
      --method qshmm \
      --qshmm ${PBSIM3_DATA}/QSHMM-RSII.model \
      --depth ${READ_DEPTH} \
      --genome ${BOTH_REF} \
      --seed ${SEED} \
      --pass-num ${PASS_NUM}} \
      --accuracy-mean 1 \
      --accuracy-min 1 \
      --difference-ratio 0:0:0 \
      --prefix ${OUTPUT_DIR}/simulated_reads || die "pbsim3 failed"

    samtools merge ${OUTPUT_DIR}/simulated_reads_*.bam -o ${OUTPUT_DIR}/merged_sim_reads.bam || die "samtools merge failed"
    rm ${OUTPUT_DIR}/simulated_reads_*.bam*

    ${PACBIO_CCS} ${OUTPUT_DIR}/merged_sim_reads.bam ${SIM_READS} || die "pacbio ccs failed"
    rm ${OUTPUT_DIR}/merged_sim_reads.bam
}

align_reads_pacbio_hifi() {
    base_ref=$1
    ref2=$2
    # minimap2 -cx map-hifi ${ref2} -t32 --secondary=no ${SIM_READS} | awk '{print $10/$11}' | datamash mean 1 sstdev 1 q1 1 median 1 q3 1 count 1 || die "minimap2 accuracy metric failed"
    # minimap2 -cx map-hifi ${base_ref} -t32 --secondary=no ${SIM_READS} | awk '{print $10/$11}' | datamash mean 1 sstdev 1 q1 1 median 1 q3 1 count 1 || die "minimap2 accuracy metric failed"
    minimap2 -ax map-hifi -Y ${base_ref} -t32 --secondary=no ${SIM_READS} > ${SAM} || die "minimap2 failed"
    samtools sort ${SAM} -o ${BAM} && samtools index ${BAM} || die "samtools sort and index failed"
    # samtools view -h -F 0x800 ${SAM} | sed -E 's/\tSA:Z:[^\t]+//g' | samtools sort -o ${BAM} && samtools index ${BAM} || die "samtools sort and index failed"
}

variant_call_sniffles() {
    base_ref=$1
    # Check if the virtual environment exists
    if [ ! -d "${SNIFFLES_VENV_PATH}" ]; then
        echo "Error: Virtual environment not found at ${SNIFFLES_VENV_PATH}"
        exit 1
    fi
    source "${SNIFFLES_VENV_PATH}/bin/activate"
    # sniffles --reference ${base_ref} --input ${BAM} --vcf ${VCF} --phase --minsvlen 50 --allow-overwrite  --no-qc || die "sniffles failed"
    sniffles --reference ${base_ref} --input ${BAM} --vcf ${VCF} --minsvlen 50 --allow-overwrite  --output-rnames || die "sniffles failed"
    deactivate
    rm -rf ${VCF}.gz && ${BGZIP} -c ${VCF} > ${VCF}.gz && ${TABIX} -p vcf ${VCF}.gz || die "bgzip and tabix failed"
}

run_svscanner() {
    base_ref=$1
    vcf=$2

    source "${SVSCANNER_VENV_PATH}/bin/activate"
    ${SVSCANNER} --out ${OUTPUT_DIR}/svscanner --vcf ${vcf} --ref ${base_ref} > ${OUTPUT_DIR}/svclass_stdout || die "sv classifier script failed"
    deactivate
}

create_output_dir
create_visor_bed 2
simulate_sv_using_visor ${BASE_REF} ${OUTPUT_DIR}/base_ref/ref_0.fa ${OUTPUT_DIR}/base_ref/visor_hack_0.bed ref_0
simulate_sv_using_visor ${BASE_REF} ${OUTPUT_DIR}/base_ref/ref_1.fa ${OUTPUT_DIR}/base_ref/visor_hack_1.bed ref_1
simulated_reads_pacbio_hifi ${OUTPUT_DIR}/base_ref/ref_0.fa ${OUTPUT_DIR}/base_ref/ref_1.fa
align_reads_pacbio_hifi ${BASE_REF}
variant_call_sniffles ${BASE_REF}
run_svscanner ${BASE_REF} ${OUTPUT_DIR}/sniffles.vcf.gz

info "success"

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

# simulated_reads_pacbio_hifi() {
#     ${PBSIM3} --strategy wgs \
#       --method qshmm \
#       --qshmm ${PBSIM3_DATA}/QSHMM-RSII.model \
#       --depth ${READ_DEPTH} \
#       --genome ${BASE_REF} \
#       --seed ${SEED} \
#       --pass-num ${PASS_NUM}} \
#       --accuracy-mean 1 \
#       --accuracy-min 1 \
#       --difference-ratio 0:0:0 \
#       --prefix ${OUTPUT_DIR}/simulated_reads || die "pbsim3 failed"

#     ${PACBIO_CCS} ${OUTPUT_DIR}/simulated_reads_0001.bam ${OUTPUT_DIR}/base_reads.fq.gz || die "pacbio ccs failed"
#     rm ${OUTPUT_DIR}/simulated_reads_0001.bam*

#     ${PBSIM3} --strategy wgs \
#       --method qshmm \
#       --qshmm ${PBSIM3_DATA}/QSHMM-RSII.model \
#       --depth ${READ_DEPTH} \
#       --genome ${SV_TREATED_REF} \
#       --seed ${SEED} \
#       --pass-num ${PASS_NUM}} \
#       --accuracy-mean 1 \
#       --accuracy-min 1 \
#       --difference-ratio 0:0:0 \
#       --prefix ${OUTPUT_DIR}/simulated_reads || die "pbsim3 failed"

#     ${PACBIO_CCS} ${OUTPUT_DIR}/simulated_reads_0001.bam ${OUTPUT_DIR}/sv_reads.fq.gz || die "pacbio ccs failed"

#     zcat ${OUTPUT_DIR}/base_reads.fq.gz ${OUTPUT_DIR}/sv_reads.fq.gz | gzip > ${SIM_READS} || die "zcat failed"
# }

# align_reads_ont() {
#     minimap2 -cx map-ont ${BASE_REF} -t32 --secondary=no ${SIM_READS} | awk '{print $10/$11}' | datamash mean 1 sstdev 1 q1 1 median 1 q3 1 count 1 || die "minimap2 accuracy metric failed"
#     minimap2 -ax map-ont ${BASE_REF} -t32 --secondary=no ${SIM_READS} > ${SAM} || die "minimap2 failed"
#     samtools sort ${SAM} -o ${BAM} && samtools index ${BAM} || die "samtools sort and index failed"
# }

# variant_call_bcftools() {
#     ${BCFTOOLS} mpileup -f ${BASE_REF} ${BAM} | ${BCFTOOLS} call -mv -Ov -o ${OUTPUT_DIR}/bcftools_mpileup.vcf
# }

# make_consensus() {
#     ${BCFTOOLS} view -e 'ALT ~ "<"' ${VCF}.gz -O b -o ${OUTPUT_DIR}/sym_filtered.vcf.gz
#     ${TABIX} ${OUTPUT_DIR}/sym_filtered.vcf.gz
#     ${BCFTOOLS} consensus -f ${BASE_REF} ${OUTPUT_DIR}/sym_filtered.vcf.gz -o ${OUTPUT_DIR}/with_variant_ref.fa 2>${OUTPUT_DIR}/bcftools.stderr || die "bcftools failed"
#     ${SAMTOOLS} faidx ${OUTPUT_DIR}/with_variant_ref.fa || die "samtools faidx failed"
# }
