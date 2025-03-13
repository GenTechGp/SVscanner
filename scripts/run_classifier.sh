#!/bin/bash

# set -x
RED='\033[0;31m' ; GREEN='\033[0;32m' ; NC='\033[0m' # No Color
die() { echo -e "${RED}$1${NC}" >&2 ; echo ; exit 1 ; } # terminate script
info() {  echo -e "${GREEN}$1${NC}" >&2 ; }
info "$(date)"

# Directories (change)
OUTPUT_DIR=$(realpath "SVtoolkit_output_sim_ref")

# Input Files (change)
SAMPLE="sim_ref"
REF_FASTA=$(realpath "test/sim_ref/base_ref/base_ref.fa")
SV_VCF=$(realpath "test/sim_ref/sniffles.vcf.gz")
STR_BED=$(realpath "test/STRchive-disease-loci.bed")

# Repeat Masker and TRF, bcftools, bgzip, tabix  programs (change if necessary)
TRF_BINARY=$(realpath "trf409.linux64")
REPEAT_MASKER=$(realpath "RepeatMasker/RepeatMasker")

BCFTOOLS=$(realpath bcftools-1.21/bcftools)
BGZIP=$(realpath htslib-1.21/bgzip)
TABIX=$(realpath htslib-1.21/tabix)

# NTHREADS=$NSLOTS   # Total number of threads
NTHREADS=$(nproc --all)
# MAX_JOBS=8          # Max number of RepeatMasker process to run in parallel 
# THREADS_PER_JOB=$((NTHREADS / MAX_JOBS)) # Number of threads allocated to each RepeatMasker job (internal)

# Parameters (change if necessary)
MIN_INTERSECT=0.05   # Min intersect between SV and repeat to be considered repetitive
MIN_COVERAGE=0.5     # Min coverage of an SV by repeat(s) to be considered repetitive
INTERVAL=0.05
DIAGRAM_LEN=100

# Python and bash scripts (keep as it is)
EXTRACT_SV_FLANKINGS=src/extract_sv.py
ANNOTATION=src/repeat_annotation.py
# Internal directories (keep as it is)
OUTPUT_SAMPLE_DIR=${OUTPUT_DIR}/${SAMPLE}
EXTRACT_SV_FLANKS_OUT=${OUTPUT_SAMPLE_DIR}/extract_sv_flanks_out
ANNOTATIONS_OUT=${OUTPUT_SAMPLE_DIR}/annotations_out
RM_TMP=${OUTPUT_SAMPLE_DIR}/RMtmp
# File Intermediates (keep as it is)
INFO_FILE=${OUTPUT_SAMPLE_DIR}/${SAMPLE}_info.tab
RM_FILE=${OUTPUT_SAMPLE_DIR}/${SAMPLE}_rm.tab
TRF_FILE=${OUTPUT_SAMPLE_DIR}/${SAMPLE}_trf.tab

# Final Outputs (change if necessary)
VIS_OUTPUT=${OUTPUT_SAMPLE_DIR}/${SAMPLE}_diagrams.txt
ANNOTATED_VCF=${OUTPUT_SAMPLE_DIR}/${SAMPLE}_annotated.vcf
RM_TSV=${OUTPUT_SAMPLE_DIR}/${SAMPLE}_annotatedRM.tsv
TRF_TSV=${OUTPUT_SAMPLE_DIR}/${SAMPLE}_annotatedTRF.tsv

check_required() {
    [ -n "$VIRTUAL_ENV" ] && info "venv ($(basename "$VIRTUAL_ENV"))  found" || die "No venv found. Please activate the venv"
    [ -z "$OUTPUT_DIR" ] && die "OUTPUT_DIR is not set"
    [ -z "$REF_FASTA" ] && die "REF_FASTA is not set"
    [ -z "$SV_VCF" ] && die "SV_VCF is not set"
    [ -z "$STR_BED" ] && die "STR_BED is not set"

    info "Output dir: ${OUTPUT_DIR}"
    info "Reference: ${REF_FASTA}"
    info "Input SV VCF: ${SV_VCF}"
    info "Input BED: ${STR_BED}"

    info "NTHREADS: ${NTHREADS}"
    
    command -v split >/dev/null 2>&1 || die "split program not found"
    command -v ${TRF_BINARY} >/dev/null 2>&1 || die "TRF binary not found"
    command -v ${REPEAT_MASKER} >/dev/null 2>&1 || die "RepeatMasker not found"
    command -v ${BCFTOOLS} >/dev/null 2>&1 || die "${BCFTOOLS} not found"
    command -v ${BGZIP} >/dev/null 2>&1 || die "${BGZIP} not found"
    command -v ${TABIX} >/dev/null 2>&1 || die "${TABIX} not found"
    command -v parallel >/dev/null 2>&1 || die "parallel not found"
}

create_output_dir() {
	test -d "${OUTPUT_DIR}" && rm -r "${OUTPUT_DIR}"
	mkdir "${OUTPUT_DIR}" || die "Failed creating ${OUTPUT_DIR}"
    mkdir -p ${RM_TMP}
}

extract_flanking_regions() {
    # 1) Extract sequence and flanking regions for variants
    info "1. Extracting structural variant sequences from VCF..."
    python3 ${EXTRACT_SV_FLANKINGS} --vcf ${SV_VCF} --ref ${REF_FASTA} --out ${EXTRACT_SV_FLANKS_OUT} --min 10 -n 1 --info ${INFO_FILE} || die "failed"
    info "done"
}

run_trf() {
    # set -x
    # 3) Run Tandem Repeat Finder and RepeatMasker - wait for both to complete
    info "3. Running Tandem Repeat Finder..."
    find "${EXTRACT_SV_FLANKS_OUT}" -name "*.fa" | parallel -j ${NTHREADS} --bar "${TRF_BINARY} {} 2 7 7 80 10 50 500 -h -ngs > {.}.dat" || die "failed"
    info "done"
}

run_repeatmasker() {
    info "3. Running RepeatMasker..."
    T2=$(date +%s)
    cd ${RM_TMP}
    # set -x
    # find "${EXTRACT_SV_FLANKS_OUT}" -name "${SAMPLE}.*.fa" | parallel -j "$MAX_JOBS" "RepeatMasker {} -pa $THREADS_PER_JOB -html -gff -dir '${EXTRACT_SV_FLANKS_OUT}'"
    # find "../" -name "${SAMPLE}.*.fa" | parallel -j "${MAX_JOBS}" "${REPEAT_MASKER} {} -pa ${THREADS_PER_JOB} -html -gff -dir '../'" || die "failed"
    # find "../" -name "${SAMPLE}.*.fa" | xargs -I {} taskset -c 0-9 "${REPEAT_MASKER}" {}  -html -gff -dir '../'
    # ls *.fasta | parallel --load 75% -j $(( $(nproc) / 4 )) repeatmasker -engine nhmmer -pa 2 {}

    # nhmmer search engine takes 2 cpus per job. hence -pa 2 means 2x2=4 cpus are used. Only runs new jobs if system load is below 75%.
    # find "../" -name "${SAMPLE}.*.fa" | parallel --load 75% -j $(( $(nproc) / 4 )) ${REPEAT_MASKER} {} -pa 2 -html -gff -dir "../" || die "failed"
    # find ${EXTRACT_SV_FLANKS_OUT} -name "*.fa"
    find ${EXTRACT_SV_FLANKS_OUT} -name "*.fa" | xargs -I {} "${REPEAT_MASKER}" {} -pa 8 -html -gff -dir ${EXTRACT_SV_FLANKS_OUT}
    cd -
    T3=$(date +%s)
    RM_TIME=$((T3 - T2))
    rm -r ${RM_TMP} || die "failed to remove ${RM_TMP}"
    info "done. RepeatMasker took ${RM_TIME} seconds"
}

process_repeatmasker_output() {
    # Process the RepeatMasker files (remove the header and 16th column)
    info "Process the RepeatMasker files (remove the header and 16th column)..."
    for rm_output in "${EXTRACT_SV_FLANKS_OUT}"/*.fa.out; do
        tail -n +4 "${rm_output}" | awk '{
            if ($16 == "*") $16 = "";
            else $16 = "";
            print $0
        }' OFS='\t' > "${rm_output}.tab" || die "failed"
    done
    # rm -rf "${EXTRACT_SV_FLANKS_OUT}/*.fa.out"
    info "done"
}

combine_split_files() {
    cat ${EXTRACT_SV_FLANKS_OUT}/*.dat > ${TRF_FILE}
    cat ${EXTRACT_SV_FLANKS_OUT}/*.out.tab > ${RM_FILE}
    # rm -rf ${EXTRACT_SV_FLANKS_OUT}/*.dat
    # rm -rf ${EXTRACT_SV_FLANKS_OUT}/*.out.tab
}

annotation() {
    info "4. Annotating..."
    test -d "${ANNOTATIONS_OUT}" && rm -r "${ANNOTATIONS_OUT}"
    python3 ${ANNOTATION} \
        --vcf ${SV_VCF}\
        --rm ${RM_FILE}\
        --trf ${TRF_FILE}\
        --info ${INFO_FILE}\
        --str ${STR_BED}\
        --out ${ANNOTATIONS_OUT}\
        --minsec ${MIN_INTERSECT}\
        --minrep ${MIN_COVERAGE}\
        --div ${INTERVAL}\
        --len ${DIAGRAM_LEN} || die "failed"
    
    COL_LIST=$(head -n 1 ${ANNOTATIONS_OUT}/vcf_annotate.tsv | cut -c2- | tr '\t' ',')
    ${BGZIP} ${ANNOTATIONS_OUT}/vcf_annotate.tsv -c > ${ANNOTATIONS_OUT}/vcf_annotate.gz || die "${BGZIP} failed"
    ${TABIX} -s1 -b2 -e2 ${ANNOTATIONS_OUT}/vcf_annotate.gz || die "${TABIX} failed"
    ${BCFTOOLS} annotate -a ${ANNOTATIONS_OUT}/vcf_annotate.gz -c ${COL_LIST} -h ${ANNOTATIONS_OUT}/vcf_header.txt ${SV_VCF} -o ${ANNOTATED_VCF} || die "${BCFTOOLS} annotate failed"

    info "done"

}

sort_and_index_vcf() {
    # Sort and Index the annotated VCF
    info "6. Sort and Index the annotated VCF..."
    ${BCFTOOLS} sort -Oz -o ${ANNOTATED_VCF}.gz ${ANNOTATED_VCF} || die "${BCFTOOLS} sort failed"
    ${BCFTOOLS} index -t ${ANNOTATED_VCF}.gz || die "${BCFTOOLS} index failed"
    
    info "done"
}

show_output_paths() {
    info "Annotation outputs dir: ${ANNOTATIONS_OUT}"
    info "SV VCF with repeats annotated: ${ANNOTATED_VCF}"
}

T0=$(date +%s)

check_required
create_output_dir
extract_flanking_regions
run_trf
run_repeatmasker
process_repeatmasker_output
combine_split_files
annotation
sort_and_index_vcf
show_output_paths

T1=$(date +%s)
ELAPSED_TIME=$((T1 - T0))
info "The SVclassifier pipeline took ${ELAPSED_TIME} seconds"

info "$(date)"
info "Success!"
exit 0