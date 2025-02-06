#!/bin/bash

# set -x
RED='\033[0;31m' ; GREEN='\033[0;32m' ; NC='\033[0m' # No Color
die() { echo -e "${RED}$1${NC}" >&2 ; echo ; exit 1 ; } # terminate script
info() {  echo -e "${GREEN}$1${NC}" >&2 ; }
info "$(date)"

# Directories (change)
OUTPUT_DIR=$(realpath "SVtoolkit_output_1")
CLASSIFIER_DIR=$(realpath "SVclassifier")

# Input Files (change)
SAMPLE="HG002_subset"
REF_FASTA=$(realpath "/genome/hg38noAlt.fa")
SV_VCF=$(realpath "test/${SAMPLE}/${SAMPLE}.vcf.gz")
STR_BED=$(realpath "test/STRchive-disease-loci.bed")

# Repeat Masker and TRF programs (change if necessary)
TRF_BINARY=$(realpath "trf409.linux64")
REPEAT_MASKER=$(realpath "RepeatMasker/RepeatMasker")

# NTHREADS=$NSLOTS   # Total number of threads
NTHREADS=$(nproc --all)
# MAX_JOBS=8          # Max number of RepeatMasker process to run in parallel 
# THREADS_PER_JOB=$((NTHREADS / MAX_JOBS)) # Number of threads allocated to each RepeatMasker job (internal)

# Parameters (change if necessary)
NUM_SPLIT=100        # Number of sequences per file #todo: this should be calculated depending on the number of parallel jobs that can be run. also each split .fasta file should have similar size approx.
MIN_INTERSECT=0.05   # Min intersect between SV and repeat to be considered repetitive
MIN_COVERAGE=0.5     # Min coverage of an SV by repeat(s) to be considered repetitive
INTERVAL=0.05
DIAGRAM_LEN=100

# Python and bash scripts (keep as it is)
EXTRACT_SV_FLANKINGS=${CLASSIFIER_DIR}/extractSVs.py
TO_FASTA=${CLASSIFIER_DIR}/toFasta.sh
ANNOTATE=${CLASSIFIER_DIR}/repeatAnnotation.py
VISUALISE=${CLASSIFIER_DIR}/repeatDiagram.py
# Internal directories (keep as it is)
OUTPUT_SAMPLE_DIR=${OUTPUT_DIR}/${SAMPLE}
SPLIT_DIR=${OUTPUT_SAMPLE_DIR}/${SAMPLE}_${NUM_SPLIT}
SPLIT_RM=${SPLIT_DIR}/RMtmp
# File Intermediates (keep as it is)
SV_TAB=${OUTPUT_SAMPLE_DIR}/variant_flanking.tab
ID_FILE=${OUTPUT_SAMPLE_DIR}/${SAMPLE}_id.tab
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
    [ -z "$CLASSIFIER_DIR" ] && die "CLASSIFIER_DIR is not set"
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
    command -v bcftools >/dev/null 2>&1 || die "bcftools not found"
    command -v parallel >/dev/null 2>&1 || die "parallel not found"
}

create_output_dir() {
	test -d "${OUTPUT_DIR}" && rm -r "${OUTPUT_DIR}"
	mkdir "${OUTPUT_DIR}" || die "Failed creating ${OUTPUT_DIR}"
    mkdir -p ${SPLIT_DIR}
    mkdir -p ${SPLIT_RM}
}

extract_flanking_regions() {
    # 1) Extract sequence and flanking regions for variants
    info "1. Extracting structural variant sequences from VCF..."
    python3 ${EXTRACT_SV_FLANKINGS} -vcf ${SV_VCF} -fa ${REF_FASTA} -out ${SV_TAB} || die "failed"
    info "done"
}

create_fasta_files() {
    # 2) Create fasta sequences
    info "2. Converting structural variant sequences to FASTA..."
    split --numeric-suffixes=1 --suffix-length=3 -l ${NUM_SPLIT} "${SV_TAB}" ${SPLIT_DIR}/${SAMPLE} || die "split failed"
    chmod +r ${SPLIT_DIR}/*
    chmod +w ${SPLIT_DIR}/*
    find "${SPLIT_DIR}" -name "${SAMPLE}.[0-9]*" | parallel -j ${NTHREADS} "
        output_fasta=\"{}.fa\"
        id_file=\"{}_id.tab\"
        bash \"${TO_FASTA}\" {} \"\$output_fasta\" \"\$id_file\"" || die "failed"
    info "done"
}

run_trf() {
    # set -x
    # 3) Run Tandem Repeat Finder and RepeatMasker - wait for both to complete
    info "3. Running Tandem Repeat Finder..."
    find "${SPLIT_DIR}" -name "${SAMPLE}.*.fa" | parallel -j ${NTHREADS} --bar "${TRF_BINARY} {} 2 7 7 80 10 50 500 -h -ngs > {.}.dat" || die "failed"
    info "done"
}

run_repeatmasker() {
    info "3. Running RepeatMasker..."
    T2=$(date +%s)
    cd ${SPLIT_RM}
    # set -x
    # find "${SPLIT_DIR}" -name "${SAMPLE}.*.fa" | parallel -j "$MAX_JOBS" "RepeatMasker {} -pa $THREADS_PER_JOB -html -gff -dir '${SPLIT_DIR}'"
    # find "../" -name "${SAMPLE}.*.fa" | parallel -j "${MAX_JOBS}" "${REPEAT_MASKER} {} -pa ${THREADS_PER_JOB} -html -gff -dir '../'" || die "failed"
    # find "../" -name "${SAMPLE}.*.fa" | xargs -I {} taskset -c 0-9 "${REPEAT_MASKER}" {}  -html -gff -dir '../'
    # ls *.fasta | parallel --load 75% -j $(( $(nproc) / 4 )) repeatmasker -engine nhmmer -pa 2 {}

    # nhmmer search engine takes 2 cpus per job. hence -pa 2 means 2x2=4 cpus are used. Only runs new jobs if system load is below 75%.
    find "../" -name "${SAMPLE}.*.fa" | parallel --load 75% -j $(( $(nproc) / 4 )) ${REPEAT_MASKER} {} -pa 2 -html -gff -dir "../" || die "failed"
    wait 
    cd -
    T3=$(date +%s)
    RM_TIME=$((T3 - T2))
    info "done. RepeatMasker took ${RM_TIME} seconds"
}

process_repeatmasker_output() {
    # Process the RepeatMasker files (remove the header and 16th column)
    info "Process the RepeatMasker files (remove the header and 16th column)..."
    for rm_output in "${SPLIT_DIR}"/*.fa.out; do
        tail -n +4 "${rm_output}" | awk '{
            if ($16 == "*") $16 = "";
            else $16 = "";
            print $0
        }' OFS='\t' > "${rm_output}.tab" || die "failed"
    done
    info "done"
}

combine_split_files() {
    # # COMBINE SPLIT FILES
    cat ${SPLIT_DIR}/*_id.tab > ${ID_FILE}
    cat ${SPLIT_DIR}/*.out.tab > ${RM_FILE}
    cat ${SPLIT_DIR}/*.dat > ${TRF_FILE}
}

run_visualiser() {
    info "4. Running visualiser..."
    python3 ${VISUALISE} -sv ${ID_FILE} -trf ${TRF_FILE} -rm ${RM_FILE} -out ${VIS_OUTPUT} -length ${DIAGRAM_LEN} -min ${MIN_INTERSECT} || die "failed"
    info "done"
}

annotate_vcf_with_repeats() {
    info "5. Annotating SV VCF with repeat information..."
    python3 ${ANNOTATE} -id ${ID_FILE} -trf ${TRF_FILE} -rm ${RM_FILE} -vcf ${ANNOTATED_VCF} -trf_tsv ${TRF_TSV} -rm_tsv ${RM_TSV} -sv_vcf ${SV_VCF} -min ${MIN_INTERSECT} -mr ${MIN_COVERAGE} -div ${INTERVAL} -strchive ${STR_BED} || die "failed"
    info "done"
}

sort_and_index_vcf() {
    # Sort and Index the annotated VCF
    info "6. Sort and Index the annotated VCF..."
    bcftools sort -Oz -o ${ANNOTATED_VCF}.gz ${ANNOTATED_VCF} || die "bcftools sort failed"
    bcftools index -t ${ANNOTATED_VCF}.gz || die "bcftools index failed"
    info "done"
}

remove_intermediates() {
    ####################### Cleanup: Remove Intermediates ##########################
    rm -r ${SPLIT_DIR} || die "failed to remove ${SPLIT_DIR}"
    filesToRemove=("${SV_TAB}" "${ID_FILE}" "${RM_FILE}" "${TRF_FILE}")
    # Loop through and remove each file
    for file in "${filesToRemove[@]}"; do
        rm -f "${file}" || die "failed to remove ${file}"
    done
}

show_output_paths() {
    info "RepeatMasker output: ${RM_TSV}"
    info "TRF output: ${TRF_TSV}"
    info "Visualisation: ${VIS_OUTPUT}"
    info "SV VCF with repeats annotated: ${ANNOTATED_VCF}"
}

T0=$(date +%s)

check_required
create_output_dir
extract_flanking_regions
create_fasta_files
run_trf
run_repeatmasker
process_repeatmasker_output
combine_split_files
run_visualiser
annotate_vcf_with_repeats
sort_and_index_vcf
remove_intermediates
show_output_paths

T1=$(date +%s)
ELAPSED_TIME=$((T1 - T0))
info "The SVclassifier pipeline took ${ELAPSED_TIME} seconds"

info "$(date)"
info "Success!"
exit 0