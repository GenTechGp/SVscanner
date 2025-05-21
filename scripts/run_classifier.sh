#!/bin/bash

# set -x
RED='\033[0;31m' ; GREEN='\033[0;32m' ; NC='\033[0m' # No Color
die() { echo -e "${RED}$1${NC}" >&2 ; echo ; exit 1 ; } # terminate script
info() {  echo -e "${GREEN}$1${NC}" >&2 ; }
info "$(date)"

# Input/Output (change)
#OUTPUT_DIR=""
#SAMPLE=""
# SV_VCF=$(realpath "test/${SAMPLE}/${SAMPLE}.vcf.gz")
#SV_VCF=""
echo "Path"
echo $PATH
echo "Pythonpath"
echo $PYTHONPATH
which python
#REF_FASTA=$(realpath "/g/data/te53/ontsv/references/hg38_reference_files/hg38.analysisSet.fa")
#REF_FASTA="/g/data/te53/variantcall/referenceresource/genome/pipeface/chm13XX.fasta"
#REF_FASTA=$(realpath "/genome/hg38.analysisSet.fa")
STR_BED=$(realpath "test/STRchive-disease-loci.bed")

##FOR SIMULATION AND TESTING
# OUTPUT_DIR=$(realpath "test/output_sim_ref")
# SAMPLE="sim_ref"
# REF_FASTA=$(realpath "test/sim_ref/base_ref/base_ref.fa")
# SV_VCF=$(realpath "test/sim_ref/sniffles.vcf.gz")

# Repeat Masker species (change)
# SPECIES="mammalia"
SPECIES="human"

# Repeat Masker and TRF, bcftools, bgzip, tabix  programs (change if necessary)
TRF_BINARY=$(realpath "trf409.linux64")
REPEAT_MASKER=$(realpath "RepeatMasker/RepeatMasker")

BCFTOOLS=$(realpath bcftools-1.21/bcftools)
BGZIP=$(realpath htslib-1.21/bgzip)
TABIX=$(realpath htslib-1.21/tabix)

# NTHREADS=$NSLOTS   # Total number of threads
NSPLIT_FILES=500
NTHREADS=$(nproc --all)
echo "Number of Threads"
echo ${NTHREADS}
MAX_JOBS=48          # Max number of RepeatMasker process to run in parallel 
THREADS_PER_JOB=$((NTHREADS / MAX_JOBS)) # Number of threads allocated to each RepeatMasker job (internal)

# Parameters (change if necessary)
MIN_SV_COVERAGE=0.05 #The minimum intersection between a repeat element and SV (aka sv_coverage) e.g. 0.05 (5%) (0 < min_sv_coverage < 1)
MIN_CLASS_SV_COVERAGE=0.25 #The minimum class sv coverage by repeat elements to be considered repetitive
MIN_TOTAL_SV_COVERAGE=0.75 #The minimum total sv coverage by repeat elements to be considered repetitive
INTERVAL=0.05
DIAGRAM_LEN=100

# Python and bash scripts (keep as it is)
EXTRACT_SV_FLANKINGS=$(realpath src/extract_sv.py)
ANNOTATION=$(realpath src/repeat_annotation.py)

# Function to show usage
usage() {
    echo "Usage: $0 --output_dir DIR --sample SAMPLE --sv_vcf FILE --ref_fasta FILE [options]"
    echo "Required arguments:"
    echo "  --output_dir DIR      Path to output directory"
    echo "  --sample SAMPLE       Sample name"
    echo "  --sv_vcf FILE         Path to SV VCF file"
    echo "  --ref_fasta FILE      Path to reference FASTA file"
    echo "Optional arguments:"
    echo "  --str_bed FILE        Path to STR BED file (default: $STR_BED)"
    echo "  --species NAME        Species name for RepeatMasker (default: $SPECIES)"
    echo "  --min_sv_coverage VAL   Minimum intersection between a repeat element and SV (aka sv_coverage) (default: $MIN_SV_COVERAGE)"
    echo "  --min_class_sv_coverage VAL Minimum class sv coverage by repeat elements to be considered repetitive (default: $MIN_CLASS_SV_COVERAGE)"
    echo "  --min_total_sv_coverage VAL Minimum total sv coverage by repeat elements to be considered repetitive (default: $MIN_TOTAL_SV_COVERAGE)"
    echo "  --interval VAL        Interval value (default: $INTERVAL)"
    echo "  --diagram_len VAL     Diagram length (default: $DIAGRAM_LEN)"
    echo "  --nsplit_files INT    Number of split files (default: $NSPLIT_FILES)"
    exit 1
}

# Function to parse arguments
parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --output_dir)
                OUTPUT_DIR=$(realpath "$2"); shift 2;;
            --sample)
                SAMPLE="$2"; shift 2;;
            --sv_vcf)
                SV_VCF=$(realpath "$2"); shift 2;;
            --ref_fasta)
                REF_FASTA=$(realpath "$2"); shift 2;;
            --str_bed)
                STR_BED=$(realpath "$2"); shift 2;;
            --species)
                SPECIES="$2"; shift 2;;
            --min_sv_coverage)
                MIN_SV_COVERAGE="$2"; shift 2;;
            --min_class_sv_coverage)
                MIN_CLASS_SV_COVERAGE="$2"; shift 2;;
            --min_total_sv_coverage)
                MIN_TOTAL_SV_COVERAGE="$2"; shift 2;;
            --interval)
                INTERVAL="$2"; shift 2;;
            --diagram_len)
                DIAGRAM_LEN="$2"; shift 2;;
	    --nsplit_files)
                NSPLIT_FILES="$2"; shift 2;;
            --help)
                usage;;
            *)
                echo "Unknown argument: $1"; usage;;
        esac
    done

    # Check required arguments
    if [[ -z "$OUTPUT_DIR" || -z "$SAMPLE" || -z "$SV_VCF" || -z "$REF_FASTA" ]]; then
        echo "Error: --output_dir, --sample, --sv_vcf and --ref_fasta are required."
        usage
    fi

    # Internal directories (keep as it is)
    EXTRACT_SV_FLANKS_OUT=${OUTPUT_DIR}/extract_sv_flanks_out
    ANNOTATIONS_OUT=${OUTPUT_DIR}/annotations_out
    RM_TMP=${OUTPUT_DIR}/RMtmp
    # File Intermediates (keep as it is)
    INFO_FILE=${OUTPUT_DIR}/${SAMPLE}_info.tab
    RM_FILE=${OUTPUT_DIR}/${SAMPLE}_rm.tab
    TRF_FILE=${OUTPUT_DIR}/${SAMPLE}_trf.tab

    # Final Outputs (change if necessary)
    VIS_OUTPUT=${OUTPUT_DIR}/${SAMPLE}_diagrams.txt
    ANNOTATED_VCF=${OUTPUT_DIR}/${SAMPLE}_annotated.vcf
    RM_TSV=${OUTPUT_DIR}/${SAMPLE}_annotatedRM.tsv
    TRF_TSV=${OUTPUT_DIR}/${SAMPLE}_annotatedTRF.tsv
}

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
	#test -d "${OUTPUT_DIR}" && rm -r "${OUTPUT_DIR}"
	mkdir -p "${OUTPUT_DIR}" || die "Failed creating ${OUTPUT_DIR}"
}

extract_flanking_regions() {
    # 1) Extract sequence and flanking regions for variants
    info "1. Extracting structural variant sequences from VCF..."
    python3 ${EXTRACT_SV_FLANKINGS} --vcf ${SV_VCF} --ref ${REF_FASTA} --out ${EXTRACT_SV_FLANKS_OUT} --min 10 -n ${NSPLIT_FILES} --info ${INFO_FILE} || die "failed"
    info "done"
}

run_trf() {
    # set -x
    # 3) Run Tandem Repeat Finder and RepeatMasker - wait for both to complete
    info "3. Running Tandem Repeat Finder..."
    T4=$(date +%s)
    find "${EXTRACT_SV_FLANKS_OUT}" -name "*.fa" | parallel -j ${NTHREADS} --bar "${TRF_BINARY} {} 2 7 7 80 10 50 500 -h -ngs > {.}.dat" || die "failed"
    T5=$(date +%s) || die "failed to get T3"
    TRF_TIME=$((T5 - T4)) || die "failed to calculate time"
    info "done. Tandem Repeat Finder took ${TRF_TIME} seconds"
}

run_repeatmasker() {
    mkdir -p ${RM_TMP}
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
    #find ${EXTRACT_SV_FLANKS_OUT} -name "*.fa" | xargs -I {} "${REPEAT_MASKER}" {} -pa $(( $(nproc) / 2 )) -html -gff -dir ${EXTRACT_SV_FLANKS_OUT} -species ${SPECIES}
    #find ${EXTRACT_SV_FLANKS_OUT} -name "*.fa" | parallel -j "${MAX_JOBS}" "${REPEAT_MASKER} {} -pa ${THREADS_PER_JOB} -html -gff -dir ${EXTRACT_SV_FLANKS_OUT} -species ${SPECIES}"
    # Check if any .fa files exist first
    shopt -s nullglob
    fa_files=(${EXTRACT_SV_FLANKS_OUT}/*.fa)
    (( ${#fa_files[@]} == 0 )) && die "No FASTA files found in ${EXTRACT_SV_FLANKS_OUT}"
    shopt -u nullglob

    # Run RepeatMasker in parallel with error sensitivity
    #find ${EXTRACT_SV_FLANKS_OUT} -name "*.fa" | parallel --halt now,fail=1 -j "${MAX_JOBS}" "${REPEAT_MASKER} {} -pa ${THREADS_PER_JOB} -html -gff -dir ${EXTRACT_SV_FLANKS_OUT} -species ${SPECIES}" || die "RepeatMasker failed"
    find ${EXTRACT_SV_FLANKS_OUT} -name "*.fa" | parallel --halt now,fail=1 -j "${MAX_JOBS}" "${REPEAT_MASKER} {} -pa ${THREADS_PER_JOB} -html -gff -dir ${EXTRACT_SV_FLANKS_OUT} -species ${SPECIES} > {}.log 2>&1" || die "RepeatMasker failed"

    cd - || die "cd - failed"
    T3=$(date +%s) || die "failed to get T3"
    RM_TIME=$((T3 - T2)) || die "failed to calculate time"
    rm -r ${RM_TMP} || die "failed to remove ${RM_TMP}"
    info "done. RepeatMasker took ${RM_TIME} seconds"
}

process_repeatmasker_output() {
    # Process the RepeatMasker files (remove the header and 16th column)
    info "Process the RepeatMasker files (remove the header and 16th column)..."
    for rm_output in "${EXTRACT_SV_FLANKS_OUT}"/*.fa.out; do
        tail -n +4 "${rm_output}" | awk '{
            $16 = "";
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
        --min_sv_coverage ${MIN_SV_COVERAGE}\
        --min_class_sv_coverage ${MIN_CLASS_SV_COVERAGE}\
        --min_total_sv_coverage ${MIN_TOTAL_SV_COVERAGE}\
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
    rm ${ANNOTATED_VCF}
    info "done"
}

show_output_paths() {
    info "Annotation outputs dir: ${ANNOTATIONS_OUT}"
    info "SV VCF with repeats annotated: ${ANNOTATED_VCF}"
}

T0=$(date +%s)

parse_args "$@"
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
