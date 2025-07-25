#!/bin/bash

VERSION="SVclassifier v0.1.0"

# set -x
die() { echo -e "$1" >&2 ; echo ; exit 1 ; } # terminate script

# Input/Output (change)
#REF=$(realpath "/g/data/te53/ontsv/references/hg38_reference_files/hg38.analysisSet.fa")
#REF="/g/data/te53/variantcall/referenceresource/genome/pipeface/chm13XX.fasta"
#REF=$(realpath "/genome/hg38.analysisSet.fa")
STR_BED=$(realpath "test/databases/STRchive-disease-loci.bed")

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

MAX_JOBS=48          # Max number of RepeatMasker process to run in parallel 
THREADS_PER_JOB=$((NTHREADS / MAX_JOBS)) # Number of threads allocated to each RepeatMasker job (internal)

# Parameters (change if necessary)
MIN_SV_COVERAGE=0.05 #The minimum intersection between a repeat element and SV (aka sv_coverage) e.g. 0.05 (5%) (0 < min_sv_coverage < 1)
MIN_CLASS_SV_COVERAGE=0.25 #The minimum class sv coverage by repeat elements to be considered repetitive
MIN_TOTAL_SV_COVERAGE=0.75 #The minimum total sv coverage by repeat elements to be considered repetitive
MAX_TRF_OVERLAP=0.1 #The maximum TRF element overlap to be considered non-overlapping (0 < max_trf_overlap < 1)
INTERVAL=0.05
DIAGRAM_LEN=100

FLAG_DELETE_TMP_FILES=1 # Flag to delete temporary files (1 = yes, 0 = no)
PREFIX="" # Prefix for output files (default: None)

# Python and bash scripts (keep as it is)
EXTRACT_SV_FLANKINGS=$(realpath src/extract_sv.py)
ANNOTATION=$(realpath src/repeat_annotation.py)
PLOT=$(realpath src/generate_plots.py)

# Function to show usage
usage() {
    echo "Usage: $0 --out DIR --vcf FILE --ref FILE [options]"
    echo "Required arguments:"
    echo "  --out DIR      Path to output directory"
    echo "  --vcf FILE     Path to SV VCF file"
    echo "  --ref FILE     Path to reference FASTA file"
    echo "Optional arguments:"
    echo "  --prefix NAME           Prefix for output files (default: None)"
    echo "  --str_bed FILE          Path to STR BED file (default: $STR_BED)"
    echo "  --species NAME          Species name for RepeatMasker (default: $SPECIES)"
    echo "  --min_sv_coverage VAL   Minimum intersection between a repeat element and SV (aka sv_coverage) (default: $MIN_SV_COVERAGE)"
    echo "  --min_class_sv_coverage VAL Minimum class sv coverage by repeat elements to be considered repetitive (default: $MIN_CLASS_SV_COVERAGE)"
    echo "  --min_total_sv_coverage VAL Minimum total sv coverage by repeat elements to be considered repetitive (default: $MIN_TOTAL_SV_COVERAGE)"
    echo "  --max_trf_overlap VAL   Maximum TRF element overlap to be considered non-overlapping (default: $MAX_TRF_OVERLAP)"
    echo "  --interval VAL          Interval value (default: $INTERVAL)"
    echo "  --diagram_len VAL       Diagram length (default: $DIAGRAM_LEN)"
    echo "  --nsplit_files INT      Number of split files (default: $NSPLIT_FILES)"
    echo "  --keep_tmp_files        Keep temporary files (default: delete)"
    echo "  --help                  Show this help message"
    echo "  --version               Show version information"
    echo "Version: $VERSION"
    exit 0
}

# Function to parse arguments
parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --out)
                OUTPUT_DIR=$(realpath "$2"); shift 2;;
            --vcf)
                VCF=$(realpath "$2"); shift 2;;
            --ref)
                REF=$(realpath "$2"); shift 2;;
            --prefix)
                PREFIX="$2_"; shift 2;;
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
            --max_trf_overlap)
                MAX_TRF_OVERLAP="$2"; shift 2;;
            --interval)
                INTERVAL="$2"; shift 2;;
            --diagram_len)
                DIAGRAM_LEN="$2"; shift 2;;
            --nsplit_files)
                NSPLIT_FILES="$2"; shift 2;;
            --keep_tmp_files)
                FLAG_DELETE_TMP_FILES=0; shift;;
            --help)
                usage;;
            --version)
                echo "$VERSION"; exit 0;;
            *)
                echo "Unknown argument: $1"; usage;;
        esac
    done

    # Check required arguments
    if [[ -z "$OUTPUT_DIR" || -z "$VCF" || -z "$REF" ]]; then
        echo "Error: --out, --vcf and --ref are required."
        usage
    fi

    # Internal directories (keep as it is)
    EXTRACT_SV_FLANKS_OUT=${OUTPUT_DIR}/${PREFIX}extract_sv_flanks_out
    ANNOTATIONS_OUT=${OUTPUT_DIR}/${PREFIX}annotations_out
    RM_TMP=${OUTPUT_DIR}/RMtmp
    # File Intermediates (keep as it is)
    INFO_FILE=${OUTPUT_DIR}/${PREFIX}info.tab
    RM_FILE=${OUTPUT_DIR}/${PREFIX}rm.tab
    TRF_FILE=${OUTPUT_DIR}/${PREFIX}trf.tab

    # Final Outputs (change if necessary)
    ANNOTATED_VCF=${OUTPUT_DIR}/${PREFIX}annotated.vcf
}

check_required() {

    [ -n "$VIRTUAL_ENV" ] && echo "venv ($(basename "$VIRTUAL_ENV"))  found" || die "No venv found. Please activate the venv"
    [ -z "$OUTPUT_DIR" ] && die "OUTPUT_DIR is not set"
    [ -z "$REF" ] && die "REF is not set"
    [ -z "$VCF" ] && die "VCF is not set"
    [ -z "$STR_BED" ] && die "STR_BED is not set"

    echo "Output dir: ${OUTPUT_DIR}"
    echo "Reference: ${REF}"
    echo "Input SV VCF: ${VCF}"
    echo "Input BED: ${STR_BED}"

    echo "Number of Threads: ${NTHREADS}"
    
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
    # error if output directory already exists
    if [[ -d "${OUTPUT_DIR}" ]]; then
        die "Output directory ${OUTPUT_DIR} already exists. Please choose a different output directory or delete the existing one."
    fi
	mkdir -p "${OUTPUT_DIR}" || die "Failed creating ${OUTPUT_DIR}"
}

extract_flanking_regions() {
    # 1) Extract sequence and flanking regions for variants
    echo "1. Extracting structural variant sequences from VCF..."
    python3 ${EXTRACT_SV_FLANKINGS} --vcf ${VCF} --ref ${REF} --out ${EXTRACT_SV_FLANKS_OUT} --min 10 -n ${NSPLIT_FILES} --info ${INFO_FILE} || die "failed"
    echo "done"
}

run_trf() {
    # set -x
    # 3) Run Tandem Repeat Finder and RepeatMasker - wait for both to complete
    echo "3. Running Tandem Repeat Finder..."
    T4=$(date +%s)
    find "${EXTRACT_SV_FLANKS_OUT}" -name "*.fa" | parallel -j ${NTHREADS} --bar "${TRF_BINARY} {} 2 7 7 80 10 50 500 -h -ngs > {.}.dat" || die "failed"
    T5=$(date +%s) || die "failed to get T3"
    TRF_TIME=$((T5 - T4)) || die "failed to calculate time"
    echo "done. Tandem Repeat Finder took ${TRF_TIME} seconds"
}

run_repeatmasker() {
    mkdir -p ${RM_TMP}
    echo "3. Running RepeatMasker..."
    T2=$(date +%s)
    cd ${RM_TMP}

    shopt -s nullglob
    fa_files=(${EXTRACT_SV_FLANKS_OUT}/*.fa)
    (( ${#fa_files[@]} == 0 )) && die "No FASTA files found in ${EXTRACT_SV_FLANKS_OUT}"
    shopt -u nullglob

    # Run RepeatMasker in parallel with error sensitivity
    find ${EXTRACT_SV_FLANKS_OUT} -name "*.fa" | parallel --halt now,fail=1 -j "${MAX_JOBS}" "${REPEAT_MASKER} {} -pa ${THREADS_PER_JOB} -html -gff -dir ${EXTRACT_SV_FLANKS_OUT} -species ${SPECIES} > {}.log 2>&1" || die "RepeatMasker failed"

    cd - || die "cd - failed"
    T3=$(date +%s) || die "failed to get T3"
    RM_TIME=$((T3 - T2)) || die "failed to calculate time"
    rm -r ${RM_TMP} || die "failed to remove ${RM_TMP}"
    echo "done. RepeatMasker took ${RM_TIME} seconds"
}

process_repeatmasker_output() {
    # Process the RepeatMasker files (remove the header and 16th column)
    echo "Process the RepeatMasker files (remove the header and 16th column)..."
    for rm_output in "${EXTRACT_SV_FLANKS_OUT}"/*.fa.out; do
        tail -n +4 "${rm_output}" | awk '{
            $16 = "";
            print $0
        }' OFS='\t' > "${rm_output}.tab" || die "failed"
    done
    # rm -rf "${EXTRACT_SV_FLANKS_OUT}/*.fa.out"
    echo "done"
}

combine_split_files() {
    cat ${EXTRACT_SV_FLANKS_OUT}/*.dat > ${TRF_FILE}
    cat ${EXTRACT_SV_FLANKS_OUT}/*.out.tab > ${RM_FILE}
    # rm -rf ${EXTRACT_SV_FLANKS_OUT}/*.dat
    # rm -rf ${EXTRACT_SV_FLANKS_OUT}/*.out.tab
}

annotation() {
    echo "4. Annotating..."
    test -d "${ANNOTATIONS_OUT}" && rm -r "${ANNOTATIONS_OUT}"
    python3 ${ANNOTATION} \
        --vcf ${VCF}\
        --rm ${RM_FILE}\
        --trf ${TRF_FILE}\
        --info ${INFO_FILE}\
        --str ${STR_BED}\
        --out ${ANNOTATIONS_OUT}\
        --min_sv_coverage ${MIN_SV_COVERAGE}\
        --min_class_sv_coverage ${MIN_CLASS_SV_COVERAGE}\
        --min_total_sv_coverage ${MIN_TOTAL_SV_COVERAGE}\
        --max_trf_overlap ${MAX_TRF_OVERLAP}\
        --div ${INTERVAL}\
        --len ${DIAGRAM_LEN} || die "failed"
    
    COL_LIST=$(head -n 1 ${ANNOTATIONS_OUT}/vcf_annotate.tsv | cut -c2- | tr '\t' ',')
    ${BGZIP} ${ANNOTATIONS_OUT}/vcf_annotate.tsv -c > ${ANNOTATIONS_OUT}/vcf_annotate.gz || die "${BGZIP} failed"
    ${TABIX} -s1 -b2 -e2 ${ANNOTATIONS_OUT}/vcf_annotate.gz || die "${TABIX} failed"
    ${BCFTOOLS} annotate -a ${ANNOTATIONS_OUT}/vcf_annotate.gz -c ${COL_LIST} -h ${ANNOTATIONS_OUT}/vcf_header.txt ${VCF} -o ${ANNOTATED_VCF} || die "${BCFTOOLS} annotate failed"

    echo "done"

}

sort_and_index_vcf() {
    # Sort and Index the annotated VCF
    echo "6. Sort and Index the annotated VCF..."
    ${BCFTOOLS} sort -Oz -o ${ANNOTATED_VCF}.gz ${ANNOTATED_VCF} || die "${BCFTOOLS} sort failed"
    ${BCFTOOLS} index -t ${ANNOTATED_VCF}.gz || die "${BCFTOOLS} index failed"
    rm ${ANNOTATED_VCF}
    echo "done"
}

plot_classifications() {
    echo "4. Plotting..."
    python3 ${PLOT} \
        --out ${ANNOTATIONS_OUT}/plots\
        --tsv ${ANNOTATIONS_OUT}/plot_annotate.tsv || die "failed"
    echo "done"
}

show_output_paths() {

    if [[ ${FLAG_DELETE_TMP_FILES} -eq 1 ]]; then
        echo "Deleting temporary files..."
        mv ${ANNOTATIONS_OUT}/plots/plots.pdf ${OUTPUT_DIR}/${PREFIX}plots.pdf || die "failed to move plots to ${OUTPUT_DIR}"
        mv ${ANNOTATIONS_OUT}/diagram.txt ${OUTPUT_DIR}/${PREFIX}diagram.txt || die "failed to move diagram.txt to ${OUTPUT_DIR}"
        mv ${ANNOTATIONS_OUT}/rm_diagram.tsv ${OUTPUT_DIR}/${PREFIX}rm_diagram.tsv || die "failed to move rm_diagram.tsv to ${OUTPUT_DIR}"
        mv ${ANNOTATIONS_OUT}/trf_diagram.tsv ${OUTPUT_DIR}/${PREFIX}trf_diagram.tsv || die "failed to move trf_diagram.tsv to ${OUTPUT_DIR}"

        rm -rf ${ANNOTATIONS_OUT} || die "failed to remove ${ANNOTATIONS_OUT}"
        rm -rf ${EXTRACT_SV_FLANKS_OUT} || die "failed to remove ${EXTRACT_SV_FLANKS_OUT}"
        # rm -f ${INFO_FILE} || die "failed to remove ${INFO_FILE}"
        # rm -f ${RM_FILE} || die "failed to remove ${RM_FILE}"
        # rm -f ${TRF_FILE} || die "failed to remove ${TRF_FILE}"
        echo "Annotation outputs dir: ${OUTPUT_DIR}"
        echo "Plots dir: ${OUTPUT_DIR}/plots"
        echo "SV VCF with repeats annotated: ${ANNOTATED_VCF}.gz"
    else
        echo "Annotation outputs dir: ${ANNOTATIONS_OUT}"
        echo "Plots dir: ${ANNOTATIONS_OUT}/plots"
        echo "SV VCF with repeats annotated: ${ANNOTATED_VCF}"
    fi

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
plot_classifications
show_output_paths

T1=$(date +%s)
ELAPSED_TIME=$((T1 - T0))
echo "The SVclassifier pipeline took ${ELAPSED_TIME} seconds"

echo "$(date)"
echo "Success!"
exit 0
