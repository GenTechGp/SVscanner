#!/bin/bash

# set -x
die() { echo -e "$1" >&2 ; echo ; exit 1 ; } # terminate script

# Repository root, resolved from this script's own location so the workflow can be
# invoked from any working directory (e.g. via a module that puts it on $PATH).
SVSCANNER_HOME=$(cd -- "$(dirname -- "$(realpath "${BASH_SOURCE[0]}")")/.." && pwd) || die "could not resolve SVscanner root"

# The VERSION file at the repository root is the single source of truth for the
# release number - never hardcode it here. scripts/check_version.sh enforces that
# this script, the README and the git tag all agree with that file.
VERSION="SVscanner v$(cat "${SVSCANNER_HOME}/VERSION" 2>/dev/null || echo unknown)"

# Input/Output (change)
#REF=$(realpath "/g/data/te53/ontsv/references/hg38_reference_files/hg38.analysisSet.fa")
#REF="/g/data/te53/variantcall/referenceresource/genome/pipeface/chm13XX.fasta"
#REF=$(realpath "/genome/hg38.analysisSet.fa")
STR_BED="${SVSCANNER_HOME}/test/databases/STRchive-disease-loci.bed"

# Repeat Masker species (change)
# SPECIES="mammalia"
SPECIES="human"

# RepeatMasker search engine (-e), e.g. rmblast, hmmer. Empty means RepeatMasker uses
# its own configured default - RMBlast in the container (see docs/docker.md).
ENGINE=""

TRF_BINARY=""
REPEAT_MASKER=""
BCFTOOLS=""
BGZIP=""

# Oldest RepeatMasker supported. 4.2.6 fixes forward-strand SVAs being split into
# fragments, and adds the FAMDB_DATA_DIR and REPEATMASKER_CACHE_DIR settings used below.
RM_MIN_VERSION="4.2.6"
RM_VERSION=""

# Dfam database (FamDB 3.0, Dfam 4.0+), see setup_dfam_library. At most one of:
#   --dfam_db    a complete database (base file dfamNN.0.h5 + partitions), used
#                instead of the one RepeatMasker is configured with
#   --dfam_extra partitions of the same Dfam release added to the configured one
# Neither: RepeatMasker's own database (in the container, Dfam 4.0 base + curated
# consensus, which is all RMBlast needs for -species human).
DFAM_DB=""
DFAM_EXTRA=""
# --uncurated: RepeatMasker -uncurated, uncurated Dfam families as well as curated ones.
# They live in their own partitions (dfamNN.uncurated.*), none of which are built in.
UNCURATED=0
# --rm_gc: RepeatMasker -gc, the background GC% its scoring matrices are chosen for. Empty
# leaves the choice to RepeatMasker, which picks per internal batch: 43% for a batch of
# several sequences, the sequence's own GC for a lone one over 2 kb (see README Notes).
RM_GC=""
STOCK_FAMDB_DIR=""        # the FamDB package (famdb.py) RepeatMasker is configured with
STOCK_FAMDB_DATA_DIR=""   # the directory of *.h5 files it reads by default

NSPLIT_FILES=500
# Under a scheduler, use the allocation rather than the whole node. Plain `nproc`
# (unlike `nproc --all`) honours the cpuset/affinity mask, so it is already
# correct on a shared node; $PBS_NCPUS is preferred where set to make it explicit.
NTHREADS=${PBS_NCPUS:-$(nproc)}
MAX_JOBS=""          # Max number of RepeatMasker processes to run in parallel (default: NTHREADS)
THREADS_PER_JOB=1    # Threads per RepeatMasker job; derived in resolve_thread_counts

# Parameters (change if necessary)
MIN_SV_COVERAGE=0.05 #The minimum intersection between a repeat element and SV (aka sv_coverage) e.g. 0.05 (5%) (0 < min_sv_coverage < 1)
MIN_CLASS_SV_COVERAGE=0.25 #The minimum class sv coverage by repeat elements to be considered repetitive
MIN_TOTAL_SV_COVERAGE=0.75 #The minimum total sv coverage by repeat elements to be considered repetitive
MAX_TRF_OVERLAP=0.1 #The maximum TRF element overlap to be considered non-overlapping (0 < max_trf_overlap < 1)
INTERVAL=0.05
DIAGRAM_LEN=100
RM_CHAIN=1 # Flag to chain RepeatMasker alignments per element for the RM chain diagram block (1 = yes, 0 = no)
RM_CHAIN_JUMP_COST=30 # Cost (in matched bases) of each jump between RepeatMasker alignment blocks when chaining them per element

FLAG_DELETE_TMP_FILES=1 # Flag to delete temporary files (1 = yes, 0 = no)
FLAG_OVERWRITE=0 # Flag to overwrite existing output files (1 = yes, 0 = no)
RESUME=0 # Flag to resume from existing TRF/RepeatMasker outputs (1 = yes, 0 = no)
POPULATE_IDS=0 # Flag to let SVscanner fill missing/duplicated VCF IDs in a copy of the VCF (1 = yes, 0 = no)
NEED_ID_REWRITE=0 # Set by check_vcf_ids when the input VCF needs its IDs rewritten
SKIP_REF_CHECK=0 # Flag to skip checking that --ref is the reference the VCF was called against (1 = yes, 0 = no)
PREFIX="" # Prefix for output files (default: None)

# Python and bash scripts (keep as it is)
CHECK_REQUIRED="${SVSCANNER_HOME}/scripts/check_required_python.sh"
EXTRACT_SV_FLANKINGS="${SVSCANNER_HOME}/src/extract_sv.py"
RM_ALIGN_CHAIN="${SVSCANNER_HOME}/src/rm_align_chain.py"
ANNOTATION="${SVSCANNER_HOME}/src/repeat_annotation.py"
PLOT="${SVSCANNER_HOME}/src/generate_plots.py"
TRACEBACK_PLOT="${SVSCANNER_HOME}/src/generate_traceback_plots.py"
VCF_ANNOTATE="${SVSCANNER_HOME}/src/annotate_vcf.py"
CHECK_VCF_IDS="${SVSCANNER_HOME}/src/check_vcf_ids.py"
CHECK_REF="${SVSCANNER_HOME}/src/check_ref.py"

# Function to show usage
usage() {
    echo "Usage: $0 --out DIR --vcf FILE --ref FILE [options]"
    echo "Required arguments:"
    echo "  --out DIR      Path to output directory"
    echo "  --vcf FILE     Path to SV VCF file"
    echo "  --ref FILE     Path to reference FASTA file"
    echo "Optional arguments:"
    echo "  --prefix NAME           Prefix for output files (default: None; e.g. Project_, Project.)"
    echo "  --str_bed FILE          Path to STR BED file (default: $STR_BED)"
    echo "  --species NAME          Species name for RepeatMasker (default: $SPECIES)"
    echo "  --engine NAME           RepeatMasker search engine, e.g. rmblast, hmmer (default: RepeatMasker's own configured default)"
    echo "  --dfam_db DIR           A complete Dfam database (FamDB 3.0: base file dfamNN.0.h5 plus partitions)"
    echo "                          to use instead of the built-in one, e.g. another Dfam release."
    echo "  --dfam_extra DIR        Dfam partitions (*.h5, no base file) to add to the built-in database, e.g."
    echo "                          dfam40.curated.hmm.* for --engine hmmer. Must be the same Dfam release."
    echo "                          Default (neither): the built-in database (Dfam 4.0 curated consensus in the container)."
    echo "  --uncurated             Use uncurated Dfam families as well as curated ones (RepeatMasker -uncurated)."
    echo "                          Needs uncurated partitions covering --species, added with --dfam_extra or present in"
    echo "                          --dfam_db: dfamNN.uncurated.consensus.* for RMBlast, dfamNN.uncurated.hmm.* for HMMER."
    echo "  --rm_gc INT             Background GC% for RepeatMasker's scoring matrices (RepeatMasker -gc, 1-100), e.g. 43"
    echo "                          so that annotations do not depend on how SVs are grouped into files (default: not set,"
    echo "                          RepeatMasker chooses; see README Notes)"
    echo "  --min_sv_coverage VAL   Minimum intersection between a repeat element and SV (aka sv_coverage) (default: $MIN_SV_COVERAGE)"
    echo "  --min_class_sv_coverage VAL Minimum class sv coverage by repeat elements to be considered repetitive (default: $MIN_CLASS_SV_COVERAGE)"
    echo "  --min_total_sv_coverage VAL Minimum total sv coverage by repeat elements to be considered repetitive (default: $MIN_TOTAL_SV_COVERAGE)"
    echo "  --max_trf_overlap VAL   Maximum TRF element overlap to be considered non-overlapping (default: $MAX_TRF_OVERLAP)"
    echo "  --interval VAL          Interval value (default: $INTERVAL)"
    echo "  --diagram_len VAL       Diagram length (default: $DIAGRAM_LEN)"
    echo "  --no_rm_chain           Do not chain RepeatMasker alignments per element: RepeatMasker runs without -a and the"
    echo "                          diagram has no RM chain block (default: chain)"
    echo "  --rm_chain_jump_cost VAL Cost, in matched bases, of each jump between RepeatMasker alignment blocks when they are"
    echo "                          chained into one alignment per element (default: $RM_CHAIN_JUMP_COST). --resume reuses the"
    echo "                          existing rm_chain.tab, so a different value needs a full run."
    echo "  --nsplit_files INT      Number of split files (default: $NSPLIT_FILES)"
    echo "  --keep_tmp_files        Keep temporary files (default: delete)"
    echo "  --overwrite             Overwrite existing output files (default: no overwrite)"
    echo "  --resume                Reuse existing info/rm/trf .tab files and skip TRF + RepeatMasker (default: full run)"
    echo "  --populate_ids          Fill missing and duplicated IDs in the VCF ID column (default: exit with an error)."
    echo "                          SVscanner joins its annotations onto the VCF by ID. The input VCF is never edited;"
    echo "                          a repaired copy is written to the output directory and annotated instead."
    echo "  --skip_ref_check        Do not check that --ref is the reference the VCF was called against (default: check,"
    echo "                          and exit with an error on a contig length or REF allele mismatch)"
    echo "  --nthread INT           Number of threads to use (default: \$PBS_NCPUS if set, else all threads available to this process)"
    echo "  --njob INT              Number of parallel jobs for RepeatMasker (default: same as --nthread)"
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
                VCF=$(realpath -e "$2" 2>/dev/null) || die "VCF file not found: $2"; shift 2;;
            --ref)
                REF=$(realpath -e "$2" 2>/dev/null) || die "Reference file not found: $2"; shift 2;;
            --prefix)
                PREFIX="$2"; shift 2;;
            --str_bed)
                STR_BED=$(realpath -e "$2" 2>/dev/null) || die "STR BED file not found: $2"; shift 2;;
            --species)
                SPECIES="$2"; shift 2;;
            --engine)
                ENGINE="$2"; shift 2;;
            --dfam_db)
                DFAM_DB=$(realpath -e "$2" 2>/dev/null) || die "Dfam database directory not found: $2"; shift 2;;
            --dfam_extra)
                DFAM_EXTRA=$(realpath -e "$2" 2>/dev/null) || die "Dfam partition directory not found: $2"; shift 2;;
            --uncurated)
                UNCURATED=1; shift;;
            --rm_gc)
                RM_GC="$2"; shift 2;;
            --dfam_dir)
                die "--dfam_dir has been removed. SVscanner now needs RepeatMasker ${RM_MIN_VERSION}+ (Dfam 4.0+, FamDB 3.0); use\n  --dfam_db DIR     for a complete Dfam database used instead of the built-in one, or\n  --dfam_extra DIR  for partitions added to the built-in database.\nThe container's built-in Dfam 4.0 covers -species human with RMBlast without either.";;
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
            --no_rm_chain)
                RM_CHAIN=0; shift;;
            --rm_chain_jump_cost)
                RM_CHAIN_JUMP_COST="$2"; shift 2;;
            --nsplit_files)
                NSPLIT_FILES="$2"; shift 2;;
            --keep_tmp_files)
                FLAG_DELETE_TMP_FILES=0; shift;;
            --overwrite)
                FLAG_OVERWRITE=1; shift;;
            --resume)
                RESUME=1; shift;;
            --populate_ids)
                POPULATE_IDS=1; shift;;
            --skip_ref_check)
                SKIP_REF_CHECK=1; shift;;
            --nthread)
                NTHREADS="$2"; shift 2;;
            --njob)
                MAX_JOBS="$2"; shift 2;;
            --help)
                usage;;
            --version)
                echo "$VERSION"; exit 0;;
            *)
                echo "Unknown argument: $1"; usage;;
        esac
    done

    # Check required arguments. --ref is only used to extract flanking regions, which --resume skips.
    if [[ -z "$OUTPUT_DIR" || -z "$VCF" ]]; then
        echo "Error: --out and --vcf are required."
        usage
    fi
    if [[ ${RESUME} -ne 1 && -z "$REF" ]]; then
        echo "Error: --ref is required for a full run."
        usage
    fi

    # --resume reuses the existing output dir; --overwrite would delete it. Disallow both.
    if [[ ${RESUME} -eq 1 && ${FLAG_OVERWRITE} -eq 1 ]]; then
        die "Error: --resume and --overwrite are mutually exclusive (--overwrite deletes the .tab files --resume needs)."
    fi
    if [[ -n "${DFAM_DB}" && -n "${DFAM_EXTRA}" ]]; then
        die "Error: --dfam_db and --dfam_extra are mutually exclusive.\n--dfam_db replaces the built-in Dfam database; put any extra partitions in that same directory."
    fi
    # RepeatMasker treats -gc 0 as not given, so 0 would silently fall back to its own choice.
    if [[ -n "${RM_GC}" ]] && ! [[ "${RM_GC}" =~ ^[0-9]+$ && "${RM_GC}" -ge 1 && "${RM_GC}" -le 100 ]]; then
        die "Error: --rm_gc must be an integer from 1 to 100 (got '${RM_GC}')."
    fi

    # Internal directories (keep as it is)
    EXTRACT_SV_FLANKS_OUT=${OUTPUT_DIR}/${PREFIX}extract_sv_flanks_out
    ANNOTATIONS_OUT=${OUTPUT_DIR}/${PREFIX}annotations_out
    RM_TMP=${OUTPUT_DIR}/RMtmp
    # Built-in + --dfam_extra partitions, as symlinks; only created with --dfam_extra.
    RM_FAMDB_DATA=${OUTPUT_DIR}/${PREFIX}rm_famdb
    # RepeatMasker's library cache (REPEATMASKER_CACHE_DIR), private to this run.
    RM_CACHE_DIR=${OUTPUT_DIR}/${PREFIX}rm_cache
    # File Intermediates (keep as it is)
    INFO_FILE=${OUTPUT_DIR}/${PREFIX}info.tab
    RM_FILE=${OUTPUT_DIR}/${PREFIX}rm.tab
    # One chained alignment per RepeatMasker element, from the .align files (see chain_rm_alignments)
    RM_CHAIN_FILE=${OUTPUT_DIR}/${PREFIX}rm_chain.tab
    TRF_FILE=${OUTPUT_DIR}/${PREFIX}trf.tab

    # Final Outputs (change if necessary)
    ANNOTATED_VCF=${OUTPUT_DIR}/${PREFIX}annotated.vcf

    # Repaired copy of the input VCF; only written when --populate_ids was needed.
    # Kept after the run: --resume has to annotate this file, not the original.
    POPULATED_VCF=${OUTPUT_DIR}/${PREFIX}ids_populated.vcf.gz
    POPULATED_IDS_MAP=${OUTPUT_DIR}/${PREFIX}ids_populated.map.tsv
}

resolve_thread_counts() {
    # One RepeatMasker process per available thread by default, matching the
    # previous behaviour on a full node (48 jobs x 1 thread).
    [[ -z "${MAX_JOBS}" ]] && MAX_JOBS=${NTHREADS}

    # Never run more parallel jobs than we have threads, and never hand
    # RepeatMasker '-pa 0', which integer division would otherwise produce
    # whenever MAX_JOBS exceeds NTHREADS.
    (( MAX_JOBS < 1 )) && MAX_JOBS=1
    (( MAX_JOBS > NTHREADS )) && MAX_JOBS=${NTHREADS}

    THREADS_PER_JOB=$(( NTHREADS / MAX_JOBS ))
    (( THREADS_PER_JOB < 1 )) && THREADS_PER_JOB=1
}

check_binary() {
    local name=$1
    shift
    local candidates=("$@")

    for cmd in "${candidates[@]}"; do
        if [[ -f "$cmd" && -x "$cmd" ]]; then
            echo "$(realpath "$cmd")"
            return 0
        fi

        local full_path
        full_path=$(command -v "$cmd" 2>/dev/null)
        if [[ -n "$full_path" ]]; then
            echo "$full_path"
            return 0
        fi
    done

    die "$name binary not found in any of: ${candidates[*]}"
}

check_repeatmasker() {
    # Version gate, then the FamDB package and data directory RepeatMasker is configured
    # with: --dfam_extra builds on that data directory, and the log should say which
    # database a run used.
    local rm_install
    rm_install=$(dirname "$(realpath "${REPEAT_MASKER}")")
    RM_VERSION=$(${REPEAT_MASKER} -v 2>/dev/null | awk '/RepeatMasker version/ { print $3; exit }')
    [[ -n "${RM_VERSION}" ]] || die "Could not read the RepeatMasker version from '${REPEAT_MASKER} -v'."
    if [[ "$(printf '%s\n%s\n' "${RM_MIN_VERSION}" "${RM_VERSION}" | sort -V | head -n1)" != "${RM_MIN_VERSION}" ]]; then
        die "RepeatMasker ${RM_VERSION} is too old: SVscanner needs ${RM_MIN_VERSION} or later.\nThe SVscanner container image ships a supported version (see docs/docker.md)."
    fi

    # RepeatMasker's own precedence for each setting: environment, then configure's value.
    cfg() {
        perl -I"${rm_install}" -MRepeatMaskerConfig \
            -e 'print $RepeatMaskerConfig::configuration->{$ARGV[0]}->{value} // ""' "$1" 2>/dev/null
    }
    STOCK_FAMDB_DIR="${FAMDB_DIR:-$(cfg FAMDB_DIR)}"
    [[ -n "${STOCK_FAMDB_DIR}" && -x "${STOCK_FAMDB_DIR}/famdb.py" ]] ||
        die "RepeatMasker ${RM_VERSION} at ${rm_install} has no FamDB configured (FAMDB_DIR is '${STOCK_FAMDB_DIR}'); -species needs it."
    # Unset, FamDB falls back to its own famdb.conf or Libraries/famdb; 'info' reports
    # whichever it resolved.
    STOCK_FAMDB_DATA_DIR="${FAMDB_DATA_DIR:-$(cfg FAMDB_DATA_DIR)}"
    if [[ -z "${STOCK_FAMDB_DATA_DIR}" ]]; then
        STOCK_FAMDB_DATA_DIR=$("${STOCK_FAMDB_DIR}/famdb.py" info 2>/dev/null |
            awk -F': *' '/^FamDB Directory/ { print $2; exit }')
    fi
    [[ -d "${STOCK_FAMDB_DATA_DIR}" ]] ||
        die "Could not find the FamDB data directory RepeatMasker uses (FamDB package: ${STOCK_FAMDB_DIR})."
    echo "RepeatMasker: ${RM_VERSION} (FamDB package: ${STOCK_FAMDB_DIR}; data: ${STOCK_FAMDB_DATA_DIR})"
}

check_required() {
    # 1. Check Python version
    ${CHECK_REQUIRED} || die "Python version check failed"

    [ -z "$OUTPUT_DIR" ] && die "OUTPUT_DIR is not set"
    [[ ${RESUME} -ne 1 && -z "$REF" ]] && die "REF is not set"
    [ -z "$VCF" ] && die "VCF is not set"
    [ -z "$STR_BED" ] && die "STR_BED is not set"

    echo "Output dir: ${OUTPUT_DIR}"
    echo "Reference: ${REF}"
    echo "Input SV VCF: ${VCF}"
    echo "Input BED: ${STR_BED}"
    echo "RepeatMasker matrix GC: ${RM_GC:-RepeatMasker default (43% for multi-sequence batches, own GC for a lone sequence >2 kb)}"

    echo "Number of Threads: ${NTHREADS}"
    echo "Number of RepeatMasker jobs: ${MAX_JOBS} (${THREADS_PER_JOB} thread(s) each)"

    # TRF and RepeatMasker are only needed for a full run; --resume reuses their .tab outputs.
    if [[ ${RESUME} -ne 1 ]]; then
        TRF_BINARY=$(check_binary "TRF" "trf" "trf409.linux64") || exit 1
        REPEAT_MASKER=$(check_binary "RepeatMasker" "RepeatMasker" "RepeatMasker/RepeatMasker") || exit 1
        parallel=$(command -v parallel) || die "parallel not found"
        echo "REPEAT_MASKER: ${REPEAT_MASKER}"
        echo "TRF_BINARY: ${TRF_BINARY}"
        echo "parallel: ${parallel}"
        check_repeatmasker

        # Fail on an empty --dfam_db/--dfam_extra before doing any work; the full checks,
        # which need FamDB, are in setup_dfam_library.
        local dir
        for dir in "${DFAM_DB}" "${DFAM_EXTRA}"; do
            [[ -n "${dir}" ]] || continue
            [[ -n "$(find "${dir}" -maxdepth 1 -name '*.h5' -print -quit)" ]] ||
                die "No FamDB files (*.h5) in ${dir}.\nDfam 4.0+ releases are at https://www.dfam.org/releases/ ."
        done
    else
        echo "Resume mode: skipping TRF + RepeatMasker (reusing existing .tab files)"
    fi

    BCFTOOLS=$(check_binary "bcftools" "bcftools" "bcftools-1.21/bcftools") || exit 1
    BGZIP=$(check_binary "bgzip" "bgzip" "htslib-1.21/bgzip") || exit 1
    echo "BCFTOOLS: ${BCFTOOLS}"
    echo "BGZIP: ${BGZIP}"

}

check_vcf_ids() {
    # SVscanner joins its annotations back onto the VCF by the ID column, so every
    # record needs a present, unique ID. Check the VCF now, before any output directory
    # is created, so a failure leaves nothing behind that blocks the retry.
    # check_vcf_ids.py exits 0 = usable, 1 = repairable, 2 = not repairable, 3 = I/O error.
    local rc=0
    echo "Checking the VCF ID column..."

    if [[ ${RESUME} -eq 1 && -s "${POPULATED_VCF}" ]]; then
        # The previous run annotated a repaired copy; the saved .tab files carry its IDs.
        echo "Resume mode: using the VCF with populated IDs from the previous run: ${POPULATED_VCF}"
        VCF="${POPULATED_VCF}"
    fi

    python3 "${CHECK_VCF_IDS}" --vcf "${VCF}" --mode check || rc=$?
    case ${rc} in
        0)
            echo "done"
            return 0;;
        1)
            if [[ ${RESUME} -eq 1 ]]; then
                die "Cannot resume: the VCF ID column is not usable and --resume cannot repair it.\nRe-run without --resume (add --populate_ids to let SVscanner fill the IDs)."
            elif [[ ${POPULATE_IDS} -eq 1 ]]; then
                echo "Warning: the VCF ID column is not usable; --populate_ids is set, so a repaired copy will be written to ${POPULATED_VCF}" >&2
                NEED_ID_REWRITE=1
                return 0
            fi
            die "The VCF ID column is not usable (see above). SVscanner annotates the VCF by ID.\nAdd every missing or duplicated ID yourself, or re-run with --populate_ids to let SVscanner fill them in a copy of the VCF (the input is never edited).";;
        2)
            die "The VCF ID column is not usable and SVscanner will not rewrite it (see above), with or without --populate_ids.";;
        *)
            die "VCF ID check failed (exit code ${rc})";;
    esac
}

check_reference() {
    # A wrong reference does not make the run fail: DEL/DUP/INV/BND sequences and every
    # flank are fetched from it and annotated as if they were right. Check it before any
    # output directory is created, so a failure leaves nothing behind that blocks the retry.
    # check_ref.py exits 0 = matches, 1 = does not match, 3 = I/O error.
    if [[ ${SKIP_REF_CHECK} -eq 1 ]]; then
        echo "Warning: --skip_ref_check is set; not checking that the reference matches the VCF" >&2
        return 0
    fi
    local rc=0
    echo "Checking the reference against the VCF..."
    python3 "${CHECK_REF}" --vcf "${VCF}" --ref "${REF}" || rc=$?
    case ${rc} in
        0)
            echo "done";;
        1)
            die "The reference (${REF}) does not match the VCF (see above).\nUse the reference the VCF was called against, or re-run with --skip_ref_check if you are sure it is right.";;
        *)
            die "Reference check failed (exit code ${rc})";;
    esac
}

populate_vcf_ids() {
    # Second stage of check_vcf_ids, once the output directory exists: write the repaired
    # copy, verify it, and annotate that from here on.
    [[ ${NEED_ID_REWRITE} -eq 1 ]] || return 0
    echo "Populating the VCF ID column in a copy of the VCF..."
    python3 "${CHECK_VCF_IDS}" --vcf "${VCF}" --mode rewrite --out "${POPULATED_VCF}" --map "${POPULATED_IDS_MAP}" || die "failed to populate VCF IDs"
    python3 "${CHECK_VCF_IDS}" --vcf "${POPULATED_VCF}" --mode check || die "the VCF with populated IDs still fails the ID check"
    VCF="${POPULATED_VCF}"
    echo "Input SV VCF is now: ${VCF}"
    echo "done"
}

create_output_dir() {
    if [[ -d "${OUTPUT_DIR}" ]]; then
        if [[ ${FLAG_OVERWRITE} -eq 1 ]]; then
            echo "Output directory ${OUTPUT_DIR} already exists. Overwriting..."
            rm -rf "${OUTPUT_DIR}" || die "Failed to remove existing output directory ${OUTPUT_DIR}"
        else
            die "Output directory ${OUTPUT_DIR} already exists. Please choose a different output directory or delete the existing one."
        fi
    fi
	mkdir -p "${OUTPUT_DIR}" || die "Failed creating ${OUTPUT_DIR}"
}

setup_dfam_library() {
    # Point RepeatMasker at the Dfam database for this run, via FAMDB_DATA_DIR (the
    # directory of FamDB *.h5 files; RepeatMasker hands it to famdb.py as -i).
    #
    # A Dfam release is one base file (dfamNN.0.h5: names and taxonomy, always needed)
    # plus partitions by component and taxon (curated/uncurated, consensus/HMM). Partitions
    # only work with the base file of their own release.
    local base_re='^dfam[0-9]+\.0\.h5$'
    local data_dir h5 name
    local supplied=()

    if [[ -z "${DFAM_DB}" && -z "${DFAM_EXTRA}" ]]; then
        echo "Dfam: built-in database (${STOCK_FAMDB_DATA_DIR})"
        return 0
    fi

    if [[ -n "${DFAM_DB}" ]]; then
        local bases=()
        for h5 in "${DFAM_DB}"/*.h5; do
            name=$(basename "${h5}")
            if [[ "${name}" =~ ${base_re} ]]; then bases+=("${name}"); else supplied+=("${name}"); fi
        done
        (( ${#bases[@]} == 1 )) || {
            (( ${#bases[@]} == 0 )) &&
                die "--dfam_db ${DFAM_DB} has no base file (dfamNN.0.h5), so it is not a complete Dfam database.\nTo add partitions to the built-in database, use --dfam_extra instead."
            die "--dfam_db ${DFAM_DB} has more than one base file (${bases[*]}); keep one Dfam release per directory."
        }
        data_dir="${DFAM_DB}"
        echo "Dfam: ${DFAM_DB} (base ${bases[0]}), replacing the built-in database"
    else
        # Partitions name their release (dfam40.curated.hmm.1.h5); they must match the
        # built-in base file, which is linked in beside them.
        local stock_base release
        stock_base=$(find "${STOCK_FAMDB_DATA_DIR}/" -maxdepth 1 -name 'dfam*.0.h5' -printf '%f\n' | head -n1)
        [[ -n "${stock_base}" ]] || die "No base file (dfamNN.0.h5) in the built-in database ${STOCK_FAMDB_DATA_DIR}."
        release="${stock_base%%.*}"
        for h5 in "${DFAM_EXTRA}"/*.h5; do
            name=$(basename "${h5}")
            [[ "${name}" =~ ${base_re} ]] &&
                die "--dfam_extra ${DFAM_EXTRA} contains ${name}, the base file of a complete Dfam database.\nUse --dfam_db for a complete database; --dfam_extra takes partitions only."
            [[ "${name}" == "${release}."* ]] ||
                die "--dfam_extra file ${name} is not from the built-in Dfam release (${release}).\nPartitions only work with their own release's base file: use --dfam_db with a complete database of that release."
            supplied+=("${name}")
        done

        rm -rf "${RM_FAMDB_DATA}" || die "failed to clear ${RM_FAMDB_DATA}"
        mkdir -p "${RM_FAMDB_DATA}" || die "failed to create ${RM_FAMDB_DATA}"
        # Built-in files first, then --dfam_extra, which overrides a same-named file.
        for h5 in "${STOCK_FAMDB_DATA_DIR}"/*.h5 "${DFAM_EXTRA}"/*.h5; do
            ln -sfn "$(realpath "${h5}")" "${RM_FAMDB_DATA}/$(basename "${h5}")" || die "failed to link ${h5}"
        done
        data_dir="${RM_FAMDB_DATA}"
        echo "Dfam: built-in database (${STOCK_FAMDB_DATA_DIR}) + ${#supplied[@]} partition(s) from ${DFAM_EXTRA}"
    fi

    # FamDB skips files it cannot use (another release's partitions, Dfam <= 3.9's FamDB
    # 2.0 format) without an error, and RepeatMasker would then annotate against less
    # than this log claims. 'info' names each partition it loaded as [file.h5]; the base
    # file holds no families and is not listed.
    local info_log="${OUTPUT_DIR}/${PREFIX}famdb_info.log"
    local unused=()
    if "${STOCK_FAMDB_DIR}/famdb.py" -i "${data_dir}" info > "${info_log}" 2>&1; then
        for name in "${supplied[@]}"; do
            grep -qF "[${name}]" "${info_log}" || unused+=("${name}")
        done
    else
        unused=("(famdb.py info failed)")
    fi
    if (( ${#unused[@]} )); then
        echo "FamDB did not load: ${unused[*]} (full output: ${info_log})" >&2
        die "The Dfam files given cannot be used by RepeatMasker ${RM_VERSION}.\nIt reads Dfam 4.0+ (FamDB 3.0) files; Dfam 3.9 and earlier are not supported."
    fi
    rm -f "${info_log}"

    # Exported so the RepeatMasker processes spawned by parallel inherit it.
    export FAMDB_DATA_DIR="${data_dir}"
    echo "FamDB files in use (FAMDB_DATA_DIR=${FAMDB_DATA_DIR}):"
    for h5 in "${FAMDB_DATA_DIR}"/*.h5; do
        echo "  $(basename "${h5}") -> $(realpath "${h5}")"
    done
}

check_resume_inputs() {
    # Verify the checkpoint files from a previous run are present so we can skip TRF + RepeatMasker.
    echo "Resume mode: validating existing checkpoint files in ${OUTPUT_DIR}..."

    [[ -d "${OUTPUT_DIR}" ]] || die "Cannot resume: output directory ${OUTPUT_DIR} does not exist."

    local missing=0
    for f in "${INFO_FILE}" "${RM_FILE}" "${TRF_FILE}"; do
        if [[ ! -s "$f" ]]; then
            echo "  missing or empty: $f" >&2
            missing=1
        fi
    done

    if [[ ${missing} -eq 1 ]]; then
        die "Cannot resume: required checkpoint file(s) above are missing or empty.\nThese are produced by a previous full run and must use the same --prefix ('${PREFIX}')."
    fi

    echo "  found: ${INFO_FILE}"
    echo "  found: ${RM_FILE}"
    echo "  found: ${TRF_FILE}"
    # Optional: runs made with --no_rm_chain or before the RepeatMasker chain step have no
    # rm_chain.tab. Their .align files are gone, so the diagram has no RM chain block.
    if (( RM_CHAIN )); then
        if [[ -s "${RM_CHAIN_FILE}" ]]; then
            echo "  found: ${RM_CHAIN_FILE}"
        else
            echo "  not found: ${RM_CHAIN_FILE} (the diagram will have no RM chain block; a full run makes it)"
        fi
    fi
    echo "done"
}

extract_flanking_regions() {
    # 1) Extract sequence and flanking regions for variants
    echo "Extracting structural variant sequences from VCF..."
    python3 ${EXTRACT_SV_FLANKINGS} --vcf ${VCF} --ref ${REF} --out ${EXTRACT_SV_FLANKS_OUT} --min 10 -n ${NSPLIT_FILES} --info ${INFO_FILE} || die "failed"
    echo "done"
}

run_trf() {
    # set -x
    # 3) Run Tandem Repeat Finder and RepeatMasker - wait for both to complete
    echo "Running Tandem Repeat Finder..."
    T4=$(date +%s)
    # TRF 4.09/4.10 truncate input paths of ~260 characters or more and exit with "not
    # found", so a deep --out stopped the run with no useful message. Run each job from
    # inside the directory and give TRF the bare file name.
    find "${EXTRACT_SV_FLANKS_OUT}" -name "*.fa" | parallel -j ${NTHREADS} "cd {//} && ${TRF_BINARY} {/} 2 7 7 80 10 50 500 -h -ngs > {/.}.dat" || die "failed"
    T5=$(date +%s) || die "failed to get T3"
    TRF_TIME=$((T5 - T4)) || die "failed to calculate time"
    echo "done. Tandem Repeat Finder took ${TRF_TIME} seconds"
}

set_rm_search_args() {
    # The options that decide which library RepeatMasker searches, shared by the warm-up
    # and the fan-out so both use the same one. Kept as an array: a species such as
    # "danio rerio" is one argument. The fan-out passes them through GNU parallel's
    # command string, which a shell parses again, so it shell-quotes them there.
    RM_SEARCH_ARGS=(-species "${SPECIES}")
    [[ -n "${ENGINE}" ]] && RM_SEARCH_ARGS+=(-e "${ENGINE}")
    (( UNCURATED )) && RM_SEARCH_ARGS+=(-uncurated)
    [[ -n "${RM_GC}" ]] && RM_SEARCH_ARGS+=(-gc "${RM_GC}")
    return 0
}

warm_repeatmasker_cache() {
    # RepeatMasker derives a per-species library from FamDB on first use and caches it.
    # Build it once, single process, before the fan-out. Concurrent builds are safe since
    # RepeatMasker 4.2.6, so this is no longer about races: it stops a run whose species
    # resolves to nothing before extraction and TRF, and it is the only place that records
    # which database and how many families the run annotated against.
    #
    # REPEATMASKER_CACHE_DIR is searched before the install's Libraries and
    # ~/.RepeatMaskerCache, so the cache is private to this run: nothing is written into a
    # read-only image or the user's shared home, and nothing from another run is reused.
    local warm_dir="${OUTPUT_DIR}/${PREFIX}rm_warmup"
    local warm_fa="${warm_dir}/warmup.fa"
    local warm_log="${warm_dir}/warmup.log"

    rm -rf "${warm_dir}" || die "failed to clear ${warm_dir}"
    mkdir -p "${warm_dir}" "${RM_CACHE_DIR}" || die "failed to create ${warm_dir} or ${RM_CACHE_DIR}"
    # Exported so the RepeatMasker processes spawned by parallel inherit it.
    export REPEATMASKER_CACHE_DIR="${RM_CACHE_DIR}"

    # Content is irrelevant - the point is to make RepeatMasker build the library.
    {
        echo ">svscanner_library_warmup"
        echo "ACGTTGCAAGCTTAGCCATGGATCCGTAAGCTTGGATCCAAGCTTGCATGCCTGCAGGTCG"
        echo "TTACGGCATTCAGGCCTAAGCTTGGCCTAGGCATTAGCCGTAAGGCTTAACCGGATTCCAG"
        echo "GGCATTAACCGGTTAAGCCTTAAGGCCATTAACCGGTTAACCGGAATTCCGGTTAACCGGA"
        echo "CCTTAAGGCCTTAAGGCCTTAAGGAATTCCGGAATTCCGGTTAACCGGTTAACCGGAATTC"
    } > "${warm_fa}" || die "failed to write ${warm_fa}"

    echo "Preparing RepeatMasker libraries for species '${SPECIES}' (single process)..."
    local T_WARM_START T_WARM_END
    T_WARM_START=$(date +%s)
    # Run from inside warm_dir: createTempDir uses cwd(), not -dir, so RepeatMasker would
    # otherwise scatter RM_<pid> directories through the caller's directory. Same reason
    # run_repeatmasker cds into RM_TMP before the fan-out.
    set_rm_search_args
    if ! ( cd "${warm_dir}" && ${REPEAT_MASKER} "${warm_fa}" -pa 1 -dir "${warm_dir}" "${RM_SEARCH_ARGS[@]}" > "${warm_log}" 2>&1 ); then
        echo "RepeatMasker reported (${warm_log}):" >&2
        tail -n 20 "${warm_log}" | sed 's/^/  /' >&2
        die "RepeatMasker could not build its libraries for species '${SPECIES}'"
    fi
    T_WARM_END=$(date +%s)

    # RepeatMasker reports the database it opened and the families it resolved only into
    # its per-chunk logs, which cleanup deletes - so a successful run leaves no record of
    # what it annotated against. This warm-up used the same library, so lift its summary
    # into the main log. 'families in ancestor taxa' is the number to watch: 1387 for
    # -species human with Dfam 4.0 curated consensus.
    echo "RepeatMasker library summary:"
    grep -E 'RepeatMasker version|Search Engine|Master RepeatMasker Database|Title[[:space:]]*:|Version[[:space:]]*:|Date[[:space:]]*:|Families[[:space:]]*:|Taxonomy ID|Including (only curated|curated and uncurated) families|families in ancestor taxa|uncurated families|Building (general|species) libraries|partition' \
        "${warm_log}" | head -n 20 | sed 's/^/  /'

    # Two ways to end up with a library that annotates nothing, both of which
    # RepeatMasker exits 0 on: famdb.py raising when a clade's descendants live in an
    # uninstalled partition (RepeatMasker ignores its exit status and carries on), and a
    # species that resolves to no families. Reporting success there is worse than
    # stopping. Both trigger on positive evidence only, so a change in RepeatMasker's
    # wording can never block a working run.
    #
    # Some of the families Dfam has for this search are in partitions that are not
    # installed (e.g. --uncurated without the uncurated partitions, or --engine hmmer with
    # only some of the curated HMM partitions). RepeatMasker warns and carries on; famdb.py
    # usually then raises too, so this is checked before the traceback for the clearer
    # message. Never fires for the built-in database under RMBlast: its curated consensus
    # partition is root-level and complete.
    if grep -q 'Some Dfam components/partitions applicable to this search are not present' "${warm_log}"; then
        echo "RepeatMasker reported missing partitions for species '${SPECIES}':" >&2
        grep -E 'Including (only curated|curated and uncurated) families|families in ancestor taxa' "${warm_log}" | sed 's/^/  /' >&2
        local hint=""
        (( UNCURATED )) && hint=" (for --uncurated: dfam40.uncurated.consensus.* for RMBlast, dfam40.uncurated.hmm.* for HMMER)"
        die "Part of the Dfam library for species '${SPECIES}' is not installed, so RepeatMasker would annotate against an incomplete library.\nAdd the missing partitions with --dfam_extra${hint}, or use --dfam_db with a complete database.\n'famdb.py check \"${SPECIES}\"' lists them.\nFull log kept at: ${warm_log}"
    fi
    if grep -q 'Traceback (most recent call last)' "${warm_log}"; then
        echo "RepeatMasker's famdb.py failed while resolving species '${SPECIES}':" >&2
        grep -hE '^[A-Za-z_.]*(Error|Exception):' "${warm_log}" | sort -u | head -n 3 | sed 's/^/  /' >&2
        die "Cannot trust the RepeatMasker library for species '${SPECIES}'.\nThis usually means '${SPECIES}' has descendant taxa in a FamDB partition that is not installed - RepeatMasker would carry on and annotate against an empty library.\nAdd the partitions covering '${SPECIES}' (--dfam_extra, or --dfam_db for a complete database) or choose a species they cover.\nFull log kept at: ${warm_log}"
    fi
    # Counts are "present/total" (e.g. 1387/1387); none present means an empty library.
    if grep -qE ' 0(/[0-9]+)? families in ancestor taxa; 0(/[0-9]+)? lineage-specific families' "${warm_log}"; then
        die "RepeatMasker resolved no families for species '${SPECIES}' - every repeat annotation would be empty.\nAdd the partitions covering '${SPECIES}' (--dfam_extra, or --dfam_db for a complete database).\nFull log kept at: ${warm_log}"
    fi
    # With -uncurated, RepeatMasker prints this note when Dfam itself has no uncurated
    # families for the species (e.g. human, danio rerio) - the totals come from the base
    # file, so adding partitions cannot change it. It then carries on curated-only.
    if (( UNCURATED )) && grep -q 'No uncurated families were found for this taxon' "${warm_log}"; then
        die "--uncurated was given, but Dfam has no uncurated families for species '${SPECIES}', so the run would be curated-only.\nDrop --uncurated: no partition can add uncurated families for this species.\nFull log kept at: ${warm_log}"
    fi

    rm -rf "${warm_dir}" || die "failed to remove ${warm_dir}"
    echo "done. Library preparation took $((T_WARM_END - T_WARM_START)) seconds"
}

run_repeatmasker() {
    mkdir -p ${RM_TMP}
    echo "Running RepeatMasker..."
    T2=$(date +%s)
    cd ${RM_TMP}

    shopt -s nullglob
    fa_files=(${EXTRACT_SV_FLANKS_OUT}/*.fa)
    (( ${#fa_files[@]} == 0 )) && die "No FASTA files found in ${EXTRACT_SV_FLANKS_OUT}"
    shopt -u nullglob

    # Run RepeatMasker in parallel with error sensitivity.
    # MAX_JOBS and THREADS_PER_JOB are set by resolve_thread_counts.
    # -a writes the alignments (.align) that chain_rm_alignments reads; it leaves .out unchanged.
    set_rm_search_args
    local align_arg=""
    (( RM_CHAIN )) && align_arg="-a"
    if ! find ${EXTRACT_SV_FLANKS_OUT} -name "*.fa" | parallel --halt now,fail=1 -j "${MAX_JOBS}" "${REPEAT_MASKER} {} -pa ${THREADS_PER_JOB} ${align_arg} -html -gff -dir ${EXTRACT_SV_FLANKS_OUT} $(printf '%q ' "${RM_SEARCH_ARGS[@]}")> {}.log 2>&1"; then
        # parallel reports only the failing command line; the reason is in that job's own
        # log. Surface it - most often a FamDB partition that does not cover --species.
        local last_log
        last_log=$(find "${EXTRACT_SV_FLANKS_OUT}" -name '*.fa.log' -printf '%T@ %p\n' 2>/dev/null |
            sort -rn | head -n1 | cut -d' ' -f2-)
        if [[ -n "${last_log}" ]]; then
            echo "RepeatMasker reported (${last_log}):" >&2
            tail -n 20 "${last_log}" | sed 's/^/  /' >&2
        fi
        echo "Remaining RepeatMasker logs: ${EXTRACT_SV_FLANKS_OUT}/*.fa.log" >&2
        die "RepeatMasker failed"
    fi

    cd - || die "cd - failed"
    T3=$(date +%s) || die "failed to get T3"
    RM_TIME=$((T3 - T2)) || die "failed to calculate time"
    rm -r ${RM_TMP} || die "failed to remove ${RM_TMP}"
    echo "done. RepeatMasker took ${RM_TIME} seconds"
}

chain_rm_alignments() {
    # RepeatMasker reports one element as several alignment blocks under one .out ID (e.g.
    # one per copy of an expanded SVA VNTR). Chain them into one alignment per element and
    # keep only the per-element table; the .align files go with the temporary files.
    echo "Chaining RepeatMasker alignments per element..."
    python3 ${RM_ALIGN_CHAIN} \
        --rm_dir ${EXTRACT_SV_FLANKS_OUT} \
        --info ${INFO_FILE} \
        --out ${RM_CHAIN_FILE} \
        --jump_cost ${RM_CHAIN_JUMP_COST} \
        --threads ${NTHREADS} || die "failed to chain RepeatMasker alignments"
    echo "done"
}

process_trf_repeatmasker_output() {
    # Process the RepeatMasker files (remove the header and 16th column)
    echo "Process the RepeatMasker files (remove the header and 16th column)..."
    for rm_output in "${EXTRACT_SV_FLANKS_OUT}"/*.fa.out; do
        tail -n +4 "${rm_output}" | awk '{
            $16 = "";
            print $0
        }' OFS='\t' > "${rm_output}.tab" || die "failed"
    done
    # rm -rf "${EXTRACT_SV_FLANKS_OUT}/*.fa.out"

    cat ${EXTRACT_SV_FLANKS_OUT}/*.dat > ${TRF_FILE}
    cat ${EXTRACT_SV_FLANKS_OUT}/*.out.tab > ${RM_FILE}
    echo "done"
}

annotation() {
    echo "Annotating..."
    test -d "${ANNOTATIONS_OUT}" && rm -r "${ANNOTATIONS_OUT}"
    # Absent when --resume picks up a run made with --no_rm_chain or before the chain step.
    local rm_chain_args=()
    (( RM_CHAIN )) && [[ -s "${RM_CHAIN_FILE}" ]] && rm_chain_args=(--rm_chain "${RM_CHAIN_FILE}")
    python3 ${ANNOTATION} \
        "${rm_chain_args[@]}" \
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
    echo "done"
}

plot_classifications() {
    echo "Plotting..."
    python3 ${PLOT} \
        --out ${ANNOTATIONS_OUT}/plots\
        --tsv ${ANNOTATIONS_OUT}/plot_annotate.tsv || die "failed"
    echo "done"
}

plot_tracebacks() {
    echo "Plotting tracebacks..."
    python3 ${TRACEBACK_PLOT} \
        --output ${ANNOTATIONS_OUT}/traceback_plots.pdf \
        --traceback ${ANNOTATIONS_OUT}/traceback.tsv || die "failed"
    echo "done"
}

apply_annotations() {
    # apply annotations to the VCF
    echo "Applying annotations to the VCF..."
    python3 ${VCF_ANNOTATE} \
        --vcf ${VCF} \
        --header ${ANNOTATIONS_OUT}/vcf_header.txt \
        --tsv ${ANNOTATIONS_OUT}/vcf_annotate.tsv \
        --mate ${ANNOTATIONS_OUT}/vcf_annotate_bnd_mate.tsv \
        --output ${ANNOTATED_VCF} || die "failed to annotate VCF"
    echo "done"
}

sort_and_index_vcf() {
    # Sort and Index the annotated VCF
    echo "Sort and Index the annotated VCF..."
    # bcftools sort spills to --temp-dir once the annotated VCF exceeds its
    # 768M in-memory buffer (routine for large cohorts). Default temp-dir is
    # /tmp, which on HPC nodes (e.g. NCI Gadi) is often the same disk PBS
    # tracks as jobfs - a quota jobs don't request by default, so runs get
    # killed mid-sort. OUTPUT_DIR is already sized for this run's output, so
    # point temp files there instead.
    ${BCFTOOLS} sort --temp-dir ${OUTPUT_DIR}/sort_tmp -Oz -o ${ANNOTATED_VCF}.gz ${ANNOTATED_VCF} || die "${BCFTOOLS} sort failed"
    ${BCFTOOLS} index -t ${ANNOTATED_VCF}.gz || die "${BCFTOOLS} index failed"
    rm ${ANNOTATED_VCF}
    echo "done"
}

show_output_paths() {

    if [[ ${FLAG_DELETE_TMP_FILES} -eq 1 ]]; then
        echo "Deleting temporary files..."
        #if ${ANNOTATIONS_OUT}/plots/plots.pdf exists, move it to the output directory
        if [[ -f ${ANNOTATIONS_OUT}/plots/plots.pdf ]]; then
            mv ${ANNOTATIONS_OUT}/plots/plots.pdf ${OUTPUT_DIR}/${PREFIX}plots.pdf || die "failed to move plots to ${OUTPUT_DIR}"
        else
            mv ${ANNOTATIONS_OUT}/plots ${OUTPUT_DIR}/${PREFIX}plots || die "failed to move plots to ${OUTPUT_DIR}"
            echo "Plots dir: ${OUTPUT_DIR}/${PREFIX}plots"
        fi

        mv ${ANNOTATIONS_OUT}/diagram.txt ${OUTPUT_DIR}/${PREFIX}diagram.txt || die "failed to move diagram.txt to ${OUTPUT_DIR}"
        mv ${ANNOTATIONS_OUT}/rm_diagram.tsv ${OUTPUT_DIR}/${PREFIX}rm_diagram.tsv || die "failed to move rm_diagram.tsv to ${OUTPUT_DIR}"
        mv ${ANNOTATIONS_OUT}/trf_diagram.tsv ${OUTPUT_DIR}/${PREFIX}trf_diagram.tsv || die "failed to move trf_diagram.tsv to ${OUTPUT_DIR}"
        mv ${ANNOTATIONS_OUT}/traceback_plots.pdf ${OUTPUT_DIR}/${PREFIX}traceback_plots.pdf # only if generated
        mv ${ANNOTATIONS_OUT}/vcf_annotate.tsv ${OUTPUT_DIR}/${PREFIX}vcf_annotate.tsv || die "failed to move vcf_annotate.tsv to ${OUTPUT_DIR}"
        mv ${ANNOTATIONS_OUT}/vcf_annotate_bnd_mate.tsv ${OUTPUT_DIR}/${PREFIX}vcf_annotate_bnd_mate.tsv || die "failed to move vcf_annotate_bnd_mate.tsv to ${OUTPUT_DIR}"

        rm -rf ${ANNOTATIONS_OUT} || die "failed to remove ${ANNOTATIONS_OUT}"
        rm -rf ${EXTRACT_SV_FLANKS_OUT} || die "failed to remove ${EXTRACT_SV_FLANKS_OUT}"
        # Only ever symlinks to the partitions (see setup_dfam_library), never the files.
        [[ -d ${RM_FAMDB_DATA} ]] && { rm -rf ${RM_FAMDB_DATA} || die "failed to remove ${RM_FAMDB_DATA}"; }
        [[ -d ${RM_CACHE_DIR} ]] && { rm -rf ${RM_CACHE_DIR} || die "failed to remove ${RM_CACHE_DIR}"; }
        # rm -f ${INFO_FILE} || die "failed to remove ${INFO_FILE}"
        # rm -f ${RM_FILE} || die "failed to remove ${RM_FILE}"
        # rm -f ${TRF_FILE} || die "failed to remove ${TRF_FILE}"
        echo "Annotation outputs dir: ${OUTPUT_DIR}"
        echo "SV VCF with repeats annotated: ${ANNOTATED_VCF}.gz"
    else
        echo "Annotation outputs dir: ${ANNOTATIONS_OUT}"
        echo "Plots dir: ${ANNOTATIONS_OUT}/plots"
        echo "SV VCF with repeats annotated: ${ANNOTATED_VCF}.gz"
    fi

}

T0=$(date +%s)
# Banner after parse_args, so that --version and --help print only their own output.
parse_args "$@"
echo "SVscanner version: ${VERSION}"
resolve_thread_counts
check_required
if [[ ${RESUME} -eq 1 ]]; then
    check_resume_inputs
    check_vcf_ids
else
    check_vcf_ids
    check_reference
    create_output_dir
    populate_vcf_ids
    setup_dfam_library
    # Before extraction: a species the installed partitions cannot cover is fatal, and
    # there is no point spending an extraction and a TRF pass to find that out.
    warm_repeatmasker_cache
    extract_flanking_regions
    run_trf
    run_repeatmasker
    if [[ ${RM_CHAIN} -eq 1 ]]; then
        chain_rm_alignments
    fi
    process_trf_repeatmasker_output
fi
annotation
plot_classifications
plot_tracebacks
apply_annotations
sort_and_index_vcf
show_output_paths

T1=$(date +%s)
ELAPSED_TIME=$((T1 - T0))
echo "The SVscanner pipeline took ${ELAPSED_TIME} seconds"
echo "SVscanner version: ${VERSION}"
echo "$(date)"
echo "Success!"
exit 0
