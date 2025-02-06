#!/bin/bash

# set -x
RED='\033[0;31m' ; GREEN='\033[0;32m' ; NC='\033[0m' # No Color
die() { echo -e "${RED}$1${NC}" >&2 ; echo ; exit 1 ; } # terminate script
info() {  echo -e "${GREEN}$1${NC}" >&2 ; }
info "$(date)"

RM_DIR=$(realpath "RepeatMasker")
TERM="40674"

# Directories (change)
OUTPUT_DIR=$(realpath "test/search_dfam/${TERM}")

create_output_dir() {
	test -d "${OUTPUT_DIR}" && rm -r "${OUTPUT_DIR}"
	mkdir "${OUTPUT_DIR}" || die "Failed creating ${OUTPUT_DIR}"
}

print_lineage() {
    python ${RM_DIR}/famdb.py lineage ${TERM} -ad || die "could not fetch lineage tree"
    python ${RM_DIR}/famdb.py lineage -ad --format totals ${TERM} || die "could not fetch the total number of records"
}

print_record_stats() {
    python ${RM_DIR}/famdb.py families --format summary ${TERM} -ad > ${OUTPUT_DIR}/summary || die "could not fetch record information"
    cat ${OUTPUT_DIR}/summary | tail -n +37 | awk -F'=' '{print $2}' > ${OUTPUT_DIR}/lens || die "could not fetch record lens"

    echo -e "count\tmin\tmax\tmode\tmean\tstddev" > ${OUTPUT_DIR}/stats
    cat ${OUTPUT_DIR}/lens | datamash count 1 min 1 max 1 mode 1 mean 1 sstdev 1 >> ${OUTPUT_DIR}/stats || die "cold not calculate record lens stats"
    cat ${OUTPUT_DIR}/stats
}

info "searching dfam database for term:${TERM}"

info "9606 Homo sapiens(0) [52] means NCBI Taxonomy ID for Homo sapiens, Scientific name, Partition number in FamDB(0), [52] repeat families assigned specifically to Homo sapiens"

create_output_dir
print_lineage
print_record_stats
