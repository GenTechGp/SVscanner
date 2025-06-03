#!/bin/bash

# set -x
RED='\033[0;31m' ; GREEN='\033[0;32m' ; NC='\033[0m' # No Color
die() { echo -e "${RED}$1${NC}" >&2 ; echo ; exit 1 ; } # terminate script
info() {  echo -e "${GREEN}$1${NC}" >&2 ; }
info "$(date)"

SEED=42

OUTPUT_DIR=$(realpath "test/test_visor")

BASE_REF=$(realpath "test/databases/test_visor/base_ref.fa")
BED_FILE=$(realpath "test/databases/test_visor/visor_hack.bed")
SV_REF_CONTROL=$(realpath "test/databases/test_visor/sv_treated_control.fa")
SV_TREATED_REF=${OUTPUT_DIR}/sv_treated_ref.fa

VISOR_VENV_PATH="/data/hiruna/SVtoolkit/VISOR-1.1.2.1/visor"

create_output_dir() {
	test -d "${OUTPUT_DIR}" && rm -r "${OUTPUT_DIR}"
    mkdir ${OUTPUT_DIR}
}

simulate_sv_using_visor() {
    source "${VISOR_VENV_PATH}/bin/activate"
    VISOR HACk -b ${BED_FILE} -g ${BASE_REF} -o ${OUTPUT_DIR}/visor_hack || die "visor hack failed"
    deactivate
    awk '/^>/ {if (seq) print seq; print; seq=""; next} {seq=seq $0} END {if (seq) print seq}' ${OUTPUT_DIR}/visor_hack/h1.fa > ${SV_TREATED_REF} || die "awk after visor failed"
    # diff -q ${SV_REF_CONTROL} ${SV_TREATED_REF} || die "diff failed, visor output is not as expected"
}

create_output_dir
simulate_sv_using_visor

info "success"