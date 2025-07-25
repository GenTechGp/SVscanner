#!/bin/bash

# set -x
RED='\033[0;31m' ; GREEN='\033[0;32m' ; NC='\033[0m' # No Color
die() { echo -e "${RED}$1${NC}" >&2 ; echo ; exit 1 ; } # terminate script
info() {  echo -e "${GREEN}$1${NC}" >&2 ; }
info "$(date)"

TRF_BINARY=$(realpath "trf409.linux64")
REPEAT_MASKER=$(realpath "RepeatMasker/RepeatMasker")

NUM_SEQ=10
SEED=100
# Directories (change)
OUTPUT_DIR=$(realpath "test/bench")

# seq_lens=("50 100 200 300 400 500 1000 5000 10000 20000 30000 40000 50000 100000 200000 300000 400000 500000 600000 700000 800000 900000 1000000")
seq_lens=("50 100 200 300 400 500")

create_output_dir() {
	test -d "${OUTPUT_DIR}" && rm -r "${OUTPUT_DIR}"
	mkdir "${OUTPUT_DIR}" || die "Failed creating ${OUTPUT_DIR}"
}

generate_seq() {
    for i in ${seq_lens} ; do
        info "running with ${i}"
        python src/random_seq_tr.py -l ${i} -n ${NUM_SEQ} -o "${OUTPUT_DIR}/${i}_${NUM_SEQ}.fa" --seed ${SEED} || die "generate seq failed"
    done
}

bench_trf() {
    for i in ${seq_lens} ; do
        info "running with ${i}"
        /usr/bin/time -v ${TRF_BINARY} "${OUTPUT_DIR}/${i}_${NUM_SEQ}.fa" 2 7 7 80 10 50 500 -h -ngs > "${OUTPUT_DIR}/${i}_${NUM_SEQ}.dat" || die "TRF failed"
    done
}

bench_rm() {
    for i in ${seq_lens} ; do
        info "running with ${i}"
        /usr/bin/time -v ${REPEAT_MASKER} "${OUTPUT_DIR}/${i}_${NUM_SEQ}.fa" -pa 2 -html -gff -dir ${OUTPUT_DIR}  || die "RM failed"
    done
}

create_output_dir
generate_seq
bench_trf
bench_rm
