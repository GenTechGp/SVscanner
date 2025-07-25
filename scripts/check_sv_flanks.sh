#!/bin/bash

# set -x
RED='\033[0;31m' ; GREEN='\033[0;32m' ; NC='\033[0m' # No Color
die() { echo -e "${RED}$1${NC}" >&2 ; echo ; exit 1 ; } # terminate script
info() {  echo -e "${GREEN}$1${NC}" >&2 ; }
info "$(date)"

SEED=42

OUTPUT_DIR=$(realpath "test/check_sv_flanks")

BASE_REF="test/sim_ref/base_ref/base_ref.fa"
SV_TREATED_REF="test/sim_ref/visor_hack/h1.fa"

FLANKS_FASTA="SVtoolkit_output_sim_ref/sim_ref/extract_sv_flanks_out/0.fa"
INS_FASTA=${OUTPUT_DIR}/ins.fasta
NON_INS_FASTA=${OUTPUT_DIR}/non_ins.fasta

creat_output_dir() {
	test -d "${OUTPUT_DIR}" && rm -r "${OUTPUT_DIR}"
    mkdir ${OUTPUT_DIR}
}

separate_fasta() {
    local fasta_file="$1"
    local ins_file="$2"
    local rest_file="$3"

    awk -v ins_file="$ins_file" -v rest_file="$rest_file" '
    BEGIN { out="" }
    /^>/ {
        if ($0 ~ /^>INS/) {
            out=ins_file
        } else {
            out=rest_file
        }
    }
    out != "" { print > out }
    ' "$fasta_file"
}


# ref $1 reads $2 
align_reads() {
    # minimap2 -cx map-ont "$1" -t32 --secondary=no "$2" -o ${OUTPUT_DIR}/1.paf
    # minimap2 -cx map-ont "$1" -t32 --secondary=no "$2" 2> /dev/null | awk '{print $1,$10,$11,$10/$11}'
    minimap2 -cx map-ont "$1" -t32 --secondary=no "$2" 2> /dev/null | awk '{
        for (i=1; i<=NF; i++) {
            if ($i ~ /^cg:/) {
                print $1,$10,$11,$10/$11, $i;
                break;
            }
        }
    }'
    # minimap2 -cx map-ont "$1" -t32 --secondary=no "$2" | awk '{print $10/$11}' | datamash mean 1 sstdev 1 q1 1 median 1 q3 1 count 1 || die "minimap2 accuracy metric failed"
    # minimap2 -ax map-ont ${BASE_REF} -t32 --secondary=no ${SIM_READS} > ${SAM} || die "minimap2 failed"
    # samtools sort ${SAM} -o ${BAM} && samtools index ${BAM} || die "samtools sort and index failed"
}

creat_output_dir
separate_fasta ${FLANKS_FASTA} ${INS_FASTA} ${NON_INS_FASTA}
align_reads ${BASE_REF} ${NON_INS_FASTA}
align_reads ${SV_TREATED_REF} ${INS_FASTA}

# make_consensus
