#!/bin/bash
#PBS -P ox63
#PBS -N svclassifier_1
#PBS -l storage=gdata/kr68+gdata/ox63+gdata/if89
#PBS -l ncpus=48
#PBS -l mem=128GB
#PBS -l walltime=20:00:00
#PBS -l wd
#PBS -V

die() { echo -e "$1" >&2 ; echo ; exit 1 ; } # terminate script
info() {  echo -e "$1" >&2 ; }
info "$(date)"

export MODULEPATH=/g/data/ox63/install/modules:$MODULEPATH
module load svclass/1.0 || die "could not load svclass module"

echo "--output_dir $OUTPUT_DIR --sample $SAMPLE --sv_vcf $SV_VCF --ref_fasta $REF"
svclass --output_dir $OUTPUT_DIR --sample $SAMPLE --sv_vcf $SV_VCF --ref_fasta $REF || die "could not run svclass"