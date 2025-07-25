#!/bin/bash
#PBS -P project
#PBS -N test
#PBS -l storage=gdata/kr68+gdata/ox63+gdata/if89
#PBS -l ncpus=48
#PBS -l mem=128GB
#PBS -l walltime=20:00:00
#PBS -l wd


die() { echo -e "$1" >&2 ; echo ; exit 1 ; } # terminate script

#export MODULEPATH=/g/data/ox63/install/modules:$MODULEPATH
module use /g/data/ox63/install/modules || die "could not use module path"
module load svclass/1.0 || die "could not load svclass module"

echo "--out $OUTPUT_DIR --vcf $VCF --ref $REF"
svclass --out $OUTPUT_DIR --vcf $VCF --ref $REF || die "could not run svclass"
