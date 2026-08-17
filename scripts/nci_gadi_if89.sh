#!/bin/bash
#PBS -P project
#PBS -N test
#PBS -l storage=gdata/if89
#PBS -l ncpus=48
#PBS -l mem=128GB
#PBS -l walltime=06:00:00
#PBS -l wd

die() { echo -e "$1" >&2 ; echo ; exit 1 ; } # terminate script

module use -a /g/data/if89/apps/modulefiles
module load SVscanner

echo "--vnv $VNV --out $OUT --vcf $VCF --ref $REF"
svscanner --version || die "could not run svscanner"
svscanner --out $OUT --vcf $VCF --ref $REF || die "could not run svscanner"

echo "svscanner finished successfully"

