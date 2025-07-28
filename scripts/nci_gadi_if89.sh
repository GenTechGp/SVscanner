#!/bin/bash
#PBS -P kr68
#PBS -N test
#PBS -l storage=gdata/kr68+gdata/ox63+gdata/if89
#PBS -l ncpus=48
#PBS -l mem=128GB
#PBS -l walltime=05:00:00
#PBS -l wd


die() { echo -e "$1" >&2 ; echo ; exit 1 ; } # terminate script

echo "--vnv $VNV --out $OUT --vcf $VCF --ref $REF"
./scripts/run_workflow_if89.sh --vnv $VNV --out $OUT --vcf $VCF --ref $REF || die "could not run svclass"