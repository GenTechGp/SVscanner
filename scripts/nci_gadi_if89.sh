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
module load RepeatMasker/4.2.0
module load parallel
module load python3/3.9.2
module load pythonlib/3.9.2
module load bcftools/1.22
module load htslib/1.22.1

echo "--vnv $VNV --out $OUT --vcf $VCF --ref $REF"
./scripts/run_workflow.sh --out $OUT --vcf $VCF --ref $REF || die "could not run workflow"
