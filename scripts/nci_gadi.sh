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

module load parallel || die "could not load parallel module"
cd "/scratch/ox63/hm4078/SVtoolkit"
source "svtools/bin/activate" || die "could not activate svtools venv"

echo "--output_dir $OUTPUT_DIR --sample $SAMPLE --sv_vcf $SV_VCF --ref_fasta $REF"
./scripts/run_classifier.sh --output_dir $OUTPUT_DIR --sample $SAMPLE --sv_vcf $SV_VCF --ref_fasta $REF || die "could not run_classifier.sh"

## qsub -v OUTPUT_DIR=/scratch/ox63/hm4078/sv_out,SAMPLE=sample_01,SV_VCF=/scratch/ox63/hm4078/SVtoolkit/test/HG002_subset_mini/HG002_subset_mini.vcf.gz /scratch/ox63/hm4078/SVtoolkit/scripts/nci_gadi.sh