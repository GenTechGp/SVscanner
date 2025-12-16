#!/bin/bash

module use -a /g/data/if89/apps/modulefiles
module load RepeatMasker/4.2.0
module load parallel
module load python3/3.9.2
module load pythonlib/3.9.2
module load bcftools/1.22
module load htslib/1.22.1
module load nextflow/25.04.6

nextflow run main.nf --vcf ../test/HG002_subset_mini/HG002_subset_mini.vcf.gz --ref /g/data/kr68/genome/hg38.analysisSet.fa
