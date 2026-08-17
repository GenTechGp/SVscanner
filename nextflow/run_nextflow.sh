#!/bin/bash

die() { echo -e "$1" >&2 ; echo ; exit 1 ; } # terminate script

module use -a /g/data/if89/apps/modulefiles
module load SVscanner || die "could not load SVscanner module"
module load nextflow/25.04.6 || die "could not load nextflow module"

# main.nf copies ${projectDir}/../src and runs that, not the module's own copy,
# so the module's svscanner wrapper is bypassed and its pinned Python env has to
# be put on PATH here. Do not use pythonlib/3.9.2: it ships pysam 0.23.3, whose
# bundled htslib resolves INFO/END against INFO/SVLEN differently (see
# get_sv_end in src/extract_sv.py). PYTHONPATH is unset for the same reason the
# module wrapper unsets it - venvs do not ignore it.
SVSCANNER_INSTALL=$(dirname "$(dirname "$(readlink -f "$(command -v svscanner)")")")
[[ -x "${SVSCANNER_INSTALL}/venv/bin/python3" ]] || die "could not locate the SVscanner venv under ${SVSCANNER_INSTALL}"
unset PYTHONPATH
export PATH="${SVSCANNER_INSTALL}/venv/bin:${PATH}"

nextflow run main.nf --vcf ../test/HG002_subset_mini/HG002_subset_mini.vcf.gz --ref /g/data/kr68/genome/hg38.analysisSet.fa
