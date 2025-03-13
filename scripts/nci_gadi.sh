#!/bin/bash
#PBS -P ox63
#PBS -N a_SVclassifier
#PBS -l ncpus=48
#PBS -l mem=128GB
#PBS -l walltime=1:00:00
#PBS -l wd

die() { echo -e "$1" >&2 ; echo ; exit 1 ; } # terminate script
info() {  echo -e "$1" >&2 ; }
info "$(date)"

module load parallel || die "could not load parallel module"
source svtools/bin/activate || die "could not activate svtools venv"

./scripts/run_classifier.sh || die "could not run_classifier.sh"
