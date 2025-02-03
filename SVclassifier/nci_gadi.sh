#!/bin/bash
#PBS -P ox63
#PBS -N a_SVclassifier
#PBS -l ncpus=16
#PBS -l mem=4GB
#PBS -l walltime=1:00:00
#PBS -l wd

RED='\033[0;31m' ; GREEN='\033[0;32m' ; NC='\033[0m' # No Color
die() { echo -e "${RED}$1${NC}" >&2 ; echo ; exit 1 ; } # terminate script
info() {  echo -e "${GREEN}$1${NC}" >&2 ; }
info "$(date)"

module load bcftools || die "could not load bcftools module"
module load parallel || die "could not load parallel module"
source svtools/bin/activate || die "could not activate svtools venv"

./SVclassifier/run_classifier.sh || die "could not run_classifier.sh"
