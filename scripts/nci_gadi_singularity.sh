#!/bin/bash
#PBS -P project
#PBS -N svscanner
#PBS -l storage=gdata/project
#PBS -l ncpus=48
#PBS -l mem=128GB
#PBS -l walltime=06:00:00
#PBS -l jobfs=10GB
#PBS -l wd

# Run SVscanner from its Singularity image on NCI Gadi (or any PBS cluster).
#
#   qsub -P [project] -l storage=gdata/[project_holding_your_data] \
#        -v SIF=/path/svscanner_X.Y.Z.sif,VCF=/path/in.vcf.gz,REF=/path/ref.fa,OUT=/path/out \
#        scripts/nci_gadi_singularity.sh
#
# Optional: PREFIX=... (output file prefix), SPECIES=... (default human).
# Pull the image beforehand on a node with internet access:
#   singularity pull svscanner_X.Y.Z.sif docker://ghcr.io/gentechgp/svscanner:X.Y.Z

die() { echo -e "$1" >&2; exit 1; }

module load singularity || die "could not load the singularity module"

[[ -n "${SIF}" && -f "${SIF}" ]] || die "SIF is not set or not found: '${SIF}' - submit with -v SIF=/path/to/svscanner.sif"
[[ -n "${VCF}" && -f "${VCF}" ]] || die "VCF is not set or not found: '${VCF}' - submit with -v VCF=/path/to/input.vcf.gz"
[[ -n "${REF}" && -f "${REF}" ]] || die "REF is not set or not found: '${REF}' - submit with -v REF=/path/to/reference.fasta"
[[ -n "${OUT}" ]] || die "OUT is not set - submit with -v OUT=/path/to/output_dir"
SPECIES=${SPECIES:-human}

# OUT is the run's own output directory (created by svscanner itself); its parent is what
# gets bound, so svscanner can create OUT under it. A bind source must already exist.
OUT_BASE=$(dirname "${OUT}")
mkdir -p "${OUT_BASE}" || die "failed to create ${OUT_BASE}"

echo "Image: ${SIF}"
echo "VCF: ${VCF}"
echo "Reference: ${REF}"
echo "Output: ${OUT}"

singularity exec \
  --bind "$(dirname "${VCF}"):/vcfdir:ro" \
  --bind "$(dirname "${REF}"):/refdir:ro" \
  --bind "${OUT_BASE}:/out" \
  "${SIF}" \
  svscanner --vcf "/vcfdir/$(basename "${VCF}")" \
            --ref "/refdir/$(basename "${REF}")" \
            --out "/out/$(basename "${OUT}")" \
            --species "${SPECIES}" \
            --nthread "${PBS_NCPUS}" \
            ${PREFIX:+--prefix "${PREFIX}"} \
  || die "svscanner failed - see the output above"

echo "svscanner finished successfully"
