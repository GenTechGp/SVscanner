# SVscanner
![Illustration](/images/svscanner_logo.png)

A workflow to annotate tandem repeats and mobile elements within structural variants (SVs) using [Tandem Repeat Finder (TRF)](https://github.com/Benson-Genomics-Lab/TRF) and [RepeatMasker (RM)](https://github.com/Dfam-consortium/RepeatMasker).

## Overview

### Inputs:
1. A structural variant (SV) file in VCF format
2. A reference genome
### Process:
1. Flanking sequences are extracted around each SV to form a query sequence (details about the [extraction process](docs/sv_extraction.md)):
    - query_sequence = left_flank + SV + right_flank
2. These query sequences are annotated using TRF and RM. RM's alignments are also chained into one alignment per element, for the diagrams.

### Outputs:
1. A VCF file with repeat annotations embedded (details about [VCF INFO tags](docs/repeat_annotation_VCF_tags.md) and [annotation process](docs/repeat_annotation_steps.md)).
2. Diagrams visualizing Repeat/SV annotations, including where an element carries sequence its consensus does not have, e.g. an SVA VNTR expansion (details about the [diagram file format](docs/Repeat-SV_diagram.md)).
3. Histograms and density plots summarizing Repeat/SV annotations.

![Illustration](/images/SVscanner_workflow.png)

### Notes:
1. STRchive dataset for genome (in extended BED format) is optional input for `repeat_annotation.py`([hg19](https://strchive.org/_astro/STRchive-disease-loci.hg19.DWACvaXd.bed), [hg38](https://strchive.org/_astro/STRchive-disease-loci.hg38.DR-UScgX.bed), [T2T-chm13](https://strchive.org/_astro/STRchive-disease-loci.T2T-chm13.Cm-HAugT.bed))
2. Usage of [run_workflow.sh](docs/Commands.md#run_workflowsh)
3. Usage of workflow components: [extract_sv.py](docs/Commands.md#extrract_svpy), [repeat_annotation.py](docs/Commands.md#repeat_annotatoinpy), [generate_plots.py](docs/Commands.md#generate_plotpy)

## Running SVscanner

SVscanner is distributed as a container image, published to the GitHub Container Registry
on every release. It carries everything the workflow needs: RepeatMasker 4.2.6 with the
Dfam 4.0 database, TRF, RMBlast, bcftools and the Python environment. Nothing needs
cloning, compiling or downloading besides the image itself and your reference genome.

Replace `X.Y.Z` below with a release number (see the
[releases](https://github.com/GenTechGp/SVscanner/releases) page). Pin an exact version
for any analysis whose results you intend to keep.

### Docker

```
docker run --rm \
  -v /path/to/reference:/ref:ro \
  -v "$PWD":/data \
  ghcr.io/gentechgp/svscanner:X.Y.Z \
  svscanner --vcf /data/input.vcf.gz --ref /ref/hg38.fa --out /data/svscanner_out
```

### Singularity / Apptainer

Convert the image to a `.sif` file once, then run that file:

```
singularity pull svscanner_X.Y.Z.sif docker://ghcr.io/gentechgp/svscanner:X.Y.Z

singularity exec \
  --bind /path/to/reference:/ref:ro \
  svscanner_X.Y.Z.sif \
  svscanner --vcf input.vcf.gz --ref /ref/hg38.fa --out svscanner_out
```

`svscanner` takes the same arguments everywhere; run `svscanner --help` for the full list,
or see [run_workflow.sh](docs/Commands.md#run_workflowsh). More detail, including Nextflow
and building the image yourself, is in [docs/docker.md](docs/docker.md).

### Quick test run

The image includes a 54-record test VCF. With a human reference genome this takes a few
minutes:

```
singularity exec --bind /path/to/reference:/ref:ro svscanner_X.Y.Z.sif \
  svscanner --vcf /opt/svscanner/test/HG002_subset_mini/HG002_subset_mini.vcf.gz \
            --ref /ref/hg38.fa --out ./svscanner_test
```

### Dfam database

The image's built-in Dfam 4.0 holds all curated consensus sequences, so with the default
RMBlast search engine `--species` alone covers human or any other species with curated
Dfam families. Beyond that:

- `--dfam_extra DIR` adds partitions of the same Dfam release, e.g. the curated HMM
  partitions needed for `--engine hmmer`;
- `--dfam_db DIR` replaces the built-in database with a complete one (base file plus
  partitions), e.g. a different Dfam release;
- `--uncurated` uses Dfam's uncurated families as well, for species with few curated ones.
  Their partitions are not built in, so add them with `--dfam_extra`.

See [Dfam databases](docs/Commands.md#dfam-databases).

### On an HPC cluster (example: NCI Gadi)

Compute nodes usually have no internet access, so pull the `.sif` once on a login or data
transfer node (on Gadi, a `copyq` job or `gadi-dm.nci.org.au`) and keep it on project
storage. Then submit jobs that run it:

```bash
#!/bin/bash
#PBS -P [your_project]
#PBS -l storage=gdata/[project_holding_your_data]
#PBS -l ncpus=48
#PBS -l mem=128GB
#PBS -l walltime=06:00:00
#PBS -l jobfs=10GB
#PBS -l wd

module load singularity

singularity exec \
  --bind "$(dirname "$VCF")":/vcfdir:ro \
  --bind "$(dirname "$REF")":/refdir:ro \
  --bind "$(dirname "$OUT")":/out \
  "$SIF" \
  svscanner --vcf "/vcfdir/$(basename "$VCF")" \
            --ref "/refdir/$(basename "$REF")" \
            --out "/out/$(basename "$OUT")" \
            --nthread "$PBS_NCPUS"
```

```
qsub -N [job_name] -v SIF=[path_to_sif],OUT=[path_to_out],VCF=[path_to_vcf],REF=[path_to_ref] [job_script]
```

A fuller version with input checks is in [scripts/nci_gadi_singularity.sh](scripts/nci_gadi_singularity.sh).

### Notes

1. **Thread and job allocation.** Inside a container, pass `--nthread` explicitly (as above):
   SVscanner otherwise falls back to `nproc`, which can report the whole host. RepeatMasker
   is run as `--njob` parallel processes of `--nthread / --njob` threads each (`-pa`),
   defaulting to one process per thread. Examples for `--nthread 48`:

    | flags | RepeatMasker processes | threads each (`-pa`) | total |
    |---|---|---|---|
    | *(none — default)* | 48 | 1 | 48 |
    | `--njob 24` | 24 | 2 | 48 |
    | `--njob 12` | 12 | 4 | 48 |
    | `--njob 8` | 8 | 6 | 48 |

    Fewer, wider processes (`--njob 12`) suit a small number of very long query sequences;
    the default suits many short ones. `--njob` is capped at `--nthread`, and each process
    always gets at least one thread.
2. Each run spends about 20 seconds preparing RepeatMasker's search library before
   annotating (`Library preparation took N seconds` in the log). The library is built
   inside the run's output directory and removed at the end, so nothing is written into the
   image or your home directory, and concurrent runs never share it.

## Releasing

The version number lives only in the [VERSION](VERSION) file, which `run_workflow.sh` reads at startup. CI and a pre-push hook both reject a release tag `vX.Y.Z` published on a commit that does not say `X.Y.Z`. See [docs/releasing.md](docs/releasing.md).

## Bug Reports

Please report/request any issues/features via [GitHub Issues](https://github.com/GenTechGp/SVscanner/issues).