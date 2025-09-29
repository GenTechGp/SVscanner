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
2. These query sequences are annotated using TRF and RM.

### Outputs:
1. A VCF file with repeat annotations embedded (details about [VCF INFO tags](docs/repeat_annotation_VCF_tags.md) and [annotation process](docs/repeat_annotation_steps.md)).
2. Diagrams visualizing Repeat/SV annotations (details about the [diagram file format](docs/Repeat-SV_diagram.md)).
3. Histograms and density plots summarizing Repeat/SV annotations.

![Illustration](/images/SVscanner_workflow.png)

### Notes:
1. STRchive dataset for genome (in extended BED format) is optional input for `repeat_annotation.py`([hg19](https://strchive.org/_astro/STRchive-disease-loci.hg19.DWACvaXd.bed), [hg38](https://strchive.org/_astro/STRchive-disease-loci.hg38.DR-UScgX.bed), [T2T-chm13](https://strchive.org/_astro/STRchive-disease-loci.T2T-chm13.Cm-HAugT.bed))
2. Usage of [run_workflow.sh](docs/Commands.md#run_workflowsh)
3. Usage of workflow components: [extract_sv.py](docs/Commands.md#extrract_svpy), [repeat_annotation.py](docs/Commands.md#repeat_annotatoinpy), [generate_plots.py](docs/Commands.md#generate_plotpy)

## Installation 

1. Clone the repository

```
git clone git@github.com:KCCGGenomeTechLab/SVscanner.git
```

2. Set up Virtual Environment and install required packages. Tested with `python 3.8` and should work with higher versions as well.

```
cd SVscanner
python3 -m venv svscanner
source svscanner/bin/activate 
pip install --upgrade pip
pip install -r requirements.txt
```

3. Follow instruction listed [here](docs/install_rm.md) to install `RepeatMasker` if not available already.

4. Check if the following tools are available. If not install them.
 - `bcftools` (v1.21 or above recommended)
 - `bgzip` (v1.21 or above recommended)
 - `tabix` (v1.21 or above recommended)
 - `trf`
 - `parallel`
 
```
./scripts/install_tools.sh [tools to be installed]
e.g. ./scripts/install_tools.sh bcftools htslib
```
*`htslib` installs both `bgzip` and `tabix`

5. Check if the workflow is working
```
./scripts/run_workflow.sh
```

## Quick example run

The following example processes the `test/HG002_subset_mini.vcf.gz` dataset (100 records).

It will take about 10 minutes. The majority of time is taken by the RepeatMasker step. Please provide the path to human genome after `--ref` argument.

```
./scripts/run_workflow.sh --vcf test/HG002_subset_mini/HG002_subset_mini.vcf.gz --ref [human genome] --out test/output
```

## SVscanner on NCI Gadi (Project if89 users)
Until `SVscanner` is available as an official NCI module, users on Gadi can run the workflow using the following steps:

1. Setup:
```
git clone git@github.com:KCCGGenomeTechLab/SVscanner.git
cd SVscanner
```

2. Submit job:
```
cd SVscanner
qsub -N [job_name] -v OUT=[path_to_out],VCF=[path_to_vcf],REF=[path_to_ref] scripts/nci_gadi_if89.sh
```

### Notes
1. Edit `scripts/nci_gadi_if89.sh` and insert your NCI project code.
2. No need to install `TRF`, `RepeatMasker`, `bcftools`, `bgzip`, `tabix`, or `GNU parallel`. These are all available as pre-installed NCI modules and are loaded by the workflow.
3. Database Dfam `3.9`; FamDB Format `2.0`; Partition `7` [dfam39_full.7.h5]: Mammalia (57 GB) is used with RepeatMasker module (`4.2.0`) [more info](https://www.dfam.org/releases/Dfam_3.9/families/FamDB/README.txt) 
4. To pass additional arguments to the workflow, edit `scripts/nci_gadi_if89.sh` as needed — it forwards parameters to `scripts/run_workflow.sh`.
5. A simple workflow runtime benchmark done on NCI Gadi ([link](docs/nci_benchmark.md))

## Bug Reports

Please report/request any issues/features via [GitHub Issues](https://github.com/KCCGGenomeTechLab/SVscanner/issues).