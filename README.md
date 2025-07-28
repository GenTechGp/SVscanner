# SVscanner

*A workflow to annotate and analyse structural variants (SV) with Repeat information.* 

## Overview

A workflow for annotating tandem repeats and mobile elements identified by Tandem Repeat Finder ([TRF](https://github.com/Benson-Genomics-Lab/TRF)) and RepeatMasker ([RM](https://github.com/Dfam-consortium/RepeatMasker)) within SVs. The workflow takes SV information (VCF) and reference genome as inputs. First, flanking sequences around the SV are extracted. Then extracted sequences are annotated using TRF and RM. The workflow has three main outputs.

1. SV VCF file with repeat information annotated.
2. Diagrams with Repeat/SV annotations ([details](docs/Repeat-SV_diagram.md)).
3. Histograms and density plots summarising Repeat/SV annotations.

![Illustration](/images/SVscanner_workflow.png)

- STRchive dataset for genome (in extended BED format) is optional input for `repeat_annotation.py`([hg19](https://strchive.org/_astro/STRchive-disease-loci.hg19.DWACvaXd.bed), [hg38](https://strchive.org/_astro/STRchive-disease-loci.hg38.DR-UScgX.bed), [T2T-chm13](https://strchive.org/_astro/STRchive-disease-loci.T2T-chm13.Cm-HAugT.bed))

- Details about [SV extraction](docs/extract_sv_and_seq_consensus_documentation.md)
- Details about [repeat annotation](docs/repeat_annotation_steps.md)
- Usage of [extract_sv.py](docs/Commands.md#extrract_svpy), [repeat_annotation.py](docs/Commands.md#repeat_annotatoinpy), [generate_plots.py](docs/Commands.md#generate_plotpy)

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

3. Install `RepeatMasker` if not available

- Follow instructions [here](docs/install_rm.md)


4. Check if the following tools are available. If not install them.
 - `bcftools` (v1.21 or above recommended)
 - `bgzip` (v1.21 or above recommended)
 - `tabix` (v1.21 or above recommended)
 - `trf`
 - `parallel`
 
```
./scripts/install_tool.sh [tools to be installed]
e.g. ./scripts/install_tool.sh bcftools htslib
```
*`htslib` installs both `bgzip` and `tabix`

5. Check if the workflow is working
```
./scripts/run_workflow.sh

```

### Quick example run

The following example processes the `test/HG002_subset_mini.vcf.gz` dataset (100 records).

It will take about 10 minutes. The majority of time is taken by the RepeatMasker step. Please provide the path to human genome after `--ref` argument.

```
./scripts/run_workflow.sh --vcf test/HG002_subset_mini/HG002_subset_mini.vcf.gz --ref [human genome] --out test/output
```

### NCI setup
TBC
