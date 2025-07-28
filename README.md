# SVscanner
![Illustration](/images/svcanner_logo.png)

This workflow annotates tandem repeats and mobile elements within structural variants (SVs) using [Tandem Repeat Finder (TRF)](https://github.com/Benson-Genomics-Lab/TRF) and [RepeatMasker (RM)](https://github.com/Dfam-consortium/RepeatMasker).

## Overview

### Inputs:
- A structural variant (SV) file in VCF format
- A reference genome
### Process:
1. Flanking sequences are extracted around each SV to form a query sequence:
    - query_sequence = left_flank + SV + right_flank
2. These query sequences are annotated using TRF and RM.

### Outputs:
1. A VCF file with repeat annotations embedded.
2. Diagrams visualizing Repeat/SV annotations ([details](docs/Repeat-SV_diagram.md)).
3. Histograms and density plots summarizing Repeat/SV annotations.

![Illustration](/images/SVscanner_workflow.png)

### Notes:
1. STRchive dataset for genome (in extended BED format) is optional input for `repeat_annotation.py`([hg19](https://strchive.org/_astro/STRchive-disease-loci.hg19.DWACvaXd.bed), [hg38](https://strchive.org/_astro/STRchive-disease-loci.hg38.DR-UScgX.bed), [T2T-chm13](https://strchive.org/_astro/STRchive-disease-loci.T2T-chm13.Cm-HAugT.bed))

2. Details about [SV extraction](docs/extract_sv_and_seq_consensus_documentation.md)
3.  Details about [repeat annotation](docs/repeat_annotation_steps.md)
4. Usage of [extract_sv.py](docs/Commands.md#extrract_svpy), [repeat_annotation.py](docs/Commands.md#repeat_annotatoinpy), [generate_plots.py](docs/Commands.md#generate_plotpy)

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
./scripts/install_tool.sh [tools to be installed]
e.g. ./scripts/install_tool.sh bcftools htslib
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

## Gadi | NCI setup
[TBC](https://nci.org.au/news-events/events/introduction-gadi-4)
